const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");
const matrix = @import("../../../matrix.zig");

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");
const float = @import("../../../float.zig");

const linalg = @import("../../../linalg.zig");

/// Performs a rank-1 update of a symmetric matrix defined as:
///
/// ```zig
/// A = alpha * x * xᵀ + A,
/// ```
///
/// where `alpha` is a numeric, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` symmetric matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.syr(layout: matrix.Layout, uplo: matrix.Uplo, n: usize, alpha: Al, x: [*]const X, incx: isize, a: [*]A, lda: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the symmetric matrix `A` is used.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * n`. On
///   return, contains the result of the operation.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
/// * `opts`: Optional parameters:
///   * `num_threads` (`usize = 0`): Number of threads to spawn:
///     * `0`: automatic. The thread count is derived from `m * n` and
///       `parallel_threshold`:
///       ```zig
///       threads = max(1, min(std.Thread.getCpuCount(), options.max_threads, (n * n) / parallel_threshold))
///       ```
///     * 1: force serial execution. parallel_threshold is ignored.
///     * N >= 2: use exactly N threads, clamped by
///       std.Thread.getCpuCount() and options.max_threads as a hard safety
///       ceiling. parallel_threshold is ignored.
///   * parallel_threshold (usize = 2_097_152 / @sizeOf(meta.Child(A))):
///     Minimum number of matrix elements (`n * n`) required to trigger
///     multithreaded execution.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is 0, or if `lda` is less
///   than `max(1, n)`.
pub fn syr(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    a: anytype,
    lda: usize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 2_097_152 / @sizeOf(meta.Child(@TypeOf(a))),
    },
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var A: type = @TypeOf(a);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.blas.syr: alpha must be a numeric, x must be a many-item pointer to numerics, and a must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ta: " ++ @typeName(A) ++ "\n");

    X = meta.Child(X);
    A = meta.Child(A);

    if (lda < int.max(1, n) or incx == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();

    // Quick return if possible.
    if (n == 0 or numeric.eq(alpha, 0))
        return;

    if (opts.num_threads == 1)
        return k_syr(eff_uplo, n, alpha, x, incx, a, lda);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, (n * n) / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);
    num_threads = int.min(num_threads, n);

    if (num_threads <= 1)
        return k_syr(eff_uplo, n, alpha, x, incx, a, lda);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);
    num_threads = int.min(num_threads, n);

    if (num_threads <= 1)
        return k_syr(eff_uplo, n, alpha, x, incx, a, lda);

    const Worker = struct {
        fn execute(
            worker_uplo: matrix.Uplo,
            worker_n: usize,
            worker_alpha: Al,
            worker_a: [*]A,
            worker_lda: usize,
            worker_x: [*]const X,
            worker_incx: isize,
            counter: *std.atomic.Value(usize),
            comptime worker_tile_size: comptime_int,
            worker_num_tiles: usize,
        ) void {
            while (true) {
                const idx = counter.fetchAdd(1, .monotonic);

                if (idx >= worker_num_tiles) // When all tiles have been assigned, break.
                    break;

                // Map 1D atomic index to 2D upper triangular coordinates (tile_i, tile_j) using triangular numbers.
                var tile_j = numeric.cast(usize, (float.sqrt(1.0 + 8.0 * numeric.cast(f64, idx)) - 1.0) / 2.0);

                while (tile_j * (tile_j + 1) / 2 > idx)
                    tile_j -= 1;

                while ((tile_j + 1) * (tile_j + 2) / 2 <= idx)
                    tile_j += 1;

                const tile_i = idx - tile_j * (tile_j + 1) / 2;

                const phys_r = if (worker_uplo == .upper) tile_i else tile_j;
                const phys_c = if (worker_uplo == .upper) tile_j else tile_i;

                const r_start = phys_r * worker_tile_size;
                const c_start = phys_c * worker_tile_size;
                const r_len = int.min(worker_tile_size, worker_n - r_start);
                const c_len = int.min(worker_tile_size, worker_n - c_start);

                if (tile_i == tile_j) {
                    // Diagonal tile
                    var local_x: [worker_tile_size]X = undefined;

                    const px = if (worker_incx == 1)
                        worker_x + r_start
                    else blk: {
                        @import("../level1/copy.zig").k_copy(
                            r_len,
                            worker_x + numeric.cast(usize, if (worker_incx > 0)
                                numeric.cast(isize, r_start) * worker_incx
                            else
                                (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incx),
                            worker_incx,
                            @as([*]X, &local_x),
                            1,
                        );

                        break :blk @as([*]const X, &local_x);
                    };

                    k_syr(
                        worker_uplo,
                        r_len,
                        worker_alpha,
                        px,
                        1,
                        worker_a + r_start + c_start * worker_lda,
                        worker_lda,
                    );
                } else {
                    const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(X, numeric.Mul(Al, X), A)) orelse 2);

                    // Off-diagonal tile
                    var local_x_r: [worker_tile_size]X = undefined;
                    var local_x_c: [worker_tile_size]X = undefined;

                    const px_r = if (worker_incx == 1)
                        worker_x + r_start
                    else blk: {
                        @import("../level1/copy.zig").k_copy(
                            r_len,
                            worker_x + numeric.cast(usize, if (worker_incx > 0)
                                numeric.cast(isize, r_start) * worker_incx
                            else
                                (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incx),
                            worker_incx,
                            @as([*]X, &local_x_r),
                            1,
                        );

                        break :blk @as([*]const X, &local_x_r);
                    };

                    const px_c = if (worker_incx == 1)
                        worker_x + c_start
                    else blk: {
                        @import("../level1/copy.zig").k_copy(
                            c_len,
                            worker_x + numeric.cast(usize, if (worker_incx > 0)
                                numeric.cast(isize, c_start) * worker_incx
                            else
                                (-numeric.cast(isize, worker_n) + numeric.cast(isize, c_start + c_len)) * worker_incx),
                            worker_incx,
                            @as([*]X, &local_x_c),
                            1,
                        );

                        break :blk @as([*]const X, &local_x_c);
                    };

                    var j: usize = 0;
                    while (j < c_len) : (j += 1) {
                        // temp = worker_alpha * px_c[j]
                        const temp = numeric.mul(
                            worker_alpha,
                            px_c[j],
                        );

                        var i: usize = 0;
                        while (i < (r_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda] += temp * px_r[i + u]
                                numeric.fmaInto(
                                    &worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda],
                                    temp,
                                    px_r[i + u],
                                    worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda],
                                );
                            }
                        }

                        while (i < r_len) : (i += 1) {
                            // worker_a[r_start + c_start * worker_lda + i + j * worker_lda] += temp * px_r[i]
                            numeric.fmaInto(
                                &worker_a[r_start + c_start * worker_lda + i + j * worker_lda],
                                temp,
                                px_r[i],
                                worker_a[r_start + c_start * worker_lda + i + j * worker_lda],
                            );
                        }
                    }
                }
            }
        }
    };

    const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(X, numeric.Mul(Al, X), A)) orelse 2);
    comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (2 * @sizeOf(X) + @sizeOf(A)));
    tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

    const k = (n + tile_size - 1) / tile_size;
    const num_tiles = k * (k + 1) / 2;

    var atomic_counter = std.atomic.Value(usize).init(0);
    var threads: [options.max_threads]std.Thread = undefined;

    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        if (std.Thread.spawn(.{}, Worker.execute, .{
            eff_uplo,
            n,
            alpha,
            a,
            lda,
            x,
            incx,
            &atomic_counter,
            tile_size,
            num_tiles,
        })) |th| {
            threads[i] = th;
            spawned_count += 1;
        } else |err| {
            spawn_err = err;
            break;
        }
    }

    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
    }

    if (spawn_err) |err|
        return err;
}

fn k_syr(uplo: matrix.Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, a: anytype, lda: usize) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (n == 0 or numeric.eq(alpha, 0))
        return;

    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    if (uplo == .upper) {
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(X, numeric.Mul(Al, X), A)) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < n) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, n - tile_i);
            var local_x: [tile_size]X = undefined;

            const px = if (incx == 1)
                x + tile_i
            else blk: {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    x + numeric.cast(usize, if (incx > 0)
                        numeric.cast(isize, tile_i) * incx
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incx),
                    incx,
                    @as([*]X, &local_x),
                    1,
                );

                break :blk @as([*]const X, &local_x);
            };

            var j: usize = tile_i;
            var jx: isize = kx + numeric.cast(isize, tile_i) * incx;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                    const temp = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    if (j < tile_i + b_len) {
                        var i: usize = 0;
                        while (i < ((j - tile_i) / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += px[i + u] * temp
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    px[i + u],
                                    temp,
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < j - tile_i) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                px[i],
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }

                        // a[j + j * lda] += px[j - tile_i] * temp
                        numeric.fmaInto(
                            &a[j + j * lda],
                            px[j - tile_i],
                            temp,
                            a[j + j * lda],
                        );
                    } else {
                        var i: usize = 0;
                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += px[i + u] * temp
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    px[i + u],
                                    temp,
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                px[i],
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }
                    }
                }

                jx += incx;
            }
        }
    } else {
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(X, numeric.Mul(Al, X), A)) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < n) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, n - tile_i);
            var local_x: [tile_size]X = undefined;

            const px = if (incx == 1)
                x + tile_i
            else blk: {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    x + numeric.cast(usize, if (incx > 0)
                        numeric.cast(isize, tile_i) * incx
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incx),
                    incx,
                    @as([*]X, &local_x),
                    1,
                );

                break :blk @as([*]const X, &local_x);
            };

            var j: usize = 0;
            var jx: isize = kx;
            while (j < tile_i + b_len) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                    const temp = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    if (j >= tile_i) {
                        // a[j + j * lda] += px[j - tile_i] * temp
                        numeric.fmaInto(
                            &a[j + j * lda],
                            px[j - tile_i],
                            temp,
                            a[j + j * lda],
                        );

                        var i: usize = j - tile_i + 1;
                        while (i < b_len and i % unroll != 0) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                px[i],
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }

                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += px[i + u] * temp
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    px[i + u],
                                    temp,
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                px[i],
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }
                    } else {
                        var i: usize = 0;
                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += px[i + u] * temp
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    px[i + u],
                                    temp,
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                px[i],
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }
                    }
                }

                jx += incx;
            }
        }
    }
}
