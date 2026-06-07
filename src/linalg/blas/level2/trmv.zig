const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;
const Diag = meta.Diag;

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");
const float = @import("../../../float.zig");

const linalg = @import("../../../linalg.zig");

/// Computes a matrix-vector product using a triangular matrix defined as:
///
/// ```zig
/// x = A * x,
/// ```
///
/// or
///
/// ```zig
/// x = Aᵀ * x,
/// ```
///
/// or
///
/// ```zig
/// x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
/// x = Aᴴ * x,
/// ```
///
/// where `x` is an `n`-element vector, and `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.trmv(layout: Layout, uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, a: [*]const A, lda: usize, x: [*]X, incx: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
///   triangular matrix:
/// * `transa` (`linalg.Transpose`): Specifies the operation to be performed on
///   `A`:
///   * `no_transpose`: `x = A * x`
///   * `transpose`: `x = Aᵀ * x`
///   * `conj_no_transpose`: `x = conj(A) * x`
///   * `conj_transpose`: `x = Aᴴ * x`
/// * `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`. On return, contains the result of the
///   operation.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `opts`: Optional parameters:
///   * `num_threads` (`usize = 0`): Number of threads to spawn:
///     * `0`: automatic. The thread count is derived from `m * n` and
///       `parallel_threshold`:
///       ```zig
///       threads = max(1, min(std.Thread.getCpuCount(), options.max_threads, ((n * n + n) / 2) / parallel_threshold))
///       ```
///     * 1: force serial execution. parallel_threshold is ignored.
///     * N >= 2: use exactly N threads, clamped by
///       std.Thread.getCpuCount() and options.max_threads as a hard safety
///       ceiling. parallel_threshold is ignored.
///   * parallel_threshold (usize = 33_554_432 / @sizeOf(meta.Child(X))):
///     Minimum number of matrix elements (`m * n`) required to trigger
///     multithreaded execution.
///   * work (`?[*]X = null`): Mutable many-item pointer, size at least `n`.
///     Required for multithreading.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, n)`, or
///   if `incx` is 0.
pub fn trmv(
    layout: Layout,
    uplo: Uplo,
    transa: linalg.Transpose,
    diag: Diag,
    n: usize,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 33_554_432 / @sizeOf(meta.Child(@TypeOf(x))),
        work: ?[*]meta.Child(@TypeOf(x)) = null,
    },
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.trmv: a must be a many-item pointer to numerics, and x must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);

    if (lda < int.max(1, n) or incx == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();
    const eff_transa = if (layout == .col_major) transa else transa.invert();
    const noconj = transa == .no_trans or transa == .trans;

    // Quick return if possible.
    if (n == 0)
        return;

    if (opts.num_threads == 1 or opts.work == null)
        return if (noconj)
            k_trmv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, true)
        else
            k_trmv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, false);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, ((n * n + n) / 2) / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);
    num_threads = int.min(num_threads, n);

    if (num_threads <= 1)
        return if (noconj)
            k_trmv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, true)
        else
            k_trmv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, false);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);
    num_threads = int.min(num_threads, n);

    if (num_threads <= 1)
        return if (noconj)
            k_trmv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, true)
        else
            k_trmv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, false);

    @import("../level1/copy.zig").k_copy(n, x, incx, opts.work.?, 1);

    const Worker = struct {
        fn execute(
            worker_uplo: Uplo,
            worker_transa: linalg.Transpose,
            worker_diag: Diag,
            worker_n: usize,
            worker_a: [*]const A,
            worker_lda: usize,
            worker_x: [*]const X,
            worker_y: [*]X,
            worker_incy: isize,
            chunk_start: usize,
            chunk_end: usize,
            worker_noconj: bool,
        ) void {
            return if (worker_noconj)
                k_trmv_nd(worker_uplo, worker_transa, worker_diag, worker_n, worker_a, worker_lda, worker_x, worker_y, worker_incy, chunk_start, chunk_end, true)
            else
                k_trmv_nd(worker_uplo, worker_transa, worker_diag, worker_n, worker_a, worker_lda, worker_x, worker_y, worker_incy, chunk_start, chunk_end, false);
        }
    };

    var threads: [options.max_threads]std.Thread = undefined;
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const start_w = numeric.cast(f64, i) / numeric.cast(f64, num_threads);
        const end_w = numeric.cast(f64, i + 1) / numeric.cast(f64, num_threads);

        const chunk_start = numeric.cast(usize, if (eff_uplo == .upper)
            numeric.cast(f64, n) * (1.0 - float.sqrt(1.0 - start_w))
        else
            numeric.cast(f64, n) * float.sqrt(start_w));
        const chunk_end = int.min(n, numeric.cast(usize, if (eff_uplo == .upper)
            numeric.cast(f64, n) * (1.0 - float.sqrt(1.0 - end_w))
        else
            numeric.cast(f64, n) * float.sqrt(end_w)));

        if (std.Thread.spawn(.{}, Worker.execute, .{
            eff_uplo,
            eff_transa,
            diag,
            n,
            a,
            lda,
            opts.work.?,
            x,
            incx,
            chunk_start,
            chunk_end,
            noconj,
        })) |th| {
            threads[i] = th;
            spawned_count += 1;
        } else |err| {
            spawn_err = err;
            break;
        }
    }

    var t: usize = 0;
    while (t < spawned_count) : (t += 1)
        threads[t].join();

    if (spawn_err) |err|
        return err;
}

pub fn k_trmv(
    uplo: Uplo,
    transa: linalg.Transpose,
    diag: Diag,
    n: usize,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
    comptime noconj: bool,
) void {
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (n == 0)
        return;

    // Set up the start points in x.
    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Form  y = alpha * A * x + y  or  y = alpha * conj(A) * x + y.
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(X, A, X)) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        if (uplo == .upper) {
            var tile_j: usize = 0;
            while (tile_j < n) : (tile_j += tile_size) {
                const c_end = int.min(tile_j + tile_size, n);

                if (tile_j > 0) {
                    var tile_i: usize = 0;
                    while (tile_i < tile_j) : (tile_i += tile_size) {
                        const r_end = int.min(tile_i + tile_size, tile_j);
                        const r_len = r_end - tile_i;
                        var local_x: [tile_size]X = undefined;

                        const px = if (incx == 1)
                            x + tile_i
                        else blk: {
                            @import("../level1/copy.zig").k_copy(
                                r_len,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_i) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + r_len)) * incx),
                                incx,
                                @as([*]X, &local_x),
                                1,
                            );

                            break :blk @as([*]X, &local_x);
                        };

                        var j: usize = tile_j;
                        while (j < c_end) : (j += 1) {
                            const temp = x[numeric.cast(usize, kx + numeric.cast(isize, j) * incx)];

                            var i: usize = 0;
                            while (i < (r_len / unroll) * unroll) : (i += unroll) {
                                inline for (0..unroll) |u| {
                                    // px[i + u] += temp * a[tile_i + i + u + j * lda]
                                    numeric.fmaInto(
                                        &px[i + u],
                                        temp,
                                        if (comptime noconj)
                                            a[tile_i + i + u + j * lda]
                                        else
                                            numeric.conj(a[tile_i + i + u + j * lda]),
                                        px[i + u],
                                    );
                                }
                            }

                            while (i < r_len) : (i += 1) {
                                // px[i] += temp * a[tile_i + i + j * lda]
                                numeric.fmaInto(
                                    &px[i],
                                    temp,
                                    if (comptime noconj)
                                        a[tile_i + i + j * lda]
                                    else
                                        numeric.conj(a[tile_i + i + j * lda]),
                                    px[i],
                                );
                            }
                        }

                        if (incx != 1) {
                            @import("../level1/copy.zig").k_copy(
                                r_len,
                                @as([*]const X, &local_x),
                                1,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_i) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + r_len)) * incx),
                                incx,
                            );
                        }
                    }
                }

                var j: usize = tile_j;
                while (j < c_end) : (j += 1) {
                    const jx = kx + numeric.cast(isize, j) * incx;

                    const temp = x[numeric.cast(usize, jx)];

                    if (incx == 1) {
                        var i: usize = tile_j;
                        while (i < (j / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // x[i + u] += a[i + u + j * lda]
                                numeric.fmaInto(
                                    &x[i + u],
                                    temp,
                                    if (comptime noconj)
                                        a[i + u + j * lda]
                                    else
                                        numeric.conj(a[i + u + j * lda]),
                                    x[i + u],
                                );
                            }
                        }

                        while (i < j) : (i += 1) {
                            // x[i] += a[i + j * lda]
                            numeric.fmaInto(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                            );
                        }
                    } else {
                        var ix = kx + numeric.cast(isize, tile_j) * incx;
                        var i: usize = tile_j;
                        while (i < j) : (i += 1) {
                            // x[ix] += temp * a[i + j * lda]
                            numeric.fmaInto(
                                &x[numeric.cast(usize, ix)],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix += incx;
                        }
                    }

                    if (diag == .non_unit) {
                        // x[jx] *= a[j + j * lda]
                        numeric.mulInto(
                            &x[numeric.cast(usize, jx)],
                            x[numeric.cast(usize, jx)],
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                        );
                    }
                }
            }
        } else {
            var ct: usize = (n + tile_size - 1) / tile_size;
            while (ct > 0) {
                ct -= 1;

                const tile_j = ct * tile_size;
                const c_end = int.min(tile_j + tile_size, n);

                if (c_end < n) {
                    var tile_i: usize = c_end;
                    while (tile_i < n) : (tile_i += tile_size) {
                        const r_end = int.min(tile_i + tile_size, n);
                        const r_len = r_end - tile_i;
                        var local_x: [tile_size]X = undefined;

                        const px = if (incx == 1)
                            x + tile_i
                        else blk: {
                            @import("../level1/copy.zig").k_copy(
                                r_len,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_i) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + r_len)) * incx),
                                incx,
                                @as([*]X, &local_x),
                                1,
                            );

                            break :blk @as([*]X, &local_x);
                        };

                        var j: usize = tile_j;
                        while (j < c_end) : (j += 1) {
                            const temp = x[numeric.cast(usize, kx + numeric.cast(isize, j) * incx)];

                            var i: usize = 0;
                            while (i < (r_len / unroll) * unroll) : (i += unroll) {
                                inline for (0..unroll) |u| {
                                    // px[i + u] += temp * a[tile_i + i + u + j * lda]
                                    numeric.fmaInto(
                                        &px[i + u],
                                        temp,
                                        if (comptime noconj)
                                            a[tile_i + i + u + j * lda]
                                        else
                                            numeric.conj(a[tile_i + i + u + j * lda]),
                                        px[i + u],
                                    );
                                }
                            }

                            while (i < r_len) : (i += 1) {
                                // px[i] += temp * a[tile_i + i + j * lda]
                                numeric.fmaInto(
                                    &px[i],
                                    temp,
                                    if (comptime noconj)
                                        a[tile_i + i + j * lda]
                                    else
                                        numeric.conj(a[tile_i + i + j * lda]),
                                    px[i],
                                );
                            }
                        }

                        if (incx != 1) {
                            @import("../level1/copy.zig").k_copy(
                                r_len,
                                @as([*]const X, &local_x),
                                1,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_i) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + r_len)) * incx),
                                incx,
                            );
                        }
                    }
                }

                var j: usize = c_end;
                while (j > tile_j) {
                    j -= 1;

                    const jx = kx + numeric.cast(isize, j) * incx;

                    const temp = x[numeric.cast(usize, jx)];

                    if (incx == 1) {
                        var i: usize = j + 1;
                        while (i < c_end and i % unroll != 0) : (i += 1) {
                            // x[i] += temp * a[i + j * lda]
                            numeric.fmaInto(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                            );
                        }

                        while (i < (c_end / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // x[i + u] += a[i + u + j * lda]
                                numeric.fmaInto(
                                    &x[i + u],
                                    temp,
                                    if (comptime noconj)
                                        a[i + u + j * lda]
                                    else
                                        numeric.conj(a[i + u + j * lda]),
                                    x[i + u],
                                );
                            }
                        }

                        while (i < c_end) : (i += 1) {
                            // x[i] += a[i + j * lda]
                            numeric.fmaInto(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                            );
                        }
                    } else {
                        var ix = kx + numeric.cast(isize, j + 1) * incx;
                        var i: usize = j + 1;
                        while (i < c_end) : (i += 1) {
                            // x[ix] += a[i + j * lda]
                            numeric.fmaInto(
                                &x[numeric.cast(usize, ix)],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix += incx;
                        }
                    }

                    if (diag == .non_unit) {
                        numeric.mulInto(
                            &x[numeric.cast(usize, jx)],
                            x[numeric.cast(usize, jx)],
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                        );
                    }
                }
            }
        }
    } else {
        // Form  x = Aᵀ * x  or  x = Aᴴ * x.
        const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(A, X))) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(meta.Accumulator(numeric.Mul(A, X))) + @sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        if (uplo == .upper) {
            var rt: usize = (n + tile_size - 1) / tile_size;
            while (rt > 0) {
                rt -= 1;

                const tile_i = rt * tile_size;
                const r_end = int.min(tile_i + tile_size, n);

                var local_sums: [tile_size]meta.Accumulator(numeric.Mul(A, X)) = undefined;
                @memset(local_sums[0 .. r_end - tile_i], numeric.zero(meta.Accumulator(numeric.Mul(A, X))));

                var i_tile: usize = 0;
                while (i_tile < tile_i) : (i_tile += tile_size) {
                    const i_end = int.min(i_tile + tile_size, tile_i);
                    const i_len = i_end - i_tile;
                    var local_x: [tile_size]X = undefined;

                    const px = if (incx == 1)
                        x + i_tile
                    else blk: {
                        @import("../level1/copy.zig").k_copy(
                            i_len,
                            x + numeric.cast(usize, if (incx > 0)
                                numeric.cast(isize, i_tile) * incx
                            else
                                (-numeric.cast(isize, n) + numeric.cast(isize, i_tile + i_len)) * incx),
                            incx,
                            @as([*]X, &local_x),
                            1,
                        );

                        break :blk @as([*]const X, &local_x);
                    };

                    var j: usize = tile_i;
                    while (j < r_end) : (j += 1) {
                        var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                        var i: usize = 0;
                        while (i < (i_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // sums[u] += a[i_tile + i + u + j * lda] * px[i + u]
                                numeric.fmaInto(
                                    &sums[u],
                                    if (comptime noconj)
                                        a[i_tile + i + u + j * lda]
                                    else
                                        numeric.conj(a[i_tile + i + u + j * lda]),
                                    px[i + u],
                                    sums[u],
                                );
                            }
                        }

                        var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                        inline for (0..unroll) |u| {
                            // temp += sums[u]
                            numeric.addInto(
                                &temp,
                                temp,
                                sums[u],
                            );
                        }

                        while (i < i_len) : (i += 1) {
                            // temp += a[i_tile + i + j * lda] * px[i]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[i_tile + i + j * lda]
                                else
                                    numeric.conj(a[i_tile + i + j * lda]),
                                px[i],
                                temp,
                            );
                        }

                        // local_sums[j - tile_i] += temp
                        numeric.addInto(
                            &local_sums[j - tile_i],
                            local_sums[j - tile_i],
                            temp,
                        );
                    }
                }

                var j: usize = r_end;
                while (j > tile_i) {
                    j -= 1;

                    const jx = kx + numeric.cast(isize, j) * incx;

                    var temp = local_sums[j - tile_i];

                    if (diag == .non_unit) {
                        // temp += a[j + j * lda] * x[jx]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[numeric.cast(usize, jx)],
                            temp,
                        );
                    } else {
                        // temp += x[jx]
                        numeric.addInto(
                            &temp,
                            temp,
                            x[numeric.cast(usize, jx)],
                        );
                    }

                    var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                    if (incx == 1) {
                        var i: usize = tile_i;
                        while (i < (j / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // sums[u] += a[i + u + j * lda] * x[i + u]
                                numeric.fmaInto(
                                    &sums[u],
                                    if (comptime noconj)
                                        a[i + u + j * lda]
                                    else
                                        numeric.conj(a[i + u + j * lda]),
                                    x[i + u],
                                    sums[u],
                                );
                            }
                        }

                        inline for (0..unroll) |u| {
                            // temp += sums[u]
                            numeric.addInto(
                                &temp,
                                temp,
                                sums[u],
                            );
                        }

                        while (i < j) : (i += 1) {
                            // temp += a[i + j * lda] * x[i]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                                temp,
                            );
                        }
                    } else {
                        var i: usize = tile_i;
                        while (i < (j / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                const ix = kx + numeric.cast(isize, i + u) * incx;

                                // sums[u] += a[i + u + j * lda] * x[ix]
                                numeric.fmaInto(
                                    &sums[u],
                                    if (comptime noconj)
                                        a[i + u + j * lda]
                                    else
                                        numeric.conj(a[i + u + j * lda]),
                                    x[numeric.cast(usize, ix)],
                                    sums[u],
                                );
                            }
                        }

                        inline for (0..unroll) |u| {
                            // temp += sums[u]
                            numeric.addInto(&temp, temp, sums[u]);
                        }

                        while (i < j) : (i += 1) {
                            const ix = kx + numeric.cast(isize, i) * incx;

                            // temp += a[i + j * lda] * x[ix]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                                temp,
                            );
                        }
                    }

                    numeric.set(&x[numeric.cast(usize, jx)], temp);
                }
            }
        } else {
            var tile_i: usize = 0;
            while (tile_i < n) : (tile_i += tile_size) {
                const r_end = int.min(tile_i + tile_size, n);

                var local_sums: [tile_size]meta.Accumulator(numeric.Mul(A, X)) = undefined;
                @memset(local_sums[0 .. r_end - tile_i], numeric.zero(meta.Accumulator(numeric.Mul(A, X))));

                var i_tile: usize = r_end;
                while (i_tile < n) : (i_tile += tile_size) {
                    const i_end = int.min(i_tile + tile_size, n);
                    const i_len = i_end - i_tile;
                    var local_x: [tile_size]X = undefined;

                    const px = if (incx == 1)
                        x + i_tile
                    else blk: {
                        @import("../level1/copy.zig").k_copy(
                            i_len,
                            x + numeric.cast(usize, if (incx > 0)
                                numeric.cast(isize, i_tile) * incx
                            else
                                (-numeric.cast(isize, n) + numeric.cast(isize, i_tile + i_len)) * incx),
                            incx,
                            @as([*]X, &local_x),
                            1,
                        );

                        break :blk @as([*]const X, &local_x);
                    };

                    var j: usize = tile_i;
                    while (j < r_end) : (j += 1) {
                        var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                        var i: usize = 0;
                        while (i < (i_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // sums[u] = a[i_tile + i + u + j * lda] * px[i + u]
                                numeric.fmaInto(
                                    &sums[u],
                                    if (comptime noconj)
                                        a[i_tile + i + u + j * lda]
                                    else
                                        numeric.conj(a[i_tile + i + u + j * lda]),
                                    px[i + u],
                                    sums[u],
                                );
                            }
                        }

                        var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                        inline for (0..unroll) |u| {
                            numeric.addInto(&temp, temp, sums[u]);
                        }

                        while (i < i_len) : (i += 1) {
                            // temp = a[i_tile + i + j * lda] * px[i]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[i_tile + i + j * lda]
                                else
                                    numeric.conj(a[i_tile + i + j * lda]),
                                px[i],
                                temp,
                            );
                        }

                        // local_sums[j - tile_i] += temp
                        numeric.addInto(
                            &local_sums[j - tile_i],
                            local_sums[j - tile_i],
                            temp,
                        );
                    }
                }

                var j: usize = tile_i;
                while (j < r_end) : (j += 1) {
                    const jx = kx + numeric.cast(isize, j) * incx;

                    var temp = local_sums[j - tile_i];

                    if (diag == .non_unit) {
                        // temp += a[j + j * lda] * x[jx]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[numeric.cast(usize, jx)],
                            temp,
                        );
                    } else {
                        // temp += x[jx]
                        numeric.addInto(
                            &temp,
                            temp,
                            x[numeric.cast(usize, jx)],
                        );
                    }

                    var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                    if (incx == 1) {
                        var i: usize = j + 1;
                        while (i < r_end and i % unroll != 0) : (i += 1) {
                            // temp += a[i + j * lda] * x[i]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                                temp,
                            );
                        }

                        while (i < (r_end / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // sums[u] += a[i + u + j * lda] * x[i + u]
                                numeric.fmaInto(
                                    &sums[u],
                                    if (comptime noconj)
                                        a[i + u + j * lda]
                                    else
                                        numeric.conj(a[i + u + j * lda]),
                                    x[i + u],
                                    sums[u],
                                );
                            }
                        }

                        inline for (0..unroll) |u| {
                            // temp += sums[u]
                            numeric.addInto(
                                &temp,
                                temp,
                                sums[u],
                            );
                        }

                        while (i < r_end) : (i += 1) {
                            // temp += a[i + j * lda] * x[i]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                                temp,
                            );
                        }
                    } else {
                        var i: usize = j + 1;
                        while (i < r_end and i % unroll != 0) : (i += 1) {
                            const ix = kx + numeric.cast(isize, i) * incx;

                            // temp += a[i + j * lda] * x[ix]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                                temp,
                            );
                        }

                        while (i < (r_end / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                const ix = kx + numeric.cast(isize, i + u) * incx;

                                // temp += a[i + u + j * lda] * x[ix]
                                numeric.fmaInto(
                                    &sums[u],
                                    if (comptime noconj)
                                        a[i + u + j * lda]
                                    else
                                        numeric.conj(a[i + u + j * lda]),
                                    x[numeric.cast(usize, ix)],
                                    sums[u],
                                );
                            }
                        }

                        inline for (0..unroll) |u| {
                            // temp += sums[u]
                            numeric.addInto(
                                &temp,
                                temp,
                                sums[u],
                            );
                        }

                        while (i < r_end) : (i += 1) {
                            const ix = kx + numeric.cast(isize, i) * incx;

                            // temp += a[i + j * lda] * x[ix]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                                temp,
                            );
                        }
                    }

                    numeric.set(&x[numeric.cast(usize, jx)], temp);
                }
            }
        }
    }

    return;
}

fn k_trmv_nd(
    uplo: Uplo,
    transa: linalg.Transpose,
    diag: Diag,
    n: usize,
    a: anytype,
    lda: usize,
    x: anytype,
    y: anytype,
    incy: isize,
    chunk_start: usize,
    chunk_end: usize,
    comptime noconj: bool,
) void {
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Quick return if possible.
    if (n == 0 or chunk_start >= chunk_end)
        return;

    // Set up the start points in y.
    const ky: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Form  y = A * x  or  y = conj(A) * x.
        const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(A, Y))) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(meta.Accumulator(numeric.Mul(A, X))) + @sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        if (uplo == .upper) {
            var tile_i: usize = chunk_start;
            while (tile_i < chunk_end) : (tile_i += tile_size) {
                const ti_end = int.min(tile_i + tile_size, chunk_end);
                const ti_len = ti_end - tile_i;

                var local_sums: [tile_size]meta.Accumulator(numeric.Mul(A, X)) = undefined;
                @memset(local_sums[0..ti_len], numeric.zero(meta.Accumulator(numeric.Mul(A, X))));

                var j: usize = tile_i;
                while (j < ti_end) : (j += 1) {
                    if (diag == .non_unit) {
                        // local_sums[j - tile_i] += a[j + j * lda] * x[j]
                        numeric.fmaInto(
                            &local_sums[j - tile_i],
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[j],
                            local_sums[j - tile_i],
                        );
                    } else {
                        // local_sums[j - tile_i] += x[j]
                        numeric.addInto(
                            &local_sums[j - tile_i],
                            local_sums[j - tile_i],
                            x[j],
                        );
                    }

                    var i: usize = tile_i;
                    // unroll?
                    while (i < j) : (i += 1) {
                        // local_sums[i - tile_i]
                        numeric.fmaInto(
                            &local_sums[i - tile_i],
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[j],
                            local_sums[i - tile_i],
                        );
                    }
                }

                var tile_j: usize = ti_end;
                while (tile_j < n) : (tile_j += tile_size) {
                    const tj_end = int.min(tile_j + tile_size, n);
                    const tj_len = tj_end - tile_j;

                    const px = x + tile_j;

                    j = 0;
                    while (j < tj_len) : (j += 1) {
                        var i: usize = 0;
                        while (i < (ti_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // local_sums[i + u] += a[tile_i + i + u + (tile_j + j) * lda] + px[j]
                                numeric.fmaInto(
                                    &local_sums[i + u],
                                    if (comptime noconj)
                                        a[tile_i + i + u + (tile_j + j) * lda]
                                    else
                                        numeric.conj(a[tile_i + i + u + (tile_j + j) * lda]),
                                    px[j],
                                    local_sums[i + u],
                                );
                            }
                        }

                        while (i < ti_len) : (i += 1) {
                            // local_sums[i] += a[tile_i + i + (tile_j + j) * lda] + px[j]
                            numeric.fmaInto(
                                &local_sums[i],
                                if (comptime noconj)
                                    a[tile_i + i + (tile_j + j) * lda]
                                else
                                    numeric.conj(a[tile_i + i + (tile_j + j) * lda]),
                                px[j],
                                local_sums[i],
                            );
                        }
                    }
                }

                var i: usize = 0;
                while (i < ti_len) : (i += 1) {
                    numeric.set(&y[numeric.cast(usize, ky + numeric.cast(isize, tile_i + i) * incy)], local_sums[i]);
                }
            }
        } else {
            var tile_i: usize = chunk_start;
            while (tile_i < chunk_end) : (tile_i += tile_size) {
                const ti_end = int.min(tile_i + tile_size, chunk_end);
                const ti_len = ti_end - tile_i;

                var local_sums: [tile_size]meta.Accumulator(numeric.Mul(A, X)) = undefined;
                @memset(local_sums[0..ti_len], numeric.zero(meta.Accumulator(numeric.Mul(A, X))));

                var tile_j: usize = 0;
                while (tile_j < tile_i) : (tile_j += tile_size) {
                    const tj_end = int.min(tile_j + tile_size, tile_i);
                    const tj_len = tj_end - tile_j;

                    const px = x + tile_j;

                    var j: usize = 0;
                    while (j < tj_len) : (j += 1) {
                        var i: usize = 0;
                        while (i < (ti_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // local_sums[i + u] += a[tile_i + i + u + (tile_j + j) * lda] * px[j]
                                numeric.fmaInto(
                                    &local_sums[i + u],
                                    if (comptime noconj)
                                        a[tile_i + i + u + (tile_j + j) * lda]
                                    else
                                        numeric.conj(a[tile_i + i + u + (tile_j + j) * lda]),
                                    px[j],
                                    local_sums[i + u],
                                );
                            }
                        }

                        while (i < ti_len) : (i += 1) {
                            // local_sums[i] += a[tile_i + i + (tile_j + j) * lda] * px[j]
                            numeric.fmaInto(
                                &local_sums[i],
                                if (comptime noconj)
                                    a[tile_i + i + (tile_j + j) * lda]
                                else
                                    numeric.conj(a[tile_i + i + (tile_j + j) * lda]),
                                px[j],
                                local_sums[i],
                            );
                        }
                    }
                }

                var j: usize = tile_i;
                while (j < ti_end) : (j += 1) {
                    if (diag == .non_unit) {
                        // local_sums[j - tile_i] += a[j + j * lda] * x[j]
                        numeric.fmaInto(
                            &local_sums[j - tile_i],
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[j],
                            local_sums[j - tile_i],
                        );
                    } else {
                        // local_sums[j - tile_i] += x[j]
                        numeric.addInto(
                            &local_sums[j - tile_i],
                            local_sums[j - tile_i],
                            x[j],
                        );
                    }

                    var i: usize = j + 1;
                    // unroll?
                    while (i < ti_end) : (i += 1) {
                        // local_sums[i - tile_i] += a[i + j * lda] * x[j]
                        numeric.fmaInto(
                            &local_sums[i - tile_i],
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[j],
                            local_sums[i - tile_i],
                        );
                    }
                }

                var i: usize = 0;
                while (i < ti_len) : (i += 1) {
                    numeric.set(&y[numeric.cast(usize, ky + numeric.cast(isize, tile_i + i) * incy)], local_sums[i]);
                }
            }
        }
    } else {
        // Form  y = Aᵀ * x  or  y = Aᴴ * x.
        const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(A, Y))) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(meta.Accumulator(numeric.Mul(A, X))) + @sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        if (uplo == .upper) {
            var tile_j: usize = chunk_start;
            while (tile_j < chunk_end) : (tile_j += tile_size) {
                const tj_end = int.min(tile_j + tile_size, chunk_end);

                var local_sums: [tile_size]meta.Accumulator(numeric.Mul(A, X)) = undefined;
                @memset(local_sums[0 .. tj_end - tile_j], numeric.zero(meta.Accumulator(numeric.Mul(A, X))));

                var j: usize = tile_j;
                while (j < tj_end) : (j += 1) {
                    if (diag == .non_unit) {
                        // local_sums[j - tile_j] += a[j + j * lda]
                        numeric.fmaInto(
                            &local_sums[j - tile_j],
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[j],
                            local_sums[j - tile_j],
                        );
                    } else {
                        // local_sums[j - tile_j] += x[j]
                        numeric.addInto(
                            &local_sums[j - tile_j],
                            local_sums[j - tile_j],
                            x[j],
                        );
                    }
                }

                var tile_i: usize = 0;
                while (tile_i < tile_j) : (tile_i += tile_size) {
                    const ti_end = int.min(tile_i + tile_size, tile_j);
                    const ti_len = ti_end - tile_i;

                    const px = x + tile_i;

                    j = tile_j;
                    while (j < tj_end) : (j += 1) {
                        var sums: [unroll]meta.Accumulator(numeric.Mul(A, Y)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, Y)))} ** unroll;

                        var i: usize = 0;
                        while (i < (ti_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // sums[u] += a[tile_i + i + u + j * lda] * px[i + u]
                                numeric.fmaInto(
                                    &sums[u],
                                    if (comptime noconj)
                                        a[tile_i + i + u + j * lda]
                                    else
                                        numeric.conj(a[tile_i + i + u + j * lda]),
                                    px[i + u],
                                    sums[u],
                                );
                            }
                        }

                        var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, Y)));

                        inline for (0..unroll) |u| {
                            // temp += sums[u]
                            numeric.addInto(
                                &temp,
                                temp,
                                sums[u],
                            );
                        }

                        while (i < ti_len) : (i += 1) {
                            // temp += a[tile_i + i + j * lda] * px[i]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[tile_i + i + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + j * lda]),
                                px[i],
                                temp,
                            );
                        }

                        // local_sums[j - tile_j] += temp
                        numeric.addInto(
                            &local_sums[j - tile_j],
                            local_sums[j - tile_j],
                            temp,
                        );
                    }
                }

                j = tile_j + 1;
                while (j < tj_end) : (j += 1) {
                    var sums: [unroll]meta.Accumulator(numeric.Mul(A, Y)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, Y)))} ** unroll;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, Y)));

                    var i: usize = tile_j;
                    while (i < j and i % unroll != 0) : (i += 1) {
                        // temp += a[i + j * lda] * x[i]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[i],
                            temp,
                        );
                    }

                    while (i < (j / unroll) * unroll) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // sums[u] += a[i + u + j * lda] * x[i + u]
                            numeric.fmaInto(
                                &sums[u],
                                if (comptime noconj)
                                    a[i + u + j * lda]
                                else
                                    numeric.conj(a[i + u + j * lda]),
                                x[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        // temp += sums[u]
                        numeric.addInto(
                            &temp,
                            temp,
                            sums[u],
                        );
                    }

                    while (i < j) : (i += 1) {
                        // temp += a[i + j * lda] * x[i]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[i],
                            temp,
                        );
                    }

                    // local_sums[j - tile_j] += temp
                    numeric.addInto(
                        &local_sums[j - tile_j],
                        local_sums[j - tile_j],
                        temp,
                    );
                }

                j = tile_j;
                while (j < tj_end) : (j += 1) {
                    numeric.set(&y[numeric.cast(usize, ky + numeric.cast(isize, j) * incy)], local_sums[j - tile_j]);
                }
            }
        } else {
            var tile_j: usize = chunk_start;
            while (tile_j < chunk_end) : (tile_j += tile_size) {
                const tj_end = int.min(tile_j + tile_size, chunk_end);

                var local_sums: [tile_size]meta.Accumulator(numeric.Mul(A, X)) = undefined;
                @memset(local_sums[0 .. tj_end - tile_j], numeric.zero(meta.Accumulator(numeric.Mul(A, X))));

                var j: usize = tile_j;
                while (j < tj_end) : (j += 1) {
                    if (diag == .non_unit) {
                        // local_sums[j - tile_j] += a[j + j * lda]
                        numeric.fmaInto(
                            &local_sums[j - tile_j],
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[j],
                            local_sums[j - tile_j],
                        );
                    } else {
                        // local_sums[j - tile_j] += x[j]
                        numeric.addInto(
                            &local_sums[j - tile_j],
                            local_sums[j - tile_j],
                            x[j],
                        );
                    }
                }

                j = tile_j;
                while (j + 1 < tj_end) : (j += 1) {
                    var sums: [unroll]meta.Accumulator(numeric.Mul(A, Y)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, Y)))} ** unroll;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, Y)));

                    var i: usize = j + 1;
                    while (i < tj_end and i % unroll != 0) : (i += 1) {
                        // temp += a[i + j * lda] * x[i]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[i],
                            temp,
                        );
                    }

                    while (i < (tj_end / unroll) * unroll) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // sums[u] += a[i + u + j * lda] * x[i + u]
                            numeric.fmaInto(
                                &sums[u],
                                if (comptime noconj)
                                    a[i + u + j * lda]
                                else
                                    numeric.conj(a[i + u + j * lda]),
                                x[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        // temp += sums[u]
                        numeric.addInto(
                            &temp,
                            temp,
                            sums[u],
                        );
                    }

                    while (i < tj_end) : (i += 1) {
                        // temp += a[i + j * lda] * x[i]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[i],
                            temp,
                        );
                    }

                    // local_sums[j - tile_j] += temp
                    numeric.addInto(
                        &local_sums[j - tile_j],
                        local_sums[j - tile_j],
                        temp,
                    );
                }

                var tile_i: usize = tj_end;
                while (tile_i < n) : (tile_i += tile_size) {
                    const ti_end = int.min(tile_i + tile_size, n);
                    const ti_len = ti_end - tile_i;

                    const px = x + tile_i;

                    j = tile_j;
                    while (j < tj_end) : (j += 1) {
                        var sums: [unroll]meta.Accumulator(numeric.Mul(A, Y)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, Y)))} ** unroll;

                        var i: usize = 0;
                        while (i < (ti_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // sums[u] += a[tile_i + i + u + j * lda] * px[i + u]
                                numeric.fmaInto(
                                    &sums[u],
                                    if (comptime noconj)
                                        a[tile_i + i + u + j * lda]
                                    else
                                        numeric.conj(a[tile_i + i + u + j * lda]),
                                    px[i + u],
                                    sums[u],
                                );
                            }
                        }

                        var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, Y)));

                        inline for (0..unroll) |u| {
                            // temp += sums[u]
                            numeric.addInto(
                                &temp,
                                temp,
                                sums[u],
                            );
                        }

                        while (i < ti_len) : (i += 1) {
                            // temp += a[tile_i + i + j * lda] * px[i]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconj)
                                    a[tile_i + i + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + j * lda]),
                                px[i],
                                temp,
                            );
                        }

                        // local_sums[j - tile_j] += temp
                        numeric.addInto(
                            &local_sums[j - tile_j],
                            local_sums[j - tile_j],
                            temp,
                        );
                    }
                }

                j = tile_j;
                while (j < tj_end) : (j += 1) {
                    numeric.set(&y[numeric.cast(usize, ky + numeric.cast(isize, j) * incy)], local_sums[j - tile_j]);
                }
            }
        }
    }
}
