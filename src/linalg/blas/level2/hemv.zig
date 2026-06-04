const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");
const float = @import("../../../float.zig");

const linalg = @import("../../../linalg.zig");

/// Computes a matrix-vector product with a Hermitian matrix defined as:
///
/// ```zig
/// y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are numerics, `x` and `y` are `n`-element vectors,
/// and `A` is an `n`-by-`n` Hermitian matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.hemv(layout: Layout, uplo: Uplo, n: usize, alpha: Al, a: [*]const A, lda: usize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of
///   the Hermitian matrix `A` is used.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `y` need not be set on input.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`. On return, contains the result of the
///   operation.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
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
///   * parallel_threshold (usize = 2_097_152 / @sizeOf(meta.Child(Y))):
///     Minimum number of matrix elements (`n * n`) required to trigger
///     multithreaded execution.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, n)`, or
///   if `incx` or `incy` is 0.
pub fn hemv(
    layout: Layout,
    uplo: Uplo,
    n: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 2_097_152 / @sizeOf(meta.Child(@TypeOf(y))),
    },
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.hemv: alpha and beta must be numerics, ap and x must be many-item pointers to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);
    Y = meta.Child(Y);

    if (lda < int.max(1, n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();
    const noconj = layout == .col_major;

    // Quick return if possible.
    if (n == 0)
        return;

    if (numeric.eq(alpha, 0)) {
        if (numeric.ne(beta, 1))
            @import("../level1/scal.zig").k_scal(n, beta, y, incy);

        return;
    }

    if (opts.num_threads == 1)
        return if (noconj)
            k_hemv(eff_uplo, n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_hemv(eff_uplo, n, alpha, a, lda, x, incx, beta, y, incy, false);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, (n * n) / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);
    num_threads = int.min(num_threads, n);

    if (num_threads <= 1)
        return if (noconj)
            k_hemv(eff_uplo, n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_hemv(eff_uplo, n, alpha, a, lda, x, incx, beta, y, incy, false);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);
    num_threads = int.min(num_threads, n);

    if (num_threads <= 1)
        return if (noconj)
            k_hemv(eff_uplo, n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_hemv(eff_uplo, n, alpha, a, lda, x, incx, beta, y, incy, false);

    // Scale y before threading because threads will only accumulate onto y.
    if (numeric.ne(beta, 1))
        @import("../level1/scal.zig").k_scal(n, beta, y, incy);

    if (numeric.eq(alpha, 0))
        return;

    const Worker = struct {
        fn execute(
            worker_uplo: Uplo,
            worker_n: usize,
            worker_alpha: Al,
            worker_a: [*]const A,
            worker_lda: usize,
            worker_x: [*]const X,
            worker_incx: isize,
            worker_y: [*]Y,
            worker_incy: isize,
            comptime worker_noconj: bool,
            counter: *std.atomic.Value(usize),
            comptime worker_tile_size: comptime_int,
            worker_num_tiles: usize,
        ) void {
            const ky: isize = if (worker_incy < 0) (-numeric.cast(isize, worker_n) + 1) * worker_incy else 0;

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
                    var local_y: [worker_tile_size]Y = .{numeric.zero(Y)} ** worker_tile_size;

                    const px = if (worker_incx == 1)
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, r_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incx)
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

                    k_hemv(
                        worker_uplo,
                        r_len,
                        worker_alpha,
                        worker_a + r_start + c_start * worker_lda,
                        worker_lda,
                        px,
                        1,
                        numeric.one(Be),
                        @as([*]Y, &local_y),
                        1,
                        worker_noconj,
                    );

                    // Flush back y to global memory.
                    var i: usize = 0;
                    while (i < r_len) : (i += 1) {
                        // y += local_y[ky + (r_start + i) * incy]
                        numeric.atomicAdd_(
                            &worker_y[numeric.cast(usize, ky + numeric.cast(isize, r_start + i) * worker_incy)],
                            local_y[i],
                        );
                    }
                } else {
                    const unroll = 2 * int.min(
                        std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2,
                        std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(if (comptime worker_noconj) A else numeric.Conj(A), X))) orelse 2,
                    );

                    // Off-diagonal tile
                    var local_x_r: [worker_tile_size]X = undefined;
                    var local_y_r: [worker_tile_size]Y = .{numeric.zero(Y)} ** worker_tile_size;
                    var local_x_c: [worker_tile_size]X = undefined;
                    var local_y_c: [worker_tile_size]Y = .{numeric.zero(Y)} ** worker_tile_size;

                    const px_r = if (worker_incx == 1)
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, r_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incx)
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
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, c_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, c_start + c_len)) * worker_incx)
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
                        // temp1 = worker_alpha * local_x_c[j]
                        const temp1 = numeric.mul(worker_alpha, px_c[j]);

                        var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime worker_noconj) A else numeric.Conj(A), X)));

                        var sums: [unroll]meta.Accumulator(numeric.Mul(if (worker_noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime worker_noconj) A else numeric.Conj(A), X)))} ** unroll;

                        var i: usize = 0;
                        while (i < (r_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // local_y_r[i + u] += temp1 * worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda]
                                numeric.fma_(
                                    &local_y_r[i + u],
                                    temp1,
                                    if (comptime worker_noconj)
                                        worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda]
                                    else
                                        numeric.conj(worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda]),
                                    local_y_r[i + u],
                                );

                                // sums[u] += conj(worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda]) * local_x_r[i + u]
                                numeric.fma_(
                                    &sums[u],
                                    if (comptime worker_noconj)
                                        numeric.conj(worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda])
                                    else
                                        worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda],
                                    px_r[i + u],
                                    sums[u],
                                );
                            }
                        }

                        inline for (0..unroll) |u| {
                            numeric.add_(&temp2, temp2, sums[u]);
                        }

                        while (i < r_len) : (i += 1) {
                            // local_y_r[i] += temp1 * worker_a[r_start + c_start * worker_lda + i + j * worker_lda]
                            numeric.fma_(
                                &local_y_r[i],
                                temp1,
                                if (comptime worker_noconj)
                                    worker_a[r_start + c_start * worker_lda + i + j * worker_lda]
                                else
                                    numeric.conj(worker_a[r_start + c_start * worker_lda + i + j * worker_lda]),
                                local_y_r[i],
                            );

                            // temp2 += conj(worker_a[r_start + c_start * worker_lda + i + j * worker_lda]) * local_x_r[i]
                            numeric.fma_(
                                &temp2,
                                if (comptime worker_noconj)
                                    numeric.conj(worker_a[r_start + c_start * worker_lda + i + j * worker_lda])
                                else
                                    worker_a[r_start + c_start * worker_lda + i + j * worker_lda],
                                px_r[i],
                                temp2,
                            );
                        }

                        numeric.fma_(&local_y_c[j], worker_alpha, temp2, local_y_c[j]);
                    }

                    // Flush y back to global memory.
                    var i: usize = 0;
                    while (i < r_len) : (i += 1) {
                        // y += local_y_r[ky + (r_start + i) * incy]
                        numeric.atomicAdd_(
                            &worker_y[numeric.cast(usize, ky + numeric.cast(isize, r_start + i) * worker_incy)],
                            local_y_r[i],
                        );
                    }

                    j = 0;
                    while (j < c_len) : (j += 1) {
                        // y += local_y_c[ky + (c_start + j) * incy]
                        numeric.atomicAdd_(
                            &worker_y[numeric.cast(usize, ky + numeric.cast(isize, c_start + j) * worker_incy)],
                            local_y_c[j],
                        );
                    }
                }
            }
        }
    };

    const unroll = 2 * int.min(
        std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2,
        std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(A, X))) orelse 2,
    );
    comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (2 * @sizeOf(X) + 2 * @sizeOf(Y) + @sizeOf(A)));
    tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

    const k = (n + tile_size - 1) / tile_size;
    const num_tiles = k * (k + 1) / 2;

    var atomic_counter = std.atomic.Value(usize).init(0);
    var threads: [options.max_threads]std.Thread = undefined;

    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        if (if (noconj)
            std.Thread.spawn(.{}, Worker.execute, .{
                eff_uplo,
                n,
                alpha,
                a,
                lda,
                x,
                incx,
                y,
                incy,
                true,
                &atomic_counter,
                tile_size,
                num_tiles,
            })
        else
            std.Thread.spawn(.{}, Worker.execute, .{
                eff_uplo,
                n,
                alpha,
                a,
                lda,
                x,
                incx,
                y,
                incy,
                false,
                &atomic_counter,
                tile_size,
                num_tiles,
            })) |th|
        {
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

fn k_hemv(uplo: Uplo, n: usize, alpha: anytype, a: anytype, lda: usize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize, comptime noconj: bool) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Quick return if possible.
    if (n == 0 or (numeric.eq(alpha, 0) and numeric.eq(beta, 1)))
        return;

    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;

    // First form  y = beta * y.
    if (numeric.ne(beta, 1))
        @import("../level1/scal.zig").k_scal(n, beta, y, incy);

    if (numeric.eq(alpha, 0))
        return;

    if (uplo == .upper) {
        // Form  y  when upper triangle of A is stored.
        const unroll = 2 * int.min(
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2,
            std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X))) orelse 2,
        );
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(Y) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < n) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, n - tile_i);
            var local_x: [tile_size]X = undefined;
            var local_y: [tile_size]Y = undefined;

            const px = if (incx == 1)
                x + numeric.cast(usize, if (incx > 0)
                    numeric.cast(isize, tile_i) * incx
                else
                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incx)
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

                break :blk @as([*]X, &local_x);
            };

            const py = if (incy == 1)
                y + numeric.cast(usize, if (incy > 0)
                    numeric.cast(isize, tile_i) * incy
                else
                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy)
            else blk: {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                    @as([*]Y, &local_y),
                    1,
                );

                break :blk @as([*]Y, &local_y);
            };

            var j: usize = tile_i;
            var jx: isize = kx + numeric.cast(isize, tile_i) * incx;
            var jy: isize = ky + numeric.cast(isize, tile_i) * incy;
            while (j < n) : (j += 1) {
                if (j < tile_i + b_len) {
                    // temp1 = alpha * px[j - tile_i]
                    const temp1 = numeric.mul(alpha, px[j - tile_i]);

                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                    var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                    var i: usize = 0;
                    while (i < ((j - tile_i) / unroll) * unroll) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // py[i + u] += temp1 * a[tile_i + i + u + j * lda]
                            numeric.fma_(
                                &py[i + u],
                                temp1,
                                if (comptime noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                py[i + u],
                            );

                            // sums[u] += conj(a[tile_i + i + u + j * lda]) * px[i + u]
                            numeric.fma_(
                                &sums[u],
                                if (comptime !noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                px[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        numeric.add_(&temp2, temp2, sums[u]);
                    }

                    while (i < j - tile_i) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fma_(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) * px[i]
                        numeric.fma_(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    // py[j - tile_i] += temp1 * re(a[j + j * lda]) + alpha * temp2
                    numeric.fma_(
                        &py[j - tile_i],
                        temp1,
                        numeric.re(a[j + j * lda]),
                        numeric.fma(
                            alpha,
                            temp2,
                            py[j - tile_i],
                        ),
                    );
                } else {
                    // temp1 = alpha * x[jx]
                    const temp1 = numeric.mul(alpha, x[numeric.cast(usize, jx)]);

                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                    var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                    var i: usize = 0;
                    while (i + unroll <= b_len) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // py[i + u] += temp1 * a[tile_i + i + u + j * lda]
                            numeric.fma_(
                                &py[i + u],
                                temp1,
                                if (comptime noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                py[i + u],
                            );

                            // sums[u] += conj(a[tile_i + i + u + j * lda]) * px[i + u]
                            numeric.fma_(
                                &sums[u],
                                if (comptime !noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                px[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        numeric.add_(&temp2, temp2, sums[u]);
                    }

                    while (i < b_len) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fma_(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) + px[i]
                        numeric.fma_(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    // y[jy] += alpha * temp2
                    numeric.fma_(
                        &y[numeric.cast(usize, jy)],
                        alpha,
                        temp2,
                        y[numeric.cast(usize, jy)],
                    );
                }

                jx += incx;
                jy += incy;
            }

            if (incy != 1) {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    @as([*]Y, &local_y),
                    1,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                );
            }
        }
    } else {
        // Form  y  when lower triangle of A is stored.
        const unroll = 2 * int.min(
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2,
            std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X))) orelse 2,
        );
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(Y) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < n) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, n - tile_i);
            var local_x: [tile_size]X = undefined;
            var local_y: [tile_size]Y = undefined;

            const px = if (incx == 1)
                x + numeric.cast(usize, if (incx > 0)
                    numeric.cast(isize, tile_i) * incx
                else
                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incx)
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

                break :blk @as([*]X, &local_x);
            };

            const py = if (incy == 1)
                y + numeric.cast(usize, if (incy > 0)
                    numeric.cast(isize, tile_i) * incy
                else
                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy)
            else blk: {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                    @as([*]Y, &local_y),
                    1,
                );

                break :blk @as([*]Y, &local_y);
            };

            var j: usize = 0;
            var jx: isize = kx;
            var jy: isize = ky;
            while (j < tile_i + b_len) : (j += 1) {
                if (j >= tile_i) {
                    // temp1 = alpha * px[j - tile_i]
                    const temp1 = numeric.mul(
                        alpha,
                        px[j - tile_i],
                    );

                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                    // py[j - tile_i] += temp1 * re(a[j + j * lda])
                    numeric.fma_(
                        &py[j - tile_i],
                        temp1,
                        numeric.re(a[j + j * lda]),
                        py[j - tile_i],
                    );

                    var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                    var i: usize = j - tile_i + 1;
                    while (i < b_len and i % unroll != 0) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fma_(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) * px[i]
                        numeric.fma_(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    while (i < (b_len / unroll) * unroll) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // py[i + u] += temp1 * a[tile_i + i + u + j * lda]
                            numeric.fma_(
                                &py[i + u],
                                temp1,
                                if (comptime noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                py[i + u],
                            );

                            // sums[u] += conj(a[tile_i + i + u + j * lda]) * px[i + u]
                            numeric.fma_(
                                &sums[u],
                                if (comptime !noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                px[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        numeric.add_(&temp2, temp2, sums[u]);
                    }

                    while (i < b_len) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fma_(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) * px[i]
                        numeric.fma_(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    // py[j - tile_i] += alpha * temp2
                    numeric.fma_(
                        &py[j - tile_i],
                        alpha,
                        temp2,
                        py[j - tile_i],
                    );
                } else {
                    // temp1 = alpha * x[jx]
                    const temp1 = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                    var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                    var i: usize = 0;
                    while (i < (b_len / unroll) * unroll) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // py[i + u] += temp1 * a[tile_i + i + u + j * lda]
                            numeric.fma_(
                                &py[i + u],
                                temp1,
                                if (comptime noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                py[i + u],
                            );

                            // sums[u] += conj(a[tile_i + i + u + j * lda]) * px[i + u]
                            numeric.fma_(
                                &sums[u],
                                if (comptime !noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                px[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        numeric.add_(&temp2, temp2, sums[u]);
                    }

                    while (i < b_len) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fma_(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) * px[i]
                        numeric.fma_(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    // y[jy] += alpha * temp2
                    numeric.fma_(
                        &y[numeric.cast(usize, jy)],
                        alpha,
                        temp2,
                        y[numeric.cast(usize, jy)],
                    );
                }

                jx += incx;
                jy += incy;
            }

            if (incy != 1) {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    @as([*]Y, &local_y),
                    1,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                );
            }
        }
    }

    return;
}
