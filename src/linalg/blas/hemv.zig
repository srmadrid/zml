const std = @import("std");
const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");

const linalg = @import("../../linalg.zig");

/// Size of tiles to use when multithreading.
const tile_size = 128;

/// Computes a matrix-vector product with a Hermitian matrix defined as:
///
/// ```zig
/// y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are numerics, `x` and `y` are `n`-element vectors,
/// and `A` is an `n`-by-`n` Hermitian matrix.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
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
///   the Hermitian packed matrix `A` is used.
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
///       threads = max(1, min(std.Thread.getCpuCount(), options.max_threads, (m * n) / parallel_threshold))
///       ```
///     * 1: force serial execution. parallel_threshold is ignored.
///     * N >= 2: use exactly N threads, clamped by
///       std.Thread.getCpuCount() and options.max_threads as a hard safety
///       ceiling. parallel_threshold is ignored.
///   * parallel_threshold (usize = 4_194_304 / @sizeOf(meta.Child(Y))):
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
        parallel_threshold: usize = 4_194_304 / @sizeOf(meta.Child(@TypeOf(y))),
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

    if (comptime options.link_cblas != null and Al == A and Al == X and Al == Y and Al == Be) {
        switch (comptime meta.numericType(Al)) {
            .complex => {
                if (comptime meta.Scalar(Al) == f32)
                    return linalg.cblas.chemv(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), &alpha, a, numeric.cast(isize, lda), x, incx, &beta, y, incy)
                else if (comptime meta.Scalar(Al) == f64)
                    return linalg.cblas.zhemv(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), &alpha, a, numeric.cast(isize, lda), x, incx, &beta, y, incy);
            },
            else => {},
        }
    }

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();
    const noconj = layout == .col_major;

    // Quick return if possible.
    if (n == 0)
        return;

    if (numeric.eq(alpha, 0)) {
        if (numeric.ne(beta, 1))
            @import("scal.zig").k_scal(n, beta, y, incy);

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
        @import("scal.zig").k_scal(n, beta, y, incy);

    if (numeric.eq(alpha, 0))
        return;

    const k = (n + tile_size - 1) / tile_size;
    const num_tiles = k * (k + 1) / 2;

    var atomic_counter = std.atomic.Value(usize).init(0);
    var threads: [options.max_threads]std.Thread = undefined;

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

                    @import("copy.zig").k_copy(
                        r_len,
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, r_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incx),
                        worker_incx,
                        @as([*]X, &local_x),
                        1,
                    );

                    k_hemv(
                        worker_uplo,
                        r_len,
                        worker_alpha,
                        worker_a + r_start + c_start * worker_lda,
                        worker_lda,
                        @as([]const X, &local_x),
                        1,
                        numeric.one(Be),
                        @as([]Y, &local_y),
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

                    @import("copy.zig").k_copy(
                        r_len,
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, r_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incx),
                        worker_incx,
                        @as([*]X, &local_x_r),
                        1,
                    );

                    @import("copy.zig").k_copy(
                        c_len,
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, c_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, c_start + c_len)) * worker_incx),
                        worker_incx,
                        @as([*]X, &local_x_c),
                        1,
                    );

                    var j: usize = 0;
                    while (j < c_len) : (j += 1) {
                        // temp1 = worker_alpha * local_x_c[j]
                        const temp1 = numeric.mul(worker_alpha, local_x_c[j]);

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
                                    local_y_r[i],
                                );

                                // sums[u] += conj(worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda]) * local_x_r[i + u]
                                numeric.fma_(
                                    &sums[u],
                                    if (comptime worker_noconj)
                                        numeric.conj(worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda])
                                    else
                                        worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda],
                                    local_x_r[i + u],
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
                                local_x_r[i],
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
        @import("scal.zig").k_scal(n, beta, y, incy);

    if (numeric.eq(alpha, 0))
        return;

    if (uplo == .upper) {
        const unroll = 2 * int.min(
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2,
            std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X))) orelse 2,
        );

        // Form  y  when upper triangle of A is stored.
        if (incx == 1 and incy == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                // temp1 = alpha * x[j]
                const temp1 = numeric.mul(
                    alpha,
                    x[j],
                );

                var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                var i: usize = 0;
                while (i < (j / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // y[i + u] += temp1 * a[i + u + j * lda]
                        numeric.fma_(
                            &y[i + u],
                            temp1,
                            if (comptime noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            y[i + u],
                        );

                        // temp2 += conj(a[i + u + j * lda]) * x[i + u]
                        numeric.fma_(
                            &temp2,
                            if (comptime !noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            x[i + u],
                            temp2,
                        );
                    }
                }

                while (i < j) : (i += 1) {
                    // y[i] += temp1 * a[i + j * lda]
                    numeric.fma_(
                        &y[i],
                        temp1,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        y[i],
                    );

                    // temp2 += conj(a[i + j * lda]) * x[i]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        x[i],
                        temp2,
                    );
                }

                // y[j] += temp1 * re(a[j + j * lda]) + alpha * temp2
                numeric.fma_(
                    &y[j],
                    temp1,
                    numeric.re(a[j + j * lda]),
                    numeric.fma(
                        alpha,
                        temp2,
                        y[j],
                    ),
                );
            }
        } else {
            var jx: isize = kx;
            var jy: isize = ky;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                // temp1 = alpha * x[jx]
                const temp1 = numeric.mul(
                    alpha,
                    x[numeric.cast(usize, jx)],
                );

                var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                var ix: isize = kx;
                var iy: isize = ky;
                var i: usize = 0;
                while (i < (j / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // y[iy + u * incy] += temp1 * a[i + u + j * lda]
                        numeric.fma_(
                            &y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                            temp1,
                            if (comptime noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                        );

                        // temp2 += conj(a[i + u + j * lda]) * x[ix + u * incx]
                        numeric.fma_(
                            &temp2,
                            if (comptime !noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                            temp2,
                        );
                    }

                    ix += numeric.cast(isize, unroll) * incx;
                    iy += numeric.cast(isize, unroll) * incy;
                }

                while (i < j) : (i += 1) {
                    // y[iy] += temp1 * a[i + j * lda]
                    numeric.fma_(
                        &y[numeric.cast(usize, iy)],
                        temp1,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        y[numeric.cast(usize, iy)],
                    );

                    // temp2 += conj(a[i + j * lda]) * x[ix]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        x[numeric.cast(usize, ix)],
                        temp2,
                    );

                    ix += incx;
                    iy += incy;
                }

                // y[jy] += temp1 * re(a[j + j * lda]) + alpha * temp2
                numeric.fma_(
                    &y[numeric.cast(usize, jy)],
                    temp1,
                    numeric.re(a[j + j * lda]),
                    numeric.fma(
                        alpha,
                        temp2,
                        y[numeric.cast(usize, jy)],
                    ),
                );

                jx += incx;
                jy += incy;
            }
        }
    } else {
        const unroll = 2 * int.min(
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2,
            std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X))) orelse 2,
        );

        // Form  y  when lower triangle of A is stored.
        if (incx == 1 and incy == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                // temp1 = alpha * x[j]
                const temp1 = numeric.mul(
                    alpha,
                    x[j],
                );

                var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                // y[j] += temp1 * re(a[j + j * lda])
                numeric.fma_(
                    &y[j],
                    temp1,
                    numeric.re(a[j + j * lda]),
                    y[j],
                );

                var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                var i: usize = j + 1;
                while (i < n and i % unroll != 0) : (i += 1) {
                    // y[i] += temp1 * a[i + j * lda]
                    numeric.fma_(
                        &y[i],
                        temp1,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        y[i],
                    );

                    // temp2 += conj(a[i + j * lda]) * x[i]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        x[i],
                        temp2,
                    );
                }

                while (i < (n / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // y[i + u] += temp1 * a[i + u + j * lda]
                        numeric.fma_(
                            &y[i + u],
                            temp1,
                            if (comptime noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            y[i + u],
                        );

                        // sums[u] += conj(a[i + u + j * lda]) * x[i + u]
                        numeric.fma_(
                            &sums[u],
                            if (comptime !noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            x[i + u],
                            sums[u],
                        );
                    }
                }

                inline for (0..unroll) |u| {
                    numeric.add_(&temp2, temp2, sums[u]);
                }

                while (i < n) : (i += 1) {
                    // y[i] += temp1 * a[i + j * lda]
                    numeric.fma_(
                        &y[i],
                        temp1,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        y[i],
                    );

                    // temp2 += conj(a[i + j * lda]) * x[i]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        x[i],
                        temp2,
                    );
                }

                // y[j] += alpha * temp2
                numeric.fma_(
                    &y[j],
                    alpha,
                    temp2,
                    y[j],
                );
            }
        } else {
            var jx: isize = kx;
            var jy: isize = ky;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                // temp1 = alpha * x[jx]
                const temp1 = numeric.mul(
                    alpha,
                    x[numeric.cast(usize, jx)],
                );

                var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                // y[jy] += temp1 * re(a[j + j * lda])
                numeric.fma_(
                    &y[numeric.cast(usize, jy)],
                    temp1,
                    numeric.re(a[j + j * lda]),
                    y[numeric.cast(usize, jy)],
                );

                var ix: isize = jx + incx;
                var iy: isize = jy + incy;
                var i: usize = j + 1;
                while (i < n and i % unroll != 0) : (i += 1) {
                    // y[iy] += temp1 * a[i + j * lda]
                    numeric.fma_(
                        &y[numeric.cast(usize, iy)],
                        temp1,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        y[numeric.cast(usize, iy)],
                    );

                    // temp2 += conj(a[i + j * lda]) * x[ix]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        x[numeric.cast(usize, ix)],
                        temp2,
                    );

                    ix += incx;
                    iy += incy;
                }

                var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                while (i < (n / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // y[iy + u * incy] += temp1 * a[i + u + j * lda]
                        numeric.fma_(
                            &y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                            temp1,
                            if (comptime noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                        );

                        // sums[u] += conj(a[i + u + j * lda]) * x[ix + u * incx]
                        numeric.fma_(
                            &sums[u],
                            if (comptime !noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                            sums[u],
                        );
                    }

                    ix += numeric.cast(isize, unroll) * incx;
                    iy += numeric.cast(isize, unroll) * incy;
                }

                inline for (0..unroll) |u| {
                    numeric.add_(&temp2, temp2, sums[u]);
                }

                while (i < n) : (i += 1) {
                    // y[iy] += temp1 * a[i + j * lda]
                    numeric.fma_(
                        &y[numeric.cast(usize, iy)],
                        temp1,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        y[numeric.cast(usize, iy)],
                    );

                    // temp2 += conj(a[i + j * lda]) * x[ix]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        x[numeric.cast(usize, ix)],
                        temp2,
                    );

                    ix += incx;
                    iy += incy;
                }

                // y[jy] += alpha * temp2
                numeric.fma_(
                    &y[numeric.cast(usize, jy)],
                    alpha,
                    temp2,
                    y[numeric.cast(usize, jy)],
                );

                jx += incx;
                jy += incy;
            }
        }
    }

    return;
}
