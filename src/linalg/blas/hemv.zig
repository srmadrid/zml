const std = @import("std");
const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

const k_gemv = @import("gemv.zig").k_gemv;

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

    var threads: [options.max_threads]std.Thread = undefined;

    const Worker = struct {
        fn execute(worker_uplo: Uplo, worker_n: usize, worker_chunk_start: usize, worker_chunk_end: usize, worker_alpha: Al, worker_a: [*]const A, worker_lda: usize, worker_x: [*]const X, worker_incx: isize, worker_beta: Be, worker_y: [*]Y, worker_incy: isize, worker_noconj: bool) void {
            // chunk = chunk_start..chunk_end
            const chunk_len = worker_chunk_end - worker_chunk_start;

            // Sub-vector pointers for x[chunk] and y[chunk].
            const y_chunk_ptr = worker_y + numeric.cast(usize, if (worker_incy > 0)
                numeric.cast(isize, worker_chunk_start) * worker_incy
            else
                (numeric.cast(isize, worker_chunk_end) - numeric.cast(isize, worker_n)) * worker_incy);
            const x_chunk_ptr = worker_x + numeric.cast(usize, if (worker_incx > 0)
                numeric.cast(isize, worker_chunk_start) * worker_incx
            else
                (numeric.cast(isize, worker_chunk_end) - numeric.cast(isize, worker_n)) * worker_incx);

            // Hemv on the diagonal block A[chunk, chunk].
            // y[chunk] = beta * y[chunk] + alpha * A[chunk, chunk_end..n] * x[chunk_end..n].
            if (worker_noconj)
                k_hemv(
                    worker_uplo,
                    chunk_len,
                    worker_alpha,
                    worker_a + worker_chunk_start + worker_chunk_start * worker_lda,
                    worker_lda,
                    x_chunk_ptr,
                    worker_incx,
                    worker_beta,
                    y_chunk_ptr,
                    worker_incy,
                    true,
                )
            else
                k_hemv(
                    worker_uplo,
                    chunk_len,
                    worker_alpha,
                    worker_a + worker_chunk_start + worker_chunk_start * worker_lda,
                    worker_lda,
                    x_chunk_ptr,
                    worker_incx,
                    worker_beta,
                    y_chunk_ptr,
                    worker_incy,
                    false,
                );

            if (worker_uplo == .upper) {
                // Gemv on the right rectangular block A[chunk, chunk_end..n].
                // y[chunk] += alpha * A[chunk, chunk_end..n] * x[chunk_end..n].
                if (worker_chunk_end < worker_n) {
                    const right_len = worker_n - worker_chunk_end;
                    const x_right_ptr = worker_x + numeric.cast(usize, if (worker_incx > 0)
                        numeric.cast(isize, worker_chunk_end) * worker_incx
                    else
                        0);

                    if (worker_noconj)
                        k_gemv(
                            .no_trans,
                            chunk_len,
                            right_len,
                            worker_alpha,
                            worker_a + worker_chunk_start + worker_chunk_end * worker_lda,
                            worker_lda,
                            x_right_ptr,
                            worker_incx,
                            numeric.one(Be),
                            y_chunk_ptr,
                            worker_incy,
                            true,
                        )
                    else
                        k_gemv(
                            .conj_no_trans,
                            chunk_len,
                            right_len,
                            worker_alpha,
                            worker_a + worker_chunk_start + worker_chunk_end * worker_lda,
                            worker_lda,
                            x_right_ptr,
                            worker_incx,
                            numeric.one(Be),
                            y_chunk_ptr,
                            worker_incy,
                            false,
                        );
                }

                // Gemv on the top rectangular block A[0..chunk_start, chunk].
                // y[chunk] += alpha * (A[0..chunk_start, chunk])ᴴ * x[0..chunk_start].
                if (worker_chunk_start > 0) {
                    const top_len = worker_chunk_start;
                    const x_top_ptr = worker_x + numeric.cast(usize, if (worker_incx > 0)
                        0
                    else
                        (numeric.cast(isize, top_len) - numeric.cast(isize, worker_n)) * worker_incx);

                    if (worker_noconj)
                        k_gemv(
                            .conj_trans,
                            top_len,
                            chunk_len,
                            worker_alpha,
                            worker_a + worker_chunk_start * worker_lda,
                            worker_lda,
                            x_top_ptr,
                            worker_incx,
                            numeric.one(Be),
                            y_chunk_ptr,
                            worker_incy,
                            false,
                        )
                    else
                        k_gemv(
                            .trans,
                            top_len,
                            chunk_len,
                            worker_alpha,
                            worker_a + worker_chunk_start * worker_lda,
                            worker_lda,
                            x_top_ptr,
                            worker_incx,
                            numeric.one(Be),
                            y_chunk_ptr,
                            worker_incy,
                            true,
                        );
                }
            } else {
                // Gemv on the left rectangular block A[chunk, 0..chunk_start].
                // y[chunk] += alpha * A[chunk, 0..chunk_start] * x[0..chunk_start].
                if (worker_chunk_start > 0) {
                    const left_len = worker_chunk_start;
                    const x_left_ptr = worker_x + numeric.cast(usize, if (worker_incx > 0)
                        0
                    else
                        (numeric.cast(isize, left_len) - numeric.cast(isize, worker_n)) * worker_incx);

                    if (worker_noconj)
                        k_gemv(
                            .no_trans,
                            chunk_len,
                            left_len,
                            worker_alpha,
                            worker_a + worker_chunk_start,
                            worker_lda,
                            x_left_ptr,
                            worker_incx,
                            numeric.one(Be),
                            y_chunk_ptr,
                            worker_incy,
                            true,
                        )
                    else
                        k_gemv(
                            .conj_no_trans,
                            chunk_len,
                            left_len,
                            worker_alpha,
                            worker_a + worker_chunk_start,
                            worker_lda,
                            x_left_ptr,
                            worker_incx,
                            numeric.one(Be),
                            y_chunk_ptr,
                            worker_incy,
                            false,
                        );
                }

                // Gemv on the bottom rectangular block A[chunk_end..n, chunk].
                // y[chunk] += alpha * (A[chunk_end..n, chunk])ᴴ * x[chunk_end..n].
                if (worker_chunk_end < worker_n) {
                    const bottom_len = worker_n - worker_chunk_end;
                    const x_bottom_ptr = worker_x + numeric.cast(usize, if (worker_incx > 0)
                        numeric.cast(isize, worker_chunk_end) * worker_incx
                    else
                        0);

                    if (worker_noconj)
                        k_gemv(
                            .conj_trans,
                            bottom_len,
                            chunk_len,
                            worker_alpha,
                            worker_a + worker_chunk_end + worker_chunk_start * worker_lda,
                            worker_lda,
                            x_bottom_ptr,
                            worker_incx,
                            numeric.one(Be),
                            y_chunk_ptr,
                            worker_incy,
                            false,
                        )
                    else
                        k_gemv(
                            .trans,
                            bottom_len,
                            chunk_len,
                            worker_alpha,
                            worker_a + worker_chunk_end + worker_chunk_start * worker_lda,
                            worker_lda,
                            x_bottom_ptr,
                            worker_incx,
                            numeric.one(Be),
                            y_chunk_ptr,
                            worker_incy,
                            true,
                        );
                }
            }
        }
    };

    const chunk_size = int.div(n, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) n else chunk_start + chunk_size;

        if (std.Thread.spawn(.{}, Worker.execute, .{
            eff_uplo,
            n,
            chunk_start,
            chunk_end,
            alpha,
            a,
            lda,
            x,
            incx,
            beta,
            y,
            incy,
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
