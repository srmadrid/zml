const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Computes a matrix-vector product with a Hermitian packed matrix defined as:
///
/// ```zig
/// y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are numerics, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.hpmv(layout: Layout, uplo: Uplo, n: usize, alpha: Al, ap: [*]const Ap, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of
///   the Hermitian band matrix `A` is used.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `ap` (`anytype`): Many-item pointer, size at least `(n * (n + 1)) / 2`.
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
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is 0.
pub fn hpmv(layout: Layout, uplo: Uplo, n: usize, alpha: anytype, ap: anytype, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize) !void {
    const Al: type = @TypeOf(alpha);
    comptime var Ap: type = @TypeOf(ap);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(Ap) or !meta.isNumeric(meta.Child(Ap)) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.hpmv: alpha and beta must be numerics, ap and x must be many-item pointers to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tap: " ++ @typeName(Ap) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    Ap = meta.Child(Ap);
    X = meta.Child(X);
    Y = meta.Child(Y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Al == Ap and Al == X and Al == Y and Al == Be) {
        switch (comptime meta.numericType(Ap)) {
            .complex => {
                if (comptime meta.Scalar(Ap) == f32)
                    return linalg.cblas.chpmv(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), &alpha, ap, x, incx, &beta, y, incy)
                else if (comptime meta.Scalar(Ap) == f64)
                    return linalg.cblas.zhpmv(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), &alpha, ap, x, incx, &beta, y, incy);
            },
            else => {},
        }
    }

    return if (layout == .col_major)
        k_hpmv(uplo, n, alpha, ap, x, incx, beta, y, incy, true)
    else
        k_hpmv(uplo.invert(), n, alpha, ap, x, incx, beta, y, incy, false);
}

fn k_hpmv(uplo: Uplo, n: usize, alpha: anytype, ap: anytype, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize, comptime noconj: bool) !void {
    const Ap: type = meta.Child(@TypeOf(ap));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (n == 0 or (numeric.eq(alpha, 0) and numeric.eq(beta, 1)))
        return;

    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;

    // First form  y = beta * y.
    if (numeric.ne(beta, 1)) {
        if (incy == 1) {
            if (numeric.eq(beta, 0)) {
                var i: usize = 0;
                while (i < n) : (i += 1) {
                    // y[i] = 0
                    numeric.set(&y[i], 0);
                }
            } else {
                var i: usize = 0;
                while (i < n) : (i += 1) {
                    // y[i] *= beta
                    numeric.mul_(
                        &y[i],
                        y[i],
                        beta,
                    );
                }
            }
        } else {
            var iy: isize = ky;
            if (numeric.eq(beta, 0)) {
                var i: usize = 0;
                while (i < n) : (i += 1) {
                    // y[iy] = 0
                    numeric.set(&y[numeric.cast(usize, iy)], 0);

                    iy += incy;
                }
            } else {
                var i: usize = 0;
                while (i < n) : (i += 1) {
                    // y[iy] *= beta
                    numeric.mul_(
                        &y[numeric.cast(usize, iy)],
                        y[numeric.cast(usize, iy)],
                        beta,
                    );

                    iy += incy;
                }
            }
        }
    }

    if (numeric.eq(alpha, 0))
        return;

    var kk: usize = 0;
    if (uplo == .upper) {
        // Form  y  when upper triangle of A is stored.
        if (incx == 1 and incy == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                // temp1 = alpha * x[j]
                const temp1 = numeric.mul(
                    alpha,
                    x[j],
                );

                var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));

                var k: usize = kk;
                var i: usize = 0;
                while (i < j) : (i += 1) {
                    // y[i] += temp1 * ap[k]
                    numeric.fma_(
                        &y[i],
                        temp1,
                        if (comptime noconj)
                            ap[k]
                        else
                            numeric.conj(ap[k]),
                        y[i],
                    );

                    // temp2 += conj(ap[k]) * x[i]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            ap[k]
                        else
                            numeric.conj(ap[k]),
                        x[i],
                        temp2,
                    );

                    k += 1;
                }

                // y[j] += temp1 * re(ap[kk + j]) + alpha * temp2
                numeric.fma_(
                    &y[j],
                    temp1,
                    numeric.re(ap[kk + j]),
                    numeric.fma(
                        alpha,
                        temp2,
                        y[j],
                    ),
                );

                kk += j + 1;
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

                var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));

                var ix: isize = kx;
                var iy: isize = ky;
                var k: isize = kk;
                while (k < kk + j) : (k += 1) {
                    // y[iy] += temp1 * ap[k]
                    numeric.fma_(
                        &y[numeric.cast(usize, iy)],
                        temp1,
                        if (comptime noconj)
                            ap[k]
                        else
                            numeric.conj(ap[k]),
                        y[numeric.cast(usize, iy)],
                    );

                    // temp2 += conj(ap[k]) * x[ix]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            ap[k]
                        else
                            numeric.conj(ap[k]),
                        x[numeric.cast(usize, ix)],
                        temp2,
                    );

                    ix += incx;
                    iy += incy;
                }

                // y[jy] += temp1 * re(ap[kk + j]) + alpha * temp2
                numeric.fma_(
                    &y[numeric.cast(usize, jy)],
                    temp1,
                    numeric.re(ap[kk + j]),
                    numeric.fma(
                        alpha,
                        temp2,
                        y[numeric.cast(usize, jy)],
                    ),
                );

                jx += incx;
                jy += incy;
                kk += j + 1;
            }
        }
    } else {
        // Form  y  when lower triangle of A is stored.
        if (incx == 1 and incy == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                // temp1 = alpha * x[j]
                const temp1 = numeric.mul(
                    alpha,
                    x[j],
                );

                var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));

                // y[j] += temp1 * re(ap[kk])
                numeric.fma_(
                    &y[j],
                    temp1,
                    numeric.re(ap[kk]),
                    y[j],
                );

                var k: usize = kk + 1;
                var i: usize = j + 1;
                while (i < n) : (i += 1) {
                    // y[i] += temp1 * ap[k]
                    numeric.fma_(
                        &y[i],
                        temp1,
                        if (comptime noconj)
                            ap[k]
                        else
                            numeric.conj(ap[k]),
                        y[i],
                    );

                    // temp2 += conj(ap[k]) * x[i]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            ap[k]
                        else
                            numeric.conj(ap[k]),
                        x[i],
                        temp2,
                    );

                    k += 1;
                }

                // y[j] += alpha * temp2
                numeric.fma_(
                    &y[j],
                    alpha,
                    temp2,
                    y[j],
                );

                kk += n - j;
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

                var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));

                // y[jy] += temp1 * re(ap[kk])
                numeric.fma_(
                    &y[numeric.cast(usize, jy)],
                    temp1,
                    numeric.re(ap[kk]),
                    y[numeric.cast(usize, jy)],
                );

                var ix: isize = jx;
                var iy: isize = jy;
                var k: usize = kk + 1;
                while (k < kk + n - j) : (k += 1) {
                    ix += incx;
                    iy += incy;

                    // y[iy] += temp1 * ap[k]
                    numeric.fma_(
                        &y[numeric.cast(usize, iy)],
                        temp1,
                        if (comptime noconj)
                            ap[k]
                        else
                            numeric.conj(ap[k]),
                        y[numeric.cast(usize, iy)],
                    );

                    // temp2 += conj(ap[k]) * x[ix]
                    numeric.fma_(
                        &temp2,
                        if (comptime !noconj)
                            ap[k]
                        else
                            numeric.conj(ap[k]),
                        x[numeric.cast(usize, ix)],
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

                jx += incx;
                jy += incy;
                kk += n - j;
            }
        }
    }

    return;
}
