const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Performs a rank-2 update of a Hermitian packed matrix.
///
/// The `hpr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
/// A = alpha * x * yᴴ + conj(alpha) * y * xᴴ + A,
/// ```
///
/// where `alpha` is a numeric, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.hpr2(layout: Layout, uplo: Uplo, n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]const Y, incy: isize, ap: [*]A) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of
///   the Hermitian packed matrix `A` is used.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `ap` (`anytype`): Mutable many-item pointer, size at least
///   `(n * (n + 1)) / 2`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is 0.
pub fn hpr2(layout: Layout, uplo: Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, y: anytype, incy: isize, ap: anytype) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var Ap: type = @TypeOf(ap);

    comptime if (!meta.isNumeric(Al) or !meta.isReal(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or !meta.isNumeric(meta.Child(Y)) or
        !meta.isManyItemPointer(Ap) or meta.isConstPointer(Ap) or !meta.isNumeric(meta.Child(Ap)))
        @compileError("zsl.linalg.blas.hpr2: alpha must be a numeric, x and y must be many-item pointers to numerics, and ap must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tap: " ++ @typeName(Ap) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);
    Ap = meta.Child(Ap);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Al == X and Al == Y and Al == Ap) {
        switch (comptime meta.numericType(X)) {
            .complex => {
                if (comptime meta.Scalar(Al) == f32)
                    return linalg.cblas.chpr2(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), alpha, x, incx, y, incy, ap)
                else if (comptime meta.Scalar(Al) == f64)
                    return linalg.cblas.zhpr2(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), alpha, x, incx, y, incy, ap);
            },
            else => {},
        }
    }

    return if (layout == .col_major)
        k_hpr2(uplo, n, alpha, x, incx, y, incy, ap, true)
    else
        k_hpr2(uplo.invert(), n, alpha, x, incx, y, incy, ap, false);
}

fn k_hpr2(uplo: Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, y: anytype, incy: isize, ap: anytype, comptime noconj: bool) !void {
    // Quick return if possible.
    if (n == 0 or numeric.eq(alpha, 0))
        return;

    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;

    var kk: usize = 0;
    if (uplo == .upper) {
        // Form  A  when upper triangle of A is stored.
        if (incx == 1 and incy == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[j], 0) or numeric.ne(y[j], 0)) {
                    // temp1 = alpha * conj(y[j])
                    const temp1 = numeric.mul(
                        alpha,
                        if (comptime !noconj)
                            y[j]
                        else
                            numeric.conj(y[j]),
                    );

                    // temp2 = conj(alpha * x[j])
                    const temp2 = numeric.conj(
                        numeric.mul(
                            alpha,
                            if (comptime noconj)
                                x[j]
                            else
                                numeric.conj(x[j]),
                        ),
                    );

                    var k: usize = kk;
                    var i: usize = 0;
                    while (i < j) : (i += 1) {
                        // ap[k] += x[i] * temp1 + y[i] * temp2
                        numeric.fma_(
                            &ap[k],
                            if (comptime noconj)
                                x[i]
                            else
                                numeric.conj(x[i]),
                            temp1,
                            numeric.fma(
                                if (comptime noconj)
                                    y[i]
                                else
                                    numeric.conj(y[i]),
                                temp2,
                                ap[k],
                            ),
                        );

                        k += 1;
                    }

                    // ap[kk + j] = re(ap[kk + j]) + re(x[j] * temp1 + y[j] * temp2)
                    numeric.fma_(
                        &ap[kk + j],
                        if (comptime !noconj)
                            numeric.im(y[j])
                        else
                            numeric.neg(numeric.im(y[j])),
                        numeric.im(temp2),
                        numeric.fma(
                            if (comptime !noconj)
                                numeric.im(x[j])
                            else
                                numeric.neg(numeric.im(x[j])),
                            numeric.im(temp1),
                            numeric.fma(
                                numeric.re(y[j]),
                                numeric.re(temp2),
                                numeric.fma(
                                    numeric.re(x[j]),
                                    numeric.re(temp1),
                                    numeric.re(ap[kk + j]),
                                ),
                            ),
                        ),
                    );
                } else {
                    numeric.set( // ap[kk + j] = re(ap[kk + j])
                        &ap[kk + j],
                        numeric.re(ap[kk + j]),
                    );
                }

                kk += j + 1;
            }
        } else {
            var jx: isize = kx;
            var jy: isize = ky;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0) or numeric.ne(y[numeric.cast(usize, jy)], 0)) {
                    // temp1 = alpha * conj(y[jy])
                    const temp1 = numeric.mul(
                        alpha,
                        if (comptime !noconj)
                            y[numeric.cast(usize, jy)]
                        else
                            numeric.conj(y[numeric.cast(usize, jy)]),
                    );

                    // temp2 = conj(alpha * x[jx])
                    const temp2 = numeric.conj(
                        numeric.mul(
                            alpha,
                            if (comptime noconj)
                                x[numeric.cast(usize, jx)]
                            else
                                numeric.conj(x[numeric.cast(usize, jx)]),
                        ),
                    );

                    var ix: isize = kx;
                    var iy: isize = ky;
                    var k: usize = kk;
                    while (k < kk + j) : (k += 1) {
                        // ap[k] += x[ix] * temp1 + y[iy] * temp2
                        numeric.fma_(
                            &ap[k],
                            if (comptime noconj)
                                x[numeric.cast(usize, ix)]
                            else
                                numeric.conj(x[numeric.cast(usize, ix)]),
                            temp1,
                            numeric.fma(
                                if (comptime noconj)
                                    y[numeric.cast(usize, iy)]
                                else
                                    numeric.conj(y[numeric.cast(usize, iy)]),
                                temp2,
                                ap[k],
                            ),
                        );

                        ix += incx;
                        iy += incy;
                    }

                    // ap[kk + j] = re(ap[kk + j]) + re(x[jx] * temp1 + y[jy] * temp2)
                    numeric.fma_(
                        &ap[kk + j],
                        if (comptime !noconj)
                            numeric.im(y[numeric.cast(usize, jy)])
                        else
                            numeric.neg(numeric.im(y[numeric.cast(usize, jy)])),
                        numeric.im(temp2),
                        numeric.fma(
                            if (comptime !noconj)
                                numeric.im(x[numeric.cast(usize, jx)])
                            else
                                numeric.neg(numeric.im(x[numeric.cast(usize, jx)])),
                            numeric.im(temp1),
                            numeric.fma(
                                numeric.re(y[numeric.cast(usize, jy)]),
                                numeric.re(temp2),
                                numeric.fma(
                                    numeric.re(x[numeric.cast(usize, jx)]),
                                    numeric.re(temp1),
                                    numeric.re(ap[kk + j]),
                                ),
                            ),
                        ),
                    );
                } else {
                    numeric.set( // ap[kk + j] = re(ap[kk + j])
                        &ap[kk + j],
                        numeric.re(ap[kk + j]),
                    );
                }

                jx += incx;
                jy += incy;
                kk += j + 1;
            }
        }
    } else {
        if (incx == 1 and incy == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[j], 0) or numeric.ne(y[j], 0)) {
                    // temp1 = alpha * conj(y[j])
                    const temp1 = numeric.mul(
                        alpha,
                        if (comptime !noconj)
                            y[j]
                        else
                            numeric.conj(y[j]),
                    );

                    // temp2 = conj(alpha * x[j])
                    const temp2 = numeric.conj(
                        numeric.mul(
                            alpha,
                            if (comptime noconj)
                                x[j]
                            else
                                numeric.conj(x[j]),
                        ),
                    );

                    // ap[kk] = re(ap[kk]) + re(x[j] * temp1 + y[j] * temp2)
                    numeric.fma_(
                        &ap[kk],
                        if (comptime !noconj)
                            numeric.im(y[j])
                        else
                            numeric.neg(numeric.im(y[j])),
                        numeric.im(temp2),
                        numeric.fma(
                            if (comptime !noconj)
                                numeric.im(x[j])
                            else
                                numeric.neg(numeric.im(x[j])),
                            numeric.im(temp1),
                            numeric.fma(
                                numeric.re(y[j]),
                                numeric.re(temp2),
                                numeric.fma(
                                    numeric.re(x[j]),
                                    numeric.re(temp1),
                                    numeric.re(ap[kk]),
                                ),
                            ),
                        ),
                    );

                    var k: usize = kk + 1;
                    var i: usize = j + 1;
                    while (i < n) : (i += 1) {
                        // ap[k] += x[i] * temp1 + y[i] * temp2
                        numeric.fma_(
                            &ap[k],
                            if (comptime noconj)
                                x[i]
                            else
                                numeric.conj(x[i]),
                            temp1,
                            numeric.fma(
                                if (comptime noconj)
                                    y[i]
                                else
                                    numeric.conj(y[i]),
                                temp2,
                                ap[k],
                            ),
                        );

                        k += 1;
                    }
                } else {
                    numeric.set( // ap[kk] = re(ap[kk])
                        &ap[kk],
                        numeric.re(ap[kk]),
                    );
                }

                kk += n - j;
            }
        } else {
            var jx: isize = kx;
            var jy: isize = ky;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0) or numeric.ne(y[numeric.cast(usize, jy)], 0)) {
                    // temp1 = alpha * conj(y[jy])
                    const temp1 = numeric.mul(
                        alpha,
                        if (comptime !noconj)
                            y[numeric.cast(usize, jy)]
                        else
                            numeric.conj(y[numeric.cast(usize, jy)]),
                    );

                    // temp2 = conj(alpha * x[jx])
                    const temp2 = numeric.conj(
                        numeric.mul(
                            alpha,
                            if (comptime noconj)
                                x[numeric.cast(usize, jx)]
                            else
                                numeric.conj(x[numeric.cast(usize, jx)]),
                        ),
                    );

                    // ap[kk] = re(ap[kk]) + re(x[jx] * temp1 + y[jy] * temp2)
                    numeric.fma_(
                        &ap[kk],
                        if (comptime !noconj)
                            numeric.im(y[numeric.cast(usize, jy)])
                        else
                            numeric.neg(numeric.im(y[numeric.cast(usize, jy)])),
                        numeric.im(temp2),
                        numeric.fma(
                            if (comptime !noconj)
                                numeric.im(x[numeric.cast(usize, jx)])
                            else
                                numeric.neg(numeric.im(x[numeric.cast(usize, jx)])),
                            numeric.im(temp1),
                            numeric.fma(
                                numeric.re(y[numeric.cast(usize, jy)]),
                                numeric.re(temp2),
                                numeric.fma(
                                    numeric.re(x[numeric.cast(usize, jx)]),
                                    numeric.re(temp1),
                                    numeric.re(ap[kk]),
                                ),
                            ),
                        ),
                    );

                    var ix: isize = jx;
                    var iy: isize = jy;
                    var k: usize = kk + 1;
                    while (k < kk + n - j) : (k += 1) {
                        ix += incx;
                        iy += incy;

                        // ap[k] += x[ix] * temp1 + y[iy] * temp2
                        numeric.fma_(
                            &ap[k],
                            if (comptime noconj)
                                x[numeric.cast(usize, ix)]
                            else
                                numeric.conj(x[numeric.cast(usize, ix)]),
                            temp1,
                            numeric.fma(
                                if (comptime noconj)
                                    y[numeric.cast(usize, iy)]
                                else
                                    numeric.conj(y[numeric.cast(usize, iy)]),
                                temp2,
                                ap[k],
                            ),
                        );
                    }
                } else {
                    numeric.set( // ap[kk] = re(ap[kk])
                        &ap[kk],
                        numeric.re(ap[kk]),
                    );
                }

                jx += incx;
                jy += incy;
                kk += n - j;
            }
        }
    }

    return;
}
