const options = @import("options");

const meta = @import("../../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

/// Performs a rank-2 update of a symmetric packed matrix.
///
/// The `spr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
/// A = alpha * x * yᵀ + alpha * y * xᵀ + A,
/// ```
///
/// where `alpha` is a numeric, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.spr2(layout: Layout, uplo: Uplo, n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]const Y, incy: isize, ap: [*]A) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of
///   the symmetric packed matrix `A` is used.
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
pub fn spr2(layout: Layout, uplo: Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, y: anytype, incy: isize, ap: anytype) !void {
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
            .float => {
                if (comptime Al == f32)
                    return linalg.cblas.sspr2(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), alpha, x, incx, y, incy, ap)
                else if (comptime Al == f64)
                    return linalg.cblas.dspr2(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), alpha, x, incx, y, incy, ap);
            },
            else => {},
        }
    }

    return if (layout == .col_major)
        k_spr2(uplo, n, alpha, x, incx, y, incy, ap)
    else
        k_spr2(uplo.invert(), n, alpha, x, incx, y, incy, ap);
}

fn k_spr2(uplo: Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, y: anytype, incy: isize, ap: anytype) void {
    // Quick return if possible.
    if (n == 0 or numeric.eq(alpha, 0))
        return;

    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;

    var kk: usize = 0;
    if (uplo == .upper) {
        if (incx == 1 and incy == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[j], 0) or numeric.ne(y[j], 0)) {
                    // temp1 = alpha * y[j]
                    const temp1 = numeric.mul(
                        alpha,
                        y[j],
                    );

                    // temp2 = alpha * x[j]
                    const temp2 = numeric.mul(
                        alpha,
                        x[j],
                    );

                    var k: usize = kk;
                    var i: usize = 0;
                    while (i < j) : (i += 1) {
                        // ap[k] += x[i] * temp1 + y[i] * temp2
                        numeric.fma_(
                            &ap[k],
                            x[i],
                            temp1,
                            numeric.fma(
                                y[i],
                                temp2,
                                ap[k],
                            ),
                        );

                        k += 1;
                    }

                    // ap[kk + j] += x[j] * temp1 + y[j] * temp2
                    numeric.fma_(
                        &ap[kk + j],
                        x[j],
                        temp1,
                        numeric.fma(
                            y[j],
                            temp2,
                            ap[kk + j],
                        ),
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
                    // temp1 = alpha * y[jy]
                    const temp1 = numeric.mul(
                        alpha,
                        y[numeric.cast(usize, jy)],
                    );

                    // temp2 = alpha * x[jx]
                    const temp2 = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    var ix: isize = kx;
                    var iy: isize = ky;
                    var k: usize = kk;
                    while (k < kk + j) : (k += 1) {
                        // ap[k] += x[ix] * temp1 + y[iy] * temp2
                        numeric.fma_(
                            &ap[k],
                            x[numeric.cast(usize, ix)],
                            temp1,
                            numeric.fma(
                                y[numeric.cast(usize, iy)],
                                temp2,
                                ap[k],
                            ),
                        );

                        ix += incx;
                        iy += incy;
                    }

                    // ap[kk + j] += x[jx] * temp1 + y[jy] * temp2
                    numeric.fma_(
                        &ap[kk + j],
                        x[numeric.cast(usize, jx)],
                        temp1,
                        numeric.fma(
                            y[numeric.cast(usize, jy)],
                            temp2,
                            ap[kk + j],
                        ),
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
                    // temp1 = alpha * y[j]
                    const temp1 = numeric.mul(
                        alpha,
                        y[j],
                    );

                    // temp2 = alpha * x[j]
                    const temp2 = numeric.mul(
                        alpha,
                        x[j],
                    );

                    // ap[kk] += x[j] * temp1 + y[j] * temp2
                    numeric.fma_(
                        &ap[kk],
                        x[j],
                        temp1,
                        numeric.fma(
                            y[j],
                            temp2,
                            ap[kk],
                        ),
                    );

                    var k: usize = kk + 1;
                    var i: usize = j + 1;
                    while (i < n) : (i += 1) {
                        // ap[k] += x[i] * temp1 + y[i] * temp2
                        numeric.fma_(
                            &ap[k],
                            x[i],
                            temp1,
                            numeric.fma(
                                y[i],
                                temp2,
                                ap[k],
                            ),
                        );

                        k += 1;
                    }
                }

                kk += n - j;
            }
        } else {
            var jx: isize = kx;
            var jy: isize = ky;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0) or numeric.ne(y[numeric.cast(usize, jy)], 0)) {
                    // temp1 = alpha * y[jy]
                    const temp1 = numeric.mul(
                        alpha,
                        y[numeric.cast(usize, jy)],
                    );

                    // temp2 = alpha * x[jx]
                    const temp2 = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    // ap[kk] += x[jx] * temp1 + y[jy] * temp2
                    numeric.fma_(
                        &ap[kk],
                        x[numeric.cast(usize, jx)],
                        temp1,
                        numeric.fma(
                            y[numeric.cast(usize, jy)],
                            temp2,
                            ap[kk],
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
                            x[numeric.cast(usize, ix)],
                            temp1,
                            numeric.fma(
                                y[numeric.cast(usize, iy)],
                                temp2,
                                ap[k],
                            ),
                        );
                    }
                }

                jx += incx;
                jy += incy;
                kk += n - j;
            }
        }
    }

    return;
}
