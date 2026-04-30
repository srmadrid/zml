const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Performs a rank-1 update of a Hermitian packed matrix defined as:
///
/// ```zig
/// A = alpha * x * xᴴ + A,
/// ```
///
/// where `alpha` is a real numeric, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.hpr(layout: Layout, uplo: Uplo, n: usize, alpha: Al, x: [*]const X, incx: isize, ap: [*]Ap) !void
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
/// * `ap` (`anytype`): Mutable many-item pointer, size at least
///   `(n * (n + 1)) / 2`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is 0.
pub fn hpr(layout: Layout, uplo: Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, ap: anytype) !void {
    const Al: type = @TypeOf(alpha);
    comptime var Ap: type = @TypeOf(ap);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(Al) or !meta.isReal(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Ap) or meta.isConstPointer(Ap) or !meta.isNumeric(meta.Child(Ap)))
        @compileError("zsl.linalg.blas.hpr: alpha must be a numeric, x must be a many-item pointer to numerics, and ap must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tap: " ++ @typeName(Ap) ++ "\n");

    Ap = meta.Child(Ap);
    X = meta.Child(X);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and meta.Real(X) == Al and X == Ap) {
        switch (comptime meta.numericType(X)) {
            .complex => {
                if (comptime meta.Scalar(X) == f32)
                    return linalg.cblas.chpr(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), alpha, x, incx, ap)
                else if (comptime meta.Scalar(X) == f64)
                    return linalg.cblas.zhpr(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), alpha, x, incx, ap);
            },
            else => {},
        }
    }

    return if (layout == .col_major)
        k_hpr(uplo, n, alpha, x, incx, ap, true)
    else
        k_hpr(uplo.invert(), n, alpha, x, incx, ap, false);
}

fn k_hpr(uplo: Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, ap: anytype, comptime noconj: bool) !void {
    // Quick return if possible.
    if (n == 0 or numeric.eq(alpha, 0))
        return;

    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    var kk: usize = 0;
    if (uplo == .upper) {
        // Form  A  when upper triangle of A is stored.
        if (incx == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[j], 0)) {
                    // temp = alpha * conj(x[j])
                    const temp = numeric.mul(
                        alpha,
                        if (comptime !noconj)
                            x[j]
                        else
                            numeric.conj(x[j]),
                    );

                    var k: usize = kk;
                    var i: usize = 0;
                    while (i < j) : (i += 1) {
                        // ap[k] += x[i] * temp
                        numeric.fma(
                            &ap[k],
                            if (comptime noconj)
                                x[i]
                            else
                                numeric.conj(x[i]),
                            temp,
                            ap[k],
                        );

                        k += 1;
                    }

                    // ap[kk + j] = re(ap[kk + j]) + re(x[j] * temp)
                    numeric.fma_(
                        &ap[kk + j],
                        if (comptime !noconj)
                            numeric.im(x[j])
                        else
                            numeric.neg(numeric.im(x[j])),
                        numeric.im(temp),
                        numeric.fma(
                            numeric.re(x[j]),
                            numeric.re(temp),
                            numeric.re(ap[kk + j]),
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
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                    // temp = alpha * conj(x[jx])
                    const temp = numeric.mul(
                        alpha,
                        if (comptime !noconj)
                            x[numeric.cast(usize, jx)]
                        else
                            numeric.conj(x[numeric.cast(usize, jx)]),
                    );

                    var ix: isize = kx;
                    var k: usize = kk;
                    while (k < kk + j) : (k += 1) {
                        // ap[k] += x[ix] * temp
                        numeric.fma(
                            &ap[k],
                            if (comptime noconj)
                                x[numeric.cast(usize, ix)]
                            else
                                numeric.conj(x[numeric.cast(usize, ix)]),
                            temp,
                            ap[k],
                        );

                        ix += incx;
                    }

                    // ap[kk + j] = re(ap[kk + j]) + re(x[jx] * temp)
                    numeric.fma_(
                        &ap[kk + j],
                        if (comptime !noconj)
                            numeric.im(x[numeric.cast(usize, jx)])
                        else
                            numeric.neg(numeric.im(x[numeric.cast(usize, jx)])),
                        numeric.im(temp),
                        numeric.fma(
                            numeric.re(x[numeric.cast(usize, jx)]),
                            numeric.re(temp),
                            numeric.re(ap[kk + j]),
                        ),
                    );
                } else {
                    numeric.set( // ap[kk + j] = re(ap[kk + j])
                        &ap[kk + j],
                        numeric.re(ap[kk + j]),
                    );
                }

                jx += incx;
                kk += j + 1;
            }
        }
    } else {
        if (incx == 1) {
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[j], 0)) {
                    // temp = alpha * conj(x[j])
                    const temp = numeric.mul(
                        alpha,
                        if (comptime !noconj)
                            x[j]
                        else
                            numeric.conj(x[j]),
                    );

                    // ap[kk] = re(ap[kk]) + re(x[j] * temp)
                    numeric.fma_(
                        &ap[kk],
                        if (comptime !noconj)
                            numeric.im(x[j])
                        else
                            numeric.neg(numeric.im(x[j])),
                        numeric.im(temp),
                        numeric.fma(
                            numeric.re(x[j]),
                            numeric.re(temp),
                            numeric.re(ap[kk]),
                        ),
                    );

                    var k: usize = kk + 1;
                    var i: usize = j + 1;
                    while (i < n) : (i += 1) {
                        // ap[k] += x[i] * temp
                        numeric.fma_(
                            &ap[k],
                            if (comptime noconj)
                                x[i]
                            else
                                numeric.conj(x[i]),
                            temp,
                            ap[k],
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
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                    // temp = alpha * conj(x[jx])
                    const temp = numeric.mul(
                        alpha,
                        if (comptime !noconj)
                            x[numeric.cast(usize, jx)]
                        else
                            numeric.conj(x[numeric.cast(usize, jx)]),
                    );

                    // ap[kk] = re(ap[kk]) + re(x[jx] * temp)
                    numeric.fma_(
                        &ap[kk],
                        if (comptime !noconj)
                            numeric.im(x[numeric.cast(usize, jx)])
                        else
                            numeric.neg(numeric.im(x[numeric.cast(usize, jx)])),
                        numeric.im(temp),
                        numeric.fma(
                            numeric.re(x[numeric.cast(usize, jx)]),
                            numeric.re(temp),
                            numeric.re(ap[kk]),
                        ),
                    );

                    var ix: isize = jx;
                    var k: usize = kk + 1;
                    while (k < kk + n - j) : (k += 1) {
                        ix += incx;

                        // ap[k] += x[ix] * temp
                        numeric.fma_(
                            &ap[k],
                            if (comptime noconj)
                                x[numeric.cast(usize, ix)]
                            else
                                numeric.conj(x[numeric.cast(usize, ix)]),
                            temp,
                            ap[k],
                        );
                    }
                } else {
                    numeric.set( // ap[kk] = re(ap[kk])
                        &ap[kk],
                        numeric.re(ap[kk]),
                    );
                }

                jx += incx;
                kk += n - j;
            }
        }
    }

    return;
}
