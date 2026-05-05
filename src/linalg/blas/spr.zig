const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Performs a rank-1 update of a symmetric packed matrix defined as:
///
/// ```zig
/// A = alpha * x * xᵀ + A,
/// ```
///
/// where `alpha` is a real numeric, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.spr(layout: Layout, uplo: Uplo, n: usize, alpha: Al, x: [*]const X, incx: isize, ap: [*]Ap) !void
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
/// * `ap` (`anytype`): Mutable many-item pointer, size at least
///   `(n * (n + 1)) / 2`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is 0.
pub fn spr(layout: Layout, uplo: Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, ap: anytype) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Ap: type = @TypeOf(ap);

    comptime if (!meta.isNumeric(Al) or !meta.isReal(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Ap) or meta.isConstPointer(Ap) or !meta.isNumeric(meta.Child(Ap)))
        @compileError("zsl.linalg.blas.spr: alpha must be a numeric, x must be a many-item pointer to numerics, and ap must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tap: " ++ @typeName(Ap) ++ "\n");

    X = meta.Child(X);
    Ap = meta.Child(Ap);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Al == X and Al == Ap) {
        switch (comptime meta.numericType(Al)) {
            .float => {
                if (comptime Al == f32)
                    return linalg.cblas.sspr(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), alpha, x, incx, ap)
                else if (comptime Al == f64)
                    return linalg.cblas.dspr(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), alpha, x, incx, ap);
            },
            else => {},
        }
    }

    return if (layout == .col_major)
        k_spr(uplo, n, alpha, x, incx, ap)
    else
        k_spr(uplo.invert(), n, alpha, x, incx, ap);
}

fn k_spr(uplo: Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, ap: anytype) void {
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
                    // temp = alpha * x[j]
                    const temp = numeric.mul(
                        alpha,
                        x[j],
                    );

                    var k: usize = kk;
                    var i: usize = 0;
                    while (i < j) : (i += 1) {
                        // ap[k] += x[i] * temp
                        numeric.fma(
                            &ap[k],
                            x[i],
                            temp,
                            ap[k],
                        );

                        k += 1;
                    }

                    // ap[kk + j] += x[j] * temp
                    numeric.fma_(
                        &ap[kk + j],
                        x[j],
                        temp,
                        ap[kk + j],
                    );
                }

                kk += j + 1;
            }
        } else {
            var jx: isize = kx;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                    // temp = alpha * x[jx]
                    const temp = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    var ix: isize = kx;
                    var k: usize = kk;
                    while (k < kk + j) : (k += 1) {
                        // ap[k] += x[ix] * temp
                        numeric.fma(
                            &ap[k],
                            x[numeric.cast(usize, ix)],
                            temp,
                            ap[k],
                        );

                        ix += incx;
                    }

                    // ap[kk + j] += x[jx] * temp
                    numeric.fma_(
                        &ap[kk + j],
                        x[numeric.cast(usize, jx)],
                        temp,
                        ap[kk + j],
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
                    // temp = alpha * x[j]
                    const temp = numeric.mul(
                        alpha,
                        x[j],
                    );

                    // ap[kk] += x[j] * temp
                    numeric.fma_(
                        &ap[kk],
                        x[j],
                        temp,
                        ap[kk],
                    );

                    var k: usize = kk + 1;
                    var i: usize = j + 1;
                    while (i < n) : (i += 1) {
                        // ap[k] += x[i] * temp
                        numeric.fma_(
                            &ap[k],
                            x[i],
                            temp,
                            ap[k],
                        );

                        k += 1;
                    }
                }

                kk += n - j;
            }
        } else {
            var jx: isize = kx;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                    // temp = alpha * x[jx]
                    const temp = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    // ap[kk] += x[jx] * temp
                    numeric.fma_(
                        &ap[kk],
                        x[numeric.cast(usize, jx)],
                        temp,
                        ap[kk],
                    );

                    var ix: isize = jx;
                    var k: usize = kk + 1;
                    while (k < kk + n - j) : (k += 1) {
                        ix += incx;

                        // ap[k] += x[ix] * temp
                        numeric.fma_(
                            &ap[k],
                            x[numeric.cast(usize, ix)],
                            temp,
                            ap[k],
                        );
                    }
                }

                jx += incx;
                kk += n - j;
            }
        }
    }

    return;
}
