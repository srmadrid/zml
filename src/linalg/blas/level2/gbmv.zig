const options = @import("options");

const meta = @import("../../../meta.zig");
const Layout = meta.Layout;

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

/// Computes a matrix-vector product with a general band matrix defined as:
///
/// ```zig
/// y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * Aᵀ * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * Aᴴ * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are numerics, `x` and `y` are vectors, and `A` is an
/// `m`-by-`n` band matrix with `kl` sub-diagonals and `ku` super-diagonals.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.gbmv(layout: Layout, transa: linalg.Transpose, m: usize, n: usize, kl: usize, ku: usize, alpha: Al, a: [*]const A, lda: usize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `transa` (`linalg.Transpose`): Specifies the operation to be performed on
///   `A`:
///   * `no_transpose`: `y = alpha * A * x + beta * y`
///   * `transpose`: `y = alpha * Aᵀ * x + beta * y`
///   * `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
///   * `conj_transpose`: `y = alpha * Aᴴ * x + beta * y`
/// * `m` (`usize`): Specifies the number of rows of the matrix `A`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `kl` (`usize`): Specifies the number of sub-diagonals of the matrix `A`.
/// * `ku` (`usize`): Specifies the number of super-diagonals of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `kl + ku + 1`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)` when `transa` is `no_transpose` or
///   `conj_no_transpose`, or `1 + (m - 1) * abs(incx)` otherwise.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `y` need not be set on input.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (m - 1) * abs(incy)` when `transa` is `no_transpose` or
///   `conj_no_transpose`, or `1 + (n - 1) * abs(incy)` otherwise. On return,
///   contains the result of the operation.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `kl + ku + 1`,
///   or if `incx` or `incy` is 0.
pub fn gbmv(layout: Layout, transa: linalg.Transpose, m: usize, n: usize, kl: usize, ku: usize, alpha: anytype, a: anytype, lda: usize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.gbmv: alpha and beta must be numerics, a and x must be many-item pointers to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);
    Y = meta.Child(Y);

    if (lda < (kl + ku + 1) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Al == A and Al == X and Al == Y and Al == Be) {
        switch (comptime meta.numericType(Al)) {
            .float => {
                if (comptime Al == f32)
                    return linalg.cblas.sgbmv(layout.toInt(c_int), transa.toInt(c_int), numeric.cast(isize, m), numeric.cast(isize, n), numeric.cast(isize, kl), numeric.cast(isize, ku), alpha, a, numeric.cast(isize, lda), x, incx, beta, y, incy)
                else if (comptime Al == f64)
                    return linalg.cblas.dgbmv(layout.toInt(c_int), transa.toInt(c_int), numeric.cast(isize, m), numeric.cast(isize, n), numeric.cast(isize, kl), numeric.cast(isize, ku), alpha, a, numeric.cast(isize, lda), x, incx, beta, y, incy);
            },
            .complex => {
                if (comptime meta.Scalar(Al) == f32)
                    return linalg.cblas.cgbmv(layout.toInt(c_int), transa.toInt(c_int), numeric.cast(isize, m), numeric.cast(isize, n), numeric.cast(isize, kl), numeric.cast(isize, ku), &alpha, a, numeric.cast(isize, lda), x, incx, &beta, y, incy)
                else if (comptime meta.Scalar(Al) == f64)
                    return linalg.cblas.zgbmv(layout.toInt(c_int), transa.toInt(c_int), numeric.cast(isize, m), numeric.cast(isize, n), numeric.cast(isize, kl), numeric.cast(isize, ku), &alpha, a, numeric.cast(isize, lda), x, incx, &beta, y, incy);
            },
            else => {},
        }
    }

    if (layout == .col_major) {
        return if (transa == .no_trans or transa == .trans)
            k_gbmv(transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_gbmv(transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, false);
    } else {
        return if (transa == .no_trans or transa == .trans)
            k_gbmv(transa.invert(), n, m, ku, kl, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_gbmv(transa.invert(), n, m, ku, kl, alpha, a, lda, x, incx, beta, y, incy, false);
    }
}

fn k_gbmv(transa: linalg.Transpose, m: usize, n: usize, kl: usize, ku: usize, alpha: anytype, a: anytype, lda: usize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize, comptime noconj: bool) void {
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (m == 0 or n == 0 or (numeric.eq(alpha, 0) and numeric.eq(beta, 1)))
        return;

    // Set lenx and leny, the lengths of the vectors x and y, and set up the
    // start points in x and y.
    var lenx: usize = undefined;
    var leny: usize = undefined;
    if (transa == .no_trans or transa == .conj_no_trans) {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    var kx: isize = if (incx < 0) (-numeric.cast(isize, lenx) + 1) * incx else 0;
    var ky: isize = if (incy < 0) (-numeric.cast(isize, leny) + 1) * incy else 0;

    // First form y = beta * y.
    if (numeric.ne(beta, 1))
        @import("../level1/scal.zig").k_scal(leny, beta, y, incy);

    if (numeric.eq(alpha, 0))
        return;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Form  y = alpha * A * x + y  or  y = alpha * conj(A) * x + y.
        var jx: isize = kx;
        var j: usize = 0;
        while (j < n) : (j += 1) {
            // temp = alpha * x[jx]
            const temp = numeric.mul(
                alpha,
                x[numeric.cast(usize, jx)],
            );

            var i: usize = j -| ku;
            if (incy == 1) {
                while (i < int.min(m, j + kl + 1)) : (i += 1) {
                    // y[i] += temp * a[i + ku + j * (lda - 1)]
                    numeric.fma_(
                        &y[i],
                        temp,
                        if (comptime noconj)
                            a[i + ku + j * (lda - 1)]
                        else
                            numeric.conj(a[i + ku + j * (lda - 1)]),
                        y[i],
                    );
                }
            } else {
                var iy: isize = ky;
                while (i < int.min(m, j + kl + 1)) : (i += 1) {
                    // y[iy] += temp * a[i + ku + j * (lda - 1)]
                    numeric.fma_(
                        &y[numeric.cast(usize, iy)],
                        temp,
                        if (comptime noconj)
                            a[i + ku + j * (lda - 1)]
                        else
                            numeric.conj(a[i + ku + j * (lda - 1)]),
                        y[numeric.cast(usize, iy)],
                    );

                    iy += incy;
                }
            }

            jx += incx;

            if (j >= ku)
                ky += incy;
        }
    } else {
        // Form  y = alpha * Aᵀ * x + y  or  y = alpha * Aᴴ * x + y.
        var jy: isize = ky;
        var j: usize = 0;
        while (j < n) : (j += 1) {
            var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

            var i: usize = j -| ku;
            if (incx == 1) {
                while (i < int.min(m, j + kl + 1)) : (i += 1) {
                    // temp += a[i + ku + j * (lda - 1)] * x[i]
                    numeric.fma_(
                        &temp,
                        if (comptime noconj)
                            a[i + ku + j * (lda - 1)]
                        else
                            numeric.conj(a[i + ku + j * (lda - 1)]),
                        x[i],
                        temp,
                    );
                }
            } else {
                var ix: isize = kx;
                while (i < int.min(m, j + kl + 1)) : (i += 1) {
                    // temp += a[i + ku + j * (lda - 1)] * x[ix]
                    numeric.fma_(
                        &temp,
                        if (comptime noconj)
                            a[i + ku + j * (lda - 1)]
                        else
                            numeric.conj(a[i + ku + j * (lda - 1)]),
                        x[numeric.cast(usize, ix)],
                        temp,
                    );

                    ix += incx;
                }
            }

            // y[jy] += alpha * temp
            numeric.fma_(
                &y[numeric.cast(usize, jy)],
                alpha,
                temp,
                y[numeric.cast(usize, jy)],
            );

            jy += incy;

            if (j >= ku)
                kx += incx;
        }
    }

    return;
}
