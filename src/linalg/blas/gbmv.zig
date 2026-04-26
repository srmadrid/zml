const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

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
/// linalg.blas.gbmv(layout: Layout, transa: linalg.Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: Al, a: [*]const A, lda: isize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
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
/// * `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
///   greater than or equal to 0.
/// * `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
///   greater than or equal to 0.
/// * `kl` (`isize`): Specifies the number of sub-diagonals of the matrix `A`.
///   Must be greater than or equal to 0.
/// * `ku` (`isize`): Specifies the number of super-diagonals of the matrix `A`.
///   Must be greater than or equal to 0.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Array, size at least `lda * n`.
/// * `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `kl + ku + 1`.
/// * `x` (`anytype`): Array, size at least `1 + (n - 1) * abs(incx)` when
///   `transa` is `no_transpose` or `conj_no_transpose`, or
///   `1 + (m - 1) * abs(incx)` otherwise.
/// * `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
///   different from 0.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `y` need not be set on input.
/// * `y` (`anytype`): Array, size at least `1 + (m - 1) * abs(incy)` when
///   `transa` is `no_transpose` or `conj_no_transpose`, or
///   `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
///   operation.
/// * `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
///   different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `m`, `n`, `kl` or `ku` are less
///   than 0, if `lda` is less than `kl + ku + 1`, or if `incx` or `incy` is 0.
pub fn gbmv(layout: Layout, transa: linalg.Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: anytype, a: anytype, lda: isize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize) !void {
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

    if (m < 0 or n < 0 or kl < 0 or ku < 0 or lda < (kl + ku + 1) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Al == A and Al == X and Al == Y and Al == Be) {
        switch (comptime meta.numericType(A)) {
            .float => {
                if (comptime A == f32)
                    return linalg.cblas.sgbmv(layout.toInt(c_int), transa.toInt(c_int), m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
                else if (comptime A == f64)
                    return linalg.cblas.dgbmv(layout.toInt(c_int), transa.toInt(c_int), m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
            },
            .complex => {
                if (comptime meta.Scalar(A) == f32)
                    return linalg.cblas.cgbmv(layout.toInt(c_int), transa.toInt(c_int), m, n, kl, ku, &alpha, a, lda, x, incx, &beta, y, incy)
                else if (comptime meta.Scalar(A) == f64)
                    return linalg.cblas.zgbmv(layout.toInt(c_int), transa.toInt(c_int), m, n, kl, ku, &alpha, a, lda, x, incx, &beta, y, incy);
            },
            else => {},
        }
    }

    if (layout == .col_major) {
        return k_gbmv(transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        return k_gbmv(transa.invert(), n, m, ku, kl, alpha, a, lda, x, incx, beta, y, incy);
    }
}

fn k_gbmv(transa: linalg.Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: anytype, a: anytype, lda: isize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize) !void {
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (m == 0 or n == 0 or (numeric.eq(alpha, 0) and numeric.eq(beta, 1)))
        return;

    const noconj: bool = transa == .no_trans or transa == .trans;

    // Set lenx and leny, the lengths of the vectors x and y, and set up the
    // start points in x and y.
    var lenx: isize = undefined;
    var leny: isize = undefined;
    if (transa == .no_trans or transa == .conj_no_trans) {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    var kx: isize = if (incx < 0) (-lenx + 1) * incx else 0;
    var ky: isize = if (incy < 0) (-leny + 1) * incy else 0;

    // First form y = beta * y.
    if (numeric.ne(beta, 1)) {
        if (incy == 1) {
            if (numeric.eq(beta, 0)) {
                var i: usize = 0;
                while (i < leny) : (i += 1) {
                    // y[i] = 0
                    numeric.set(&y[i], 0);
                }
            } else {
                var i: usize = 0;
                while (i < leny) : (i += 1) {
                    // y[i] *= beta
                    numeric.mul_(&y[i], y[i], beta);
                }
            }
        } else {
            var iy: isize = ky;
            if (numeric.eq(beta, 0)) {
                var i: usize = 0;
                while (i < leny) : (i += 1) {
                    // y[iy] = 0
                    numeric.set(&y[numeric.cast(usize, iy)], 0);

                    iy += incy;
                }
            } else {
                var i: usize = 0;
                while (i < leny) : (i += 1) {
                    // y[iy] *= beta
                    numeric.mul_(&y[numeric.cast(usize, iy)], y[numeric.cast(usize, iy)], beta);

                    iy += incy;
                }
            }
        }
    }

    if (numeric.eq(alpha, 0))
        return;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Form  y = alpha * A * x + y  or  y = alpha * conj(A) * x + y.
        var jx: isize = kx;
        if (incy == 1) {
            var j: isize = 0;
            while (j < n) : (j += 1) {
                // temp = alpha * x[jx]
                const temp = numeric.mul(
                    alpha,
                    x[numeric.cast(usize, jx)],
                );

                const k: isize = ku - j;
                if (noconj) {
                    var i: isize = int.max(0, j - ku);
                    while (i < int.min(m, j + kl + 1)) : (i += 1) {
                        // y[i] += temp * a[k + i + j * lda]
                        numeric.add_(
                            &y[numeric.cast(usize, i)],
                            y[numeric.cast(usize, i)],
                            numeric.mul(
                                temp,
                                a[numeric.cast(usize, k + i + j * lda)],
                            ),
                        );
                    }
                } else {
                    var i: isize = int.max(0, j - ku);
                    while (i < int.min(m, j + kl + 1)) : (i += 1) {
                        // y[i] += temp * conj(a[k + i + j * lda])
                        numeric.add_(
                            &y[numeric.cast(usize, i)],
                            y[numeric.cast(usize, i)],
                            numeric.mul(
                                temp,
                                numeric.conj(a[numeric.cast(usize, k + i + j * lda)]),
                            ),
                        );
                    }
                }

                jx += incx;
            }
        } else {
            var j: isize = 0;
            while (j < n) : (j += 1) {
                // temp = alpha * x[jx]
                const temp = numeric.mul(
                    alpha,
                    x[numeric.cast(usize, jx)],
                );

                const k: isize = ku - j;
                if (noconj) {
                    var iy: isize = ky;
                    var i: isize = int.max(0, j - ku);
                    while (i < int.min(m, j + kl + 1)) : (i += 1) {
                        // y[iy] += temp * a[k + i + j * lda]
                        numeric.add_(
                            &y[numeric.cast(usize, iy)],
                            y[numeric.cast(usize, iy)],
                            numeric.mul(
                                temp,
                                a[numeric.cast(usize, k + i + j * lda)],
                            ),
                        );
                        iy += incy;
                    }
                } else {
                    var iy: isize = ky;
                    var i: isize = int.max(0, j - ku);
                    while (i < int.min(m, j + kl + 1)) : (i += 1) {
                        // y[iy] += temp * a[k + i + j * lda]
                        numeric.add_(
                            &y[numeric.cast(usize, iy)],
                            y[numeric.cast(usize, iy)],
                            numeric.mul(
                                temp,
                                numeric.conj(a[numeric.cast(usize, k + i + j * lda)]),
                            ),
                        );

                        iy += incy;
                    }
                }

                jx += incx;

                if (j >= ku)
                    ky += incy;
            }
        }
    } else {
        // Form  y = alpha * Aᵀ * x + y  or  y = alpha * Aᴴ * x + y.
        var jy: isize = ky;
        if (incx == 1) {
            var j: isize = 0;
            while (j < n) : (j += 1) {
                var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                const k: isize = ku - j;
                if (noconj) {
                    var i: isize = int.max(0, j - ku);
                    while (i < int.min(m, j + kl + 1)) : (i += 1) {
                        // temp += a[k + i + j * lda] * x[i]
                        numeric.add_(
                            &temp,
                            temp,
                            numeric.mul(
                                a[numeric.cast(usize, k + i + j * lda)],
                                x[numeric.cast(usize, i)],
                            ),
                        );
                    }
                } else {
                    var i: isize = int.max(0, j - ku);
                    while (i < int.min(m, j + kl + 1)) : (i += 1) {
                        // temp += conj(a[k + i + j * lda]) * x[i]
                        numeric.add_(
                            &temp,
                            temp,
                            numeric.mul(
                                numeric.conj(a[numeric.cast(usize, k + i + j * lda)]),
                                x[numeric.cast(usize, i)],
                            ),
                        );
                    }
                }

                // y[jy] += alpha * temp
                numeric.add_(
                    &y[numeric.cast(usize, jy)],
                    y[numeric.cast(usize, jy)],
                    numeric.mul(
                        alpha,
                        temp,
                    ),
                );

                jy += incy;
            }
        } else {
            var j: isize = 0;
            while (j < n) : (j += 1) {
                var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                const k: isize = ku - j;
                if (noconj) {
                    var ix: isize = kx;
                    var i: isize = int.max(0, j - ku);
                    while (i < int.min(m, j + kl + 1)) : (i += 1) {
                        // temp += a[k + i + j * lda] * x[ix]
                        numeric.add_(
                            &temp,
                            temp,
                            numeric.mul(
                                a[numeric.cast(usize, k + i + j * lda)],
                                x[numeric.cast(usize, ix)],
                            ),
                        );

                        ix += incx;
                    }
                } else {
                    var ix: isize = kx;
                    var i: isize = int.max(0, j - ku);
                    while (i < int.min(m, j + kl + 1)) : (i += 1) {
                        // temp += conj(a[k + i + j * lda]) * x[ix]
                        numeric.add_(
                            &temp,
                            temp,
                            numeric.mul(
                                numeric.conj(a[numeric.cast(usize, k + i + j * lda)]),
                                x[numeric.cast(usize, ix)],
                            ),
                        );

                        ix += incx;
                    }
                }

                // y[jy] += alpha * temp
                numeric.add_(
                    &y[numeric.cast(usize, jy)],
                    y[numeric.cast(usize, jy)],
                    numeric.mul(
                        alpha,
                        temp,
                    ),
                );

                jy += incy;

                if (j >= ku)
                    kx += incx;
            }
        }
    }

    return;
}
