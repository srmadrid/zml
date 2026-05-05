const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;
const Diag = meta.Diag;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Solves a system of linear equations whose coefficients are in a triangular
/// band matrix defined as:
///
/// ```zig
/// A * x = b,
/// ```
///
/// or
///
/// ```zig
/// conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
/// Aᵀ * x = b,
/// ```
///
/// or
///
/// ```zig
/// Aᴴ * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, and `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular band matrix, with `k + 1` diagonals. The
/// function does not test for singularity or near-singularity, such tests must
/// be performed before calling this function.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.tbsv(layout: Layout, uplo: Uplo, transa: Transpose, diag: Diag, n: usize, k: usize, a: [*]const A, lda: usize, x: [*]X, incx: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
///   triangular matrix.
/// * `transa` (`linalg.Transpose`): Specifies the system of equations to be
///   solved:
///   * `no_transpose`: `A * x = b`
///   * `transpose`: `Aᵀ * x = b`
///   * `conj_no_transpose`: `conj(A) * x = b`
///   * `conj_transpose`: `Aᴴ * x = b`
/// * `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `k` (`usize`): Specifies the number of super-diagonals or sub-diagonals of
///   the matrix `A`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `k + 1`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`. Before entry, the incremented array `x` must
///   contain the n-element right-hand side vector `b`. On return, it contains
///   the solution vector `x`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `k + 1`, or if
///   `incx` is 0.
pub fn tbsv(layout: Layout, uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, k: usize, a: anytype, lda: usize, x: anytype, incx: isize) !void {
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.tbsv: a must be a many-item pointer to numerics, and x must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);

    if (lda < (k + 1) or incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and A == X) {
        switch (comptime meta.numericType(A)) {
            .float => {
                if (comptime A == f32)
                    return linalg.cblas.stbsv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), numeric.cast(isize, k), a, numeric.cast(isize, lda), x, incx)
                else if (comptime A == f64)
                    return linalg.cblas.dtbsv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), numeric.cast(isize, k), a, numeric.cast(isize, lda), x, incx);
            },
            .complex => {
                if (comptime meta.Scalar(A) == f32)
                    return linalg.cblas.ctbsv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), numeric.cast(isize, k), a, numeric.cast(isize, lda), x, incx)
                else if (comptime meta.Scalar(A) == f64)
                    return linalg.cblas.ztbsv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), numeric.cast(isize, k), a, numeric.cast(isize, lda), x, incx);
            },
            else => {},
        }
    }

    if (layout == .col_major) {
        return if (transa == .no_trans or transa == .trans)
            k_tbsv(uplo, transa, diag, n, k, a, lda, x, incx, true)
        else
            k_tbsv(uplo, transa, diag, n, k, a, lda, x, incx, false);
    } else {
        return if (transa == .no_trans or transa == .trans)
            k_tbsv(uplo.invert(), transa.invert(), diag, n, k, a, lda, x, incx, true)
        else
            k_tbsv(uplo.invert(), transa.invert(), diag, n, k, a, lda, x, incx, false);
    }
}

fn k_tbsv(uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, k: usize, a: anytype, lda: usize, x: anytype, incx: isize, comptime noconj: bool) void {
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (n == 0)
        return;

    const nounit: bool = diag == .non_unit;

    var kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        if (uplo == .upper) {
            if (incx == 1) {
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    if (numeric.ne(x[j], 0)) {
                        if (nounit) {
                            // x[j] /= a[k + j * lda]
                            numeric.div_(
                                &x[j],
                                x[j],
                                if (comptime noconj)
                                    a[k + j * lda]
                                else
                                    numeric.conj(a[k + j * lda]),
                            );
                        }

                        const temp = x[j];

                        var i: usize = j;
                        while (i > j -| k) {
                            i -= 1;

                            // x[i] -= temp * a[((k + i) - j) + j * lda]
                            numeric.fma_(
                                &x[i],
                                numeric.neg(temp),
                                if (comptime noconj)
                                    a[((k + i) - j) + j * lda]
                                else
                                    numeric.conj(a[((k + i) - j) + j * lda]),
                                x[i],
                            );
                        }
                    }
                }
            } else {
                kx += (n - 1) * incx;

                var jx: isize = kx;
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    kx -= incx;

                    if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                        if (nounit) {
                            // x[jx] /= a[k + j * lda]
                            numeric.div_(
                                &x[numeric.cast(usize, jx)],
                                x[numeric.cast(usize, jx)],
                                if (comptime noconj)
                                    a[k + j * lda]
                                else
                                    numeric.conj(a[k + j * lda]),
                            );
                        }

                        const temp = x[numeric.cast(usize, jx)];

                        var ix: isize = kx;
                        var i: usize = j;
                        while (i > j -| k) {
                            i -= 1;

                            // x[ix] -= temp * a[((k + i) - j) + j * lda]
                            numeric.fma_(
                                &x[numeric.cast(usize, ix)],
                                numeric.neg(temp),
                                if (comptime noconj)
                                    a[((k + i) - j) + j * lda]
                                else
                                    numeric.conj(a[((k + i) - j) + j * lda]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix -= incx;
                        }
                    }

                    jx -= incx;
                }
            }
        } else {
            if (incx == 1) {
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    if (numeric.ne(x[j], 0)) {
                        if (nounit) {
                            // x[j] /= a[0 + j * lda]
                            numeric.div_(
                                &x[j],
                                x[j],
                                if (comptime noconj)
                                    a[0 + j * lda]
                                else
                                    numeric.conj(a[0 + j * lda]),
                            );
                        }

                        const temp = x[j];

                        var i: usize = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            // x[i] -= temp * a[(i - j) + j * lda]
                            numeric.fma_(
                                &x[i],
                                numeric.neg(temp),
                                if (comptime noconj)
                                    a[(i - j) + j * lda]
                                else
                                    numeric.conj(a[(i - j) + j * lda]),
                                x[i],
                            );
                        }
                    }
                }
            } else {
                var jx: isize = kx;
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    kx += incx;

                    if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                        var ix: isize = kx;

                        if (nounit) {
                            // x[jx] /= a[0 + j * lda]
                            numeric.div_(
                                &x[numeric.cast(usize, jx)],
                                x[numeric.cast(usize, jx)],
                                if (comptime noconj)
                                    a[0 + j * lda]
                                else
                                    numeric.conj(a[0 + j * lda]),
                            );
                        }

                        const temp = x[numeric.cast(usize, jx)];

                        var i: usize = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            // x[ix] -= temp * a[(i - j) + j * lda]
                            numeric.fma_(
                                &x[numeric.cast(usize, ix)],
                                numeric.neg(temp),
                                if (comptime noconj)
                                    a[(i - j) + j * lda]
                                else
                                    numeric.conj(a[(i - j) + j * lda]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix += incx;
                        }
                    }

                    jx += incx;
                }
            }
        }
    } else {
        if (uplo == .upper) {
            if (incx == 1) {
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));
                    numeric.set(&temp, x[j]);

                    var i: usize = j -| k;
                    while (i < j) : (i += 1) {
                        // temp -= a[k + i + j * (lda - 1)] * x[i]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[k + i + j * (lda - 1)]
                            else
                                numeric.conj(a[k + i + j * (lda - 1)]),
                            numeric.neg(x[i]),
                            temp,
                        );
                    }

                    if (nounit) {
                        // temp /= a[k + j * lda]
                        numeric.div_(
                            &temp,
                            temp,
                            if (comptime noconj)
                                a[k + j * lda]
                            else
                                numeric.conj(a[k + j * lda]),
                        );
                    }

                    // x[j] = temp
                    numeric.set(
                        &x[j],
                        temp,
                    );
                }
            } else {
                var jx: isize = kx;
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));
                    numeric.set(&temp, x[numeric.cast(usize, jx)]);

                    var ix: isize = kx;
                    var i: usize = j -| k;
                    while (i < j) : (i += 1) {
                        // temp -= a[k + i + j * (lda - 1)] * x[ix]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[k + i + j * (lda - 1)]
                            else
                                numeric.conj(a[k + i + j * (lda - 1)]),
                            numeric.neg(x[numeric.cast(usize, ix)]),
                            temp,
                        );

                        ix += incx;
                    }

                    if (nounit) {
                        // temp /= a[k + j * lda]
                        numeric.div_(
                            &temp,
                            temp,
                            if (comptime noconj)
                                a[k + j * lda]
                            else
                                numeric.conj(a[k + j * lda]),
                        );
                    }

                    numeric.set(&x[numeric.cast(usize, jx)], temp);

                    jx += incx;

                    if (j >= k) {
                        kx += incx;
                    }
                }
            }
        } else {
            if (incx == 1) {
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));
                    numeric.set(&temp, x[j]);

                    var i: usize = int.min(n - 1, j + k);
                    while (i > j) {
                        i -= 1;

                        // temp -= a[(1 + i - j) + j * lda] * x[i + 1]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[(1 + i - j) + j * lda]
                            else
                                numeric.conj(a[(1 + i - j) + j * lda]),
                            numeric.neg(x[i + 1]),
                            temp,
                        );
                    }

                    if (nounit) {
                        // temp /= a[0 + j * lda]
                        numeric.div_(
                            &temp,
                            temp,
                            if (comptime noconj)
                                a[0 + j * lda]
                            else
                                numeric.conj(a[0 + j * lda]),
                        );
                    }

                    numeric.set(&x[j], temp);
                }
            } else {
                kx += (n - 1) * incx;
                var jx: isize = kx;
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));
                    numeric.set(&temp, x[numeric.cast(usize, jx)]);

                    var ix: isize = kx;
                    var i: usize = int.min(n - 1, j + k);
                    while (i > j) {
                        i -= 1;

                        // temp -= a[(1 + i - j) + j * lda] * x[ix]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[(1 + i - j) + j * lda]
                            else
                                numeric.conj(a[(1 + i - j) + j * lda]),
                            numeric.neg(x[numeric.cast(usize, ix)]),
                            temp,
                        );

                        ix -= incx;
                    }

                    if (nounit) {
                        // temp /= a[0 + j * lda]
                        numeric.div_(
                            &temp,
                            temp,
                            if (comptime noconj)
                                a[0 + j * lda]
                            else
                                numeric.conj(a[0 + j * lda]),
                        );
                    }

                    numeric.set(&x[numeric.cast(usize, jx)], temp);

                    jx -= incx;

                    if (n - 1 - j >= k) {
                        kx -= incx;
                    }
                }
            }
        }
    }

    return;
}
