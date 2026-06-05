const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;
const Diag = meta.Diag;

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");
const float = @import("../../../float.zig");

const linalg = @import("../../../linalg.zig");

/// Computes a matrix-vector product using a triangular matrix defined as:
///
/// ```zig
/// x = A * x,
/// ```
///
/// or
///
/// ```zig
/// x = Aᵀ * x,
/// ```
///
/// or
///
/// ```zig
/// x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
/// x = Aᴴ * x,
/// ```
///
/// where `x` is an `n`-element vector, and `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.trmv(layout: Layout, uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, a: [*]const A, lda: usize, x: [*]X, incx: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
///   triangular matrix:
/// * `transa` (`linalg.Transpose`): Specifies the operation to be performed on
///   `A`:
///   * `no_transpose`: `y = A * x`
///   * `transpose`: `x = Aᵀ * x`
///   * `conj_no_transpose`: `x = conj(A) * x`
///   * `conj_transpose`: `x = Aᴴ * x`
/// * `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`. On return, contains the result of the
///   operation.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, n)`, or
///   if `incx` is 0.
pub fn trmv(
    layout: Layout,
    uplo: Uplo,
    transa: linalg.Transpose,
    diag: Diag,
    n: usize,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.trmv: a must be a many-item pointer to numerics, and x must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);

    if (lda < int.max(1, n) or incx == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();
    const eff_transa = if (layout == .col_major) transa else transa.invert();
    const noconj = transa == .no_trans or transa == .trans;
    // const no_trans = eff_transa == .no_trans or eff_transa == .conj_no_trans;

    // Quick return if possible.
    if (n == 0)
        return;

    if (noconj)
        return k_trmv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, true)
    else
        return k_trmv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, false);
}

fn k_trmv(uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, a: anytype, lda: usize, x: anytype, incx: isize, comptime noconj: bool) void {
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (n == 0)
        return;

    // Set up the start points in x.
    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Form  y = alpha * A * x + y  or  y = alpha * conj(A) * x + y.
        if (uplo == .upper) {
            if (incx == 1) {
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    if (numeric.ne(x[j], 0)) {
                        const temp = x[j];

                        var i: usize = 0;
                        while (i < j) : (i += 1) {
                            // x[i] += temp * a[i + j * lda]
                            numeric.fma_(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                            );
                        }

                        if (diag == .non_unit) {
                            // x[j] *= a[j + j * lda]
                            numeric.mul_(
                                &x[j],
                                x[j],
                                if (comptime noconj)
                                    a[j + j * lda]
                                else
                                    numeric.conj(a[j + j * lda]),
                            );
                        }
                    }
                }
            } else {
                var jx: isize = kx;
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                        const temp = x[numeric.cast(usize, jx)];

                        var ix: isize = kx;
                        var i: usize = 0;
                        while (i < j) : (i += 1) {
                            // x[ix] += temp * a[i + j * lda]
                            numeric.fma_(
                                &x[numeric.cast(usize, ix)],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix += incx;
                        }

                        if (diag == .non_unit) {
                            // x[jx] *= a[j + j * lda]
                            numeric.mul_(
                                &x[numeric.cast(usize, jx)],
                                x[numeric.cast(usize, jx)],
                                if (comptime noconj)
                                    a[j + j * lda]
                                else
                                    numeric.conj(a[j + j * lda]),
                            );
                        }
                    }

                    jx += incx;
                }
            }
        } else {
            if (incx == 1) {
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    if (numeric.ne(x[j], 0)) {
                        const temp = x[j];

                        var i: usize = n;
                        while (i > j + 1) {
                            i -= 1;

                            // x[i] += temp * a[i + j * lda]
                            numeric.fma_(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                            );
                        }

                        if (diag == .non_unit) {
                            numeric.mul_( // x[j] *= a[j + j * lda]
                                &x[j],
                                x[j],
                                if (comptime noconj)
                                    a[j + j * lda]
                                else
                                    numeric.conj(a[j + j * lda]),
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

                    if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                        const temp = x[numeric.cast(usize, jx)];

                        var ix: isize = kx;
                        var i: usize = n;
                        while (i > j + 1) {
                            i -= 1;

                            // x[ix] += temp * a[i + j * lda]
                            numeric.fma_(
                                &x[numeric.cast(usize, ix)],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix -= incx;
                        }

                        if (diag == .non_unit) {
                            numeric.mul_( // x[jx] *= a[j + j * lda]
                                &x[numeric.cast(usize, jx)],
                                x[numeric.cast(usize, jx)],
                                if (comptime noconj)
                                    a[j + j * lda]
                                else
                                    numeric.conj(a[j + j * lda]),
                            );
                        }
                    }

                    jx -= incx;
                }
            }
        }
    } else {
        if (uplo == .upper) {
            if (incx == 1) {
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                    if (diag == .non_unit) {
                        // temp += a[j + j * lda] * x[j]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[j],
                            temp,
                        );
                    } else {
                        // temp += x[j]
                        numeric.add_(
                            &temp,
                            temp,
                            x[j],
                        );
                    }

                    var i: usize = j;
                    while (i > 0) {
                        i -= 1;

                        // temp += a[i + j * lda] * x[i]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[i],
                            temp,
                        );
                    }

                    // x[j] = temp
                    numeric.set(&x[j], temp);
                }
            } else {
                var jx: isize = kx + (n - 1) * incx;
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                    if (diag == .non_unit) {
                        // temp += a[j + j * lda] * x[jx]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[numeric.cast(usize, jx)],
                            temp,
                        );
                    } else {
                        // temp += x[jx]
                        numeric.add_(
                            &temp,
                            temp,
                            x[numeric.cast(usize, jx)],
                        );
                    }

                    var ix: isize = jx;
                    var i: usize = j;
                    while (i > 0) {
                        ix -= incx;
                        i -= 1;

                        // temp += a[i + j * lda] * x[ix]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[numeric.cast(usize, ix)],
                            temp,
                        );
                    }

                    // x[jx] = temp
                    numeric.set(&x[numeric.cast(usize, jx)], temp);

                    jx -= incx;
                }
            }
        } else {
            if (incx == 1) {
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                    if (diag == .non_unit) {
                        // temp += a[j + j * lda] * x[j]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[j],
                            temp,
                        );
                    } else {
                        // temp += x[j]
                        numeric.add_(
                            &temp,
                            temp,
                            x[j],
                        );
                    }

                    var i: usize = j + 1;
                    while (i < n) : (i += 1) {
                        numeric.add_( // temp += a[i + j * lda] * x[i]
                            &temp,
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[i],
                            temp,
                        );
                    }

                    // x[j] = temp
                    numeric.set(&x[j], temp);
                }
            } else {
                var jx: isize = kx;
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                    if (diag == .non_unit) {
                        // temp += a[j + j * lda] * x[jx]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                            x[numeric.cast(usize, jx)],
                            temp,
                        );
                    } else {
                        // temp += x[jx]
                        numeric.add_(
                            &temp,
                            temp,
                            x[numeric.cast(usize, jx)],
                        );
                    }

                    var ix: usize = jx;
                    var i: usize = j + 1;
                    while (i < n) : (i += 1) {
                        ix += incx;

                        // temp += a[i + j * lda] * x[ix]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                a[i + j * lda]
                            else
                                numeric.conj(a[i + j * lda]),
                            x[numeric.cast(usize, ix)],
                            temp,
                        );
                    }

                    // x[jx] = temp
                    numeric.set(&x[numeric.cast(usize, jx)], temp);

                    jx += incx;
                }
            }
        }
    }

    return;
}
