const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Performs a symmetric rank-`k` update defined as:
///
/// ```zig
///     C = alpha * A * Aᵀ + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * Aᵀ * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are numerics, `A` is an `n`-by-`k` or `k`-by-`n`
/// matrix, and `C` is an `n`-by-`n` symmetric matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.syrk(layout: matrix.Layout, uplo: matrix.Uplo, transa: linal.blas.Transpose, n: usize, k: usize, alpha: Al, a: [*]const A, lda: usize, beta: Be, c: [*]C, ldc: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the symmetric matrix `C` is used.
/// * `transa` (`Transpose`): Specifies whether the conjugate transpose of
///   matrix `A` appears on the left or right in the operation:
///   * `no_trans`: `C = alpha * A * Aᵀ + beta * C`.
///   * `trans`: `C = alpha * Aᵀ * A + beta * C`.
///   * `conj_no_trans`: `C = alpha * conj(A) * Aᴴ + beta * C`.
///   * `conj_trans`: `C = alpha * Aᴴ * conj(A) + beta * C`.
/// * `n` (`usize`): Specifies the number of rows or columns of the matrix `A`
///   and the size of the matrix `C`.
/// * `k` (`usize`): Specifies the number of columns or rows of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * ka`, where ka` is
///   `k` when `transa` is `no_trans`, or `n` when `transa` is `conj_trans`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)` when
///   `transa` is `no_trans`, or `max(1, k)` when `transa` is `conj_trans`.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `c` need not be set on input.
/// * `c` (`anytype`): Many-item pointer, size at least `ldc * n`.
/// * `ldc` (`usize`): Specifies the leading dimension of `c` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, n)` or
///   `max(1, k)`, or if `ldc` is less than `max(1, n)`.
pub fn syrk(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    transa: linalg.blas.Transpose,
    n: usize,
    k: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    beta: anytype,
    c: anytype,
    ldc: usize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    const Be: type = @TypeOf(beta);
    comptime var C: type = @TypeOf(c);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(C) or meta.isConstPointer(C) or !meta.isNumeric(meta.Child(C)))
        @compileError("zsl.linalg.blas.syrk: alpha and beta must be numerics, a must be a many-item pointer to numerics, and c must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\tc: " ++ @typeName(C) ++ "\n");

    A = meta.Child(A);
    C = meta.Child(C);

    const nrowa: usize = if (transa == .no_trans) n else k;
    const ncola: usize = if (transa == .no_trans) k else n;

    if (layout == .col_major) {
        if (transa == .conj_no_trans or transa == .conj_trans or lda < int.max(1, nrowa) or ldc < int.max(1, n))
            return linalg.blas.Error.InvalidArgument;

        // Quick return if possible.
        if (n == 0)
            return;

        if (numeric.eq(alpha, 0)) {
            if (numeric.ne(beta, 1)) {
                for (0..n) |j| {
                    linalg.blas.scal(
                        if (uplo == .upper)
                            j + 1
                        else
                            n - j,
                        beta,
                        c + if (uplo == .upper)
                            j * ldc
                        else
                            j + j * ldc,
                        1,
                    ) catch unreachable;
                }
            }

            return;
        }

        return k_syrk(uplo, transa, n, k, alpha, a, lda, beta, c, ldc);
    } else {
        if (transa == .conj_no_trans or transa == .conj_trans or lda < int.max(1, ncola) or ldc < int.max(1, n))
            return linalg.blas.Error.InvalidArgument;

        // Quick return if possible.
        if (n == 0)
            return;

        if (numeric.eq(alpha, 0)) {
            if (numeric.ne(beta, 1)) {
                for (0..n) |i| {
                    linalg.blas.scal(
                        if (uplo == .upper)
                            n - i
                        else
                            i + 1,
                        beta,
                        c +
                            if (uplo == .upper)
                                i + i * ldc
                            else
                                i * ldc,
                        1,
                    ) catch unreachable;
                }
            }

            return;
        }

        return k_syrk(uplo.invert(), transa.invert(), n, k, alpha, a, lda, beta, c, ldc);
    }
}

fn k_syrk(
    uplo: matrix.Uplo,
    transa: linalg.blas.Transpose,
    n: usize,
    k: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    beta: anytype,
    c: anytype,
    ldc: usize,
) !void {
    const A: type = meta.Child(@TypeOf(a));

    // First form  C = beta * C.
    if (numeric.ne(beta, 1)) {
        for (0..n) |j| {
            linalg.blas.scal(
                if (uplo == .upper)
                    j + 1
                else
                    n - j,
                beta,
                c + if (uplo == .upper)
                    j * ldc
                else
                    j + j * ldc,
                1,
            ) catch unreachable;
        }
    }

    if (transa == .no_trans) {
        if (uplo == .upper) {
            for (0..n) |j| {
                for (0..k) |l| {
                    if (numeric.ne(a[j + l * lda], 0)) {
                        // temp = alpha * a[j + l * lda]
                        const temp = numeric.mul(alpha, a[j + l * lda]);

                        for (0..j + 1) |i| {
                            // c[i + j * ldc] += temp * a[i + l * lda]
                            numeric.fmaInto(
                                &c[i + j * ldc],
                                temp,
                                a[i + l * lda],
                                c[i + j * ldc],
                            );
                        }
                    }
                }
            }
        } else {
            for (0..n) |j| {
                for (0..k) |l| {
                    if (numeric.ne(a[j + l * lda], 0)) {
                        // temp = alpha * a[j + l * lda]
                        const temp = numeric.mul(
                            alpha,
                            a[j + l * lda],
                        );

                        for (j..n) |i| {
                            // c[i + j * ldc] += temp * a[i + l * lda]
                            numeric.fmaInto(
                                &c[i + j * ldc],
                                temp,
                                a[i + l * lda],
                                c[i + j * ldc],
                            );
                        }
                    }
                }
            }
        }
    } else {
        if (uplo == .upper) {
            for (0..n) |j| {
                for (0..j + 1) |i| {
                    var temp = numeric.cast(meta.Accumulator(numeric.Mul(A, A)), 0);

                    for (0..k) |l| {
                        // temp += a[l + i * lda] * a[l + j * lda]
                        numeric.fmaInto(
                            &temp,
                            a[l + i * lda],
                            a[l + j * lda],
                            temp,
                        );
                    }

                    // c[i + j * ldc] += alpha * temp
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        alpha,
                        temp,
                        c[i + j * ldc],
                    );
                }
            }
        } else {
            for (0..n) |j| {
                for (j..n) |i| {
                    var temp = numeric.cast(meta.Accumulator(numeric.Mul(A, A)), 0);

                    for (0..k) |l| {
                        // temp += a[l + i * lda] * a[l + j * lda]
                        numeric.fmaInto(
                            &temp,
                            a[l + i * lda],
                            a[l + j * lda],
                            temp,
                        );
                    }

                    // c[i + j * ldc] += alpha * temp
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        alpha,
                        temp,
                        c[i + j * ldc],
                    );
                }
            }
        }
    }

    return;
}
