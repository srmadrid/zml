const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Performs a symmetric rank-`2k` update defined as:
///
/// ```zig
///     C = alpha * A * Bᵀ + alpha * B * Aᵀ + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * Aᵀ * B + alpha * Bᵀ * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are numerics, `A` and `B` are `n`-by-`k` or
/// `k`-by-`n` matrices, and `C` is an `n`-by-`n` symmetric matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.syr2k(layout: matrix.Layout, uplo: matrix.Uplo, transab: linalg.blas.Transpose, n: usize, k: usize, alpha: Al, a: [*]const A, lda: usize, b: [*]const B, ldb: usize, beta: Be, c: [*]C, ldc: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the symmetric matrix `C` is used.
/// * `transab` (`Transpose`): Specifies whether the conjugate transpose of
///   matrix `A` appears on the left or right in the operation:
///   * `no_trans`: `C = alpha * A * Aᵀ + beta * C`.
///   * `conj_trans`: `C = alpha * Aᵀ * A + beta * C`.
/// * `n` (`usize`): Specifies the number of rows or columns of the matrix `A`
///   and the size of the matrix `C`.
/// * `k` (`usize`): Specifies the number of columns or rows of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * ka`, where ka` is
///   `k` when `transab` is `no_trans`, or `n` when `transab` is `conj_trans`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)` when
///   `transab` is `no_trans`, or `max(1, k)` when `transab` is `conj_trans`.
/// * `b` (`anytype`): Many-item pointer, size at least `ldb * kb`, where kb` is
///   `k` when `transab` is `no_trans`, or `n` when `transab` is `conj_trans`.
/// * `ldb` (`usize`): Specifies the leading dimension of `b` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)` when
///   `transab` is `no_trans`, or `max(1, k)` when `transab` is `conj_trans`.
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
/// * `linalg.blas.Error.InvalidArgument`: If `lda` or `ldb` is less than
///   `max(1, n)` or `max(1, k)`, or if `ldc` is less than `max(1, n)`.
pub fn syr2k(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    transab: linalg.blas.Transpose,
    n: usize,
    k: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
    beta: anytype,
    c: anytype,
    ldc: usize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);
    const Be: type = @TypeOf(beta);
    comptime var C: type = @TypeOf(c);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(B) or !meta.isNumeric(meta.Child(B)) or
        !meta.isManyItemPointer(C) or meta.isConstPointer(C) or !meta.isNumeric(meta.Child(C)))
        @compileError("zsl.linalg.blas.syr2k: alpha and beta must be numerics, a and b must be many-item pointers to numerics, and c must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tb: " ++ @typeName(B) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\tc: " ++ @typeName(C) ++ "\n");

    A = meta.Child(A);
    B = meta.Child(B);
    C = meta.Child(C);

    const nrowab: usize = if (transab == .no_trans) n else k;
    const ncolab: usize = if (transab == .no_trans) k else n;

    if (layout == .col_major) {
        if (transab == .conj_no_trans or transab == .conj_trans or lda < int.max(1, nrowab) or ldb < int.max(1, nrowab) or ldc < int.max(1, n))
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

        return k_syr2k(uplo, transab, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        if (transab == .conj_no_trans or transab == .conj_trans or lda < int.max(1, ncolab) or ldb < int.max(1, ncolab) or ldc < int.max(1, n))
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

        return k_syr2k(uplo.invert(), transab.invert(), n, k, alpha, b, ldb, a, lda, beta, c, ldc);
    }
}

fn k_syr2k(
    uplo: matrix.Uplo,
    transab: linalg.blas.Transpose,
    n: usize,
    k: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
    beta: anytype,
    c: anytype,
    ldc: usize,
) !void {
    const A: type = meta.Child(@TypeOf(a));
    const B: type = meta.Child(@TypeOf(b));

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

    if (transab == .no_trans) {
        if (uplo == .upper) {
            for (0..n) |j| {
                for (0..k) |l| {
                    if (numeric.ne(a[j + l * lda], 0) or numeric.ne(b[j + l * ldb], 0)) {
                        // temp1 = alpha * b[j + l * ldb]
                        const temp1 = numeric.mul(alpha, b[j + l * ldb]);

                        // temp2 = alpha * a[j + l * lda]
                        const temp2 = numeric.mul(alpha, a[j + l * lda]);

                        for (0..j + 1) |i| {
                            // c[i + j * ldc] += a[i + l * lda] * temp1 + b[i + l * ldb] * temp2
                            numeric.fmaInto(
                                &c[i + j * ldc],
                                a[i + l * lda],
                                temp1,
                                numeric.fma(
                                    b[i + l * ldb],
                                    temp2,
                                    c[i + j * ldc],
                                ),
                            );
                        }
                    }
                }
            }
        } else {
            for (0..n) |j| {
                for (0..k) |l| {
                    if (numeric.ne(a[j + l * lda], 0) or numeric.ne(b[j + l * ldb], 0)) {
                        // temp1 = alpha * b[j + l * ldb]
                        const temp1 = numeric.mul(alpha, b[j + l * ldb]);

                        // temp2 = alpha * a[j + l * lda]
                        const temp2 = numeric.mul(alpha, a[j + l * lda]);

                        for (j..n) |i| {
                            // c[i + j * ldc] += a[i + l * lda] * temp1 + b[i + l * ldb] * temp2
                            numeric.fmaInto(
                                &c[i + j * ldc],
                                a[i + l * lda],
                                temp1,
                                numeric.fma(
                                    b[i + l * ldb],
                                    temp2,
                                    c[i + j * ldc],
                                ),
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
                    var temp1 = numeric.zero(meta.Accumulator(numeric.Mul(A, B)));
                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(A, B)));

                    for (0..k) |l| {
                        // temp1 += a[l + i * lda] * b[l + j * ldb]
                        numeric.fmaInto(
                            &temp1,
                            a[l + i * lda],
                            b[l + j * ldb],
                            temp1,
                        );

                        // temp2 += b[l + i * ldb] * a[l + j * lda]
                        numeric.fmaInto(
                            &temp2,
                            b[l + i * ldb],
                            a[l + j * lda],
                            temp2,
                        );
                    }

                    // c[i + j * ldc] += alpha * temp1 + alpha * temp2
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        alpha,
                        temp1,
                        numeric.fma(
                            alpha,
                            temp2,
                            c[i + j * ldc],
                        ),
                    );
                }
            }
        } else {
            for (0..n) |j| {
                for (j..n) |i| {
                    var temp1 = numeric.zero(meta.Accumulator(numeric.Mul(A, B)));
                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(A, B)));

                    for (0..k) |l| {
                        // temp1 += a[l + i * lda] * b[l + j * ldb]
                        numeric.fmaInto(
                            &temp1,
                            a[l + i * lda],
                            b[l + j * ldb],
                            temp1,
                        );

                        // temp2 += b[l + i * ldb] * a[l + j * lda]
                        numeric.fmaInto(
                            &temp2,
                            b[l + i * ldb],
                            a[l + j * lda],
                            temp2,
                        );
                    }

                    // c[i + j * ldc] += alpha * temp1 + alpha * temp2
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        alpha,
                        temp1,
                        numeric.fma(
                            alpha,
                            temp2,
                            c[i + j * ldc],
                        ),
                    );
                }
            }
        }
    }

    return;
}
