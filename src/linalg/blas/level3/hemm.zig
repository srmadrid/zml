const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Computes a matrix-matrix product where one input matrix is Hermitian defined
/// as:
///
/// ```zig
///     C = alpha * A * B + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * B * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are numerics, `A` is an `m`-by-`m` or `n`-by-`n`
/// Hermitian matrix, and `B` and `C` are `m`-by-`n` general matrices.
///
/// ## Signature
/// ```zig
/// linalg.blas.hemm(layout: matrix.Layout, side: linalg.blas.Side, uplo: matrix.Uplo, m: usize, n: usize, alpha: Al, a: [*]const A, lda: usize, b: [*]const B, ldb: usize, beta: Be, c: [*]C, ldc: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `side` (`linalg.blas.Side`): Specifies whether the Hermitian matrix `A`
///   appears on the left or right in the operation:
///   * `left`: `C = alpha * A * B + beta * C`.
///   * `right`: `C = alpha * B * A + beta * C`.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the Hermitian matrix `A` is used.
/// * `m` (`usize`): Specifies the number of rows of the matrix `C`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `C`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * ka`, where `ka`
///   is `m` if `side` is `left`, or `n` if `side` is `right`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `side` is `left`, or `max(1, n)` when `side` is `right`.
/// * `b` (`anytype`): Many-item pointer, size at least `ldb * kb`, where `kb`
///   is `n` if `layout` is `col_major`, or `m` if `layout` is `row_major`.
/// * `ldb` (`usize`): Specifies the leading dimension of `b` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `c` need not be set on input.
/// * `c` (`anytype`): Many-item pointer, size at least `ldc * kc`, where `kc`
///   is `n` when `layout` is `col_major`, or `m` when `layout` is `row_major`.
/// * `ldc` (`usize`): Specifies the leading dimension of `c` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, m)` or
///   `max(1, n)`, if `ldb` is less than `max(1, m)` or `max(1, n)`, or if `ldc`
///   is less than `max(1, m)` or `max(1, n)`.
pub fn hemm(
    layout: matrix.Layout,
    side: linalg.blas.Side,
    uplo: matrix.Uplo,
    m: usize,
    n: usize,
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
        @compileError("zsl.linalg.blas.hemm: alpha and beta must be numerics, a and b must be many-item pointers to numerics, and c must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tb: " ++ @typeName(B) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\tc: " ++ @typeName(C) ++ "\n");

    A = meta.Child(A);
    B = meta.Child(B);
    C = meta.Child(C);

    const nrowa: usize = if (side == .left) m else n;

    if (layout == .col_major) {
        if (lda < int.max(1, nrowa) or ldb < int.max(1, m) or ldc < int.max(1, m))
            return linalg.blas.Error.InvalidArgument;

        // Quick return if possible.
        if (m == 0 or n == 0)
            return;

        if (numeric.eq(alpha, 0)) {
            if (numeric.ne(beta, 1)) {
                for (0..n) |j| {
                    linalg.blas.scal(m, beta, c + j * ldc, 1) catch unreachable;
                }
            }
            return;
        }

        return k_hemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        if (lda < int.max(1, nrowa) or ldb < int.max(1, n) or ldc < int.max(1, n))
            return linalg.blas.Error.InvalidArgument;

        // Quick return if possible.
        if (m == 0 or n == 0)
            return;

        if (numeric.eq(alpha, 0)) {
            if (numeric.ne(beta, 1)) {
                for (0..m) |i| {
                    linalg.blas.scal(n, beta, c + i * ldc, 1) catch unreachable;
                }
            }
            return;
        }

        return k_hemm(side.invert(), uplo.invert(), n, m, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

fn k_hemm(
    side: linalg.blas.Side,
    uplo: matrix.Uplo,
    m: usize,
    n: usize,
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
            linalg.blas.scal(m, beta, c + j * ldc, 1) catch unreachable;
        }
    }

    if (side == .left) {
        if (uplo == .upper) {
            for (0..n) |j| {
                for (0..m) |i| {
                    // temp1 = alpha * b[i + j * ldb]
                    const temp1 = numeric.mul(alpha, b[i + j * ldb]);
                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(B, A)));

                    for (0..i) |k| {
                        // c[k + j * ldc] += temp1 * a[k + i * lda]
                        numeric.fmaInto(
                            &c[k + j * ldc],
                            temp1,
                            a[k + i * lda],
                            c[k + j * ldc],
                        );

                        // temp2 += b[k + j * ldb] * conj(a[k + i * lda])
                        numeric.fmaInto(
                            &temp2,
                            b[k + j * ldb],
                            numeric.conj(a[k + i * lda]),
                            temp2,
                        );
                    }

                    // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        temp1,
                        numeric.re(a[i + i * lda]),
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
                var i: usize = m;
                while (i > 0) {
                    i -= 1;

                    // temp1 = alpha * b[i + j * ldb]
                    const temp1 = numeric.mul(alpha, b[i + j * ldb]);
                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(B, A)));

                    for (i + 1..m) |k| {
                        // c[k + j * ldc] += temp1 * a[k + i * lda]
                        numeric.fmaInto(
                            &c[k + j * ldc],
                            temp1,
                            a[k + i * lda],
                            c[k + j * ldc],
                        );

                        // temp2 += b[k + j * ldb] * conj(a[k + i * lda])
                        numeric.fmaInto(
                            &temp2,
                            b[k + j * ldb],
                            numeric.conj(a[k + i * lda]),
                            temp2,
                        );
                    }

                    // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        temp1,
                        numeric.re(a[i + i * lda]),
                        numeric.fma(
                            alpha,
                            temp2,
                            c[i + j * ldc],
                        ),
                    );
                }
            }
        }
    } else {
        for (0..n) |j| {
            // temp1 = alpha * re(a[j + j * lda])
            var temp1 = numeric.mul(alpha, numeric.re(a[j + j * lda]));

            for (0..m) |i| {
                // c[i + j * ldc] += temp1 * b[i + j * ldb]
                numeric.fmaInto(
                    &c[i + j * ldc],
                    temp1,
                    b[i + j * ldb],
                    c[i + j * ldc],
                );
            }

            for (0..j) |k| {
                if (uplo == .upper) {
                    // temp1 = alpha * a[k + j * lda]
                    numeric.mulInto(
                        &temp1,
                        alpha,
                        a[k + j * lda],
                    );
                } else {
                    // temp1 = alpha * conj(a[j + k * lda])
                    numeric.mulInto(
                        &temp1,
                        alpha,
                        numeric.conj(a[j + k * lda]),
                    );
                }

                for (0..m) |i| {
                    // c[i + j * ldc] += temp1 * b[i + k * ldb]
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        temp1,
                        b[i + k * ldb],
                        c[i + j * ldc],
                    );
                }
            }

            for (j + 1..n) |k| {
                if (uplo == .upper) {
                    // temp1 = alpha * conj(a[j + k * lda])
                    numeric.mulInto(
                        &temp1,
                        alpha,
                        numeric.conj(a[j + k * lda]),
                    );
                } else {
                    // temp1 = alpha * a[k + j * lda]
                    numeric.mulInto(
                        &temp1,
                        alpha,
                        a[k + j * lda],
                    );
                }

                for (0..m) |i| {
                    // c[i + j * ldc] += temp1 * b[i + k * ldb]
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        temp1,
                        b[i + k * ldb],
                        c[i + j * ldc],
                    );
                }
            }
        }
    }

    return;
}
