const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Computes a matrix-matrix product where one input matrix is triangular
/// defined as:
///
/// ```zig
///     B = alpha * op(A) * B,
/// ```
///
/// or
///
/// ```zig
///     B = alpha * B * op(A),
/// ```
///
/// where `alpha` is a numeric, `op(A)` is `A`, `Aᵀ`, `conj(A)` oe `Aᴴ`, op(A)`
/// is an `m`-by-`m` or `n`-by-`n` matrix, and `B` is an `m`-by-`n` matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.trmm(layout: matrix.Layout, side: linalg.blas.Side, uplo: matrix.Uplo, transa: linalg.blas.Transpose, diag: matrix.Diag, m: usize, n: usize, alpha: Al, a: [*]const A, lda: usize, b: [*]B, ldb: usize) !void
/// ```
///
/// ## Parameters
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `side` (`linalg.blas.Side`): Specifies whether the triangular matrix `A`
///   appears on the left or right in the operation:
///   * `left`: `B = alpha * op(A) * B`.
///   * `right`: `B = alpha * B * op(A)`.
/// * `uplo` (`matrix.Uplo`): Specifies whether the matrix `A` is upper or lower
///   triangular.
/// * `transa` (`linalg.blas.Transpose`): Specifies the operation to be
///   performed on `A`:
///   * `no_transpose`: `B = alpha * A * B` or `B = alpha * B * A`
///   * `transpose`: `B = alpha * Aᵀ * B` or `B = alpha * B * Aᵀ`
///   * `conj_no_transpose`: `B = alpha * conj(A) * B` or `B = alpha * B * conj(A)`
///   * `conj_transpose`: `B = alpha * Aᴴ * B` or `B = alpha * B * Aᴴ`
/// * `diag` (`matrix.Diag`): Specifies whether the matrix `A` is unit
///   triangular.
/// * `m` (`usize`): Specifies the number of rows of the matrix `B`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `B`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * ka`, where `ka`
///   is `m` if `side` is `left`, or `n` if `side` is `right`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `side` is `left`, or `max(1, n)` when `side` is `right`.
/// * `b` (`anytype`): Many-item pointer, size at least `ldb * kb`, where `kb`
///   is `n` when `layout` is `col_major`, or `m` when `layout` is `row_major`.
/// * `ldb` (`usize`): Specifies the leading dimension of `b` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, m)` or
///   `max(1, n)`, or if `ldb` is less than `max(1, m)` or `max(1, n)`.
pub fn trmm(
    layout: matrix.Layout,
    side: linalg.blas.Side,
    uplo: matrix.Uplo,
    transa: linalg.blas.Transpose,
    diag: matrix.Diag,
    m: usize,
    n: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(B) or meta.isConstPointer(B) or !meta.isNumeric(meta.Child(B)))
        @compileError("zsl.linalg.blas.trmm: alpha must be a numeric, a must be a many-item pointer to numerics, and b must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tb: " ++ @typeName(B) ++ "\n");

    A = meta.Child(A);
    B = meta.Child(B);

    const noconja = transa == .no_trans or transa == .trans;
    const nrowa: usize = if (side == .left) m else n;

    if (layout == .col_major) {
        if (lda < int.max(1, nrowa) or ldb < int.max(1, m))
            return linalg.blas.Error.InvalidArgument;

        // Quick return if possible.
        if (m == 0 or n == 0)
            return;

        if (numeric.eq(alpha, 0)) {
            for (0..n) |j| {
                linalg.blas.scal(m, alpha, b + j * ldb, 1) catch unreachable;
            }

            return;
        }

        return if (noconja)
            k_trmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, true)
        else
            k_trmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, false);
    } else {
        if (lda < int.max(1, nrowa) or ldb < int.max(1, n))
            return linalg.blas.Error.InvalidArgument;

        // Quick return if possible.
        if (m == 0 or n == 0)
            return;

        if (numeric.eq(alpha, 0)) {
            for (0..m) |i| {
                linalg.blas.scal(n, alpha, b + i * ldb, 1) catch unreachable;
            }
            return;
        }

        return if (noconja)
            k_trmm(side.invert(), uplo.invert(), transa, diag, n, m, alpha, a, lda, b, ldb, true)
        else
            k_trmm(side.invert(), uplo.invert(), transa, diag, n, m, alpha, a, lda, b, ldb, false);
    }
}

fn k_trmm(
    side: linalg.blas.Side,
    uplo: matrix.Uplo,
    transa: linalg.blas.Transpose,
    diag: matrix.Diag,
    m: usize,
    n: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
    comptime noconja: bool,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const B: type = meta.Child(@TypeOf(b));

    if (side == .left) {
        if (transa == .no_trans or transa == .conj_no_trans) {
            if (uplo == .upper) {
                for (0..n) |j| {
                    for (0..m) |k| {
                        if (numeric.ne(b[k + j * ldb], 0)) {
                            const temp = numeric.mul(
                                alpha,
                                b[k + j * ldb],
                            );

                            for (0..k) |i| {
                                // b[i + j * ldb] += temp * a[i + k * lda]
                                numeric.fmaInto(
                                    &b[i + j * ldb],
                                    temp,
                                    if (comptime noconja)
                                        a[i + k * lda]
                                    else
                                        numeric.conj(a[i + k * lda]),
                                    b[i + j * ldb],
                                );
                            }

                            if (diag == .non_unit) {
                                // b[k + j * ldb] = temp * a[k + k * lda]
                                numeric.mulInto(
                                    &b[k + j * ldb],
                                    temp,
                                    if (comptime noconja)
                                        a[k + k * lda]
                                    else
                                        numeric.conj(a[k + k * lda]),
                                );
                            } else {
                                numeric.set(&b[k + j * ldb], temp);
                            }
                        }
                    }
                }
            } else {
                for (0..n) |j| {
                    var k: usize = m;
                    while (k > 0) {
                        k -= 1;

                        if (numeric.ne(b[k + j * ldb], 0)) {
                            const temp = numeric.mul(
                                alpha,
                                b[k + j * ldb],
                            );

                            numeric.set(&b[k + j * ldb], temp);

                            if (diag == .non_unit) {
                                // b[k + j * ldb] *= a[k + k * lda]
                                numeric.mulInto(
                                    &b[k + j * ldb],
                                    b[k + j * ldb],
                                    if (comptime noconja)
                                        a[k + k * lda]
                                    else
                                        numeric.conj(a[k + k * lda]),
                                );
                            }

                            for (k + 1..m) |i| {
                                // b[i + j * ldb] += temp * a[i + k * lda]
                                numeric.fmaInto(
                                    &b[i + j * ldb],
                                    temp,
                                    if (comptime noconja)
                                        a[i + k * lda]
                                    else
                                        numeric.conj(a[i + k * lda]),
                                    b[i + j * ldb],
                                );
                            }
                        }
                    }
                }
            }
        } else {
            if (uplo == .upper) {
                for (0..n) |j| {
                    var i: usize = m;
                    while (i > 0) {
                        i -= 1;

                        var temp = numeric.cast(meta.Accumulator(numeric.Mul(A, B)), b[i + j * ldb]);

                        if (diag == .non_unit) {
                            // temp *= a[i + i * lda]
                            numeric.mulInto(
                                &temp,
                                if (comptime noconja)
                                    a[i + i * lda]
                                else
                                    numeric.conj(a[i + i * lda]),
                                temp,
                            );
                        }

                        for (0..i) |k| {
                            // temp += a[k + i * lda] * b[k + j * ldb]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconja)
                                    a[k + i * lda]
                                else
                                    numeric.conj(a[k + i * lda]),
                                b[k + j * ldb],
                                temp,
                            );
                        }

                        // b[i + j * ldb] = alpha * temp
                        numeric.mulInto(
                            &b[i + j * ldb],
                            alpha,
                            temp,
                        );
                    }
                }
            } else {
                for (0..n) |j| {
                    for (0..m) |i| {
                        var temp = numeric.cast(meta.Accumulator(numeric.Mul(A, B)), b[i + j * ldb]);

                        if (diag == .non_unit) {
                            // temp *= a[i + i * lda]
                            numeric.mulInto(
                                &temp,
                                temp,
                                if (comptime noconja)
                                    a[i + i * lda]
                                else
                                    numeric.conj(a[i + i * lda]),
                            );
                        }

                        for (i + 1..m) |k| {
                            // temp += a[k + i * lda] * b[k + j * ldb]
                            numeric.fmaInto(
                                &temp,
                                if (comptime noconja)
                                    a[k + i * lda]
                                else
                                    numeric.conj(a[k + i * lda]),
                                b[k + j * ldb],
                                temp,
                            );
                        }

                        // b[i + j * ldb] = alpha * temp
                        numeric.mulInto(
                            &b[i + j * ldb],
                            alpha,
                            temp,
                        );
                    }
                }
            }
        }
    } else {
        if (transa == .no_trans or transa == .conj_no_trans) {
            if (uplo == .upper) {
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.cast(numeric.Mul(Al, A), alpha);

                    if (diag == .non_unit) {
                        // temp *= a[j + j * lda]
                        numeric.mulInto(
                            &temp,
                            temp,
                            if (comptime noconja)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                        );
                    }

                    for (0..m) |i| {
                        // b[i + j * ldb] *= temp
                        numeric.mulInto(
                            &b[i + j * ldb],
                            b[i + j * ldb],
                            temp,
                        );
                    }

                    for (0..j) |k| {
                        if (numeric.ne(a[k + j * lda], 0)) {
                            // temp = alpha * a[k + j * lda]
                            numeric.mulInto(
                                &temp,
                                alpha,
                                if (comptime noconja)
                                    a[k + j * lda]
                                else
                                    numeric.conj(a[k + j * lda]),
                            );

                            for (0..m) |i| {
                                // b[i + j * ldb] += temp * b[i + k * ldb]
                                numeric.fmaInto(
                                    &b[i + j * ldb],
                                    temp,
                                    b[i + k * ldb],
                                    b[i + j * ldb],
                                );
                            }
                        }
                    }
                }
            } else {
                for (0..n) |j| {
                    var temp = numeric.cast(numeric.Mul(Al, A), alpha);

                    if (diag == .non_unit) {
                        // temp *= a[j + j * lda]
                        numeric.mulInto(
                            &temp,
                            temp,
                            if (comptime noconja)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                        );
                    }

                    for (0..m) |i| {
                        // b[i + j * ldb] *= temp
                        numeric.mulInto(
                            &b[i + j * ldb],
                            b[i + j * ldb],
                            temp,
                        );
                    }

                    for (j + 1..n) |k| {
                        if (numeric.ne(a[k + j * lda], 0)) {
                            // temp = alpha * a[k + j * lda]
                            numeric.mulInto(
                                &temp,
                                alpha,
                                if (comptime noconja)
                                    a[k + j * lda]
                                else
                                    numeric.conj(a[k + j * lda]),
                            );

                            for (0..m) |i| {
                                // b[i + j * ldb] += temp * b[i + k * ldb]
                                numeric.fmaInto(
                                    &b[i + j * ldb],
                                    temp,
                                    b[i + k * ldb],
                                    b[i + j * ldb],
                                );
                            }
                        }
                    }
                }
            }
        } else {
            if (uplo == .upper) {
                for (0..n) |k| {
                    for (0..k) |j| {
                        if (numeric.ne(a[j + k * lda], 0)) {
                            var temp = numeric.cast(numeric.Mul(Al, A), alpha);

                            // temp *= a[j + k * lda]
                            numeric.mulInto(
                                &temp,
                                temp,
                                if (comptime noconja)
                                    a[j + k * lda]
                                else
                                    numeric.conj(a[j + k * lda]),
                            );

                            for (0..m) |i| {
                                // b[i + j * ldb] += temp * b[i + k * ldb]
                                numeric.fmaInto(
                                    &b[i + j * ldb],
                                    temp,
                                    b[i + k * ldb],
                                    b[i + j * ldb],
                                );
                            }
                        }
                    }

                    var temp = numeric.cast(numeric.Mul(Al, A), alpha);

                    if (diag == .non_unit) {
                        // temp *= a[k + k * lda]
                        numeric.mulInto(
                            &temp,
                            temp,
                            if (comptime noconja)
                                a[k + k * lda]
                            else
                                numeric.conj(a[k + k * lda]),
                        );
                    }

                    if (numeric.ne(temp, 1)) {
                        for (0..m) |i| {
                            // b[i + k * ldb] *= temp
                            numeric.mulInto(
                                &b[i + k * ldb],
                                b[i + k * ldb],
                                temp,
                            );
                        }
                    }
                }
            } else {
                var k: usize = n;
                while (k > 0) {
                    k -= 1;

                    for (k + 1..n) |j| {
                        if (numeric.ne(a[j + k * lda], 0)) {
                            var temp = numeric.cast(numeric.Mul(Al, A), alpha);

                            // temp *= a[j + k * lda]
                            numeric.mulInto(
                                &temp,
                                temp,
                                if (comptime noconja)
                                    a[j + k * lda]
                                else
                                    numeric.conj(a[j + k * lda]),
                            );

                            for (0..m) |i| {
                                // b[i + j * ldb] += temp * b[i + k * ldb]
                                numeric.fmaInto(
                                    &b[i + j * ldb],
                                    temp,
                                    b[i + k * ldb],
                                    b[i + j * ldb],
                                );
                            }
                        }
                    }

                    var temp = numeric.cast(numeric.Mul(Al, A), alpha);

                    if (diag == .non_unit) {
                        // temp *= a[k + k * lda]
                        numeric.mulInto(
                            &temp,
                            temp,
                            if (comptime noconja)
                                a[k + k * lda]
                            else
                                numeric.conj(a[k + k * lda]),
                        );
                    }

                    if (numeric.ne(temp, 1)) {
                        for (0..m) |i| {
                            // b[i + k * ldb] *= temp
                            numeric.mulInto(
                                &b[i + k * ldb],
                                b[i + k * ldb],
                                temp,
                            );
                        }
                    }
                }
            }
        }
    }

    return;
}
