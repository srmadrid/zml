const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Computes a matrix-matrix product with general matrices defined as:
///
/// ```zig
/// C = alpha * op(A) * op(B) + beta * C,
/// ```
///
/// where `alpha` and `beta` are numerics, `op(X)` is `X`, `Xᵀ`, `conj(X)` or
/// `Xᴴ`, op(A)` is an `m`-by-`k` matrix, `op(B)` is a `k`-by-`n`, and `C` is an
/// `m`-by-`n` matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.gemm(layout: matrix.Layout, transa: linalg.blas.Transpose, transb: linalg.blas.Transpose, m: usize, n: usize, k: usize, alpha: Al, a: [*]const A, lda: usize, b: [*]const B, ldb: usize, beta: Be, c: [*]C, ldc: usize) !void
/// ```
///
/// ## Parameters
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `transa` (`linalg.blas.Transpose`): Specifies the operation to be
///   performed on `A`:
///   * `no_transpose`: `C = alpha * A * op(B) + beta * C`
///   * `transpose`: `C = alpha * Aᵀ * op(B) + beta * C`
///   * `conj_no_transpose`: `C = alpha * conj(A) * op(B) + beta * C`
///   * `conj_transpose`: `C = alpha * Aᴴ * op(B) + beta * C`
/// * `transa` (`linalg.blas.Transpose`): Specifies the operation to be
///   performed on `B`:
///   * `no_transpose`: `C = alpha * op(A) * B + beta * C`
///   * `transpose`: `C = alpha * op(A) * Bᵀ + beta * C`
///   * `conj_no_transpose`: `C = alpha * op(A) * conj(B) + beta * C`
///   * `conj_transpose`: `C = alpha * op(A) * Bᴴ + beta * C`
/// * `m` (`usize`): Specifies the number of rows of the matrices `op(A)` and
///   `C`.
/// * `n` (`usize`): Specifies the number of columns of the matrices `op(B)` and
///   `C`.
/// * `k` (`usize`): Specifies the number of columns of the matrix `op(A)` and
///   the number of rows of the matrix `op(B)`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least:
///
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `lda * k`                                       | `lda * m`                                 |
/// | `order = row_major` | `lda * m`                                       | `lda * k`                                 |
///
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, m)`                                     | `max(1, k)`                               |
/// | `order = row_major` | `max(1, k)`                                     | `max(1, m)`                               |
///
/// * `b` (`anytype`): Many-item pointer, size at least:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `ldb * n`                                       | `ldb * k`                                 |
/// | `order = row_major` | `ldb * k`                                       | `ldb * n`                                 |
///
/// * `ldb` (`usize`): Specifies the leading dimension of `b` as declared in the
///   calling (sub)program. Must be greater than or equal to:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, k)`                                     | `max(1, n)`                               |
/// | `order = row_major` | `max(1, n)`                                     | `max(1, k)`                               |
///
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `c` need not be set on input.
/// * `c` (`anytype`): Many-item pointer, size at least `ldc * l`, where `l` is
///   `n` when `layout` is `col_major`, or `m` when `layout` is `row_major`.
/// * `ldc` (`usize`): Specifies the leading dimension of `c` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, n)` or
///   `max(1, k)`, if `ldb` is less than `max(1, k)` or `max(1, n)`, or if `ldc`
///   is less than `max(1, m)` or `max(1, n)`.
pub fn gemm(
    layout: matrix.Layout,
    transa: linalg.blas.Transpose,
    transb: linalg.blas.Transpose,
    m: usize,
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
        @compileError("zsl.linalg.blas.gemm: alpha and beta must be numerics, a and b must be many-item pointers to numerics, and c must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tb: " ++ @typeName(B) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\tc: " ++ @typeName(C) ++ "\n");

    A = meta.Child(A);
    B = meta.Child(B);
    C = meta.Child(C);

    const notransa = transa == .no_trans or transa == .conj_no_trans;
    const notransb = transb == .no_trans or transb == .conj_no_trans;
    const noconja = transa == .no_trans or transa == .trans;
    const noconjb = transb == .no_trans or transb == .trans;

    if (layout == .col_major) {
        const nrowa: usize = if (notransa) m else k;
        const nrowb: usize = if (notransb) k else n;

        if (lda < int.max(1, nrowa) or ldb < int.max(1, nrowb) or ldc < int.max(1, m))
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

        switch (noconja) {
            true => switch (noconjb) {
                true => return k_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, true, true),
                false => return k_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, true, false),
            },
            false => switch (noconjb) {
                true => return k_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, false, true),
                false => return k_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, false, false),
            },
        }
    } else {
        const nrowa: usize = if (notransa) k else m;
        const nrowb: usize = if (notransb) n else k;

        if (lda < int.max(1, nrowa) or ldb < int.max(1, nrowb) or ldc < int.max(1, n))
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

        switch (noconja) {
            true => switch (noconjb) {
                true => return k_gemm(transb, transa, n, m, k, alpha, b, ldb, a, lda, beta, c, ldc, true, true),
                false => return k_gemm(transb, transa, n, m, k, alpha, b, ldb, a, lda, beta, c, ldc, false, true),
            },
            false => switch (noconjb) {
                true => return k_gemm(transb, transa, n, m, k, alpha, b, ldb, a, lda, beta, c, ldc, true, false),
                false => return k_gemm(transb, transa, n, m, k, alpha, b, ldb, a, lda, beta, c, ldc, false, false),
            },
        }
    }
}

fn k_gemm(
    transa: linalg.blas.Transpose,
    transb: linalg.blas.Transpose,
    m: usize,
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
    comptime noconja: bool,
    comptime noconjb: bool,
) !void {
    const A: type = meta.Child(@TypeOf(a));
    const B: type = meta.Child(@TypeOf(b));

    // First form  C = beta * C.
    if (numeric.ne(beta, 1)) {
        for (0..n) |j| {
            linalg.blas.scal(m, beta, c + j * ldc, 1) catch unreachable;
        }
    }

    if (transb == .no_trans or transb == .conj_no_trans) {
        if (transa == .no_trans or transa == .conj_no_trans) {
            for (0..n) |j| {
                for (0..k) |l| {
                    // temp = alpha * b[l + j * ldb]
                    const temp = numeric.mul(
                        alpha,
                        if (comptime noconjb)
                            b[l + j * ldb]
                        else
                            numeric.conj(b[l + j * ldb]),
                    );

                    for (0..m) |i| {
                        // c[i + j * ldc] += temp * a[i + l * lda]
                        numeric.fmaInto(
                            &c[i + j * ldc],
                            temp,
                            if (comptime noconja)
                                a[i + l * lda]
                            else
                                numeric.conj(a[i + l * lda]),
                            c[i + j * ldc],
                        );
                    }
                }
            }
        } else {
            for (0..n) |j| {
                for (0..m) |i| {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, B)));

                    for (0..k) |l| {
                        // temp += a[l + i * lda] * b[l + j * ldb]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconja)
                                a[l + i * lda]
                            else
                                numeric.conj(a[l + i * lda]),
                            if (comptime noconjb)
                                b[l + j * ldb]
                            else
                                numeric.conj(b[l + j * ldb]),
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
    } else {
        if (transa == .no_trans or transa == .conj_no_trans) {
            for (0..n) |j| {
                for (0..k) |l| {
                    // temp = alpha * b[j + l * ldb]
                    const temp = numeric.mul(
                        alpha,
                        if (comptime noconjb)
                            b[j + l * ldb]
                        else
                            numeric.conj(b[j + l * ldb]),
                    );

                    for (0..m) |i| {
                        // c[i + j * ldc] += temp * a[i + l * lda]
                        numeric.fmaInto(
                            &c[i + j * ldc],
                            temp,
                            if (comptime noconja)
                                a[i + l * lda]
                            else
                                numeric.conj(a[i + l * lda]),
                            c[i + j * ldc],
                        );
                    }
                }
            }
        } else {
            for (0..n) |j| {
                for (0..m) |i| {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, B)));

                    for (0..k) |l| {
                        // temp += a[l + i * lda] * b[j + l * ldb]
                        numeric.fmaInto(
                            &temp,
                            if (comptime noconja)
                                a[l + i * lda]
                            else
                                numeric.conj(a[l + i * lda]),
                            if (comptime noconjb)
                                b[j + l * ldb]
                            else
                                numeric.conj(b[j + l * ldb]),
                            temp,
                        );
                    }

                    // c[i + j * ldc] += alpha * temp
                    numeric.fmaInto(
                        &c[i + j * ldc],
                        temp,
                        alpha,
                        c[i + j * ldc],
                    );
                }
            }
        }
    }

    return;
}
