const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Solves a system of linear equations with a Cholesky-factored symmetric or
/// Hermitian positive-definite coefficient matrix, defined as:
///
/// The `potrs` routine solves for `X` the system of linear equations
/// `A * X = B` with a symmetric positive-definite or, for complex data,
/// Hermitian positive-definite matrix `A`, given the Cholesky factorization of
/// `A`:
///
/// ```zig
/// A * X = B,
/// ```
///
/// where `A` is the Cholesky factorization of a symmetric or Hermitian
/// positive-definite `n × n` matrix `A`, as computed by `linalg.lapack.potrf`,
/// `B` is an `n × nrhs` matrix of right-hand sides, and `X` is an `n × nrhs`
/// matrix of solutions.
///
/// ## Signature
/// ```zig
/// linalg.lapack.potrs(layout: matrix.Layout, transa: linalg.blas.Transpose, n: usize, nrhs: usize, a: [*]const A, lda: usize, b: [*]B, ldb: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the matrix `A` is used, and which factorization has been computed:
///   * `upper`: `A = Uᵀ * U` or `A = Uᴴ * U`
///   * `lower`: `A = L * Lᵀ` or `A = L * Lᴴ`
/// * `n` (`usize`): Specifies the size of the matrix `A`, and the number of
///   rows of the matrices `B` and `X`.
/// * `nrhs` (`usize`): Specifies the number of columns of the matrices `B` and
///   `X`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
/// * `b` (`anytype`): Mutable many-item pointer, size at least `ldb * k`, where
///   `k` is `nrhs` when `layout` is `col_major`, or `n` when `layout` is
///   `row_major`. On return, contains the solution matrix `X`.
/// * `ldb` (`usize`): Specifies the leading dimension of `b` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)` when
///   `layout` is `col_major`, or `max(1, nrhs)` when `layout` is `row_major`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.lapack.Error.InvalidArgument`: If `lda` is less than `max(1, n)`,
///   or if `ldb` is less than `max(1, n)` or `max(1, nrhs)`.
pub fn potrs(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    nrhs: usize,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
) !void {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);

    comptime if (!meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(B) or meta.isConstPointer(B) or !meta.isNumeric(meta.Child(B)))
        @compileError("zsl.linalg.lapack.potrs: a must be a many-item pointer to numerics, and b must be a mutable many-item pointer to numerics, got \n\ta: " ++
            @typeName(A) ++ "\n\tb: " ++ @typeName(B) ++ "\n");

    if (lda < int.max(1, n) or ldb < int.max(1, if (layout == .col_major) n else nrhs))
        return linalg.lapack.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return;

    if (uplo == .upper) {
        // Solve  A * X = B, where  A = Uᵀ * U  or  A = Uᴴ * U.

        // Solve  Uᵀ * X = B  or  Uᴴ * X = B, overwriting B with X.
        linalg.blas.trsm(
            layout,
            .left,
            .upper,
            .conj_trans,
            .non_unit,
            n,
            nrhs,
            1,
            a,
            lda,
            b,
            ldb,
        ) catch unreachable;

        // Solve  U * X = B, overwriting B with X.
        linalg.blas.trsm(
            layout,
            .left,
            .upper,
            .no_trans,
            .non_unit,
            n,
            nrhs,
            1,
            a,
            lda,
            b,
            ldb,
        ) catch unreachable;
    } else {
        // Solve  A * X = B, where  A = L * Lᵀ  or  A = L * Lᴴ.

        // Solve  L * X = B, overwriting B with X.
        linalg.blas.trsm(
            layout,
            .left,
            .lower,
            .no_trans,
            .non_unit,
            n,
            nrhs,
            1,
            a,
            lda,
            b,
            ldb,
        ) catch unreachable;

        // Solve  Lᵀ * X = B  or  Lᴴ * X = B, overwriting B with X.
        linalg.blas.trsm(
            layout,
            .left,
            .lower,
            .conj_trans,
            .non_unit,
            n,
            nrhs,
            1,
            a,
            lda,
            b,
            ldb,
        ) catch unreachable;
    }

    return;
}
