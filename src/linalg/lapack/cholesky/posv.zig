const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");

/// Solves a system of linear equations with a symmetric or Hermitian
/// positive-definite coefficient matrix, defined as:
///
/// ```zig
/// A * X = B,
/// ```
///
/// where `A` is a symmetric or Hermitian positive-definite `n × n` matrix, `B`
/// is an `n × nrhs` matrix of right-hand sides, and `X` is an `n × nrhs` matrix
/// of solutions.
///
/// ## Signature
/// ```zig
/// linalg.lapackposv(layout: matrix.Layout, uplo: matrix.Uplo, n: usize, nrhs: usize, a: [*]A, lda: usize, b: [*]B, ldb: usize) !usize
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
pub fn posv(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    nrhs: usize,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
) !i32 {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);

    comptime if (!meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(B) or meta.isConstPointer(B) or !meta.isNumeric(meta.Child(B)))
        @compileError("zsl.linalg.lapack.posv: a must be a many-item pointer to numerics, and b must be a mutable many-item pointer to numerics, got \n\ta: " ++
            @typeName(A) ++ "\n\tb: " ++ @typeName(B) ++ "\n");

    if (lda < int.max(1, n) or ldb < int.max(1, if (layout == .col_major) n else nrhs))
        return linalg.lapack.Error.InvalidArgument;

    var info = int.highest(usize);

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return info;

    // Compute the Cholesky factorization  A = Uᵀ * U  or  A = Uᴴ * U  or  A = L * Lᵀ  or  A = L * Lᴴ.
    info = linalg.lapack.potrf(
        layout,
        uplo,
        n,
        a,
        lda,
    ) catch unreachable;

    if (info == int.highest(usize)) {
        // Solve the system  A * X = B, overwriting B with X.
        linalg.lapack.potrs(
            layout,
            uplo,
            n,
            nrhs,
            a,
            lda,
            b,
            ldb,
        ) catch unreachable;
    }

    return info;
}
