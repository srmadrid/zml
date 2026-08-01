const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");

/// Solves a system of linear equations with an LU-factored square coefficient
/// matrix, defined as:
///
/// ```zig
/// A * X = B,
/// ```
///
/// or
///
/// ```zig
/// Aᵀ * X = B,
/// ```
///
/// or
///
/// ```zig
/// conj(A) * X = B,
/// ```
///
/// or
///
/// ```zig
/// Aᴴ * X = B,
/// ```
///
/// where `A` is the LU factorization of a general `n × n` matrix `A`, as
/// computed by `linalg.lapack.getrf`, `B` is an `n × nrhs` matrix of right-hand
/// sides, and `X` is an `n × nrhs` matrix of solutions.
///
/// ## Signature
/// ```zig
/// linalg.lapack.getrs(layout: matrix.Layout, transa: linalg.blas.Transpose, n: usize, nrhs: usize, a: [*]const A, lda: usize, ipiv: [*]const usize, b: [*]B, ldb: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `transa` (`linalg.blas.Transpose`): Specifies the operation to be
///   performed on `A`:
///   * `no_transpose`: `A * X  = B`
///   * `transpose`: `Aᵀ * X = B`
///   * `conj_no_transpose`: `conj(A) * X = B`
///   * `conj_transpose`: `Aᴴ * X = B`
/// * `n` (`usize`): Specifies the size of the matrix `A`, and the number of
///   rows of the matrices `B` and `X`.
/// * `nrhs` (`usize`): Specifies the number of columns of the matrices `B` and
///   `X`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
/// * `ipiv` (`[*]const usize`): Mutable many-item pointer, size at least
///   `max(1, min(m, n))`.
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
pub fn getrs(
    layout: matrix.Layout,
    transa: linalg.blas.Transpose,
    n: usize,
    nrhs: usize,
    a: anytype,
    lda: usize,
    ipiv: [*]const usize,
    b: anytype,
    ldb: usize,
) !void {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);

    comptime if (!meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(B) or meta.isConstPointer(B) or !meta.isNumeric(meta.Child(B)))
        @compileError("zsl.linalg.lapack.getrs: a must be a many-item pointer to numerics, and b must be a mutable many-item pointer to numerics, got \n\ta: " ++
            @typeName(A) ++ "\n\tb: " ++ @typeName(B) ++ "\n");

    if (lda < int.max(1, n) or ldb < int.max(1, if (layout == .col_major) n else nrhs))
        return linalg.lapack.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Solve  A * X = B  or  conj(A) * X = B.

        // Apply row interchanges to the right hand sides.
        linalg.lapack.laswp(
            layout,
            nrhs,
            b,
            ldb,
            0,
            n,
            ipiv,
            1,
        ) catch unreachable;

        // Solve  L * X = B  or  conj(L) * X = B, overwriting B with X.
        linalg.blas.trsm(
            layout,
            .left,
            .lower,
            transa,
            .unit,
            n,
            nrhs,
            1,
            a,
            lda,
            b,
            ldb,
        ) catch unreachable;

        // Solve  U * X = B  or  conj(U) * X = B, overwriting B with X.
        linalg.blas.trsm(
            layout,
            .left,
            .upper,
            transa,
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
        // Solve  Aᵀ * X = B  or  Aᴴ * X = B.

        // Solve  Uᵀ * X = B  or  Uᴴ * X = B, overwriting B with X.
        linalg.blas.trsm(
            layout,
            .left,
            .upper,
            transa,
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
            transa,
            .unit,
            n,
            nrhs,
            1,
            a,
            lda,
            b,
            ldb,
        ) catch unreachable;

        // Apply row interchanges to the solution vectors.
        linalg.lapack.laswp(
            layout,
            nrhs,
            b,
            ldb,
            0,
            n,
            ipiv,
            -1,
        ) catch unreachable;
    }

    return;
}
