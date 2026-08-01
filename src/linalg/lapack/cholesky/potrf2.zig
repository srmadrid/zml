const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Computes Cholesky factorization of a real symmetric or complex Hermitian
/// positive definite matrix, defined as:
///
/// ```zig
/// A = Uᵀ * U,
/// ```
///
/// or
///
/// ```zig
/// A = L * Lᵀ,
/// ```
///
/// or
///
/// ```zig
/// A = Uᴴ * U,
/// ```
///
/// or
///
/// ```zig
/// A = L * Lᴴ,
/// ```
///
/// where `A` is an `n × n` symmetric or Hermitian positive difinite matrix, `U`
/// is an `n × n` upper triangular matrix, and `L` is an `n × n` lower
/// triangular matrix.
///
/// This is the recursive version of the algorithm. It divides the matrix into
/// four submatrices:
///
/// ```zig
///     ┌          ┐
///     │ A₁₁  A₁₂ │
/// A = │ A₂₁  A₂₂ │
///     └          ┘
/// ```
///
/// where `A₁₁` is `n₁ × n₁` and `A₂₂` is `n₂ × n₂` with `n₁ = n / 2`, and
/// `n₂ = n - n₁`. The function calls itself to factor `A₁₁`,  update and scale
/// `A₂₁` or `A₁₂`, and update `A₂₂`, and then it calls itself to factor `A₂₂`.
///
/// ## Signature
/// ```zig
/// linalg.lapack.potrf2(layout: matrix.Layout, uplo: matrix.Uplo, n: usize, a: [*]A, lda: usize) !usize
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the matrix `A` is used, and which factorization is computed:
///   * `upper`: `A = Uᵀ * U` or `A = Uᴴ * U`
///   * `lower`: `A = L * Lᵀ` or `A = L * Lᴴ`
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `c` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// ## Returns
/// `usize`: `int.highest(usize)` if successful, or `i` if `uᵢᵢ` is exactly
/// zero.
///
/// ## Errors
/// * `linalg.lapack.Error.InvalidArgument`: If `lda` is less than `max(1, n)`.
pub fn potrf2(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    a: anytype,
    lda: usize,
) !usize {
    const A: type = @TypeOf(a);

    comptime if (!meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.lapack.potrf2: a must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n");

    if (lda < int.max(1, n))
        return linalg.lapack.Error.InvalidArgument;

    return k_potrf2(layout, uplo, n, a, lda);
}

fn k_potrf2(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    a: anytype,
    lda: usize,
) usize {
    // Quick return if possible.
    if (n == 0)
        return int.highest(usize);

    if (n == 1) {
        // Test for non-positive-definiteness.
        const ajj = numeric.re(a[0]);
        if (numeric.le(ajj, 0) or numeric.isNan(ajj))
            return 0;

        // Factor.
        numeric.sqrtInto(&a[0], ajj);
    } else {
        // Use recursive code.
        const n1 = int.div(n, 2);
        const n2 = n - n1;

        // Factor A₁₁.
        var iinfo = k_potrf2(
            layout,
            uplo,
            n1,
            a,
            lda,
        );

        if (iinfo != int.highest(usize))
            return iinfo;

        if (uplo == .upper) {
            // Compute the Cholesky factorization  A = Uᵀ * U  or  A = Uᴴ * U.

            // Update and scale A₁₂.
            linalg.blas.trsm(
                layout,
                .left,
                .upper,
                .conj_trans,
                .non_unit,
                n1,
                n2,
                1,
                a,
                lda,
                a + utils.index(layout, 0, n1, lda),
                lda,
            ) catch unreachable;

            // Update and factor A₂₂.
            linalg.blas.herk(
                layout,
                uplo,
                .conj_trans,
                n2,
                n1,
                -1,
                a + utils.index(layout, 0, n1, lda),
                lda,
                1,
                a + utils.index(layout, n1, n1, lda),
                lda,
            ) catch unreachable;

            iinfo = k_potrf2(
                layout,
                uplo,
                n2,
                a + utils.index(layout, n1, n1, lda),
                lda,
            ) catch unreachable;

            if (iinfo != int.highest(usize))
                return iinfo + n1;
        } else {
            // Compute the Cholesky factorization  A = L * Lᵀ  or  A = L * Lᴴ.

            // Update and scale A₂₁.
            linalg.blas.trsm(
                layout,
                .right,
                .lower,
                .conj_trans,
                .non_unit,
                n2,
                n1,
                1,
                a,
                lda,
                a + utils.index(layout, n1, 0, lda),
                lda,
            ) catch unreachable;

            // Update and factor A₂₂.
            linalg.blas.herk(
                layout,
                uplo,
                .no_trans,
                n2,
                n1,
                -1,
                a + utils.index(layout, n1, 0, lda),
                lda,
                1,
                a + utils.index(layout, n1, n1, lda),
                lda,
            ) catch unreachable;

            iinfo = k_potrf2(
                layout,
                uplo,
                n2,
                a + utils.index(layout, n1, n1, lda),
                lda,
            ) catch unreachable;

            if (iinfo != int.highest(usize))
                return iinfo + n1;
        }
    }

    return int.highest(usize);
}
