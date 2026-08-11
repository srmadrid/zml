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
/// A = Uᵀ U,
/// ```
///
/// or
///
/// ```zig
/// A = L Lᵀ,
/// ```
///
/// or
///
/// ```zig
/// A = Uᴴ U,
/// ```
///
/// or
///
/// ```zig
/// A = L Lᴴ,
/// ```
///
/// where `A` is an `n × n` symmetric or Hermitian positive difinite matrix, `U`
/// is an `n × n` upper triangular matrix, and `L` is an `n × n` lower
/// triangular matrix.
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
///   * `upper`: `A = Uᵀ U` or `A = Uᴴ U`
///   * `lower`: `A = L Lᵀ` or `A = L Lᴴ`
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
pub fn potrf(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    a: anytype,
    lda: usize,
) !usize {
    const A: type = @TypeOf(a);

    comptime if (!meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.lapack.potrf: a must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n");

    if (lda < int.max(1, n))
        return linalg.lapack.Error.InvalidArgument;

    return k_potrf(layout, uplo, n, a, lda);
}

pub fn k_potrf(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    a: anytype,
    lda: usize,
) !usize {
    var info = int.highest(usize);

    // Quick return if possible.
    if (n == 0)
        return info;

    const nb = 64; // Adapt to cache size

    if (nb <= 1 or nb >= n) {
        // Use unblocked code.
        info = linalg.lapack.potrf2(
            layout,
            uplo,
            n,
            a,
            lda,
        ) catch unreachable;
    } else {
        // Use blocked code.
        if (uplo == .upper) {
            // Compute the Cholesky factorization `A = Uᵀ U` or `A = Uᴴ U`.
            var j: usize = 0;
            while (j < n) : (j += nb) {
                const jb: usize = int.min(n - j, nb);

                // Update and factorize the current diagonal block and test for non-positive-definiteness.
                linalg.blas.herk(
                    layout,
                    .upper,
                    .conj_trans,
                    jb,
                    j,
                    -1,
                    a + utils.index(layout, 0, j, lda),
                    lda,
                    1,
                    a + utils.index(layout, j, j, lda),
                    lda,
                ) catch unreachable;

                info = linalg.lapack.potrf2(
                    layout,
                    .upper,
                    jb,
                    a + utils.index(layout, j, j, lda),
                    lda,
                ) catch unreachable;

                if (info != int.highest(usize))
                    return info + j;

                if (j + jb < n) {
                    // Compute the current block row.
                    linalg.blas.gemm(
                        layout,
                        .conj_trans,
                        .no_trans,
                        jb,
                        n - j - jb,
                        j,
                        -1,
                        a + utils.index(layout, 0, j, lda),
                        lda,
                        a + utils.index(layout, 0, j + jb, lda),
                        lda,
                        1,
                        a + utils.index(layout, j, j + jb, lda),
                        lda,
                    ) catch unreachable;

                    linalg.blas.trsm(
                        layout,
                        .left,
                        .upper,
                        .conj_trans,
                        .non_unit,
                        jb,
                        n - j - jb,
                        1,
                        a + utils.index(layout, j, j, lda),
                        lda,
                        a + utils.index(layout, j, j + jb, lda),
                        lda,
                    ) catch unreachable;
                }
            }
        } else {
            // Compute the Cholesky factorization `A = L Lᵀ` or `A = L Lᴴ`.
            var j: usize = 0;
            while (j < n) : (j += nb) {
                const jb: usize = int.min(n - j, nb);

                // Update and factorize the current diagonal block and test for non-positive-definiteness.
                linalg.blas.herk(
                    layout,
                    .lower,
                    .no_trans,
                    jb,
                    j,
                    -1,
                    a + utils.index(layout, j, 0, lda),
                    lda,
                    1,
                    a + utils.index(layout, j, j, lda),
                    lda,
                ) catch unreachable;

                info = linalg.lapack.potrf2(
                    layout,
                    .lower,
                    jb,
                    a + utils.index(layout, j, j, lda),
                    lda,
                ) catch unreachable;

                if (info != int.highest(usize))
                    return info + j;

                if (j + jb < n) {
                    // Compute the current block column.
                    linalg.blas.gemm(
                        layout,
                        .no_trans,
                        .conj_trans,
                        n - j - jb,
                        jb,
                        j,
                        -1,
                        a + utils.index(layout, j + jb, 0, lda),
                        lda,
                        a + utils.index(layout, j, 0, lda),
                        lda,
                        1,
                        a + utils.index(layout, j + jb, j, lda),
                        lda,
                    ) catch unreachable;

                    linalg.blas.trsm(
                        layout,
                        .right,
                        .lower,
                        .conj_trans,
                        .non_unit,
                        n - j - jb,
                        jb,
                        1,
                        a + utils.index(layout, j, j, lda),
                        lda,
                        a + utils.index(layout, j + jb, j, lda),
                        lda,
                    ) catch unreachable;
                }
            }
        }
    }

    return info;
}
