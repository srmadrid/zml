const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Computes the LU factorization of a general matrix, defined as:
///
/// ```zig
/// A = P * L * U,
/// ```
///
/// where `A` is an `m × n` matrix, `P` is an `m × m` permutation matrix, `L` is
/// an `m × min(m, n)` unit lower triangular matrix, and `U` is an
/// `min(m, n) × n` upper triangular matrix.
///
/// ## Signature
/// ```zig
/// linalg.lapack.getrf2(layout: matrix.Layout, m: usize, n: usize, a: [*]A, lda: usize, ipiv: [*]usize) !usize
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `m` (`usize`): Specifies the number of rows of the matrix `A`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * k`, where
///   `k` is `n` when `layout` is `col_major`, or `m` when `layout` is
///   `row_major`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
/// * `ipiv` (`[*]usize`): Mutable many-item pointer, size at least
///   `max(1, min(m, n))`. On return contains the pivot indices, i.e., for
///   `0 <= i < min(m, n)`, row `i` of the matrix was interchanged with row
///   `ipiv[i]`.
///
/// ## Returns
/// `usize`: `int.highest(usize)` if successful, or `i` if `uᵢᵢ` is exactly
/// zero.
///
/// ## Errors
/// * `linalg.lapack.Error.InvalidArgument`: If `lda` is less than `max(1, m)`
///   or `max(1, n)`.
pub fn getrf(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    a: anytype,
    lda: usize,
    ipiv: [*]usize,
) !usize {
    const A: type = @TypeOf(a);

    comptime if (!meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.lapack.getrf: a must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n");

    if (lda < int.max(1, if (layout == .col_major) m else n))
        return linalg.lapack.Error.InvalidArgument;

    return k_getrf(layout, m, n, a, lda, ipiv);
}

fn k_getrf(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    a: anytype,
    lda: usize,
    ipiv: [*]usize,
) usize {
    var info = int.highest(usize);

    // Quick return if possible.
    if (m == 0 or n == 0)
        return info;

    const nb = 64; // Adapt to cache size

    if (nb <= 1 or nb >= int.min(m, n)) {
        // Use unblocked code.
        info = linalg.lapack.getrf2(
            layout,
            m,
            n,
            a,
            lda,
            ipiv,
        ) catch unreachable;
    } else {
        // Use blocked code.
        var j: usize = 0;
        while (j < int.min(m, n)) : (j += nb) {
            const jb: usize = int.min(int.min(m, n) - j, nb);

            // Factor diagonal and subdiagonal blocks and test for exact singularity.
            const iinfo: usize = linalg.lapack.getrf2(
                layout,
                m - j,
                jb,
                a + utils.index(layout, j, j, lda),
                lda,
                ipiv + j,
            ) catch unreachable;

            // Adjust info and the pivot indices.
            if (info == int.highest(usize) and iinfo != int.highest(usize))
                info = iinfo + j;

            for (j..int.min(m, j + jb)) |i| {
                ipiv[i] += j;
            }

            // Apply interchanges to columns 0..j.
            linalg.lapack.laswp(
                layout,
                j,
                a,
                lda,
                j,
                j + jb,
                ipiv,
                1,
            ) catch unreachable;

            if (j + jb < n) {
                // Apply interchanges to columns j + jb..n.
                linalg.lapack.laswp(
                    layout,
                    n - j - jb,
                    a + utils.index(layout, 0, j + jb, lda),
                    lda,
                    j,
                    j + jb,
                    ipiv,
                    1,
                ) catch unreachable;

                // Compute block row of U.
                linalg.blas.trsm(
                    layout,
                    .left,
                    .lower,
                    .no_trans,
                    .unit,
                    jb,
                    n - j - jb,
                    1,
                    a + utils.index(layout, j, j, lda),
                    lda,
                    a + utils.index(layout, j, j + jb, lda),
                    lda,
                ) catch unreachable;

                if (j + jb < m) {
                    // Update trailing submatrix.
                    linalg.blas.gemm(
                        layout,
                        .no_trans,
                        .no_trans,
                        m - j - jb,
                        n - j - jb,
                        jb,
                        -1,
                        a + utils.index(layout, j + jb, j, lda),
                        lda,
                        a + utils.index(layout, j, j + jb, lda),
                        lda,
                        1,
                        a + utils.index(layout, j + jb, j + jb, lda),
                        lda,
                    ) catch unreachable;
                }
            }
        }
    }

    return info;
}
