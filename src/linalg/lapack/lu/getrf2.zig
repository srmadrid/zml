const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Computes the LU factorization of a general matrix, defined as:
///
/// ```zig
/// A = P L U,
/// ```
///
/// where `A` is an `m × n` matrix, `P` is an `m × m` permutation matrix, `L` is
/// an `m × min(m, n)` unit lower triangular matrix, and `U` is an
/// `min(m, n) × n` upper triangular matrix.
///
/// This is the recursive version of the algorithm. It divides the matrix into
/// four submatrices:
///
/// ```zig
///     ┌           ┐
///     │ A₁₁   A₁₂ │
/// A = │ A₂₁   A₂₂ │
///     └           ┘
/// ```
///
/// where `A₁₁` is `n₁ × n₁` and `A₂₂` is `n₂ × n₂` with `n₁ = min(m, n)`, and
/// `n₂ = n - n₁`. The function calls itself to factor `[ A₁₁   A₁₂ ]`,  does
/// the swaps on `[ A₁₂   A₂₂ ]ᵀ`, solves `A₁₂`, updates `A₂₂`, and then it
/// calls itself to factor `A₂₂` and do the swaps on `A₂₁`.
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
pub fn getrf2(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    a: anytype,
    lda: usize,
    ipiv: [*]usize,
) !usize {
    const A: type = @TypeOf(a);

    comptime if (!meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.lapack.getrf2: a must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n");

    if (lda < int.max(1, if (layout == .col_major) m else n))
        return linalg.lapack.Error.InvalidArgument;

    return k_getrf2(layout, m, n, a, lda, ipiv);
}

fn k_getrf2(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    a: anytype,
    lda: usize,
    ipiv: [*]usize,
) usize {
    const A: type = meta.Child(@TypeOf(a));

    var info = int.highest(usize);

    // Quick return if possible.
    if (m == 0 or n == 0)
        return info;

    if (m == 1) {
        // Use unblocked code for one row case. Just need to handle `ipiv` and `info`.
        ipiv[0] = 0;

        if (numeric.eq(a[0], 0))
            info = 1;
    } else if (n == 1) {
        // Use unblocked code for one column case.

        // Find pivot and test for singularity.
        const ii = linalg.blas.iamax(
            m,
            a,
            numeric.cast(isize, utils.col_ld(layout, lda)),
        ) catch unreachable;

        ipiv[0] = ii;

        if (numeric.ne(a[utils.index(layout, ii, 0, lda)], 0)) {
            // Apply the interchange.
            if (ii != 0) {
                const temp = a[0];
                a[0] = a[utils.index(layout, ii, 0, lda)];
                a[utils.index(layout, ii, 0, lda)] = temp;
            }

            // Compute elements `1..m` of the column.
            if (numeric.ge(numeric.abs1(a[0]), numeric.smallest(meta.Real(A)))) {
                linalg.blas.scal(
                    m - 1,
                    numeric.div(1, a[0]),
                    a + utils.index(layout, 1, 0, lda),
                    numeric.cast(isize, utils.col_ld(layout, lda)),
                ) catch unreachable;
            } else {
                for (1..m) |i| {
                    // `aᵢ₀ /= a₀`
                    numeric.divInto(
                        &a[utils.index(layout, i, 0, lda)],
                        a[utils.index(layout, i, 0, lda)],
                        a[0],
                    );
                }
            }
        } else {
            info = 0;
        }
    } else {
        // Use recursive code.
        const n1 = int.div(int.min(m, n), 2);
        const n2 = n - n1;

        //         ┌     ┐
        //         │ A₁₁ │
        // Factor  │ A₂₁ │
        //         └     ┘
        var iinfo = k_getrf2(
            layout,
            m,
            n1,
            a,
            lda,
            ipiv,
        );

        if (info == int.highest(usize) and iinfo != int.highest(usize))
            info = iinfo;

        //                        ┌     ┐
        //                        │ A₁₂ │
        // Apply interchanges to  │ A₂₂ │
        //                        └     ┘
        linalg.lapack.laswp(
            layout,
            n2,
            a + utils.index(layout, 0, n1, lda),
            lda,
            0,
            n1,
            ipiv,
            1,
        ) catch unreachable;

        // Solve `A₁₂`.
        linalg.blas.trsm(
            layout,
            .left,
            .lower,
            .no_trans,
            .unit,
            n1,
            n2,
            1,
            a,
            lda,
            a + utils.index(layout, 0, n1, lda),
            lda,
        ) catch unreachable;

        // Update `A₂₂`.
        linalg.blas.gemm(
            layout,
            .no_trans,
            .no_trans,
            m - n1,
            n2,
            n1,
            -1,
            a + utils.index(layout, n1, 0, lda),
            lda,
            a + utils.index(layout, 0, n1, lda),
            lda,
            1,
            a + utils.index(layout, n1, n1, lda),
            lda,
        ) catch unreachable;

        // Factor `A₂₂`.
        iinfo = k_getrf2(
            layout,
            m - n1,
            n2,
            a + utils.index(layout, n1, n1, lda),
            lda,
            ipiv + n1,
        );

        // Adjust info and the pivot indices.
        if (info == int.highest(usize) and iinfo != int.highest(usize))
            info = iinfo + n1;

        for (n1..int.min(m, n)) |i| {
            ipiv[i] += n1;
        }

        // Apply interchanges to `A₂₁`.
        linalg.lapack.laswp(
            layout,
            n1,
            a,
            lda,
            n1,
            int.min(m, n),
            ipiv,
            1,
        ) catch unreachable;
    }

    return info;
}
