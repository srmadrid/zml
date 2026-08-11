const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Computes the QR factorization of a general matrix, defined as:
///
/// ```zig
/// A = Q R,
/// ```
///
/// where `A` is an `m × n` matrix, `Q` is an `m × m` unitary matrix, and `R` is
/// an `m × n` upper triangular matrix.
///
/// The function does not form the matrix `Q` explicitly. Instead, `Q` is
/// represented as a product of `min(m, n)` elementary reflectors:
///
/// ```zig
/// Q = H₁ H₂ ⋯ Hₖ,
/// ```
///
/// where `k = min(m, n)`. Each `Hᵢ` has the form
///
/// ```zig
/// Hᵢ = 𝕀 − τ v vᴴ,
/// ```
///
/// where `τ` is a numeric stored in `tau[i]`, and `v` is a vector with
/// `v[1..i - 1] = 0` and `vᵢ = 1`. On return, `v[i + 1..m]` is stored in
/// `a[i + 1..m, i]`.
///
/// ## Signature
/// ```zig
/// linalg.lapack.geqr2(layout: matrix.Layout, m: usize, n: usize, a: [*]A, lda: usize, tau: [*]Ta, work: [*]W) !usize
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
/// * `tau` (`anytype`): Mutable many-item pointer, size at least
///   `max(1, min(m, n))`.
/// * `work` (`anytype`): Mutable many-item pointer, size at least `max(1, n)`.
///
/// ## Returns
/// `usize`: `int.highest(usize)` if successful, or `i` if `uᵢᵢ` was an illegal
/// value.
///
/// ## Errors
/// * `linalg.lapack.Error.InvalidArgument`: If `lda` is less than `max(1, m)`
///   or `max(1, n)`.
pub fn geqr2(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    a: anytype,
    lda: usize,
    tau: anytype,
    work: anytype,
) !void {
    const A: type = @TypeOf(a);
    const Ta: type = @TypeOf(tau);

    comptime if (!meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(Ta) or meta.isConstPointer(Ta) or !meta.isNumeric(meta.Child(Ta)))
        @compileError("zsl.linalg.lapack.larfg: a and tau must be mutable many-item pointers to numerics, got \n\ta: " ++ @typeName(A) ++ "\n\ttau: " ++ @typeName(Ta) ++ "\n");

    if (lda < int.max(1, if (layout == .col_major) m else n))
        return linalg.lapack.Error.InvalidArgument;

    for (0..int.min(m, n)) |i| {
        // Generate elementary reflector `Hᵢ` to annihilate `A[i + 1..m, i]`.
        linalg.lapack.larfg(
            m - i,
            &a[utils.index(layout, i, i, lda)],
            a + utils.index(layout, int.min(i + 1, m - 1), i, lda),
            utils.col_ld(layout, lda),
            &tau[i],
        ) catch unreachable;

        if (i < n - 1) {
            // Apply `Hᵢ` to `A[i..m, i + 1..n]` from the left.
            linalg.lapack.larf1f(
                layout,
                .left,
                m - i,
                n - i - 1,
                a + utils.index(layout, i, i, lda),
                utils.col_ld(layout, lda),
                numeric.conj(tau[i]),
                a + utils.index(layout, i, i + 1, lda),
                lda,
                work,
            ) catch unreachable;
        }
    }
}
