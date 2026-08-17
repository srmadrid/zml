const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Forms the triangular factor `T` of a block reflector `H` of size `n`, which
/// is defined as a product of `k` elementary reflectors.
///
/// If `direction` is `forward`:
///
/// ```zig
/// H = H₁ H₂ ⋯ Hₖ,
/// ```
///
/// and `T` is upper triangular; If `direction` is `backward`:
///
/// ```zig
/// H = Hₖ ⋯ H₂ H₁,
/// ```
///
/// and `T` is lower triangular.
///
/// If `storev` is `columnwise`, the vector which defines the elementary
/// reflector `Hᵢ` is stored in the `i`-th column of the matrix `V`, and
///
/// ```zig
/// H = 𝕀 - V T Vᴴ.
/// ```
///
/// If `storev` is `rowwise`, the vector which defines the elementary reflector
/// `Hᵢ` is stored in the `i`-th row of the matrix `V`, and
///
/// ```zig
/// H = 𝕀 - Vᴴ T V.
/// ```
///
/// The shape of the matrix `V` and the storage of the vectors which define the
/// `Hᵢ` is best illustrated by the following example with `n = 5` and `k = 3`.
/// The elements equal to 1 are not stored. If `direction` is `forward` and
/// `storev` is `columnwise`:
///
/// ```zig
///     ┌              ┐
///     │  1           │
///     │ V₁    1      │
/// V = │ V₁   V₂    1 │,
///     │ V₁   V₂   V₃ │
///     │ V₁   V₂   V₃ │
///     └              ┘
/// ```
///
/// if `direction` is `forward` and `storev` is `rowwise`:
///
/// ```zig
///     ┌                        ┐
///     │  1   V₁   V₁   V₁   V₁ │
/// V = │       1   V₂   V₂   V₂ │,
///     │            1   V₃   V₃ │
///     └                        ┘
/// ```
///
/// if `direction` is `backward` and `storev` is `columnwise`:
///
/// ```zig
///     ┌              ┐
///     │ V₁   V₂   V₃ │
///     │ V₁   V₂   V₃ │
/// V = │  1   V₂   V₃ │,
///     │       1   V₃ │
///     │            1 │
///     └              ┘
/// ```
///
/// and if `direction` is `backward` and `storev` is `rowwise`:
///
/// ```zig
///     ┌                        ┐
///     │ V₁   V₁    1           │
/// V = │ V₂   V₂   V₂    1      │.
///     │ V₃   V₃   V₃   V₃    1 │
///     └                        ┘
/// ```
///
/// ## Signature
/// ```zig
/// linalg.lapack.larft(layout: matrix.Layout, direct: linalg.lapack.Direction, storev: linalg.lapack.Storage, n: usize, k: usize, v: [*]const V, ldv: usize, tau: [*]const Ta, t: [*]T, ldt: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `direct` (`linalg.lapack.Direction`): Specifies how `H` is formed from a
///   product of elementary reflectors:
///   * `forward`: `H = H₁ H₂ ⋯ Hₖ`.
///   * `backward`: `H = Hₖ ⋯ H₂ H₁`.
/// * `storev` (`linal.lapack.Storage`): Specifies how the vectors which define
///   the elementary reflectors are stored.
/// * `n` (`usize`): Specifies the size of the block reflector `H`.
/// * `k` (`usize`): Specifies the size of the triangular factor `T` and the
///   number of elementary reflectors.
/// * `v` (`anytype`): Many-item pointer, size at least `ldv * n`.
/// * `ldv` (`usize`): Specifies the leading dimension of `v` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, k)` when
///   `storage` is `columnwise`, or `max(1, n)` when `storage` is `rowwise`.
/// * `tau` (`anytype`): Many-item pointer, size at least `k`.
/// * `t` (`anytype`): Mutable many-item pointer, size at least `ldt * k`.
/// * `ldt` (`usize`): Specifies the leading dimension of `y` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, k)`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.lapack.Error.InvalidArgument`: If `ldv` is less than `max(1, k)`
///   or `max(1, n)`, or if `ldt` is less than `max(1, k)`.
pub fn larft(
    layout: matrix.Layout,
    direct: linalg.lapack.Direction,
    storev: linalg.lapack.Storage,
    n: usize,
    k: usize,
    v: anytype,
    ldv: usize,
    tau: anytype,
    t: anytype,
    ldt: usize,
) !void {
    const V: type = @TypeOf(v);
    const Ta: type = @TypeOf(tau);
    const T: type = @TypeOf(t);

    comptime if (!meta.isPointer(V) or !meta.isNumeric(meta.Child(V)) or
        !meta.isPointer(Ta) or !meta.isNumeric(meta.Child(Ta)) or
        !meta.isPointer(Ta) or meta.isConstPointer(Ta) or !meta.isNumeric(meta.Child(Ta)))
        @compileError("zsl.linalg.lapack.larft: v and tau must be many-item pointers to numerics, and t must be a mutable many-item pointer to numerics, got \n\tv: " ++ @typeName(V) ++ "\n\ttau: " ++ @typeName(Ta) ++ "\n\tt: " ++ @typeName(T) ++ "\n");

    if (ldv < int.max(1, if (storev == .columnwise) k else n) or ldt < int.max(1, k))
        return linalg.lapack.Error.InvalidArgument;

    return k_larft(layout, direct, storev, n, k, v, ldv, tau, t, ldt);
}

pub fn k_larft(
    layout: matrix.Layout,
    direct: linalg.lapack.Direction,
    storev: linalg.lapack.Storage,
    n: usize,
    k: usize,
    v: anytype,
    ldv: usize,
    tau: anytype,
    t: anytype,
    ldt: usize,
) void {
    // Quick return if possible.
    if (n == 0 or k == 0)
        return;

    if (n == 1 or k == 1) {
        numeric.set(&t[utils.index(layout, 0, 0, ldt)], tau[0]);

        return;
    }

    const l = k / 2;

    // Determine what kind of q we need to compute:
    // * QR: when we have forward direction in column storage.
    // * LQ: when we have forward direction in row storage.
    // * QL: when we have backward direction in column storage.
    // * RQ: when we have backward direction in row storage.
    const qr = direct == .forward and storev == .columnwise;
    const lq = direct == .forward and storev == .rowwise;
    const ql = direct == .backward and storev == .columnwise;

    if (qr) {
        // Break `V` apart into 6 components:
        //
        // ```zig
        //     ┌           ┐
        //     │ V₁₁     0 │
        // V = │ V₂₁   V₂₂ │,
        //     │ V₃₁   V₃₂ │
        //     └           ┘
        // ```
        //
        // where `V₁₁ ∈ ℂˡ×ˡ` unit lower triangular, `V₂₁ ∈ ℂᵏ⁻ˡ×ˡ`,
        // `V₃₁ ∈ ℂⁿ⁻ᵏ×ˡ`, `V₂₂ ∈ ℂᵏ⁻ˡ×ᵏ⁻ˡ` unit lower triangular, and
        // `V₃₂ ∈ ℂⁿ⁻ᵏ×ᵏ⁻ˡ`, with `l = ⌊ᵏ/₂⌋`. We will construct the matrix `T`
        // as:
        //
        // ```zig
        //     ┌           ┐
        //     │ T₁₁   T₁₂ │
        // T = │   0   T₂₂ │,
        //     └           ┘
        // ```
        //
        // where `T` is the triangular factor obtained from block reflectors,
        // `T₁₁ ∈ ℂˡ×ˡ` upper triangular, `T₂₂ ∈ ℂᵏ⁻ˡ×ᵏ⁻ˡ` upper triangular, and
        // `T₁₂ ∈ ℂˡ×ᵏ⁻ˡ`. To motivate the structure, assume we have already
        // computed `T₁₁` and `T₂₂`. Then collect the associated reflectors in
        // `V₁` and `V₂`. Consider the product:
        //
        // ```zig
        // (𝕀 - V₁ T₁₁ V₁ᴴ) (𝕀 - V₂ T₂₂ V₂ᴴ)
        //     = 𝕀 - V₁ T₁₁ V₁ᴴ - V₂ T₂₂ V₂ᴴ + V₁ T₁₁ V₁ᴴ V₂ T₂₂ V₂ᴴ.
        // ```
        //
        // Let `T₁₂ = -T₁₁ V₁ᴴ V₂ T₂₂`. Then, we can define the matrix `V` as:
        //
        // ```zig
        //     ┌         ┐
        // V = │ V₁   V₂ │,
        //     └         ┘
        // ```
        //
        // so, our product is equivalent to the matrix product:
        //
        // ```zig
        // 𝕀 - V T Vᴴ.
        // ```
        //
        // This means we can compute `T₁₁` and `T₂₂`, then use this information
        // to compute `T₁₂`.

        // Compute `T₁₁` recursively.
        k_larft(
            layout,
            direct,
            storev,
            n,
            l,
            v,
            ldv,
            tau,
            t,
            ldt,
        );

        // Compute `T₂₂` recursively.
        k_larft(
            layout,
            direct,
            storev,
            n - l,
            k - l,
            v + utils.index(layout, l, l, ldv),
            ldv,
            tau + l,
            t + utils.index(layout, l, l, ldt),
            ldt,
        );

        // Compute `T₁₂ = V₂₁ᴴ`.
        for (0..l) |j| {
            for (0..k - l) |i| {
                numeric.conjInto(
                    &t[utils.index(layout, j, l + i, ldt)],
                    v[utils.index(layout, l + i, j, ldv)],
                );
            }
        }

        // `T₁₂ = T₁₂ V₂₂`.
        linalg.blas.trmm(
            layout,
            .right,
            .lower,
            .no_trans,
            .unit,
            l,
            k - l,
            1,
            v + utils.index(layout, l, l, ldv),
            ldv,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;

        // `T₁₂ = V₃₁ᴴ V₃₂ + T₁₂`. We assume `k ≤ n`, and gemm will do nothing
        // if `n = k`.
        linalg.blas.gemm(
            layout,
            .conj_trans,
            .no_trans,
            l,
            k - l,
            n - k,
            1,
            v + utils.index(layout, k, 0, ldv),
            ldv,
            v + utils.index(layout, k, l, ldv),
            ldv,
            1,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;

        // At this point, we have that `T₁₂ = V₁ᴴ V₂`. All that is left to do is
        // to pre- and post-multiply by `-T₁₁` and `T₂₂`, respectively.

        // `T₁₂ = -T₁₁ * T₁₂`.
        linalg.blas.trmm(
            layout,
            .left,
            .upper,
            .no_trans,
            .non_unit,
            l,
            k - l,
            -1,
            t,
            ldt,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;

        // `T₁₂ = T₁₂ * T₂₂`.
        linalg.blas.trmm(
            layout,
            .right,
            .upper,
            .no_trans,
            .non_unit,
            l,
            k - l,
            1,
            t + utils.index(layout, l, l, ldt),
            ldt,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;
    } else if (lq) {
        // Break `V` apart into 6 components:
        //
        // ```zig
        //     ┌                 ┐
        //     │ V₁₁   V₁₂   V₁₃ │
        // V = │   0   V₂₂   V₂₃ │,
        //     └                 ┘
        // ```
        //
        // where ´V₁₁ ∈ ℂˡ×ˡ´ unit upper triangular, ´V₁₂ ∈ ℂˡ×ᵏ⁻ˡ´,
        // ´V₁₃ ∈ ℂˡ×ⁿ⁻ᵏ´, `V₂₂ ∈ ℂᵏ⁻ˡ×ᵏ⁻ˡ` unit upper triangular, and
        // `V₂₃ ∈ ℂᵏ⁻ˡ×ⁿ⁻ᵏ`, with `l = ⌊ᵏ/₂⌋`. We will construct the matrix `T`
        // as:
        //
        // ```zig
        //     ┌           ┐
        //     │ T₁₁   T₁₂ │
        // T = │   0   T₂₂ │,
        //     └           ┘
        // ```
        //
        // where `T` is the triangular factor obtained from block reflectors,
        // `T₁₁ ∈ ℂˡ×ˡ` upper triangular, `T₂₂ ∈ ℂᵏ⁻ˡ×ᵏ⁻ˡ` upper triangular, and
        // `T₁₂ ∈ ℂˡ×ᵏ⁻ˡ`. To motivate the structure, assume we have already
        // computed `T₁₁` and `T₂₂`. Then collect the associated reflectors in
        // `V₁` and `V₂`. Consider the product:
        //
        // ```zig
        // (𝕀 - V₁ᴴ T₁₁ V₁) (𝕀 - V₂ᴴ T₂₂ V₂)
        //     = 𝕀 - V₁ᴴ T₁₁ V₁ - V₂ᴴ T₂₂ V₂ + V₁ᴴ T₁₁ V₁ V₂ᴴ T₂₂ V₂.
        // ```
        //
        // Let `T₁₂ = -T₁₁ V₁ V₂ᴴ T₂₂`. Then, we can define the matrix `V` as:
        //
        // ```zig
        //     ┌    ┐
        //     │ V₁ │
        // V = │ V₂ │,
        //     └    ┘
        // ```
        //
        // so, our product is equivalent to the matrix product:
        //
        // ```zig
        // 𝕀 - Vᴴ T V.
        // ```
        //
        // This means we can compute `T₁₁` and `T₂₂`, then use this information
        // to compute `T₁₂`.

        // Compute `T₁₁` recursively.
        k_larft(
            layout,
            direct,
            storev,
            n,
            l,
            v,
            ldv,
            tau,
            t,
            ldt,
        );

        // Compute `T₂₂` recursively.
        try larft(
            layout,
            direct,
            storev,
            n - l,
            k - l,
            v + utils.index(layout, l, l, ldv),
            ldv,
            tau + l,
            t + utils.index(layout, l, l, ldt),
            ldt,
        );

        // Compute `T₁₂ = V₁₂`.
        linalg.lapack.lacpy(
            layout,
            .full,
            l,
            k - l,
            v + utils.index(layout, 0, l, ldv),
            ldv,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;

        // `T₁₂ = T₁₂ V₂₂ᴴ`.
        linalg.blas.trmm(
            layout,
            .right,
            .upper,
            .conj_trans,
            .unit,
            l,
            k - l,
            1,
            v + utils.index(layout, l, l, ldv),
            ldv,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;

        // `T₁₂ = V₁₃ V₂₃ᴴ + T₁₂`. We assume `k ≤ n`, and gemm will do nothing
        // if `n = k`.
        linalg.blas.gemm(
            layout,
            .no_trans,
            .conj_trans,
            l,
            k - l,
            n - k,
            1,
            v + utils.index(layout, 0, k, ldv),
            ldv,
            v + utils.index(layout, l, k, ldv),
            ldv,
            1,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;

        // At this point, we have that `T₁₂ = V₁ V₂ᴴ`. All that is left is to do
        // is to pre- and post-multiply by `-T₁₁` and `T₂₂`, respectively.

        // `T₁₂ = -T₁₁ T₁₂`.
        linalg.blas.trmm(
            layout,
            .left,
            .upper,
            .no_trans,
            .non_unit,
            l,
            k - l,
            -1,
            t,
            ldt,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;

        // `T₁₂ = T₁₂ T₂₂`.
        linalg.blas.trmm(
            layout,
            .right,
            .upper,
            .no_trans,
            .non_unit,
            l,
            k - l,
            1,
            t + utils.index(layout, l, l, ldt),
            ldt,
            t + utils.index(layout, 0, l, ldt),
            ldt,
        ) catch unreachable;
    } else if (ql) {
        // Break `V` apart into 6 components:
        //
        // ```zig
        //     ┌           ┐
        //     │ V₁₁   V₁₂ │
        // V = │ V₂₁   V₂₂ │,
        //     │   0   V₃₂ │
        //     └           ┘
        // ```
        //
        // where `V₁₁ ∈ ℂⁿ⁻ᵏ×ᵏ⁻ˡ`, `V₂₁ ∈ ℂᵏ⁻ˡ×ᵏ⁻ˡ` unit upper triangular,
        // `V₁₂ ∈ ℂⁿ⁻ᵏ×ˡ`, `V₂₂ ∈ ℂᵏ⁻ˡ×ˡ, and `V₃₂ ∈ ℂˡ×ˡ` unit upper
        // triangular, with `l = ⌊ᵏ/₂⌋`. We will construct the matrix `T` as:
        //
        // ```zig
        //     ┌           ┐
        //     │ T₁₁     0 │
        // T = │ T₂₁   T₂₂ │,
        //     └           ┘
        // ```
        //
        // where `T` is the triangular factor obtained from block reflectors,
        // `T₁₁ ∈ ℂᵏ⁻ˡ×ᵏ⁻ˡ` lower triangular, `T₂₂ ∈ ℂˡ×ˡ` lower triangular, and
        // `T₂₁ ∈ ℂᵏ⁻ˡ×ˡ`. To motivate the structure, assume we have already
        // computed `T₁₁` and `T₂₂`. then collect the associated reflectors in
        // `V₁` and `V₂`. Consider the product:
        //
        // ```zig
        // (𝕀 - V₂ T₂₂ V₂ᴴ) (𝕀 - V₁ T₁₁ V₁ᴴ)
        //     = 𝕀 - V₂ T₂₂ V₂ᴴ - V₁ T₁₁ V₁ᴴ + V₂ T₂₂ V₂ᴴ V₁ T₁₁ V₁ᴴ.
        // ```
        //
        // Let `T₂₁ = -T₂₂ V₂ᴴ V₁ T₁₁`. Then, we can define the matrix `V` as:
        //
        // ```zig
        //     ┌         ┐
        // V = │ V₁   V₂ │,
        //     └         ┘
        // ```
        //
        // so, our product is equivalent to the matrix product:
        //
        // ```zig
        // 𝕀 - V T Vᴴ.
        // ```
        //
        // This means we can compute `T₁₁` and `T₂₂`, then use this information
        // to compute `T₂₁`.

        // Compute `T₁₁` recursively.
        k_larft(
            layout,
            direct,
            storev,
            n - l,
            k - l,
            v,
            ldv,
            tau,
            t,
            ldt,
        );

        // Compute T₂₂ recursively
        k_larft(
            layout,
            direct,
            storev,
            n,
            l,
            v + utils.index(layout, 0, k - l, ldv),
            ldv,
            tau + k - l,
            t + utils.index(layout, k - l, k - l, ldt),
            ldt,
        );

        // Compute `T₂₁ = V₂₂ᴴ`.
        for (0..k - l) |j| {
            for (0..l) |i| {
                numeric.conjInto(
                    &t[utils.index(layout, k - l + i, j, ldt)],
                    v[utils.index(layout, n - k + j, k - l + i, ldv)],
                );
            }
        }

        // `T₂₁ = T₂₁ V₂₁`.
        linalg.blas.trmm(
            layout,
            .right,
            .upper,
            .no_trans,
            .unit,
            l,
            k - l,
            1,
            v + utils.index(layout, n - k, 0, ldv),
            ldv,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;

        // `T₂₁ = V₂₂ᴴ V₂₁ + T₂₁`. We assume `k ≤ n`, and gemm will do nothing
        // if `n = k`.
        linalg.blas.gemm(
            layout,
            .conj_trans,
            .no_trans,
            l,
            k - l,
            n - k,
            1,
            v + utils.index(layout, 0, k - l, ldv),
            ldv,
            v,
            ldv,
            1,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;

        // At this point, we have that `T₂₁ = V₂ᴴ V₁`. All that is left is to do
        // is to pre- and post-multiply by `-T₂₂` and `T₁₁` respectively.

        // `T₂₁ = -T₂₂ T₂₁`.
        linalg.blas.trmm(
            layout,
            .left,
            .lower,
            .no_trans,
            .non_unit,
            l,
            k - l,
            -1,
            t + utils.index(layout, k - l, k - l, ldt),
            ldt,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;

        // `T₂₁ = T₂₁ T₁₁`.
        linalg.blas.trmm(
            layout,
            .right,
            .lower,
            .no_trans,
            .non_unit,
            l,
            k - l,
            1,
            t,
            ldt,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;
    } else { // RQ
        // Break `V` apart into 6 components:
        //
        // ```zig
        //     ┌                 ┐
        //     │ V₁₁   V₁₂     0 │
        // V = │ V₂₁   V₂₂   V₂₃ │,
        //     └                 ┘
        // ```
        //
        // `V₁₁ ∈ ℂᵏ⁻ˡ×ⁿ⁻ᵏ`, `V₁₂ ∈ ℂᵏ⁻ˡ×ᵏ⁻ˡ` unit lower triangular,
        // `V₂₁ ∈ ℂˡ×ⁿ⁻ᵏ`, `V₂₂ ∈ ℂˡ×ᵏ⁻ˡ`, and `V₂₃ ∈ ℂˡ×ˡ` unit lower
        // triangular, with `l = ⌊ᵏ/₂⌋`. We will construct the matrix `T` as:
        //
        // ```zig
        //     ┌           ┐
        //     │ T₁₁     0 │
        // T = │ T₂₁   T₂₂ │,
        //     └           ┘
        // ```
        //
        // where `T` is the triangular factor obtained from block reflectors,
        // `T₁₁ ∈ ℂᵏ⁻ˡ×ᵏ⁻ˡ` lower triangular, `T₂₂ ∈ ℂˡ×ˡ` lower triangular, and
        // `T₂₁ ∈ ℂᵏ⁻ˡ×ˡ`. To motivate the structure, assume we have already
        // computed `T₁₁` and `T₂₂`. then collect the associated reflectors in
        // `V₁` and `V₂`. Consider the product:
        //
        // ```zig
        // (𝕀 - V₂ᴴ T₂₂ V₂) (𝕀 - V₁ᴴ T₁₁ V₁)
        //     = 𝕀 - V₂ᴴ T₂₂ V₂ - V₁ᴴ T₁₁ V₁ + V₂ᴴ T₂₂ V₂ V₁ᴴ T₁₁ V₁.
        // ```
        //
        // Let `T₂₁ = -T₂₂ V₂ V₁ᴴ T₁₁`. Then, we can define the matrix `V` as:
        //
        // ```zig
        //     ┌    ┐
        //     │ V₁ │
        // V = │ V₂ │,
        //     └    ┘
        // ```
        //
        // so, our product is equivalent to the matrix product:
        //
        // ```zig
        // 𝕀 - Vᴴ T V.
        // ```
        //
        // This means, we can compute `T₁₁` and `T₂₂`, then use this information
        // to compute `T₂₁`.

        // Compute `T₁₁` recursively.
        k_larft(
            layout,
            direct,
            storev,
            n - l,
            k - l,
            v,
            ldv,
            tau,
            t,
            ldt,
        );

        // Compute `T₂₂` recursively.
        k_larft(
            layout,
            direct,
            storev,
            n,
            l,
            v + utils.index(layout, k - l, 0, ldv),
            ldv,
            tau + k - l,
            t + utils.index(layout, k - l, k - l, ldt),
            ldt,
        );

        // Compute `T₂₁ = V₂₂`.
        linalg.lapack.lacpy(
            layout,
            .full,
            l,
            k - l,
            v + utils.index(layout, k - l, n - k, ldv),
            ldv,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;

        // `T₂₁ = T₂₁ V₁₂ᴴ`.
        linalg.blas.trmm(
            layout,
            .right,
            .lower,
            .conj_trans,
            .unit,
            l,
            k - l,
            1,
            v + utils.index(layout, 0, n - k, ldv),
            ldv,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;

        // `T₂₁ = V₂₁ V₁₁ᴴ + T₂₁`. We assume `k ≤ n`, and gemm will do nothing
        // if `n = k`.
        linalg.blas.gemm(
            layout,
            .no_trans,
            .conj_trans,
            l,
            k - l,
            n - k,
            1,
            v + utils.index(layout, k - l, 0, ldv),
            ldv,
            v,
            ldv,
            1,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;

        // At this point, we have that `T₂₁ = V₂ V₁ᴴ`. All that is left to do is
        // to pre- and post-multiply by `-T₂₂` and `T₁₁` respectively.

        // `T₂₁ = -T₂₂ T₂₁`.
        linalg.blas.trmm(
            layout,
            .left,
            .lower,
            .no_trans,
            .non_unit,
            l,
            k - l,
            -1,
            t + utils.index(layout, k - l, k - l, ldt),
            ldt,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;

        // `T₂₁ = T₂₁ T₁₁`.
        linalg.blas.trmm(
            layout,
            .right,
            .lower,
            .no_trans,
            .non_unit,
            l,
            k - l,
            1,
            t,
            ldt,
            t + utils.index(layout, k - l, 0, ldt),
            ldt,
        ) catch unreachable;
    }
}
