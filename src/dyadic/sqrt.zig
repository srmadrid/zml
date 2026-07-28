const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const int = @import("../int.zig");

/// Returns the square root `√x` of a float.
///
/// ## Signature
/// ```zig
/// dyadic.sqrt(x: X) X
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the square root of.
///
/// ## Returns
/// `@TypeOf(x)`: The square root of `x`.
///
/// ## Method
/// The algorithm computes the exactly rounded square root of `x: Dyadic(N, E)`
/// by pulling the entire operation out of floating-point land and into the
/// integer domain [1]:
/// 1. Alignment Shift: To find the square root of a normalized input
///    `x = m × 2ᵉ`, we first shift the mantissa `m` left by a specific amount
///    `a`. We choose `a` to be either `N` or `N - 1` so that the new exponent
///    `(e - a)` is even. This scales the mantissa up to `n = m × 2ᵃ`, dropping
///    it into the bounds of `[2²ᴺ⁻², 2²ᴺ)`. By doing this shift, we
///    mathematically guarantee that our integer square root `⌊√n⌋` will land
///    right back in the normalized `N`-bit range `[2ᴺ⁻¹, 2ᴺ)`. This completely
///    eliminates the need for any post-normalization shifts [1].
/// 2. Integer Convergence: Next, we calculate the integer square root
///    `q = ⌊√n⌋` using an integer adaptation of Newton's method (Heron's
///    method), `qₖ₊₁ = ⌊(qₖ + n / qₖ) / 2⌋`. By seeding the loop with an
///    overestimation, the natural truncation of integer division guarantees
///    that the sequence will converge monotonically [2].
/// 3. Exact Rounding: Because our working values `n` and `q` are pure integers,
///    falling exactly on a `0.5` tie (`√n = q + 0.5`) is mathematically
///    impossible [1]. We only need to round up if the true root is strictly
///    greater than `q + 0.5`. Algebraically, this reduces to a simple check: we
///    round up if and only if the remainder `n − q² > q`.
///
/// ## References
/// [1] Muller, J.-M., et al. (2018). Handbook of Floating-Point Arithmetic (2nd ed.). Birkhäuser. https://doi.org/10.1007/978-3-319-76526-6
/// [2] Warren, H. S. (2013). Hacker's Delight (2nd ed.). Addison-Wesley Professional. ISBN: 978-0321842688
pub fn sqrt(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or meta.numericType(X) != .dyadic)
        @compileError("zsl.dyadic.sqrt: x must be a dyadic, got \n\tx: " ++ @typeName(X) ++ "\n");

    // NaN check
    if (x.isNan())
        return .nan;

    // Negative check
    if (!x.positive) {
        if (x.isZero())
            return x;

        return .nan;
    }

    // Inf check
    if (x.isInf())
        return .inf;

    // Zero check
    if (x.isZero())
        return x;

    const mantissa_bits = @typeInfo(X.Mantissa).int.bits;

    // For finite positive x = m * 2^e:
    //   sqrt(x) = sqrt(m << a) * 2^((e - a) / 2)
    // where a ∈ {N, N - 1} is chosen so (e - a) is even. Either choice puts
    // m = m << a in a range where isqrt(m) lands in [2^(N-1), 2^N), exactly
    // the N-bit normalized range, so no post-sqrt renormalization is needed.
    const e_is_odd: bool = @mod(x.exponent, 2) != 0;
    const shift_amount: u16 = if (e_is_odd != (comptime (mantissa_bits & 1) != 0)) mantissa_bits - 1 else mantissa_bits;

    const n: X.WideMantissa = numeric.cast(X.WideMantissa, x.mantissa) << @intCast(shift_amount);

    // Integer square root via Newton's method. Starting from an overestimation,
    // the iteration converges monotonically downward to floor(sqrt(n)) in
    // O(log(N)) steps.
    const q: X.WideMantissa = blk: {
        const msb_pos: u16 = 2 * mantissa_bits - 1 - @clz(n);
        var k: X.WideMantissa = @as(X.WideMantissa, 1) << @intCast(msb_pos / 2 + 1);
        while (true) {
            const next_k: X.WideMantissa = (k + n / k) / 2;

            if (next_k >= k)
                break :blk k;

            k = next_k;
        }
    };

    // For integer n, sqrt(n) is never exactly q + 0.5 (that would require
    // n = q^2 + q + 1/4, impossible for integer n), so there are no ties. The
    // decision reduces to:
    //   rem > q  iff  n ≥ q^2 + q + 1  iff  sqrt(n) > q + 0.5  iff  round up
    const rem: X.WideMantissa = n - q * q;

    var mantissa: X.Mantissa = @truncate(q);
    var result_exp: X.WideExponent = @divExact(
        numeric.cast(X.WideExponent, x.exponent) - numeric.cast(X.WideExponent, shift_amount),
        2,
    );

    if (rem > q) {
        const inc = @addWithOverflow(mantissa, 1);
        mantissa = inc[0];
        if (inc[1] != 0) {
            // q was 2^N - 1; rounding overflowed to 2^N. Renormalize.
            mantissa = @as(X.Mantissa, 1) << (mantissa_bits - 1);
            result_exp +|= 1;
        }
    }

    // Check for overflow
    if (result_exp >= int.highest(X.Exponent))
        return .{
            .mantissa = 0,
            .exponent = int.highest(X.Exponent),
            .positive = true,
        };

    // Check for underflow
    if (result_exp <= int.lowest(X.Exponent))
        return .{
            .mantissa = 0,
            .exponent = int.lowest(X.Exponent),
            .positive = true,
        };

    return .{
        .mantissa = mantissa,
        .exponent = numeric.cast(X.Exponent, result_exp),
        .positive = true,
    };
}
