const std = @import("std");

const types = @import("../types.zig");
const Cmp = types.Cmp;
const int = @import("../int.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Compares a `Rational` with another lower or equal precision numeric type for
/// ordering.
///
/// Signature
/// ---------
/// ```zig
/// fn cmp(x: X, y: Y) Cmp
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// Returns
/// -------
/// `Cmp`:
/// The comparison result.
pub fn cmp(x: anytype, y: anytype) Cmp {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .rational and types.numericType(Y) == .rational) and
        !(types.numericType(X) == .rational and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .rational and types.numericType(Y) == .float) and
        !(types.numericType(X) == .rational and types.numericType(Y) == .int) and
        !(types.numericType(X) == .rational and types.numericType(Y) == .bool) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .rational) and
        !(types.numericType(X) == .float and types.numericType(Y) == .rational) and
        !(types.numericType(X) == .int and types.numericType(Y) == .rational) and
        !(types.numericType(X) == .bool and types.numericType(Y) == .rational))
        @compileError("rational.cmp requires x or y to be a rational type, the other must be a rational, integer, float, int or bool type, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(X)) {
        .rational => switch (comptime types.numericType(Y)) {
            .rational => {
                if (x.num.positive != y.num.positive)
                    return if (!x.num.positive and y.num.positive) .lt else .gt;

                if (x.num.size == 0 or y.num.size == 0) {
                    if (x.num.size == 0 and y.num.size == 0) return .eq;

                    const cmp_abs: Cmp = if (x.num.size == 0) .lt else .gt;
                    return if (x.num.positive) cmp_abs else cmp_abs.invert();
                }

                const bit_size_nx: u32 = x.num.size * 32 - @clz(x.num.limbs[x.num.size - 1]);
                const bit_size_dx: u32 = x.den.size * 32 - @clz(x.den.limbs[x.den.size - 1]);
                const bit_size_ny: u32 = y.num.size * 32 - @clz(y.num.limbs[y.num.size - 1]);
                const bit_size_dy: u32 = y.den.size * 32 - @clz(y.den.limbs[y.den.size - 1]);

                const left_size: u32 = bit_size_nx + bit_size_dy;
                const right_size: u32 = bit_size_ny + bit_size_dx;

                if (left_size > right_size + 1) {
                    return if (x.num.positive) .gt else .lt;
                } else if (right_size > left_size + 1) {
                    return if (x.num.positive) .lt else .gt;
                }

                // a/b ? c/d  <=>  a*d ? c*b
                var carry: i128 = 0;
                //var last_rem: i128 = 0;
                var any_nonzero: bool = false;

                var k: u32 = 0;
                while (k < int.max(x.num.size + y.den.size, y.num.size + x.den.size)) : (k += 1) {
                    var acc: i128 = carry;

                    // acc += sum_{i + j = k} x.num[i] * y.den[j]
                    const i_min_ad: u32 = if (k >= y.den.size) k - (y.den.size - 1) else 0;
                    const i_max_ad: u32 = @min(k, x.num.size - 1);
                    var i: u32 = i_min_ad;
                    while (i <= i_max_ad) : (i += 1) {
                        const j: u32 = k - i;
                        acc += types.scast(i128, types.scast(u128, x.num.limbs[i]) * types.scast(u128, y.den.limbs[j]));
                    }

                    // acc -= sum_{i + j = k} y.num[i] * x.den[j]
                    const i_min_cb: u32 = if (k >= x.den.size) k - (x.den.size - 1) else 0;
                    const i_max_cb: u32 = @min(k, y.num.size - 1);
                    var t: u32 = i_min_cb;
                    while (t <= i_max_cb) : (t += 1) {
                        const j: u32 = k - t;
                        acc -= types.scast(i128, types.scast(u128, y.num.limbs[t]) * types.scast(u128, x.den.limbs[j]));
                    }

                    // Normalize: acc = carry_next * 2^32 + rem, with rem in [0, B)
                    // Use truncating division toward zero; then fix-up to make rem non-negative.
                    var carry_next: i128 = int.div(acc, 1 << 32);
                    var rem: i128 = acc - carry_next * @as(i128, 1 << 32);
                    if (rem < 0) {
                        rem += @as(i128, 1 << 32);
                        carry_next -= 1;
                    }

                    carry = carry_next;
                    any_nonzero = any_nonzero or (rem != 0); // track if any remainder limb is non-zero
                }

                // Decide sign of a*d - c*b
                var cmp_mag: Cmp = .eq;
                if (carry > 0) {
                    cmp_mag = .gt;
                } else if (carry < 0) {
                    cmp_mag = .lt;
                } else if (any_nonzero) {
                    cmp_mag = .gt; // any non-zero remainder means left > right since carry == 0
                } else {
                    cmp_mag = .eq;
                }

                // Adjust for overall sign (we handled opposite-sign early)
                return if (x.num.positive) cmp_mag else cmp_mag.invert();
            },
            .integer => return cmp(x, y.asRational()),
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return cmp(x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return cmp(x, ty[0]);
            },
            .bool => return cmp(x, types.cast(Rational, y, .{}) catch unreachable),
            else => unreachable,
        },
        .integer => switch (comptime types.numericType(Y)) {
            .rational => return cmp(x.asRational(), y),
            else => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .rational => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                return cmp(tx[0], y);
            },
            else => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .rational => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                return cmp(tx[0], y);
            },
            else => unreachable,
        },
        .bool => switch (comptime types.numericType(Y)) {
            .rational => return cmp(types.cast(Rational, x, .{}) catch unreachable, y),
            else => unreachable,
        },
        else => unreachable,
    }
}
