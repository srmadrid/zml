//! Namespace for dyadic operations.

const dyadic = @This();

const std = @import("std");

const meta = @import("meta.zig");

const numeric = @import("numeric.zig");

const int = @import("int.zig");

/// Arbitrary-precision dyadic type.
pub fn Dyadic(mantissa_bits: u16, exponent_bits: u16) type {
    if (mantissa_bits == 0 or exponent_bits == 0 or
        mantissa_bits >= int.highest(u16) / 2 or exponent_bits >= int.highest(u16) / 2)
        @compileError(std.fmt.comptimePrint(
            "zsl.Dyadic: both mantissa_bits and exponent_bits must be non-zero and less than {}, got\n\tmantissa_bits: {}\n\texponent_bits: {}\n",
            .{ int.highest(u16) / 2, mantissa_bits, exponent_bits },
        ));

    return struct {
        mantissa: Mantissa,
        exponent: Exponent,
        positive: bool,

        // Type signature
        pub const is_numeric = true;

        pub const Accumulator = blk: {
            if (mantissa_bits <= 16)
                break :blk Dyadic(32, exponent_bits)
            else
                break :blk Dyadic(mantissa_bits, exponent_bits);
        };

        pub const Mantissa = @Int(.unsigned, mantissa_bits);
        pub const WideMantissa = @Int(.unsigned, 2 * mantissa_bits);
        pub const Exponent = @Int(.signed, exponent_bits);
        pub const WideExponent = @Int(.signed, 2 * exponent_bits);

        // Constants
        pub const inf: Dyadic(mantissa_bits, exponent_bits) = .{ .mantissa = 0, .exponent = int.highest(Exponent), .positive = true };
        pub const nan: Dyadic(mantissa_bits, exponent_bits) = .{ .mantissa = 1, .exponent = int.highest(Exponent), .positive = true };
        pub const zero: Dyadic(mantissa_bits, exponent_bits) = .{ .mantissa = 0, .exponent = int.lowest(Exponent), .positive = true };
        pub const one: Dyadic(mantissa_bits, exponent_bits) = .{ .mantissa = @as(Mantissa, 1) << (mantissa_bits - 1), .exponent = -numeric.cast(Exponent, mantissa_bits - 1), .positive = true };
        pub const two: Dyadic(mantissa_bits, exponent_bits) = .{ .mantissa = @as(Mantissa, 1) << (mantissa_bits - 1), .exponent = -numeric.cast(Exponent, mantissa_bits), .positive = true };

        pub fn isInf(self: Dyadic(mantissa_bits, exponent_bits)) bool {
            return self.exponent == int.highest(Exponent) and self.mantissa == 0;
        }

        pub fn isNan(self: Dyadic(mantissa_bits, exponent_bits)) bool {
            return self.exponent == int.highest(Exponent) and self.mantissa != 0;
        }

        pub fn isZero(self: Dyadic(mantissa_bits, exponent_bits)) bool {
            return self.exponent == int.lowest(Exponent) and self.mantissa == 0;
        }

        /// Initializes a dyadic from any numeric value.
        ///
        /// ## Arguments
        /// * `value` (`anytype`): The value to set the dyadic to. Must be a
        ///   numeric.
        ///
        /// ## Returns
        /// `Dyadic(mantissa_bits, exponent_bits)`: The new dyadic.
        pub fn initValue(value: anytype) Dyadic(mantissa_bits, exponent_bits) {
            const V: type = @TypeOf(value);

            comptime if (!meta.isNumeric(V))
                @compileError("zsl.Dyadic(mantissa_bits, exponent_bits).initValue: value must be a numeric, got \n\tvalue: " ++ @typeName(V) ++ "\n");

            switch (comptime meta.numericType(V)) {
                .bool => return if (value) .one else .zero,
                .int => {
                    if (value == 0)
                        return .zero;

                    const UV = @Int(.unsigned, int.bitSize(value));
                    const abs_value = numeric.cast(UV, @abs(value));

                    const msb_pos_value: u16 = @typeInfo(UV).int.bits - 1 - @clz(abs_value);
                    const msb_pos_result: u16 = mantissa_bits - 1;

                    var mantissa: Mantissa = undefined;
                    var exponent: WideExponent = 0;
                    if (msb_pos_value > msb_pos_result) {
                        var shift: u16 = msb_pos_value - msb_pos_result;
                        const shift_minus_1 = shift - 1;

                        const shifted: UV = abs_value >> @intCast(shift);

                        const round_bit = (abs_value >> @intCast(shift_minus_1)) & 1;
                        const sticky = if (@ctz(abs_value) < shift_minus_1) @as(u1, 1) else 0;

                        var rounded: UV = shifted;
                        if (round_bit == 1 and (sticky == 1 or (shifted & 1) != 0))
                            rounded += 1;

                        if (comptime @typeInfo(UV).int.bits > mantissa_bits) {
                            const overflow = comptime 1 << mantissa_bits;
                            if (rounded == overflow) {
                                rounded >>= 1;
                                shift += 1;
                            }
                        }

                        mantissa = numeric.cast(Mantissa, rounded);
                        exponent +|= numeric.cast(WideExponent, shift);
                    } else {
                        const shift: u16 = msb_pos_result - msb_pos_value;
                        mantissa = numeric.cast(Mantissa, abs_value) << @intCast(shift);
                        exponent -|= numeric.cast(WideExponent, shift);
                    }

                    // Check for overflow
                    if (exponent >= int.highest(Exponent))
                        return .{
                            .mantissa = 0,
                            .exponent = int.highest(Exponent),
                            .positive = value >= 0,
                        };

                    // Check for underflow
                    if (exponent <= int.lowest(Exponent))
                        return .{
                            .mantissa = 0,
                            .exponent = int.lowest(Exponent),
                            .positive = value >= 0,
                        };

                    return .{
                        .mantissa = mantissa,
                        .exponent = numeric.cast(Exponent, exponent),
                        .positive = value >= 0,
                    };
                },
                .float => {
                    comptime if (V == f80)
                        @compileError("zsl.Dyadic(mantissa_bits, exponent_bits).initValue: f80 not yet supported");

                    const Bits = @Int(.unsigned, @typeInfo(V).float.bits);
                    const bits: Bits = @bitCast(value);

                    const f_mantissa_bits = comptime std.math.floatMantissaBits(V);
                    const f_exponent_bits = comptime std.math.floatExponentBits(V);
                    const bias = comptime (1 << (f_exponent_bits - 1)) - 1;
                    const exp_mask = comptime (1 << f_exponent_bits) - 1;
                    const mant_mask = comptime (1 << f_mantissa_bits) - 1;

                    const positive = (bits >> (f_mantissa_bits + f_exponent_bits)) == 0;
                    const biased_exp = (bits >> f_mantissa_bits) & exp_mask;
                    const frac = bits & mant_mask;

                    //NaN / Inf check
                    if (biased_exp == exp_mask) {
                        if (frac != 0)
                            return .nan;

                        return if (positive) .inf else neg(.inf);
                    }

                    // Zero check
                    if (biased_exp == 0 and frac == 0)
                        return .{
                            .mantissa = 0,
                            .exponent = int.lowest(Exponent),
                            .positive = positive,
                        };

                    const ExponentRaw = @Int(.signed, int.max(@typeInfo(WideExponent).int.bits, 32));

                    var raw_m: Bits = undefined;
                    var raw_e: ExponentRaw = undefined;
                    var msb_pos_value: u16 = undefined;

                    if (biased_exp == 0) {
                        // Subnormal: non-normalized
                        raw_m = frac;
                        raw_e = 1 - bias - f_mantissa_bits;
                        msb_pos_value = @typeInfo(Bits).int.bits - 1 - @clz(raw_m);
                    } else {
                        // Normal: implicitly normalized
                        raw_m = frac | (@as(Bits, 1) << f_mantissa_bits);
                        raw_e = numeric.cast(ExponentRaw, biased_exp) - bias - f_mantissa_bits;
                        msb_pos_value = f_mantissa_bits;
                    }

                    const msb_pos_result: u16 = mantissa_bits - 1;
                    var mantissa: Mantissa = undefined;
                    var exponent: ExponentRaw = raw_e;

                    if (msb_pos_value > msb_pos_result) {
                        var shift: u16 = msb_pos_value - msb_pos_result;
                        const shift_minus_1 = shift - 1;

                        const shifted: Bits = raw_m >> @intCast(shift);

                        const round_bit = (raw_m >> @intCast(shift_minus_1)) & 1;
                        const sticky = if (@ctz(raw_m) < shift_minus_1) @as(u1, 1) else 0;

                        var rounded: Bits = shifted;
                        if (round_bit == 1 and (sticky == 1 or (shifted & 1) != 0))
                            rounded += 1;

                        if (comptime @typeInfo(Bits).int.bits > mantissa_bits) {
                            const overflow = comptime 1 << mantissa_bits;
                            if (rounded == overflow) {
                                rounded >>= 1;
                                shift += 1;
                            }
                        }

                        mantissa = numeric.cast(Mantissa, rounded);
                        exponent += shift;
                    } else if (msb_pos_value < msb_pos_result) {
                        const shift: u16 = msb_pos_result - msb_pos_value;
                        const ShiftM = std.math.Log2Int(Mantissa);
                        mantissa = numeric.cast(Mantissa, raw_m) << @as(ShiftM, @intCast(shift));
                        exponent -= shift;
                    } else {
                        mantissa = numeric.cast(Mantissa, raw_m);
                    }

                    // Check for overflow
                    if (exponent >= int.highest(Exponent))
                        return .{
                            .mantissa = 0,
                            .exponent = int.highest(Exponent),
                            .positive = positive,
                        };

                    // Check for underflow
                    if (exponent <= int.lowest(Exponent))
                        return .{
                            .mantissa = 0,
                            .exponent = int.lowest(Exponent),
                            .positive = positive,
                        };

                    return .{
                        .mantissa = mantissa,
                        .exponent = numeric.cast(Exponent, exponent),
                        .positive = positive,
                    };
                },
                .dyadic => @compileError("zsl.Dyadic(mantissa_bits, exponent_bits).initValue: dyadics not yet supported"),
                .complex => return initValue(value.re),
                .custom => return numeric.cast(Dyadic(mantissa_bits, exponent_bits), value),
            }
        }

        /// Normalizes the dyadic in place, ensuring the mantissa is
        /// left-aligned (i.e., the most significant bit is set).
        ///
        /// ## Arguments
        /// * `self` (`*Dyadic(mantissa_bits, exponent_bits)`): A pointer to the
        ///   dyadic to normalize.
        ///
        /// ## Returns
        /// `void`
        pub fn normalize(self: *Dyadic(mantissa_bits, exponent_bits)) void {
            if (self.mantissa == 0) {
                self.exponent = int.lowest(Exponent);
                return;
            }

            const lz = @clz(self.mantissa);
            self.mantissa <<= lz;
            self.exponent -|= numeric.cast(Exponent, lz);

            if (self.exponent == int.lowest(Exponent))
                self.mantissa = 0;

            return;
        }

        pub fn add(x: Dyadic(mantissa_bits, exponent_bits), y: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            // NaN check
            if (x.isNan() or y.isNan())
                return .nan;

            // Infinity check
            if (x.isInf()) {
                if (y.isInf()) {
                    if (x.positive == y.positive)
                        return x
                    else
                        return .nan;
                } else {
                    return x;
                }
            } else if (y.isInf()) {
                return y;
            }

            // Zero check
            if (x.isZero())
                return if (y.isZero())
                    .{ .mantissa = 0, .exponent = int.lowest(Exponent), .positive = x.positive or y.positive }
                else
                    y
            else if (y.isZero())
                return x;

            // Addition or subtraction
            const order_abs: std.math.Order = if (x.exponent != y.exponent)
                int.order(x.exponent, y.exponent)
            else
                int.order(x.mantissa, y.mantissa);

            if (x.positive == y.positive) {
                var result: Dyadic(mantissa_bits, exponent_bits) =
                    if (order_abs == .gt)
                        _addAbs(x, y)
                    else
                        _addAbs(y, x);
                result.positive = x.positive;
                return result;
            }

            if (order_abs == .eq)
                return .zero;

            var result: Dyadic(mantissa_bits, exponent_bits) =
                if (order_abs == .gt)
                    _subAbs(x, y)
                else
                    _subAbs(y, x);
            result.positive = if (order_abs == .gt) x.positive else y.positive;
            return result;
        }

        pub fn sub(x: Dyadic(mantissa_bits, exponent_bits), y: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            return x.add(y.neg());
        }

        fn _addAbs(x: Dyadic(mantissa_bits, exponent_bits), y: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            // |x| >= |y|, so exponent difference is non-negative
            const exp_diff: WideExponent = numeric.cast(WideExponent, x.exponent) - numeric.cast(WideExponent, y.exponent);

            const x_wide: WideMantissa = numeric.cast(WideMantissa, x.mantissa) << @intCast(mantissa_bits);
            const y_wide: WideMantissa = numeric.cast(WideMantissa, y.mantissa) << @intCast(mantissa_bits);

            var y_shifted: WideMantissa = undefined;
            var sticky: u1 = 0;
            if (exp_diff >= 2 * mantissa_bits) {
                y_shifted = 0;
                sticky = 1;
            } else {
                y_shifted = y_wide >> @intCast(exp_diff);

                if (@ctz(y_wide) < exp_diff)
                    sticky = 1;
            }

            const sum_ov = @addWithOverflow(x_wide, y_shifted);
            var sum: WideMantissa = sum_ov[0];
            const carry: u1 = sum_ov[1];
            var exponent: Exponent = x.exponent;
            if (carry != 0) {
                if ((sum & 1) != 0)
                    sticky = 1;

                sum >>= 1;
                sum |= (@as(WideMantissa, 1) << (mantissa_bits * 2 - 1));
                exponent +|= 1;
            }

            const remainder_mask = comptime ((1 << mantissa_bits) - 1);
            const halfway_mask = comptime (1 << (mantissa_bits - 1));

            var mantissa: Mantissa = @truncate(sum >> @intCast(mantissa_bits));
            const remainder: WideMantissa = sum & remainder_mask;
            const round_up = (remainder > halfway_mask) or (remainder == halfway_mask and (sticky == 1 or (mantissa & 1) == 1));

            if (round_up) {
                const round = @addWithOverflow(mantissa, 1);
                mantissa = round[0];

                if (round[1] != 0) {
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    exponent +|= 1;
                }
            }

            // Check for overflow
            if (exponent == int.highest(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = true,
                };

            return .{
                .mantissa = mantissa,
                .exponent = exponent,
                .positive = true,
            };
        }

        fn _subAbs(x: Dyadic(mantissa_bits, exponent_bits), y: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            // |x| >= |y|, so exponent difference is non-negative
            const exp_diff: WideExponent = numeric.cast(WideExponent, x.exponent) - numeric.cast(WideExponent, y.exponent);

            const x_wide: WideMantissa = numeric.cast(WideMantissa, x.mantissa) << @intCast(mantissa_bits);
            const y_wide: WideMantissa = numeric.cast(WideMantissa, y.mantissa) << @intCast(mantissa_bits);

            var y_shifted: WideMantissa = undefined;
            var sticky: u1 = 0;
            if (exp_diff >= 2 * mantissa_bits) {
                y_shifted = 0;
                sticky = 1;
            } else {
                y_shifted = y_wide >> @intCast(exp_diff);

                if (@ctz(y_wide) < exp_diff)
                    sticky = 1;
            }

            var diff: WideMantissa = x_wide - y_shifted;
            if (sticky == 1)
                diff -= 1;

            if (diff == 0)
                return .zero;

            const lz = @clz(diff);

            var exponent: WideExponent = numeric.cast(WideExponent, x.exponent) - numeric.cast(WideExponent, lz);

            const remainder_mask = comptime ((1 << mantissa_bits) - 1);
            const halfway_mask = comptime (1 << (mantissa_bits - 1));

            const shifted_diff = diff << @intCast(lz);

            var mantissa: Mantissa = @truncate(shifted_diff >> @intCast(mantissa_bits));
            const remainder: WideMantissa = shifted_diff & remainder_mask;
            const round_up = (remainder > halfway_mask) or (remainder == halfway_mask and (sticky == 1 or (mantissa & 1) == 1));

            if (round_up) {
                const round = @addWithOverflow(mantissa, 1);
                mantissa = round[0];

                if (round[1] != 0) {
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    exponent +|= 1;
                }
            }

            // Check for underflow
            if (exponent <= int.lowest(Exponent))
                return .zero;

            return .{
                .mantissa = mantissa,
                .exponent = numeric.cast(Exponent, exponent),
                .positive = true,
            };
        }

        pub fn mul(x: Dyadic(mantissa_bits, exponent_bits), y: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            // NaN check
            if (x.isNan() or y.isNan())
                return .nan;

            const result_positive = x.positive == y.positive;

            // Infinity check
            if (x.isInf()) {
                if (y.isZero())
                    return .nan; // Inf * 0 = NaN

                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = result_positive,
                };
            } else if (y.isInf()) {
                if (x.isZero())
                    return .nan; // 0 * Inf = NaN

                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = result_positive,
                };
            }

            // Zero check
            if (x.isZero() or y.isZero())
                return .{
                    .mantissa = 0,
                    .exponent = int.lowest(Exponent),
                    .positive = result_positive,
                };

            // Multiplication
            var product: WideMantissa = numeric.cast(WideMantissa, x.mantissa) * numeric.cast(WideMantissa, y.mantissa);

            // Base exponent assumes the MSB landed at 2 * mantissa_bits - 1
            var exponent: WideExponent = numeric.cast(WideExponent, x.exponent) + numeric.cast(WideExponent, y.exponent) + numeric.cast(WideExponent, mantissa_bits);

            // Normalize: If the MSB is at 2 * mantissa_bits - 2, shift left to standardize it
            if ((product >> (mantissa_bits * 2 - 1)) == 0) {
                product <<= 1;
                exponent -= 1;
            }

            const halfway_mask = comptime (1 << (mantissa_bits - 1));
            const remainder_mask = comptime ((1 << mantissa_bits) - 1);

            var mantissa: Mantissa = @truncate(product >> @intCast(mantissa_bits));
            const remainder: WideMantissa = product & remainder_mask;

            // Round to nearest, tie to even
            if (remainder > halfway_mask or (remainder == halfway_mask and (mantissa & 1) == 1)) {
                const round = @addWithOverflow(mantissa, 1);
                mantissa = round[0];

                if (round[1] != 0) {
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    exponent += 1;
                }
            }

            // Check for overflow
            if (exponent >= int.highest(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = result_positive,
                };

            // Check for underflow
            if (exponent <= int.lowest(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.lowest(Exponent),
                    .positive = result_positive,
                };

            return .{
                .mantissa = mantissa,
                .exponent = numeric.cast(Exponent, exponent),
                .positive = result_positive,
            };
        }

        pub fn fma(x: Dyadic(mantissa_bits, exponent_bits), y: Dyadic(mantissa_bits, exponent_bits), z: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            // NaN check
            if (x.isNan() or y.isNan() or z.isNan())
                return .nan;

            const xy_positive = x.positive == y.positive;

            // Infinity check
            if (x.isInf() or y.isInf()) {
                if ((x.isInf() and y.isZero()) or (x.isZero() and y.isInf()))
                    return .nan; // Inf * 0 or 0 * Inf = NaN regardless of z

                if (z.isInf() and z.positive != xy_positive)
                    return .nan; // (±Inf) + (∓Inf) = NaN

                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = xy_positive,
                };
            } else if (z.isInf()) {
                return z;
            }

            // Zero check
            if (x.isZero() or y.isZero()) {
                if (z.isZero())
                    return .{
                        .mantissa = 0,
                        .exponent = int.lowest(Exponent),
                        .positive = xy_positive or z.positive,
                    };

                return z;
            }

            if (z.isZero())
                return x.mul(y);

            // Fused multiply-add (single rounding)
            var p_mant: WideMantissa = numeric.cast(WideMantissa, x.mantissa) * numeric.cast(WideMantissa, y.mantissa);
            var p_exp: WideExponent = numeric.cast(WideExponent, x.exponent) + numeric.cast(WideExponent, y.exponent);

            // Check for MSB position
            if ((p_mant >> (2 * mantissa_bits - 1)) == 0) {
                p_mant <<= 1;
                p_exp -= 1;
            }

            const z_mant: WideMantissa = numeric.cast(WideMantissa, z.mantissa) << @intCast(mantissa_bits);
            const z_exp: WideExponent = numeric.cast(WideExponent, z.exponent) - numeric.cast(WideExponent, mantissa_bits);

            // Order by magnitude
            const p_larger = if (p_exp != z_exp) p_exp > z_exp else p_mant >= z_mant;
            const larger_mant: WideMantissa = if (p_larger) p_mant else z_mant;
            const larger_exp: WideExponent = if (p_larger) p_exp else z_exp;
            const larger_pos: bool = if (p_larger) xy_positive else z.positive;
            const smaller_mant: WideMantissa = if (p_larger) z_mant else p_mant;
            const smaller_exp: WideExponent = if (p_larger) z_exp else p_exp;
            const smaller_pos: bool = if (p_larger) z.positive else xy_positive;

            // Align smaller to larger
            const exp_diff: WideExponent = larger_exp - smaller_exp;
            var smaller_shifted: WideMantissa = smaller_mant;

            // frac holds the exact bits shifted out of smaller_mant, aligned to the MSB
            var frac: WideMantissa = 0;
            var sticky: u1 = 0;

            if (exp_diff >= 2 * mantissa_bits) {
                smaller_shifted = 0;
                sticky = 1; // smaller_mant is strictly nonzero here
            } else if (exp_diff > 0) {
                const shift: std.math.Log2Int(WideMantissa) = @intCast(exp_diff);
                smaller_shifted = smaller_mant >> shift;

                const shifted_out = smaller_mant & ((@as(WideMantissa, 1) << shift) - 1);
                if (shifted_out != 0) {
                    frac = shifted_out << @intCast(2 * mantissa_bits - shift);
                }
            }

            // Add (same sign) or subtract (opposite sign)
            var result_mant: WideMantissa = undefined;
            var result_exp: WideExponent = larger_exp;

            if (larger_pos == smaller_pos) {
                const sum_ov = @addWithOverflow(larger_mant, smaller_shifted);
                result_mant = sum_ov[0];
                if (sum_ov[1] != 0) {
                    // Shift right to accommodate carry, pushing LSB into frac
                    sticky |= @intCast(frac & 1);
                    frac = (frac >> 1) | (@as(WideMantissa, @intCast(result_mant & 1)) << @intCast(2 * mantissa_bits - 1));
                    result_mant = (result_mant >> 1) | (@as(WideMantissa, 1) << @intCast(2 * mantissa_bits - 1));
                    result_exp +|= 1;
                }
            } else {
                var diff: WideMantissa = larger_mant - smaller_shifted;

                if (frac != 0 or sticky == 1) {
                    diff -= 1; // Borrow 1 from the integer part

                    // Compute the new fractional part
                    frac = if (sticky == 1) ~frac else ~frac +% 1;
                }

                if (diff == 0 and frac == 0) {
                    return .{
                        .mantissa = 0,
                        .exponent = int.lowest(Exponent),
                        .positive = true, // Exact zero is +0.0
                    };
                }

                // Renormalize
                if (diff == 0) {
                    // Massive cancellation: the result lives entirely in the fractional part
                    const lz = @clz(frac);
                    result_mant = frac << @intCast(lz);
                    result_exp -|= numeric.cast(WideExponent, lz + 1);
                    frac = 0;
                } else {
                    const lz = @clz(diff);
                    diff <<= @intCast(lz);
                    if (lz > 0) {
                        // Pull the top lz bits of frac up into the bottom of diff
                        diff |= frac >> @intCast(2 * mantissa_bits - lz);
                        frac <<= @intCast(lz);
                    }

                    result_mant = diff;
                    result_exp -|= numeric.cast(WideExponent, lz);
                }
            }

            // Restore sticky
            if (frac != 0)
                sticky = 1;

            const remainder_mask = comptime ((1 << mantissa_bits) - 1);
            const halfway_mask = comptime (1 << (mantissa_bits - 1));

            var mantissa: Mantissa = @truncate(result_mant >> @intCast(mantissa_bits));
            const remainder: WideMantissa = result_mant & remainder_mask;
            var final_exp: WideExponent = result_exp +| numeric.cast(WideExponent, mantissa_bits);
            const round_up = (remainder > halfway_mask) or (remainder == halfway_mask and (sticky == 1 or (mantissa & 1) == 1));

            if (round_up) {
                const inc = @addWithOverflow(mantissa, 1);
                mantissa = inc[0];
                if (inc[1] != 0) {
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    final_exp +|= 1;
                }
            }

            // Check for overflow
            if (final_exp >= int.highest(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = larger_pos,
                };

            // Check for underflow
            if (final_exp <= int.lowest(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.lowest(Exponent),
                    .positive = larger_pos,
                };

            return .{
                .mantissa = mantissa,
                .exponent = numeric.cast(Exponent, final_exp),
                .positive = larger_pos,
            };
        }

        pub fn div(x: Dyadic(mantissa_bits, exponent_bits), y: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            // NaN check
            if (x.isNan() or y.isNan())
                return .nan;

            const result_positive = x.positive == y.positive;

            // Infinity check
            if (x.isInf()) {
                if (y.isInf())
                    return .nan;

                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = result_positive,
                };
            } else if (y.isInf()) {
                return .{
                    .mantissa = 0,
                    .exponent = int.lowest(Exponent),
                    .positive = result_positive,
                };
            }

            // Zero check
            if (y.isZero()) {
                if (x.isZero())
                    return .nan; // 0 / 0 = NaN

                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = result_positive,
                };
            } else if (x.isZero()) {
                return .{
                    .mantissa = 0,
                    .exponent = int.lowest(Exponent),
                    .positive = result_positive,
                };
            }

            // Division
            const x_wide: WideMantissa = numeric.cast(WideMantissa, x.mantissa) << @intCast(mantissa_bits);
            const y_wide: WideMantissa = numeric.cast(WideMantissa, y.mantissa);

            const quotient: WideMantissa = x_wide / y_wide;
            const remainder: WideMantissa = x_wide % y_wide;

            var exponent: WideExponent = numeric.cast(WideExponent, x.exponent) - numeric.cast(WideExponent, y.exponent) - numeric.cast(WideExponent, mantissa_bits);
            var mantissa: Mantissa = undefined;
            var round_up = false;

            const q_msb_mask = comptime (1 << mantissa_bits);

            if ((quotient & q_msb_mask) != 0) {
                // MSB is at mantissa_bits. We must shift right by 1 to normalize
                mantissa = @truncate(quotient >> 1);
                exponent +|= 1;

                // Dropped bit is 1 AND (sticky remainder exists OR tie-to-even)
                round_up = (quotient & 1) == 1 and (remainder > 0 or (mantissa & 1) == 1);
            } else {
                // MSB is at mantissa_bits - 1. Fits perfectly
                mantissa = @truncate(quotient);
                const rem_doubled = remainder << 1;

                // > 0.5 OR (exactly 0.5 AND tie-to-even)
                round_up = (rem_doubled > y_wide) or (rem_doubled == y_wide and (mantissa & 1) == 1);
            }

            if (round_up) {
                const round = @addWithOverflow(mantissa, 1);
                mantissa = round[0];

                if (round[1] != 0) {
                    // Mantissa overflowed during rounding (e.g., 1111 + 1 = 10000)
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    exponent +|= 1;
                }
            }

            // Check for overflow
            if (exponent >= int.highest(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.highest(Exponent),
                    .positive = result_positive,
                };

            // Check for underflow
            if (exponent <= int.lowest(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.lowest(Exponent),
                    .positive = result_positive,
                };

            return .{
                .mantissa = mantissa,
                .exponent = numeric.cast(Exponent, exponent),
                .positive = result_positive,
            };
        }

        pub fn order(x: Dyadic(mantissa_bits, exponent_bits), y: Dyadic(mantissa_bits, exponent_bits)) std.math.Order {
            // NaN check
            if (x.isNan())
                return if (y.isNan()) .eq else .gt;

            if (y.isNan()) return .lt;

            // Zero check
            if (x.isZero() and y.isZero())
                return .eq;

            // Sign check
            if (x.positive != y.positive)
                return if (x.positive) .gt else .lt;

            // Signs equal, evaluate magnitude
            var mag_order: std.math.Order = .eq;

            if (x.isInf()) {
                mag_order = if (y.isInf()) .eq else .gt;
            } else if (y.isInf()) {
                mag_order = .lt;
            } else {
                mag_order = if (x.exponent != y.exponent)
                    int.order(x.exponent, y.exponent)
                else
                    int.order(x.mantissa, y.mantissa);
            }

            return if (x.positive) mag_order else mag_order.invert();
        }

        pub fn abs(self: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            return .{
                .mantissa = self.mantissa,
                .exponent = self.exponent,
                .positive = true,
            };
        }

        pub fn neg(self: Dyadic(mantissa_bits, exponent_bits)) Dyadic(mantissa_bits, exponent_bits) {
            return .{
                .mantissa = self.mantissa,
                .exponent = self.exponent,
                .positive = !self.positive,
            };
        }

        pub fn toFloat(self: Dyadic(mantissa_bits, exponent_bits), comptime Float: type) Float {
            comptime if (!meta.isNumeric(Float) or meta.numericType(Float) != .float)
                @compileError("zsl.Dyadic(mantissa_bits, exponent_bits).toFloat: Float must be a float type, got \n\tFloat = " ++ @typeName(Float) ++ "\n");

            comptime if (Float == f80)
                @compileError("zsl.Dyadic(mantissa_bits, exponent_bits).toFloat: f80 not yet supported");

            if (self.isNan())
                return std.math.nan(Float);

            if (self.isInf())
                return if (self.positive) std.math.inf(Float) else -std.math.inf(Float);

            if (self.isZero())
                return if (self.positive) 0.0 else -0.0;

            const f_mantissa_bits = comptime std.math.floatMantissaBits(Float);
            const f_exponent_bits = comptime std.math.floatExponentBits(Float);
            const bias = comptime (1 << (f_exponent_bits - 1)) - 1;
            const max_biased = comptime (1 << f_exponent_bits) - 1;
            const frac_mask = comptime ((1 << f_mantissa_bits) - 1);

            const Bits = @Int(.unsigned, @typeInfo(Float).float.bits);
            const ExponentRaw = @Int(.signed, int.max(@typeInfo(WideExponent).int.bits, 32));

            const raw_e: ExponentRaw = numeric.cast(ExponentRaw, self.exponent) +| (mantissa_bits - 1);
            var biased_exp: ExponentRaw = raw_e +| bias;

            var right_shift: ExponentRaw = (numeric.cast(ExponentRaw, mantissa_bits) - 1) - f_mantissa_bits;

            if (biased_exp >= max_biased)
                return if (self.positive) std.math.inf(Float) else -std.math.inf(Float);

            if (biased_exp <= 0) {
                // Handle subnormals: exponent gets pinned to 0, and we shift right by the deficit
                right_shift += (1 - biased_exp);
                biased_exp = 0;
            }

            const Wide = @Int(.unsigned, int.max(mantissa_bits, f_mantissa_bits + 2));
            var m: Wide = undefined;

            if (right_shift < 0) {
                // Left shift, no rounding possible
                m = numeric.cast(Wide, self.mantissa) << @intCast(-right_shift);
            } else if (right_shift == 0) {
                // Exact alignment
                m = numeric.cast(Wide, self.mantissa);
            } else {
                const rs: u16 = numeric.cast(u16, right_shift);
                if (rs > mantissa_bits) {
                    // Completely shifted out of bounds
                    m = 0;
                } else if (rs == mantissa_bits) {
                    // Round bit is the MSB (which is 1). Sticky is the rest
                    const sticky_mask = comptime (1 << (mantissa_bits - 1)) - 1;
                    m = if ((self.mantissa & sticky_mask) != 0) 1 else 0;
                } else {
                    const shift_minus_1 = rs - 1;
                    const shifted = self.mantissa >> @intCast(rs);

                    const round_bit = (self.mantissa >> @intCast(shift_minus_1)) & 1;
                    const sticky = if (@ctz(self.mantissa) < shift_minus_1) @as(u1, 1) else 0;

                    m = numeric.cast(Wide, shifted);
                    if (round_bit == 1 and (sticky == 1 or (shifted & 1) != 0))
                        m += 1;
                }
            }

            // Renormalization
            if (biased_exp > 0) {
                const normal_threshold = comptime 1 << (f_mantissa_bits + 1);
                if (m >= normal_threshold) {
                    m >>= 1;
                    biased_exp += 1;
                    if (biased_exp >= max_biased)
                        return if (self.positive) std.math.inf(Float) else -std.math.inf(Float);
                }
            } else {
                const subnormal_threshold = comptime 1 << f_mantissa_bits;
                if (m >= subnormal_threshold) {
                    // Subnormal rounded up completely into a normal number
                    biased_exp = 1;
                }
            }

            const frac: Bits = numeric.cast(Bits, m & frac_mask);
            const exp_field: Bits = numeric.cast(Bits, biased_exp);
            const sign_bit: Bits = if (self.positive) 0 else 1;

            const float_bits: Bits = (sign_bit << (f_mantissa_bits + f_exponent_bits)) | (exp_field << f_mantissa_bits) | frac;

            return @bitCast(float_bits);
        }
    };
}

pub const Coerce = @import("dyadic/coerce.zig").Coerce;

pub const pi = @import("dyadic/pi.zig").pi;
pub const tau = @import("dyadic/tau.zig").tau;
pub const e = @import("dyadic/e.zig").e;

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.Add: at least one of X or Y must be a dyadic type, the other must be a bool, an int, a float or a dyadic type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return dyadic.Coerce(X, Y);
}

/// Performs addition between two operands of dyadic, float, int or bool types,
/// where at least one operand must be of dyadic type. The result type is
/// determined by coercing the operand types, and the operation is performed by
/// casting both operands to the result type, then adding them.
///
/// ## Signature
/// ```zig
/// dyadic.add(x: X, y: Y) dyadic.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `dyadic.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
///
/// ## Method
/// The algorithm implements exact floating-point addition and subtraction of
/// `x: Dyadic(N, E)` and `y: Dyadic(N, E)`, ensuring exact rounding and
/// graceful handling of catastrophic cancellation [1]:
/// 1. Magnitude Ordering: The operands are evaluated by absolute magnitude. The
///    operation is delegated to either absolute addition (`_addAbs`) or
///    absolute subtraction (`_subAbs`) such that the larger magnitude operand
///    is always processed first.
/// 2. Alignment & Sticky Bit: To add or subtract mantissas, their exponents
///    must align. Both mantissas are upcast to a `WideMantissa` representation
///    (shifted left by `N` bits) to act as guard digits. The smaller operand's
///    wide mantissa is right-shifted by the exponent difference. Any non-zero
///    bits shifted entirely out of the representable range are logically OR'd
///    into a single sticky bit. This bit acts as a mathematical proxy for all
///    lost precision, strictly required for exact rounding [2].
/// 3. Execution & Renormalization:
///    * Addition: The aligned mantissas are added. If the sum overflows
///      (produces a carry bit), the result is right-shifted by one, the sticky
///      bit is updated, the most significant bit is restored, and the exponent
///      is incremented.
///    * Subtraction: The smaller aligned mantissa is subtracted from the
///      larger. If the sticky bit is set, a `1` is borrowed (subtracted) from
///      the difference. Because subtraction of similar magnitudes can result in
///      massive precision loss (catastrophic cancellation), the result is
///      rapidly renormalized by counting leading zeros, shifting the mantissa
///      left, and decrementing the exponent accordingly [1].
/// 4. Rounding: The result is rounded using the standard "round to nearest,
///    ties to even" mode. This is achieved by inspecting the remainder against
///    a halfway mask alongside the sticky bit. If rounding causes the mantissa
///    to overflow, a final post-normalization shift is applied [1].
///
/// ## References
/// [1] Muller, J.-M., et al. (2018). Handbook of Floating-Point Arithmetic (2nd ed.). Birkhäuser. https://doi.org/10.1007/978-3-319-76526-6
/// [2] Goldberg, D. (1991). What Every Computer Scientist Should Know About Floating-Point Arithmetic. ACM Computing Surveys, 23(1), 5-48. https://doi.org/10.1145/103162.103163
pub fn add(x: anytype, y: anytype) dyadic.Add(@TypeOf(x), @TypeOf(y)) {
    const R: type = Add(@TypeOf(x), @TypeOf(y));

    return R.add(numeric.cast(R, x), numeric.cast(R, y));
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.Sub: at least one of X or Y must be a dyadic type, the other must be a bool, an int, a float or a dyadic type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return dyadic.Coerce(X, Y);
}

/// Performs subtraction between two operands of dyadic, float, int or bool
/// types, where at least one operand must be of dyadic type. The result type is
/// determined by coercing the operand types, and the operation is performed by
/// casting both operands to the result type, then subtracting them.
///
/// ## Signature
/// ```zig
/// dyadic.sub(x: X, y: Y) dyadic.Sub(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `dyadic.Sub(@TypeOf(x), @TypeOf(y))`: The result of the subtraction.
///
/// ## Method
/// The algorithm implements exact floating-point subtraction of
/// `x: Dyadic(N, E)` and `y: Dyadic(N, E)` by algebraically reducing it to
/// addition: `x - y = x + (-y)`. Because floating-point negation is an exact,
/// bitwise operation (simply flipping the sign flag) that incurs no loss of
/// precision, this reduction is mathematically lossless [1].
///
/// ## References
/// [1] Muller, J.-M., et al. (2018). Handbook of Floating-Point Arithmetic (2nd ed.). Birkhäuser. https://doi.org/10.1007/978-3-319-76526-6
pub fn sub(x: anytype, y: anytype) Sub(@TypeOf(x), @TypeOf(y)) {
    const R: type = Sub(@TypeOf(x), @TypeOf(y));

    return R.sub(numeric.cast(R, x), numeric.cast(R, y));
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.Mul: at least one of X or Y must be a dyadic type, the other must be a bool, an int, a float or a dyadic type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return dyadic.Coerce(X, Y);
}

/// Performs multiplication between two operands of dyadic, float, int or bool
/// types, where at least one operand must be of dyadic type. The result type is
/// determined by coercing the operand types, and the operation is performed by
/// casting both operands to the result type, then multiplying them.
///
/// ## Signature
/// ```zig
/// dyadic.mul(x: X, y: Y) dyadic.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `dyadic.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
///
/// ## Method
/// The algorithm implements exact floating-point multiplication of
/// `x: Dyadic(N, E)` and `y: Dyadic(N, E)` by computing the full-width product
/// of the integer mantissas and strictly normalizing [1]:
/// 1. Full-Width Multiplication: The `N`-bit normalized mantissas are
///    multiplied to form a complete, exact `2N`-bit wide product. Because the
///    input mantissas mathematically fall within the normalized bounds of
///    `[2ᴺ⁻¹, 2ᴺ)`, their algebraic product is strictly bounded within
///    `[2²ᴺ⁻², 2²ᴺ)`. The sign of the result is determined independently.
/// 2. Renormalization: Due to the bounds established in step 1, the most
///    significant bit (MSB) of the product is mathematically guaranteed to land
///    at either index `2N - 1` or `2N - 2`. If it lands at `2N - 2`, a `1`-bit
///    left shift renormalizes the mantissa, and the working exponent is
///    decremented by `1`. No iterative shifting or leading-zero counting is
///    required [1].
/// 3. Exact Rounding: Because the operation computes the exact `2N`-bit
///    product, no precision is lost during the calculation. The lower `N` bits
///    act as an exact remainder. This remainder is checked against a halfway
///    mask to apply standard "round to nearest, ties to even" logic before
///    truncation.
///
/// ## References
/// [1] Muller, J.-M., et al. (2018). Handbook of Floating-Point Arithmetic (2nd ed.). Birkhäuser. https://doi.org/10.1007/978-3-319-76526-6
pub fn mul(x: anytype, y: anytype) Mul(@TypeOf(x), @TypeOf(y)) {
    const R: type = Mul(@TypeOf(x), @TypeOf(y));

    return R.mul(numeric.cast(R, x), numeric.cast(R, y));
}

pub fn Fma(comptime X: type, comptime Y: type, comptime Z: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or !meta.isNumeric(Z) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or !meta.numericType(Z).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic and meta.numericType(Z) != .dyadic))
        @compileError("zsl.dyadic.Fma: at least one of X, Y or Z must be a dyadic type, the others must be bool, int, float or dyadic types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tZ = " ++ @typeName(Z) ++ "\n");

    return dyadic.Coerce(X, numeric.Coerce(Y, Z));
}

/// Performs fused multiplication and addition (x * y + z) between three
/// operands of dyadic, float, int or bool types, where at least one operand
/// must be of dyadic type. The result type is determined by coercing the
/// operand types, and the operation is performed by casting all three operands
/// to the result type, then performing the fused operation.
///
/// ## Signature
/// ```zig
/// dyadic.fma(x: X, y: Y, z: Z) dyadic.Fma(X, Y, Z)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left multiplication operand.
/// * `y` (`anytype`): The right multiplication operand.
/// * `z` (`anytype`): The addition operand.
///
/// ## Returns
/// `dyadic.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z))`: The result of the fused
/// multiplication and addition.
///
/// ## Method
/// The algorithm implements exact fused multiply-add for `x: Dyadic(N, E)`,
/// `y: Dyadic(N, E)`, and `z: Dyadic(N, E)` computing the infinitely precise
/// mathematical result with only a single final rounding step to eliminate
/// double-rounding errors [1]:
/// 1. Exact Product: The `N`-bit mantissas of `x` and `y` are multiplied to
///    form an exact `2N`-bit product. Because the input mantissas
///    mathematically fall within `[2ᴺ⁻¹, 2ᴺ)`, the unrounded algebraic product
///    is strictly bounded within `[2²ᴺ⁻², 2²ᴺ)`. No truncation or rounding
///    occurs here.
/// 2. Alignment & Fractional Tracking: The exact product and the addend `z` are
///    ordered by absolute magnitude. To align the exponents, the smaller
///    operand is right-shifted. Any exact bits shifted out of the primary `2N`
///    working width are perfectly preserved in a frac register. Bits shifted
///    entirely beyond `frac` are logically OR'd into a single sticky bit. This
///    architecture evaluates the true sum without requiring hardware registers
///    of infinite width [1].
/// 3. Execution & Renormalization: The aligned operands are added or
///    subtracted. If a subtraction of similar magnitudes triggers massive
///    precision loss (catastrophic cancellation), the result is rapidly
///    renormalized by counting leading zeros. The exact bits stored in the
///    `frac` register are pulled back up into the primary mantissa to fully
///    recover the mathematically exact difference [1].
/// 4. Single Exact Rounding: The single, final rounding step evaluates the
///    primary remainder alongside the `frac` and sticky bits. This applies
///    standard "round to nearest, ties to even" logic to the exact algebraic
///    outcome before standardizing the exponent.
///
/// ## References
/// [1] Muller, J.-M., et al. (2018). Handbook of Floating-Point Arithmetic (2nd ed.). Birkhäuser. https://doi.org/10.1007/978-3-319-76526-6
pub fn fma(x: anytype, y: anytype, z: anytype) dyadic.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    const R: type = dyadic.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z));

    return R.fma(numeric.cast(R, x), numeric.cast(R, y), numeric.cast(R, z));
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.Div: at least one of X or Y must be a dyadic type, the other must be a bool, an int, a float or a dyadic type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return dyadic.Coerce(X, Y);
}

/// Performs division between two operands of dyadic, float, int or bool types,
/// where at least one operand must be of dyadic type. The result type is
/// determined by coercing the operand types, and the operation is performed by
/// casting both operands to the result type, then dividing them.
///
/// ## Signature
/// ```zig
/// dyadic.div(x: X, y: Y) dyadic.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `dyadic.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
///
/// ## Method
/// The algorithm implements exact floating-point division of
/// `x: Dyadic(N, E)` and `y: Dyadic(N, E)` by leveraging exact integer division
/// of the scaled mantissas [1]:
/// 1. Scaled Integer Division: The dividend mantissa `x` is shifted left by `N`
///    bits, producing a `2N`-bit wide integer. This is divided by the `N`-bit
///    divisor mantissa `y` to yield an exact integer quotient and remainder.
///    Because the normalized inputs mathematically fall within the bounds of
///    `[2ᴺ⁻¹, 2ᴺ)`, the resulting quotient is strictly bounded within
///    `[2ᴺ⁻¹, 2ᴺ⁺¹)`.
/// 2. Reormalization: Due to the bounds established in step 1, the most
///    significant bit (MSB) of the quotient is mathematically guaranteed to
///    land at either index `N` or `N - 1`. If it lands at `N`, a `1`-bit right
///    shift renormalizes the mantissa, and the working exponent is incremented
///    by `1`.
/// 3. Exact Rounding: The exact mathematical remainder is used to apply
///    standard "round to nearest, ties to even" logic without any loss of
///    floating-point precision. If the quotient required a right shift, the
///    dropped bit and the presence of a non-zero remainder (acting as a sticky
///    bit) dictate the rounding. If no shift was required, the remainder is
///    algebraically doubled and compared directly against the divisor to
///    determine if the fractional part exceeds exactly `0.5` [1].
///
/// ## References
/// [1] Muller, J.-M., et al. (2018). Handbook of Floating-Point Arithmetic (2nd ed.). Birkhäuser. https://doi.org/10.1007/978-3-319-76526-6
pub fn div(x: anytype, y: anytype) Div(@TypeOf(x), @TypeOf(y)) {
    const R: type = Div(@TypeOf(x), @TypeOf(y));

    return numeric.cast(R, x).div(numeric.cast(R, y));
}

/// Compares two operands of dyadic, float, int or bool types, where at least
/// one operand must be of dyadic type, for ordering. The operation is performed
/// by casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// dyadic.order(x: X, y: Y) std.math.Order
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `std.math.Order`: The result of the comparison.
pub fn order(x: anytype, y: anytype) std.math.Order {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.order: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = dyadic.Coerce(X, Y);

    return C.order(numeric.cast(C, x), numeric.cast(C, y));
}

/// Compares two operands of dyadic, float, int or bool types, where at least
/// one operand must be of dyadic type, for equality. The operation is performed
/// by casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// dyadic.eq(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are equal, `false` otherwise.
pub fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.eq: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return dyadic.order(x, y) == .eq;
}

/// Compares two operands of dyadic, float, int or bool types, where at least
/// one operand must be of dyadic type, for inequality. The operation is
/// performed by casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// dyadic.ne(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are not equal, `false` otherwise.
pub fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.ne: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return dyadic.order(x, y) != .eq;
}

/// Compares two operands of dyadic, float, int or bool types, where at least
/// one operand must be of dyadic type, for for less-than ordering. The
/// operation is performed by casting both operands to the coerced type, then
/// comparing them.
///
/// ## Signature
/// ```zig
/// dyadic.lt(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is less than `y`, `false` otherwise.
pub fn lt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.lt: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return dyadic.order(x, y) == .lt;
}

/// Compares two operands of dyadic, float, int or bool types, where at least
/// one operand must be of dyadic type, for less-than or equal ordering. The
/// operation is performed by casting both operands to the coerced type, then
/// comparing them.
///
/// ## Signature
/// ```zig
/// dyadic.le(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is less than or equal to `y`, `false` otherwise.
pub fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.le: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const order_xy = dyadic.order(x, y);
    return order_xy == .eq or order_xy == .lt;
}

/// Compares two operands of dyadic, float, int or bool types, where at least
/// one operand must be of dyadic type, for greater-than ordering. The operation
/// is performed by casting both operands to the coerced type, then comparing
/// them.
///
/// ## Signature
/// ```zig
/// dyadic.gt(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is greater than `y`, `false` otherwise.
pub fn gt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.gt: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return dyadic.order(x, y) == .gt;
}

/// Compares two operands of dyadic, float, int or bool types, where at least
/// one operand must be of dyadic type, for greater-than or equality ordering.
/// The operation is performed by casting both operands to the coerced type,
/// then comparing them.
///
/// ## Signature
/// ```zig
/// dyadic.ge(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is greater than or equal to `y`, `false` otherwise.
pub fn ge(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.dyadic) or !meta.numericType(Y).le(.dyadic) or
        (meta.numericType(X) != .dyadic and meta.numericType(Y) != .dyadic))
        @compileError("zsl.dyadic.le: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const order_xy = dyadic.order(x, y);
    return order_xy == .eq or order_xy == .gt;
}

pub fn isNan(x: anytype) bool {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or meta.numericType(X) != .dyadic)
        @compileError("zsl.dyadic.isNan: x must be a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n");

    return x.isNan();
}

pub fn abs(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or meta.numericType(X) != .dyadic)
        @compileError("zsl.dyadic.abs: x must be a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n");

    return x.abs();
}

pub const sign = @import("dyadic/sign.zig").sign;

pub const sqrt = @import("dyadic/sqrt.zig").sqrt;
