//! Namespace for dyadic operations.

const std = @import("std");

const types = @import("types.zig");
const Cmp = types.Cmp;
const int = @import("int.zig");
const ops = @import("ops.zig");

/// Arbitrary-precision (chosen at compile time)
/// dyadic type.
pub fn Dyadic(mantissa_bits: u16, exponent_bits: u16) type {
    if (mantissa_bits == 0 or exponent_bits == 0 or
        mantissa_bits >= int.maxVal(u16) / 2 or exponent_bits >= int.maxVal(u16) / 2)
        @compileError(std.fmt.comptimePrint(
            "zml.Dyadic: both mantissa_bits and exponent_bits must be non-zero and less than {}, got\n\tmantissa_bits: {}\n\texponent_bits: {}\n",
            .{ int.maxVal(u16) / 2, mantissa_bits, exponent_bits },
        ));

    return packed struct {
        mantissa: Mantissa,
        exponent: Exponent,
        positive: bool,

        /// Type signature
        pub const is_numeric = true;
        pub const is_dyadic = true;
        pub const is_real_type = true;
        pub const is_signed = true;

        pub const Mantissa = std.meta.Int(.unsigned, mantissa_bits);
        pub const Exponent = std.meta.Int(.signed, exponent_bits);

        pub const inf: Dyadic(mantissa_bits, exponent_bits) = .{
            .mantissa = 0,
            .exponent = int.maxVal(Exponent),
            .positive = true,
        };
        pub const nan: Dyadic(mantissa_bits, exponent_bits) = .{
            .mantissa = 1,
            .exponent = int.maxVal(Exponent),
            .positive = true,
        };
        pub const zero: Dyadic(mantissa_bits, exponent_bits) = .{
            .mantissa = 0,
            .exponent = int.minVal(Exponent),
            .positive = true,
        };
        pub const one: Dyadic(mantissa_bits, exponent_bits) = .{
            .mantissa = @as(Mantissa, 1) << (mantissa_bits - 1),
            .exponent = int.minVal(Exponent) +| 1,
            .positive = true,
        };
        pub const negOne: Dyadic(mantissa_bits, exponent_bits) = .{
            .mantissa = @as(Mantissa, 1) << (mantissa_bits - 1),
            .exponent = int.minVal(Exponent) +| 1,
            .positive = false,
        };

        /// Initializes a new dyadic with the specified value.
        ///
        /// Signature
        /// ---------
        /// ```zig
        /// fn initSet(value: V) Dyadic(mantissa_bits, exponent_bits)
        /// ```
        ///
        /// Parameters
        /// ----------
        /// `value` (`anytype`):
        /// The value to set the dyadic to. Must be a numeric type or a string.
        ///
        /// Returns
        /// -------
        /// `Dyadic(mantissa_bits, exponent_bits)`:
        /// The newly initialized dyadic.
        pub fn initSet(value: anytype) Dyadic(mantissa_bits, exponent_bits) {
            _ = value;
            return .{ .mantissa = 0, .exponent = 0, .positive = true }; // TODO: implement
        }

        /// Sets the value of the dyadic in-place.
        ///
        /// Signature
        /// ---------
        /// ```zig
        /// fn set(value: V) void
        /// ```
        ///
        /// Parameters
        /// ----------
        /// `self` (`*Dyadic(mantissa_bits, exponent_bits)`):
        /// A pointer to the dyadic to set.
        ///
        /// `value` (`anytype`):
        /// The value to set the dyadic to. Must be a numeric type or a string.
        ///
        /// Returns
        /// -------
        /// `void`
        pub fn set(self: *Dyadic(mantissa_bits, exponent_bits), value: anytype) void {
            // TODO: implement
            _ = value;
            self.* = .{ .mantissa = 0, .exponent = 0, .positive = true };
        }

        pub fn isInf(self: Dyadic(mantissa_bits, exponent_bits)) bool {
            return self.exponent == int.maxVal(Exponent) and self.mantissa == 0;
        }

        pub fn isNan(self: Dyadic(mantissa_bits, exponent_bits)) bool {
            return self.exponent == int.maxVal(Exponent) and self.mantissa != 0;
        }

        pub fn isZero(self: Dyadic(mantissa_bits, exponent_bits)) bool {
            return self.exponent == int.minVal(Exponent) and self.mantissa == 0;
        }

        /// Normalizes the dyadic in place, ensuring the mantissa is
        /// left-aligned (i.e., the most significant bit is set).
        ///
        /// Parameters
        /// ----------
        /// `self` (`*Dyadic(mantissa_bits, exponent_bits)`):
        /// A pointer to the dyadic to normalize.
        ///
        /// Returns
        /// -------
        /// `void`
        pub fn normalize(self: *Dyadic(mantissa_bits, exponent_bits)) void {
            if (self.mantissa == 0) {
                self.exponent = int.minVal(Exponent);
                return;
            }

            const lz = @clz(self.mantissa);
            self.mantissa <<= lz;
            self.exponent -= types.scast(Exponent, lz);

            return;
        }

        pub fn add(
            x: Dyadic(mantissa_bits, exponent_bits),
            y: Dyadic(mantissa_bits, exponent_bits),
        ) Dyadic(mantissa_bits, exponent_bits) {
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
                return y
            else if (y.isZero())
                return x;

            // Addition or subtraction
            const cmp_abs: types.Cmp = x.abs().cmp(y.abs());
            if (x.positive == y.positive) {
                var result: Dyadic(mantissa_bits, exponent_bits) = if (cmp_abs == .gt)
                    _addAbs(x, y)
                else
                    _addAbs(y, x);
                result.positive = x.positive;
                return result;
            }

            if (cmp_abs == .eq)
                return .zero;

            var result: Dyadic(mantissa_bits, exponent_bits) = if (cmp_abs == .gt)
                _subAbs(x, y)
            else
                _subAbs(y, x);
            result.positive = if (cmp_abs == .gt) x.positive else y.positive;
            return result;
        }

        pub fn sub(
            x: Dyadic(mantissa_bits, exponent_bits),
            y: Dyadic(mantissa_bits, exponent_bits),
        ) Dyadic(mantissa_bits, exponent_bits) {
            return x.add(y.neg());
        }

        fn _addAbs(
            x: Dyadic(mantissa_bits, exponent_bits),
            y: Dyadic(mantissa_bits, exponent_bits),
        ) Dyadic(mantissa_bits, exponent_bits) {
            const Wide = std.meta.Int(.unsigned, mantissa_bits * 2);

            // |x| >= |y|, so exponent difference is non-negative
            const exp_diff: Exponent = x.exponent - y.exponent;

            const x_wide: Wide = types.scast(Wide, x.mantissa) << @intCast(mantissa_bits);
            const y_wide: Wide = types.scast(Wide, y.mantissa) << @intCast(mantissa_bits);

            var y_shifted: Wide = undefined;
            var sticky: u1 = 0;
            if (exp_diff >= 2 * mantissa_bits + 2) {
                y_shifted = 0;

                if (y.mantissa != 0)
                    sticky = 1;
            } else {
                y_shifted = y_wide >> @intCast(exp_diff);

                if ((y_wide & ((@as(Wide, 1) << @intCast(exp_diff)) - 1)) != 0)
                    sticky = 1;
            }

            const sum_ov = @addWithOverflow(x_wide, y_shifted);
            var sum: Wide = sum_ov[0];
            const carry: u1 = sum_ov[1];
            var exponent: Exponent = x.exponent;
            if (carry != 0) {
                if ((sum & 1) != 0)
                    sticky = 1;

                sum >>= 1;
                sum |= (@as(Wide, 1) << (mantissa_bits * 2 - 1));
                exponent += 1;
            }

            var mantissa: Mantissa = @truncate(sum >> @intCast(mantissa_bits));
            const remainder: Wide = sum & ((@as(Wide, 1) << @intCast(mantissa_bits)) - 1);
            const halfway: Wide = @as(Wide, 1) << @intCast(mantissa_bits - 1);
            var round_up = false;
            if (remainder > halfway) {
                round_up = true;
            } else if (remainder == halfway) {
                // Tie case
                if (sticky == 1 or (mantissa & 1) == 1)
                    round_up = true;
            }

            if (round_up) {
                const round = @addWithOverflow(mantissa, 1);
                mantissa = round[0];

                if (round[1] != 0) {
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    exponent +|= 1;
                }
            }

            // Check for overflow
            if (exponent == int.maxVal(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.maxVal(Exponent),
                    .positive = true,
                };

            return .{
                .mantissa = mantissa,
                .exponent = exponent,
                .positive = true,
            };
        }

        fn _subAbs(
            x: Dyadic(mantissa_bits, exponent_bits),
            y: Dyadic(mantissa_bits, exponent_bits),
        ) Dyadic(mantissa_bits, exponent_bits) {
            const Wide = std.meta.Int(.unsigned, mantissa_bits * 2);

            // |x| >= |y|, so exponent difference is non-negative
            const exp_diff: Exponent = x.exponent - y.exponent;

            const x_wide: Wide = types.scast(Wide, x.mantissa) << @intCast(mantissa_bits);
            const y_wide: Wide = types.scast(Wide, y.mantissa) << @intCast(mantissa_bits);

            var y_shifted: Wide = undefined;
            var sticky: u1 = 0;
            if (exp_diff >= mantissa_bits * 2) {
                y_shifted = 0;

                if (y.mantissa != 0)
                    sticky = 1;
            } else {
                y_shifted = y_wide >> @intCast(exp_diff);

                if ((y_wide & ((@as(Wide, 1) << @intCast(exp_diff)) - 1)) != 0)
                    sticky = 1;
            }

            var diff: Wide = x_wide - y_shifted;
            if (sticky == 1)
                diff -= 1;

            if (diff == 0)
                return .zero;

            const lz = @clz(diff);

            var exponent: Exponent = x.exponent -| types.scast(Exponent, lz);
            var mantissa: Mantissa = @truncate((diff << @intCast(lz)) >> @intCast(mantissa_bits));
            const remainder: Wide = (diff << @intCast(lz)) & ((@as(Wide, 1) << @intCast(mantissa_bits)) - 1);
            const halfway: Wide = @as(Wide, 1) << @intCast(mantissa_bits - 1);
            var round_up = false;
            if (remainder > halfway) {
                round_up = true;
            } else if (remainder == halfway) {
                // Tie case
                if (sticky == 1 or (mantissa & 1) == 1)
                    round_up = true;
            }

            if (round_up) {
                const round = @addWithOverflow(mantissa, 1);
                mantissa = round[0];

                if (round[1] != 0) {
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    exponent +|= 1;
                }
            }

            // Check for underflow
            if (exponent == int.minVal(Exponent))
                return .zero;

            return .{
                .mantissa = mantissa,
                .exponent = exponent,
                .positive = true,
            };
        }

        pub fn mul(
            x: Dyadic(mantissa_bits, exponent_bits),
            y: Dyadic(mantissa_bits, exponent_bits),
        ) Dyadic(mantissa_bits, exponent_bits) {
            // NaN check
            if (x.isNan() or y.isNan())
                return .nan;

            // Infinity check
            if (x.isInf()) {
                if (y.isZero())
                    return .nan // Inf * 0 = NaN
                else
                    return .{
                        .mantissa = 0,
                        .exponent = int.maxVal(Exponent),
                        .positive = x.positive == y.positive,
                    };
            } else if (y.isInf()) {
                if (x.isZero())
                    return .nan // 0 * Inf = NaN
                else
                    return .{
                        .mantissa = 0,
                        .exponent = int.maxVal(Exponent),
                        .positive = x.positive == y.positive,
                    };
            }

            // Zero check
            if (x.isZero() or y.isZero())
                return .{
                    .mantissa = 0,
                    .exponent = int.minVal(Exponent),
                    .positive = x.positive == y.positive,
                };

            // Multiplication
            const Wide = std.meta.Int(.unsigned, mantissa_bits * 2);

            // MSB will be at bit position mantissa_bits * 2 - 1 or mantissa_bits * 2 - 2
            const product: Wide = types.scast(Wide, x.mantissa) * types.scast(Wide, y.mantissa);
            var exponent: Exponent = x.exponent +| y.exponent -| types.scast(Exponent, mantissa_bits - 1);

            var mantissa: Mantissa = undefined;
            var remainder: Wide = undefined;
            var halfway: Wide = undefined;
            if ((product & (@as(Wide, 1) << @intCast(mantissa_bits * 2 - 1))) != 0) {
                // MSB at mantissa_bits * 2 - 1
                mantissa = @truncate(product >> @intCast(mantissa_bits));
                remainder = product & ((@as(Wide, 1) << @intCast(mantissa_bits)) - 1);
                halfway = @as(Wide, 1) << @intCast(mantissa_bits - 1);
                exponent +|= 1;
            } else {
                // MSB at mantissa_bits * 2 - 2
                mantissa = @truncate(product >> @intCast(mantissa_bits - 1));
                remainder = product & ((@as(Wide, 1) << @intCast(mantissa_bits - 1)) - 1);
                halfway = @as(Wide, 1) << @intCast(mantissa_bits - 2);
            }

            if (remainder > halfway or (remainder == halfway and (mantissa & 1) == 1)) {
                const round = @addWithOverflow(mantissa, 1);
                mantissa = round[0];

                if (round[1] != 0) {
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    exponent +|= 1;
                }
            }

            // Check for overflow
            if (exponent == int.maxVal(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.maxVal(Exponent),
                    .positive = x.positive == y.positive,
                };

            // Check for underflow
            if (exponent == int.minVal(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.minVal(Exponent),
                    .positive = x.positive == y.positive,
                };

            return .{
                .mantissa = mantissa,
                .exponent = exponent,
                .positive = x.positive == y.positive,
            };
        }

        pub fn div(
            x: Dyadic(mantissa_bits, exponent_bits),
            y: Dyadic(mantissa_bits, exponent_bits),
        ) Dyadic(mantissa_bits, exponent_bits) {
            // NaN check
            if (x.isNan() or y.isNan())
                return .nan;

            // Infinity check
            if (x.isInf()) {
                if (y.isInf())
                    return .nan // Inf / Inf = NaN
                else
                    return .{
                        .mantissa = 0,
                        .exponent = int.maxVal(Exponent),
                        .positive = x.positive == y.positive,
                    };
            } else if (y.isInf()) {
                // Finite / Inf = 0
                return .{
                    .mantissa = 0,
                    .exponent = int.minVal(Exponent),
                    .positive = x.positive == y.positive,
                };
            }

            // Zero check
            if (y.isZero()) {
                if (x.isZero())
                    return .nan // 0 / 0 = NaN
                else
                    return .{
                        .mantissa = 0,
                        .exponent = int.maxVal(Exponent),
                        .positive = x.positive == y.positive,
                    };
            } else if (x.isZero()) {
                return .{
                    .mantissa = 0,
                    .exponent = int.minVal(Exponent),
                    .positive = x.positive == y.positive,
                };
            }

            // Division
            const Wide = std.meta.Int(.unsigned, mantissa_bits * 2);

            const x_wide: Wide = types.scast(Wide, x.mantissa) << @intCast(mantissa_bits);
            const y_wide: Wide = types.scast(Wide, y.mantissa);
            const quotient: Wide = x_wide / y_wide;
            const remainder: Wide = x_wide % y_wide;

            var exponent: Exponent = x.exponent -| y.exponent +| types.scast(Exponent, mantissa_bits - 1);
            var mantissa: Mantissa = undefined;
            var round_remainder: Wide = undefined;
            var halfway: Wide = undefined;
            if ((quotient & (@as(Wide, 1) << @intCast(mantissa_bits))) != 0) {
                // MSB at mantissa_bits
                mantissa = @truncate(quotient >> 1);
                round_remainder = (quotient & 1) * y_wide + remainder;
                halfway = y_wide;
                exponent +|= 1;
            } else {
                // MSB at mantissa_bits - 1
                mantissa = @truncate(quotient);
                round_remainder = remainder;
                halfway = y_wide >> 1;
            }

            if (round_remainder > halfway or (round_remainder == halfway and (mantissa & 1) == 1)) {
                const round = @addWithOverflow(mantissa, 1);
                mantissa = round[0];

                if (round[1] != 0) {
                    mantissa = @as(Mantissa, 1) << (mantissa_bits - 1);
                    exponent +|= 1;
                }
            }

            // Check for overflow
            if (exponent == int.maxVal(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.maxVal(Exponent),
                    .positive = x.positive == y.positive,
                };

            // Check for underflow
            if (exponent == int.minVal(Exponent))
                return .{
                    .mantissa = 0,
                    .exponent = int.minVal(Exponent),
                    .positive = x.positive == y.positive,
                };

            return .{
                .mantissa = mantissa,
                .exponent = exponent,
                .positive = x.positive == y.positive,
            };
        }

        pub fn cmp(
            x: Dyadic(mantissa_bits, exponent_bits),
            y: Dyadic(mantissa_bits, exponent_bits),
        ) Cmp {
            // NaN check
            if (x.isNan()) {
                if (y.isNan())
                    return .eq
                else
                    return .gt;
            } else if (y.isNan()) {
                return .lt;
            }

            // Zero check (ignoring sign)
            if (x.isZero() and y.isZero())
                return .eq;

            // Sign check
            if (x.positive != y.positive) {
                if (x.positive)
                    return .gt
                else
                    return .lt;
            }

            // Infinity check
            if (x.isInf()) {
                if (y.isInf())
                    return .eq
                else
                    return if (x.positive) .gt else .lt;
            } else if (y.isInf()) {
                return if (x.positive) .lt else .gt;
            }

            // Exponent comparison
            if (x.exponent > y.exponent)
                return if (x.positive) .gt else .lt
            else if (x.exponent < y.exponent)
                return if (x.positive) .lt else .gt;

            // Mantissa comparison
            if (x.mantissa > y.mantissa)
                return if (x.positive) .gt else .lt
            else if (x.mantissa < y.mantissa)
                return if (x.positive) .lt else .gt;

            return .eq;
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
    };
}

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.dyadic) or !types.numericType(Y).le(.dyadic) or
        (types.numericType(X) != .dyadic and types.numericType(Y) != .dyadic))
        @compileError("zml.dyadic.add: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub fn add(x: anytype, y: anytype) Add(@TypeOf(x), @TypeOf(y)) {
    const R: type = Add(@TypeOf(x), @TypeOf(y));

    return types.scast(R, x).add(types.scast(R, y));
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.dyadic) or !types.numericType(Y).le(.dyadic) or
        (types.numericType(X) != .dyadic and types.numericType(Y) != .dyadic))
        @compileError("zml.dyadic.sub: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub fn sub(x: anytype, y: anytype) Sub(@TypeOf(x), @TypeOf(y)) {
    const R: type = Sub(@TypeOf(x), @TypeOf(y));

    return types.scast(R, x).sub(types.scast(R, y));
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.dyadic) or !types.numericType(Y).le(.dyadic) or
        (types.numericType(X) != .dyadic and types.numericType(Y) != .dyadic))
        @compileError("zml.dyadic.mul: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub fn mul(x: anytype, y: anytype) Mul(@TypeOf(x), @TypeOf(y)) {
    const R: type = Mul(@TypeOf(x), @TypeOf(y));

    return types.scast(R, x).mul(types.scast(R, y));
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.dyadic) or !types.numericType(Y).le(.dyadic) or
        (types.numericType(X) != .dyadic and types.numericType(Y) != .dyadic))
        @compileError("zml.dyadic.div: at least one of x or y must be a dyadic, the other must be a bool, an int, a float or a dyadic, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub fn div(x: anytype, y: anytype) Div(@TypeOf(x), @TypeOf(y)) {
    const R: type = Div(@TypeOf(x), @TypeOf(y));

    return types.scast(R, x).div(types.scast(R, y));
}

pub const sign = @import("dyadic/sign.zig").sign;
