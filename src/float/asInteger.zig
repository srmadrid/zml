const std = @import("std");
const builtin = @import("builtin");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");

/// Only for internal use.
///
/// Converts a float type to a `Integer` as a view. After calling, you must execute:
/// ```zig
/// var i = asInteger(x);
/// i[0].limbs = &i[1];
/// ```
pub fn asInteger(x: anytype) !t: {
    const X: type = @TypeOf(x);
    var size = 0;
    if (X == comptime_float) {
        if (x == 0)
            size = 1
        else {
            size = std.math.log2(int.abs(types.scast(comptime_int, float.trunc(x)))) / 32 + 1;
        }
    } else {
        switch (@typeInfo(X).float.bits) {
            16 => size = 1,
            32 => size = 4,
            64 => size = 32,
            80 => size = 512,
            128 => size = 512,
            else => unreachable,
        }
    }

    break :t struct { integer.Integer, [size]u32 };
} {
    const X: type = @TypeOf(x);

    comptime if (types.numericType(X) != .float)
        @compileError("float.asInteger requires x to be a float type, got " ++ @typeName(X));

    const v: X = float.trunc(x);

    if (comptime X == comptime_float) {
        if (x == 0) {
            return .{
                .{
                    .limbs = undefined,
                    .size = 0,
                    ._llen = 1,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                },
                .{0},
            };
        }

        return @import("../int/asInteger.zig").asInteger(comptime types.scast(comptime_int, v));
    }

    if (!std.math.isFinite(v))
        return integer.Error.NotFinite;

    const bits: u16 = @typeInfo(X).float.bits;
    const size: u32 = switch (bits) {
        16 => 1,
        32 => 4,
        64 => 32,
        80 => 512,
        128 => 512,
        else => unreachable,
    };

    var limbs: [size]u32 = undefined;

    if (v == 0) {
        return .{
            .{
                .limbs = undefined,
                .size = 0,
                ._llen = size,
                .positive = true,
                .flags = .{ .owns_data = false, .writable = false },
            },
            limbs,
        };
    }

    var actual_size: u32 = 0;
    switch (bits) {
        16 => {
            const uvalue: u16 = @bitCast(v);
            const exponent: i16 = types.scast(i16, (uvalue & 0x7C00) >> 10) - 15;
            var mantissa: u32 = types.scast(u32, uvalue & 0x3FF);

            if (exponent == -15) {
                return; // Zero or subnormal
            } else {
                mantissa |= 0x400;
            }

            if (exponent > 10) {
                mantissa <<= @as(u5, @intCast(exponent - 10));
            } else {
                mantissa >>= @as(u5, @intCast(10 - exponent));
            }

            // Mantissa is now the integer value
            if (mantissa != 0) {
                limbs[0] = mantissa;
                actual_size = 1;
            }
        },
        32 => {
            const uvalue: u32 = @bitCast(v);
            const exponent: i32 = types.scast(i32, (uvalue & 0x7F800000) >> 23) - 127;
            var mantissa: u160 = @intCast(uvalue & 0x7FFFFF);

            if (exponent == -127) {
                return integer.Error.NotFinite; // Zero or subnormal
            } else {
                mantissa |= 0x800000;
            }

            if (exponent > 23) {
                mantissa <<= @as(u7, @intCast(exponent - 23));
            } else if (exponent >= 0) {
                mantissa >>= @as(u7, @intCast(23 - exponent));
            } else {
                mantissa = 0;
            }

            // Mantissa is now the integer value
            while (mantissa != 0) {
                limbs[actual_size] = @truncate(mantissa & 0xFFFFFFFF);
                actual_size += 1;
                mantissa >>= 32;
            }
        },
        64 => {
            const uvalue: u64 = @bitCast(v);
            const exponent: i32 = types.scast(i32, (uvalue & 0x7FF0000000000000) >> 52) - 1023;
            var mantissa: u1088 = @intCast(uvalue & 0xFFFFFFFFFFFFF);

            if (exponent == -1023) {
                return integer.Error.NotFinite; // Zero or subnormal
            } else {
                mantissa |= 0x10000000000000;
            }

            if (exponent > 52) {
                mantissa <<= @as(u11, @intCast(exponent - 52));
            } else if (exponent >= 0) {
                mantissa >>= @as(u11, @intCast(52 - exponent));
            } else {
                mantissa = 0;
            }

            // Mantissa is now the integer value
            while (mantissa != 0) {
                limbs[actual_size] = @truncate(mantissa & 0xFFFFFFFF);
                actual_size += 1;
                mantissa >>= 32;
            }
        },
        80 => {
            return asInteger(types.scast(f128, x));
        },
        128 => {
            const uvalue: u128 = @bitCast(v);
            const exponent: i32 = types.scast(i32, (uvalue & 0x7FFF0000000000000000000000000000) >> 112) - 16383;
            var mantissa: u16512 = @intCast(uvalue & 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF);

            if (exponent == -16383) {
                return integer.Error.NotFinite; // Zero or subnormal
            } else {
                mantissa |= 0x10000000000000000000000000000;
            }

            if (exponent > 112) {
                mantissa <<= @as(u15, @intCast(exponent - 112));
            } else if (exponent >= 0) {
                mantissa >>= @as(u15, @intCast(112 - exponent));
            } else {
                mantissa = 0;
            }

            // Mantissa is now the integer value
            while (mantissa != 0) {
                limbs[actual_size] = @truncate(mantissa & 0xFFFFFFFF);
                actual_size += 1;
                mantissa >>= 32;
            }
        },
        else => unreachable,
    }

    return .{
        .{
            .limbs = undefined,
            .size = actual_size,
            ._llen = size,
            .positive = x >= 0,
            .flags = .{ .owns_data = false, .writable = false },
        },
        limbs,
    };
}
