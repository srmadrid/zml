const std = @import("std");
const builtin = @import("builtin");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

/// Only for internal use.
///
/// Converts a float type to a `Complex` as a view. After calling, you must execute:
/// ```zig
/// var c = asComplex(x);
/// c[0].re.num.limbs = &c[1][0];
/// c[0].re.den.limbs = &c[1][1];
/// ```
pub fn asComplex(x: anytype) !t: {
    const X: type = @TypeOf(x);

    if (!types.isNumeric(X) or types.numericType(X) != .float)
        @compileError("zml.float.asComplex: x must be a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    var size = 0;
    if (X == comptime_float) {
        size = 512;
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

    break :t struct { complex.Complex(rational.Rational), [2][size]u32 };
} {
    const X: type = @TypeOf(x);

    comptime if (types.numericType(X) != .float)
        @compileError("float.asComplex requires x to be a float type, got " ++ @typeName(X));

    const size: u32 = if (X == comptime_float)
        512
    else switch (@typeInfo(X).float.bits) {
        16 => 1,
        32 => 4,
        64 => 32,
        80 => 512,
        128 => 512,
        else => unreachable,
    };

    const re = try @import("asRational.zig").asRational(x);

    return .{
        .{
            .re = .{
                .num = .{
                    .limbs = undefined,
                    .size = re[0].num.size,
                    ._llen = size,
                    .positive = re[0].num.positive,
                    .flags = .{ .owns_data = false, .writable = false },
                },
                .den = .{
                    .limbs = undefined,
                    .size = re[0].den.size,
                    ._llen = size,
                    .positive = re[0].den.positive,
                    .flags = .{ .owns_data = false, .writable = false },
                },
                .flags = .{ .owns_data = false, .writable = false },
            },
            .im = constants.zero(rational.Rational, .{}) catch unreachable,
            .flags = .{ .owns_data = false, .writable = false },
        },
        .{ re[1][0], re[1][1] },
    };
}
