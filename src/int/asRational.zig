const std = @import("std");
const builtin = @import("builtin");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

/// Only for internal use.
///
/// Converts an int type to a `Rational` as a view. After calling, you must execute:
/// ```zig
/// var r = asRational(x);
/// r[0].num.limbs = &r[1];
/// ```
pub fn asRational(x: anytype) t: {
    const X: type = @TypeOf(x);
    var size = 0;
    if (X == comptime_int) {
        if (x == 0)
            size = 1
        else
            size = std.math.log2(int.abs(x)) / 32 + 1;
    } else {
        switch (@typeInfo(X).int.bits) {
            8, 16, 32 => size = 1,
            64 => size = 2,
            128 => size = 4,
            else => unreachable,
        }
    }

    break :t struct { rational.Rational, [size]u32 };
} {
    const X: type = @TypeOf(x);

    comptime if (types.numericType(X) != .int)
        @compileError("int.asRational requires x to be an int type, got " ++ @typeName(X));

    const num = @import("asInteger.zig").asInteger(x);

    return .{
        .{
            .num = num[0],
            .den = constants.one(integer.Integer, .{}) catch unreachable,
            .flags = .{ .owns_data = false, .writable = false },
        },
        num[1],
    };
}
