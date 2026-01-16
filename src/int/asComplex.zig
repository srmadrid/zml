const std = @import("std");
const builtin = @import("builtin");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

/// Only for internal use.
///
/// Converts an int type to a `Complex` as a view. After calling, you must execute:
/// ```zig
/// var c = asComplex(x);
/// c[0].re.num.limbs = &c[1];
/// ```
pub fn asComplex(x: anytype) t: {
    const X: type = @TypeOf(x);

    if (!types.isNumeric(X) or types.numericType(X) != .int)
        @compileError("zml.int.asComplex: x must be an int, got \n\tx: " ++ @typeName(X) ++ "\n");

    var size = 0;
    if (X == comptime_int) {
        if (x == 0)
            size = 1
        else
            size = std.math.log2(int.abs(x)) / 32 + 1;
    } else {
        size = int.max(1, @typeInfo(X).int.bits / 32);
    }

    break :t struct { complex.Complex(rational.Rational), [size]u32 };
} {
    const renum = @import("asInteger.zig").asInteger(x);

    return .{
        .{
            .re = .{
                .num = renum[0],
                .den = constants.one(integer.Integer, .{}) catch unreachable,
                .flags = .{ .owns_data = false, .writable = false },
            },
            .im = constants.zero(rational.Rational, .{}) catch unreachable,
            .flags = .{ .owns_data = false, .writable = false },
        },
        renum[1],
    };
}
