const std = @import("std");
const builtin = @import("builtin");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

/// Only for internal use.
///
/// Converts a float type to a `Rational` as a view. After calling, you must execute:
/// ```zig
/// var r = asRational(x);
/// r[0].num.limbs = &r[1][0];
/// r[0].den.limbs = &r[1][1];
/// ```
pub fn asRational(x: anytype) !t: {
    const X: type = @TypeOf(x);

    if (!types.isNumeric(X) or types.numericType(X) != .float)
        @compileError("zml.float.asRational: x must be a float, got \n\tx: " ++ @typeName(X) ++ "\n");

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

    break :t struct { rational.Rational, [2][size]u32 };
} {
    @compileError("asRational is not implemented yet.");
}
