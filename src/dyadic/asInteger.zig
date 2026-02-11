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

    if (!types.isNumeric(X) or types.numericType(X) != .float)
        @compileError("zml.float.asInteger: x must be a float, got \n\tx: " ++ @typeName(X) ++ "\n");

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
    @compileError("asInteger is not implemented yet.");
}
