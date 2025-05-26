const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn cospi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.cospi: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    if (float.abs(xx) < std.math.floatEps(@TypeOf(xx)))
        return 1;

    const y: @TypeOf(xx) = float.abs(xx - 2 * float.round(0.5 * xx));
    if (y <= 0.25) {
        return float.cos(std.math.pi * y);
    } else if (y == 0.5) {
        return 0;
    } else if (y <= 0.75) {
        return float.sin(std.math.pi * (0.5 - y));
    } else {
        return -float.cos(std.math.pi * (1 - y));
    }
}
