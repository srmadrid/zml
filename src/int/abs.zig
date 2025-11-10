const types = @import("../types.zig");

pub inline fn abs(x: anytype) @TypeOf(x) {
    comptime if (types.numericType(@TypeOf(x)) != .int)
        @compileError("x must be an int");

    if (comptime @TypeOf(x) == comptime_int) {
        return @abs(x);
    }

    return @bitCast(@abs(x));
}
