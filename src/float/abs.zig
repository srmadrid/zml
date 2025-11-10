const types = @import("../types.zig");

pub inline fn abs(x: anytype) @TypeOf(x) {
    comptime if (types.numericType(@TypeOf(x)) != .float)
        @compileError("x must be a float");

    return @abs(x);
}
