const types = @import("../types.zig");

pub inline fn abs(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!types.isNumeric(X) or types.numericType(X) != .int)
        @compileError("zml.int.maxVal: x must be an int, got \n\tx: " ++ @typeName(X) ++ "\n");

    if (comptime X == comptime_int)
        return @abs(x);

    return @bitCast(@abs(x));
}
