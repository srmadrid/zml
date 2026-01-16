const types = @import("../types.zig");

pub inline fn abs(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!types.isNumeric(X) or types.numericType(X) != .float)
        @compileError("zml.float.abs: x must be a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return @abs(x);
}
