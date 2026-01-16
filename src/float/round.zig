const types = @import("../types.zig");

pub inline fn round(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!types.isNumeric(X) or types.numericType(X) != .float)
        @compileError("zml.float.round: x must be a float, got \n\tx: " ++ @typeName(@TypeOf(x)) ++ "\n");

    return @round(x);
}
