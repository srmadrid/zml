const types = @import("../types.zig");

pub inline fn ceil(x: anytype) @TypeOf(x) {
    comptime if (types.numericType(@TypeOf(x)) != .float)
        @compileError("float.ceil: x must be a float, got " ++ @typeName(@TypeOf(x)));

    return @ceil(x);
}
