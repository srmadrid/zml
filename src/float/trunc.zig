const types = @import("../types.zig");

pub inline fn trunc(x: anytype) @TypeOf(x) {
    comptime if (types.numericType(@TypeOf(x)) != .float)
        @compileError("float.trunc: x must be a float, got " ++ @typeName(@TypeOf(x)));

    return @trunc(x);
}
