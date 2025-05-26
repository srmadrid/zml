const types = @import("../types.zig");

pub inline fn round(x: anytype) @TypeOf(x) {
    comptime if (types.numericType(@TypeOf(x)) != .float)
        @compileError("float.round: x must be a float, got " ++ @typeName(@TypeOf(x)));

    switch (types.numericType(@TypeOf(x))) {
        .float => {
            return @round(x);
        },
        else => unreachable,
    }
}
