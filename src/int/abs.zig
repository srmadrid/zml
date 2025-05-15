const types = @import("../types.zig");

pub inline fn abs(x: anytype) @TypeOf(x) {
    comptime if (types.numericType(@TypeOf(x)) != .int)
        @compileError("x must be an int");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return @bitCast(@abs(x));
        },
        else => unreachable,
    }
}
