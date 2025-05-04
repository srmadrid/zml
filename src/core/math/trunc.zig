const types = @import("../types.zig");

pub inline fn trunc(x: anytype) @TypeOf(x) {
    comptime if (types.numericType(@TypeOf(x)) != .float)
        @compileError("x must be a float");

    switch (types.numericType(@TypeOf(x))) {
        .float => {
            return @trunc(x);
        },
        else => unreachable,
    }
}
