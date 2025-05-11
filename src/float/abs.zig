const types = @import("../types.zig");

pub inline fn abs(x: anytype) @TypeOf(x) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return @bitCast(@abs(x));
        },
        .float => {
            return @abs(x);
        },
        else => unreachable,
    }
}
