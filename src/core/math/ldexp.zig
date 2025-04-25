const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn ldexp(x: anytype, n: i32) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return ldexp(cast(EnsureFloat(@TypeOf(x)), x, .{}), n);
        },
        .float => {
            if (!std.math.isFinite(x) or x == 0)
                return x + x;

            return math.scalbn(x, n);
        },
        else => unreachable,
    }
}
