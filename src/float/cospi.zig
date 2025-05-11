const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn cospi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return cospi(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (float.abs(x) < std.math.floatEps(@TypeOf(x)))
                return 1;

            const xx: @TypeOf(x) = float.abs(x - 2 * float.round(0.5 * x));
            if (xx <= 0.25) {
                return float.cos(std.math.pi * xx);
            } else if (xx == 0.5) {
                return 0;
            } else if (xx <= 0.75) {
                return float.sin(std.math.pi * (0.5 - xx));
            } else {
                return -float.cos(std.math.pi * (1 - xx));
            }
        },
        else => unreachable,
    }
}
