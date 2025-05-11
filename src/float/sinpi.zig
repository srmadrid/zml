const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn sinpi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return sinpi(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (float.abs(x) < std.math.floatEps(@TypeOf(x))) {
                const ret: @TypeOf(x) = std.math.pi * x;

                if (float.abs(ret) < std.math.floatMin(@TypeOf(x))) {
                    const vret: @TypeOf(x) = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }

                return ret;
            }

            const y: @TypeOf(x) = x - 2 * float.round(0.5 * x);
            const absy: @TypeOf(x) = float.abs(y);
            if (absy == 0 or absy == 1) {
                return float.copysign(@as(@TypeOf(x), 0), x);
            } else if (absy <= 0.25) {
                return float.sin(std.math.pi * y);
            } else if (absy <= 0.75) {
                return float.copysign(float.cos(std.math.pi * (0.5 - absy)), y);
            } else {
                return float.copysign(float.sin(std.math.pi * (1 - absy)), y);
            }
        },
        else => unreachable,
    }
}
