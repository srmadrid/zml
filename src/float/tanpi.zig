const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn tanpi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return tanpi(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (float.abs(x) < std.math.floatEps(@TypeOf(x))) {
                const ret: @TypeOf(x) = std.math.pi * x;

                if (float.abs(ret) < std.math.floatMin(@TypeOf(x))) {
                    const vret: @TypeOf(x) = x * x;
                    std.mem.doNotOptimizeAway(vret);
                }

                return ret;
            }

            var y: @TypeOf(x) = x - 2 * float.round(0.5 * x);
            var absy: @TypeOf(x) = float.abs(y);
            if (absy == 0) {
                // For even integers, return +0 if positive and -0 if negative (so
                // matching sinpi(x)/cospi(x)).
                return float.copysign(@as(@TypeOf(x), 0), x);
            } else if (absy == 1) {
                // For odd integers, return -0 if positive and +0 if negative (so
                // matching sinpi(x)/cospi(x)).
                return float.copysign(@as(@TypeOf(x), 0), -x);
            } else if (absy == 0.5) {
                // Return infinity with positive sign for an even integer + 0.5
                // and negative sign for an odd integer + 0.5 (so matching
                // sinpi(x)/cospi(x)).
                return 1 / float.copysign(@as(@TypeOf(x), 0), y);
            } else if (absy > 0.5) {
                // Now we only care about the value of X mod 1, not mod 2.
                y -= float.copysign(@as(@TypeOf(x), 1), y);
                absy = float.abs(y);
            }

            if (absy <= 0.25) {
                return float.tan(std.math.pi * y);
            } else {
                return float.copysign(1 / float.tan(std.math.pi * (0.5 - absy)), y);
            }
        },
        else => unreachable,
    }
}
