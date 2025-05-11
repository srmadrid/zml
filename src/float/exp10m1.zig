const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub fn exp10m1(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return exp10m1(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (x >= -0.5 and x <= 0.5) {
                const ret: @TypeOf(x) = float.expm1(std.math.ln10 * x);

                if (float.abs(ret) < std.math.floatMin(@TypeOf(x))) {
                    const vret: @TypeOf(x) = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }

                return ret;
            } else if (x > std.math.floatMantissaBits(@TypeOf(x)) / 3 + 2) {
                const ret: @TypeOf(x) = float.exp10(x);
                return ret;
            } else if (x < -std.math.floatMantissaBits(@TypeOf(x)) / 3 - 2) {
                return -1;
            } else {
                return float.exp10(x) - 1;
            }
        },
        else => unreachable,
    }
}
