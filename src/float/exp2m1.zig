const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub fn exp2m1(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return exp2m1(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (x >= -1 and x <= 1) {
                const ret: @TypeOf(x) = float.expm1(std.math.ln2 * x);

                if (float.abs(ret) < std.math.floatMin(@TypeOf(x))) {
                    const vret: @TypeOf(x) = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }

                return ret;
            } else if (x > std.math.floatMantissaBits(@TypeOf(x)) + 2) {
                const ret: @TypeOf(x) = float.exp2(x);
                return ret;
            } else if (x < -std.math.floatMantissaBits(@TypeOf(x)) - 2) {
                return -1;
            } else {
                return float.exp2(x) - 1;
            }
        },
        else => unreachable,
    }
}
