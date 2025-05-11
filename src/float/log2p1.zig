const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub fn log2p1(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return log2p1(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (float.abs(x) < std.math.floatEps(@TypeOf(x)) / 4) {
                // Avoid spurious underflows when log1p underflows but log2p1
                // should not.
                const ret: @TypeOf(x) = std.math.log2e * x;

                if (float.abs(ret) < std.math.floatMin(@TypeOf(x))) {
                    const vret: @TypeOf(x) = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }

                return ret;
            }
            return std.math.log2e * float.log1p(x);
        },
        else => unreachable,
    }
}
