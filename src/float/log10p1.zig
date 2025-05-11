const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub fn log10p1(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return log10p1(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (float.abs(x) < std.math.floatEps(@TypeOf(x)) / 4) {
                // Ensure appropriate underflows (a wider range than for log1p,
                // with potential for zero results from nonzero arguments, in
                // which case errno should be set based on the result with any
                // excess range and precision removed) even if the result of
                // multiplying by log10e is exact.
                const ret: @TypeOf(x) = std.math.log10e * x;

                if (float.abs(ret) < std.math.floatMin(@TypeOf(x))) {
                    const vret: @TypeOf(x) = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }

                return ret;
            }

            return std.math.log10e * float.log1p(x);
        },
        else => unreachable,
    }
}
