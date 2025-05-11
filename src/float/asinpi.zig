const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn asinpi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return asinpi(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (float.abs(x) > 1) {
                @branchHint(.unlikely);

                return (x - x) / (x - x);
            }

            const ret: @TypeOf(x) = float.asin(x) / float.pi(@TypeOf(x));

            if (float.abs(ret) < std.math.floatMin(@TypeOf(x))) {
                const vret: @TypeOf(x) = ret * ret;
                std.mem.doNotOptimizeAway(vret);
            }

            // Ensure that rounding away from zero for both asin and the
            // division cannot yield a return value from asinpi with absolute
            // value greater than 0.5.
            return if (float.abs(ret) > 0.5) float.copysign(@as(@TypeOf(x), 0.5), ret) else ret;
        },
        else => unreachable,
    }
}
