const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn asinpi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.asinpi: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    if (float.abs(xx) > 1) {
        @branchHint(.unlikely);

        return (xx - xx) / (xx - xx);
    }

    const ret: @TypeOf(xx) = float.asin(xx) / float.pi(@TypeOf(xx));

    if (float.abs(ret) < std.math.floatMin(@TypeOf(xx))) {
        const vret: @TypeOf(xx) = ret * ret;
        std.mem.doNotOptimizeAway(vret);
    }

    // Ensure that rounding away from zero for both asin and the
    // division cannot yield a return value from asinpi with absolute
    // value greater than 0.5.
    return if (float.abs(ret) > 0.5) float.copysign(@as(@TypeOf(xx), 0.5), ret) else ret;
}
