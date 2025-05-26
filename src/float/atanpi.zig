const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn atanpi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.atanpi: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    const ret: @TypeOf(xx) = float.atan(xx) / std.math.pi;
    if (!std.math.isNan(x)) {
        @branchHint(.likely);

        if (float.abs(ret) < std.math.floatMin(@TypeOf(xx))) {
            const vret: @TypeOf(xx) = ret * ret;
            std.mem.doNotOptimizeAway(vret);
        }
    }

    // Ensure that rounding away from zero for both atan and the
    // division cannot yield a return value from atanpi with absolute
    // value greater than 0.5.
    return if (float.abs(ret) > 0.5) float.copysign(@as(@TypeOf(xx), 0.5), ret) else ret;
}
