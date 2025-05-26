const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn exp10m1(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.exp10m1: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    if (xx >= -0.5 and xx <= 0.5) {
        const ret: @TypeOf(xx) = float.expm1(std.math.ln10 * xx);

        if (float.abs(ret) < std.math.floatMin(@TypeOf(xx))) {
            const vret: @TypeOf(xx) = ret * ret;
            std.mem.doNotOptimizeAway(vret);
        }

        return ret;
    } else if (xx > std.math.floatMantissaBits(@TypeOf(xx)) / 3 + 2) {
        const ret: @TypeOf(xx) = float.exp10(xx);
        return ret;
    } else if (xx < -std.math.floatMantissaBits(@TypeOf(xx)) / 3 - 2) {
        return -1;
    } else {
        return float.exp10(xx) - 1;
    }
}
