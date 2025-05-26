const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn exp2m1(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.exp2m1: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    if (xx >= -1 and xx <= 1) {
        const ret: @TypeOf(xx) = float.expm1(std.math.ln2 * xx);

        if (float.abs(ret) < std.math.floatMin(@TypeOf(xx))) {
            const vret: @TypeOf(xx) = ret * ret;
            std.mem.doNotOptimizeAway(vret);
        }

        return ret;
    } else if (xx > std.math.floatMantissaBits(@TypeOf(xx)) + 2) {
        return float.exp2(xx);
    } else if (xx < -std.math.floatMantissaBits(@TypeOf(xx)) - 2) {
        return -1;
    } else {
        return float.exp2(xx) - 1;
    }
}
