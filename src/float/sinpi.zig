const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn sinpi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.sinpi: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    if (float.abs(xx) < std.math.floatEps(@TypeOf(xx))) {
        const ret: @TypeOf(xx) = std.math.pi * xx;

        if (float.abs(ret) < std.math.floatMin(@TypeOf(xx))) {
            const vret: @TypeOf(xx) = ret * ret;
            std.mem.doNotOptimizeAway(vret);
        }

        return ret;
    }

    const y: @TypeOf(xx) = xx - 2 * float.round(0.5 * xx);
    const absy: @TypeOf(xx) = float.abs(y);
    if (absy == 0 or absy == 1) {
        return float.copysign(@as(@TypeOf(xx), 0), xx);
    } else if (absy <= 0.25) {
        return float.sin(std.math.pi * y);
    } else if (absy <= 0.75) {
        return float.copysign(float.cos(std.math.pi * (0.5 - absy)), y);
    } else {
        return float.copysign(float.sin(std.math.pi * (1 - absy)), y);
    }
}
