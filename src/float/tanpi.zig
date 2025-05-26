const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn tanpi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.tanpi: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    if (float.abs(xx) < std.math.floatEps(@TypeOf(xx))) {
        const ret: @TypeOf(xx) = std.math.pi * xx;

        if (float.abs(ret) < std.math.floatMin(@TypeOf(xx))) {
            const vret: @TypeOf(xx) = ret * ret;
            std.mem.doNotOptimizeAway(vret);
        }

        return ret;
    }

    var y: @TypeOf(xx) = xx - 2 * float.round(0.5 * xx);
    var absy: @TypeOf(xx) = float.abs(y);
    if (absy == 0) {
        // For even integers, return +0 if positive and -0 if negative (so
        // matching sinpi(xx)/cospi(xx)).
        return float.copysign(@as(@TypeOf(xx), 0), xx);
    } else if (absy == 1) {
        // For odd integers, return -0 if positive and +0 if negative (so
        // matching sinpi(xx)/cospi(xx)).
        return float.copysign(@as(@TypeOf(xx), 0), -xx);
    } else if (absy == 0.5) {
        // Return infinity with positive sign for an even integer + 0.5
        // and negative sign for an odd integer + 0.5 (so matching
        // sinpi(xx)/cospi(xx)).
        return 1 / float.copysign(@as(@TypeOf(xx), 0), y);
    } else if (absy > 0.5) {
        // Now we only care about the value of X mod 1, not mod 2.
        y -= float.copysign(@as(@TypeOf(xx), 1), y);
        absy = float.abs(y);
    }

    if (absy <= 0.25) {
        return float.tan(std.math.pi * y);
    } else {
        return float.copysign(1 / float.tan(std.math.pi * (0.5 - absy)), y);
    }
}
