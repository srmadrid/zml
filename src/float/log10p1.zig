const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn log10p1(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.log10p1: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    if (float.abs(xx) < std.math.floatEps(@TypeOf(xx)) / 4) {
        // Ensure appropriate underflows (a wider range than for log1p,
        // with potential for zero results from nonzero arguments, in
        // which case errno should be set based on the result with any
        // excess range and precision removed) even if the result of
        // multiplying by log10e is exact.
        const ret: @TypeOf(xx) = std.math.log10e * xx;

        if (float.abs(ret) < std.math.floatMin(@TypeOf(xx))) {
            const vret: @TypeOf(xx) = ret * ret;
            std.mem.doNotOptimizeAway(vret);
        }

        return ret;
    }

    return std.math.log10e * float.log1p(xx);
}
