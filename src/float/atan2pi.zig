const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn atan2pi(y: anytype, x: anytype) EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))) {
    comptime if (types.numericType(@TypeOf(y)) != .int and types.numericType(@TypeOf(y)) != .float)
        @compileError("float.atan2pi: y must be an int or float, got " ++ @typeName(@TypeOf(y)));

    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.atan2pi: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const yy: EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))) = scast(EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))), y);
    const xx: EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))) = scast(EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))), x);

    const ret: @TypeOf(yy) = float.atan2(yy, xx) / std.math.pi;
    if (!std.math.isNan(ret)) {
        @branchHint(.likely);

        if (float.abs(ret) < std.math.floatMin(@TypeOf(yy))) {
            const vret: @TypeOf(yy) = ret * ret;
            std.mem.doNotOptimizeAway(vret);
        }
    }

    // Ensure that rounding away from zero for both atan2 and the
    // division cannot yield a return value from atan2pi with absolute
    // value greater than 1.
    return if (float.abs(ret) > 1) float.copysign(@as(@TypeOf(yy), 1), ret) else ret;
}
