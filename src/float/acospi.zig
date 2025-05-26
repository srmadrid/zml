const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub fn acospi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.acospi: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    const xx: EnsureFloat(@TypeOf(x)) = scast(EnsureFloat(@TypeOf(x)), x);

    if (float.abs(xx) > 1) {
        @branchHint(.unlikely);
        return (xx - xx) / (xx - xx);
    }

    const ret: @TypeOf(xx) = float.acos(xx) / float.pi(@TypeOf(xx));
    // Ensure that rounding upward for both acos and the division cannot
    // yield a return value from acospi greater than 1.
    return if (ret > 1) 1 else ret;
}
