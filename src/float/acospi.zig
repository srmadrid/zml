const types = @import("../types.zig");
const float = @import("../float.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn acospi(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return acospi(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            if (float.abs(x) > 1) {
                @branchHint(.unlikely);
                return (x - x) / (x - x);
            }

            const ret: @TypeOf(x) = float.acos(x) / float.pi(@TypeOf(x));
            // Ensure that rounding upward for both acos and the division cannot
            // yield a return value from acospi greater than 1.
            return if (ret > 1) 1 else ret;
        },
        else => unreachable,
    }
}
