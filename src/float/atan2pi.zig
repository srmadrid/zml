const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn atan2pi(y: anytype, x: anytype) EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))) {
    comptime if (!types.isFixedPrecision(@TypeOf(y)) or types.isComplex(@TypeOf(y)))
        @compileError("y must be an int or float");

    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(y))) {
        .int => {
            switch (types.numericType(@TypeOf(x))) {
                .int => return atan2pi(cast(EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))), y, .{}), cast(EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))), x, .{})),
                .float => return atan2pi(cast(Coerce(@TypeOf(y), @TypeOf(x)), y, .{}), cast(Coerce(@TypeOf(y), @TypeOf(x)), x, .{})),
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(x))) {
                .int => return atan2pi(cast(Coerce(@TypeOf(y), @TypeOf(x)), y, .{}), cast(Coerce(@TypeOf(y), @TypeOf(x)), x, .{})),
                .float => {
                    const yy: Coerce(@TypeOf(y), @TypeOf(x)) = cast(Coerce(@TypeOf(y), @TypeOf(x)), y, .{});
                    const xx: Coerce(@TypeOf(y), @TypeOf(x)) = cast(Coerce(@TypeOf(y), @TypeOf(x)), x, .{});

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
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}
