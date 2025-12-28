const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

pub fn tanh(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.tanh: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    const d: @TypeOf(z.re) = float.cosh(2.0 * z.re) + float.cos(2.0 * z.im);
    return .{
        .re = float.cosh(2.0 * z.re) / d,
        .im = float.sinh(2.0 * z.im) / d,
    };
}
