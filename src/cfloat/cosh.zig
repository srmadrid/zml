const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

pub fn cosh(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.cosh: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    return .{
        .re = float.cosh(z.re) * float.cos(z.im),
        .im = float.sinh(z.re) * float.sin(z.im),
    };
}
