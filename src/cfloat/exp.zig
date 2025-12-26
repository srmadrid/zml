const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

pub fn exp(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.acos: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    const r: @TypeOf(z.re) = float.exp(z.re);
    return .{
        .re = r * float.cos(z.im),
        .im = r * float.sin(z.im),
    };
}
