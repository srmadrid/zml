const std = @import("std");

const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");

pub fn acosh(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.acosh: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    return cfloat.acos(z).mulImag(1.0);
}
