const std = @import("std");

const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");

pub fn asinh(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.asinh: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    return cfloat.asin(z.mulImag(1.0)).mulImag(-1.0);
}
