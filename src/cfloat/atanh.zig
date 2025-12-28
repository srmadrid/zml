const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");

pub fn atanh(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.atanh: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    return cfloat.atan(z.mulImag(1.0)).mulImag(-1.0);
}
