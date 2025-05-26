const types = @import("../types.zig");
const Scalar = types.Scalar;

pub fn abs2(z: anytype) Scalar(@TypeOf(z)) {
    comptime if (types.numericType(@TypeOf(z)) == .cfloat)
        @compileError("cfloat.abs2: z must be an int or float, got " ++ @typeName(@TypeOf(z)));

    return z.re * z.re + z.im * z.im;
}
