const types = @import("../types.zig");
const Scalar = types.Scalar;
const float = @import("../float.zig");

pub fn abs(z: anytype) Scalar(@TypeOf(z)) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.abs: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    return float.hypot(z.re, z.im);
}
