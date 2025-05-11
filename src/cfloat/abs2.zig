const types = @import("../types.zig");
const Scalar = types.Scalar;

pub fn abs2(z: anytype) Scalar(@TypeOf(z)) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return z.re * z.re + z.im * z.im;
}
