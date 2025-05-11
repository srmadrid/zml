const types = @import("../types.zig");
const float = @import("../float.zig");
const Scalar = types.Scalar;

pub fn abs(z: anytype) Scalar(@TypeOf(z)) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return float.hypot(z.re, z.im);
}
