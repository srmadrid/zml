const std = @import("std");
const types = @import("../../types.zig");
const math = @import("../../math.zig");
const Scalar = types.Scalar;

pub fn logabs(z: anytype) Scalar(@TypeOf(z)) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    const rabs = @abs(z.re);
    const iabs = @abs(z.im);
    var max: @TypeOf(z.re, z.im) = undefined;
    var u: @TypeOf(z.re, z.im) = undefined;

    if (rabs >= iabs) {
        max = rabs;
        u = iabs / rabs;
    } else {
        max = iabs;
        u = rabs / iabs;
    }

    return math.log(max) + math.log1p(u * u) / 2;
}
