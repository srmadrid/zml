const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const Scalar = types.Scalar;
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const Cfloat = cfloat.Cfloat;

pub fn log(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.log: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz = types.scast(Cfloat(types.EnsureFloat(types.Scalar(@TypeOf(z)))), z);

    var p: types.Scalar(@TypeOf(zz)) = cfloat.abs(zz);
    p = float.log(p);
    return .{
        .re = p,
        .im = float.atan2(zz.im, zz.re),
    };
}
