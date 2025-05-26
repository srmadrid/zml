const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn asin(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.asin: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    if (std.math.isNan(zz.re) or std.math.isNan(zz.im)) {
        if (zz.re == 0) {
            return zz;
        } else if (std.math.isInf(zz.re) or std.math.isInf(zz.im)) {
            return .{
                .re = std.math.nan(Scalar(@TypeOf(zz))),
                .im = float.copysign(std.math.inf(Scalar(@TypeOf(zz))), zz.im),
            };
        } else {
            return .{
                .re = std.math.nan(Scalar(@TypeOf(zz))),
                .im = std.math.nan(Scalar(@TypeOf(zz))),
            };
        }
    } else {
        var y: @TypeOf(zz) = .{
            .re = -zz.im,
            .im = zz.re,
        };

        y = cfloat.asinh(y);

        return .{
            .re = y.im,
            .im = -y.re,
        };
    }
}
