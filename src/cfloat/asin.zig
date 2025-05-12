const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn asin(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return asin(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return asin(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            if (std.math.isNan(z.re) or std.math.isNan(z.im)) {
                if (z.re == 0) {
                    return z;
                } else if (std.math.isInf(z.re) or std.math.isInf(z.im)) {
                    return .{
                        .re = std.math.nan(Scalar(@TypeOf(z))),
                        .im = float.copysign(std.math.inf(Scalar(@TypeOf(z))), z.im),
                    };
                } else {
                    return .{
                        .re = std.math.nan(Scalar(@TypeOf(z))),
                        .im = std.math.nan(Scalar(@TypeOf(z))),
                    };
                }
            } else {
                var y: @TypeOf(z) = .{
                    .re = -z.im,
                    .im = z.re,
                };

                y = cfloat.asinh(y);

                return .{
                    .re = y.im,
                    .im = -y.re,
                };
            }
        },
        else => unreachable,
    }
}
