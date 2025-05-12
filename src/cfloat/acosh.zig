const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const asinh = @import("asinh.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn acosh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return acosh(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return acosh(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            if (rcls <= classify.INFINITE or icls <= classify.INFINITE) {
                var res: @TypeOf(z) = undefined;
                if (icls == classify.INFINITE) {
                    res.re = std.math.inf(Scalar(@TypeOf(z)));

                    if (rcls == classify.NAN) {
                        res.im = std.math.nan(Scalar(@TypeOf(z)));
                    } else {
                        res.im = float.copysign(if (rcls == classify.INFINITE) (if (z.re < 0) float.pi(Scalar(@TypeOf(z))) - float.pi_4(Scalar(@TypeOf(z))) else float.pi_4(Scalar(@TypeOf(z)))) else float.pi_2(Scalar(@TypeOf(z))), z.im);
                    }
                } else if (rcls == classify.INFINITE) {
                    res.re = std.math.inf(Scalar(@TypeOf(z)));

                    if (icls >= classify.ZERO) {
                        res.im = float.copysign(if (std.math.signbit(z.re)) float.pi(Scalar(@TypeOf(z))) else 0, z.im);
                    } else {
                        res.im = std.math.nan(Scalar(@TypeOf(z)));
                    }
                } else {
                    res.re = std.math.nan(Scalar(@TypeOf(z)));
                    if (rcls == classify.ZERO) {
                        res.im = float.pi_2(Scalar(@TypeOf(z)));
                    } else {
                        res.im = std.math.nan(Scalar(@TypeOf(z)));
                    }
                }

                return res;
            } else if (rcls == classify.ZERO and icls == classify.ZERO) {
                return .{
                    .re = 0,
                    .im = float.copysign(float.pi_2(Scalar(@TypeOf(z))), z.im),
                };
            } else {
                var y: @TypeOf(z) = .{
                    .re = -z.im,
                    .im = z.re,
                };

                y = asinh.kernel_asinh(y, 1);

                var res: @TypeOf(z) = undefined;
                if (std.math.signbit(z.im)) {
                    res.re = y.re;
                    res.im = -y.im;
                } else {
                    res.re = -y.re;
                    res.im = y.im;
                }

                return res;
            }
        },
        else => unreachable,
    }
}
