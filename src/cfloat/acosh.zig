const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const asinh = @import("asinh.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn acosh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.acosh: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    const rcls: u32 = classify.classify(zz.re);
    const icls: u32 = classify.classify(zz.im);

    if (rcls <= classify.INFINITE or icls <= classify.INFINITE) {
        var res: @TypeOf(zz) = undefined;
        if (icls == classify.INFINITE) {
            res.re = std.math.inf(Scalar(@TypeOf(zz)));

            if (rcls == classify.NAN) {
                res.im = std.math.nan(Scalar(@TypeOf(zz)));
            } else {
                res.im = float.copysign(if (rcls == classify.INFINITE) (if (zz.re < 0) float.pi(Scalar(@TypeOf(zz))) - float.pi_4(Scalar(@TypeOf(zz))) else float.pi_4(Scalar(@TypeOf(zz)))) else float.pi_2(Scalar(@TypeOf(zz))), zz.im);
            }
        } else if (rcls == classify.INFINITE) {
            res.re = std.math.inf(Scalar(@TypeOf(zz)));

            if (icls >= classify.ZERO) {
                res.im = float.copysign(if (std.math.signbit(zz.re)) float.pi(Scalar(@TypeOf(zz))) else 0, zz.im);
            } else {
                res.im = std.math.nan(Scalar(@TypeOf(zz)));
            }
        } else {
            res.re = std.math.nan(Scalar(@TypeOf(zz)));
            if (rcls == classify.ZERO) {
                res.im = float.pi_2(Scalar(@TypeOf(zz)));
            } else {
                res.im = std.math.nan(Scalar(@TypeOf(zz)));
            }
        }

        return res;
    } else if (rcls == classify.ZERO and icls == classify.ZERO) {
        return .{
            .re = 0,
            .im = float.copysign(float.pi_2(Scalar(@TypeOf(zz))), zz.im),
        };
    } else {
        var y: @TypeOf(zz) = .{
            .re = -zz.im,
            .im = zz.re,
        };

        y = asinh.kernel_asinh(y, 1);

        var res: @TypeOf(zz) = undefined;
        if (std.math.signbit(zz.im)) {
            res.re = y.re;
            res.im = -y.im;
        } else {
            res.re = -y.re;
            res.im = y.im;
        }

        return res;
    }
}
