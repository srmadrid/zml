const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn tanh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.tanh: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    if (!std.math.isFinite(zz.re) or !std.math.isFinite(zz.im)) {
        @branchHint(.unlikely);
        var res: @TypeOf(zz) = undefined;
        if (std.math.isInf(zz.re)) {
            res.re = float.copysign(@as(Scalar(@TypeOf(zz)), 1), zz.re);
            if (std.math.isFinite(zz.im) and float.abs(zz.im) > 1) {
                const tmp = float.sincos(zz.im);
                res.im = float.copysign(@as(Scalar(@TypeOf(zz)), 0), tmp.sinx * tmp.cosx);
            } else {
                res.im = float.copysign(@as(Scalar(@TypeOf(zz)), 0), zz.im);
            }
        } else if (zz.im == 0) {
            res = zz;
        } else {
            if (zz.re == 0) {
                res.re = zz.re;
            } else {
                res.re = std.math.nan(Scalar(@TypeOf(zz)));
            }

            res.im = std.math.nan(Scalar(@TypeOf(zz)));
        }

        return res;
    } else {
        const t: i32 = scast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(zz))) - 1) * float.ln2(Scalar(@TypeOf(zz))) / 2);

        // tanh(zz+iy) = (sinh(2zz) + i*sin(2y))/(cosh(2zz) + cos(2y))
        // = (sinh(zz)*cosh(zz) + i*sin(y)*cos(y))/(sinh(zz)^2 + cos(y)^2).
        var sinizz: Scalar(@TypeOf(zz)) = undefined;
        var cosizz: Scalar(@TypeOf(zz)) = undefined;
        if (float.abs(zz.im) > std.math.floatMin(Scalar(@TypeOf(zz)))) {
            @branchHint(.likely);
            const tmp = float.sincos(zz.im);
            sinizz = tmp.sinx;
            cosizz = tmp.cosx;
        } else {
            sinizz = zz.im;
            cosizz = 1;
        }

        var res: @TypeOf(zz) = undefined;
        if (float.abs(zz.re) > scast(Scalar(@TypeOf(zz)), t)) {
            // Avoid intermediate overflow when the imaginary part of
            // the result may be subnormal.  Ignoring negligible terms,
            // the real part is +/- 1, the imaginary part is
            // sin(y)*cos(y)/sinh(zz)^2 = 4*sin(y)*cos(y)/exp(2zz).
            const exp_2t: Scalar(@TypeOf(zz)) = float.exp(2 * scast(Scalar(@TypeOf(zz)), t));

            res.re = float.copysign(@as(Scalar(@TypeOf(zz)), 1), zz.re);
            res.im = 4 * sinizz * cosizz;
            const zzzz: Scalar(@TypeOf(zz)) = float.abs(zz.re) - scast(Scalar(@TypeOf(zz)), t);
            res.im /= exp_2t;
            if (zzzz > scast(Scalar(@TypeOf(zz)), t)) {
                // Underflow (original real part of zz has absolute value > 2t).
                res.im /= exp_2t;
            } else {
                res.im /= float.exp(2 * zzzz);
            }
        } else {
            var sinhrzz: Scalar(@TypeOf(zz)) = undefined;
            var coshrzz: Scalar(@TypeOf(zz)) = undefined;
            if (float.abs(zz.re) > std.math.floatMin(Scalar(@TypeOf(zz)))) {
                sinhrzz = float.sinh(zz.re);
                coshrzz = float.cosh(zz.re);
            } else {
                sinhrzz = zz.re;
                coshrzz = 1;
            }

            var den: Scalar(@TypeOf(zz)) = undefined;
            if (float.abs(sinhrzz) > float.abs(cosizz) * std.math.floatEps(Scalar(@TypeOf(zz)))) {
                den = sinhrzz * sinhrzz + cosizz * cosizz;
            } else {
                den = cosizz * cosizz;
            }

            res.re = sinhrzz * coshrzz / den;
            res.im = sinizz * cosizz / den;
        }

        if (float.abs(res.re) < std.math.floatMin(Scalar(@TypeOf(zz)))) {
            const vresr: Scalar(@TypeOf(zz)) = res.re * res.re;
            std.mem.doNotOptimizeAway(vresr);
        }

        if (float.abs(res.im) < std.math.floatMin(Scalar(@TypeOf(zz)))) {
            const vresi: Scalar(@TypeOf(zz)) = res.im * res.im;
            std.mem.doNotOptimizeAway(vresi);
        }

        return res;
    }
}
