const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn tan(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.tan: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    if (!std.math.isFinite(zz.re) or !std.math.isFinite(zz.im)) {
        @branchHint(.unlikely);
        var res: @TypeOf(zz) = undefined;
        if (std.math.isInf(zz.im)) {
            if (std.math.isFinite(zz.re) and float.abs(zz.re) > 1) {
                const tmp = float.sincos(zz.re);
                const sinrzz: Scalar(@TypeOf(zz)) = tmp.sinx;
                const cosrzz: Scalar(@TypeOf(zz)) = tmp.cosx;

                res.re = float.copysign(@as(Scalar(@TypeOf(zz)), 0), sinrzz * cosrzz);
            } else {
                res.re = float.copysign(@as(Scalar(@TypeOf(zz)), 0), zz.re);
            }

            res.im = float.copysign(@as(Scalar(@TypeOf(zz)), 1), zz.im);
        } else if (zz.re == 0) {
            res = zz;
        } else {
            res.re = std.math.nan(Scalar(@TypeOf(zz)));
            if (zz.im == 0) {
                res.im = zz.im;
            } else {
                res.im = std.math.nan(Scalar(@TypeOf(zz)));
            }
        }

        return res;
    } else {
        const t: i32 = scast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(zz))) - 1) * float.ln2(Scalar(@TypeOf(zz))) / 2);

        // tan(zz+iy) = (sin(2zz) + i*sinh(2y))/(cos(2zz) + cosh(2y))
        // = (sin(zz)*cos(zz) + i*sinh(y)*cosh(y)/(cos(zz)^2 + sinh(y)^2).
        var sinrzz: Scalar(@TypeOf(zz)) = undefined;
        var cosrzz: Scalar(@TypeOf(zz)) = undefined;
        if (float.abs(zz.re) > std.math.floatMin(Scalar(@TypeOf(zz)))) {
            @branchHint(.likely);
            const tmp = float.sincos(zz.re);
            sinrzz = tmp.sinx;
            cosrzz = tmp.cosx;
        } else {
            sinrzz = zz.re;
            cosrzz = 1;
        }

        var res: @TypeOf(zz) = undefined;
        if (float.abs(zz.im) > scast(Scalar(@TypeOf(zz)), t)) {
            // Avoid intermediate overflow when the real part of the
            // result may be subnormal.  Ignoring negligible terms, the
            // imaginary part is +/- 1, the real part is
            // sin(zz)*cos(zz)/sinh(y)^2 = 4*sin(zz)*cos(zz)/exp(2y).
            const exp_2t: Scalar(@TypeOf(zz)) = float.exp(2 * scast(Scalar(@TypeOf(zz)), t));

            res.im = float.copysign(@as(Scalar(@TypeOf(zz)), 1), zz.im);
            res.re = 4 * sinrzz * cosrzz;
            const zzzz: Scalar(@TypeOf(zz)) = float.abs(zz.im) - scast(Scalar(@TypeOf(zz)), t);
            res.re /= exp_2t;
            if (zzzz > scast(Scalar(@TypeOf(zz)), t)) {
                // Underflow (original imaginary part of zz has absolute
                // value > 2t).
                res.re /= exp_2t;
            } else {
                res.re /= float.exp(2 * zzzz);
            }
        } else {
            var sinhizz: Scalar(@TypeOf(zz)) = undefined;
            var coshizz: Scalar(@TypeOf(zz)) = undefined;
            if (float.abs(zz.im) > std.math.floatMin(Scalar(@TypeOf(zz)))) {
                sinhizz = float.sinh(zz.im);
                coshizz = float.cosh(zz.im);
            } else {
                sinhizz = zz.im;
                coshizz = 1;
            }

            var den: Scalar(@TypeOf(zz)) = undefined;
            if (float.abs(sinhizz) > float.abs(cosrzz) * std.math.floatEps(Scalar(@TypeOf(zz)))) {
                den = cosrzz * cosrzz + sinhizz * sinhizz;
            } else {
                den = cosrzz * cosrzz;
            }

            res.re = sinrzz * cosrzz / den;
            res.im = sinhizz * coshizz / den;
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
