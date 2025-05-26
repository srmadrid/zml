const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const x2y2m1 = @import("../float/x2y2m1.zig").x2y2m1;
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn atanh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.atanh: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    const rcls: u32 = classify.classify(zz.re);
    const icls: u32 = classify.classify(zz.im);

    if (rcls <= classify.INFINITE or icls <= classify.INFINITE) {
        @branchHint(.unlikely);
        var res: @TypeOf(zz) = undefined;
        if (icls == classify.INFINITE) {
            res.re = float.copysign(@as(Scalar(@TypeOf(zz)), 0), zz.re);
            res.im = float.copysign(float.pi_2(Scalar(@TypeOf(zz))), zz.im);
        } else if (rcls == classify.INFINITE or rcls == classify.ZERO) {
            res.re = float.copysign(@as(Scalar(@TypeOf(zz)), 0), zz.re);
            if (icls >= classify.ZERO) {
                res.im = float.copysign(float.pi_2(Scalar(@TypeOf(zz))), zz.im);
            } else {
                res.im = std.math.nan(Scalar(@TypeOf(zz)));
            }
        } else {
            res.re = std.math.nan(Scalar(@TypeOf(zz)));
            res.im = std.math.nan(Scalar(@TypeOf(zz)));
        }

        return res;
    } else if (rcls == classify.ZERO and icls == classify.ZERO) {
        @branchHint(.unlikely);
        return zz;
    } else {
        var res: @TypeOf(zz) = undefined;
        if (float.abs(zz.re) >= 16 / std.math.floatEps(Scalar(@TypeOf(zz))) or float.abs(zz.im) >= 16 / std.math.floatEps(Scalar(@TypeOf(zz)))) {
            res.im = float.copysign(float.pi_2(Scalar(@TypeOf(zz))), zz.im);
            if (float.abs(zz.im) <= 1) {
                res.re = 1 / zz.re;
            } else if (float.abs(zz.re) <= 1) {
                res.re = zz.re / zz.im / zz.im;
            } else {
                const h: Scalar(@TypeOf(zz)) = float.hypot(zz.re / 2, zz.im / 2);
                res.re = zz.re / h / h / 4;
            }
        } else {
            if (float.abs(zz.re) == 1 and float.abs(zz.im) < std.math.floatEps(Scalar(@TypeOf(zz))) * std.math.floatEps(Scalar(@TypeOf(zz)))) {
                res.re = (float.copysign(@as(Scalar(@TypeOf(zz)), 0.5), zz.re) * (float.ln2(Scalar(@TypeOf(zz))) - float.log(float.abs(zz.im))));
            } else {
                var I2: Scalar(@TypeOf(zz)) = 0;
                if (float.abs(zz.im) >= std.math.floatEps(Scalar(@TypeOf(zz))) * std.math.floatEps(Scalar(@TypeOf(zz))))
                    I2 = zz.im * zz.im;

                var num: Scalar(@TypeOf(zz)) = 1 + zz.re;
                num = I2 + num * num;

                var den: Scalar(@TypeOf(zz)) = 1 - zz.re;
                den = I2 + den * den;

                const f: Scalar(@TypeOf(zz)) = num / den;
                if (f < 0.5) {
                    res.re = 0.25 * float.log(f);
                } else {
                    num = 4 * zz.re;
                    res.re = 0.25 * float.log1p(num / den);
                }
            }
            var absr: Scalar(@TypeOf(zz)) = float.abs(zz.re);
            var absi: Scalar(@TypeOf(zz)) = float.abs(zz.im);
            if (absr < absi) {
                const t: Scalar(@TypeOf(zz)) = absr;
                absr = absi;
                absi = t;
            }

            var den: Scalar(@TypeOf(zz)) = undefined;
            if (absi < std.math.floatEps(Scalar(@TypeOf(zz))) / 2) {
                den = (1 - absr) * (1 + absr);
                if (den == 0)
                    den = 0;
            } else if (absr >= 1) {
                den = (1 - absr) * (1 + absr) - absi * absi;
            } else if (absr >= 0.75 or absi >= 0.5) {
                den = -x2y2m1(absr, absi);
            } else {
                den = (1 - absr) * (1 + absr) - absi * absi;
            }

            res.im = 0.5 * float.atan2(2 * zz.im, den);
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
