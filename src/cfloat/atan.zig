const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const x2y2m1 = @import("../float/x2y2m1.zig").x2y2m1;
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn atan(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.atan: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    const rcls: u32 = classify.classify(zz.re);
    const icls: u32 = classify.classify(zz.im);

    if (rcls <= classify.INFINITE or icls <= classify.INFINITE) {
        @branchHint(.unlikely);
        var res: @TypeOf(zz) = undefined;
        if (rcls == classify.INFINITE) {
            res.re = float.copysign(float.pi_2(Scalar(@TypeOf(zz))), zz.re);
            res.im = float.copysign(@as(Scalar(@TypeOf(zz)), 0), zz.im);
        } else if (icls == classify.INFINITE) {
            if (rcls >= classify.ZERO) {
                res.re = float.copysign(float.pi_2(Scalar(@TypeOf(zz))), zz.re);
            } else {
                res.re = std.math.nan(Scalar(@TypeOf(zz)));
            }

            res.im = float.copysign(@as(Scalar(@TypeOf(zz)), 0), zz.im);
        } else if (icls == classify.ZERO or icls == classify.INFINITE) {
            res.re = std.math.nan(Scalar(@TypeOf(zz)));
            res.im = float.copysign(@as(Scalar(@TypeOf(zz)), 0), zz.im);
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
            res.re = float.copysign(float.pi_2(Scalar(@TypeOf(zz))), zz.re);
            if (float.abs(zz.re) <= 1) {
                res.im = 1 / zz.im;
            } else if (float.abs(zz.im) <= 1) {
                res.im = zz.im / zz.re / zz.re;
            } else {
                const h: Scalar(@TypeOf(zz)) = float.hypot(zz.re / 2, zz.im / 2);
                res.im = zz.im / h / h / 4;
            }
        } else {
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

            res.re = 0.5 * float.atan2(2 * zz.re, den);

            if (float.abs(zz.im) == 1 and float.abs(zz.re) < std.math.floatEps(Scalar(@TypeOf(zz))) * std.math.floatEps(Scalar(@TypeOf(zz)))) {
                res.im = (float.copysign(@as(Scalar(@TypeOf(zz)), 0.5), zz.im) * (float.ln2(Scalar(@TypeOf(zz))) - float.log(float.abs(zz.re))));
            } else {
                var r2: Scalar(@TypeOf(zz)) = 0;
                if (float.abs(zz.re) >= std.math.floatEps(Scalar(@TypeOf(zz))) * std.math.floatEps(Scalar(@TypeOf(zz))))
                    r2 = zz.re * zz.re;

                var num: Scalar(@TypeOf(zz)) = zz.im + 1;
                num = r2 + num * num;

                den = zz.im - 1;
                den = r2 + den * den;

                const f: Scalar(@TypeOf(zz)) = num / den;
                if (f < 0.5) {
                    res.im = 0.25 * float.log(f);
                } else {
                    num = 4 * zz.im;
                    res.im = 0.25 * float.log1p(num / den);
                }
            }
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
