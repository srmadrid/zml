const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const x2y2m1 = @import("../float/x2y2m1.zig").x2y2m1;
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn atan(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return atan(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return atan(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            if (rcls <= classify.INFINITE or icls <= classify.INFINITE) {
                @branchHint(.unlikely);
                var res: @TypeOf(z) = undefined;
                if (rcls == classify.INFINITE) {
                    res.re = float.copysign(float.pi_2(Scalar(@TypeOf(z))), z.re);
                    res.im = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.im);
                } else if (icls == classify.INFINITE) {
                    if (rcls >= classify.ZERO) {
                        res.re = float.copysign(float.pi_2(Scalar(@TypeOf(z))), z.re);
                    } else {
                        res.re = std.math.nan(Scalar(@TypeOf(z)));
                    }

                    res.im = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.im);
                } else if (icls == classify.ZERO or icls == classify.INFINITE) {
                    res.re = std.math.nan(Scalar(@TypeOf(z)));
                    res.im = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.im);
                } else {
                    res.re = std.math.nan(Scalar(@TypeOf(z)));
                    res.im = std.math.nan(Scalar(@TypeOf(z)));
                }

                return res;
            } else if (rcls == classify.ZERO and icls == classify.ZERO) {
                @branchHint(.unlikely);
                return z;
            } else {
                var res: @TypeOf(z) = undefined;
                if (float.abs(z.re) >= 16 / std.math.floatEps(Scalar(@TypeOf(z))) or float.abs(z.im) >= 16 / std.math.floatEps(Scalar(@TypeOf(z)))) {
                    res.re = float.copysign(float.pi_2(Scalar(@TypeOf(z))), z.re);
                    if (float.abs(z.re) <= 1) {
                        res.im = 1 / z.im;
                    } else if (float.abs(z.im) <= 1) {
                        res.im = z.im / z.re / z.re;
                    } else {
                        const h: Scalar(@TypeOf(z)) = float.hypot(z.re / 2, z.im / 2);
                        res.im = z.im / h / h / 4;
                    }
                } else {
                    var absr: Scalar(@TypeOf(z)) = float.abs(z.re);
                    var absi: Scalar(@TypeOf(z)) = float.abs(z.im);
                    if (absr < absi) {
                        const t: Scalar(@TypeOf(z)) = absr;
                        absr = absi;
                        absi = t;
                    }

                    var den: Scalar(@TypeOf(z)) = undefined;
                    if (absi < std.math.floatEps(Scalar(@TypeOf(z))) / 2) {
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

                    res.re = 0.5 * float.atan2(2 * z.re, den);

                    if (float.abs(z.im) == 1 and float.abs(z.re) < std.math.floatEps(Scalar(@TypeOf(z))) * std.math.floatEps(Scalar(@TypeOf(z)))) {
                        res.im = (float.copysign(@as(Scalar(@TypeOf(z)), 0.5), z.im) * (float.ln2(Scalar(@TypeOf(z))) - float.log(float.abs(z.re))));
                    } else {
                        var r2: Scalar(@TypeOf(z)) = 0;
                        if (float.abs(z.re) >= std.math.floatEps(Scalar(@TypeOf(z))) * std.math.floatEps(Scalar(@TypeOf(z))))
                            r2 = z.re * z.re;

                        var num: Scalar(@TypeOf(z)) = z.im + 1;
                        num = r2 + num * num;

                        den = z.im - 1;
                        den = r2 + den * den;

                        const f: Scalar(@TypeOf(z)) = num / den;
                        if (f < 0.5) {
                            res.im = 0.25 * float.log(f);
                        } else {
                            num = 4 * z.im;
                            res.im = 0.25 * float.log1p(num / den);
                        }
                    }
                }

                if (float.abs(res.re) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                    const vresr: Scalar(@TypeOf(z)) = res.re * res.re;
                    std.mem.doNotOptimizeAway(vresr);
                }

                if (float.abs(res.im) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                    const vresi: Scalar(@TypeOf(z)) = res.im * res.im;
                    std.mem.doNotOptimizeAway(vresi);
                }

                return res;
            }
        },
        else => unreachable,
    }
}
