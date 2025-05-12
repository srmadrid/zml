const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const x2y2m1 = @import("../float/x2y2m1.zig").x2y2m1;
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn atanh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return atanh(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return atanh(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            if (rcls <= classify.INFINITE or icls <= classify.INFINITE) {
                @branchHint(.unlikely);
                var res: @TypeOf(z) = undefined;
                if (icls == classify.INFINITE) {
                    res.re = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.re);
                    res.im = float.copysign(float.pi_2(Scalar(@TypeOf(z))), z.im);
                } else if (rcls == classify.INFINITE or rcls == classify.ZERO) {
                    res.re = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.re);
                    if (icls >= classify.ZERO) {
                        res.im = float.copysign(float.pi_2(Scalar(@TypeOf(z))), z.im);
                    } else {
                        res.im = std.math.nan(Scalar(@TypeOf(z)));
                    }
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
                    res.im = float.copysign(float.pi_2(Scalar(@TypeOf(z))), z.im);
                    if (float.abs(z.im) <= 1) {
                        res.re = 1 / z.re;
                    } else if (float.abs(z.re) <= 1) {
                        res.re = z.re / z.im / z.im;
                    } else {
                        const h: Scalar(@TypeOf(z)) = float.hypot(z.re / 2, z.im / 2);
                        res.re = z.re / h / h / 4;
                    }
                } else {
                    if (float.abs(z.re) == 1 and float.abs(z.im) < std.math.floatEps(Scalar(@TypeOf(z))) * std.math.floatEps(Scalar(@TypeOf(z)))) {
                        res.re = (float.copysign(@as(Scalar(@TypeOf(z)), 0.5), z.re) * (float.ln2(Scalar(@TypeOf(z))) - float.log(float.abs(z.im))));
                    } else {
                        var I2: Scalar(@TypeOf(z)) = 0;
                        if (float.abs(z.im) >= std.math.floatEps(Scalar(@TypeOf(z))) * std.math.floatEps(Scalar(@TypeOf(z))))
                            I2 = z.im * z.im;

                        var num: Scalar(@TypeOf(z)) = 1 + z.re;
                        num = I2 + num * num;

                        var den: Scalar(@TypeOf(z)) = 1 - z.re;
                        den = I2 + den * den;

                        const f: Scalar(@TypeOf(z)) = num / den;
                        if (f < 0.5) {
                            res.re = 0.25 * float.log(f);
                        } else {
                            num = 4 * z.re;
                            res.re = 0.25 * float.log1p(num / den);
                        }
                    }
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

                    res.im = 0.5 * float.atan2(2 * z.im, den);
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
