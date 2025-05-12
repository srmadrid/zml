const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn asinh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return asinh(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return asinh(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            if (rcls <= classify.INFINITE or icls <= classify.INFINITE) {
                var res: @TypeOf(z) = undefined;
                if (icls == classify.INFINITE) {
                    res.re = float.copysign(std.math.inf(Scalar(@TypeOf(z))), z.re);

                    if (rcls == classify.NAN) {
                        res.im = std.math.nan(Scalar(@TypeOf(z)));
                    } else {
                        res.im = float.copysign(if (rcls >= classify.ZERO) float.pi_2(Scalar(@TypeOf(z))) else float.pi_4(Scalar(@TypeOf(z))), z.im);
                    }
                } else if (rcls <= classify.INFINITE) {
                    res.re = z.re;
                    if ((rcls == classify.INFINITE and icls >= classify.ZERO) or (rcls == classify.NAN and icls == classify.ZERO)) {
                        res.im = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.im);
                    } else {
                        res.im = std.math.nan(Scalar(@TypeOf(z)));
                    }
                } else {
                    res.re = std.math.nan(Scalar(@TypeOf(z)));
                    res.im = std.math.nan(Scalar(@TypeOf(z)));
                }

                return res;
            } else if (rcls == classify.ZERO and icls == classify.ZERO) {
                return z;
            } else {
                return kernel_asinh(z, 0);
            }
        },
        else => unreachable,
    }
}

pub fn kernel_asinh(x: anytype, adj: i32) @TypeOf(x) {
    // Avoid cancellation by reducing to the first quadrant.
    const rx: Scalar(@TypeOf(x)) = float.abs(x.re);
    const ix: Scalar(@TypeOf(x)) = float.abs(x.im);

    var res: @TypeOf(x) = undefined;
    if (rx >= 1 / std.math.floatEps(Scalar(@TypeOf(x))) or ix >= 1 / std.math.floatEps(Scalar(@TypeOf(x)))) {
        // For large x in the first quadrant, x + csqrt (1 + x * x)
        // is sufficiently close to 2 * x to make no significant
        // difference to the result; avoid possible overflow from
        // the squaring and addition.
        var y: @TypeOf(x) = .{
            .re = rx,
            .im = ix,
        };

        if (adj != 0) {
            const t: Scalar(@TypeOf(x)) = y.re;
            y.re = float.copysign(y.im, x.im);
            y.im = t;
        }

        res = cfloat.log(y);
        res.re += float.ln2(Scalar(@TypeOf(x)));
    } else if (rx >= 0.5 and ix < std.math.floatEps(Scalar(@TypeOf(x))) / 8) {
        const s: Scalar(@TypeOf(x)) = float.hypot(1, rx);

        res.re = float.log(rx + s);
        if (adj != 0) {
            res.im = float.atan2(s, x.im);
        } else {
            res.im = float.atan2(ix, s);
        }
    } else if (rx < std.math.floatEps(Scalar(@TypeOf(x))) / 8 and ix >= 1.5) {
        const s: Scalar(@TypeOf(x)) = float.sqrt((ix + 1) * (ix - 1));

        res.re = float.log(ix + s);
        if (adj != 0) {
            res.im = float.atan2(rx, float.copysign(s, x.im));
        } else {
            res.im = float.atan2(s, rx);
        }
    } else if (ix > 1 and ix < 1.5 and rx < 0.5) {
        if (rx < std.math.floatEps(Scalar(@TypeOf(x))) * std.math.floatEps(Scalar(@TypeOf(x)))) {
            const ix2m1: Scalar(@TypeOf(x)) = (ix + 1) * (ix - 1);
            const s: Scalar(@TypeOf(x)) = float.sqrt(ix2m1);

            res.re = float.log1p(2 * (ix2m1 + ix * s)) / 2;
            if (adj != 0) {
                res.im = float.atan2(rx, float.copysign(s, x.im));
            } else {
                res.im = float.atan2(s, rx);
            }
        } else {
            const ix2m1: Scalar(@TypeOf(x)) = (ix + 1) * (ix - 1);
            const rx2: Scalar(@TypeOf(x)) = rx * rx;
            const f: Scalar(@TypeOf(x)) = rx2 * (2 + rx2 + 2 * ix * ix);
            const d: Scalar(@TypeOf(x)) = float.sqrt(ix2m1 * ix2m1 + f);
            const dp: Scalar(@TypeOf(x)) = d + ix2m1;
            const dm: Scalar(@TypeOf(x)) = f / dp;
            const r1: Scalar(@TypeOf(x)) = float.sqrt((dm + rx2) / 2);
            const r2: Scalar(@TypeOf(x)) = rx * ix / r1;

            res.re = float.log1p(rx2 + dp + 2 * (rx * r1 + ix * r2)) / 2;
            if (adj != 0) {
                res.im = float.atan2(rx + r1, float.copysign(ix + r2, x.im));
            } else {
                res.im = float.atan2(ix + r2, rx + r1);
            }
        }
    } else if (ix == 1 and rx < 0.5) {
        if (rx < std.math.floatEps(Scalar(@TypeOf(x))) / 8) {
            res.re = float.log1p(2 * (rx + float.sqrt(rx))) / 2;
            if (adj != 0) {
                res.im = float.atan2(float.sqrt(rx), float.copysign(@as(Scalar(@TypeOf(x)), 1), x.im));
            } else {
                res.im = float.atan2(1, float.sqrt(rx));
            }
        } else {
            const d: Scalar(@TypeOf(x)) = rx * float.sqrt(4 + rx * rx);
            const s1: Scalar(@TypeOf(x)) = float.sqrt((d + rx * rx) / 2);
            const s2: Scalar(@TypeOf(x)) = float.sqrt((d - rx * rx) / 2);

            res.re = float.log1p(rx * rx + d + 2 * (rx * s1 + s2)) / 2;
            if (adj != 0) {
                res.im = float.atan2(rx + s1, float.copysign(1 + s2, x.im));
            } else {
                res.im = float.atan2(1 + s2, rx + s1);
            }
        }
    } else if (ix < 1 and rx < 0.5) {
        if (ix >= std.math.floatEps(Scalar(@TypeOf(x)))) {
            if (rx < std.math.floatEps(Scalar(@TypeOf(x))) * std.math.floatEps(Scalar(@TypeOf(x)))) {
                const onemix2: Scalar(@TypeOf(x)) = (1 + ix) * (1 - ix);
                const s: Scalar(@TypeOf(x)) = float.sqrt(onemix2);

                res.re = float.log1p(2 * rx / s) / 2;
                if (adj != 0) {
                    res.im = float.atan2(s, x.im);
                } else {
                    res.im = float.atan2(ix, s);
                }
            } else {
                const onemix2: Scalar(@TypeOf(x)) = (1 + ix) * (1 - ix);
                const rx2: Scalar(@TypeOf(x)) = rx * rx;
                const f: Scalar(@TypeOf(x)) = rx2 * (2 + rx2 + 2 * ix * ix);
                const d: Scalar(@TypeOf(x)) = float.sqrt(onemix2 * onemix2 + f);
                const dp: Scalar(@TypeOf(x)) = d + onemix2;
                const dm: Scalar(@TypeOf(x)) = f / dp;
                const r1: Scalar(@TypeOf(x)) = float.sqrt((dp + rx2) / 2);
                const r2: Scalar(@TypeOf(x)) = rx * ix / r1;

                res.re = float.log1p(rx2 + dm + 2 * (rx * r1 + ix * r2)) / 2;
                if (adj != 0) {
                    res.im = float.atan2(rx + r1, float.copysign(ix + r2, x.im));
                } else {
                    res.im = float.atan2(ix + r2, rx + r1);
                }
            }
        } else {
            const s: Scalar(@TypeOf(x)) = float.hypot(1, rx);

            res.re = float.log1p(2 * rx * (rx + s)) / 2;
            if (adj != 0) {
                res.im = float.atan2(s, x.im);
            } else {
                res.im = float.atan2(ix, s);
            }
        }

        if (res.re < std.math.floatMin(Scalar(@TypeOf(x)))) {
            const vresr: Scalar(@TypeOf(x)) = res.re * res.re;
            std.mem.doNotOptimizeAway(vresr);
        }
    } else {
        var y: @TypeOf(x) = .{
            .re = (rx - ix) * (rx + ix) + 1,
            .im = 2 * rx * ix,
        };

        y = cfloat.sqrt(y);

        y.re += rx;
        y.im += ix;

        if (adj != 0) {
            const t: Scalar(@TypeOf(x)) = y.re;
            y.re = float.copysign(y.im, x.im);
            y.im = t;
        }

        res = cfloat.log(y);
    }

    // Give results the correct sign for the original argument.
    return .{
        .re = float.copysign(res.re, x.re),
        .im = float.copysign(res.im, if (adj != 0) 1 else x.im),
    };
}
