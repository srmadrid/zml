const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn tanh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return tanh(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return tanh(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            if (!std.math.isFinite(z.re) or !std.math.isFinite(z.im)) {
                @branchHint(.unlikely);
                var res: @TypeOf(z) = undefined;
                if (std.math.isInf(z.re)) {
                    res.re = float.copysign(@as(Scalar(@TypeOf(z)), 1), z.re);
                    if (std.math.isFinite(z.im) and float.abs(z.im) > 1) {
                        const tmp = float.sincos(z.im);
                        res.im = float.copysign(@as(Scalar(@TypeOf(z)), 0), tmp.sinx * tmp.cosx);
                    } else {
                        res.im = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.im);
                    }
                } else if (z.im == 0) {
                    res = z;
                } else {
                    if (z.re == 0) {
                        res.re = z.re;
                    } else {
                        res.re = std.math.nan(Scalar(@TypeOf(z)));
                    }

                    res.im = std.math.nan(Scalar(@TypeOf(z)));
                }

                return res;
            } else {
                const t: i32 = cast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(z))) - 1) * float.ln2(Scalar(@TypeOf(z))) / 2, .{});

                // tanh(z+iy) = (sinh(2z) + i*sin(2y))/(cosh(2z) + cos(2y))
                // = (sinh(z)*cosh(z) + i*sin(y)*cos(y))/(sinh(z)^2 + cos(y)^2).
                var siniz: Scalar(@TypeOf(z)) = undefined;
                var cosiz: Scalar(@TypeOf(z)) = undefined;
                if (float.abs(z.im) > std.math.floatMin(Scalar(@TypeOf(z)))) {
                    @branchHint(.likely);
                    const tmp = float.sincos(z.im);
                    siniz = tmp.sinx;
                    cosiz = tmp.cosx;
                } else {
                    siniz = z.im;
                    cosiz = 1;
                }

                var res: @TypeOf(z) = undefined;
                if (float.abs(z.re) > cast(Scalar(@TypeOf(z)), t, .{})) {
                    // Avoid intermediate overflow when the imaginary part of
                    // the result may be subnormal.  Ignoring negligible terms,
                    // the real part is +/- 1, the imaginary part is
                    // sin(y)*cos(y)/sinh(z)^2 = 4*sin(y)*cos(y)/exp(2z).
                    const exp_2t: Scalar(@TypeOf(z)) = float.exp(2 * cast(Scalar(@TypeOf(z)), t, .{}));

                    res.re = float.copysign(@as(Scalar(@TypeOf(z)), 1), z.re);
                    res.im = 4 * siniz * cosiz;
                    const zz: Scalar(@TypeOf(z)) = float.abs(z.re) - cast(Scalar(@TypeOf(z)), t, .{});
                    res.im /= exp_2t;
                    if (zz > cast(Scalar(@TypeOf(z)), t, .{})) {
                        // Underflow (original real part of z has absolute value > 2t).
                        res.im /= exp_2t;
                    } else {
                        res.im /= float.exp(2 * zz);
                    }
                } else {
                    var sinhrz: Scalar(@TypeOf(z)) = undefined;
                    var coshrz: Scalar(@TypeOf(z)) = undefined;
                    if (float.abs(z.re) > std.math.floatMin(Scalar(@TypeOf(z)))) {
                        sinhrz = float.sinh(z.re);
                        coshrz = float.cosh(z.re);
                    } else {
                        sinhrz = z.re;
                        coshrz = 1;
                    }

                    var den: Scalar(@TypeOf(z)) = undefined;
                    if (float.abs(sinhrz) > float.abs(cosiz) * std.math.floatEps(Scalar(@TypeOf(z)))) {
                        den = sinhrz * sinhrz + cosiz * cosiz;
                    } else {
                        den = cosiz * cosiz;
                    }

                    res.re = sinhrz * coshrz / den;
                    res.im = siniz * cosiz / den;
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
