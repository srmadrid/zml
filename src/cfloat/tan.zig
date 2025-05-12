const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn tan(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return tan(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return tan(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            if (!std.math.isFinite(z.re) or !std.math.isFinite(z.im)) {
                @branchHint(.unlikely);
                var res: @TypeOf(z) = undefined;
                if (std.math.isInf(z.im)) {
                    if (std.math.isFinite(z.re) and float.abs(z.re) > 1) {
                        const tmp = float.sincos(z.re);
                        const sinrz: Scalar(@TypeOf(z)) = tmp.sinx;
                        const cosrz: Scalar(@TypeOf(z)) = tmp.cosx;

                        res.re = float.copysign(@as(Scalar(@TypeOf(z)), 0), sinrz * cosrz);
                    } else {
                        res.re = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.re);
                    }

                    res.im = float.copysign(@as(Scalar(@TypeOf(z)), 1), z.im);
                } else if (z.re == 0) {
                    res = z;
                } else {
                    res.re = std.math.nan(Scalar(@TypeOf(z)));
                    if (z.im == 0) {
                        res.im = z.im;
                    } else {
                        res.im = std.math.nan(Scalar(@TypeOf(z)));
                    }
                }

                return res;
            } else {
                const t: i32 = cast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(z))) - 1) * float.ln2(Scalar(@TypeOf(z))) / 2, .{});

                // tan(z+iy) = (sin(2z) + i*sinh(2y))/(cos(2z) + cosh(2y))
                // = (sin(z)*cos(z) + i*sinh(y)*cosh(y)/(cos(z)^2 + sinh(y)^2).
                var sinrz: Scalar(@TypeOf(z)) = undefined;
                var cosrz: Scalar(@TypeOf(z)) = undefined;
                if (float.abs(z.re) > std.math.floatMin(Scalar(@TypeOf(z)))) {
                    @branchHint(.likely);
                    const tmp = float.sincos(z.re);
                    sinrz = tmp.sinx;
                    cosrz = tmp.cosx;
                } else {
                    sinrz = z.re;
                    cosrz = 1;
                }

                var res: @TypeOf(z) = undefined;
                if (float.abs(z.im) > cast(Scalar(@TypeOf(z)), t, .{})) {
                    // Avoid intermediate overflow when the real part of the
                    // result may be subnormal.  Ignoring negligible terms, the
                    // imaginary part is +/- 1, the real part is
                    // sin(z)*cos(z)/sinh(y)^2 = 4*sin(z)*cos(z)/exp(2y).
                    const exp_2t: Scalar(@TypeOf(z)) = float.exp(2 * cast(Scalar(@TypeOf(z)), t, .{}));

                    res.im = float.copysign(@as(Scalar(@TypeOf(z)), 1), z.im);
                    res.re = 4 * sinrz * cosrz;
                    const zz: Scalar(@TypeOf(z)) = float.abs(z.im) - cast(Scalar(@TypeOf(z)), t, .{});
                    res.re /= exp_2t;
                    if (zz > cast(Scalar(@TypeOf(z)), t, .{})) {
                        // Underflow (original imaginary part of z has absolute
                        // value > 2t).
                        res.re /= exp_2t;
                    } else {
                        res.re /= float.exp(2 * zz);
                    }
                } else {
                    var sinhiz: Scalar(@TypeOf(z)) = undefined;
                    var coshiz: Scalar(@TypeOf(z)) = undefined;
                    if (float.abs(z.im) > std.math.floatMin(Scalar(@TypeOf(z)))) {
                        sinhiz = float.sinh(z.im);
                        coshiz = float.cosh(z.im);
                    } else {
                        sinhiz = z.im;
                        coshiz = 1;
                    }

                    var den: Scalar(@TypeOf(z)) = undefined;
                    if (float.abs(sinhiz) > float.abs(cosrz) * std.math.floatEps(Scalar(@TypeOf(z)))) {
                        den = cosrz * cosrz + sinhiz * sinhiz;
                    } else {
                        den = cosrz * cosrz;
                    }

                    res.re = sinrz * cosrz / den;
                    res.im = sinhiz * coshiz / den;
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
