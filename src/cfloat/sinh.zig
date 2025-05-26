const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn sinh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.sinh: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    const negate: bool = std.math.signbit(zz.re);
    const rcls: u32 = classify.classify(zz.re);
    const icls: u32 = classify.classify(zz.im);

    const x: @TypeOf(zz) = .{
        .re = float.abs(zz.re),
        .im = zz.im,
    };

    if (rcls >= classify.ZERO) {
        @branchHint(.likely);
        // Real part is finite.
        if (icls >= classify.ZERO) {
            @branchHint(.likely);
            var retval: @TypeOf(zz) = undefined;
            // Imaginary part is finite.
            const t = scast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(zz))) - 1) * float.ln2(Scalar(@TypeOf(zz))));

            var sinix: Scalar(@TypeOf(zz)) = undefined;
            var cosix: Scalar(@TypeOf(zz)) = undefined;
            if (float.abs(x.im) > std.math.floatMin(Scalar(@TypeOf(zz)))) {
                @branchHint(.likely);
                const tmp = float.sincos(x.im);
                sinix = tmp.sinx;
                cosix = tmp.cosx;
            } else {
                sinix = x.im;
                cosix = 1;
            }

            if (negate)
                cosix = -cosix;

            if (float.abs(x.re) > scast(Scalar(@TypeOf(zz)), t)) {
                const exp_t: Scalar(@TypeOf(zz)) = float.exp(scast(Scalar(@TypeOf(zz)), t));
                var rx: Scalar(@TypeOf(zz)) = float.abs(x.re);
                if (std.math.signbit(x.re))
                    cosix = -cosix;

                rx -= scast(Scalar(@TypeOf(zz)), t);
                sinix *= exp_t / 2;
                cosix *= exp_t / 2;
                if (rx > scast(Scalar(@TypeOf(zz)), t)) {
                    rx -= scast(Scalar(@TypeOf(zz)), t);
                    sinix *= exp_t;
                    cosix *= exp_t;
                }
                if (rx > scast(Scalar(@TypeOf(zz)), t)) {
                    // Overflow (original real part of x > 3t).
                    retval.re = std.math.floatMax(Scalar(@TypeOf(zz))) * cosix;
                    retval.im = std.math.floatMax(Scalar(@TypeOf(zz))) * sinix;
                } else {
                    const exp_val: Scalar(@TypeOf(zz)) = float.exp(rx);
                    retval.re = exp_val * cosix;
                    retval.im = exp_val * sinix;
                }
            } else {
                retval.re = float.sinh(x.re) * cosix;
                retval.im = float.cosh(x.re) * sinix;
            }

            if (float.abs(retval.re) < std.math.floatMin(Scalar(@TypeOf(zz)))) {
                const vretvalr: Scalar(@TypeOf(zz)) = retval.re * retval.re;
                std.mem.doNotOptimizeAway(vretvalr);
            }

            if (float.abs(retval.im) < std.math.floatMin(Scalar(@TypeOf(zz)))) {
                const vretvali: Scalar(@TypeOf(zz)) = retval.im * retval.im;
                std.mem.doNotOptimizeAway(vretvali);
            }

            return retval;
        } else {
            if (rcls == classify.ZERO) {
                // Real part is 0.
                return .{
                    .re = float.copysign(@as(Scalar(@TypeOf(zz)), 0), @as(Scalar(@TypeOf(zz)), if (negate) -1 else 1)),
                    .im = x.im - x.im,
                };
            } else {
                return .{
                    .re = std.math.nan(Scalar(@TypeOf(zz))),
                    .im = std.math.nan(Scalar(@TypeOf(zz))),
                };
            }
        }
    } else if (rcls == classify.INFINITE) {
        var retval: @TypeOf(zz) = undefined;
        // Real part is infinite.
        if (icls > classify.ZERO) {
            @branchHint(.likely);
            // Imaginary part is finite.
            var sinix: Scalar(@TypeOf(zz)) = undefined;
            var cosix: Scalar(@TypeOf(zz)) = undefined;
            if (float.abs(x.im) > std.math.floatMin(Scalar(@TypeOf(zz)))) {
                @branchHint(.likely);
                const tmp = float.sincos(x.im);
                sinix = tmp.sinx;
                cosix = tmp.cosx;
            } else {
                sinix = x.im;
                cosix = 1;
            }

            retval.re = float.copysign(std.math.inf(Scalar(@TypeOf(zz))), cosix);
            retval.im = float.copysign(std.math.inf(Scalar(@TypeOf(zz))), sinix);

            if (negate)
                retval.re = -retval.re;
        } else if (icls == classify.ZERO) {
            // Imaginary part is 0.
            retval.re = if (negate) -std.math.inf(Scalar(@TypeOf(zz))) else std.math.inf(Scalar(@TypeOf(zz)));
            retval.im = x.im;
        } else {
            retval.re = std.math.inf(Scalar(@TypeOf(zz)));
            retval.im = x.im - x.im;
        }

        return retval;
    } else {
        return .{
            .re = std.math.nan(Scalar(@TypeOf(zz))),
            .im = if (x.im == 0) x.im else std.math.nan(Scalar(@TypeOf(zz))),
        };
    }
}
