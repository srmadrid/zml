const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn sin(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.sin: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    const negate: bool = std.math.signbit(zz.re);
    const rcls: u32 = classify.classify(zz.re);
    const icls: u32 = classify.classify(zz.im);

    const x: @TypeOf(zz) = .{ .re = float.abs(zz.re), .im = zz.im };

    if (icls >= classify.ZERO) {
        @branchHint(.likely);
        // Imaginary part is finite.
        if (rcls >= classify.ZERO) {
            @branchHint(.likely);
            // Real part is finite.
            const t: i32 = scast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(zz))) - 1) * float.ln2(Scalar(@TypeOf(zz))));

            var sinix: Scalar(@TypeOf(zz)) = undefined;
            var cosix: Scalar(@TypeOf(zz)) = undefined;
            if (x.re > std.math.floatMin(Scalar(@TypeOf(zz)))) {
                @branchHint(.likely);
                const tmp = float.sincos(x.re);
                sinix = tmp.sinx;
                cosix = tmp.cosx;
            } else {
                sinix = x.re;
                cosix = 1;
            }

            if (negate)
                sinix = -sinix;

            var retval: @TypeOf(zz) = undefined;
            if (float.abs(x.im) > scast(Scalar(@TypeOf(zz)), t)) {
                const exp_t: Scalar(@TypeOf(zz)) = float.exp(scast(Scalar(@TypeOf(zz)), t));
                var ix: Scalar(@TypeOf(zz)) = float.abs(x.im);
                if (std.math.signbit(x.im))
                    cosix = -cosix;
                ix -= scast(Scalar(@TypeOf(zz)), t);
                sinix *= exp_t / 2;
                cosix *= exp_t / 2;
                if (ix > scast(Scalar(@TypeOf(zz)), t)) {
                    ix -= scast(Scalar(@TypeOf(zz)), t);
                    sinix *= exp_t;
                    cosix *= exp_t;
                }

                if (ix > scast(Scalar(@TypeOf(zz)), t)) {
                    // Overflow (original imaginary part of x > 3t).
                    retval.re = std.math.floatMax(Scalar(@TypeOf(zz))) * sinix;
                    retval.im = std.math.floatMax(Scalar(@TypeOf(zz))) * cosix;
                } else {
                    const exp_val: Scalar(@TypeOf(zz)) = float.exp(ix);
                    retval.re = exp_val * sinix;
                    retval.im = exp_val * cosix;
                }
            } else {
                retval.re = float.cosh(x.im) * sinix;
                retval.im = float.sinh(x.im) * cosix;
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
            if (icls == classify.ZERO) {
                // Imaginary part is 0.
                return .{
                    .re = x.re - x.re,
                    .im = x.im,
                };
            } else {
                return .{
                    .re = std.math.nan(Scalar(@TypeOf(zz))),
                    .im = std.math.nan(Scalar(@TypeOf(zz))),
                };
            }
        }
    } else if (icls == classify.INFINITE) {
        // Imaginary part is infinite.
        if (rcls == classify.ZERO) {
            // Real part is 0.
            return .{
                .re = float.copysign(@as(Scalar(@TypeOf(zz)), 0), @as(Scalar(@TypeOf(zz)), if (negate) -1 else 1)),
                .im = x.im,
            };
        } else if (rcls > classify.ZERO) {
            // Real part is finite.
            var sinix: Scalar(@TypeOf(zz)) = undefined;
            var cosix: Scalar(@TypeOf(zz)) = undefined;
            if (x.re > std.math.floatMin(Scalar(@TypeOf(zz)))) {
                @branchHint(.likely);
                const tmp = float.sincos(x.re);
                sinix = tmp.sinx;
                cosix = tmp.cosx;
            } else {
                sinix = x.re;
                cosix = 1;
            }

            var retval: @TypeOf(zz) = .{
                .re = float.copysign(std.math.inf(Scalar(@TypeOf(zz))), sinix),
                .im = float.copysign(std.math.inf(Scalar(@TypeOf(zz))), cosix),
            };

            if (negate)
                retval.re = -retval.re;

            if (std.math.signbit(x.im))
                retval.im = -retval.im;

            return retval;
        } else {
            return .{
                .re = x.re - x.re,
                .im = std.math.inf(Scalar(@TypeOf(zz))),
            };
        }
    } else {
        return .{
            .re = if (rcls == classify.ZERO) float.copysign(@as(Scalar(@TypeOf(zz)), 0), @as(Scalar(@TypeOf(zz)), if (negate) -1 else 1)) else std.math.nan(Scalar(@TypeOf(zz))),
            .im = std.math.nan(Scalar(@TypeOf(zz))),
        };
    }
}
