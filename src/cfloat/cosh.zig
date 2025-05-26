const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn cosh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.cosh: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    const rcls: u32 = classify.classify(zz.re);
    const icls: u32 = classify.classify(zz.im);

    if (rcls >= classify.ZERO) {
        @branchHint(.likely);
        // Real part is finite.
        if (icls >= classify.ZERO) {
            // Imaginary part is finite.
            const t: i32 = scast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(zz))) - 1) * float.ln2(Scalar(@TypeOf(zz))));

            var sinizz: Scalar(@TypeOf(zz)) = undefined;
            var cosizz: Scalar(@TypeOf(zz)) = undefined;
            if (float.abs(zz.im) > std.math.floatMin(Scalar(@TypeOf(zz)))) {
                @branchHint(.likely);
                const tmp = float.sincos(zz.im);
                sinizz = tmp.sinx;
                cosizz = tmp.cosx;
            } else {
                sinizz = zz.im;
                cosizz = 1;
            }

            var retval: @TypeOf(zz) = undefined;
            if (float.abs(zz.re) > scast(Scalar(@TypeOf(zz)), t)) {
                const exp_t: Scalar(@TypeOf(zz)) = float.exp(scast(Scalar(@TypeOf(zz)), t));
                var rzz: Scalar(@TypeOf(zz)) = float.abs(zz.re);
                if (std.math.signbit(zz.re))
                    sinizz = -sinizz;
                rzz -= scast(Scalar(@TypeOf(zz)), t);
                sinizz *= exp_t / 2;
                cosizz *= exp_t / 2;
                if (rzz > scast(Scalar(@TypeOf(zz)), t)) {
                    rzz -= scast(Scalar(@TypeOf(zz)), t);
                    sinizz *= exp_t;
                    cosizz *= exp_t;
                }

                if (rzz > scast(Scalar(@TypeOf(zz)), t)) {
                    // Overflow (original real part of zz > 3t).
                    retval.re = std.math.floatMax(Scalar(@TypeOf(zz))) * cosizz;
                    retval.im = std.math.floatMax(Scalar(@TypeOf(zz))) * sinizz;
                } else {
                    const exp_val: Scalar(@TypeOf(zz)) = float.exp(rzz);
                    retval.re = exp_val * cosizz;
                    retval.im = exp_val * sinizz;
                }
            } else {
                retval.re = float.cosh(zz.re) * cosizz;
                retval.im = float.sinh(zz.re) * sinizz;
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
            return .{
                .re = zz.im - zz.im,
                .im = if (zz.re == 0) 0 else std.math.nan(Scalar(@TypeOf(zz))),
            };
        }
    } else if (rcls == classify.INFINITE) {
        // Real part is infinite.
        if (icls > classify.ZERO) {
            @branchHint(.likely);
            // Imaginary part is finite.
            var sinizz: Scalar(@TypeOf(zz)) = undefined;
            var cosizz: Scalar(@TypeOf(zz)) = undefined;
            if (float.abs(zz.im) > std.math.floatMin(Scalar(@TypeOf(zz)))) {
                @branchHint(.likely);
                const tmp = float.sincos(zz.im);
                sinizz = tmp.sinx;
                cosizz = tmp.cosx;
            } else {
                sinizz = zz.im;
                cosizz = 1;
            }

            return .{
                .re = float.copysign(std.math.inf(Scalar(@TypeOf(zz))), cosizz),
                .im = (float.copysign(std.math.inf(Scalar(@TypeOf(zz))), sinizz) * float.copysign(@as(Scalar(@TypeOf(zz)), 1), zz.re)),
            };
        } else if (icls == classify.ZERO) {
            // Imaginary part is 0.
            return .{
                .re = std.math.inf(Scalar(@TypeOf(zz))),
                .im = zz.im * float.copysign(@as(Scalar(@TypeOf(zz)), 1), zz.re),
            };
        } else {
            return .{
                .re = std.math.inf(Scalar(@TypeOf(zz))),
                .im = zz.im - zz.im,
            };
        }
    } else {
        return .{
            .re = std.math.nan(Scalar(@TypeOf(zz))),
            .im = if (zz.im == 0) zz.im else std.math.nan(Scalar(@TypeOf(zz))),
        };
    }
}
