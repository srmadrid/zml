const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const scast = types.scast;

pub fn exp(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.acos: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    const rcls: u32 = classify.classify(z.re);
    const icls: u32 = classify.classify(z.im);

    if (rcls >= classify.ZERO) {
        @branchHint(.likely);
        // Real part is finite.
        if (icls >= classify.ZERO) {
            @branchHint(.likely);
            // Imaginary part is finite.
            const t: Scalar(@TypeOf(z)) = scast(Scalar(@TypeOf(z)), scast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(z))) - 1) * float.ln2(Scalar(@TypeOf(z)))));
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

            var zt: Scalar(@TypeOf(z)) = z.re;
            if (zt > t) {
                const exp_t: Scalar(@TypeOf(z)) = float.exp(t);
                zt -= t;
                siniz *= exp_t;
                cosiz *= exp_t;
                if (zt > t) {
                    zt -= t;
                    siniz *= exp_t;
                    cosiz *= exp_t;
                }
            }

            if (zt > t) {
                // Overflow (original real part of z > 3t).
                const retre: Scalar(@TypeOf(z)) = std.math.floatMax(Scalar(@TypeOf(z))) * cosiz;
                if (float.abs(retre) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                    const vretre: Scalar(@TypeOf(z)) = retre * retre;
                    std.mem.doNotOptimizeAway(vretre);
                }
                const retim: Scalar(@TypeOf(z)) = std.math.floatMax(Scalar(@TypeOf(z))) * siniz;
                if (float.abs(retim) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                    const vretim: Scalar(@TypeOf(z)) = retim * retim;
                    std.mem.doNotOptimizeAway(vretim);
                }

                return .{
                    .re = retre,
                    .im = retim,
                };
            } else {
                const exp_val: Scalar(@TypeOf(z)) = float.exp(zt);

                const retre: Scalar(@TypeOf(z)) = exp_val * cosiz;
                if (float.abs(retre) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                    const vretre: Scalar(@TypeOf(z)) = retre * retre;
                    std.mem.doNotOptimizeAway(vretre);
                }
                const retim: Scalar(@TypeOf(z)) = exp_val * siniz;
                if (float.abs(retim) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                    const vretim: Scalar(@TypeOf(z)) = retim * retim;
                    std.mem.doNotOptimizeAway(vretim);
                }

                return .{
                    .re = retre,
                    .im = retim,
                };
            }
        } else {
            // If the imaginary part is +-inf or NaN and the real part
            // is not +-inf the result is NaN + iNaN.
            return .{
                .re = std.math.nan(Scalar(@TypeOf(z))),
                .im = std.math.nan(Scalar(@TypeOf(z))),
            };
        }
    } else if (rcls == classify.INFINITE) {
        @branchHint(.likely);
        // Real part is infinite.
        if (icls >= classify.ZERO) {
            @branchHint(.likely);
            // Imaginary part is finite.
            const value: Scalar(@TypeOf(z)) = if (std.math.signbit(z.re)) 0 else std.math.inf(Scalar(@TypeOf(z)));

            if (icls == classify.ZERO) {
                // Imaginary part is 0.0.
                return .{
                    .re = value,
                    .im = z.im,
                };
            } else {
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

                return .{
                    .re = float.copysign(value, cosiz),
                    .im = float.copysign(value, siniz),
                };
            }
        } else if (!std.math.signbit(z.re)) {
            return .{
                .re = std.math.floatMax(Scalar(@TypeOf(z))),
                .im = z.im - z.im,
            };
        } else {
            return .{
                .re = 0,
                .im = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.im),
            };
        }
    } else {
        // If the real part is NaN the result is NaN + iNaN unless the
        // imaginary part is zero.
        if (icls == classify.ZERO) {
            return .{
                .re = std.math.nan(Scalar(@TypeOf(z))),
                .im = z.im,
            };
        } else {
            return .{
                .re = std.math.nan(Scalar(@TypeOf(z))),
                .im = std.math.nan(Scalar(@TypeOf(z))),
            };
        }
    }
}
