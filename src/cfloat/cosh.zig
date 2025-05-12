const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn cosh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return cosh(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return cosh(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            if (rcls >= classify.ZERO) {
                @branchHint(.likely);
                // Real part is finite.
                if (icls >= classify.ZERO) {
                    // Imaginary part is finite.
                    const t: i32 = cast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(z))) - 1) * float.ln2(Scalar(@TypeOf(z))), .{});

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

                    var retval: @TypeOf(z) = undefined;
                    if (float.abs(z.re) > cast(Scalar(@TypeOf(z)), t, .{})) {
                        const exp_t: Scalar(@TypeOf(z)) = float.exp(cast(Scalar(@TypeOf(z)), t, .{}));
                        var rz: Scalar(@TypeOf(z)) = float.abs(z.re);
                        if (std.math.signbit(z.re))
                            siniz = -siniz;
                        rz -= cast(Scalar(@TypeOf(z)), t, .{});
                        siniz *= exp_t / 2;
                        cosiz *= exp_t / 2;
                        if (rz > cast(Scalar(@TypeOf(z)), t, .{})) {
                            rz -= cast(Scalar(@TypeOf(z)), t, .{});
                            siniz *= exp_t;
                            cosiz *= exp_t;
                        }

                        if (rz > cast(Scalar(@TypeOf(z)), t, .{})) {
                            // Overflow (original real part of z > 3t).
                            retval.re = std.math.floatMax(Scalar(@TypeOf(z))) * cosiz;
                            retval.im = std.math.floatMax(Scalar(@TypeOf(z))) * siniz;
                        } else {
                            const exp_val: Scalar(@TypeOf(z)) = float.exp(rz);
                            retval.re = exp_val * cosiz;
                            retval.im = exp_val * siniz;
                        }
                    } else {
                        retval.re = float.cosh(z.re) * cosiz;
                        retval.im = float.sinh(z.re) * siniz;
                    }

                    if (float.abs(retval.re) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                        const vretvalr: Scalar(@TypeOf(z)) = retval.re * retval.re;
                        std.mem.doNotOptimizeAway(vretvalr);
                    }

                    if (float.abs(retval.im) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                        const vretvali: Scalar(@TypeOf(z)) = retval.im * retval.im;
                        std.mem.doNotOptimizeAway(vretvali);
                    }

                    return retval;
                } else {
                    return .{
                        .re = z.im - z.im,
                        .im = if (z.re == 0) 0 else std.math.nan(Scalar(@TypeOf(z))),
                    };
                }
            } else if (rcls == classify.INFINITE) {
                // Real part is infinite.
                if (icls > classify.ZERO) {
                    @branchHint(.likely);
                    // Imaginary part is finite.
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
                        .re = float.copysign(std.math.inf(Scalar(@TypeOf(z))), cosiz),
                        .im = (float.copysign(std.math.inf(Scalar(@TypeOf(z))), siniz) * float.copysign(@as(Scalar(@TypeOf(z)), 1), z.re)),
                    };
                } else if (icls == classify.ZERO) {
                    // Imaginary part is 0.
                    return .{
                        .re = std.math.inf(Scalar(@TypeOf(z))),
                        .im = z.im * float.copysign(@as(Scalar(@TypeOf(z)), 1), z.re),
                    };
                } else {
                    return .{
                        .re = std.math.inf(Scalar(@TypeOf(z))),
                        .im = z.im - z.im,
                    };
                }
            } else {
                return .{
                    .re = std.math.nan(Scalar(@TypeOf(z))),
                    .im = if (z.im == 0) z.im else std.math.nan(Scalar(@TypeOf(z))),
                };
            }
        },
        else => unreachable,
    }
}
