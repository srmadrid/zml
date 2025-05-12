const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn sinh(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return sinh(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return sinh(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const negate: bool = std.math.signbit(z.re);
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            const x: @TypeOf(z) = .{
                .re = float.abs(z.re),
                .im = z.im,
            };

            if (rcls >= classify.ZERO) {
                @branchHint(.likely);
                // Real part is finite.
                if (icls >= classify.ZERO) {
                    @branchHint(.likely);
                    var retval: @TypeOf(z) = undefined;
                    // Imaginary part is finite.
                    const t = cast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(z))) - 1) * float.ln2(Scalar(@TypeOf(z))), .{});

                    var sinix: Scalar(@TypeOf(z)) = undefined;
                    var cosix: Scalar(@TypeOf(z)) = undefined;
                    if (float.abs(x.im) > std.math.floatMin(Scalar(@TypeOf(z)))) {
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

                    if (float.abs(x.re) > cast(Scalar(@TypeOf(z)), t, .{})) {
                        const exp_t: Scalar(@TypeOf(z)) = float.exp(cast(Scalar(@TypeOf(z)), t, .{}));
                        var rx: Scalar(@TypeOf(z)) = float.abs(x.re);
                        if (std.math.signbit(x.re))
                            cosix = -cosix;

                        rx -= cast(Scalar(@TypeOf(z)), t, .{});
                        sinix *= exp_t / 2;
                        cosix *= exp_t / 2;
                        if (rx > cast(Scalar(@TypeOf(z)), t, .{})) {
                            rx -= cast(Scalar(@TypeOf(z)), t, .{});
                            sinix *= exp_t;
                            cosix *= exp_t;
                        }
                        if (rx > cast(Scalar(@TypeOf(z)), t, .{})) {
                            // Overflow (original real part of x > 3t).
                            retval.re = std.math.floatMax(Scalar(@TypeOf(z))) * cosix;
                            retval.im = std.math.floatMax(Scalar(@TypeOf(z))) * sinix;
                        } else {
                            const exp_val: Scalar(@TypeOf(z)) = float.exp(rx);
                            retval.re = exp_val * cosix;
                            retval.im = exp_val * sinix;
                        }
                    } else {
                        retval.re = float.sinh(x.re) * cosix;
                        retval.im = float.cosh(x.re) * sinix;
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
                    if (rcls == classify.ZERO) {
                        // Real part is 0.
                        return .{
                            .re = float.copysign(@as(Scalar(@TypeOf(z)), 0), @as(Scalar(@TypeOf(z)), if (negate) -1 else 1)),
                            .im = x.im - x.im,
                        };
                    } else {
                        return .{
                            .re = std.math.nan(Scalar(@TypeOf(z))),
                            .im = std.math.nan(Scalar(@TypeOf(z))),
                        };
                    }
                }
            } else if (rcls == classify.INFINITE) {
                var retval: @TypeOf(z) = undefined;
                // Real part is infinite.
                if (icls > classify.ZERO) {
                    @branchHint(.likely);
                    // Imaginary part is finite.
                    var sinix: Scalar(@TypeOf(z)) = undefined;
                    var cosix: Scalar(@TypeOf(z)) = undefined;
                    if (float.abs(x.im) > std.math.floatMin(Scalar(@TypeOf(z)))) {
                        @branchHint(.likely);
                        const tmp = float.sincos(x.im);
                        sinix = tmp.sinx;
                        cosix = tmp.cosx;
                    } else {
                        sinix = x.im;
                        cosix = 1;
                    }

                    retval.re = float.copysign(std.math.inf(Scalar(@TypeOf(z))), cosix);
                    retval.im = float.copysign(std.math.inf(Scalar(@TypeOf(z))), sinix);

                    if (negate)
                        retval.re = -retval.re;
                } else if (icls == classify.ZERO) {
                    // Imaginary part is 0.
                    retval.re = if (negate) -std.math.inf(Scalar(@TypeOf(z))) else std.math.inf(Scalar(@TypeOf(z)));
                    retval.im = x.im;
                } else {
                    retval.re = std.math.inf(Scalar(@TypeOf(z)));
                    retval.im = x.im - x.im;
                }

                return retval;
            } else {
                return .{
                    .re = std.math.nan(Scalar(@TypeOf(z))),
                    .im = if (x.im == 0) x.im else std.math.nan(Scalar(@TypeOf(z))),
                };
            }
        },
        else => unreachable,
    }
}
