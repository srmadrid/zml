const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn sin(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return sin(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return sin(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const negate: bool = std.math.signbit(z.re);
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            const x: @TypeOf(z) = .{ .re = float.abs(z.re), .im = z.im };

            if (icls >= classify.ZERO) {
                @branchHint(.likely);
                // Imaginary part is finite.
                if (rcls >= classify.ZERO) {
                    @branchHint(.likely);
                    // Real part is finite.
                    const t: i32 = cast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(z))) - 1) * float.ln2(Scalar(@TypeOf(z))), .{});

                    var sinix: Scalar(@TypeOf(z)) = undefined;
                    var cosix: Scalar(@TypeOf(z)) = undefined;
                    if (x.re > std.math.floatMin(Scalar(@TypeOf(z)))) {
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

                    var retval: @TypeOf(z) = undefined;
                    if (float.abs(x.im) > cast(Scalar(@TypeOf(z)), t, .{})) {
                        const exp_t: Scalar(@TypeOf(z)) = float.exp(cast(Scalar(@TypeOf(z)), t, .{}));
                        var ix: Scalar(@TypeOf(z)) = float.abs(x.im);
                        if (std.math.signbit(x.im))
                            cosix = -cosix;
                        ix -= cast(Scalar(@TypeOf(z)), t, .{});
                        sinix *= exp_t / 2;
                        cosix *= exp_t / 2;
                        if (ix > cast(Scalar(@TypeOf(z)), t, .{})) {
                            ix -= cast(Scalar(@TypeOf(z)), t, .{});
                            sinix *= exp_t;
                            cosix *= exp_t;
                        }

                        if (ix > cast(Scalar(@TypeOf(z)), t, .{})) {
                            // Overflow (original imaginary part of x > 3t).
                            retval.re = std.math.floatMax(Scalar(@TypeOf(z))) * sinix;
                            retval.im = std.math.floatMax(Scalar(@TypeOf(z))) * cosix;
                        } else {
                            const exp_val: Scalar(@TypeOf(z)) = float.exp(ix);
                            retval.re = exp_val * sinix;
                            retval.im = exp_val * cosix;
                        }
                    } else {
                        retval.re = float.cosh(x.im) * sinix;
                        retval.im = float.sinh(x.im) * cosix;
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
                    if (icls == classify.ZERO) {
                        // Imaginary part is 0.
                        return .{
                            .re = x.re - x.re,
                            .im = x.im,
                        };
                    } else {
                        return .{
                            .re = std.math.nan(Scalar(@TypeOf(z))),
                            .im = std.math.nan(Scalar(@TypeOf(z))),
                        };
                    }
                }
            } else if (icls == classify.INFINITE) {
                // Imaginary part is infinite.
                if (rcls == classify.ZERO) {
                    // Real part is 0.
                    return .{
                        .re = float.copysign(@as(Scalar(@TypeOf(z)), 0), @as(Scalar(@TypeOf(z)), if (negate) -1 else 1)),
                        .im = x.im,
                    };
                } else if (rcls > classify.ZERO) {
                    // Real part is finite.
                    var sinix: Scalar(@TypeOf(z)) = undefined;
                    var cosix: Scalar(@TypeOf(z)) = undefined;
                    if (x.re > std.math.floatMin(Scalar(@TypeOf(z)))) {
                        @branchHint(.likely);
                        const tmp = float.sincos(x.re);
                        sinix = tmp.sinx;
                        cosix = tmp.cosx;
                    } else {
                        sinix = x.re;
                        cosix = 1;
                    }

                    var retval: @TypeOf(z) = .{
                        .re = float.copysign(std.math.inf(Scalar(@TypeOf(z))), sinix),
                        .im = float.copysign(std.math.inf(Scalar(@TypeOf(z))), cosix),
                    };

                    if (negate)
                        retval.re = -retval.re;

                    if (std.math.signbit(x.im))
                        retval.im = -retval.im;

                    return retval;
                } else {
                    return .{
                        .re = x.re - x.re,
                        .im = std.math.inf(Scalar(@TypeOf(z))),
                    };
                }
            } else {
                return .{
                    .re = if (rcls == classify.ZERO) float.copysign(@as(Scalar(@TypeOf(z)), 0), @as(Scalar(@TypeOf(z)), if (negate) -1 else 1)) else std.math.nan(Scalar(@TypeOf(z))),
                    .im = std.math.nan(Scalar(@TypeOf(z))),
                };
            }
        },
        else => unreachable,
    }
}
