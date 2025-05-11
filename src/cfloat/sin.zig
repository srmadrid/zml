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

            const zz: @TypeOf(z) = .{ .re = float.abs(z.re), .im = z.im };

            if (icls >= classify.ZERO) {
                @branchHint(.likely);
                // Imaginary part is finite.
                if (rcls >= classify.ZERO) {
                    @branchHint(.likely);
                    // Real part is finite.
                    const t: i32 = cast(i32, (std.math.floatExponentMax(Scalar(@TypeOf(z))) - 1) * float.ln2(Scalar(@TypeOf(z))), .{});

                    var sinizz: Scalar(@TypeOf(z)) = undefined;
                    var cosizz: Scalar(@TypeOf(z)) = undefined;
                    if (zz.re > std.math.floatMin(Scalar(@TypeOf(z)))) {
                        @branchHint(.likely);
                        const tmp = float.sincos(zz.re);
                        sinizz = tmp.sinx;
                        cosizz = tmp.cosx;
                    } else {
                        sinizz = zz.re;
                        cosizz = 1;
                    }

                    if (negate)
                        sinizz = -sinizz;

                    var retval: @TypeOf(z) = undefined;
                    if (float.abs(zz.im) > t) {
                        const ezzp_t: Scalar(@TypeOf(z)) = float.exp(cast(Scalar(@TypeOf(z)), t, .{}));
                        var izz: Scalar(@TypeOf(z)) = float.abs(zz.im);
                        if (std.math.signbit(zz.im))
                            cosizz = -cosizz;
                        izz -= t;
                        sinizz *= ezzp_t / 2;
                        cosizz *= ezzp_t / 2;
                        if (izz > t) {
                            izz -= t;
                            sinizz *= ezzp_t;
                            cosizz *= ezzp_t;
                        }

                        if (izz > t) {
                            // Overflow (original imaginary part of zz > 3t).
                            retval.re = std.math.floatMax(Scalar(@TypeOf(z))) * sinizz;
                            retval.im = std.math.floatMax(Scalar(@TypeOf(z))) * cosizz;
                        } else {
                            const ezzp_val: Scalar(@TypeOf(z)) = float.exp(izz);
                            retval.re = ezzp_val * sinizz;
                            retval.im = ezzp_val * cosizz;
                        }
                    } else {
                        retval.re = float.cosh(zz.im) * sinizz;
                        retval.im = float.sinh(zz.im) * cosizz;
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
                            .re = zz.re - zz.re,
                            .im = zz.im,
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
                        .im = zz.im,
                    };
                } else if (rcls > classify.ZERO) {
                    // Real part is finite.
                    var sinizz: Scalar(@TypeOf(z)) = undefined;
                    var cosizz: Scalar(@TypeOf(z)) = undefined;
                    if (zz.re > std.math.floatMin(Scalar(@TypeOf(z)))) {
                        @branchHint(.likely);
                        const tmp = float.sincos(zz.re);
                        sinizz = tmp.sinx;
                        cosizz = tmp.cosx;
                    } else {
                        sinizz = zz.re;
                        cosizz = 1;
                    }

                    var retval: @TypeOf(z) = .{
                        .re = float.copysign(std.math.inf(Scalar(@TypeOf(z))), sinizz),
                        .im = float.copysign(std.math.inf(Scalar(@TypeOf(z))), cosizz),
                    };

                    if (negate)
                        retval.re = -retval.re;

                    if (std.math.signbit(zz.im))
                        retval.im = -retval.im;

                    return retval;
                } else {
                    return .{
                        .re = zz.re - zz.re,
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
