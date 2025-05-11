const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;

pub fn sqrt(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            if (z >= 0) {
                return .{
                    .re = float.sqrt(z),
                    .im = 0,
                };
            } else {
                return .{
                    .re = 0,
                    .im = float.sqrt(-z),
                };
            }
        },
        .float => {
            if (z >= 0) {
                return .{
                    .re = float.sqrt(z),
                    .im = 0,
                };
            } else {
                return .{
                    .re = 0,
                    .im = float.sqrt(-z),
                };
            }
        },
        .cfloat => {
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            if (rcls <= classify.INFINITE or icls <= classify.INFINITE) {
                @branchHint(.unlikely);
                if (icls == classify.INFINITE) {
                    return .{
                        .re = std.math.inf(Scalar(@TypeOf(z))),
                        .im = z.im,
                    };
                } else if (rcls == classify.INFINITE) {
                    if (z.re < 0) {
                        return .{
                            .re = if (icls == classify.NAN) std.math.nan(Scalar(@TypeOf(z))) else 0,
                            .im = float.copysign(std.math.inf(Scalar(@TypeOf(z))), z.im),
                        };
                    } else {
                        return .{
                            .re = z.re,
                            .im = if (icls == classify.NAN) std.math.nan(Scalar(@TypeOf(z))) else float.copysign(@as(Scalar(@TypeOf(z)), 0), z.im),
                        };
                    }
                } else {
                    return .{
                        .re = std.math.nan(Scalar(@TypeOf(z))),
                        .im = std.math.nan(Scalar(@TypeOf(z))),
                    };
                }
            } else {
                if (icls == classify.ZERO) {
                    @branchHint(.unlikely);
                    if (z.re < 0) {
                        return .{
                            .re = 0,
                            .im = float.copysign(float.sqrt(-z.re), z.im),
                        };
                    } else {
                        return .{
                            .re = float.abs(float.sqrt(z.re)),
                            .im = float.copysign(@as(Scalar(@TypeOf(z)), 0), z.im),
                        };
                    }
                } else if (rcls == classify.ZERO) {
                    @branchHint(.unlikely);
                    var r: Scalar(@TypeOf(z)) = undefined;
                    if (float.abs(z.im) >= 2 * std.math.floatMin(Scalar(@TypeOf(z)))) {
                        r = float.sqrt(0.5 * float.abs(z.im));
                    } else {
                        r = 0.5 * float.sqrt(2 * float.abs(z.im));
                    }

                    return .{
                        .re = r,
                        .im = float.copysign(r, z.im),
                    };
                } else {
                    var scale: i32 = 0;
                    var x = z;

                    if (float.abs(x.re) > std.math.floatMax(Scalar(@TypeOf(z))) / 4) {
                        scale = 1;
                        x.re = float.scalbn(x.re, -2 * scale);
                        x.im = float.scalbn(x.im, -2 * scale);
                    } else if (float.abs(x.im) > std.math.floatMax(Scalar(@TypeOf(z))) / 4) {
                        scale = 1;
                        if (float.abs(x.re) >= 4 * std.math.floatMin(Scalar(@TypeOf(z)))) {
                            x.re = float.scalbn(x.re, -2 * scale);
                        } else {
                            x.re = 0;
                        }
                        x.im = float.scalbn(x.im, -2 * scale);
                    } else if (float.abs(x.re) < 2 * std.math.floatMin(Scalar(@TypeOf(z))) and float.abs(x.im) < 2 * std.math.floatMin(Scalar(@TypeOf(z)))) {
                        scale = -((std.math.floatMantissaBits(Scalar(@TypeOf(z))) + 1) / 2);
                        x.re = float.scalbn(x.re, -2 * scale);
                        x.im = float.scalbn(x.im, -2 * scale);
                    }

                    const d: Scalar(@TypeOf(z)) = float.hypot(x.re, x.im);
                    var r: Scalar(@TypeOf(z)) = undefined;
                    var s: Scalar(@TypeOf(z)) = undefined;
                    // Use the identity   2  Re res  Im res = Im x to avoid cancellation error in  d +/- Re x.
                    if (x.re > 0) {
                        r = float.sqrt(0.5 * (d + x.re));
                        if (scale == 1 and float.abs(x.im) < 1) {
                            // Avoid possible intermediate underflow.
                            s = x.im / r;
                            r = float.scalbn(r, scale);
                            scale = 0;
                        } else {
                            s = 0.5 * (x.im / r);
                        }
                    } else {
                        s = float.sqrt(0.5 * (d - x.re));
                        if (scale == 1 and float.abs(x.im) < 1) {
                            // Avoid possible intermediate underflow.
                            r = float.abs(x.im / s);
                            s = float.scalbn(s, scale);
                            scale = 0;
                        } else {
                            r = float.abs(0.5 * (x.im / s));
                        }
                    }

                    if (scale != 0) {
                        r = float.scalbn(r, scale);
                        s = float.scalbn(s, scale);
                    }

                    if (float.abs(r) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                        const vr: Scalar(@TypeOf(z)) = r * r;
                        std.mem.doNotOptimizeAway(vr);
                    }

                    if (float.abs(s) < std.math.floatMin(Scalar(@TypeOf(z)))) {
                        const vs: Scalar(@TypeOf(z)) = s * s;
                        std.mem.doNotOptimizeAway(vs);
                    }

                    return .{
                        .re = r,
                        .im = float.copysign(s, x.im),
                    };
                }
            }
        },
        else => unreachable,
    }
}
