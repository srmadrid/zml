const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const x2y2m1 = @import("../float/x2y2m1.zig").x2y2m1;
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn log10(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return log10(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return log10(cast(Cfloat(Scalar(@TypeOf(z))), z, .{}));
        },
        .cfloat => {
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            if (rcls == classify.ZERO and icls == classify.ZERO) {
                @branchHint(.unlikely);
                // Real and imaginary part are 0.
                return .{
                    .re = -1 / float.abs(z.re),
                    .im = float.copysign(@as(Scalar(@TypeOf(z)), if (std.math.signbit(z.re)) float.log10e(Scalar(@TypeOf(z))) else 0), z.im),
                };
            } else if (rcls != classify.NAN and icls != classify.NAN) {
                @branchHint(.likely);
                // Neither real nor imaginary part is NaN.
                var absr: Scalar(@TypeOf(z)) = float.abs(z.re);
                var absi: Scalar(@TypeOf(z)) = float.abs(z.im);
                var scale: i32 = 0;

                if (absr < absi) {
                    const t: Scalar(@TypeOf(z)) = absr;
                    absr = absi;
                    absi = t;
                }

                if (absr > std.math.floatMax(Scalar(@TypeOf(z))) / 2) {
                    scale = -1;
                    absr = float.scalbn(absr, scale);
                    absi = if (absi >= std.math.floatMin(Scalar(@TypeOf(z))) * 2) float.scalbn(absi, scale) else 0;
                } else if (absr < std.math.floatMin(Scalar(@TypeOf(z))) and absi < std.math.floatMin(Scalar(@TypeOf(z)))) {
                    scale = std.math.floatMantissaBits(Scalar(@TypeOf(z)));
                    absr = float.scalbn(absr, scale);
                    absi = float.scalbn(absi, scale);
                }

                var result: Cfloat(Scalar(@TypeOf(z))) = undefined;
                if (absr == 1 and scale == 0) {
                    result.re = (float.log1p(absi * absi) * (@as(Scalar(@TypeOf(z)), float.log10e(Scalar(@TypeOf(z)))) / 2));

                    if (result.re < std.math.floatMin(Scalar(@TypeOf(z)))) {
                        const v: Scalar(@TypeOf(z)) = result.re * result.re;
                        std.mem.doNotOptimizeAway(v);
                    }
                } else if (absr > 1 and absr < 2 and absi < 1 and scale == 0) {
                    var d2m1: Scalar(@TypeOf(z)) = (absr - 1) * (absr + 1);
                    if (absi >= std.math.floatEps(Scalar(@TypeOf(z))))
                        d2m1 += absi * absi;
                    result.re = float.log1p(d2m1) * (@as(Scalar(@TypeOf(z)), float.log10e(Scalar(@TypeOf(z)))) / 2);
                } else if (absr < 1 and absr >= 0.5 and absi < std.math.floatEps(Scalar(@TypeOf(z))) / 2 and scale == 0) {
                    const d2m1: Scalar(@TypeOf(z)) = (absr - 1) * (absr + 1);
                    result.re = float.log1p(d2m1) * (@as(Scalar(@TypeOf(z)), float.log10e(Scalar(@TypeOf(z)))) / 2);
                } else if (absr < 1 and absr >= 0.5 and scale == 0 and absr * absr + absi * absi >= 0.5) {
                    const d2m1: Scalar(@TypeOf(z)) = x2y2m1(absr, absi);
                    result.re = float.log1p(d2m1) * (@as(Scalar(@TypeOf(z)), float.log10e(Scalar(@TypeOf(z)))) / 2);
                } else {
                    const log10_2 = 0.3010299956639811952137388947244930267682;
                    const d: Scalar(@TypeOf(z)) = float.hypot(absr, absi);
                    result.re = float.log10(d) - cast(Scalar(@TypeOf(z)), scale, .{}) * log10_2;
                }

                result.im = float.log10e(Scalar(@TypeOf(z))) * float.atan2(z.im, z.re);

                return result;
            } else {
                return .{
                    .re = if (rcls == classify.INFINITE or icls == classify.INFINITE) std.math.inf(Scalar(@TypeOf(z))) else std.math.nan(Scalar(@TypeOf(z))),
                    .im = std.math.nan(Scalar(@TypeOf(z))),
                };
            }
        },
        else => unreachable,
    }
}
