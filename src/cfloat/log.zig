const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("../float/classify.zig");
const x2y2m1 = @import("../float/x2y2m1.zig").x2y2m1;
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn log(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.log: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    const rcls: u32 = classify.classify(zz.re);
    const icls: u32 = classify.classify(zz.im);

    if (rcls == classify.ZERO and icls == classify.ZERO) {
        @branchHint(.unlikely);
        // Real and imaginary part are 0.
        return .{
            .re = -1 / float.abs(zz.re),
            .im = float.copysign(@as(Scalar(@TypeOf(zz)), if (std.math.signbit(zz.re)) float.pi(Scalar(@TypeOf(zz))) else 0), zz.im),
        };
    } else if (rcls != classify.NAN and icls != classify.NAN) {
        @branchHint(.likely);
        // Neither real nor imaginary part is NaN.
        var absr: Scalar(@TypeOf(zz)) = float.abs(zz.re);
        var absi: Scalar(@TypeOf(zz)) = float.abs(zz.im);
        var scale: i32 = 0;

        if (absr < absi) {
            const t: Scalar(@TypeOf(zz)) = absr;
            absr = absi;
            absi = t;
        }

        if (absr > std.math.floatMax(Scalar(@TypeOf(zz))) / 2) {
            scale = -1;
            absr = float.scalbn(absr, scale);
            absi = if (absi >= std.math.floatMin(Scalar(@TypeOf(zz))) * 2) float.scalbn(absi, scale) else 0;
        } else if (absr < std.math.floatMin(Scalar(@TypeOf(zz))) and absi < std.math.floatMin(Scalar(@TypeOf(zz)))) {
            scale = std.math.floatMantissaBits(Scalar(@TypeOf(zz)));
            absr = float.scalbn(absr, scale);
            absi = float.scalbn(absi, scale);
        }

        var result: @TypeOf(zz) = undefined;
        if (absr == 1 and scale == 0) {
            result.re = float.log1p(absi * absi) / 2;

            if (result.re < std.math.floatMin(Scalar(@TypeOf(zz)))) {
                const v: Scalar(@TypeOf(zz)) = result.re * result.re;
                std.mem.doNotOptimizeAway(v);
            }
        } else if (absr > 1 and absr < 2 and absi < 1 and scale == 0) {
            var d2m1: Scalar(@TypeOf(zz)) = (absr - 1) * (absr + 1);
            if (absi >= std.math.floatEps(Scalar(@TypeOf(zz))))
                d2m1 += absi * absi;
            result.re = float.log1p(d2m1) / 2;
        } else if (absr < 1 and absr >= 0.5 and absi < std.math.floatEps(Scalar(@TypeOf(zz))) / 2 and scale == 0) {
            const d2m1: Scalar(@TypeOf(zz)) = (absr - 1) * (absr + 1);
            result.re = float.log1p(d2m1) / 2;
        } else if (absr < 1 and absr >= 0.5 and scale == 0 and absr * absr + absi * absi >= 0.5) {
            const d2m1: Scalar(@TypeOf(zz)) = x2y2m1(absr, absi);
            result.re = float.log1p(d2m1) / 2;
        } else {
            const d: Scalar(@TypeOf(zz)) = float.hypot(absr, absi);
            result.re = float.log(d) - scast(Scalar(@TypeOf(zz)), scale) * float.ln2(Scalar(@TypeOf(zz)));
        }

        result.im = float.atan2(zz.im, zz.re);

        return result;
    } else {
        return .{
            .re = if (rcls == classify.INFINITE or icls == classify.INFINITE) std.math.inf(Scalar(@TypeOf(zz))) else std.math.nan(Scalar(@TypeOf(zz))),
            .im = std.math.nan(Scalar(@TypeOf(zz))),
        };
    }
}
