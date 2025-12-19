const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const Scalar = types.Scalar;
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const Cfloat = cfloat.Cfloat;

pub fn sqrt(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.sqrt: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    switch (comptime types.numericType(@TypeOf(z))) {
        .int, .float => {
            if (z >= 0)
                return .{
                    .re = float.sqrt(z),
                    .im = 0,
                }
            else
                return .{
                    .re = 0,
                    .im = float.sqrt(-z),
                };
        },
        .cfloat => {
            if (z.re == 0.0 and z.im == 0.0)
                return .{
                    .re = 0,
                    .im = z.im,
                };

            if (std.math.isInf(z.im))
                return .{
                    .re = std.math.inf(@TypeOf(z.re)),
                    .im = z.im,
                };

            if (std.math.isNan(z.re))
                return .{
                    .re = std.math.nan(@TypeOf(z.re)),
                    .im = std.math.nan(@TypeOf(z.im)),
                };

            if (std.math.isInf(z.re)) {
                if (std.math.signbit(z.re))
                    return .{
                        .re = float.abs(z.im - z.im),
                        .im = float.copysign(z.re, z.im),
                    }
                else
                    return .{
                        .re = z.re,
                        .im = float.copysign(z.im - z.im, z.im),
                    };
            }

            switch (comptime @TypeOf(z.re)) {
                f16 => {
                    if (z.re >= 0) {
                        const t: f32 = float.sqrt((types.scast(f32, z.re) + float.hypot(z.re, z.im)) * 0.5);
                        return .{
                            .re = types.scast(f16, t),
                            .im = types.scast(f16, types.scast(f32, z.im) / (2.0 * t)),
                        };
                    } else {
                        const t: f32 = float.sqrt((types.scast(f32, -z.re) + float.hypot(z.re, z.im)) * 0.5);
                        return .{
                            .re = types.scast(f16, types.scast(f32, float.abs(z.im)) / (2.0 * t)),
                            .im = types.scast(f16, float.copysign(t, z.im)),
                        };
                    }
                },
                f32 => {
                    if (z.re >= 0) {
                        const t: f64 = float.sqrt((types.scast(f64, z.re) + float.hypot(types.scast(f64, z.re), types.scast(f64, z.im))) * 0.5);
                        return .{
                            .re = types.scast(f32, t),
                            .im = types.scast(f32, types.scast(f64, z.im) / (2.0 * t)),
                        };
                    } else {
                        const t: f64 = float.sqrt((types.scast(f64, -z.re) + float.hypot(types.scast(f64, z.re), types.scast(f64, z.im))) * 0.5);
                        return .{
                            .re = types.scast(f32, types.scast(f64, float.abs(z.im)) / (2.0 * t)),
                            .im = types.scast(f32, float.copysign(t, z.im)),
                        };
                    }
                },
                f64, f80, f128 => {
                    var a: @TypeOf(z.re) = z.re;
                    var b: @TypeOf(z.im) = z.im;

                    var scale: @TypeOf(z.re) = 1.0;

                    const abs_a = float.abs(a);
                    const abs_b = float.abs(b);

                    const thresh_max = std.math.floatMax(@TypeOf(z.re)) / 2.414213562373095048801688724209698; // max / 1 + sqrt(2)
                    const thresh_min = std.math.floatMin(@TypeOf(z.re)) * 2.414213562373095048801688724209698; // min * (1 + sqrt(2))
                    if (abs_a >= thresh_max or abs_b >= thresh_max) {
                        a *= 0.25;
                        b *= 0.25;
                        scale = 2.0; // sqrt(z/4) * 2 = sqrt(z)
                    } else if (abs_a <= thresh_min and abs_b <= thresh_min) {
                        a *= 4.0;
                        b *= 4.0;
                        scale = 0.5; // sqrt(z*4) * 0.5 = sqrt(z)

                    }

                    const t: @TypeOf(z.re) = float.sqrt((float.abs(a) + float.hypot(a, b)) * 0.5);

                    var result: @TypeOf(z) = undefined;
                    if (a >= 0) {
                        result.re = t * scale;
                        result.im = (b * scale) / (2.0 * t);
                    } else {
                        result.re = (float.abs(b) * scale) / (2.0 * t);
                        result.im = float.copysign(t * scale, b);
                    }

                    return result;
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}
