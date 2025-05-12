const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const classify = @import("../float/classify.zig");
const asinh = @import("asinh.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn acos(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return acos(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return acos(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const rcls: u32 = classify.classify(z.re);
            const icls: u32 = classify.classify(z.im);

            if (rcls <= classify.INFINITE or icls <= classify.INFINITE or (rcls == classify.ZERO and icls == classify.ZERO)) {
                const y = cfloat.asin(z);

                var res: Scalar(@TypeOf(z)) = float.pi_2(Scalar(@TypeOf(z))) - y.re;

                if (res == 0)
                    res = 0;

                return .{
                    .re = res,
                    .im = -y.im,
                };
            } else {
                var y: @TypeOf(z) = .{
                    .re = -z.im,
                    .im = z.re,
                };

                y = asinh.kernel_asinh(y, 1);

                return .{
                    .re = y.im,
                    .im = y.re,
                };
            }
        },
        else => unreachable,
    }
}
