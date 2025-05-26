const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const classify = @import("../float/classify.zig");
const asinh = @import("asinh.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const scast = types.scast;

pub fn acos(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.acos: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    const zz: Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) = scast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z);

    const rcls: u32 = classify.classify(zz.re);
    const icls: u32 = classify.classify(zz.im);

    if (rcls <= classify.INFINITE or icls <= classify.INFINITE or (rcls == classify.ZERO and icls == classify.ZERO)) {
        const y = cfloat.asin(zz);

        var res: Scalar(@TypeOf(zz)) = float.pi_2(Scalar(@TypeOf(zz))) - y.re;

        if (res == 0)
            res = 0;

        return .{
            .re = res,
            .im = -y.im,
        };
    } else {
        var y: @TypeOf(zz) = .{
            .re = -zz.im,
            .im = zz.re,
        };

        y = asinh.kernel_asinh(y, 1);

        return .{
            .re = y.im,
            .im = y.re,
        };
    }
}
