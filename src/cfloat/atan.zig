const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

pub fn atan(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.atan: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    if (z.re == 0.0 and z.im > 1.0)
        return .{
            .re = std.math.inf(@TypeOf(z.re)),
            .im = std.math.inf(@TypeOf(z.im)),
        };

    const x2: @TypeOf(z.re) = z.re * z.re;
    var a: @TypeOf(z.re) = 1.0 - x2 - (z.im * z.im);
    if (a == 0.0)
        return .{
            .re = std.math.inf(@TypeOf(z.re)),
            .im = std.math.inf(@TypeOf(z.im)),
        };

    const w: @TypeOf(z.re) = float.atan2(2.0 * z.re, a) * 0.5;
    var t: @TypeOf(z.re) = w / float.pi(@TypeOf(z.re));
    if (t >= 0.0)
        t += 0.5
    else
        t -= 0.5;

    const i: i128 = types.scast(i128, t);
    t = types.scast(@TypeOf(z.re), i);
    t = switch (@TypeOf(z.re)) {
        f16, f32 => ((w -
            t * 3.140625) -
            t * 9.67502593994140625e-4) -
            t * 1.509957990978376432e-7,
        f64 => ((w -
            t * 3.14159265160560607910) -
            t * 1.98418714791870343106e-9) -
            t * 1.14423774522196636802e-17,
        f80, f128 => ((w -
            t * 3.14159265358979323829596852490908531763125) -
            t * 1.6667485837041756656403424829301998703007e-19) -
            t * 1.8830410776607851167459095484560349402753e-39,
        else => unreachable,
    };
    var r: @TypeOf(z) = undefined;
    r.re = t;

    t = z.im - 1.0;
    a = x2 + (t * t);
    if (a == 0.0)
        return .{
            .re = std.math.inf(@TypeOf(z.re)),
            .im = std.math.inf(@TypeOf(z.im)),
        };

    t = z.im + 1.0;
    a = (x2 + (t * t)) / a;
    r.im = 0.25 * float.log(a);
    return r;
}
