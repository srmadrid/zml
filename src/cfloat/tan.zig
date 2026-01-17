const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

pub fn tan(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or !types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.tan: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    var d: @TypeOf(z.re) = float.cos(2.0 * z.re) + float.cosh(2.0 * z.im);

    if (float.abs(d) < 0.25) {
        var x: @TypeOf(z.re) = float.abs(2.0 * z.re);
        var y: @TypeOf(z.im) = float.abs(2.0 * z.im);

        var t: @TypeOf(z.re) = x / float.pi(@TypeOf(z.re));
        if (t >= 0.0)
            t += 0.5
        else
            t -= 0.5;

        const i: i128 = types.scast(i128, t);
        t = types.scast(@TypeOf(z.re), i);
        t = switch (@TypeOf(z.re)) {
            f16, f32 => ((x -
                t * 3.140625) -
                t * 9.67502593994140625e-4) -
                t * 1.509957990978376432e-7,
            f64 => ((x -
                t * 3.14159265160560607910) -
                t * 1.98418714791870343106e-9) -
                t * 1.14423774522196636802e-17,
            f80, f128 => ((x -
                t * 3.14159265358979323829596852490908531763125) -
                t * 1.6667485837041756656403424829301998703007e-19) -
                t * 1.8830410776607851167459095484560349402753e-39,
            else => unreachable,
        };
        x = t;

        x = x * x;
        y = y * y;
        var x2: @TypeOf(z.re) = 1.0;
        var y2: @TypeOf(z.im) = 1.0;
        var f: @TypeOf(z.re) = 1.0;
        var rn: @TypeOf(z.re) = 0.0;
        d = 0.0;
        while (true) {
            rn += 1.0;
            f *= rn;
            rn += 1.0;
            f *= rn;
            x2 *= x;
            y2 *= y;
            t = y2 + x2;
            t /= f;
            d += t;

            rn += 1.0;
            f *= rn;
            rn += 1.0;
            f *= rn;
            x2 *= x;
            y2 *= y;
            t = y2 - x2;
            t /= f;
            d += t;

            if (float.abs(t / d) <= std.math.floatEps(@TypeOf(z.re)))
                break;
        }
    }

    if (d == 0.0)
        return .{
            .re = std.math.inf(@TypeOf(z.re)),
            .im = std.math.inf(@TypeOf(z.im)),
        };

    return .{
        .re = float.sin(2.0 * z.re) / d,
        .im = float.sinh(2.0 * z.im) / d,
    };
}
