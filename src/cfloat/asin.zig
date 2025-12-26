const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

pub fn asin(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.asin: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    if (z.im == 0.0) {
        return if (float.abs(z.re) > 1.0)
            .{
                .re = 1.570796326794896619231321691639751442098585,
                .im = 0.0,
            }
        else
            .{
                .re = float.asin(z.re),
                .im = 0.0,
            };
    }

    var b: @TypeOf(z.re) = cfloat.abs(z);
    if (float.abs(b) < 0.125) {
        const z2: @TypeOf(z) = .{
            .re = (z.re - z.im) * (z.re + z.im),
            .im = (2.0 * z.re * z.im),
        };
        var cn: @TypeOf(z.re) = 1.0;
        var n: @TypeOf(z.re) = 1.0;
        var ca: @TypeOf(z) = z;
        var sum: @TypeOf(z) = z;
        while (true) {
            var ct: @TypeOf(z) = z2.mul(ca);
            ca = ct;

            cn *= n;
            n += 1.0;
            cn /= n;
            n += 1.0;
            b = cn / n;

            ct = ct.mulReal(b);
            sum = sum.add(ct);
            b = cfloat.abs(ct);

            if (b <= std.math.floatEps(@TypeOf(z.re)))
                break;
        }

        return sum;
    }

    const ct: @TypeOf(z) = z.mulImag(1.0);
    var zz: @TypeOf(z) = .{
        .re = 1.0 - (z.re - z.im) * (z.re + z.im),
        .im = -(2.0 * z.re * z.im),
    };
    const z2: @TypeOf(z) = cfloat.sqrt(zz);
    zz = ct.add(z2);
    zz = cfloat.log(zz);
    return zz.mulImag(-1.0);
}
