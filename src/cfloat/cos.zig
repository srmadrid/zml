const types = @import("../types.zig");
const float = @import("../float.zig");

pub fn cos(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.cos: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    var s: @TypeOf(z.re) = undefined;
    var c: @TypeOf(z.re) = undefined;
    if (float.abs(z.re) <= 0.5) {
        s = float.sinh(z.im);
        c = float.cosh(z.im);
    } else {
        var e: @TypeOf(z.re) = float.exp(z.im);
        const e_inv: @TypeOf(z.re) = 0.5 / e;
        e = 0.5 * e;
        s = e - e_inv;
        c = e + e_inv;
    }

    return .{
        .re = float.cos(z.re) * c,
        .im = float.sin(z.re) * s,
    };
}
