const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn ldexp(x: anytype, exp: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return cast(f16, std.math.ldexp(cast(f32, x), cast(i32, exp))),
        f32 => return std.math.ldexp(x, cast(i32, exp)),
        f64 => return std.math.ldexp(x, cast(i32, exp)),
        f80 => return cast(f80, std.math.ldexp(cast(f64, x), cast(i32, exp))),
        f128 => return cast(f128, std.math.ldexp(cast(f64, x), cast(i32, exp))),
        else => @compileError("x must be a float"),
    }
}
