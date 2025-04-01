const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn gamma(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return cast(f16, std.math.gamma(f32, cast(f32, x))),
        f32 => return std.math.gamma(f32, x),
        f64 => return std.math.gamma(f64, x),
        f80 => return cast(f80, std.math.gamma(f64, cast(f64, x))),
        f128 => return cast(f128, std.math.gamma(f64, cast(f64, x))),
        else => @compileError("x must be a float"),
    }
}
