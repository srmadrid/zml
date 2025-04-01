const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn asin(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return cast(f16, std.math.asin(cast(f32, x))),
        f32 => return std.math.asin(x),
        f64 => return std.math.asin(x),
        f80 => return cast(f80, std.math.asin(cast(f64, x))),
        f128 => return cast(f128, std.math.asin(cast(f64, x))),
        else => @compileError("x must be a float"),
    }
}
