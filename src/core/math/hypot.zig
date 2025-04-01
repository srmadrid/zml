const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn hypot(x: anytype, y: anytype) @TypeOf(x, y) {
    switch (@TypeOf(y)) {
        f16 => return cast(f16, std.math.hypot(cast(f32, x), cast(f32, y))),
        f32 => return std.math.hypot(x, y),
        f64 => return std.math.hypot(x, y),
        f80 => return cast(f80, std.math.hypot(cast(f64, x), cast(f64, y))),
        f128 => return cast(f128, std.math.hypot(cast(f64, x), cast(f64, y))),
        else => @compileError("x must be a float"),
    }
}
