const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn asinh(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return cast(f16, std.math.asinh(cast(f32, x))),
        f32 => return std.math.asinh(x),
        f64 => return std.math.asinh(x),
        f80 => return cast(f80, std.math.asinh(cast(f64, x))),
        f128 => return cast(f128, std.math.asinh(cast(f64, x))),
        else => @compileError("x must be a float"),
    }
}
