const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn pow(left: anytype, right: anytype) @TypeOf(left, right) {
    switch (@TypeOf(right)) {
        f16 => return cast(f16, std.math.pow(f32, cast(f32, left), cast(f32, right))),
        f32 => return std.math.pow(f32, left, right),
        f64 => return std.math.pow(f32, left, right),
        f80 => return cast(f80, std.math.pow(f64, cast(f64, left), cast(f64, right))),
        f128 => return cast(f128, std.math.pow(f64, cast(f64, left), cast(f64, right))),
        else => @compileError("left must be a float"),
    }
}
