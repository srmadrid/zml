const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn ceil(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @ceil(x),
        f32 => return @ceil(x),
        f64 => return @ceil(x),
        f80 => return @ceil(x),
        f128 => return @ceil(x),
        else => @compileError("x must be a float"),
    }
}
