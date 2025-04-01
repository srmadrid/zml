const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn floor(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @floor(x),
        f32 => return @floor(x),
        f64 => return @floor(x),
        f80 => return @floor(x),
        f128 => return @floor(x),
        else => @compileError("x must be a float"),
    }
}
