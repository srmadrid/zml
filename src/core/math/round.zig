const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn round(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @round(x),
        f32 => return @round(x),
        f64 => return @round(x),
        f80 => return @round(x),
        f128 => return @round(x),
        else => @compileError("x must be a float"),
    }
}
