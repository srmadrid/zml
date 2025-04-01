const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn trunc(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @trunc(x),
        f32 => return @trunc(x),
        f64 => return @trunc(x),
        f80 => return @trunc(x),
        f128 => return @trunc(x),
        else => @compileError("x must be a float"),
    }
}
