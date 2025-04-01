const cast = @import("../types.zig").cast;
const std = @import("std");

pub fn remainder(left: anytype, right: anytype) @TypeOf(left, right) {
    switch (@TypeOf(right)) {
        f16 => return @rem(left, right),
        f32 => return @rem(left, right),
        f64 => return @rem(left, right),
        f80 => return @rem(left, right),
        f128 => return @rem(left, right),
        else => @compileError("left must be a float"),
    }
}
