const std = @import("std");

const types = @import("../types.zig");

const check_aliasing = @import("check_aliasing.zig").check_aliasing;

pub fn check_aliasing_alloc(allocator: std.mem.Allocator, o: anytype, x: anytype) !@TypeOf(x) {
    switch (comptime types.numericType(@TypeOf(x))) {
        .bool, .int, .float, .cfloat => return x,
        .integer, .rational, .real, .complex => return if (check_aliasing(o, x))
            x.copy(allocator)
        else
            x,
        else => unreachable,
    }
}
