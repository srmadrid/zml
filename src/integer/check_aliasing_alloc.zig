const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");

/// Checks if `x` is aliased with `o`, and if so, returns a copy of `x`
/// allocated with `allocator`. Otherwise, returns `x`.
pub fn check_aliasing_alloc(allocator: std.mem.Allocator, o: *const integer.Integer, x: anytype) !@TypeOf(x) {
    switch (comptime types.numericType(@TypeOf(x))) {
        .bool, .int, .float, .dyadic, .cfloat => return x,
        .integer => return if (@import("check_aliasing.zig").check_aliasing(o, x))
            x.copy(allocator)
        else
            x,
        else => unreachable,
    }
}
