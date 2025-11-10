const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

pub fn check_aliasing_alloc(allocator: std.mem.Allocator, o: *const integer.Integer, x: anytype) !@TypeOf(x) {
    var tx: @TypeOf(x) = undefined;
    switch (comptime types.numericType(@TypeOf(x))) {
        .bool, .int, .float, .cfloat => {
            tx = x;
        },
        .integer => {
            if (@import("check_aliasing.zig").check_aliasing(o, x)) {
                tx = try x.copy(allocator);
            } else {
                tx = x;
                tx.flags.owns_data = false;
            }
        },
        else => unreachable,
    }

    return tx;
}
