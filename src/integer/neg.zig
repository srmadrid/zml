const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn neg(allocator: ?std.mem.Allocator, x: Integer) !Integer {
    if (allocator) |a| {
        var result: Integer = try x.copy(a);
        result.positive = !result.positive;
        return result;
    }

    var result: Integer = x;
    result.positive = !result.positive;
    result.flags.owns_data = false;
    return result;
}
