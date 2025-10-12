const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn abs(allocator: ?std.mem.Allocator, x: Integer) !Integer {
    if (allocator) |a| {
        var result: Integer = try x.copy(a);
        result.positive = true;
        return result;
    }

    var result: Integer = x;
    result.positive = true;
    result.flags.owns_data = false;
    return result;
}
