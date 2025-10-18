const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

pub fn neg(allocator: ?std.mem.Allocator, x: Rational) !Rational {
    if (allocator) |a| {
        var result: Rational = try x.copy(a);
        result.num.positive = !result.num.positive;
        return result;
    }

    var result: Rational = x;
    result.num.positive = !result.num.positive;
    result.flags.owns_data = false;
    return result;
}
