const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const integer_check_aliasing = @import("../integer/check_aliasing.zig").check_aliasing;

pub fn check_aliasing(o: *const rational.Rational, x: anytype) bool {
    switch (comptime types.numericType(@TypeOf(x))) {
        .bool, .int, .float, .cfloat => return false,
        .integer => return integer_check_aliasing(&o.num, x) or integer_check_aliasing(&o.den, x),
        .rational => return integer_check_aliasing(&o.num, x.num) or integer_check_aliasing(&o.num, x.den) or
            integer_check_aliasing(&o.den, x.num) or integer_check_aliasing(&o.den, x.den),
        else => unreachable,
    }
}
