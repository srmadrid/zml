const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

pub fn check_aliasing(o: *const integer.Integer, x: anytype) bool {
    switch (comptime types.numericType(@TypeOf(x))) {
        .bool, .int, .float, .cfloat => return false,
        .integer => return o.limbs == x.limbs,
        else => unreachable,
    }
}
