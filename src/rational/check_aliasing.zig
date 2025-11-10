const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

pub fn check_aliasing(o: *const rational.Rational, x: anytype) bool {
    switch (comptime types.numericType(@TypeOf(x))) {
        .bool, .int, .float, .cfloat => return false,
        .integer => return o.num.limbs == x.limbs or o.den.limbs == x.limbs,
        .rational => return o.num.limbs == x.num.limbs or o.num.limbs == x.den.limbs or
            o.den.limbs == x.num.limbs or o.den.limbs == x.den.limbs,
        else => unreachable,
    }
}
