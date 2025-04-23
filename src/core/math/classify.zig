const std = @import("std");

pub const NAN: u32 = 0;
pub const INFINITE: u32 = 1;
pub const ZERO: u32 = 2;
pub const SUBNORMAL: u32 = 3;
pub const NORMAL: u32 = 4;

pub inline fn classify(x: anytype) u32 {
    if (std.math.isNan(x)) {
        return NAN;
    }

    if (std.math.isInf(x)) {
        return INFINITE;
    }

    if (std.math.isPositiveZero(x) or std.math.isNegativeZero(x)) {
        return ZERO;
    }

    return NORMAL;
}
