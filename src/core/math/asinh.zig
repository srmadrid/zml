const std = @import("std");

// TODO: implement acosh for f128
pub fn asinh128(x: f128) f128 {
    return @floatCast(std.math.asinh(@as(f64, @floatCast(x))));
}
