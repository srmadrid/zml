const std = @import("std");

// TODO: implement pow for f128
pub fn pow128(x: f128, y: f128) f128 {
    return @floatCast(std.math.pow(f64, @as(f64, @floatCast(x)), @as(f64, @floatCast(y))));
}
