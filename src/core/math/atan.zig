const std = @import("std");

// TODO: implement atan for f128
pub fn atan128(x: f128) f128 {
    return @floatCast(std.math.atan(@as(f64, @floatCast(x))));
}
