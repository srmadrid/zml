const std = @import("std");

// TODO: implement tanh for f128
pub fn tanh128(x: f128) f128 {
    return @floatCast(std.math.tanh(@as(f64, @floatCast(x))));
}
