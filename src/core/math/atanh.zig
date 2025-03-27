const std = @import("std");

// TODO: implement atanh for f128
pub fn atanh128(x: f128) f128 {
    return @floatCast(std.math.atanh(@as(f64, @floatCast(x))));
}
