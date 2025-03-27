const std = @import("std");

// TODO: implement acosh for f128
pub fn acosh128(x: f128) f128 {
    return @floatCast(std.math.acosh(@as(f64, @floatCast(x))));
}
