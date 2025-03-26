const std = @import("std");

// TODO: implement cosh for f128
pub fn sinh128(x: f128) f128 {
    return @floatCast(std.math.sinh(@as(f64, @floatCast(x))));
}
