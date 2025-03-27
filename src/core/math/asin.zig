const std = @import("std");

// TODO: implement asin for f128
pub fn asin128(x: f128) f128 {
    return @floatCast(std.math.asin(@as(f64, @floatCast(x))));
}
