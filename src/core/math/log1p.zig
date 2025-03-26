const std = @import("std");

// TODO: implement log for f128
pub fn log1p128(x: f128) f128 {
    return @floatCast(std.math.log1p(@as(f64, @floatCast(x))));
}
