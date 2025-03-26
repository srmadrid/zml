const std = @import("std");
const atan128 = @import("atan.zig").atan128;

// TODO: implement atan2 for f128
pub fn atan2_128(y: f128, x: f128) f128 {
    return @floatCast(std.math.atan2(@as(f64, @floatCast(y)), @as(f64, @floatCast(x))));
}
