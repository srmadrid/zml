const std = @import("std");

// TODO: implement acos for f128
pub fn acos128(x: f128) f128 {
    return @floatCast(std.math.acos(@as(f64, @floatCast(x))));
}
