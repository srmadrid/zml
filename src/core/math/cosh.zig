const std = @import("std");

// TODO: implement cosh for f128
pub fn cosh128(x: f128) f128 {
    return @floatCast(std.math.cosh(@as(f64, @floatCast(x))));
}
