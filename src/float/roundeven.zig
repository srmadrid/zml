const std = @import("std");
const float = @import("../float.zig");

// Round x to nearest integer value in floating-point format, rounding halfway
// cases to even.  If the input is non finite the result is unspecified.
pub inline fn roundeven_finite(x: f64) f64 {
    if (std.math.isInf(x))
        unreachable;

    var y: f64 = float.round(x);
    if (float.abs(x - y) == 0.5) {
        const u: u64 = @bitCast(y);
        const v: u64 = @bitCast(y - float.copysign(@as(f64, 1), x));
        if (@ctz(v) > @ctz(u))
            y = @bitCast(v);
    }
    return y;
}
