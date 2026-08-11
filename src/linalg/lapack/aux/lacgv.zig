const numeric = @import("../../../numeric.zig");

pub fn lacgv(n: usize, x: anytype, incx: isize) !void {
    if (incx == 1) {
        for (0..n) |i| {
            numeric.conjInto(&x[i], x[i]);
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        for (0..n) |_| {
            numeric.conjInto(&x[numeric.cast(usize, ix)], x[numeric.cast(usize, ix)]);

            ix += incx;
        }
    }
}
