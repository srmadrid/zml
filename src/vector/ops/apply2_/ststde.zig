const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    if (y.inc == 1) {
        inline for (0..@TypeOf(x).len) |i| {
            op_(&o.data[i], x.data[i], y.data[i]);
        }
    } else {
        const iy: isize = if (y.inc < 0) (-numeric.cast(isize, @TypeOf(x).len) + 1) * y.inc else 0;
        inline for (0..@TypeOf(x).len) |i| {
            op_(&o.data[i], x.data[i], y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]);
        }
    }
}
