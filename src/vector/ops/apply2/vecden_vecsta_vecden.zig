const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    if (o.inc == 1 and y.inc == 1) {
        inline for (0..@TypeOf(x).len) |i| {
            opInto(&o.data[i], x.data[i], y.data[i]);
        }
    } else {
        const io: isize = if (o.inc < 0) (-numeric.cast(isize, @TypeOf(x).len) + 1) * o.inc else 0;
        const iy: isize = if (y.inc < 0) (-numeric.cast(isize, @TypeOf(x).len) + 1) * y.inc else 0;
        inline for (0..@TypeOf(x).len) |i| {
            opInto(&o.data[numeric.cast(usize, io + numeric.cast(isize, i) * o.inc)], x.data[i], y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]);
        }
    }
}
