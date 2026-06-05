const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    if (o.inc == 1) {
        inline for (0..@TypeOf(x).len) |i| {
            op_(&o.data[i], x.data[i], y.data[i]);
        }
    } else {
        const io: isize = if (o.inc < 0) (-numeric.cast(isize, @TypeOf(x).len) + 1) * o.inc else 0;
        inline for (0..@TypeOf(x).len) |i| {
            op_(&o.data[numeric.cast(usize, io + numeric.cast(isize, i) * o.inc)], x.data[i], y.data[i]);
        }
    }
}
