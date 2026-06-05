const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    if (x.inc == 1 and y.inc == 1) {
        inline for (0..meta.Child(@TypeOf(o)).len) |i| {
            op_(&o.data[i], x.data[i], y.data[i]);
        }
    } else {
        const ix: isize = if (x.inc < 0) (-numeric.cast(isize, meta.Child(@TypeOf(o)).len) + 1) * x.inc else 0;
        const iy: isize = if (y.inc < 0) (-numeric.cast(isize, meta.Child(@TypeOf(o)).len) + 1) * y.inc else 0;
        inline for (0..meta.Child(@TypeOf(o)).len) |i| {
            op_(&o.data[i], x.data[numeric.cast(usize, ix + numeric.cast(isize, i) * x.inc)], y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]);
        }
    }
}
