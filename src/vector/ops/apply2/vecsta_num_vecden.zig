const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    if (y.inc == 1) {
        inline for (0..meta.Child(@TypeOf(o)).len) |i| {
            opInto(&o.data[i], x, y.data[i]);
        }
    } else {
        const iy: isize = if (y.inc < 0) (-numeric.cast(isize, meta.Child(@TypeOf(o)).len) + 1) * y.inc else 0;
        inline for (0..meta.Child(@TypeOf(o)).len) |i| {
            opInto(&o.data[i], x, y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]);
        }
    }
}
