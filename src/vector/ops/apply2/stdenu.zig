const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    if (x.inc == 1) {
        inline for (0..meta.Child(@TypeOf(o)).len) |i| {
            opInto(&o.data[i], x.data[i], y);
        }
    } else {
        const ix: isize = if (x.inc < 0) (-numeric.cast(isize, meta.Child(@TypeOf(o)).len) + 1) * x.inc else 0;
        inline for (0..meta.Child(@TypeOf(o)).len) |i| {
            opInto(&o.data[i], x.data[numeric.cast(usize, ix + numeric.cast(isize, i) * x.inc)], y);
        }
    }
}
