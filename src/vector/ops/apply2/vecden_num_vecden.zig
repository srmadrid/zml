const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    if (o.inc == 1 and y.inc == 1) {
        var i: usize = 0;
        while (i < o.len) : (i += 1) {
            opInto(&o.data[i], x, y.data[i]);
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var iy: isize = if (y.inc < 0) (-numeric.cast(isize, y.len) + 1) * y.inc else 0;
        var i: usize = 0;
        while (i < o.len) : (i += 1) {
            opInto(&o.data[numeric.cast(usize, io)], x, y.data[numeric.cast(usize, iy)]);

            io += o.inc;
            iy += y.inc;
        }
    }
}
