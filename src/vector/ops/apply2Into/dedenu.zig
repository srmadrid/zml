const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    var i: usize = 0;
    if (o.inc == 1 and x.inc == 1) {
        while (i < o.len) : (i += 1) {
            opInto(&o.data[i], x.data[i], y);
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
        while (i < o.len) : (i += 1) {
            opInto(&o.data[numeric.cast(usize, io)], x.data[numeric.cast(usize, ix)], y);

            io += o.inc;
            ix += x.inc;
        }
    }
}
