const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    var i: usize = 0;
    if (x.inc == 1) {
        var iy: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                numeric.set(&o.data[i], x.data[i]);
            }

            op_(&o.data[i], x.data[i], y.data[iy]);

            i += 1;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            numeric.set(&o.data[i], x.data[i]);
        }
    } else {
        var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
        var iy: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                numeric.set(&o.data[i], x.data[numeric.cast(usize, ix)]);

                ix += x.inc;
            }

            op_(&o.data[i], x.data[numeric.cast(usize, ix)], y.data[iy]);

            i += 1;
            ix += x.inc;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            numeric.set(&o.data[i], x.data[numeric.cast(usize, ix)]);

            ix += x.inc;
        }
    }
}
