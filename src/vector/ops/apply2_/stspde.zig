const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    var i: usize = 0;
    if (y.inc == 1) {
        var ix: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                if (comptime op_ == numeric.add_)
                    numeric.set(&o.data[i], y.data[i])
                else
                    numeric.set(&o.data[i], numeric.neg(y.data[i]));
            }

            op_(&o.data[i], x.data[ix], y.data[i]);

            i += 1;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            if (comptime op_ == numeric.add_)
                numeric.set(&o.data[i], y.data[i])
            else
                numeric.set(&o.data[i], numeric.neg(y.data[i]));
        }
    } else {
        var ix: usize = 0;
        var iy: isize = if (y.inc < 0) (-numeric.cast(isize, y.len) + 1) * y.inc else 0;

        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                if (comptime op_ == numeric.add_)
                    numeric.set(&o.data[i], y.data[numeric.cast(usize, iy)])
                else
                    numeric.set(&o.data[i], numeric.neg(y.data[numeric.cast(usize, iy)]));

                iy += y.inc;
            }

            op_(&o.data[i], x.data[ix], y.data[numeric.cast(usize, iy)]);

            i += 1;
            iy += y.inc;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            if (comptime op_ == numeric.add_)
                numeric.set(&o.data[i], y.data[numeric.cast(usize, iy)])
            else
                numeric.set(&o.data[i], numeric.neg(y.data[numeric.cast(usize, iy)]));

            iy += y.inc;
        }
    }
}
