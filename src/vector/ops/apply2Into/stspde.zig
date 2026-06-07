const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    var i: usize = 0;
    if (y.inc == 1) {
        var ix: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                if (comptime opInto == numeric.addInto)
                    numeric.set(&o.data[i], y.data[i])
                else
                    numeric.set(&o.data[i], numeric.neg(y.data[i]));
            }

            opInto(&o.data[i], x.data[ix], y.data[i]);

            i += 1;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            if (comptime opInto == numeric.addInto)
                numeric.set(&o.data[i], y.data[i])
            else
                numeric.set(&o.data[i], numeric.neg(y.data[i]));
        }
    } else {
        var ix: usize = 0;
        var iy: isize = if (y.inc < 0) (-numeric.cast(isize, y.len) + 1) * y.inc else 0;

        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                if (comptime opInto == numeric.addInto)
                    numeric.set(&o.data[i], y.data[numeric.cast(usize, iy)])
                else
                    numeric.set(&o.data[i], numeric.neg(y.data[numeric.cast(usize, iy)]));

                iy += y.inc;
            }

            opInto(&o.data[i], x.data[ix], y.data[numeric.cast(usize, iy)]);

            i += 1;
            iy += y.inc;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            if (comptime opInto == numeric.addInto)
                numeric.set(&o.data[i], y.data[numeric.cast(usize, iy)])
            else
                numeric.set(&o.data[i], numeric.neg(y.data[numeric.cast(usize, iy)]));

            iy += y.inc;
        }
    }
}
