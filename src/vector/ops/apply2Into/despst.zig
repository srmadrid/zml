const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    var i: usize = 0;
    if (o.inc == 1) {
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

        while (i < o.len) : (i += 1) {
            if (comptime opInto == numeric.addInto)
                numeric.set(&o.data[i], y.data[i])
            else
                numeric.set(&o.data[i], numeric.neg(y.data[i]));
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var ix: usize = 0;

        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                if (comptime opInto == numeric.addInto)
                    numeric.set(&o.data[numeric.cast(usize, io)], y.data[i])
                else
                    numeric.set(&o.data[numeric.cast(usize, io)], numeric.neg(y.data[i]));

                io += o.inc;
            }

            opInto(&o.data[numeric.cast(usize, io)], x.data[ix], y.data[i]);

            i += 1;
            io += o.inc;
        }

        while (i < o.len) : (i += 1) {
            if (comptime opInto == numeric.addInto)
                numeric.set(&o.data[numeric.cast(usize, io)], y.data[i])
            else
                numeric.set(&o.data[numeric.cast(usize, io)], numeric.neg(y.data[i]));

            io += o.inc;
        }
    }
}
