const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));

    var i: usize = 0;
    if (o.inc == 1) {
        var ix: usize = 0;
        var iy: usize = 0;

        while (i < o.len) : (i += 1) {
            if (ix < x.nnz and x.idx[ix] == i) {
                if (iy < y.nnz and y.idx[iy] == i) {
                    opInto(&o.data[i], x.data[ix], y.data[iy]);

                    ix += 1;
                    iy += 1;
                } else {
                    numeric.set(&o.data[i], x.data[ix]);

                    ix += 1;
                }
            } else {
                if (iy < y.nnz and y.idx[iy] == i) {
                    if (comptime opInto == numeric.addInto)
                        numeric.set(&o.data[i], y.data[iy])
                    else
                        numeric.set(&o.data[i], numeric.neg(y.data[iy]));

                    iy += 1;
                } else {
                    o.data[i] = numeric.zero(meta.Numeric(O));
                }
            }
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var ix: usize = 0;
        var iy: usize = 0;

        while (i < o.len) : (i += 1) {
            if (ix < x.nnz and x.idx[ix] == i) {
                if (iy < y.nnz and y.idx[iy] == i) {
                    opInto(&o.data[numeric.cast(usize, io)], x.data[ix], y.data[iy]);

                    ix += 1;
                    iy += 1;
                } else {
                    numeric.set(&o.data[numeric.cast(usize, io)], x.data[ix]);

                    ix += 1;
                }
            } else {
                if (iy < y.nnz and y.idx[iy] == i) {
                    if (comptime opInto == numeric.addInto)
                        numeric.set(&o.data[numeric.cast(usize, io)], y.data[iy])
                    else
                        numeric.set(&o.data[numeric.cast(usize, io)], numeric.neg(y.data[iy]));

                    iy += 1;
                } else {
                    o.data[numeric.cast(usize, io)] = numeric.zero(meta.Numeric(O));
                }
            }

            io += o.inc;
        }
    }
}
