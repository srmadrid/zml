const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

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
                    opInto(&o.data[i], x.data[ix], numeric.zero(meta.Numeric(Y)));

                    ix += 1;
                }
            } else {
                if (iy < y.nnz and y.idx[iy] == i) {
                    opInto(&o.data[i], numeric.zero(meta.Numeric(X)), y.data[iy]);

                    iy += 1;
                } else {
                    opInto(&o.data[i], numeric.zero(meta.Numeric(X)), numeric.zero(meta.Numeric(Y)));
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
                    opInto(&o.data[numeric.cast(usize, io)], x.data[ix], numeric.zero(meta.Numeric(Y)));

                    ix += 1;
                }
            } else {
                if (iy < y.nnz and y.idx[iy] == i) {
                    opInto(&o.data[numeric.cast(usize, io)], numeric.zero(meta.Numeric(X)), y.data[iy]);

                    iy += 1;
                } else {
                    opInto(&o.data[numeric.cast(usize, io)], numeric.zero(meta.Numeric(X)), numeric.zero(meta.Numeric(Y)));
                }
            }

            io += o.inc;
        }
    }
}
