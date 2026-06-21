const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X: type = @TypeOf(x);

    if (o.inc == 1 and y.inc == 1) {
        var ix: usize = 0;
        var i: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                opInto(&o.data[i], numeric.zero(meta.Numeric(X)), y.data[i]);
            }

            opInto(&o.data[i], x.data[ix], y.data[i]);

            i += 1;
        }

        while (i < o.len) : (i += 1) {
            opInto(&o.data[i], numeric.zero(meta.Numeric(X)), y.data[i]);
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var ix: usize = 0;
        var iy: isize = if (y.inc < 0) (-numeric.cast(isize, y.len) + 1) * y.inc else 0;
        var i: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                opInto(&o.data[numeric.cast(usize, io)], numeric.zero(meta.Numeric(X)), y.data[numeric.cast(usize, iy)]);

                io += o.inc;
                iy += y.inc;
            }

            opInto(&o.data[numeric.cast(usize, io)], x.data[ix], y.data[numeric.cast(usize, iy)]);

            i += 1;
            io += o.inc;
            iy += y.inc;
        }

        while (i < o.len) : (i += 1) {
            opInto(&o.data[numeric.cast(usize, io)], numeric.zero(meta.Numeric(X)), y.data[numeric.cast(usize, iy)]);

            io += o.inc;
            iy += y.inc;
        }
    }
}
