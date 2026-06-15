const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const Y = @TypeOf(y);

    var i: usize = 0;
    if (o.inc == 1) {
        var iy: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                opInto(&o.data[i], x.data[i], numeric.zero(meta.Numeric(Y)));
            }

            opInto(&o.data[i], x.data[i], y.data[iy]);

            i += 1;
        }

        while (i < o.len) : (i += 1) {
            opInto(&o.data[i], x.data[i], numeric.zero(meta.Numeric(Y)));
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var iy: usize = 0;

        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                opInto(&o.data[numeric.cast(usize, io)], x.data[i], numeric.zero(meta.Numeric(Y)));

                io += o.inc;
            }

            opInto(&o.data[numeric.cast(usize, io)], x.data[i], y.data[iy]);

            i += 1;
            io += o.inc;
        }

        while (i < o.len) : (i += 1) {
            opInto(&o.data[numeric.cast(usize, io)], x.data[i], numeric.zero(meta.Numeric(Y)));

            io += o.inc;
        }
    }
}
