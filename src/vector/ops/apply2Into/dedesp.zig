const std = @import("std");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    const aliased = (comptime O == X) and std.meta.eql(o.*, x);

    if (aliased) {
        if (o.inc == 1) {
            var iy: usize = 0;
            while (iy < y.nnz) : (iy += 1) {
                opInto(&o.data[y.idx[iy]], x.data[y.idx[iy]], y.data[iy]);
            }
        } else {
            const io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
            var iy: usize = 0;
            while (iy < y.nnz) : (iy += 1) {
                const idx = numeric.cast(usize, io + numeric.cast(isize, y.idx[iy]) * o.inc);
                opInto(&o.data[idx], x.data[idx], y.data[iy]);
            }
        }
    } else {
        var i: usize = 0;
        if (o.inc == 1 and x.inc == 1) {
            var iy: usize = 0;
            while (iy < y.nnz) : (iy += 1) {
                while (i < y.idx[iy]) : (i += 1) {
                    numeric.set(&o.data[i], x.data[i]);
                }

                opInto(&o.data[i], x.data[i], y.data[iy]);

                i += 1;
            }

            while (i < o.len) : (i += 1) {
                numeric.set(&o.data[i], x.data[i]);
            }
        } else {
            var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
            var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
            var iy: usize = 0;

            while (iy < y.nnz) : (iy += 1) {
                while (i < y.idx[iy]) : (i += 1) {
                    numeric.set(&o.data[numeric.cast(usize, io)], x.data[numeric.cast(usize, ix)]);

                    io += o.inc;
                    ix += x.inc;
                }

                opInto(&o.data[numeric.cast(usize, io)], x.data[numeric.cast(usize, ix)], y.data[iy]);

                i += 1;
                io += o.inc;
                ix += x.inc;
            }

            while (i < o.len) : (i += 1) {
                numeric.set(&o.data[numeric.cast(usize, io)], x.data[numeric.cast(usize, ix)]);

                io += o.inc;
                ix += x.inc;
            }
        }
    }
}
