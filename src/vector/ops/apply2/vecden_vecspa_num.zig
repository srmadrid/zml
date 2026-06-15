const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X = @TypeOf(x);

    var i: usize = 0;
    if (o.inc == 1) {
        var ix: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                opInto(&o.data[i], numeric.zero(meta.Numeric(X)), y);
            }

            opInto(&o.data[i], x.data[ix], y);

            i += 1;
        }

        while (i < o.len) : (i += 1) {
            opInto(&o.data[i], numeric.zero(meta.Numeric(X)), y);
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var ix: usize = 0;

        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                opInto(&o.data[numeric.cast(usize, io)], numeric.zero(meta.Numeric(X)), y);

                io += o.inc;
            }

            opInto(&o.data[numeric.cast(usize, io)], x.data[ix], y);

            i += 1;
            io += o.inc;
        }

        while (i < o.len) : (i += 1) {
            opInto(&o.data[numeric.cast(usize, io)], numeric.zero(meta.Numeric(X)), y);

            io += o.inc;
        }
    }
}
