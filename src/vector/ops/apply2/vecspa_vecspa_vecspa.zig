const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    var ix: usize = 0;
    var iy: usize = 0;
    var io: usize = 0;
    while (ix < x.nnz and iy < y.nnz) {
        if (x.idx[ix] == y.idx[iy]) {
            opInto(&o.data[io], x.data[ix], y.data[iy]);

            o.idx[io] = x.idx[ix];
            ix += 1;
            iy += 1;
            io += 1;
        } else if (x.idx[ix] < y.idx[iy]) {
            opInto(&o.data[io], x.data[ix], numeric.zero(meta.Numeric(Y)));

            o.idx[io] = x.idx[ix];
            ix += 1;
            io += 1;
        } else {
            opInto(&o.data[io], numeric.zero(meta.Numeric(X)), y.data[iy]);

            o.idx[io] = y.idx[iy];
            iy += 1;
            io += 1;
        }
    }

    while (ix < x.nnz) {
        opInto(&o.data[io], x.data[ix], numeric.zero(meta.Numeric(Y)));

        o.idx[io] = x.idx[ix];
        ix += 1;
        io += 1;
    }

    while (iy < y.nnz) {
        opInto(&o.data[io], numeric.zero(meta.Numeric(X)), y.data[iy]);

        o.idx[io] = y.idx[iy];
        iy += 1;
        io += 1;
    }

    o.nnz = io;
}
