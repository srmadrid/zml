const numeric = @import("../../../numeric.zig");
const vector = @import("../../../vector.zig");

const int = @import("../../../int.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) !void {
    var nnz: usize = 0;
    var cx: usize = 0;
    var cy: usize = 0;
    while (cx < x.nnz and cy < y.nnz) {
        if (x.idx[cx] == y.idx[cy]) {
            nnz += 1;
            cx += 1;
            cy += 1;
        } else if (x.idx[cx] < y.idx[cy]) {
            nnz += 1;
            cx += 1;
        } else {
            nnz += 1;
            cy += 1;
        }
    }

    nnz += (x.nnz - cx) + (y.nnz - cy);

    if (o._dlen < nnz or o._ilen < nnz)
        return vector.Error.DimensionMismatch;

    o.nnz = nnz;

    var ix: usize = x.nnz;
    var iy: usize = y.nnz;
    var io: usize = nnz;
    while (ix > 0 and iy > 0) {
        const xi = ix - 1;
        const yi = iy - 1;
        const out_i = io - 1;

        if (x.idx[xi] == y.idx[yi]) {
            opInto(&o.data[out_i], x.data[xi], y.data[yi]);
            o.idx[out_i] = x.idx[xi];
            ix -= 1;
            iy -= 1;
            io -= 1;
        } else if (x.idx[xi] > y.idx[yi]) {
            numeric.set(&o.data[out_i], x.data[xi]);

            o.idx[out_i] = x.idx[xi];
            ix -= 1;
            io -= 1;
        } else {
            if (comptime opInto == numeric.addInto)
                numeric.set(&o.data[out_i], y.data[yi])
            else
                numeric.set(&o.data[out_i], numeric.neg(y.data[yi]));

            o.idx[out_i] = y.idx[yi];
            iy -= 1;
            io -= 1;
        }
    }

    while (ix > 0) {
        const xi = ix - 1;
        const oi = io - 1;

        numeric.set(&o.data[oi], x.data[xi]);

        o.idx[oi] = x.idx[xi];
        ix -= 1;
        io -= 1;
    }

    while (iy > 0) {
        const yi = iy - 1;
        const oi = io - 1;

        if (comptime opInto == numeric.addInto)
            numeric.set(&o.data[oi], y.data[yi])
        else
            numeric.set(&o.data[oi], numeric.neg(y.data[yi]));

        o.idx[oi] = y.idx[yi];
        iy -= 1;
        io -= 1;
    }
}
