const vector = @import("../../../vector.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) !void {
    if (o._dlen < y.nnz or o._ilen < y.nnz)
        return vector.Error.DimensionMismatch;

    o.nnz = y.nnz;

    var i: usize = 0;
    while (i < o.nnz) : (i += 1) {
        opInto(&o.data[i], x, y.data[i]);

        o.idx[i] = y.idx[i];
    }
}
