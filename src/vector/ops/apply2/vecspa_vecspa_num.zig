const vector = @import("../../../vector.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) !void {
    if (o._dlen < x.nnz or o._ilen < x.nnz)
        return vector.Error.DimensionMismatch;

    o.nnz = x.nnz;

    var i: usize = 0;
    while (i < x.nnz) : (i += 1) {
        opInto(&o.data[i], x.data[i], y);

        o.idx[i] = x.idx[i];
    }
}
