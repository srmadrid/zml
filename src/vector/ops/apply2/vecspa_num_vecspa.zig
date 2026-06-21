pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    o.nnz = y.nnz;

    var i: usize = 0;
    while (i < o.nnz) : (i += 1) {
        opInto(&o.data[i], x, y.data[i]);

        o.idx[i] = y.idx[i];
    }
}
