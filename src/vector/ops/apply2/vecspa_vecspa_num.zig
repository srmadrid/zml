pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    o.nnz = x.nnz;

    var i: usize = 0;
    while (i < x.nnz) : (i += 1) {
        opInto(&o.data[i], x.data[i], y);

        o.idx[i] = x.idx[i];
    }
}
