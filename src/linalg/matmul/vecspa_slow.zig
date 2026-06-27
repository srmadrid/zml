const meta = @import("../../meta.zig");
const int = @import("../../int.zig");
const numeric = @import("../../numeric.zig");
const linalg = @import("../../linalg.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    if (comptime meta.isVector(X)) {
        const x_len = if (comptime meta.isStaticVector(X)) X.len else x.len;

        var nnz: usize = 0;

        var j: usize = 0;
        while (j < o.len) : (j += 1) {
            var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

            var k: usize = 0;
            while (k < x_len) : (k += 1) {
                numeric.fmaInto(
                    &sum,
                    x.get(k) catch unreachable,
                    y.get(k, j) catch unreachable,
                    sum,
                );
            }

            if (numeric.eq(sum, 0))
                continue;

            numeric.set(&o.data[nnz], sum);

            if (o.flags.noconj)
                numeric.set(&o.data[nnz], sum)
            else
                numeric.conjInto(&o.data[nnz], sum);

            o.idx[nnz] = j;
            nnz += 1;
        }

        o.nnz = nnz;
    } else {
        const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

        var nnz: usize = 0;

        var i: usize = 0;
        while (i < o.len) : (i += 1) {
            var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

            var k: usize = 0;
            while (k < x_cols) : (k += 1) {
                numeric.fmaInto(
                    &sum,
                    x.get(i, k) catch unreachable,
                    y.get(k) catch unreachable,
                    sum,
                );
            }

            if (numeric.eq(sum, 0))
                continue;

            numeric.set(&o.data[nnz], sum);

            if (o.flags.noconj)
                numeric.set(&o.data[nnz], sum)
            else
                numeric.conjInto(&o.data[nnz], sum);

            o.idx[nnz] = i;
            nnz += 1;
        }

        o.nnz = nnz;
    }
}
