const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const X: type = @TypeOf(x);
    const O: type = meta.Child(@TypeOf(o));

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    var nnz: usize = 0;
    o.ptr[0] = 0;

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                var sum = numeric.cast(meta.Accumulator(meta.Numeric(O)), 0);

                var k: usize = 0;
                while (k < x_cols) : (k += 1) {
                    numeric.fmaInto(
                        &sum,
                        x.get(i, k) catch unreachable,
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

                o.idx[nnz] = i;
                nnz += 1;
            }

            o.ptr[j + 1] = nnz;
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                var sum = numeric.cast(meta.Accumulator(meta.Numeric(O)), 0);

                var k: usize = 0;
                while (k < x_cols) : (k += 1) {
                    numeric.fmaInto(
                        &sum,
                        x.get(i, k) catch unreachable,
                        y.get(k, j) catch unreachable,
                        sum,
                    );
                }

                if (numeric.eq(sum, 0))
                    continue;

                numeric.set(&o.data[nnz], sum);

                o.idx[nnz] = j;
                nnz += 1;
            }

            o.ptr[i + 1] = nnz;
        }
    }

    o.nnz = nnz;
}
