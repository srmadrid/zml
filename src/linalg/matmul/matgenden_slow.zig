const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                var k: usize = 0;
                while (k < x_cols) : (k += 1) {
                    numeric.fmaInto(
                        &sum,
                        x.get(i, k) catch unreachable,
                        y.get(k, j) catch unreachable,
                        sum,
                    );
                }

                numeric.set(&o.data[o._index(i, j)], sum);

                if (o.flags.noconj)
                    numeric.set(&o.data[o._index(i, j)], sum)
                else
                    numeric.conjInto(&o.data[o._index(i, j)], sum);
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                var k: usize = 0;
                while (k < x_cols) : (k += 1) {
                    numeric.fmaInto(
                        &sum,
                        x.get(i, k) catch unreachable,
                        y.get(k, j) catch unreachable,
                        sum,
                    );
                }

                numeric.set(&o.data[o._index(i, j)], sum);

                if (o.flags.noconj)
                    numeric.set(&o.data[o._index(i, j)], sum)
                else
                    numeric.conjInto(&o.data[o._index(i, j)], sum);
            }
        }
    }
}
