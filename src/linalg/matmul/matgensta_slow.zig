const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

pub fn matmulInto(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    if (comptime meta.layoutOf(O) == .col_major) {
        inline for (0..O.cols) |j| {
            inline for (0..O.rows) |i| {
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

                numeric.set(&o.data[O._index(i, j)], sum);
            }
        }
    } else {
        inline for (0..O.rows) |i| {
            inline for (0..O.cols) |j| {
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

                numeric.set(&o.data[O._index(i, j)], sum);
            }
        }
    }
}
