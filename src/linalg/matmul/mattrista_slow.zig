const meta = @import("../../meta.zig");
const int = @import("../../int.zig");
const numeric = @import("../../numeric.zig");

pub fn matmulInto(o: anytype, x: anytype, y: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    if (comptime meta.layoutOf(O) == .col_major) {
        inline for (0..O.cols) |j| {
            if (comptime meta.uploOf(O) == .upper) {
                inline for (0..(comptime int.min(if (meta.diagOf(O) == .non_unit) j + 1 else j, O.rows))) |i| {
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
            } else {
                inline for ((comptime int.min(if (meta.diagOf(O) == .non_unit) j else j + 1, O.rows))..O.rows) |i| {
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
    } else {
        inline for (0..O.rows) |i| {
            if (comptime meta.uploOf(O) == .lower) {
                inline for (0..(comptime int.min(if (meta.diagOf(O) == .non_unit) i + 1 else i, O.cols))) |j| {
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
            } else {
                inline for ((comptime int.min(if (meta.diagOf(O) == .non_unit) i else i + 1, O.cols))..O.cols) |j| {
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
}
