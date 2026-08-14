const meta = @import("../../meta.zig");
const int = @import("../../int.zig");
const numeric = @import("../../numeric.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var i: usize = 0;
                while (i < int.min(if (comptime meta.diagOf(O) == .non_unit) j + 1 else j, o.rows)) : (i += 1) {
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

                    numeric.set(&o.data[o._index(i, j)], sum);

                    if (o.flags.noconj)
                        numeric.set(&o.data[o._index(i, j)], sum)
                    else
                        numeric.conjInto(&o.data[o._index(i, j)], sum);
                }
            } else {
                var i: usize = int.min(if (comptime meta.diagOf(O) == .non_unit) j else j + 1, o.rows);
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

                    numeric.set(&o.data[o._index(i, j)], sum);

                    if (o.flags.noconj)
                        numeric.set(&o.data[o._index(i, j)], sum)
                    else
                        numeric.conjInto(&o.data[o._index(i, j)], sum);
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(O) == .lower) {
                var j: usize = 0;
                while (j < int.min(if (comptime meta.diagOf(O) == .non_unit) i + 1 else i, o.cols)) : (j += 1) {
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

                    numeric.set(&o.data[o._index(i, j)], sum);

                    if (o.flags.noconj)
                        numeric.set(&o.data[o._index(i, j)], sum)
                    else
                        numeric.conjInto(&o.data[o._index(i, j)], sum);
                }
            } else {
                var j: usize = int.min(if (comptime meta.diagOf(O) == .non_unit) i else i + 1, o.cols);
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

                    numeric.set(&o.data[o._index(i, j)], sum);

                    if (o.flags.noconj)
                        numeric.set(&o.data[o._index(i, j)], sum)
                    else
                        numeric.conjInto(&o.data[o._index(i, j)], sum);
                }
            }
        }
    }
}
