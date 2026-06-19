const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

pub fn matmulInto(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var i: usize = 0;
                while (i < j) : (i += 1) {
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
                }

                var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                var k: usize = 0;
                while (k < x_cols) : (k += 1) {
                    numeric.fmaInto(
                        &sum,
                        x.get(j, k) catch unreachable,
                        y.get(k, j) catch unreachable,
                        sum,
                    );
                }

                numeric.set(&o.data[o._index(j, j)], sum);
            } else {
                var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                var k: usize = 0;
                while (k < x_cols) : (k += 1) {
                    numeric.fmaInto(
                        &sum,
                        x.get(j, k) catch unreachable,
                        y.get(k, j) catch unreachable,
                        sum,
                    );
                }

                numeric.set(&o.data[o._index(j, j)], sum);

                var i: usize = j + 1;
                while (i < o.rows) : (i += 1) {
                    sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                    k = 0;
                    while (k < x_cols) : (k += 1) {
                        numeric.fmaInto(
                            &sum,
                            x.get(i, k) catch unreachable,
                            y.get(k, j) catch unreachable,
                            sum,
                        );
                    }

                    numeric.set(&o.data[o._index(i, j)], sum);
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(O) == .lower) {
                var j: usize = 0;
                while (j < i) : (j += 1) {
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
                }

                var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                var k: usize = 0;
                while (k < x_cols) : (k += 1) {
                    numeric.fmaInto(
                        &sum,
                        x.get(i, k) catch unreachable,
                        y.get(k, i) catch unreachable,
                        sum,
                    );
                }

                numeric.set(&o.data[o._index(i, i)], sum);
            } else {
                var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                var k: usize = 0;
                while (k < x_cols) : (k += 1) {
                    numeric.fmaInto(
                        &sum,
                        x.get(i, k) catch unreachable,
                        y.get(k, i) catch unreachable,
                        sum,
                    );
                }

                numeric.set(&o.data[o._index(i, i)], sum);

                var j: usize = i + 1;
                while (j < o.cols) : (j += 1) {
                    sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                    k = 0;
                    while (k < x_cols) : (k += 1) {
                        numeric.fmaInto(
                            &sum,
                            x.get(i, k) catch unreachable,
                            y.get(k, j) catch unreachable,
                            sum,
                        );
                    }

                    numeric.set(&o.data[o._index(i, j)], sum);
                }
            }
        }
    }
}
