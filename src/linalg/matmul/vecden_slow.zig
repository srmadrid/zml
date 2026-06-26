const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    if (comptime meta.isVector(X)) {
        const x_len = if (comptime meta.isStaticVector(X)) X.len else x.len;

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

            numeric.set(&o.data[o._index(j)], sum);
        }
    } else {
        const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

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

            numeric.set(&o.data[o._index(i)], sum);
        }
    }
}
