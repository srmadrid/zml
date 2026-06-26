const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    if (comptime meta.isVector(X)) {
        const x_len = if (comptime meta.isStaticVector(X)) X.len else x.len;

        inline for (0..O.len) |j| {
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

            numeric.set(&o.data[j], sum);
        }
    } else {
        const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

        inline for (0..O.len) |i| {
            var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

            for (0..x_cols) |k| {
                numeric.fmaInto(
                    &sum,
                    x.get(i, k) catch unreachable,
                    y.get(k) catch unreachable, // Y is the vector
                    sum,
                );
            }

            numeric.set(&o.data[i], sum);
        }
    }
}
