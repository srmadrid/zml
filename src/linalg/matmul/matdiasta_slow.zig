const meta = @import("../../meta.zig");

const int = @import("../../int.zig");

const numeric = @import("../../numeric.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    inline for (0..(comptime int.min(O.rows, O.cols))) |i| {
        var sum = numeric.cast(meta.Accumulator(meta.Numeric(O)), 0);

        var k: usize = 0;
        while (k < x_cols) : (k += 1) {
            numeric.fmaInto(
                &sum,
                x.get(i, k) catch unreachable,
                y.get(k, i) catch unreachable,
                sum,
            );
        }

        numeric.set(&o.data[i], sum);
    }
}
