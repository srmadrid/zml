const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");
const linalg = @import("../../linalg.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O: type = meta.Child(@TypeOf(x));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    var nnz: usize = 0;

    const is_col_major = comptime blk: {
        if (meta.isMatrix(X) and meta.layoutOf(X) != null)
            break :blk meta.layoutOf(X).? == .col_major
        else if (meta.isMatrix(Y) and meta.layoutOf(Y) != null)
            break :blk meta.layoutOf(Y).? == .col_major;

        break :blk true;
    };

    if (comptime is_col_major) {
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

                if (numeric.eq(sum, 0))
                    continue;

                numeric.set(&o.data[nnz], sum);

                o.ridx[nnz] = i;
                o.cidx[nnz] = j;
                nnz += 1;
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

                if (numeric.eq(sum, 0))
                    continue;

                numeric.set(&o.data[nnz], sum);

                o.ridx[nnz] = i;
                o.cidx[nnz] = j;
                nnz += 1;
            }
        }
    }

    o.nnz = nnz;
}
