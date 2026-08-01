const linalg = @import("../../../../linalg.zig");
const matrix = @import("../../../../matrix.zig");
const meta = @import("../../../../meta.zig");
const numeric = @import("../../../../numeric.zig");

pub fn factorIntoUnchecked(o: anytype, x: anytype) !void {
    const X: type = @TypeOf(x);

    if (comptime X.storage_layout == .col_major) {
        inline for (0..X.cols) |j| {
            inline for (0..X.rows) |i| {
                numeric.set(
                    &o.data[i + j * X.rows],
                    x.data[X._index(i, j)],
                );
            }
        }
    } else {
        inline for (0..X.rows) |i| {
            inline for (0..X.cols) |j| {
                numeric.set(
                    &o.data[i * X.cols + j],
                    x.data[X._index(i, j)],
                );
            }
        }
    }

    const info = linalg.lapack.getrf(
        X.storage_layout.?,
        X.rows,
        X.cols,
        o.data,
        if (comptime X.storage_layout == .col_major) X.rows else X.cols,
        o.pivots,
    ) catch unreachable;

    if (info != numeric.highest(usize))
        return linalg.Error.SingularMatrix;

    return;
}
