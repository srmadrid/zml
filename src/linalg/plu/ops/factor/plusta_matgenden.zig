const linalg = @import("../../../../linalg.zig");
const matrix = @import("../../../../matrix.zig");
const meta = @import("../../../../meta.zig");
const numeric = @import("../../../../numeric.zig");

pub fn factorIntoUnchecked(o: anytype, x: anytype) !void {
    const O: type = meta.Child(@TypeOf(o));

    if (comptime O.storage_layout == .col_major) {
        inline for (0..O.cols) |j| {
            inline for (0..O.rows) |i| {
                numeric.set(
                    &o.data[i + j * O.rows],
                    x.data[x._index(i, j)],
                );
            }
        }
    } else {
        inline for (0..O.rows) |i| {
            inline for (0..O.cols) |j| {
                numeric.set(
                    &o.data[i * O.cols + j],
                    x.data[x._index(i, j)],
                );
            }
        }
    }

    const info = linalg.lapack.getrf(
        O.storage_layout,
        O.rows,
        O.cols,
        @as([*]O.Numeric, &o.data),
        if (comptime O.storage_layout == .col_major) O.rows else O.cols,
        @as([*]usize, &o.pivots),
    ) catch unreachable;

    if (info != numeric.highest(usize))
        return linalg.Error.SingularMatrix;

    return;
}
