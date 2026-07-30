const linalg = @import("../../../../linalg.zig");
const matrix = @import("../../../../matrix.zig");
const meta = @import("../../../../meta.zig");
const numeric = @import("../../../../numeric.zig");

pub fn factorIntoUnchecked(o: anytype, x: anytype) !void {
    const O: type = meta.Child(@TypeOf(o));

    if (comptime O.storage_layout == .col_major) {
        if (x.flags.noconj) {
            for (0..o.cols) |j| {
                for (0..o.rows) |i| {
                    numeric.set(
                        &o.data[i + j * o.ld],
                        x.data[x._index(i, j)],
                    );
                }
            }
        } else {
            for (0..o.cols) |j| {
                for (0..o.rows) |i| {
                    numeric.conjInto(
                        &o.data[i + j * o.ld],
                        x.data[x._index(i, j)],
                    );
                }
            }
        }
    } else {
        if (x.flags.noconj) {
            for (0..o.rows) |i| {
                for (0..o.cols) |j| {
                    numeric.set(
                        &o.data[i * o.ld + j],
                        x.data[x._index(i, j)],
                    );
                }
            }
        } else {
            for (0..o.rows) |i| {
                for (0..o.cols) |j| {
                    numeric.conjInto(
                        &o.data[i * o.ld + j],
                        x.data[x._index(i, j)],
                    );
                }
            }
        }
    }

    const info = linalg.lapack.getrf(
        O.storage_layout,
        o.rows,
        o.cols,
        o.data,
        o.ld,
        o.pivots,
    ) catch unreachable;

    if (info != numeric.highest(usize))
        return linalg.Error.SingularMatrix;

    return;
}
