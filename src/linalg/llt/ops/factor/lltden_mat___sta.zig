const linalg = @import("../../../../linalg.zig");
const matrix = @import("../../../../matrix.zig");
const meta = @import("../../../../meta.zig");
const numeric = @import("../../../../numeric.zig");

pub fn factorIntoUnchecked(o: anytype, x: anytype) !void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    if (comptime O.storage_layout == .col_major) {
        if (comptime meta.uploOf(X) == .lower) {
            inline for (0..X.cols) |j| {
                inline for (j..X.rows) |i| {
                    numeric.set(
                        &o.data[i + j * o.ld],
                        x.data[X._index(i, j)],
                    );
                }
            }
        } else {
            inline for (0..X.cols) |j| {
                inline for (j..X.rows) |i| {
                    numeric.conjInto(
                        &o.data[i + j * o.ld],
                        x.data[X._index(j, i)],
                    );
                }
            }
        }
    } else {
        if (comptime meta.uploOf(X) == .lower) {
            inline for (0..X.rows) |i| {
                inline for (0..i + 1) |j| {
                    numeric.set(
                        &o.data[i * o.ld + j],
                        x.data[X._index(i, j)],
                    );
                }
            }
        } else {
            inline for (0..X.rows) |i| {
                inline for (0..i + 1) |j| {
                    numeric.conjInto(
                        &o.data[i * o.ld + j],
                        x.data[X._index(j, i)],
                    );
                }
            }
        }
    }

    const info = linalg.lapack.potrf(
        O.storage_layout,
        .lower,
        o.rows,
        o.data,
        o.ld,
    ) catch unreachable;

    if (info != numeric.highest(usize))
        return linalg.Error.IndefiniteMatrix;

    return;
}
