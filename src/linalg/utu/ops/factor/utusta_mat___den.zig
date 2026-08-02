const linalg = @import("../../../../linalg.zig");
const matrix = @import("../../../../matrix.zig");
const meta = @import("../../../../meta.zig");
const numeric = @import("../../../../numeric.zig");

pub fn factorIntoUnchecked(o: anytype, x: anytype) !void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    if (comptime O.storage_layout == .col_major) {
        if ((comptime meta.uploOf(X) == .upper) == x.flags.noconj) {
            inline for (0..O.cols) |j| {
                inline for (0..j + 1) |i| {
                    numeric.set(
                        &o.data[i + j * O.size],
                        x.data[
                            if (comptime meta.uploOf(X) == .upper)
                                x._index(i, j)
                            else
                                x._index(j, i)
                        ],
                    );
                }
            }
        } else {
            inline for (0..O.cols) |j| {
                inline for (0..j + 1) |i| {
                    numeric.conjInto(
                        &o.data[i + j * O.size],
                        x.data[
                            if (comptime meta.uploOf(X) == .upper)
                                x._index(i, j)
                            else
                                x._index(j, i)
                        ],
                    );
                }
            }
        }
    } else {
        if ((comptime meta.uploOf(X) == .upper) == x.flags.noconj) {
            inline for (0..O.rows) |i| {
                inline for (i..O.cols) |j| {
                    numeric.set(
                        &o.data[i * O.size + j],
                        x.data[
                            if (comptime meta.uploOf(X) == .upper)
                                x._index(i, j)
                            else
                                x._index(j, i)
                        ],
                    );
                }
            }
        } else {
            inline for (0..O.rows) |i| {
                inline for (i..O.cols) |j| {
                    numeric.conjInto(
                        &o.data[i * O.size + j],
                        x.data[
                            if (comptime meta.uploOf(X) == .upper)
                                x._index(i, j)
                            else
                                x._index(j, i)
                        ],
                    );
                }
            }
        }
    }

    const info = linalg.lapack.potrf(
        O.storage_layout,
        .upper,
        O.size,
        @as([*]O.Numeric, &o.data),
        O.size,
    ) catch unreachable;

    if (info != numeric.highest(usize))
        return linalg.Error.IndefiniteMatrix;

    return;
}
