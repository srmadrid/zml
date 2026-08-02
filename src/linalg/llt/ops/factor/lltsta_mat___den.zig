const linalg = @import("../../../../linalg.zig");
const matrix = @import("../../../../matrix.zig");
const meta = @import("../../../../meta.zig");
const numeric = @import("../../../../numeric.zig");

pub fn factorIntoUnchecked(o: anytype, x: anytype) !void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    if (comptime O.storage_layout == .col_major) {
        if ((comptime meta.uploOf(X) == .lower) == x.flags.noconj) {
            inline for (0..O.cols) |j| {
                inline for (j..O.rows) |i| {
                    numeric.set(
                        &o.data[i + j * O.size],
                        x.data[
                            if (comptime meta.uploOf(X) == .lower)
                                x._index(i, j)
                            else
                                x._index(j, i)
                        ],
                    );
                }
            }
        } else {
            inline for (0..O.cols) |j| {
                inline for (j..O.rows) |i| {
                    numeric.conjInto(
                        &o.data[i + j * O.size],
                        x.data[
                            if (comptime meta.uploOf(X) == .lower)
                                x._index(i, j)
                            else
                                x._index(j, i)
                        ],
                    );
                }
            }
        }
    } else {
        if ((comptime meta.uploOf(X) == .lower) == x.flags.noconj) {
            inline for (0..O.rows) |i| {
                inline for (0..i + 1) |j| {
                    numeric.set(
                        &o.data[i * O.size + j],
                        x.data[
                            if (comptime meta.uploOf(X) == .lower)
                                x._index(i, j)
                            else
                                x._index(j, i)
                        ],
                    );
                }
            }
        } else {
            inline for (0..O.rows) |i| {
                inline for (0..i + 1) |j| {
                    numeric.conjInto(
                        &o.data[i * O.size + j],
                        x.data[
                            if (comptime meta.uploOf(X) == .lower)
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
        .lower,
        O.size,
        @as([*]O.Numeric, &o.data),
        O.size,
    ) catch unreachable;

    if (info != numeric.highest(usize))
        return linalg.Error.IndefiniteMatrix;

    return;
}
