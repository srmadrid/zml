const int = @import("../../../int.zig");
const meta = @import("../../../meta.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        inline for (0..O.cols) |j| {
            if (comptime meta.uploOf(O) == .upper) {
                inline for (0..(comptime int.min(if (meta.diagOf(O) == .non_unit) j + 1 else j, O.rows))) |i| {
                    opInto(
                        &o.data[O._index(i, j)],
                        if (comptime meta.isMatrix(X))
                            x.get(i, j) catch unreachable
                        else
                            x,
                        if (comptime meta.isMatrix(Y))
                            y.get(i, j) catch unreachable
                        else
                            y,
                    );
                }
            } else {
                inline for ((comptime int.min(if (meta.diagOf(O) == .non_unit) j else j + 1, O.rows))..O.rows) |i| {
                    opInto(
                        &o.data[O._index(i, j)],
                        if (comptime meta.isMatrix(X))
                            x.get(i, j) catch unreachable
                        else
                            x,
                        if (comptime meta.isMatrix(Y))
                            y.get(i, j) catch unreachable
                        else
                            y,
                    );
                }
            }
        }
    } else {
        inline for (0..O.rows) |i| {
            if (comptime meta.uploOf(O) == .lower) {
                inline for (0..(comptime int.min(if (meta.diagOf(O) == .non_unit) i + 1 else i, O.cols))) |j| {
                    opInto(
                        &o.data[O._index(i, j)],
                        if (comptime meta.isMatrix(X))
                            x.get(i, j) catch unreachable
                        else
                            x,
                        if (comptime meta.isMatrix(Y))
                            y.get(i, j) catch unreachable
                        else
                            y,
                    );
                }
            } else {
                inline for ((comptime int.min(if (meta.diagOf(O) == .non_unit) i else i + 1, O.cols))..O.cols) |j| {
                    opInto(
                        &o.data[O._index(i, j)],
                        if (comptime meta.isMatrix(X))
                            x.get(i, j) catch unreachable
                        else
                            x,
                        if (comptime meta.isMatrix(Y))
                            y.get(i, j) catch unreachable
                        else
                            y,
                    );
                }
            }
        }
    }
}
