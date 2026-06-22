const meta = @import("../../../meta.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        inline for (0..O.cols) |j| {
            inline for (0..O.rows) |i| {
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
    } else {
        inline for (0..O.rows) |i| {
            inline for (0..O.cols) |j| {
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
