const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(meta.Child(@TypeOf(o))) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime meta.isMatrix(X))
                        x.get(i, j) catch unreachable
                    else
                        x,
                    if (comptime meta.isMatrix(Y))
                        y.get(i, j) catch unreachable
                    else
                        y,
                );

                if (!o.flags.noconj)
                    numeric.conjInto(&o.data[o._index(i, j)], o.data[o._index(i, j)]);
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime meta.isMatrix(X))
                        x.get(i, j) catch unreachable
                    else
                        x,
                    if (comptime meta.isMatrix(Y))
                        y.get(i, j) catch unreachable
                    else
                        y,
                );

                if (!o.flags.noconj)
                    numeric.conjInto(&o.data[o._index(i, j)], o.data[o._index(i, j)]);
            }
        }
    }
}
