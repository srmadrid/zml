const meta = @import("../../../meta.zig");

const int = @import("../../../int.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var i: usize = 0;
                while (i < int.min(if (comptime meta.diagOf(O) == .non_unit) j + 1 else j, o.rows)) : (i += 1) {
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
                }
            } else {
                var i: usize = int.min(if (comptime meta.diagOf(O) == .non_unit) j else j + 1, o.rows);
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
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(O) == .lower) {
                var j: usize = 0;
                while (j < int.min(if (comptime meta.diagOf(O) == .non_unit) i + 1 else i, o.cols)) : (j += 1) {
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
                }
            } else {
                var j: usize = int.min(if (comptime meta.diagOf(O) == .non_unit) i else i + 1, o.cols);
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
                }
            }
        }
    }
}
