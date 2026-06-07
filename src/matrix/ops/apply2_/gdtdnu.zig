const meta = @import("../../../meta.zig");

const int = @import("../../../int.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            if (comptime meta.uploOf(X) == .upper) {
                while (i < int.min(j, o.rows)) : (i += 1) {
                    op_(&o.data[o._index(i, j)], x.data[x._index(i, j)], y);
                }
            } else {
                while (i < int.min(j, o.rows)) : (i += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }

            if (j < o.rows) {
                if (comptime meta.diagOf(X) == .unit) {
                    if (comptime op_ == numeric.mulInto)
                        numeric.set(&o.data[o._index(j, j)], y)
                    else
                        op_(&o.data[o._index(j, j)], numeric.one(meta.Numeric(X)), y);
                } else {
                    op_(&o.data[o._index(j, j)], x.data[x._index(j, j)], y);
                }
            }

            i = int.min(j + 1, o.rows);
            if (comptime meta.uploOf(X) == .lower) {
                while (i < o.rows) : (i += 1) {
                    op_(&o.data[o._index(i, j)], x.data[x._index(i, j)], y);
                }
            } else {
                while (i < o.rows) : (i += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            if (comptime meta.uploOf(X) == .lower) {
                while (j < int.min(i, o.cols)) : (j += 1) {
                    op_(&o.data[o._index(i, j)], x.data[x._index(i, j)], y);
                }
            } else {
                while (j < int.min(i, o.cols)) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }

            if (i < o.cols) {
                if (comptime meta.diagOf(X) == .unit) {
                    if (comptime op_ == numeric.mulInto)
                        numeric.set(&o.data[o._index(i, i)], y)
                    else
                        op_(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), y);
                } else {
                    op_(&o.data[o._index(i, i)], x.data[x._index(i, i)], y);
                }
            }

            j = int.min(i + 1, o.cols);
            if (comptime meta.uploOf(X) == .upper) {
                while (j < o.cols) : (j += 1) {
                    op_(&o.data[o._index(i, j)], x.data[x._index(i, j)], y);
                }
            } else {
                while (j < o.cols) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }
        }
    }
}
