const meta = @import("../../../meta.zig");

const int = @import("../../../int.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            if (comptime meta.uploOf(Y) == .upper) {
                while (i < int.min(j, o.rows)) : (i += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }
            } else {
                while (i < int.min(j, o.rows)) : (i += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }

            if (j < o.rows) {
                if (comptime meta.diagOf(Y) == .unit)
                    op_(&o.data[o._index(j, j)], x.data[j], numeric.one(meta.Numeric(Y)))
                else
                    op_(&o.data[o._index(j, j)], x.data[j], y.data[y._index(j, j)]);
            }

            i = int.min(j + 1, o.rows);
            if (comptime meta.uploOf(Y) == .lower) {
                while (i < o.rows) : (i += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
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
            if (comptime meta.uploOf(Y) == .lower) {
                while (j < int.min(i, o.cols)) : (j += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }
            } else {
                while (j < int.min(i, o.cols)) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }

            if (i < o.cols) {
                if (comptime meta.diagOf(Y) == .unit)
                    op_(&o.data[o._index(i, i)], x.data[i], numeric.one(meta.Numeric(Y)))
                else
                    op_(&o.data[o._index(i, i)], x.data[i], y.data[y._index(i, i)]);
            }

            j = int.min(i + 1, o.cols);
            if (comptime meta.uploOf(Y) == .upper) {
                while (j < o.cols) : (j += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }
            } else {
                while (j < o.cols) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }
        }
    }
}
