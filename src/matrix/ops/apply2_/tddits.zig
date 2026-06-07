const std = @import("std");

const meta = @import("../../../meta.zig");
const int = @import("../../../int.zig");
const numeric = @import("../../../numeric.zig");
const matrix = @import("../../../matrix.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var i: usize = 0;
                while (i < int.min(j, o.rows)) : (i += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }

                if (j < o.rows) {
                    if (comptime meta.diagOf(Y) == .unit)
                        op_(&o.data[o._index(j, j)], x.data[j], numeric.one(meta.Numeric(Y)))
                    else
                        numeric.set(&o.data[o._index(j, j)], x.data[j]);
                }
            } else {
                if (j < o.rows) {
                    if (comptime meta.diagOf(Y) == .unit)
                        op_(&o.data[o._index(j, j)], x.data[j], numeric.one(meta.Numeric(Y)))
                    else
                        numeric.set(&o.data[o._index(j, j)], x.data[j]);
                }

                var i: usize = int.min(j + 1, o.rows);
                while (i < o.rows) : (i += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(O) == .lower) {
                var j: usize = 0;
                while (j < int.min(i, o.cols)) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }

                if (i < o.cols) {
                    if (comptime meta.diagOf(Y) == .unit)
                        op_(&o.data[o._index(i, i)], x.data[i], numeric.one(meta.Numeric(Y)))
                    else
                        numeric.set(&o.data[o._index(i, i)], x.data[i]);
                }
            } else {
                if (i < o.cols) {
                    if (comptime meta.diagOf(Y) == .unit)
                        op_(&o.data[o._index(i, i)], x.data[i], numeric.one(meta.Numeric(Y)))
                    else
                        numeric.set(&o.data[o._index(i, i)], x.data[i]);
                }

                var j: usize = int.min(i + 1, o.cols);
                while (j < o.cols) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            }
        }
    }

    if (comptime meta.layoutOf(Y) == .col_major) {
        var j: usize = 0;
        while (j < y.cols) : (j += 1) {
            var p: usize = y.ptr[j];
            while (p < y.ptr[j + 1]) : (p += 1) {
                if (y.idx[p] == j) {
                    op_(&o.data[o._index(y.idx[p], j)], x.data[j], y.data[p]);
                } else {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(y.idx[p], j)], y.data[p])
                    else
                        numeric.set(&o.data[o._index(y.idx[p], j)], numeric.neg(y.data[p]));
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < y.rows) : (i += 1) {
            var p: usize = y.ptr[i];
            while (p < y.ptr[i + 1]) : (p += 1) {
                if (i == y.idx[p]) {
                    op_(&o.data[o._index(i, y.idx[p])], x.data[i], y.data[p]);
                } else {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, y.idx[p])], y.data[p])
                    else
                        numeric.set(&o.data[o._index(i, y.idx[p])], numeric.neg(y.data[p]));
                }
            }
        }
    }
}
