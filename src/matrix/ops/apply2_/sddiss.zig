const std = @import("std");

const meta = @import("../../../meta.zig");
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
                while (i < j) : (i += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }

                numeric.set(&o.data[o._index(j, j)], x.data[j]);
            } else {
                numeric.set(&o.data[o._index(j, j)], x.data[j]);

                var i: usize = j + 1;
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
                while (j < i) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }

                numeric.set(&o.data[o._index(i, i)], x.data[i]);
            } else {
                numeric.set(&o.data[o._index(i, i)], x.data[i]);

                var j: usize = i + 1;
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
                if (comptime meta.uploOf(O) == meta.uploOf(Y)) {
                    if (y.idx[p] == j) {
                        op_(&o.data[o._index(y.idx[p], j)], x.data[j], y.data[p]);
                    } else {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(y.idx[p], j)], y.data[p])
                        else
                            numeric.set(&o.data[o._index(y.idx[p], j)], numeric.neg(y.data[p]));
                    }
                } else {
                    if (y.idx[p] == j) {
                        op_(&o.data[o._index(j, y.idx[p])], x.data[j], y.data[p]);
                    } else {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(j, y.idx[p])], y.data[p])
                        else
                            numeric.set(&o.data[o._index(j, y.idx[p])], numeric.neg(y.data[p]));
                    }
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < y.rows) : (i += 1) {
            var p: usize = y.ptr[i];
            while (p < y.ptr[i + 1]) : (p += 1) {
                if (comptime meta.uploOf(O) == meta.uploOf(Y)) {
                    if (i == y.idx[p]) {
                        op_(&o.data[o._index(i, y.idx[p])], x.data[i], y.data[p]);
                    } else {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, y.idx[p])], y.data[p])
                        else
                            numeric.set(&o.data[o._index(i, y.idx[p])], numeric.neg(y.data[p]));
                    }
                } else {
                    if (i == y.idx[p]) {
                        op_(&o.data[o._index(y.idx[p], i)], x.data[i], y.data[p]);
                    } else {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(y.idx[p], i)], y.data[p])
                        else
                            numeric.set(&o.data[o._index(y.idx[p], i)], numeric.neg(y.data[p]));
                    }
                }
            }
        }
    }
}
