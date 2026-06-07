const std = @import("std");

const meta = @import("../../../meta.zig");

const int = @import("../../../int.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const Y: type = @TypeOf(y);

    const aliased = (comptime O == Y) and std.meta.eql(o.*, y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                if ((comptime op_ == numeric.subInto) or !aliased) {
                    var i: usize = 0;
                    while (i < int.min(j, o.rows)) : (i += 1) {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                    }
                }

                if (j < o.rows) {
                    if (comptime meta.diagOf(Y) == .unit)
                        op_(&o.data[o._index(j, j)], x.data[j], numeric.one(meta.Numeric(Y)))
                    else
                        op_(&o.data[o._index(j, j)], x.data[j], y.data[y._index(j, j)]);
                }
            } else {
                if (j < o.rows) {
                    if (comptime meta.diagOf(Y) == .unit)
                        op_(&o.data[o._index(j, j)], x.data[j], numeric.one(meta.Numeric(Y)))
                    else
                        op_(&o.data[o._index(j, j)], x.data[j], y.data[y._index(j, j)]);
                }

                if ((comptime op_ == numeric.subInto) or !aliased) {
                    var i: usize = int.min(j + 1, o.rows);
                    while (i < o.rows) : (i += 1) {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                    }
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(O) == .lower) {
                if ((comptime op_ == numeric.subInto) or !aliased) {
                    var j: usize = 0;
                    while (j < int.min(i, o.cols)) : (j += 1) {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                    }
                }

                if (i < o.cols) {
                    if (comptime meta.diagOf(Y) == .unit)
                        op_(&o.data[o._index(i, i)], x.data[i], numeric.one(meta.Numeric(Y)))
                    else
                        op_(&o.data[o._index(i, i)], x.data[i], y.data[y._index(i, i)]);
                }
            } else {
                if (i < o.cols) {
                    if (comptime meta.diagOf(Y) == .unit)
                        op_(&o.data[o._index(i, i)], x.data[i], numeric.one(meta.Numeric(Y)))
                    else
                        op_(&o.data[o._index(i, i)], x.data[i], y.data[y._index(i, i)]);
                }

                if ((comptime op_ == numeric.subInto) or !aliased) {
                    var j: usize = int.min(i + 1, o.cols);
                    while (j < o.cols) : (j += 1) {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                    }
                }
            }
        }
    }
}
