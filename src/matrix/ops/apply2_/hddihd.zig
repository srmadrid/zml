const std = @import("std");

const meta = @import("../../../meta.zig");

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
                    while (i < j) : (i += 1) {
                        const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], ty)
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
                    }
                }

                op_(&o.data[o._index(j, j)], x.data[j], y.data[y._index(j, j)]);
            } else {
                op_(&o.data[o._index(j, j)], x.data[j], y.data[y._index(j, j)]);

                if ((comptime op_ == numeric.subInto) or !aliased) {
                    var i: usize = j + 1;
                    while (i < o.rows) : (i += 1) {
                        const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], ty)
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
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
                    while (j < i) : (j += 1) {
                        const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], ty)
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
                    }
                }

                op_(&o.data[o._index(i, i)], x.data[i], y.data[y._index(i, i)]);
            } else {
                op_(&o.data[o._index(i, i)], x.data[i], y.data[y._index(i, i)]);

                if ((comptime op_ == numeric.subInto) or !aliased) {
                    var j: usize = i + 1;
                    while (j < o.cols) : (j += 1) {
                        const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], ty)
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
                    }
                }
            }
        }
    }
}
