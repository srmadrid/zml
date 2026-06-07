const std = @import("std");

const meta = @import("../../../meta.zig");
const int = @import("../../../int.zig");
const numeric = @import("../../../numeric.zig");
const matrix = @import("../../../matrix.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const aliased = (comptime O == Y) and std.meta.eql(o.*, y);

    if ((comptime op_ == numeric.subInto) or !aliased) {
        if (comptime meta.layoutOf(O) == .col_major) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                var i: usize = 0;
                while (i < o.rows) : (i += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }
            }
        } else {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                var j: usize = 0;
                while (j < o.cols) : (j += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }
            }
        }
    }

    if (comptime meta.diagOf(X) == .unit) {
        var i: usize = 0;
        while (i < int.min(o.rows, o.cols)) : (i += 1) {
            if ((comptime op_ == numeric.addInto) or !aliased)
                op_(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), y.data[y._index(i, i)])
            else
                op_(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), numeric.neg(y.data[y._index(i, i)]));
        }
    }

    if (comptime meta.layoutOf(X) == .col_major) {
        var j: usize = 0;
        while (j < x.cols) : (j += 1) {
            var p: usize = x.ptr[j];
            while (p < x.ptr[j + 1]) : (p += 1) {
                if ((comptime op_ == numeric.addInto) or !aliased)
                    op_(&o.data[o._index(x.idx[p], j)], x.data[p], y.data[y._index(x.idx[p], j)])
                else
                    op_(&o.data[o._index(x.idx[p], j)], x.data[p], numeric.neg(y.data[y._index(x.idx[p], j)]));
            }
        }
    } else {
        var i: usize = 0;
        while (i < x.rows) : (i += 1) {
            var p: usize = x.ptr[i];
            while (p < x.ptr[i + 1]) : (p += 1) {
                if ((comptime op_ == numeric.addInto) or !aliased)
                    op_(&o.data[o._index(i, x.idx[p])], x.data[p], y.data[y._index(i, x.idx[p])])
                else
                    op_(&o.data[o._index(i, x.idx[p])], x.data[p], numeric.neg(y.data[y._index(i, x.idx[p])]));
            }
        }
    }
}
