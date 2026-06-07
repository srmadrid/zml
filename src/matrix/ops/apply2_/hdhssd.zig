const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const matrix = @import("../../../matrix.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var i: usize = 0;
                while (i <= j) : (i += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], ty)
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
                }
            } else {
                var i: usize = j;
                while (i < o.rows) : (i += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], ty)
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(O) == .lower) {
                var j: usize = 0;
                while (j <= i) : (j += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], ty)
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
                }
            } else {
                var j: usize = i;
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

    if (comptime meta.layoutOf(X) == .col_major) {
        var j: usize = 0;
        while (j < x.cols) : (j += 1) {
            var p: usize = x.ptr[j];
            while (p < x.ptr[j + 1]) : (p += 1) {
                if (comptime meta.uploOf(O) == meta.uploOf(X)) {
                    const ty = if (comptime meta.uploOf(Y) == meta.uploOf(O))
                        y.data[y._index(x.idx[p], j)]
                    else
                        numeric.conj(y.data[y._index(j, x.idx[p])]);

                    op_(&o.data[o._index(x.idx[p], j)], x.data[p], ty);
                } else {
                    const ty = if (comptime meta.uploOf(Y) == meta.uploOf(O))
                        y.data[y._index(j, x.idx[p])]
                    else
                        numeric.conj(y.data[y._index(x.idx[p], j)]);

                    op_(&o.data[o._index(j, x.idx[p])], numeric.conj(x.data[p]), ty);
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < x.rows) : (i += 1) {
            var p: usize = x.ptr[i];
            while (p < x.ptr[i + 1]) : (p += 1) {
                if (comptime meta.uploOf(O) == meta.uploOf(X)) {
                    const ty = if (comptime meta.uploOf(Y) == meta.uploOf(O))
                        y.data[y._index(i, x.idx[p])]
                    else
                        numeric.conj(y.data[y._index(x.idx[p], i)]);

                    op_(&o.data[o._index(i, x.idx[p])], x.data[p], ty);
                } else {
                    const ty = if (comptime meta.uploOf(Y) == meta.uploOf(O))
                        y.data[y._index(x.idx[p], i)]
                    else
                        numeric.conj(y.data[y._index(i, x.idx[p])]);

                    op_(&o.data[o._index(x.idx[p], i)], numeric.conj(x.data[p]), ty);
                }
            }
        }
    }
}
