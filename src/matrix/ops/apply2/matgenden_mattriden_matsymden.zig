const std = @import("std");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            if (comptime meta.uploOf(X) == .upper) {
                while (i < j) : (i += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                    opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], ty);
                }
            } else {
                while (i < j) : (i += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                    opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), ty);
                }
            }

            if (comptime meta.diagOf(X) == .unit)
                opInto(&o.data[o._index(j, j)], numeric.one(meta.Numeric(X)), y.data[y._index(j, j)])
            else
                opInto(&o.data[o._index(j, j)], x.data[x._index(j, j)], y.data[y._index(j, j)]);

            i = j + 1;
            if (comptime meta.uploOf(X) == .lower) {
                while (i < o.rows) : (i += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                    opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], ty);
                }
            } else {
                while (i < o.rows) : (i += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                    opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), ty);
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            if (comptime meta.uploOf(X) == .lower) {
                while (j < i) : (j += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                    opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], ty);
                }
            } else {
                while (j < i) : (j += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                    opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), ty);
                }
            }

            if (comptime meta.diagOf(X) == .unit)
                opInto(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), y.data[y._index(i, i)])
            else
                opInto(&o.data[o._index(i, i)], x.data[x._index(i, i)], y.data[y._index(i, i)]);

            j = i + 1;
            if (comptime meta.uploOf(X) == .upper) {
                while (j < o.cols) : (j += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                    opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], ty);
                }
            } else {
                while (j < o.cols) : (j += 1) {
                    const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                    opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), ty);
                }
            }
        }
    }
}
