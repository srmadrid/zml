const std = @import("std");

const meta = @import("../../../meta.zig");

const int = @import("../../../int.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < int.min(j, o.rows)) : (i += 1) {
                opInto(&o.data[o._index(j, j)], numeric.zero(meta.Numeric(X)), y.data[y._index(j, j)]);
            }

            if (j < o.rows) {
                opInto(&o.data[o._index(j, j)], x.data[j], y.data[y._index(j, j)]);
            }

            i = int.min(j + 1, o.rows);
            while (i < o.rows) : (i += 1) {
                opInto(&o.data[o._index(j, j)], numeric.zero(meta.Numeric(X)), y.data[y._index(j, j)]);
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < int.min(i, o.cols)) : (j += 1) {
                opInto(&o.data[o._index(j, j)], numeric.zero(meta.Numeric(X)), y.data[y._index(j, j)]);
            }

            if (i < o.cols) {
                opInto(&o.data[o._index(i, i)], x.data[i], y.data[y._index(i, i)]);
            }

            j = int.min(i + 1, o.cols);
            while (j < o.cols) : (j += 1) {
                opInto(&o.data[o._index(j, j)], numeric.zero(meta.Numeric(X)), y.data[y._index(j, j)]);
            }
        }
    }
}
