const meta = @import("../../../meta.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X: type = @TypeOf(x);

    if (comptime meta.layoutOf(meta.Child(@TypeOf(o))) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < j) : (i += 1) {
                const tx = if (comptime meta.uploOf(X) == .upper) x.data[x._index(i, j)] else x.data[x._index(j, i)];
                opInto(&o.data[o._index(i, j)], tx, y.data[y._index(i, j)]);
            }

            opInto(&o.data[o._index(j, j)], x.data[x._index(j, j)], y.data[y._index(j, j)]);

            i = j + 1;
            while (i < o.rows) : (i += 1) {
                const tx = if (comptime meta.uploOf(X) == .lower) x.data[x._index(i, j)] else x.data[x._index(j, i)];
                opInto(&o.data[o._index(i, j)], tx, y.data[y._index(i, j)]);
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < i) : (j += 1) {
                const tx = if (comptime meta.uploOf(X) == .lower) x.data[x._index(i, j)] else x.data[x._index(j, i)];
                opInto(&o.data[o._index(i, j)], tx, y.data[y._index(i, j)]);
            }

            opInto(&o.data[o._index(i, i)], x.data[x._index(i, i)], y.data[y._index(i, i)]);

            j = i + 1;
            while (j < o.cols) : (j += 1) {
                const tx = if (comptime meta.uploOf(X) == .upper) x.data[x._index(i, j)] else x.data[x._index(j, i)];
                opInto(&o.data[o._index(i, j)], tx, y.data[y._index(i, j)]);
            }
        }
    }
}
