const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < j) : (i += 1) {
                const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                opInto(&o.data[o._index(i, j)], x, ty);
            }

            opInto(&o.data[o._index(j, j)], x, y.data[y._index(j, j)]);

            i = j + 1;
            while (i < o.rows) : (i += 1) {
                const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                opInto(&o.data[o._index(i, j)], x, ty);
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < i) : (j += 1) {
                const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                opInto(&o.data[o._index(i, j)], x, ty);
            }

            opInto(&o.data[o._index(i, i)], x, y.data[y._index(i, i)]);

            j = i + 1;
            while (j < o.cols) : (j += 1) {
                const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else numeric.conj(y.data[y._index(j, i)]);
                opInto(&o.data[o._index(i, j)], x, ty);
            }
        }
    }
}
