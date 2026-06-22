const int = @import("../../../int.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < int.min(j, o.rows)) : (i += 1) {
                opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), y);
            }

            if (j < o.rows) {
                opInto(&o.data[o._index(j, j)], x.data[j], y);
            }

            i = int.min(j + 1, o.rows);
            while (i < o.rows) : (i += 1) {
                opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), y);
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < int.min(i, o.cols)) : (j += 1) {
                opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), y);
            }

            if (i < o.cols) {
                opInto(&o.data[o._index(i, i)], x.data[i], y);
            }

            j = int.min(i + 1, o.cols);
            while (j < o.cols) : (j += 1) {
                opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), y);
            }
        }
    }
}
