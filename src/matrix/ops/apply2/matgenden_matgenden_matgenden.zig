const meta = @import("../../../meta.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    if (comptime meta.layoutOf(meta.Child(@TypeOf(o))) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
            }
        }
    }
}
