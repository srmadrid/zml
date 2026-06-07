const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
            }

            if (j < o.rows) {
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(j, j)], y.data[j])
                else
                    numeric.set(&o.data[o._index(j, j)], numeric.neg(y.data[j]));
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
            }

            if (i < o.cols) {
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(i, i)], y.data[i])
                else
                    numeric.set(&o.data[o._index(i, i)], numeric.neg(y.data[i]));
            }
        }
    }

    var k: usize = 0;
    while (k < x.rows) : (k += 1) {
        const i = if (x.direction == .forward) k else x.data[k];
        const j = if (x.direction == .forward) x.data[k] else k;

        if (i == j)
            op_(&o.data[o._index(i, j)], numeric.one(meta.Numeric(X)), y.data[i])
        else
            o.data[o._index(i, j)] = numeric.one(meta.Numeric(O));
    }
}
