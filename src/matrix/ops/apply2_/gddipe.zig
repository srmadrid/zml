const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const Y = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
            }

            if (j < o.rows) {
                numeric.set(&o.data[o._index(j, j)], x.data[j]);
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
                numeric.set(&o.data[o._index(i, i)], x.data[i]);
            }
        }
    }

    var k: usize = 0;
    while (k < y.rows) : (k += 1) {
        const i = if (y.direction == .forward) k else y.data[k];
        const j = if (y.direction == .forward) y.data[k] else k;

        if (i == j) {
            op_(&o.data[o._index(i, j)], x.data[i], numeric.one(meta.Numeric(Y)));
        } else {
            if (comptime op_ == numeric.addInto)
                o.data[o._index(i, j)] = numeric.one(meta.Numeric(O))
            else
                o.data[o._index(i, j)] = numeric.neg(numeric.one(meta.Numeric(O)));
        }
    }
}
