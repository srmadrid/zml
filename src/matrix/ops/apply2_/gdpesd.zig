const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < j) : (i += 1) {
                const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(i, j)], ty)
                else
                    numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
            }

            if (comptime op_ == numeric.addInto)
                numeric.set(&o.data[o._index(j, j)], y.data[y._index(j, j)])
            else
                numeric.set(&o.data[o._index(j, j)], numeric.neg(y.data[y._index(j, j)]));

            i = j + 1;
            while (i < o.rows) : (i += 1) {
                const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(i, j)], ty)
                else
                    numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < i) : (j += 1) {
                const ty = if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(i, j)], ty)
                else
                    numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
            }

            if (comptime op_ == numeric.addInto)
                numeric.set(&o.data[o._index(i, i)], y.data[y._index(i, i)])
            else
                numeric.set(&o.data[o._index(i, i)], numeric.neg(y.data[y._index(i, i)]));

            j = i + 1;
            while (j < o.cols) : (j += 1) {
                const ty = if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else y.data[y._index(j, i)];
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(i, j)], ty)
                else
                    numeric.set(&o.data[o._index(i, j)], numeric.neg(ty));
            }
        }
    }

    var k: usize = 0;
    while (k < x.rows) : (k += 1) {
        const i = if (x.direction == .forward) k else x.data[k];
        const j = if (x.direction == .forward) x.data[k] else k;

        const ty = if (i == j)
            y.data[y._index(i, i)]
        else if (i < j)
            (if (comptime meta.uploOf(Y) == .upper)
                y.data[y._index(i, j)]
            else
                y.data[y._index(j, i)])
        else
            (if (comptime meta.uploOf(Y) == .lower)
                y.data[y._index(i, j)]
            else
                y.data[y._index(j, i)]);

        op_(&o.data[o._index(i, j)], numeric.one(meta.Numeric(X)), ty);
    }
}
