const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(Y) == .upper) {
                var i: usize = 0;
                while (i < j) : (i += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }

                const ty = if (comptime meta.diagOf(Y) == .unit) numeric.one(meta.Numeric(Y)) else y.data[y._index(j, j)];
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(j, j)], ty)
                else
                    numeric.set(&o.data[o._index(j, j)], numeric.neg(ty));

                i = j + 1;
                while (i < o.rows) : (i += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            } else {
                var i: usize = 0;
                while (i < j) : (i += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }

                const ty = if (comptime meta.diagOf(Y) == .unit) numeric.one(meta.Numeric(Y)) else y.data[y._index(j, j)];
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(j, j)], ty)
                else
                    numeric.set(&o.data[o._index(j, j)], numeric.neg(ty));

                i = j + 1;
                while (i < o.rows) : (i += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(Y) == .lower) {
                var j: usize = 0;
                while (j < i) : (j += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }

                const ty = if (comptime meta.diagOf(Y) == .unit) numeric.one(meta.Numeric(Y)) else y.data[y._index(i, i)];
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(i, i)], ty)
                else
                    numeric.set(&o.data[o._index(i, i)], numeric.neg(ty));

                j = i + 1;
                while (j < o.cols) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }
            } else {
                var j: usize = 0;
                while (j < i) : (j += 1) {
                    o.data[o._index(i, j)] = numeric.zero(meta.Numeric(O));
                }

                const ty = if (comptime meta.diagOf(Y) == .unit) numeric.one(meta.Numeric(Y)) else y.data[y._index(i, i)];
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(i, i)], ty)
                else
                    numeric.set(&o.data[o._index(i, i)], numeric.neg(ty));

                j = i + 1;
                while (j < o.cols) : (j += 1) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                    else
                        numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                }
            }
        }
    }

    var k: usize = 0;
    while (k < x.rows) : (k += 1) {
        const i = if (x.direction == .forward) k else x.data[k];
        const j = if (x.direction == .forward) x.data[k] else k;

        const ty = if (i == j)
            (if (comptime meta.diagOf(Y) == .unit) numeric.one(meta.Numeric(Y)) else y.data[y._index(i, i)])
        else if (i < j)
            (if (comptime meta.uploOf(Y) == .upper) y.data[y._index(i, j)] else numeric.zero(meta.Numeric(Y)))
        else
            (if (comptime meta.uploOf(Y) == .lower) y.data[y._index(i, j)] else numeric.zero(meta.Numeric(Y)));

        op_(&o.data[o._index(i, j)], numeric.one(meta.Numeric(X)), ty);
    }
}
