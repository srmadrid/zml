const meta = @import("../../../meta.zig");

const int = @import("../../../int.zig");

const numeric = @import("../../../numeric.zig");

const matrix = @import("../../../matrix.zig");

const utils = @import("utils.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    o.setAll(numeric.zero(meta.Numeric(O)));

    if (comptime meta.layoutOf(X) == meta.layoutOf(Y)) {
        var outer: usize = 0;
        while (outer < if (comptime meta.layoutOf(X) == .col_major) x.cols else x.rows) : (outer += 1) {
            var px = x.ptr[outer];
            var py = y.ptr[outer];
            while (px < x.ptr[outer + 1] and py < y.ptr[outer + 1]) {
                if (x.idx[px] == y.idx[py]) {
                    const i_o = if (comptime meta.layoutOf(X) == .col_major) x.idx[px] else outer;
                    const j_o = if (comptime meta.layoutOf(X) == .col_major) outer else x.idx[px];

                    op_(&o.data[o._index(i_o, j_o)], x.data[px], y.data[py]);

                    px += 1;
                    py += 1;
                } else if (x.idx[px] < y.idx[py]) {
                    const i_o = if (comptime meta.layoutOf(X) == .col_major) x.idx[px] else outer;
                    const j_o = if (comptime meta.layoutOf(X) == .col_major) outer else x.idx[px];

                    if (i_o == j_o and comptime meta.diagOf(Y) == .unit) {
                        op_(&o.data[o._index(i_o, j_o)], x.data[px], numeric.one(meta.Numeric(Y)));
                    } else {
                        numeric.set(&o.data[o._index(i_o, j_o)], x.data[px]);
                    }

                    px += 1;
                } else {
                    const i_o = if (comptime meta.layoutOf(Y) == .col_major) y.idx[py] else outer;
                    const j_o = if (comptime meta.layoutOf(Y) == .col_major) outer else y.idx[py];

                    if (i_o == j_o and comptime meta.diagOf(X) == .unit) {
                        op_(&o.data[o._index(i_o, j_o)], numeric.one(meta.Numeric(X)), y.data[py]);
                    } else {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i_o, j_o)], y.data[py])
                        else
                            numeric.set(&o.data[o._index(i_o, j_o)], numeric.neg(y.data[py]));
                    }

                    py += 1;
                }
            }

            while (px < x.ptr[outer + 1]) : (px += 1) {
                const i_o = if (comptime meta.layoutOf(X) == .col_major) x.idx[px] else outer;
                const j_o = if (comptime meta.layoutOf(X) == .col_major) outer else x.idx[px];

                if (i_o == j_o and comptime meta.diagOf(Y) == .unit) {
                    op_(&o.data[o._index(i_o, j_o)], x.data[px], numeric.one(meta.Numeric(Y)));
                } else {
                    numeric.set(&o.data[o._index(i_o, j_o)], x.data[px]);
                }
            }

            while (py < y.ptr[outer + 1]) : (py += 1) {
                const i_o = if (comptime meta.layoutOf(Y) == .col_major) y.idx[py] else outer;
                const j_o = if (comptime meta.layoutOf(Y) == .col_major) outer else y.idx[py];

                if (i_o == j_o and comptime meta.diagOf(X) == .unit) {
                    op_(&o.data[o._index(i_o, j_o)], numeric.one(meta.Numeric(X)), y.data[py]);
                } else {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i_o, j_o)], y.data[py])
                    else
                        numeric.set(&o.data[o._index(i_o, j_o)], numeric.neg(y.data[py]));
                }
            }
        }
    } else {
        var idx_x_outer: usize = 0;
        while (idx_x_outer < if (comptime meta.layoutOf(X) == .col_major) x.cols else x.rows) : (idx_x_outer += 1) {
            var px = x.ptr[idx_x_outer];
            while (px < x.ptr[idx_x_outer + 1]) : (px += 1) {
                const i_o = if (comptime meta.layoutOf(X) == .col_major) x.idx[px] else idx_x_outer;
                const j_o = if (comptime meta.layoutOf(X) == .col_major) idx_x_outer else x.idx[px];

                if (i_o == j_o and comptime meta.diagOf(Y) == .unit) {
                    op_(&o.data[o._index(i_o, j_o)], x.data[px], numeric.one(meta.Numeric(Y)));
                } else if (utils.searchSparse(y, i_o, j_o)) |ty| {
                    op_(&o.data[o._index(i_o, j_o)], x.data[px], ty);
                } else {
                    numeric.set(&o.data[o._index(i_o, j_o)], x.data[px]);
                }
            }
        }

        var idx_y_outer: usize = 0;
        while (idx_y_outer < if (comptime meta.layoutOf(Y) == .col_major) y.cols else y.rows) : (idx_y_outer += 1) {
            var py = y.ptr[idx_y_outer];
            while (py < y.ptr[idx_y_outer + 1]) : (py += 1) {
                const i_o = if (comptime meta.layoutOf(Y) == .col_major) y.idx[py] else idx_y_outer;
                const j_o = if (comptime meta.layoutOf(Y) == .col_major) idx_y_outer else y.idx[py];

                if (i_o == j_o and comptime meta.diagOf(X) == .unit) {
                    op_(&o.data[o._index(i_o, j_o)], numeric.one(meta.Numeric(X)), y.data[py]);
                } else if (utils.searchSparse(x, i_o, j_o) == null) {
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i_o, j_o)], y.data[py])
                    else
                        numeric.set(&o.data[o._index(i_o, j_o)], numeric.neg(y.data[py]));
                }
            }
        }
    }

    var idx: usize = 0;
    while (idx < int.min(o.rows, o.cols)) : (idx += 1) {
        if (utils.searchSparse(x, idx, idx) == null and utils.searchSparse(y, idx, idx) == null) {
            if (comptime meta.diagOf(X) == .unit or meta.diagOf(Y) == .unit) {
                const tx = if (comptime meta.diagOf(X) == .unit) numeric.one(meta.Numeric(X)) else numeric.zero(meta.Numeric(X));
                const ty = if (comptime meta.diagOf(Y) == .unit) numeric.one(meta.Numeric(Y)) else numeric.zero(meta.Numeric(Y));
                op_(&o.data[o._index(idx, idx)], tx, ty);
            }
        }
    }
}
