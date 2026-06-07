const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const matrix = @import("../../../matrix.zig");

const utils = @import("utils.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    o.setAll(numeric.zero(meta.Numeric(O)));

    if (comptime meta.layoutOf(X) == meta.layoutOf(Y) and meta.uploOf(X) == meta.uploOf(Y)) {
        var outer: usize = 0;
        while (outer < x.rows) : (outer += 1) {
            var px = x.ptr[outer];
            var py = y.ptr[outer];
            while (px < x.ptr[outer + 1] and py < y.ptr[outer + 1]) {
                if (x.idx[px] == y.idx[py]) {
                    const i_stored = if (comptime meta.layoutOf(X) == .col_major) x.idx[px] else outer;
                    const j_stored = if (comptime meta.layoutOf(X) == .col_major) outer else x.idx[px];
                    const i_o = if (comptime meta.uploOf(O) == meta.uploOf(X)) i_stored else j_stored;
                    const j_o = if (comptime meta.uploOf(O) == meta.uploOf(X)) j_stored else i_stored;

                    const tx = if (comptime meta.uploOf(O) == meta.uploOf(X)) x.data[px] else numeric.conj(x.data[px]);
                    const ty = if (comptime meta.uploOf(O) == meta.uploOf(Y)) y.data[py] else numeric.conj(y.data[py]);
                    op_(&o.data[o._index(i_o, j_o)], tx, ty);

                    px += 1;
                    py += 1;
                } else if (x.idx[px] < y.idx[py]) {
                    const i_stored = if (comptime meta.layoutOf(X) == .col_major) x.idx[px] else outer;
                    const j_stored = if (comptime meta.layoutOf(X) == .col_major) outer else x.idx[px];
                    const i_o = if (comptime meta.uploOf(O) == meta.uploOf(X)) i_stored else j_stored;
                    const j_o = if (comptime meta.uploOf(O) == meta.uploOf(X)) j_stored else i_stored;

                    const tx = if (comptime meta.uploOf(O) == meta.uploOf(X)) x.data[px] else numeric.conj(x.data[px]);
                    numeric.set(&o.data[o._index(i_o, j_o)], tx);

                    px += 1;
                } else {
                    const i_stored = if (comptime meta.layoutOf(Y) == .col_major) y.idx[py] else outer;
                    const j_stored = if (comptime meta.layoutOf(Y) == .col_major) outer else y.idx[py];
                    const i_o = if (comptime meta.uploOf(O) == meta.uploOf(Y)) i_stored else j_stored;
                    const j_o = if (comptime meta.uploOf(O) == meta.uploOf(Y)) j_stored else i_stored;

                    const ty = if (comptime meta.uploOf(O) == meta.uploOf(Y)) y.data[py] else numeric.conj(y.data[py]);
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i_o, j_o)], ty)
                    else
                        numeric.set(&o.data[o._index(i_o, j_o)], numeric.neg(ty));

                    py += 1;
                }
            }

            while (px < x.ptr[outer + 1]) : (px += 1) {
                const i_stored = if (comptime meta.layoutOf(X) == .col_major) x.idx[px] else outer;
                const j_stored = if (comptime meta.layoutOf(X) == .col_major) outer else x.idx[px];
                const i_o = if (comptime meta.uploOf(O) == meta.uploOf(X)) i_stored else j_stored;
                const j_o = if (comptime meta.uploOf(O) == meta.uploOf(X)) j_stored else i_stored;

                const tx = if (comptime meta.uploOf(O) == meta.uploOf(X)) x.data[px] else numeric.conj(x.data[px]);
                numeric.set(&o.data[o._index(i_o, j_o)], tx);
            }

            while (py < y.ptr[outer + 1]) : (py += 1) {
                const i_stored = if (comptime meta.layoutOf(Y) == .col_major) y.idx[py] else outer;
                const j_stored = if (comptime meta.layoutOf(Y) == .col_major) outer else y.idx[py];
                const i_o = if (comptime meta.uploOf(O) == meta.uploOf(Y)) i_stored else j_stored;
                const j_o = if (comptime meta.uploOf(O) == meta.uploOf(Y)) j_stored else i_stored;

                const ty = if (comptime meta.uploOf(O) == meta.uploOf(Y)) y.data[py] else numeric.conj(y.data[py]);
                if (comptime op_ == numeric.addInto)
                    numeric.set(&o.data[o._index(i_o, j_o)], ty)
                else
                    numeric.set(&o.data[o._index(i_o, j_o)], numeric.neg(ty));
            }
        }
    } else {
        var idx_x_outer: usize = 0;
        while (idx_x_outer < x.rows) : (idx_x_outer += 1) {
            var px = x.ptr[idx_x_outer];
            while (px < x.ptr[idx_x_outer + 1]) : (px += 1) {
                const i_stored = if (comptime meta.layoutOf(X) == .col_major) x.idx[px] else idx_x_outer;
                const j_stored = if (comptime meta.layoutOf(X) == .col_major) idx_x_outer else x.idx[px];

                const r_y = if (comptime meta.uploOf(Y) == .upper) (if (i_stored < j_stored) i_stored else j_stored) else (if (i_stored > j_stored) i_stored else j_stored);
                const c_y = if (comptime meta.uploOf(Y) == .upper) (if (i_stored > j_stored) i_stored else j_stored) else (if (i_stored < j_stored) i_stored else j_stored);

                const i_o = if (comptime meta.uploOf(O) == meta.uploOf(X)) i_stored else j_stored;
                const j_o = if (comptime meta.uploOf(O) == meta.uploOf(X)) j_stored else i_stored;

                const tx = if (comptime meta.uploOf(O) == meta.uploOf(X)) x.data[px] else numeric.conj(x.data[px]);
                if (utils.searchSparse(y, r_y, c_y)) |ty| {
                    op_(&o.data[o._index(i_o, j_o)], tx, if (comptime meta.uploOf(O) == meta.uploOf(Y)) ty else numeric.conj(ty));
                } else {
                    numeric.set(&o.data[o._index(i_o, j_o)], tx);
                }
            }
        }

        var idx_y_outer: usize = 0;
        while (idx_y_outer < y.rows) : (idx_y_outer += 1) {
            var py = y.ptr[idx_y_outer];
            while (py < y.ptr[idx_y_outer + 1]) : (py += 1) {
                const i_stored = if (comptime meta.layoutOf(Y) == .col_major) y.idx[py] else idx_y_outer;
                const j_stored = if (comptime meta.layoutOf(Y) == .col_major) idx_y_outer else y.idx[py];

                const r_x = if (comptime meta.uploOf(X) == .upper) (if (i_stored < j_stored) i_stored else j_stored) else (if (i_stored > j_stored) i_stored else j_stored);
                const c_x = if (comptime meta.uploOf(X) == .upper) (if (i_stored > j_stored) i_stored else j_stored) else (if (i_stored < j_stored) i_stored else j_stored);

                const i_o = if (comptime meta.uploOf(O) == meta.uploOf(Y)) i_stored else j_stored;
                const j_o = if (comptime meta.uploOf(O) == meta.uploOf(Y)) j_stored else i_stored;

                if (utils.searchSparse(x, r_x, c_x) == null) {
                    const ty = if (comptime meta.uploOf(O) == meta.uploOf(Y)) y.data[py] else numeric.conj(y.data[py]);
                    if (comptime op_ == numeric.addInto)
                        numeric.set(&o.data[o._index(i_o, j_o)], ty)
                    else
                        numeric.set(&o.data[o._index(i_o, j_o)], numeric.neg(ty));
                }
            }
        }
    }
}
