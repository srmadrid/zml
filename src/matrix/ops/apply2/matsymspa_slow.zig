const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const int = @import("../../../int.zig");
const matrix = @import("../../../matrix.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const O: type = meta.Child(@TypeOf(o));

    const x_nnz = switch (comptime meta.matrixKind(X)) {
        .symmetric => x.nnz,
        .diagonal => o.rows,
        .numeric => 0,
        else => unreachable,
    };

    const y_nnz = switch (comptime meta.matrixKind(Y)) {
        .symmetric => y.nnz,
        .diagonal => o.rows,
        .numeric => 0,
        else => unreachable,
    };

    if (o._dlen < x_nnz + y_nnz or o._ilen < x_nnz + y_nnz)
        return matrix.Error.InsufficientSpace;

    var nnz: usize = 0;
    o.ptr[0] = 0;

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var i: usize = 0;
                while (i <= j) : (i += 1) {
                    const val_x = if (comptime meta.isMatrix(X)) x.get(i, j) catch unreachable else x;
                    const val_y = if (comptime meta.isMatrix(Y)) y.get(i, j) catch unreachable else y;

                    if ((if (comptime !meta.isNumeric(X)) numeric.eq(val_x, 0) else true) and
                        (if (comptime !meta.isNumeric(Y)) numeric.eq(val_y, 0) else true))
                        continue;

                    opInto(&o.data[nnz], val_x, val_y);

                    o.idx[nnz] = i;
                    nnz += 1;
                }
            } else {
                var i: usize = j;
                while (i < o.rows) : (i += 1) {
                    const val_x = if (comptime meta.isMatrix(X)) x.get(i, j) catch unreachable else x;
                    const val_y = if (comptime meta.isMatrix(Y)) y.get(i, j) catch unreachable else y;

                    if ((if (comptime !meta.isNumeric(X)) numeric.eq(val_x, 0) else true) and
                        (if (comptime !meta.isNumeric(Y)) numeric.eq(val_y, 0) else true))
                        continue;

                    opInto(&o.data[nnz], val_x, val_y);

                    o.idx[nnz] = i;
                    nnz += 1;
                }
            }

            o.ptr[j + 1] = nnz;
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var j: usize = i;
                while (j < o.cols) : (j += 1) {
                    const val_x = if (comptime meta.isMatrix(X)) x.get(i, j) catch unreachable else x;
                    const val_y = if (comptime meta.isMatrix(Y)) y.get(i, j) catch unreachable else y;

                    if ((if (comptime !meta.isNumeric(X)) numeric.eq(val_x, 0) else true) and
                        (if (comptime !meta.isNumeric(Y)) numeric.eq(val_y, 0) else true))
                        continue;

                    opInto(&o.data[nnz], val_x, val_y);

                    o.idx[nnz] = j;
                    nnz += 1;
                }
            } else {
                var j: usize = 0;
                while (j <= i) : (j += 1) {
                    const val_x = if (comptime meta.isMatrix(X)) x.get(i, j) catch unreachable else x;
                    const val_y = if (comptime meta.isMatrix(Y)) y.get(i, j) catch unreachable else y;

                    if ((if (comptime !meta.isNumeric(X)) numeric.eq(val_x, 0) else true) and
                        (if (comptime !meta.isNumeric(Y)) numeric.eq(val_y, 0) else true))
                        continue;

                    opInto(&o.data[nnz], val_x, val_y);

                    o.idx[nnz] = j;
                    nnz += 1;
                }
            }

            o.ptr[i + 1] = nnz;
        }
    }

    o.nnz = nnz;
}
