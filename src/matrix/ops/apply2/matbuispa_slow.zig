const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    var nnz: usize = 0;

    const is_col_major = comptime blk: {
        if (meta.isMatrix(X) and meta.layoutOf(X) != null)
            break :blk meta.layoutOf(X).? == .col_major
        else if (meta.isMatrix(Y) and meta.layoutOf(Y) != null)
            break :blk meta.layoutOf(Y).? == .col_major;

        break :blk true;
    };

    if (comptime is_col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                const val_x = if (comptime meta.isMatrix(X)) x.get(i, j) catch unreachable else x;
                const val_y = if (comptime meta.isMatrix(Y)) y.get(i, j) catch unreachable else y;

                if ((if (comptime !meta.isNumeric(X)) numeric.eq(val_x, 0) else true) and
                    (if (comptime !meta.isNumeric(Y)) numeric.eq(val_y, 0) else true))
                    continue;

                opInto(&o.data[nnz], val_x, val_y);

                o.ridx[nnz] = i;
                o.cidx[nnz] = j;
                nnz += 1;
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                const val_x = if (comptime meta.isMatrix(X)) x.get(i, j) catch unreachable else x;
                const val_y = if (comptime meta.isMatrix(Y)) y.get(i, j) catch unreachable else y;

                if ((if (comptime !meta.isNumeric(X)) numeric.eq(val_x, 0) else true) and
                    (if (comptime !meta.isNumeric(Y)) numeric.eq(val_y, 0) else true))
                    continue;

                opInto(&o.data[nnz], val_x, val_y);

                o.ridx[nnz] = i;
                o.cidx[nnz] = j;
                nnz += 1;
            }
        }
    }

    o.nnz = nnz;
}
