const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");
const matrix = @import("../../matrix.zig");
const linalg = @import("../../linalg.zig");

pub fn matmulInto(o: anytype, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const O: type = meta.Child(@TypeOf(o));

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    const x_nnz = switch (comptime meta.matrixKind(X)) {
        .triangular => if (comptime meta.diagOf(X) == .non_unit) x.nnz else x.nnz + int.min(o.rows, o.cols),
        .diagonal => int.min(o.rows, o.cols),
        .numeric => 0,
        else => unreachable,
    };

    const y_nnz = switch (comptime meta.matrixKind(Y)) {
        .triangular => if (comptime meta.diagOf(Y) == .non_unit) y.nnz else y.nnz + int.min(o.rows, o.cols),
        .diagonal => int.min(o.rows, o.cols),
        .numeric => 0,
        else => unreachable,
    };

    if (o._dlen < int.min(o.rows * o.cols, x_nnz * y_nnz) or
        o._ilen < int.min(o.rows * o.cols, x_nnz * y_nnz))
        return linalg.Error.InsufficientSpace;

    var nnz: usize = 0;
    o.ptr[0] = 0;

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var i: usize = 0;
                while (i < int.min(if (comptime meta.uploOf(O).? == .non_unit) j + 1 else j, o.rows)) : (i += 1) {
                    var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                    var k: usize = 0;
                    while (k < x_cols) : (k += 1) {
                        numeric.fmaInto(
                            &sum,
                            x.get(i, k) catch unreachable,
                            y.get(k, j) catch unreachable,
                            sum,
                        );
                    }

                    if (numeric.eq(sum, 0))
                        continue;

                    numeric.set(&o.data[nnz], sum);

                    o.idx[nnz] = i;
                    nnz += 1;
                }
            } else {
                var i: usize = int.min(if (comptime meta.uploOf(O).? == .non_unit) j else j + 1, o.rows);
                while (i < o.rows) : (i += 1) {
                    var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                    var k: usize = 0;
                    while (k < x_cols) : (k += 1) {
                        numeric.fmaInto(
                            &sum,
                            x.get(i, k) catch unreachable,
                            y.get(k, j) catch unreachable,
                            sum,
                        );
                    }

                    if (numeric.eq(sum, 0))
                        continue;

                    numeric.set(&o.data[nnz], sum);

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
                var j: usize = int.min(if (comptime meta.uploOf(O).? == .non_unit) i else i + 1, o.cols);
                while (j < o.cols) : (j += 1) {
                    var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                    var k: usize = 0;
                    while (k < x_cols) : (k += 1) {
                        numeric.fmaInto(
                            &sum,
                            x.get(i, k) catch unreachable,
                            y.get(k, j) catch unreachable,
                            sum,
                        );
                    }

                    if (numeric.eq(sum, 0))
                        continue;

                    numeric.set(&o.data[nnz], sum);

                    o.idx[nnz] = j;
                    nnz += 1;
                }
            } else {
                var j: usize = 0;
                while (j < int.min(if (comptime meta.uploOf(O).? == .non_unit) i + 1 else i, o.cols)) : (j += 1) {
                    var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

                    var k: usize = 0;
                    while (k < x_cols) : (k += 1) {
                        numeric.fmaInto(
                            &sum,
                            x.get(i, k) catch unreachable,
                            y.get(k, j) catch unreachable,
                            sum,
                        );
                    }

                    if (numeric.eq(sum, 0))
                        continue;

                    numeric.set(&o.data[nnz], sum);

                    o.idx[nnz] = j;
                    nnz += 1;
                }
            }

            o.ptr[i + 1] = nnz;
        }
    }

    o.nnz = nnz;
}
