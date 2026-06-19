const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

pub fn matmulInto(o: anytype, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const O: type = meta.Child(@TypeOf(o));

    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    const x_nnz = switch (comptime meta.matrixKind(X)) {
        .general => x.nnz,
        .symmetric, .hermitian => 2 * x.nnz,
        .triangular => if (comptime meta.diagOf(X) == .non_unit) x.nnz else x.nnz + int.min(o.rows, o.cols),
        .diagonal, .permutation => int.min(o.rows, o.cols),
        .builder => unreachable,
        .numeric => 0,
    };

    const y_nnz = switch (comptime meta.matrixKind(Y)) {
        .general => y.nnz,
        .symmetric, .hermitian => 2 * y.nnz,
        .triangular => if (comptime meta.diagOf(Y) == .non_unit) y.nnz else y.nnz + int.min(o.rows, o.cols),
        .diagonal, .permutation => int.min(o.rows, o.cols),
        .numeric => 0,
        .builder => unreachable,
    };

    if (o._dlen < int.min(o.rows * o.cols, x_nnz * y_nnz) or
        o._ilen < int.min(o.rows * o.cols, x_nnz * y_nnz))
        return linalg.Error.InsufficientSpace;

    var nnz: usize = 0;
    o.ptr[0] = 0;

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
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

            o.ptr[j + 1] = nnz;
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
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

            o.ptr[i + 1] = nnz;
        }
    }

    o.nnz = nnz;
}
