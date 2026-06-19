const meta = @import("../../meta.zig");
const int = @import("../../int.zig");
const numeric = @import("../../numeric.zig");
const linalg = @import("../../linalg.zig");

pub fn matmulInto(o: anytype, x: anytype, y: anytype) !void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    if (comptime meta.isVector(X)) {
        const x_len = if (comptime meta.isStaticVector(X)) X.len else x.len;

        const y_nnz = switch (comptime meta.matrixKind(Y)) {
            .general => y.nnz,
            .symmetric, .hermitian => 2 * y.nnz,
            .triangular => if (comptime meta.diagOf(Y).? == .non_unit) y.nnz else y.nnz + int.min(y.rows, y.cols),
            .diagonal, .permutation => int.min(y.rows, y.cols),
            .builder => unreachable,
            .numeric => 0,
        };

        if (o._dlen < int.min(o.len, x.nnz * y_nnz) or
            o._ilen < int.min(o.len, x.nnz * y_nnz))
            return linalg.Error.InsufficientSpace;

        var nnz: usize = 0;

        var j: usize = 0;
        while (j < o.len) : (j += 1) {
            var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

            var k: usize = 0;
            while (k < x_len) : (k += 1) {
                numeric.fmaInto(
                    &sum,
                    x.get(k) catch unreachable,
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

        o.nnz = nnz;
    } else {
        const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

        const x_nnz = switch (comptime meta.matrixKind(X)) {
            .general => x.nnz,
            .symmetric, .hermitian => 2 * x.nnz,
            .triangular => if (comptime meta.diagOf(X).? == .non_unit) x.nnz else x.nnz + int.min(x.rows, x.cols),
            .diagonal, .permutation => int.min(x.rows, x.cols),
            .builder => unreachable,
            .numeric => 0,
        };

        if (o._dlen < int.min(o.len, x_nnz * y.nnz) or
            o._ilen < int.min(o.len, x_nnz * y.nnz))
            return linalg.Error.InsufficientSpace;

        var nnz: usize = 0;

        var i: usize = 0;
        while (i < o.len) : (i += 1) {
            var sum = numeric.zero(meta.Accumulator(meta.Numeric(O)));

            var k: usize = 0;
            while (k < x_cols) : (k += 1) {
                numeric.fmaInto(
                    &sum,
                    x.get(i, k) catch unreachable,
                    y.get(k) catch unreachable,
                    sum,
                );
            }

            if (numeric.eq(sum, 0))
                continue;

            numeric.set(&o.data[nnz], sum);

            o.idx[nnz] = i;
            nnz += 1;
        }

        o.nnz = nnz;
    }
}
