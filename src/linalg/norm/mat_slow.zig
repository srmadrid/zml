const std = @import("std");

const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

const linalg = @import("../../linalg.zig");

pub fn norm(x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) linalg.Norm(@TypeOf(x), order) {
    const X = @TypeOf(x);
    const R = linalg.Norm(X, order);

    const rows = if (comptime meta.isStaticMatrix(X)) X.rows else x.rows;
    const cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    switch (comptime order) {
        .l1 => {
            var max_col_sum = numeric.zero(R);

            for (0..cols) |j| {
                var col_sum = numeric.zero(meta.Accumulator(R));

                for (0..rows) |i| {
                    // sum += abs(x[i, j])
                    numeric.addInto(
                        &col_sum,
                        col_sum,
                        numeric.abs(x.get(i, j) catch unreachable),
                    );
                }

                max_col_sum = numeric.max(max_col_sum, numeric.cast(R, col_sum));
            }

            return max_col_sum;
        },
        .l2 => unreachable,
        .frobenius => {
            var sum = numeric.zero(meta.Accumulator(R));

            for (0..cols) |j| {
                for (0..rows) |i| {
                    // sum += abs(x[i, j])²
                    numeric.addInto(
                        &sum,
                        sum,
                        numeric.abs2(x.get(i, j) catch unreachable),
                    );
                }
            }

            return numeric.cast(R, numeric.sqrt(sum));
        },
        .inf => {
            var max_row_sum = numeric.zero(R);

            for (0..rows) |i| {
                var row_sum = numeric.zero(meta.Accumulator(R));

                for (0..cols) |j| {
                    // sum += abs(x[i, j])
                    numeric.addInto(
                        &row_sum,
                        row_sum,
                        numeric.abs(x.get(i, j) catch unreachable),
                    );
                }

                max_row_sum = numeric.max(max_row_sum, numeric.cast(R, row_sum));
            }

            return max_row_sum;
        },
        .p => unreachable,
    }
}

pub fn normAlloc(allocator: std.mem.Allocator, x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) !linalg.Norm(@TypeOf(x), order) {
    const X = @TypeOf(x);
    const R = linalg.Norm(X, order);

    const rows = if (comptime meta.isStaticMatrix(X)) X.rows else x.rows;
    const cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    switch (comptime order) {
        .l1 => {
            if (comptime meta.layoutOf(X) == .col_major)
                return norm(x, order);

            const col_sums = if (comptime meta.isStaticMatrix(X))
                @as([cols]meta.Accumulator(R), undefined)
            else
                try allocator.alloc(meta.Accumulator(R), cols);
            defer if (comptime !meta.isStaticMatrix(X)) allocator.free(col_sums);
            @memset(col_sums[0..cols], numeric.zero(meta.Accumulator(R)));

            for (0..rows) |i| {
                for (0..cols) |j| {
                    // col_sums[j] += abs(x[i])
                    numeric.addInto(
                        &col_sums[j],
                        col_sums[j],
                        numeric.abs(x.get(i, j) catch unreachable),
                    );
                }
            }

            var max_col_sum = numeric.zero(R);

            for (0..cols) |j| {
                max_col_sum = numeric.max(max_col_sum, numeric.cast(R, col_sums[j]));
            }

            return max_col_sum;
        },
        .l2 => @compileError("zsl.linalg.normAlloc: not implemented yet for \n\tX = " ++ @typeName(X) ++ "\n\torder = " ++ std.fmt.comptimePrint("{any}", .{order}) ++ "\n"),
        .frobenius => return norm(x, order),
        .inf => {
            if (comptime meta.layoutOf(X) == .row_major)
                return norm(x, order);

            const row_sums = if (comptime meta.isStaticMatrix(X))
                @as([rows]meta.Accumulator(R), undefined)
            else
                try allocator.alloc(meta.Accumulator(R), rows);
            defer if (comptime !meta.isStaticMatrix(X)) allocator.free(row_sums);
            @memset(row_sums[0..rows], numeric.zero(meta.Accumulator(R)));

            for (0..cols) |j| {
                for (0..rows) |i| {
                    // row_sums[i] += abs(x[i])
                    numeric.addInto(
                        &row_sums[i],
                        row_sums[i],
                        numeric.abs(x.get(i, j) catch unreachable),
                    );
                }
            }

            var max_row_sum = numeric.zero(R);

            for (0..rows) |i| {
                max_row_sum = numeric.max(max_row_sum, numeric.cast(R, row_sums[i]));
            }

            return max_row_sum;
        },
        .p => unreachable,
    }
}
