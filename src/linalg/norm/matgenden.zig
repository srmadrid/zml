const std = @import("std");

const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

const linalg = @import("../../linalg.zig");

pub fn norm(x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) linalg.Norm(@TypeOf(x), order) {
    const X = @TypeOf(x);
    const R = linalg.Norm(X, order);

    switch (comptime order) {
        .l1 => {
            var max_col_sum = numeric.zero(R);

            for (0..x.cols) |j| {
                var col_sum = numeric.zero(meta.Accumulator(R));

                for (0..x.rows) |i| {
                    // sum += abs(x[i, j])
                    numeric.addInto(
                        &col_sum,
                        col_sum,
                        numeric.abs(x.data[x._index(i, j)]),
                    );
                }

                max_col_sum = numeric.max(max_col_sum, numeric.cast(R, col_sum));
            }

            return max_col_sum;
        },
        .l2 => unreachable,
        .frobenius => {
            var sum = numeric.zero(meta.Accumulator(R));

            if (comptime meta.layoutOf(X) == .col_major) {
                for (0..x.cols) |j| {
                    for (0..x.rows) |i| {
                        // sum += abs(x[i, j])²
                        numeric.addInto(
                            &sum,
                            sum,
                            numeric.abs2(x.data[i + j * x.ld]),
                        );
                    }
                }
            } else {
                for (0..x.rows) |i| {
                    for (0..x.cols) |j| {
                        // sum += abs(x[i, j])²
                        numeric.addInto(
                            &sum,
                            sum,
                            numeric.abs2(x.data[i * x.ld + j]),
                        );
                    }
                }
            }

            return numeric.cast(R, numeric.sqrt(sum));
        },
        .inf => {
            var max_row_sum = numeric.zero(R);

            for (0..x.rows) |i| {
                var row_sum = numeric.zero(meta.Accumulator(R));

                for (0..x.cols) |j| {
                    // sum += abs(x[i, j])
                    numeric.addInto(
                        &row_sum,
                        row_sum,
                        numeric.abs(x.data[x._index(i, j)]),
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

    switch (comptime order) {
        .l1 => {
            if (comptime meta.layoutOf(X) == .col_major)
                return norm(x, order);

            const col_sums = try allocator.alloc(meta.Accumulator(R), x.cols);
            defer allocator.free(col_sums);
            @memset(col_sums, numeric.zero(meta.Accumulator(R)));

            for (0..x.rows) |i| {
                for (0..x.cols) |j| {
                    // col_sums[j] += abs(x[i])
                    numeric.addInto(
                        &col_sums[j],
                        col_sums[j],
                        numeric.abs(x.data[i * x.ld + j]),
                    );
                }
            }

            var max_col_sum = numeric.zero(R);

            for (0..x.cols) |j| {
                max_col_sum = numeric.max(max_col_sum, numeric.cast(R, col_sums[j]));
            }

            return max_col_sum;
        },
        .l2 => @compileError("zsl.linalg.normAlloc: not implemented yet for \n\tX = " ++ @typeName(X) ++ "\n\torder = " ++ std.fmt.comptimePrint("{any}", .{order}) ++ "\n"),
        .frobenius => return norm(x, order),
        .inf => {
            if (comptime meta.layoutOf(X) == .row_major)
                return norm(x, order);

            const row_sums = try allocator.alloc(meta.Accumulator(R), x.rows);
            defer allocator.free(row_sums);
            @memset(row_sums, numeric.zero(meta.Accumulator(R)));

            for (0..x.cols) |j| {
                for (0..x.rows) |i| {
                    // row_sums[i] += abs(x[i])
                    numeric.addInto(
                        &row_sums[i],
                        row_sums[i],
                        numeric.abs(x.data[i + j * x.ld]),
                    );
                }
            }

            var max_row_sum = numeric.zero(R);

            for (0..x.rows) |i| {
                max_row_sum = numeric.max(max_row_sum, numeric.cast(R, row_sums[i]));
            }

            return max_row_sum;
        },
        .p => unreachable,
    }
}
