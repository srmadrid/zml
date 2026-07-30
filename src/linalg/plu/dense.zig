const std = @import("std");

const int = @import("../../int.zig");
const linalg = @import("../../linalg.zig");
const matrix = @import("../../matrix.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn Dense(N: type, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.linalg.plu.Dense: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        pivots: [*]usize,
        rows: usize,
        cols: usize,
        ld: usize,

        // Type signatures
        pub const is_plu = true;
        pub const is_dense = true;
        pub const storage_layout: matrix.Layout = layout;

        // Numeric type
        pub const Numeric = N;

        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) !linalg.plu.Dense(N, layout) {
            if (rows == 0 or cols == 0)
                return linalg.Error.ZeroDimension;

            const data = try allocator.alloc(N, rows * cols);
            errdefer allocator.free(data);

            return .{
                .data = data.ptr,
                .pivots = (try allocator.alloc(usize, int.min(rows, cols))).ptr,
                .rows = rows,
                .cols = cols,
                .ld = if (comptime layout == .col_major) rows else cols,
            };
        }

        pub fn deinit(self: *linalg.plu.Dense(N, layout), allocator: std.mem.Allocator) void {
            allocator.free(self.data[0 .. self.rows * self.cols]);
            allocator.free(self.pivots[0..int.min(self.rows, self.cols)]);

            self.* = undefined;
        }

        pub fn p(self: linalg.plu.Dense(N, layout), allocator: std.mem.Allocator) !matrix.permutation.Sparse(N, .forward) {
            var _p: matrix.permutation.Sparse(N, .forward) = try .init(allocator, self.rows);
            errdefer _p.deinit(allocator);

            for (0..self.rows) |i| {
                _p.idx[i] = i;
            }

            var i: usize = int.min(self.rows, self.cols);
            while (i > 0) {
                i -= 1;

                const temp = _p.idx[i];
                _p.idx[i] = _p.idx[self.pivots[i]];
                _p.idx[self.pivots[i]] = temp;
            }

            return _p;
        }

        pub fn l(self: linalg.plu.Dense(N, layout)) matrix.triangular.Dense(N, .lower, .unit, layout) {
            return .{
                .data = self.data,
                .rows = self.rows,
                .cols = int.min(self.rows, self.cols),
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn u(self: linalg.plu.Dense(N, layout)) matrix.triangular.Dense(N, .upper, .non_unit, layout) {
            return .{
                .data = self.data,
                .rows = int.min(self.rows, self.cols),
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }
    };
}
