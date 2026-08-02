const std = @import("std");

const int = @import("../../int.zig");
const linalg = @import("../../linalg.zig");
const matrix = @import("../../matrix.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn Dense(N: type, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.linalg.utu.Dense: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        rows: usize,
        cols: usize,
        ld: usize,

        // Type signatures
        pub const is_utu = true;
        pub const is_dense = true;
        pub const storage_layout: matrix.Layout = layout;

        // Numeric type
        pub const Numeric = N;

        pub fn init(allocator: std.mem.Allocator, size: usize) !linalg.utu.Dense(N, layout) {
            if (size == 0)
                return linalg.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(N, size * size)).ptr,
                .rows = size,
                .cols = size,
                .ld = size,
            };
        }

        pub fn deinit(self: *linalg.utu.Dense(N, layout), allocator: std.mem.Allocator) void {
            allocator.free(self.data[0 .. self.rows * self.cols]);

            self.* = undefined;
        }

        pub fn u(self: linalg.utu.Dense(N, layout)) matrix.triangular.Dense(N, .upper, .non_unit, layout) {
            return .{
                .data = self.data,
                .rows = self.rows,
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }
    };
}
