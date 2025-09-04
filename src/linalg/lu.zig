const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const Numeric = types.Numeric;

const matrix = @import("../matrix.zig");
const General = matrix.General;
const Triangular = matrix.Triangular;
const Permutation = matrix.Permutation;

pub fn LU(T: type, order: Order) type {
    return struct {
        l: Triangular(T, .lower, .unit, order),
        u: Triangular(T, .upper, .non_unit, order),

        pub fn init(lu: anytype, m: u32, n: u32) LU(Numeric(@TypeOf(lu)), order) {
            return .{
                .l = .{
                    .data = lu,
                    .rows = m,
                    .cols = n,
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = true },
                },
                .u = .{
                    .data = lu,
                    .rows = m,
                    .cols = n,
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = false },
                },
            };
        }

        pub fn deinit(self: *LU(T, order), allocator: std.mem.Allocator) void {
            self.l.deinit(allocator);

            self.* = undefined;
        }
    };
}

pub fn PLU(T: type, order: Order) type {
    return struct {
        p: Permutation(T),
        l: Triangular(T, .lower, .unit, order),
        u: Triangular(T, .upper, .non_unit, order),

        pub fn init(p: [*]u32, lu: anytype, m: u32, n: u32) PLU(Numeric(@TypeOf(lu)), order) {
            return .{
                .p = .{
                    .data = p,
                    .size = m,
                    .flags = .{ .owns_data = true },
                },
                .l = .{
                    .data = lu,
                    .rows = m,
                    .cols = n,
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = true },
                },
                .u = .{
                    .data = lu,
                    .rows = m,
                    .cols = n,
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = false },
                },
            };
        }

        pub fn deinit(self: *PLU(T, order), allocator: std.mem.Allocator) void {
            self.p.deinit(allocator);
            self.l.deinit(allocator);

            self.* = undefined;
        }
    };
}

pub fn PLUQ(T: type, order: Order) type {
    return struct {
        p: Permutation(T),
        l: Triangular(T, .lower, .unit, order),
        u: Triangular(T, .upper, .non_unit, order),
        q: Permutation(T),

        pub fn init(p: [*]u32, lu: anytype, q: [*]u32, m: u32, n: u32) PLUQ(Numeric(@TypeOf(lu)), order) {
            return .{
                .p = .{
                    .data = p,
                    .size = m,
                    .flags = .{ .owns_data = true },
                },
                .l = .{
                    .data = lu,
                    .rows = m,
                    .cols = n,
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = true },
                },
                .u = .{
                    .data = lu,
                    .rows = m,
                    .cols = n,
                    .ld = if (order == .col_major) m else n,
                    .flags = .{ .owns_data = false },
                },
                .q = .{
                    .data = q,
                    .size = n,
                    .flags = .{ .owns_data = true },
                },
            };
        }

        pub fn deinit(self: *PLUQ(T, order), allocator: std.mem.Allocator) void {
            self.p.deinit(allocator);
            self.l.deinit(allocator);
            self.q.deinit(allocator);

            self.* = undefined;
        }
    };
}

pub fn pluq() void {
    // Make sure to convert jpiv from column permutation to row permutation
    // (Permutation matrix type assumes row permutations and, in case whe
    // multiplicate from the right, like A * P, matmul already takes care of
    // converting it to column permutation internally, i.e., doing the inverse
    // permutation). Also from 1 to 0-based indexing.
}
