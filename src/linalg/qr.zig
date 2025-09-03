const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const Numeric = types.Numeric;
const Coerce = types.Coerce;

const matrix = @import("../matrix.zig");
const General = matrix.General;
const Triangular = matrix.Triangular;
const Permutation = matrix.Permutation;

pub fn QR(T: type, order: Order) type {
    return struct {
        q: General(T, order),
        r: Triangular(T, .upper, .non_unit, order),

        pub fn init(q: anytype, r: anytype, m: u32, n: u32, k: u32) QR(Numeric(Coerce(@TypeOf(q), @TypeOf(r))), order) {
            return .{
                .q = .{
                    .data = q,
                    .rows = m,
                    .cols = k,
                    .ld = if (order == .col_major) m else k,
                    .flags = .{ .owns_data = true },
                },
                .r = .{
                    .data = r,
                    .rows = k,
                    .cols = n,
                    .ld = if (order == .col_major) k else n,
                    .flags = .{ .owns_data = true },
                },
            };
        }

        pub fn deinit(self: *QR(T, order), allocator: std.mem.Allocator) void {
            self.q.deinit(allocator);
            self.r.deinit(allocator);

            self.* = undefined;
        }
    };
}

pub fn QRP(T: type, order: Order) type {
    return struct {
        q: General(T, order),
        r: Triangular(T, .upper, .non_unit, order),
        p: Permutation(T),

        pub fn init(q: anytype, r: anytype, p: anytype, m: u32, n: u32, k: u32) QRP(Numeric(Coerce(@TypeOf(q), @TypeOf(r))), order) {
            return .{
                .q = .{
                    .data = q,
                    .rows = m,
                    .cols = k,
                    .ld = if (order == .col_major) m else k,
                    .flags = .{ .owns_data = true },
                },
                .r = .{
                    .data = r,
                    .rows = k,
                    .cols = n,
                    .ld = if (order == .col_major) k else n,
                    .flags = .{ .owns_data = true },
                },
                .p = .{
                    .data = p,
                    .size = n,
                    .flags = .{ .owns_data = true },
                },
            };
        }

        pub fn deinit(self: *QRP(T, order), allocator: std.mem.Allocator) void {
            self.q.deinit(allocator);
            self.r.deinit(allocator);
            self.p.deinit(allocator);

            self.* = undefined;
        }
    };
}
