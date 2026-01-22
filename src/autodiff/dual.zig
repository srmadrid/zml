const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");

pub fn isDual(comptime T: type) bool {
    return @hasDecl(T, "is_dual");
}

pub fn Dual(comptime T: type) type {
    if (comptime !types.isNumeric(T))
        @compileError("zml.autodiff: T must be a numeric type, got \n\tT: " ++ @typeName(T) ++ "\n");

    return struct {
        val: T,
        eps: T,

        // Type signature
        pub const is_numeric = {};
        pub const is_dual = {};
        pub const is_allocated = types.isAllocated(T);

        /// Scalar type
        pub const Scalar = T;

        pub const empty: Dual(T) = .{
            .val = undefined,
            .eps = undefined,
        };

        pub const zero = if (types.isAllocated(T)) _zeroAllocated else _zero;

        pub fn _zero() Dual(T) {
            return .{
                .val = constants.zero(T, .{}) catch unreachable,
                .eps = constants.zero(T, .{}) catch unreachable,
            };
        }

        fn _zeroAllocated(allocator: ?std.mem.Allocator) !Dual(T) {
            return .{
                .val = try constants.zero(T, .{ .allocator = allocator }),
                .eps = try constants.zero(T, .{ .allocator = allocator }),
            };
        }

        pub fn deinit(self: *Dual(T), allocator: std.mem.Allocator) void {
            if (comptime types.isAllocated(T)) {
                ops.deinit(&self.val, .{ .allocator = allocator });
                ops.deinit(&self.eps, .{ .allocator = allocator });
            }
        }

        pub const Add = _Add;
        pub const add = if (types.isAllocated(T)) _addAllocated else _add;
    };
}

fn _Add(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zml.Dual(T).add: at least one of x or y must be a dual, the other must be a numeric, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) types.Scalar(X) else X;
    const SY: type = if (isDual(Y)) types.Scalar(Y) else Y;

    return Dual(ops.Add(SX, SY));
}

fn _add(x: anytype, y: anytype) _Add(@TypeOf(x), @TypeOf(y)) {
    // if (comptime isDual(@TypeOf(x))) {
    //     if (comptime isDual(@TypeOf(y))) {
    //         return .{
    //             .val = ops.add(x.val, y.val, .{}) catch unreachable,
    //             .eps = ops.add(x.eps, y.eps, .{}) catch unreachable,
    //         };
    //     } else {
    //         return .{
    //             .val = ops.add(x.val, y, .{}) catch unreachable,
    //             .eps = x.eps,
    //         };
    //     }
    // } else {
    //     return .{
    //         .val = ops.add(x, y.val, .{}) catch unreachable,
    //         .eps = y.eps,
    //     };
    // }
    return .empty;
}

fn _addAllocated(allocator: std.mem.Allocator, x: anytype, y: anytype) !_Add(@TypeOf(x), @TypeOf(y)) {
    // const X: type = @TypeOf(x);
    // const Y: type = @TypeOf(y);
    // const R: type = _Add(X, Y);

    // if (comptime isDual(X)) {
    //     if (comptime isDual(Y)) {
    //         var val = try ops.add(allocator, x.val, y.val, .{ .allocator = allocator });
    //         errdefer ops.deinit(&val, .{ .allocator = allocator });

    //         return .{
    //             .val = val,
    //             .eps = try ops.add(allocator, x.eps, y.eps, .{ .allocator = allocator }),
    //         };
    //     } else {
    //         var val = try ops.add(allocator, x.val, y, .{ .allocator = allocator });
    //         errdefer ops.deinit(&val, .{ .allocator = allocator });

    //         return .{
    //             .val = val,
    //             .eps = try types.cast(types.Scalar(R), x.eps, .{ .allocator = allocator }),
    //         };
    //     }
    // } else {
    //     var val = try ops.add(allocator, x, y.val, .{ .allocator = allocator });
    //     errdefer ops.deinit(&val, .{ .allocator = allocator });

    //     return .{
    //         .val = val,
    //         .eps = try types.cast(types.Scalar(R), y.eps, .{ .allocator = allocator }),
    //     };
    // }
    _ = allocator;
    return .empty;
}
