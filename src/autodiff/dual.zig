const std = @import("std");

const types = @import("../types.zig");
const numeric = @import("../numeric.zig");
const constants = @import("../constants.zig");

pub fn isDual(comptime T: type) bool {
    return @hasDecl(T, "is_dual") and T.is_dual;
}

pub fn Dual(comptime T: type) type {
    if (comptime !types.isNumeric(T))
        @compileError("zml.autodiff.Dual: T must be a numeric type, got \n\tT: " ++ @typeName(T) ++ "\n");

    return struct {
        val: T,
        eps: T,

        // Type signature
        pub const is_custom = true;
        pub const is_numeric = true;
        pub const is_dual = true;
        pub const is_allocated = types.isAllocated(T);
        pub const is_complex_type = types.isComplexType(T);

        /// Scalar type
        pub const Scalar = T;

        pub const empty: Dual(T) = .{
            .val = undefined,
            .eps = undefined,
        };

        pub fn deinit(self: *Dual(T), allocator: std.mem.Allocator) void {
            if (comptime types.isAllocated(T)) {
                numeric.deinit(&self.val, .{ .allocator = allocator });
                numeric.deinit(&self.eps, .{ .allocator = allocator });
            }
        }

        // Constants
        pub const zero = if (types.isAllocated(T)) _zeroAllocated else _zero;

        fn _zero() Dual(T) {
            return .{
                .val = constants.zero(T, .{}) catch unreachable,
                .eps = constants.zero(T, .{}) catch unreachable,
            };
        }

        fn _zeroAllocated(allocator: ?std.mem.Allocator) !Dual(T) {
            var val = try constants.zero(T, .{ .allocator = allocator });
            errdefer numeric.deinit(&val, .{ .allocator = allocator });

            return .{
                .val = val,
                .eps = try constants.zero(T, .{ .allocator = allocator }),
            };
        }

        // Basic operations
        pub const Abs = _Abs;
        pub const abs = if (types.isAllocated(T)) _absAllocated else _abs;

        pub const Add = _Add;
        pub const add = if (types.isAllocated(T)) _addAllocated else _add;

        pub fn fromFloat(x: anytype) Dual(T) {
            return .{
                .val = types.scast(T, x),
                .eps = constants.zero(T, .{}) catch unreachable,
            };
        }

        pub fn toCfloat(self: Dual(T), comptime Cfloat: type) Cfloat {
            return types.scast(Cfloat, self.val);
        }
    };
}

fn _Abs(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !isDual(X))
        @compileError("zml.Dual(T).abs: X must be a dual type, got\n\tX: " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Abs(types.Scalar(X)));
}

fn _abs(x: anytype) _Abs(@TypeOf(x)) {
    return if (comptime types.isComplexType(types.Scalar(@TypeOf(x))))
        .{
            .val = numeric.abs(x.val, .{}) catch unreachable,
            .eps = numeric.div(
                numeric.re(
                    numeric.mul(
                        x.val,
                        numeric.conj(x.eps, .{}) catch unreachable,
                        .{},
                    ) catch unreachable,
                    .{},
                ) catch unreachable,
                numeric.abs(x.val, .{}) catch unreachable,
                .{},
            ) catch unreachable,
        }
    else
        .{
            .val = numeric.abs(x.val, .{}) catch unreachable,
            .eps = numeric.mul(
                numeric.sgn(x.val, .{}) catch unreachable,
                x.eps,
                .{},
            ) catch unreachable,
        };
}

fn _absAllocated(allocator: std.mem.Allocator, x: anytype) !_Abs(@TypeOf(x)) {
    var val = try numeric.abs(allocator, x.val, .{ .allocator = allocator });
    errdefer numeric.deinit(&val, .{ .allocator = allocator });

    return .{
        .val = val,
        .eps = try numeric.mul(
            numeric.sgn(x.val, .{}) catch unreachable,
            x.eps,
            .{ .allocator = allocator },
        ),
    };
}

fn _Add(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zml.Dual(T).add: at least one of x or y must be a dual, the other must be a numeric, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) types.Scalar(X) else X;
    const SY: type = if (isDual(Y)) types.Scalar(Y) else Y;

    return Dual(numeric.Add(SX, SY));
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
