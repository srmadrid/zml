const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const scast = types.scast;
const cast = types.cast;
const needsAllocator = types.needsAllocator;
const validateContext = types.validateContext;

const ops = @import("../ops.zig");
const int = @import("../int.zig");

const array = @import("../array.zig");
const max_dimensions = array.max_dimensions;
const Order = array.Order;
const Flags = array.Flags;
const Range = array.Range;

const strided = @import("strided.zig");
const Strided = strided.Strided;

const general = @import("dense/general.zig");
const symmetric = @import("dense/symmetric.zig");
const hermitian = @import("dense/hermitian.zig");
const triangular = @import("dense/triangular.zig");
const diagonal = @import("dense/diagonal.zig");
const banded = @import("dense/banded.zig");
const tridiagonal = @import("dense/tridiagonal.zig");

pub const Kind = union(enum) {
    general: General,
    symmetric: Symmetric,
    hermitian: Hermitian,
    triangular: Triangular,
    diagonal: Diagonal,
    banded: Banded,
    tridiagonal: Tridiagonal,

    pub const General = struct {};

    pub const Symmetric = packed struct {
        upper: bool = true,
    };

    pub const Hermitian = packed struct {
        upper: bool = true,
    };

    pub const Triangular = packed struct {
        upper: bool = true,
        unit: bool = false,
    };

    pub const Diagonal = struct {};

    pub const Banded = packed struct {
        lower: u16 = 0,
        upper: u16 = 0,
    };

    pub const Tridiagonal = struct {
        superdiagonal_offset: u32 = 0,
    };
};

pub fn Dense(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("Dense requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: []T, // Make [*]T?
        ndim: u32,
        shape: [max_dimensions]u32,
        strides: [max_dimensions]u32,
        size: u32,
        base: ?*anyopaque,
        flags: Flags = .{},
        kind: Kind = .{ .general = .{} },

        pub const empty: Dense(T) = .{
            .data = &.{},
            .ndim = 0,
            .shape = .{0} ** max_dimensions,
            .base = null,
            .flags = .{},
            .kind = .{},
        };

        pub fn init(
            allocator: std.mem.Allocator,
            shape: []const u32,
            options: struct {
                order: Order = .col_major,
                kind: Kind = .{ .general = .{} },
            },
        ) !Dense(T) {
            if (options.kind != .general) {
                if (shape.len > 2)
                    return array.Error.TooManyDimensions;
            } else if (shape.len > max_dimensions) {
                return array.Error.TooManyDimensions;
            }

            if (comptime !types.isComplex(T)) {
                if (options.kind == .hermitian)
                    return array.Error.InvalidFlags;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return array.Error.ZeroDimension;
                }
            }

            return switch (options.kind) {
                .general => general.init(T, allocator, shape, options.order),
                .symmetric => symmetric.init(T, allocator, shape, options.order, options.kind.symmetric.upper),
                .hermitian => hermitian.init(T, allocator, shape, options.order, options.kind.hermitian.upper),
                .triangular => triangular.init(T, allocator, shape, options.order, options.kind.triangular.upper, options.kind.triangular.unit),
                .diagonal => diagonal.init(T, allocator, shape),
                .banded => banded.init(T, allocator, shape, options.order, options.kind.banded.lower, options.kind.banded.upper),
                .tridiagonal => tridiagonal.init(T, allocator, shape),
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            shape: []const u32,
            value: anytype,
            options: struct {
                order: Order = .col_major,
                kind: Kind = .{ .general = .{} },
            },
            ctx: anytype,
        ) !Dense(T) {
            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            if (options.kind != .general) {
                if (shape.len > 2)
                    return array.Error.TooManyDimensions;
            } else if (shape.len > max_dimensions) {
                return array.Error.TooManyDimensions;
            }

            if (comptime !types.isComplex(T)) {
                if (options.kind == .hermitian)
                    return array.Error.InvalidFlags;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return array.Error.ZeroDimension;
                }
            }

            return switch (options.kind) {
                .general => general.full(T, allocator, shape, value, options.order, ctx),
                .symmetric => symmetric.full(T, allocator, shape, value, options.order, options.kind.symmetric.upper, ctx),
                .hermitian => hermitian.full(T, allocator, shape, value, options.order, options.kind.hermitian.upper, ctx),
                .triangular => triangular.full(T, allocator, shape, value, options.order, options.kind.triangular.upper, options.kind.triangular.unit, ctx),
                .diagonal => diagonal.full(T, allocator, shape, value, ctx),
                .banded => banded.full(T, allocator, shape, value, options.order, options.kind.banded.lower, options.kind.banded.upper, ctx),
                .tridiagonal => tridiagonal.full(T, allocator, shape, value, options.order, ctx),
            };
        }

        pub fn arange(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            step: anytype,
            ctx: anytype,
        ) !Dense(T) {
            comptime if (types.isComplex(T))
                @compileError("array.arange does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.arange: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.arange: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.arange: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            comptime if (!types.isNumeric(@TypeOf(step)))
                @compileError("array.arange: step must be numeric, got " ++ @typeName(@TypeOf(step)));

            comptime if (types.isComplex(@TypeOf(step)))
                @compileError("array.arange: step cannot be complex, got " ++ @typeName(@TypeOf(step)));

            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            return general.arange(T, allocator, start, stop, step, ctx);
        }

        pub fn linspace(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            options: struct {
                num: u32 = 50,
                endpoint: bool = true,
                retstep: ?*T = null,
            },
            ctx: anytype,
        ) !Dense(T) {
            comptime if (types.isComplex(T))
                @compileError("array.linspace does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.linspace: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.linspace: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.linspace: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            if (options.num == 0)
                return array.Error.ZeroDimension;

            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            return general.linspace(T, allocator, start, stop, options.num, options.endpoint, options.retstep, ctx);
        }

        pub fn logspace(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            base: anytype,
            options: struct {
                num: u32 = 50,
                endpoint: bool = true,
            },
            ctx: anytype,
        ) !Dense(T) {
            comptime if (types.isComplex(T))
                @compileError("array.logspace does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.logspace: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.logspace: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.logspace: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            if (options.num == 0)
                return array.Error.ZeroDimension;

            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            return general.logspace(T, allocator, start, stop, base, options.num, options.endpoint, ctx);
        }

        /// Cleans up the array by deinitializing its elements. If the array holds
        /// a fixed precision type this is does not do anything.
        pub fn cleanup(self: *Dense(T), ctx: anytype) void {
            _ = self;
            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            if (comptime types.isArbitraryPrecision(T)) {
                @compileError("Cleanup not implemented yet");
            }
        }

        /// Deinitializes the array, freeing its data if it owns it.
        ///
        /// If the array holds an arbitrary precision type, it will not free the
        /// elements; call `cleanup` first to do that.
        pub fn deinit(self: *Dense(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data);
            }

            self.* = undefined;
        }

        pub fn set(self: *const Dense(T), position: []const u32, value: T) !void {
            switch (self.kind) {
                .general => try general.set(T, self, position, value),
                .symmetric => try symmetric.set(T, self, position, value),
                .hermitian => try hermitian.set(T, self, position, value),
                .triangular => try triangular.set(T, self, position, value),
                .diagonal => try diagonal.set(T, self, position, value),
                .banded => try banded.set(T, self, position, value),
                .tridiagonal => try tridiagonal.set(T, self, position, value),
            }
        }

        pub fn get(self: *const Dense(T), position: []const u32) !T {
            return switch (self.kind) {
                .general => general.get(T, self, position),
                .symmetric => symmetric.get(T, self, position),
                .hermitian => hermitian.get(T, self, position),
                .triangular => triangular.get(T, self, position),
                .diagonal => diagonal.get(T, self, position),
                .banded => banded.get(T, self, position),
                .tridiagonal => tridiagonal.get(T, self, position),
            };
        }

        pub fn reshape(self: *Dense(T), shape: []const u32) !Strided(T) {
            if (shape.len > max_dimensions) {
                return array.Error.TooManyDimensions;
            }

            if (shape.len == 0) {
                return array.Error.ZeroDimension;
            }

            return switch (self.kind) {
                .general => general.reshape(T, self, shape),
                else => array.Error.InvalidKind, // Reshape is only supported for general arrays.
            };
        }

        pub fn ravel(self: *Dense(T)) !Strided(T) {
            if (self.ndim == 0) {
                return array.Error.ZeroDimension;
            }

            return switch (self.kind) {
                .general => general.ravel(T, self),
                else => array.Error.InvalidKind, // Ravel is only supported for general arrays.
            };
        }

        pub fn transpose(self: *Dense(T), axes: ?[]const u32) !Dense(T) {
            const axes_: []const u32 =
                axes orelse
                array.trivialReversePermutation(self.ndim)[0..self.ndim];

            if (axes_.len == 0) {
                return array.Error.ZeroDimension;
            }

            if (axes_.len != self.ndim) {
                return array.Error.DimensionMismatch;
            }

            if (!array.isPermutation(self.ndim, axes_)) {
                return array.Error.InvalidAxes; // axes must be a valid permutation of [0, ..., ndim - 1]
            }

            switch (self.kind) {
                .general => return general.transpose(T, self, axes_),
                .symmetric => return symmetric.transpose(T, self, axes_),
                .hermitian => return hermitian.transpose(T, self, axes_),
                .triangular => return triangular.transpose(T, self, axes_),
                .diagonal => return diagonal.transpose(T, self, axes_),
                .banded => return banded.transpose(T, self, axes_),
                .tridiagonal => return tridiagonal.transpose(T, self, axes_),
            }
        }

        pub fn slice(self: *const Dense(T), ranges: []const Range) !Strided(T) {
            if (ranges.len == 0 or ranges.len > self.ndim) {
                return error.DimensionMismatch;
            }

            switch (self.kind) {
                .general => return general.slice(T, self, ranges),
                .symmetric => return symmetric.slice(T, self, ranges),
                .hermitian => return hermitian.slice(T, self, ranges),
                .triangular => return triangular.slice(T, self, ranges),
                .diagonal => return diagonal.slice(T, self, ranges),
                .banded => return banded.slice(T, self, ranges),
                .tridiagonal => return tridiagonal.slice(T, self, ranges),
            }
        }

        pub fn broadcast(self: *const Dense(T), shape: []const u32) !Strided(T) {
            if (shape.len > max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len < self.ndim) {
                return array.Error.TooLittleDimensions;
            }

            switch (self.kind) {
                .general => return general.broadcast(T, self, shape),
                else => return array.Error.InvalidKind, // Broadcasting is not supported for symmetric, hermitian, triangular, diagonal, banded, or tridiagonal arrays.
            }
        }
    };
}

pub inline fn apply1(
    comptime T: type,
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    order: array.Order,
    ctx: anytype,
) !Dense(ReturnType1(op, T)) {
    var result: Dense(ReturnType1(op, T)) = try .init(ReturnType1(op, T), allocator, x.shape[0..x.ndim], order);
    errdefer result.deinit(allocator);

    var j: u32 = 0;
    if (result.flags.order == x.flags.order) {
        // Trivial loop
        //errdefer cleanup(ReturnType1(op, T), allocator, result.data[0..j]);

        const opinfo = @typeInfo(@TypeOf(op));
        for (0..result.size) |i| {
            if (comptime opinfo.@"fn".params.len == 1) {
                result.data[i] = op(x.data[i]);
            } else if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = try op(x.data[i], ctx);
            }

            j += 1;
        }
    } else {
        // Different order, but same shape
        const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
        errdefer strided.cleanup(ReturnType1(op, T), &result, j, iteration_order, ctx);
        const axis: u32 = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
        var iterr: array.Iterator(ReturnType1(op, T)) = .init(&result);
        var iterx: array.Iterator(T) = .init(&x);
        const opinfo = @typeInfo(@TypeOf(op));
        for (0..result.size) |_| {
            if (comptime opinfo.@"fn".params.len == 1) {
                result.data[iterr.index] = op(x.data[iterx.index]);
            } else if (comptime opinfo.@"fn".params.len == 2) {
                result.data[iterr.index] = try op(x.data[iterx.index], ctx);
            }

            _ = iterr.nextAO(axis, iteration_order);
            _ = iterx.nextAO(axis, iteration_order);
        }
    }

    return result;
}

pub fn apply1_(
    comptime O: type,
    o: anytype,
    comptime X: type,
    x: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    if (comptime !types.isDense(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        // Only run the function once if x is a scalar
        const opinfo = @typeInfo(@TypeOf(op_));
        if (comptime opinfo.@"fn".params.len == 2) {
            op_(&o.data[0], x);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            try op_(&o.data[0], x, ctx);
        }

        for (1..o.size) |i| {
            try ops.set(&o.data[i], x, ctx);
        }

        return;
    }

    var xx: Dense(X) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
        if (o.flags.order == x.flags.order) {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (comptime opinfo.@"fn".params.len == 2) {
                    op_(&o.data[i], x.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try op_(&o.data[i], x.data[i], ctx);
                }
            }

            return;
        } else {
            // Different order, but same shape
            xx = x;
        }
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });
        if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
            return array.Error.NotBroadcastable;
        }

        //xx = try broadcast(X, &x, bct.shape[0..bct.ndim]);
    }

    const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
    const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
    var itero: array.Iterator(O) = .init(o);
    var iterx: array.Iterator(X) = .init(&xx);
    const opinfo = @typeInfo(@TypeOf(op_));
    for (0..o.size) |_| {
        if (comptime opinfo.@"fn".params.len == 2) {
            op_(&o.data[itero.index], xx.data[iterx.index]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            try op_(&o.data[itero.index], xx.data[iterx.index], ctx);
        }

        _ = itero.nextAO(axis, iteration_order);
        _ = iterx.nextAO(axis, iteration_order);
    }

    return;
}

// Apply2 outline:
// 1. check for equality of shapes
//     - if equal, check for order (either equal loop, or one loops backwards)
// 2. if not equal, check for broadcasting
//     - if broadcasting is possible, apply broadcasting (efficient inner loop, like in innerLoop2RLSS, maybe moove all these inner loop functions to a separate file?, instead of DD being here in dense, and DS and SS in strided)
//     - if broadcasting is not possible, return error
pub fn apply2(
    allocator: std.mem.Allocator,
    comptime X: type,
    x: anytype,
    comptime Y: type,
    y: anytype,
    comptime op: anytype,
    order: array.Order,
    ctx: anytype,
) !Dense(ReturnType2(op, X, Y)) {
    if (comptime !types.isDense(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        var result: Dense(ReturnType2(op, X, Y)) = try .init(ReturnType2(op, X, Y), allocator, y.shape[0..y.ndim], order);
        errdefer result.deinit(allocator);

        var j: u32 = 0;
        if (result.flags.order == y.flags.order) {
            // Trivial loop
            //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

            const opinfo = @typeInfo(@TypeOf(op));
            for (0..result.size) |i| {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x, y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x, y.data[i], ctx);
                }

                j += 1;
            }
        } else {
            const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
            const axis: u32 = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
            errdefer strided.cleanup(ReturnType2(op, X, Y), &result, j, iteration_order, ctx);
            var iterr: array.Iterator(ReturnType2(op, X, Y)) = .init(&result);
            var itery: array.Iterator(Y) = .init(&y);
            const opinfo = @typeInfo(@TypeOf(op));
            for (0..result.size) |_| {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[iterr.index] = op(x, y.data[itery.index]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[iterr.index] = try op(x, y.data[itery.index], ctx);
                }

                j += 1;
                _ = iterr.nextAO(axis, iteration_order);
                _ = itery.nextAO(axis, iteration_order);
            }
        }

        return result;
    } else if (comptime !types.isDense(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        var result: Dense(ReturnType2(op, X, Y)) = try .init(ReturnType2(op, X, Y), allocator, x.shape[0..x.ndim], order);
        errdefer result.deinit(allocator);

        var j: u32 = 0;
        if (result.flags.order == x.flags.order) {
            //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

            const opinfo = @typeInfo(@TypeOf(op));
            for (0..result.size) |i| {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], y, ctx);
                }

                j += 1;
            }
        } else {
            const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
            const axis: u32 = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
            errdefer strided.cleanup(ReturnType2(op, X, Y), &result, j, iteration_order, ctx);
            var iterr: array.Iterator(ReturnType2(op, X, Y)) = .init(&result);
            var iterx: array.Iterator(X) = .init(&x);
            const opinfo = @typeInfo(@TypeOf(op));
            for (0..result.size) |_| {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[iterr.index] = op(x.data[iterx.index], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[iterr.index] = try op(x.data[iterx.index], y, ctx);
                }

                j += 1;
                _ = iterr.nextAO(axis, iteration_order);
                _ = iterx.nextAO(axis, iteration_order);
            }
        }

        return result;
    }

    var xx: Dense(X) = undefined;
    var yy: Dense(Y) = undefined;
    if (std.mem.eql(u32, x.shape[0..x.ndim], y.shape[0..y.ndim])) {
        if (order == x.flags.order and order == y.flags.order) {
            // Trivial loop
            var result: Dense(ReturnType2(op, X, Y)) = try .init(ReturnType2(op, X, Y), allocator, x.shape[0..x.ndim], order);
            errdefer result.deinit(allocator);

            var j: u32 = 0;
            //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);
            const opinfo = @typeInfo(@TypeOf(op));
            for (0..result.size) |i| {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], y.data[i], ctx);
                }

                j += 1;
            }

            return result;
        } else {
            // Different order, but same shape
            xx = x;
            yy = y;
        }
    } else {
        //const bct = try array.broadcastShapes(&.{ x.shape[0..x.ndim], y.shape[0..y.ndim] });
        //xx = try broadcast(X, &x, bct.shape[0..bct.ndim]);
        //yy = try broadcast(Y, &y, bct.shape[0..bct.ndim]);
    }

    var result: Dense(ReturnType2(op, X, Y)) = try .init(ReturnType2(op, X, Y), allocator, xx.shape[0..xx.ndim], order);
    errdefer result.deinit(allocator);

    var j: u32 = 0;
    const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
    const axis: u32 = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
    errdefer strided.cleanup(ReturnType2(op, X, Y), &result, j, iteration_order, ctx);
    var iterr: array.Iterator(ReturnType2(op, X, Y)) = .init(&result);
    var iterx: array.Iterator(X) = .init(&xx);
    var itery: array.Iterator(Y) = .init(&yy);
    const opinfo = @typeInfo(@TypeOf(op));
    for (0..result.size) |_| {
        if (comptime opinfo.@"fn".params.len == 2) {
            result.data[iterr.index] = op(x.data[iterx.index], y.data[itery.index]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            result.data[iterr.index] = try op(x.data[iterx.index], y.data[itery.index], ctx);
        }

        j += 1;
        _ = iterr.nextAO(axis, iteration_order);
        _ = iterx.nextAO(axis, iteration_order);
        _ = itery.nextAO(axis, iteration_order);
    }

    return result;
}

pub fn apply2_(
    comptime O: type,
    o: anytype,
    comptime X: type,
    x: anytype,
    comptime Y: type,
    y: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    if (comptime !types.isDense(@TypeOf(x)) and !types.isSlice(@TypeOf(x)) and
        !types.isDense(@TypeOf(y)) and !types.isSlice(@TypeOf(y)))
    {
        const opinfo = @typeInfo(@TypeOf(op_));
        if (comptime opinfo.@"fn".params.len == 3) {
            op_(&o.data[0], x, y);
        } else if (comptime opinfo.@"fn".params.len == 4) {
            try op_(&o.data[0], x, y, ctx);
        }

        for (1..o.size) |i| {
            try ops.set(&o.data[i], o.data[0], ctx);
        }

        return;
    } else if (comptime !types.isDense(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        if (std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim]) and
            o.flags.order == y.flags.order)
        {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op_(&o.data[i], x, y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op_(&o.data[i], x, y.data[i], ctx);
                }
            }

            return;
        }

        // Different order, but same shape
        var yy: Dense(Y) = undefined;
        if (std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim])) {
            yy = y;
        } else {
            const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], y.shape[0..y.ndim] });
            if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
                return array.Error.NotBroadcastable;
            }

            //yy = try broadcast(Y, &y, bct.shape[0..bct.ndim]);
        }

        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        var itery: array.Iterator(Y) = .init(&yy);
        const opinfo = @typeInfo(@TypeOf(op_));
        for (0..o.size) |_| {
            if (comptime opinfo.@"fn".params.len == 3) {
                op_(&o.data[itero.index], x, yy.data[itery.index]);
            } else if (comptime opinfo.@"fn".params.len == 4) {
                try op_(&o.data[itero.index], x, yy.data[itery.index], ctx);
            }

            _ = itero.nextAO(axis, iteration_order);
            _ = itery.nextAO(axis, iteration_order);
        }

        return;
    } else if (comptime !types.isDense(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
            o.flags.order == x.flags.order)
        {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op_(&o.data[i], x.data[i], y);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op_(&o.data[i], x.data[i], y, ctx);
                }
            }

            return;
        }

        // Different order, but same shape
        var xx: Dense(X) = undefined;
        if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
            xx = x;
        } else {
            const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });
            if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
                return array.Error.NotBroadcastable;
            }

            //xx = try broadcast(X, &x, bct.shape[0..bct.ndim]);
        }

        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        var iterx: array.Iterator(X) = .init(&xx);
        const opinfo = @typeInfo(@TypeOf(op_));
        for (0..o.size) |_| {
            if (comptime opinfo.@"fn".params.len == 3) {
                op_(&o.data[itero.index], xx.data[iterx.index], y);
            } else if (comptime opinfo.@"fn".params.len == 4) {
                try op_(&o.data[itero.index], xx.data[iterx.index], y, ctx);
            }

            _ = itero.nextAO(axis, iteration_order);
            _ = iterx.nextAO(axis, iteration_order);
        }

        return;
    }

    var xx: Dense(X) = undefined;
    var yy: Dense(Y) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
        std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim]))
    {
        if (o.flags.order == x.flags.order and o.flags.order == y.flags.order) {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op_(&o.data[i], x.data[i], y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op_(&o.data[i], x.data[i], y.data[i], ctx);
                }
            }

            return;
        } else {
            // Different order, but same shape
            xx = x;
            yy = y;
        }
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim], y.shape[0..y.ndim] });
        if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
            return array.Error.NotBroadcastable;
        }

        //xx = try broadcast(X, &x, bct.shape[0..bct.ndim]);
        //yy = try broadcast(Y, &y, bct.shape[0..bct.ndim]);
    }

    const iteration_order: array.IterationOrder = o.flags.order.resolve3(xx.flags.order, yy.flags.order).toIterationOrder();
    const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
    var itero: array.Iterator(O) = .init(o);
    var iterx: array.Iterator(X) = .init(&xx);
    var itery: array.Iterator(Y) = .init(&yy);
    const opinfo = @typeInfo(@TypeOf(op_));
    for (0..o.size) |_| {
        if (comptime opinfo.@"fn".params.len == 3) {
            op_(&o.data[itero.index], xx.data[iterx.index], yy.data[itery.index]);
        } else if (comptime opinfo.@"fn".params.len == 4) {
            try op_(&o.data[itero.index], xx.data[iterx.index], yy.data[itery.index], ctx);
        }

        _ = itero.nextAO(axis, iteration_order);
        _ = iterx.nextAO(axis, iteration_order);
        _ = itery.nextAO(axis, iteration_order);
    }

    return;
}
