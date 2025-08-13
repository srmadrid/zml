const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;

const ops = @import("../ops.zig");
const int = @import("../int.zig");

const array = @import("../array.zig");
const Order = types.Order;
const Flags = array.Flags;
const Range = array.Range;

const dense = @import("dense.zig");
const Dense = dense.Dense;

pub fn Strided(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("Strided requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        ndim: u32,
        shape: [array.max_dimensions]u32,
        strides: [array.max_dimensions]i32,
        size: u32,
        offset: u32,
        base: ?*const anyopaque,
        flags: Flags = .{},

        pub const empty: Strided(T) = .{
            .data = &.{},
            .ndim = 0,
            .shape = .{0} ** array.max_dimensions,
            .strides = .{0} ** array.max_dimensions,
            .size = 0,
            .offset = 0,
            .base = null,
            .flags = .{},
        };

        pub fn get(self: *const Strided(T), position: []const u32) !T {
            try self._checkPosition(position);

            return self.data[self._index(position)];
        }

        pub inline fn at(self: *const Strided(T), position: []const u32) T {
            // Unchecked version of get. Assumes position is valid.
            return self.data[self._index(position)];
        }

        pub fn set(self: *const Strided(T), position: []const u32, value: T) !void {
            try self._checkPosition(position);

            self.data[self._index(position)] = value;
        }

        pub inline fn put(self: *const Strided(T), position: []const u32, value: T) void {
            // Unchecked version of set. Assumes position is valid.
            self.data[self._index(position)] = value;
        }

        pub inline fn reshape(self: *const Strided(T), shape: []const u32) !Strided(T) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len == 0)
                return array.Error.ZeroDimension;

            var new_size: u32 = 1;
            var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            if (shape.len > 0) {
                for (0..shape.len) |i| {
                    const idx: u32 = if (self.flags.order == .row_major) shape.len - i - 1 else i;

                    new_strides[idx] = types.scast(i32, new_size);
                    new_size *= shape[idx];

                    new_shape[i] = shape[i];
                }
            }

            if (new_size != self.size) {
                return error.DimensionMismatch;
            }

            return Strided(T){
                .data = self.data,
                .ndim = shape.len,
                .shape = new_shape,
                .strides = new_strides,
                .size = new_size,
                .offset = self.offset,
                .base = if (self.flags.owns_data) self else self.base,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }

        pub fn transpose(self: *const Strided(T), axes: ?[]const u32) !Strided(T) {
            const axes_: []const u32 =
                axes orelse
                array.trivialReversePermutation(self.ndim)[0..self.ndim];

            if (axes_.len == 0)
                return array.Error.ZeroDimension;

            if (axes_.len != self.ndim)
                return array.Error.DimensionMismatch;

            if (!array.isPermutation(self.ndim, axes_))
                return array.Error.InvalidAxes; // axes must be a valid permutation of [0, ..., ndim - 1]

            var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            var size: u32 = 1;

            var i: u32 = 0;
            while (i < self.ndim) : (i += 1) {
                const idx: u32 = axes_[i];

                new_shape[i] = self.shape[idx];
                new_strides[i] = self.metadata.strided.strides[idx];
                size *= new_shape[i];
            }

            return Strided(T){
                .data = self.data,
                .ndim = self.ndim,
                .shape = new_shape,
                .strides = new_strides,
                .size = size,
                .offset = self.metadata.strided.offset,
                .base = if (self.flags.owns_data) self else self.base,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }

        pub inline fn broadcast(self: *const Strided(T), shape: []const u32) !Strided(T) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len < self.ndim)
                return array.Error.TooLittleDimensions;

            var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            var size: u32 = 1;

            var i: i32 = types.scast(i32, shape.len - 1);
            const diff: i32 = types.scast(i32, shape.len - self.ndim);
            while (i >= 0) : (i -= 1) {
                if (i - diff >= 0) {
                    new_shape[types.scast(u32, i)] = int.max(self.shape[types.scast(u32, i - diff)], shape[types.scast(u32, i)]);
                    strides[types.scast(u32, i)] = self.metadata.strided.strides[types.scast(u32, i - diff)];
                } else {
                    new_shape[types.scast(u32, i)] = shape[types.scast(u32, i)];
                    strides[types.scast(u32, i)] = 0; // No stride for the new dimensions.
                }

                size *= new_shape[types.scast(u32, i)];
            }

            return Strided(T){
                .data = self.data,
                .ndim = shape.len,
                .shape = new_shape,
                .strides = strides,
                .size = size,
                .offset = self.metadata.strided.offset,
                .base = if (self.flags.owns_data) self else self.base,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }

        pub fn slice(self: *const Strided(T), ranges: []const Range) !Strided(T) {
            if (ranges.len == 0 or ranges.len > self.ndim)
                return error.DimensionMismatch;

            var ndim: u32 = self.ndim;
            var size: u32 = 1;
            var shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            var offset: u32 = self.metadata.strided.offset;

            var i: u32 = 0;
            var j: u32 = 0;
            while (i < self.ndim) {
                const stride: i32 = self.metadata.strided.strides[i];

                if (i >= ranges.len) {
                    shape[j] = self.shape[i];
                    strides[j] = stride;
                    size *= self.shape[i];
                    j += 1;
                    i += 1;
                    continue;
                }

                var range: Range = ranges[i];
                if (range.start != int.maxVal(u32) and range.start == range.stop) {
                    return array.Error.InvalidRange;
                } else if (range.step > 0) {
                    if (range.start != int.maxVal(u32) and range.start >= self.shape[i] or
                        (range.stop != int.maxVal(u32) and range.stop > self.shape[i]))
                        return array.Error.RangeOutOfBounds;
                } else if (range.step < 0) {
                    if ((range.stop != int.maxVal(u32) and range.stop >= self.shape[i]) or
                        (range.start != int.maxVal(u32) and range.start > self.shape[i]))
                        return array.Error.RangeOutOfBounds;
                }

                var len_adjustment: u32 = 0;
                if (range.step > 0) {
                    if (range.start == int.maxVal(u32)) {
                        range.start = 0;
                    }

                    if (range.stop == int.maxVal(u32)) {
                        range.stop = self.shape[i];
                    }
                } else if (range.step < 0) {
                    if (range.start == int.maxVal(u32)) {
                        range.start = self.shape[i] - 1;
                    }

                    if (range.stop == int.maxVal(u32)) {
                        range.stop = 0;
                        len_adjustment = 1;
                    }
                }

                const len: u32 = range.len() + len_adjustment;
                if (len == 1) {
                    ndim -= 1;
                } else {
                    shape[j] = len;
                    strides[j] = stride * range.step;
                    size *= len;
                    j += 1;
                }

                if (stride < 0) {
                    offset -= range.start * types.scast(u32, int.abs(stride));
                } else {
                    offset += range.start * types.scast(u32, stride);
                }

                i += 1;
            }

            return Strided(T){
                .data = self.data,
                .ndim = ndim,
                .shape = shape,
                .strides = strides,
                .size = size,
                .offset = offset,
                .base = if (self.flags.owns_data) self else self.base,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }

        inline fn _index(
            self: *const Strided(T),
            position: []const u32,
        ) u32 {
            var idx: u32 = self.offset;
            var i: u32 = 0;
            while (i < position.len) : (i += 1) {
                const stride: i32 = self.strides[i];
                if (stride < 0) {
                    idx -= position[i] * types.scast(u32, int.abs(stride));
                } else {
                    idx += position[i] * types.scast(u32, stride);
                }
            }

            return idx;
        }

        inline fn _checkPosition(self: *const Strided(T), position: []const u32) !void {
            if (position.len > self.ndim)
                return array.Error.DimensionMismatch;

            var i: u32 = 0;
            while (i < position.len) : (i += 1) {
                if (position[i] >= self.shape[i]) {
                    return array.Error.PositionOutOfBounds;
                }
            }
        }
    };
}

pub inline fn cleanup(
    comptime T: type,
    arr: *const Strided(T),
    num_elements: u32,
    iteration_order: array.IterationOrder,
    ctx: anytype,
) void {
    comptime if (types.isArbitraryPrecision(T)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    comptime if (!types.needsAllocator(T))
        return;

    var iter: array.Iterator(T) = .init(arr);
    var axis: u32 = 0;
    if (iteration_order == .right_to_left) {
        axis = arr.ndim - 1;
    }
    for (0..num_elements) |_| {
        ops.deinit(arr.data[iter.index], ctx);

        _ = iter.nextAO(axis, iteration_order);
    }
}

fn loop1(
    result: anytype,
    x: anytype,
    comptime op: anytype,
    depth: u32,
    order: types.IterationOrder,
    ir: u32,
    ix: i32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (order == .left_to_right) 0 else result.ndim - 1;

        var jr: u32 = ir;
        var jx: i32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime opinfo.@"fn".params.len == 1) {
                result.data[jr] = op(x.data[types.scast(u32, jx)]);
            } else if (comptime opinfo.@"fn".params.len == 2) {
                result.data[jr] = try op(x.data[types.scast(u32, jx)], ctx);
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    } else {
        const idx: u32 = if (order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: u32 = ir;
        var jx: i32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            try loop1(
                result,
                x,
                op,
                depth - 1,
                order,
                jr,
                jx,
                ctx,
            );

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    }
}

pub fn apply1(
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    opts: struct {
        order: ?Order = null,
    },
    ctx: anytype,
) !Dense(ReturnType1(op, Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    var result: Dense(ReturnType1(op, X)) = try .init(allocator, x.shape[0..x.ndim], .{ .order = opts.order orelse x.flags.order });
    errdefer result.deinit(allocator);

    try loop1(
        &result,
        &x,
        op,
        result.ndim - 1,
        result.flags.order.toIterationOrder(),
        0,
        types.scast(i32, x.offset),
        ctx,
    );

    return result;
}

pub fn apply1_(
    comptime O: type,
    o: anytype,
    comptime X: type,
    x: anytype,
    comptime op_to: anytype,
    ctx: anytype,
) !void {
    if (comptime !types.isStrided(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        const first: u32 = itero.index;

        const opinfo = @typeInfo(@TypeOf(op_to));
        if (comptime opinfo.@"fn".params.len == 2) {
            op_to(&o.data[first], x);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            try op_to(&o.data[first], x, ctx);
        }
        _ = itero.nextAO(axis, iteration_order);

        for (1..o.size) |_| {
            try ops.set(&o.data[itero.index], o.data[first], ctx);

            _ = itero.nextAO(axis, iteration_order);
        }

        return;
    }

    const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });
    if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
        return array.Error.NotBroadcastable;
    }

    const xx: Strided(X) = try x.broadcast(bct.shape[0..bct.ndim]);

    const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
    const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
    var itero: array.Iterator(O) = .init(o);
    var iterx: array.Iterator(X) = .init(&xx);

    const opinfo = @typeInfo(@TypeOf(op_to));
    for (0..o.size) |_| {
        if (comptime opinfo.@"fn".params.len == 2) {
            op_to(&o.data[itero.index], xx.data[iterx.index]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            try op_to(&o.data[itero.index], xx.data[iterx.index], ctx);
        }

        _ = itero.nextAO(axis, iteration_order);
        _ = iterx.nextAO(axis, iteration_order);
    }

    return;
}

pub fn apply2(
    allocator: std.mem.Allocator,
    comptime X: type,
    x: anytype,
    comptime Y: type,
    y: anytype,
    comptime op: anytype,
    order: array.Order,
    ctx: anytype,
) !Strided(ReturnType2(op, X, Y)) {
    if (comptime !types.isStrided(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        var result: Strided(ReturnType2(op, X, Y)) = try dense.init(ReturnType2(op, X, Y), allocator, y.shape[0..y.ndim], order);
        errdefer result.deinit(allocator);

        var j: u32 = 0;
        const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
        const axis: u32 = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
        errdefer cleanup(ReturnType2(op, X, Y), &result, j, iteration_order, ctx);
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

        return result;
    } else if (comptime !types.isStrided(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        var result: Strided(ReturnType2(op, X, Y)) = try dense.init(ReturnType2(op, X, Y), allocator, x.shape[0..y.ndim], order);
        errdefer result.deinit(allocator);

        var j: u32 = 0;
        const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
        const axis: u32 = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
        errdefer cleanup(ReturnType2(op, X, Y), &result, j, iteration_order, ctx);
        var iterr: array.Iterator(ReturnType2(op, X, Y)) = .init(&result);
        var iterx = array.Iterator(X).init(&x);
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

        return result;
    }

    const bct = try array.broadcastShapes(&.{ x.shape[0..x.ndim], y.shape[0..y.ndim] });
    const xx: Strided(X) = if (x.flags.storage == .dense)
        try x.broadcast(X, &x, bct.shape[0..bct.ndim])
    else
        undefined;
    // else
    //     try broadcast(X, &x, bct.shape[0..bct.ndim]);
    const yy: Strided(Y) = if (y.flags.storage == .dense)
        try y.broadcast(Y, &y, bct.shape[0..bct.ndim])
    else
        undefined;
    // else
    //     try broadcast(Y, &y, bct.shape[0..bct.ndim]);

    var result: Strided(ReturnType2(op, X, Y)) = try dense.init(ReturnType2(op, X, Y), allocator, bct.shape[0..bct.ndim], order);
    errdefer result.deinit(allocator);

    var j: u32 = 0;
    const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
    const axis: u32 = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
    errdefer cleanup(ReturnType2(op, X, Y), &result, j, iteration_order, ctx);
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
    comptime op_to: anytype,
    ctx: anytype,
) !void {
    if (comptime !types.isStrided(@TypeOf(x)) and !types.isSlice(@TypeOf(x)) and
        !types.isStrided(@TypeOf(y)) and !types.isSlice(@TypeOf(y)))
    {
        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        const first: u32 = itero.index;
        const opinfo = @typeInfo(@TypeOf(op_to));
        if (comptime opinfo.@"fn".params.len == 3) {
            op_to(&o.data[first], x, y);
        } else if (comptime opinfo.@"fn".params.len == 4) {
            try op_to(&o.data[first], x, y, ctx);
        }
        _ = itero.nextAO(axis, iteration_order);

        for (1..o.size) |_| {
            try ops.set(&o.data[itero.index], o.data[first], ctx);

            _ = itero.nextAO(axis, iteration_order);
        }

        return;
    } else if (comptime !types.isStrided(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        var yy: Strided(Y) = undefined;
        if (std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim])) {
            yy = y;
        } else {
            const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], y.shape[0..y.ndim] });
            if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
                return array.Error.NotBroadcastable;
            }

            yy = try y.broadcast(Y, &y, bct.shape[0..bct.ndim]);
        }

        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        var itery: array.Iterator(Y) = .init(&yy);
        const opinfo = @typeInfo(@TypeOf(op_to));
        for (0..o.size) |_| {
            if (comptime opinfo.@"fn".params.len == 3) {
                op_to(&o.data[itero.index], x, yy.data[itery.index]);
            } else if (comptime opinfo.@"fn".params.len == 4) {
                try op_to(&o.data[itero.index], x, yy.data[itery.index], ctx);
            }

            _ = itero.nextAO(axis, iteration_order);
            _ = itery.nextAO(axis, iteration_order);
        }

        return;
    } else if (comptime !types.isStrided(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        var xx: Strided(X) = undefined;
        if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
            xx = x;
        } else {
            const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });
            if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
                return array.Error.NotBroadcastable;
            }

            xx = try x.broadcast(X, &x, bct.shape[0..bct.ndim]);
        }

        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        var iterx: array.Iterator(X) = .init(&xx);
        const opinfo = @typeInfo(@TypeOf(op_to));
        for (0..o.size) |_| {
            if (comptime opinfo.@"fn".params.len == 3) {
                op_to(&o.data[itero.index], xx.data[iterx.index], y);
            } else if (comptime opinfo.@"fn".params.len == 4) {
                try op_to(&o.data[itero.index], xx.data[iterx.index], y, ctx);
            }

            _ = itero.nextAO(axis, iteration_order);
            _ = iterx.nextAO(axis, iteration_order);
        }

        return;
    }

    var xx: Strided(X) = undefined;
    var yy: Strided(Y) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
        std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim]))
    {
        xx = x;
        yy = y;
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim], y.shape[0..y.ndim] });
        if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
            return array.Error.NotBroadcastable;
        }

        xx = try x.broadcast(X, &x, bct.shape[0..bct.ndim]);
        yy = try y.broadcast(Y, &y, bct.shape[0..bct.ndim]);
    }

    const iteration_order: array.IterationOrder = o.flags.order.resolve3(xx.flags.order, yy.flags.order).toIterationOrder();
    const axis: u32 = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
    var itero: array.Iterator(O) = .init(o);
    var iterx: array.Iterator(X) = .init(&xx);
    var itery: array.Iterator(Y) = .init(&yy);
    const opinfo = @typeInfo(@TypeOf(op_to));
    for (0..o.size) |_| {
        if (comptime opinfo.@"fn".params.len == 3) {
            op_to(&o.data[itero.index], xx.data[iterx.index], yy.data[itery.index]);
        } else if (comptime opinfo.@"fn".params.len == 4) {
            try op_to(&o.data[itero.index], xx.data[iterx.index], yy.data[itery.index], ctx);
        }

        _ = itero.nextAO(axis, iteration_order);
        _ = iterx.nextAO(axis, iteration_order);
        _ = itery.nextAO(axis, iteration_order);
    }

    return;
}
