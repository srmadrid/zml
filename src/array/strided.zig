const std = @import("std");

const types = @import("../types.zig");
const Scalar = types.Scalar;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const needsAllocator = types.needsAllocator;
const scast = types.scast;
const cast = types.cast;
const validateContext = types.validateContext;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");
const Array = array.Array;
const Range = array.Range;

const dense = @import("dense.zig");

pub inline fn index(strides: [array.max_dimensions]isize, offset: usize, position: []const usize) usize {
    var idx: usize = offset;
    for (0..position.len) |i| {
        const stride: isize = strides[i];
        if (stride < 0) {
            idx -= position[i] * scast(usize, int.abs(stride));
        } else {
            idx += position[i] * scast(usize, stride);
        }
    }

    return idx;
}

pub inline fn checkPosition(ndim: usize, shape: [array.max_dimensions]usize, position: []const usize) !void {
    if (position.len != ndim) {
        return array.Error.DimensionMismatch;
    }

    for (0..position.len) |i| {
        if (position[i] >= shape[i]) {
            return array.Error.PositionOutOfBounds;
        }
    }
}

pub inline fn cleanup(
    comptime T: type,
    arr: *const Array(T),
    num_elements: usize,
    iteration_order: array.IterationOrder,
    ctx: anytype,
) void {
    comptime if (types.isArbitraryPrecision(T)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    comptime if (!needsAllocator(T))
        return;

    var iter: array.Iterator(T) = .init(arr);
    var axis: usize = 0;
    if (iteration_order == .right_to_left) {
        axis = arr.ndim - 1;
    }
    for (0..num_elements) |_| {
        ops.deinit(arr.data[iter.index], ctx);

        _ = iter.nextAO(axis, iteration_order);
    }
}

pub inline fn set(comptime T: type, arr: *Array(T), position: []const usize, value: T) !void {
    try checkPosition(arr.ndim, arr.shape, position);

    arr.data[index(arr.metadata.strided.strides, arr.metadata.strided.offset, position)] = value;
}

pub fn get(comptime T: type, arr: Array(T), position: []const usize) !*T {
    try checkPosition(arr.ndim, arr.shape, position);

    return &arr.data[index(arr.metadata.strided.strides, arr.metadata.strided.offset, position)];
}

pub inline fn slice(comptime T: type, arr: *const Array(T), ranges: []const Range) !Array(T) {
    var ndim: usize = arr.ndim;
    var size: usize = 1;
    var shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    var strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
    var offset: usize = arr.metadata.strided.offset;

    var i: usize = 0;
    var j: usize = 0;
    while (i < arr.ndim) {
        const stride: isize = arr.metadata.strided.strides[i];

        if (i >= ranges.len) {
            shape[j] = arr.shape[i];
            strides[j] = stride;
            size *= arr.shape[i];
            j += 1;
            i += 1;
            continue;
        }

        var range: Range = ranges[i];
        if (range.start != int.maxVal(usize) and range.start == range.stop) {
            return array.Error.InvalidRange;
        } else if (range.step > 0) {
            if (range.start != int.maxVal(usize) and range.start >= arr.shape[i] or
                (range.stop != int.maxVal(usize) and range.stop > arr.shape[i]))
                return array.Error.RangeOutOfBounds;
        } else if (range.step < 0) {
            if ((range.stop != int.maxVal(usize) and range.stop >= arr.shape[i]) or
                (range.start != int.maxVal(usize) and range.start > arr.shape[i]))
                return array.Error.RangeOutOfBounds;
        }

        var len_adjustment: usize = 0;
        if (range.step > 0) {
            if (range.start == int.maxVal(usize)) {
                range.start = 0;
            }

            if (range.stop == int.maxVal(usize)) {
                range.stop = arr.shape[i];
            }
        } else if (range.step < 0) {
            if (range.start == int.maxVal(usize)) {
                range.start = arr.shape[i] - 1;
            }

            if (range.stop == int.maxVal(usize)) {
                range.stop = 0;
                len_adjustment = 1;
            }
        }

        const len: usize = range.len() + len_adjustment;
        if (len == 1) {
            ndim -= 1;
        } else {
            shape[j] = len;
            strides[j] = stride * range.step;
            size *= len;
            j += 1;
        }

        if (stride < 0) {
            offset -= range.start * scast(usize, int.abs(stride));
        } else {
            offset += range.start * scast(usize, stride);
        }

        i += 1;
    }

    return Array(T){
        .data = arr.data,
        .ndim = ndim,
        .shape = shape,
        .size = size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided,
            .owns_data = false,
        },
        .metadata = .{ .strided = .{
            .strides = strides,
            .offset = offset,
        } },
    };
}

pub inline fn reshape(comptime T: type, arr: *const Array(T), shape: []const usize) !Array(T) {
    var new_size: usize = 1;
    var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    var new_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
    if (shape.len > 0) {
        for (0..shape.len) |i| {
            const idx: usize = if (arr.flags.order == .row_major) shape.len - i - 1 else i;

            new_strides[idx] = scast(isize, new_size);
            new_size *= shape[idx];

            new_shape[i] = shape[i];
        }
    }

    if (new_size != arr.size) {
        return error.DimensionMismatch;
    }

    return Array(T){
        .data = arr.data,
        .ndim = shape.len,
        .shape = new_shape,
        .size = new_size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided,
            .owns_data = false,
        },
        .metadata = .{
            .strided = .{
                .strides = new_strides,
                .offset = 0, // No offset in reshaping.
            },
        },
    };
}

pub inline fn broadcast(
    comptime T: type,
    arr: *const Array(T),
    shape: []const usize,
) !Array(T) {
    var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    var strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
    var size: usize = 1;

    var i: isize = scast(isize, shape.len - 1);
    const diff: isize = scast(isize, shape.len - arr.ndim);
    while (i >= 0) : (i -= 1) {
        if (i - diff >= 0) {
            new_shape[scast(usize, i)] = try ops.max(arr.shape[scast(usize, i - diff)], shape[scast(usize, i)], .{});
            strides[scast(usize, i)] = arr.metadata.strided.strides[scast(usize, i - diff)];
        } else {
            new_shape[scast(usize, i)] = shape[scast(usize, i)];
            strides[scast(usize, i)] = 0; // No stride for the new dimensions.
        }

        size *= new_shape[scast(usize, i)];
    }

    return Array(T){
        .data = arr.data,
        .ndim = shape.len,
        .shape = new_shape,
        .size = size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided,
            .owns_data = false,
        },
        .metadata = .{
            .strided = .{
                .strides = strides,
                .offset = arr.metadata.strided.offset, // No change in offset during broadcast.
            },
        },
    };
}

pub fn transpose(
    comptime T: type,
    self: *const Array(T),
    axes: []const usize,
) !Array(T) {
    var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    var new_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
    var size: usize = 1;

    for (0..self.ndim) |i| {
        const idx: usize = axes[i];

        new_shape[i] = self.shape[idx];
        new_strides[i] = self.metadata.strided.strides[idx];
        size *= new_shape[i];
    }

    return Array(T){
        .data = self.data,
        .ndim = self.ndim,
        .shape = new_shape,
        .size = size,
        .base = if (self.flags.owns_data) self else self.base,
        .flags = .{
            .order = self.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided,
            .owns_data = false,
        },
        .metadata = .{
            .strided = .{
                .strides = new_strides,
                .offset = self.metadata.strided.offset, // No change in offset during transpose.
            },
        },
    };
}

pub inline fn apply1(
    comptime T: type,
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    order: array.Order,
    ctx: anytype,
) !Array(ReturnType1(op, T)) {
    var result: Array(ReturnType1(op, T)) = try dense.init(ReturnType1(op, T), allocator, x.shape[0..x.ndim], order);
    errdefer result.deinit(allocator);

    var j: usize = 0;
    if (result.flags.order == x.flags.order) {
        errdefer dense.cleanup(ReturnType1(op, T), result.data[0..j], ctx);
        var iter = array.Iterator(T).init(&x);

        const opinfo = @typeInfo(@TypeOf(op));
        for (0..result.size) |i| {
            if (comptime opinfo.@"fn".params.len == 1) {
                result.data[i] = op(x.data[iter.index]);
            } else if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = try op(x.data[iter.index], ctx);
            }

            j += 1;
            _ = iter.next();
        }
    } else {
        const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
        errdefer cleanup(ReturnType1(op, T), &result, j, iteration_order, ctx);
        const axis: usize = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
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
    comptime op_to: anytype,
    ctx: anytype,
) !void {
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: usize = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        const first: usize = itero.index;

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
    if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
        return array.Error.NotBroadcastable;
    }

    const xx: Array(X) = try x.broadcast(bct.shape[0..bct.ndim]);

    const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
    const axis: usize = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
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
) !Array(ReturnType2(op, X, Y)) {
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        var result: Array(ReturnType2(op, X, Y)) = try dense.init(ReturnType2(op, X, Y), allocator, y.shape[0..y.ndim], order);
        errdefer result.deinit(allocator);

        var j: usize = 0;
        const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
        const axis: usize = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
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
    } else if (comptime !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        var result: Array(ReturnType2(op, X, Y)) = try dense.init(ReturnType2(op, X, Y), allocator, x.shape[0..y.ndim], order);
        errdefer result.deinit(allocator);

        var j: usize = 0;
        const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
        const axis: usize = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
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
    const xx: Array(X) = if (x.flags.storage == .dense)
        try dense.broadcast(X, &x, bct.shape[0..bct.ndim])
    else
        try broadcast(X, &x, bct.shape[0..bct.ndim]);
    const yy: Array(Y) = if (y.flags.storage == .dense)
        try dense.broadcast(Y, &y, bct.shape[0..bct.ndim])
    else
        try broadcast(Y, &y, bct.shape[0..bct.ndim]);

    var result: Array(ReturnType2(op, X, Y)) = try dense.init(ReturnType2(op, X, Y), allocator, bct.shape[0..bct.ndim], order);
    errdefer result.deinit(allocator);

    var j: usize = 0;
    const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
    const axis: usize = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
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
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x)) and
        !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y)))
    {
        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: usize = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        const first: usize = itero.index;
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
    } else if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        var yy: Array(Y) = undefined;
        if (std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim])) {
            yy = y;
        } else {
            const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], y.shape[0..y.ndim] });
            if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
                return array.Error.NotBroadcastable;
            }

            yy = try broadcast(Y, &y, bct.shape[0..bct.ndim]);
        }

        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: usize = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
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
    } else if (comptime !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        var xx: Array(X) = undefined;
        if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
            xx = x;
        } else {
            const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });
            if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
                return array.Error.NotBroadcastable;
            }

            xx = try broadcast(X, &x, bct.shape[0..bct.ndim]);
        }

        const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
        const axis: usize = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
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

    var xx: Array(X) = undefined;
    var yy: Array(Y) = undefined;
    if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
        std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim]))
    {
        xx = x;
        yy = y;
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim], y.shape[0..y.ndim] });
        if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
            return array.Error.NotBroadcastable;
        }

        xx = try broadcast(X, &x, bct.shape[0..bct.ndim]);
        yy = try broadcast(Y, &y, bct.shape[0..bct.ndim]);
    }

    const iteration_order: array.IterationOrder = o.flags.order.resolve3(xx.flags.order, yy.flags.order).toIterationOrder();
    const axis: usize = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
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
