const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const scast = types.scast;
const cast = types.cast;
const needsAllocator = types.needsAllocator;
const validateContext = types.validateContext;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");
const Array = array.Array;
const Range = array.Range;

const strided = @import("strided.zig");

pub inline fn index(strides: [array.max_dimensions]usize, position: []const usize) usize {
    var idx: usize = 0;
    for (0..position.len) |i| {
        idx += position[i] * strides[i];
    }

    return idx;
}

pub inline fn checkPosition(ndim: usize, shape: [array.max_dimensions]usize, position: []const usize) !void {
    if (position.len > ndim) {
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
    data: []T,
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

    comptime if (!types.isArbitraryPrecision(T))
        return;

    for (data) |*value| {
        ops.deinit(value, ctx);
    }
}

pub inline fn init(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const usize,
    order: array.Order,
) !Array(T) {
    var size: usize = 1;
    var shapes: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    var strides: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    if (shape.len > 0) {
        for (0..shape.len) |i| {
            const idx: usize = if (order == .row_major) shape.len - i - 1 else i;

            if (shape[idx] == 1 and shape[idx] == 1) {
                strides[idx] = 0; // No stride for the new dimensions.
            } else {
                strides[idx] = size;
            }

            size *= shape[idx];

            shapes[i] = shape[i];
        }
    }

    return Array(T){
        .data = try allocator.alloc(T, size),
        .ndim = shape.len,
        .shape = shapes,
        .size = size,
        .base = null,
        .flags = .{
            .order = order,
            .storage = .dense,
            .owns_data = true,
        },
        .metadata = .{ .dense = .{
            .strides = strides,
        } },
    };
}

pub inline fn full(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const usize,
    value: anytype,
    order: array.Order,
    ctx: anytype,
) !Array(T) {
    var arr = try init(T, allocator, shape, order);
    errdefer arr.deinit(allocator);

    arr.data[0] = try ops.init(T, ctx);

    var j: usize = 1;
    errdefer cleanup(T, arr.data[0..j], ctx);

    try ops.set(&arr.data[0], value, ctx);

    for (1..arr.size) |i| {
        arr.data[i] = ops.copy(arr.data[0], ctx);

        j += 1;
    }

    return arr;
}

pub inline fn arange(
    comptime T: type,
    allocator: std.mem.Allocator,
    start: anytype,
    stop: anytype,
    step: anytype,
    ctx: anytype,
) !Array(T) {
    // Refactor to avoid unnecessary casts, like in linspace.

    const positive_step: bool = try ops.gt(step, 0, .{});
    if (try ops.eq(step, 0, .{}) or
        (try ops.lt(stop, start, .{}) and positive_step) or
        (try ops.gt(stop, start, .{}) and !positive_step))
    {
        return array.Error.InvalidRange;
    }

    var start_casted: T = try ops.init(T, ctx);
    errdefer ops.deinit(&start_casted, ctx);
    try ops.set(&start_casted, start, ctx);

    var stop_casted: T = try ops.init(T, ctx);
    errdefer ops.deinit(&stop_casted, ctx);
    try ops.set(&stop_casted, stop, ctx);

    var step_casted: T = try ops.init(T, ctx);
    errdefer ops.deinit(&step_casted, ctx);
    try ops.set(&step_casted, step, ctx);

    var diff: T = if (positive_step)
        try ops.sub(stop_casted, start_casted, ctx)
    else
        try ops.sub(start_casted, stop_casted, ctx);
    errdefer ops.deinit(&diff, ctx);

    var len_T: T = try ops.div(diff, step_casted, ctx);
    errdefer ops.deinit(&len_T, ctx);
    try ops.abs_(&len_T, len_T, ctx);
    try ops.ceil_(&len_T, len_T, ctx);
    const len: usize = scast(usize, len_T);

    if (len == 0) {
        return array.Error.InvalidRange;
    }

    var arr: Array(T) = try init(T, allocator, &.{len}, .row_major);
    errdefer arr.deinit(allocator);
    switch (len) {
        1 => {
            arr.data[0] = start_casted;

            ops.deinit(&stop_casted, ctx);
            ops.deinit(&step_casted, ctx);
            ops.deinit(&diff, ctx);
            ops.deinit(&len_T, ctx);

            return arr;
        },
        2 => {
            arr.data[0] = start_casted;
            try ops.add_(&step_casted, step_casted, arr.data[0], ctx);
            arr.data[1] = step_casted;

            ops.deinit(&stop_casted, ctx);
            ops.deinit(&diff, ctx);
            ops.deinit(&len_T, ctx);

            return arr;
        },
        3 => {
            arr.data[0] = start_casted;
            try ops.add_(&stop_casted, arr.data[0], step_casted, ctx);
            arr.data[1] = stop_casted;
            try ops.add_(&step_casted, arr.data[1], step_casted, ctx);
            arr.data[2] = step_casted;

            ops.deinit(&diff, ctx);
            ops.deinit(&len_T, ctx);

            return arr;
        },
        4 => {
            arr.data[0] = start_casted;
            try ops.add_(&diff, arr.data[0], step_casted, ctx);
            arr.data[1] = diff;
            try ops.add_(&stop_casted, arr.data[1], step_casted, ctx);
            arr.data[2] = stop_casted;
            try ops.add_(&step_casted, arr.data[2], step_casted, ctx);
            arr.data[3] = step_casted;

            ops.deinit(&len_T, ctx);

            return arr;
        },
        else => {
            arr.data[0] = start_casted;
            try ops.add_(&len_T, arr.data[0], step_casted, ctx);
            arr.data[1] = len_T;
            try ops.add_(&diff, arr.data[1], step_casted, ctx);
            arr.data[2] = diff;
            try ops.add_(&stop_casted, arr.data[2], step_casted, ctx);
            arr.data[3] = stop_casted;
        },
    }

    var j: usize = 4;
    errdefer cleanup(T, arr.data[4..j], ctx);

    for (4..len - 1) |i| {
        arr.data[i] = try ops.add(arr.data[i - 1], step_casted, ctx);
        j += 1;
    }

    try ops.add_(&step_casted, step_casted, arr.data[len - 2], ctx);
    arr.data[len - 1] = step_casted;

    return arr;
}

pub inline fn linspace(
    comptime T: type,
    allocator: std.mem.Allocator,
    start: anytype,
    stop: anytype,
    num: usize,
    endpoint: bool,
    retstep: ?*T,
    ctx: anytype,
) !Array(T) {
    var arr: Array(T) = try init(T, allocator, &.{num}, .row_major);
    errdefer arr.deinit(allocator);

    if (num == 1) {
        arr.data[0] = try ops.init(T, ctx);
        errdefer ops.deinit(&arr.data[0], ctx);
        try ops.set(&arr.data[0], start, ctx);

        return arr;
    } else if (num == 2 and endpoint) {
        arr.data[0] = try ops.init(T, ctx);
        errdefer ops.deinit(&arr.data[0], ctx);
        try ops.set(&arr.data[0], start, ctx);
        arr.data[1] = try ops.init(T, ctx);
        errdefer ops.deinit(&arr.data[1], ctx);
        try ops.set(&arr.data[1], stop, ctx);

        return arr;
    } else if (num == 2 and !endpoint) {
        arr.data[0] = try ops.init(T, ctx);
        errdefer ops.deinit(&arr.data[0], ctx);
        try ops.set(&arr.data[0], start, ctx);
        arr.data[1] = try ops.init(T, ctx);
        errdefer ops.deinit(&arr.data[1], ctx);
        try ops.set(&arr.data[1], stop, ctx);
        try ops.add_(&arr.data[1], arr.data[1], arr.data[0], ctx);
        try ops.div_(&arr.data[1], arr.data[1], 2, ctx);

        return arr;
    }

    var start_casted: T = try ops.init(T, ctx);
    errdefer ops.deinit(&start_casted, ctx);
    try ops.set(&start_casted, start, ctx);

    var stop_casted: T = try ops.init(T, ctx);
    errdefer ops.deinit(&stop_casted, ctx);
    try ops.set(&stop_casted, stop, ctx);

    var step: T = try ops.sub(stop_casted, start_casted, ctx);
    errdefer ops.deinit(&step, ctx);
    if (endpoint) {
        try ops.div_(&step, step, num - 1, ctx);
    } else {
        try ops.div_(&step, step, num, ctx);
    }

    if (retstep) |*s| {
        try ops.set(s, step, ctx);
    }

    if (num == 3 and endpoint) {
        arr.data[0] = start_casted;
        try ops.add_(&step, step, arr.data[0], ctx);
        arr.data[1] = step;
        arr.data[2] = stop_casted;

        return arr;
    } else if (num == 3 and !endpoint) {
        arr.data[0] = start_casted;
        try ops.sub_(&stop_casted, stop_casted, step, ctx);
        arr.data[2] = stop_casted;
        try ops.add_(&step, step, arr.data[0], ctx);
        arr.data[1] = step;

        return arr;
    }

    arr.data[0] = start_casted;

    var j: usize = 1;
    errdefer cleanup(T, allocator, arr.data[1..j]);
    for (1..num - 2) |i| {
        arr.data[i] = try ops.add(arr.data[i - 1], step, ctx);

        j += 1;
    }

    if (endpoint) {
        try ops.add_(&step, step, arr.data[num - 3], ctx);
        arr.data[num - 2] = step;
        arr.data[num - 1] = stop_casted;
    } else {
        try ops.sub_(&stop_casted, stop_casted, step, ctx);
        try ops.add_(&step, step, arr.data[num - 3], ctx);
        arr.data[num - 2] = step;
        arr.data[num - 1] = stop_casted;
    }

    return arr;
}

pub inline fn logspace(
    comptime T: type,
    allocator: std.mem.Allocator,
    start: anytype,
    stop: anytype,
    num: usize,
    base: anytype,
    endpoint: bool,
    ctx: anytype,
) !Array(T) {
    var arr: Array(T) = try linspace(T, allocator, start, stop, num, endpoint, null, ctx);
    errdefer arr.deinit(allocator);

    try ops.pow_(&arr, base, arr, ctx);

    return arr;
}

pub inline fn set(comptime T: type, arr: *Array(T), position: []const usize, value: T) !void {
    try checkPosition(arr.ndim, arr.shape, position);

    arr.data[index(arr.metadata.dense.strides, position)] = value;
}

pub inline fn get(comptime T: type, arr: Array(T), position: []const usize) !*T {
    try checkPosition(arr.ndim, arr.shape, position);

    return &arr.data[index(arr.metadata.dense.strides, position)];
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
            .storage = .strided, // Should be dense?
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

pub fn ravel(comptime T: type, arr: *const Array(T)) !Array(T) {
    return Array(T){
        .data = arr.data,
        .ndim = 1,
        .shape = .{arr.size} ++ .{0} ** (array.max_dimensions - 1),
        .size = arr.size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order,
            .storage = .dense,
            .owns_data = false,
        },
        .metadata = .{ .dense = .{
            .strides = .{1} ++ .{0} ** (array.max_dimensions - 1),
        } },
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
        new_strides[i] = scast(isize, self.metadata.dense.strides[idx]);
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
                .offset = 0, // No offset in transposing.
            },
        },
    };
}

pub inline fn slice(comptime T: type, arr: *const Array(T), ranges: []const Range) !Array(T) {
    var ndim: usize = arr.ndim;
    var size: usize = 1;
    var shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    var strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
    var offset: usize = 0;

    var i: usize = 0;
    var j: usize = 0;
    while (i < arr.ndim) {
        const stride: isize = scast(isize, arr.metadata.dense.strides[i]);

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
            strides[scast(usize, i)] = scast(isize, arr.metadata.dense.strides[scast(usize, i - diff)]);
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
                .offset = 0, // No offset in broadcasting.
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
    var result: Array(ReturnType1(op, T)) = try init(ReturnType1(op, T), allocator, x.shape[0..x.ndim], order);
    errdefer result.deinit(allocator);

    var j: usize = 0;
    if (result.flags.order == x.flags.order) {
        // Trivial loop
        errdefer cleanup(ReturnType1(op, T), allocator, result.data[0..j]);

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
    comptime op_: anytype,
    ctx: anytype,
) !void {
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
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

    var xx: Array(X) = undefined;
    if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
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
        if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
            return array.Error.NotBroadcastable;
        }

        xx = try broadcast(X, &x, bct.shape[0..bct.ndim]);
    }

    const iteration_order: array.IterationOrder = o.flags.order.toIterationOrder();
    const axis: usize = if (iteration_order == .right_to_left) o.ndim - 1 else 0;
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
) !Array(ReturnType2(op, X, Y)) {
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        var result: Array(ReturnType2(op, X, Y)) = try init(ReturnType2(op, X, Y), allocator, y.shape[0..y.ndim], order);
        errdefer result.deinit(allocator);

        var j: usize = 0;
        if (result.flags.order == y.flags.order) {
            // Trivial loop
            errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

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
            const axis: usize = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
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
    } else if (comptime !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        var result: Array(ReturnType2(op, X, Y)) = try init(ReturnType2(op, X, Y), allocator, x.shape[0..x.ndim], order);
        errdefer result.deinit(allocator);

        var j: usize = 0;
        if (result.flags.order == x.flags.order) {
            errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

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
            const axis: usize = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
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

    var xx: Array(X) = undefined;
    var yy: Array(Y) = undefined;
    if (std.mem.eql(usize, x.shape[0..x.ndim], y.shape[0..y.ndim])) {
        if (order == x.flags.order and order == y.flags.order) {
            // Trivial loop
            var result: Array(ReturnType2(op, X, Y)) = try init(ReturnType2(op, X, Y), allocator, x.shape[0..x.ndim], order);
            errdefer result.deinit(allocator);

            var j: usize = 0;
            errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);
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
        const bct = try array.broadcastShapes(&.{ x.shape[0..x.ndim], y.shape[0..y.ndim] });
        xx = try broadcast(X, &x, bct.shape[0..bct.ndim]);
        yy = try broadcast(Y, &y, bct.shape[0..bct.ndim]);
    }

    var result: Array(ReturnType2(op, X, Y)) = try init(ReturnType2(op, X, Y), allocator, xx.shape[0..xx.ndim], order);
    errdefer result.deinit(allocator);

    var j: usize = 0;
    const iteration_order: array.IterationOrder = result.flags.order.toIterationOrder();
    const axis: usize = if (iteration_order == .right_to_left) result.ndim - 1 else 0;
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
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x)) and
        !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y)))
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
    } else if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        if (std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim]) and
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
    } else if (comptime !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
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

    var xx: Array(X) = undefined;
    var yy: Array(Y) = undefined;
    if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
        std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim]))
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
