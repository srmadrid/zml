const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const scast = types.scast;
const cast = types.cast;
const needsAllocator = types.needsAllocator;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");
const Array = array.Array;
const Range = array.Range;

const strided = @import("strided.zig");

pub inline fn index(strides: [array.maxDimensions]usize, position: []const usize) usize {
    var idx: usize = 0;
    for (0..position.len) |i| {
        idx += position[i] * strides[i];
    }

    return idx;
}

pub inline fn checkPosition(ndim: usize, shape: [array.maxDimensions]usize, position: []const usize) !void {
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
    allocator: ?std.mem.Allocator,
    data: []T,
) void {
    comptime if (!needsAllocator(T))
        return;

    for (data) |*value| {
        ops.deinit(value, .{ .allocator = allocator.? });
    }
}

pub inline fn init(
    allocator: std.mem.Allocator,
    comptime T: type,
    shape: []const usize,
    order: array.Order,
) !Array(T) {
    var size: usize = 1;
    var shapes: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    var strides: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    if (shape.len > 0) {
        for (0..shape.len) |i| {
            const idx: usize = if (order == .rowMajor) shape.len - i - 1 else i;

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
            .ownsData = true,
            .writeable = true,
        },
        .metadata = .{ .dense = .{
            .strides = strides,
        } },
    };
}

pub inline fn full(
    allocator: std.mem.Allocator,
    comptime T: type,
    shape: []const usize,
    value: anytype,
    order: array.Order,
) !Array(T) {
    var arr = try init(allocator, T, shape, order);

    const value_casted: T = try cast(T, value, .{ .allocator = allocator, .copy = true });
    arr.data[0] = value_casted;
    for (1..arr.size) |i| {
        arr.data[i] = ops.copy(value_casted, .{ .allocator = allocator });
    }

    return arr;
}

pub inline fn arange(
    allocator: std.mem.Allocator,
    comptime T: type,
    start: anytype,
    stop: anytype,
    step: anytype,
    writeable: bool,
) !Array(T) {
    // Refactor to avoid unnecessary casts, like in linspace.
    _ = types.numericType(@TypeOf(start));
    _ = types.numericType(@TypeOf(stop));
    _ = types.numericType(@TypeOf(step));

    const positive_step: bool = try ops.gt(step, 0, .{});
    if (try ops.eq(step, 0, .{}) or
        (try ops.lt(stop, start, .{}) and positive_step) or
        (try ops.gt(stop, start, .{}) and !positive_step))
    {
        return array.Error.InvalidRange;
    }

    var start_casted: T = try cast(T, start, .{ .allocator = allocator, .copy = true });
    errdefer ops.deinit(&start_casted, .{ .allocator = allocator });
    var stop_casted: T = try cast(T, stop, .{ .allocator = allocator, .copy = true });
    errdefer ops.deinit(&stop_casted, .{ .allocator = allocator });
    var step_casted: T = try cast(T, step, .{ .allocator = allocator, .copy = true });
    errdefer ops.deinit(&step_casted, .{ .allocator = allocator });

    var diff: T = if (positive_step)
        try ops.sub(stop_casted, start_casted, .{ .allocator = allocator })
    else
        try ops.sub(start_casted, stop_casted, .{ .allocator = allocator });
    errdefer ops.deinit(&diff, .{ .allocator = allocator });

    var len_T: T = try ops.div(diff, step_casted, .{ .allocator = allocator });
    errdefer ops.deinit(&len_T, .{ .allocator = allocator });
    try ops.abs_(&len_T, .{ .allocator = allocator });
    try ops.ceil_(&len_T, .{ .allocator = allocator });
    const len: usize = try cast(usize, len_T, .{ .allocator = allocator });

    if (len == 0) {
        ops.deinit(&start_casted, .{ .allocator = allocator });
        ops.deinit(&stop_casted, .{ .allocator = allocator });
        ops.deinit(&step_casted, .{ .allocator = allocator });
        ops.deinit(&diff, .{ .allocator = allocator });
        ops.deinit(&len_T, .{ .allocator = allocator });

        return array.Error.InvalidRange;
    }

    var arr: Array(T) = try init(allocator, T, &.{len}, .rowMajor);
    errdefer arr.deinit(allocator);
    arr.flags.writeable = writeable;
    switch (len) {
        1 => {
            arr.data[0] = start_casted;

            ops.deinit(&stop_casted, .{ .allocator = allocator });
            ops.deinit(&step_casted, .{ .allocator = allocator });
            ops.deinit(&diff, .{ .allocator = allocator });
            ops.deinit(&len_T, .{ .allocator = allocator });

            return arr;
        },
        2 => {
            arr.data[0] = start_casted;
            try ops.add_(&step_casted, arr.data[0], .{ .allocator = allocator });
            arr.data[1] = step_casted;

            ops.deinit(&stop_casted, .{ .allocator = allocator });
            ops.deinit(&diff, .{ .allocator = allocator });
            ops.deinit(&len_T, .{ .allocator = allocator });

            return arr;
        },
        3 => {
            arr.data[0] = start_casted;
            try ops.add_to(&stop_casted, arr.data[0], step_casted, .{ .allocator = allocator });
            arr.data[1] = stop_casted;
            try ops.add_to(&step_casted, arr.data[1], step_casted, .{ .allocator = allocator });
            arr.data[2] = step_casted;

            ops.deinit(&diff, .{ .allocator = allocator });
            ops.deinit(&len_T, .{ .allocator = allocator });

            return arr;
        },
        4 => {
            arr.data[0] = start_casted;
            try ops.add_to(&diff, arr.data[0], step_casted, .{ .allocator = allocator });
            arr.data[1] = diff;
            try ops.add_to(&stop_casted, arr.data[1], step_casted, .{ .allocator = allocator });
            arr.data[2] = stop_casted;
            try ops.add_to(&step_casted, arr.data[2], step_casted, .{ .allocator = allocator });
            arr.data[3] = step_casted;

            ops.deinit(&len_T, .{ .allocator = allocator });

            return arr;
        },
        else => {
            arr.data[0] = start_casted;
            try ops.add_to(&len_T, arr.data[0], step_casted, .{ .allocator = allocator });
            arr.data[1] = len_T;
            try ops.add_to(&diff, arr.data[1], step_casted, .{ .allocator = allocator });
            arr.data[2] = diff;
            try ops.add_to(&stop_casted, arr.data[2], step_casted, .{ .allocator = allocator });
            arr.data[3] = stop_casted;
        },
    }

    var j: usize = 4;
    errdefer cleanup(T, arr.data[4..j], allocator);

    for (4..len - 1) |i| {
        arr.data[i] = try ops.add(arr.data[i - 1], step_casted, .{ .allocator = allocator });
        j += 1;
    }

    try ops.add_(&step_casted, arr.data[len - 2], .{ .allocator = allocator });
    arr.data[len - 1] = step_casted;

    return arr;
}

pub inline fn linspace(
    allocator: std.mem.Allocator,
    comptime T: type,
    start: anytype,
    stop: anytype,
    num: usize,
    writeable: bool,
    endpoint: bool,
) !Array(T) {
    _ = types.numericType(@TypeOf(start));
    _ = types.numericType(@TypeOf(stop));

    var arr: Array(T) = try init(allocator, T, &.{num}, .rowMajor);
    errdefer arr.deinit(allocator);
    arr.flags.writeable = writeable;

    if (num == 1) {
        arr.data[0] = try cast(T, start, .{ .allocator = allocator, .copy = true });

        return arr;
    } else if (num == 2 and endpoint) {
        arr.data[0] = try cast(T, start, .{ .allocator = allocator, .copy = true });
        errdefer ops.deinit(&arr.data[0], .{ .allocator = allocator });
        arr.data[1] = try cast(T, stop, .{ .allocator = allocator, .copy = true });

        return arr;
    } else if (num == 2 and !endpoint) {
        arr.data[0] = try cast(T, start, .{ .allocator = allocator, .copy = true });
        errdefer ops.deinit(&arr.data[0], .{ .allocator = allocator });
        arr.data[1] = try cast(T, stop, .{ .allocator = allocator, .copy = true });
        errdefer ops.deinit(&arr.data[1], .{ .allocator = allocator });
        try ops.add_(&arr.data[1], arr.data[0], .{ .allocator = allocator });
        try ops.div_(&arr.data[1], 2, .{ .allocator = allocator });

        return arr;
    }

    var start_casted: T = try cast(T, start, .{ .allocator = allocator, .copy = true });
    errdefer ops.deinit(&start_casted, .{ .allocator = allocator });
    var stop_casted: T = try cast(T, stop, .{ .allocator = allocator, .copy = true });
    errdefer ops.deinit(&stop_casted, .{ .allocator = allocator });
    var step: T = try ops.sub(stop_casted, start_casted, .{ .allocator = allocator });
    errdefer ops.deinit(&step, .{ .allocator = allocator });
    if (endpoint) {
        try ops.div_(&step, num - 1, .{ .allocator = allocator });
    } else {
        try ops.div_(&step, num, .{ .allocator = allocator });
    }

    if (num == 3 and endpoint) {
        arr.data[0] = start_casted;
        try ops.add_(&step, arr.data[0], .{ .allocator = allocator });
        arr.data[1] = step;
        arr.data[2] = stop_casted;
    } else if (num == 3 and !endpoint) {
        arr.data[0] = start_casted;
        try ops.add_(&step, arr.data[0], .{ .allocator = allocator });
        arr.data[1] = step;
        try ops.add_to(&stop_casted, arr.data[1], step, .{ .allocator = allocator });
        arr.data[2] = stop_casted;

        return arr;
    }

    arr.data[0] = start_casted;

    var j: usize = 1;
    errdefer cleanup(T, allocator, arr.data[1..j]);
    for (1..num - 2) |i| {
        arr.data[i] = try ops.add(arr.data[i - 1], step, .{ .allocator = allocator });

        j += 1;
    }

    if (endpoint) {
        try ops.add_to(&step, arr.data[num - 3], step, .{ .allocator = allocator });
        arr.data[num - 2] = step;
        arr.data[num - 1] = stop_casted;
    } else {
        try ops.sub_(&stop_casted, step, .{ .allocator = allocator });
        try ops.add_(&step, arr.data[num - 3], .{ .allocator = allocator });
        arr.data[num - 2] = step;
        arr.data[num - 1] = stop_casted;
    }

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
    var new_shape: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    var new_strides: [array.maxDimensions]isize = .{0} ** array.maxDimensions;
    if (shape.len > 0) {
        for (0..shape.len) |i| {
            const idx: usize = if (arr.flags.order == .rowMajor) shape.len - i - 1 else i;

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
        .base = if (arr.flags.ownsData) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided, // Should be dense?
            .ownsData = false,
            .writeable = arr.flags.writeable,
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
        .shape = .{arr.size} ++ .{0} ** (array.maxDimensions - 1),
        .size = arr.size,
        .base = if (arr.flags.ownsData) arr else arr.base,
        .flags = .{
            .order = arr.flags.order,
            .storage = .dense,
            .ownsData = false,
            .writeable = arr.flags.writeable,
        },
        .metadata = .{ .dense = .{
            .strides = .{1} ++ .{0} ** (array.maxDimensions - 1),
        } },
    };
}

pub fn transpose(
    comptime T: type,
    self: *const Array(T),
    axes: []const usize,
) !Array(T) {
    var new_shape: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    var new_strides: [array.maxDimensions]isize = .{0} ** array.maxDimensions;
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
        .base = if (self.flags.ownsData) self else self.base,
        .flags = .{
            .order = self.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided,
            .ownsData = false,
            .writeable = self.flags.writeable,
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
    var shape: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    var strides: [array.maxDimensions]isize = .{0} ** array.maxDimensions;
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
        if (range.start != int.max(usize) and range.start == range.stop) {
            return array.Error.InvalidRange;
        } else if (range.step > 0) {
            if (range.start != int.max(usize) and range.start >= arr.shape[i] or
                (range.stop != int.max(usize) and range.stop > arr.shape[i]))
                return array.Error.RangeOutOfBounds;
        } else if (range.step < 0) {
            if ((range.stop != int.max(usize) and range.stop >= arr.shape[i]) or
                (range.start != int.max(usize) and range.start > arr.shape[i]))
                return array.Error.RangeOutOfBounds;
        }

        var len_adjustment: usize = 0;
        if (range.step > 0) {
            if (range.start == int.max(usize)) {
                range.start = 0;
            }

            if (range.stop == int.max(usize)) {
                range.stop = arr.shape[i];
            }
        } else if (range.step < 0) {
            if (range.start == int.max(usize)) {
                range.start = arr.shape[i] - 1;
            }

            if (range.stop == int.max(usize)) {
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
        .base = if (arr.flags.ownsData) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided,
            .ownsData = false,
            .writeable = arr.flags.writeable,
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
    var new_shape: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    var strides: [array.maxDimensions]isize = .{0} ** array.maxDimensions;
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
        .base = if (arr.flags.ownsData) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided,
            .ownsData = false,
            .writeable = arr.flags.writeable,
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
    allocator: std.mem.Allocator,
    comptime T: type,
    arr: Array(T),
    comptime op: anytype,
    writeable: bool,
) !Array(ReturnType1(op, T)) {
    var newarr: Array(ReturnType1(op, T)) = try init(allocator, ReturnType1(op, T), arr.shape[0..arr.ndim], arr.flags.order);
    errdefer newarr.deinit(allocator);
    newarr.flags.writeable = writeable;

    var j: usize = 0;
    errdefer cleanup(ReturnType1(op, T), allocator, newarr.data[0..j]);

    const opinfo = @typeInfo(@TypeOf(op));
    for (0..newarr.size) |i| {
        if (opinfo.@"fn".params.len == 1) {
            newarr.data[i] = op(arr.data[i]);
        } else if (opinfo.@"fn".params.len == 2) {
            newarr.data[i] = try op(arr.data[i], .{ .allocator = allocator });
        }

        j += 1;
    }

    return newarr;
}

pub inline fn apply1_(
    comptime T: type,
    arr: *Array(T),
    comptime op_: anytype,
    allocator: ?std.mem.Allocator,
) !void {
    const opinfo = @typeInfo(@TypeOf(op_));
    for (0..arr.size) |i| {
        if (opinfo.@"fn".params.len == 1) {
            op_(&arr.data[i]);
        } else if (opinfo.@"fn".params.len == 2) {
            try op_(&arr.data[i], .{ .allocator = allocator });
        }
    }

    return;
}

pub fn apply1_to(
    comptime O: type,
    o: anytype,
    comptime X: type,
    x: anytype,
    comptime op_to: anytype,
    allocator: ?std.mem.Allocator,
) !void {
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        // Only run the function once if x is a scalar
        const opinfo = @typeInfo(@TypeOf(op_to));
        if (opinfo.@"fn".params.len == 2) {
            op_to(&o.data[0], x);
        } else if (opinfo.@"fn".params.len == 3) {
            try op_to(&o.data[0], x, .{ .allocator = allocator });
        }

        for (1..o.size) |i| {
            try ops.set(&o.data[i], x, .{ .allocator = allocator });
        }

        return;
    }

    var xx: Array(X) = undefined;
    if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
        if (o.flags.order == x.flags.order) {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_to));
            for (0..o.size) |i| {
                if (opinfo.@"fn".params.len == 2) {
                    op_to(&o.data[i], x.data[i]);
                } else if (opinfo.@"fn".params.len == 3) {
                    try op_to(&o.data[i], x.data[i], .{ .allocator = allocator });
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

    const iterationOrder: array.IterationOrder = if (o.flags.order == .rowMajor) .rightToLeft else .leftToRight;
    const axis: usize = if (o.flags.order == .rowMajor) o.ndim - 1 else 0;
    var itero: array.Iterator(O) = .init(o);
    var iterx = array.Iterator(X).init(&xx);
    const opinfo = @typeInfo(@TypeOf(op_to));
    for (0..o.size) |_| {
        if (opinfo.@"fn".params.len == 2) {
            op_to(&o.data[itero.index], xx.data[iterx.index]);
        } else if (opinfo.@"fn".params.len == 3) {
            try op_to(&o.data[itero.index], xx.data[iterx.index], .{ .allocator = allocator });
        }

        _ = itero.nextAO(axis, iterationOrder);
        _ = iterx.nextAO(axis, iterationOrder);
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
    options: struct {
        order: ?array.Order = null,
        writeable: bool = true,
    },
) !Array(ReturnType2(op, X, Y)) {
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        var newarr: Array(ReturnType2(op, X, Y)) = try init(allocator, ReturnType2(op, X, Y), y.shape[0..y.ndim], options.order orelse y.flags.order);
        errdefer newarr.deinit(allocator);
        newarr.flags.writeable = options.writeable;

        var j: usize = 0;
        if (newarr.flags.order == y.flags.order) {
            errdefer cleanup(ReturnType2(op, X, Y), allocator, newarr.data[0..j]);

            const opinfo = @typeInfo(@TypeOf(op));
            for (0..newarr.size) |i| {
                if (opinfo.@"fn".params.len == 2) {
                    newarr.data[i] = op(x, y.data[i]);
                } else if (opinfo.@"fn".params.len == 3) {
                    newarr.data[i] = try op(x, y.data[i], .{ .allocator = allocator });
                }

                j += 1;
            }
        } else {
            const iterationOrder: array.IterationOrder = if (newarr.flags.order == .rowMajor) .rightToLeft else .leftToRight;
            const axis: usize = if (newarr.flags.order == .rowMajor) newarr.ndim - 1 else 0;
            errdefer strided.cleanup(ReturnType2(op, X, Y), allocator, &newarr, j, iterationOrder);
            var itero: array.Iterator(ReturnType2(op, X, Y)) = .init(&newarr);
            var itery = array.Iterator(Y).init(&y);
            const opinfo = @typeInfo(@TypeOf(op));
            for (0..newarr.size) |_| {
                if (opinfo.@"fn".params.len == 2) {
                    newarr.data[itero.index] = op(x, y.data[itery.index]);
                } else if (opinfo.@"fn".params.len == 3) {
                    newarr.data[itero.index] = try op(x, y.data[itery.index], .{ .allocator = allocator });
                }

                j += 1;
                _ = itero.nextAO(axis, iterationOrder);
                _ = itery.nextAO(axis, iterationOrder);
            }
        }

        return newarr;
    } else if (comptime !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        var newarr: Array(ReturnType2(op, X, Y)) = try init(allocator, ReturnType2(op, X, Y), x.shape[0..x.ndim], options.order orelse x.flags.order);
        errdefer newarr.deinit(allocator);
        newarr.flags.writeable = options.writeable;

        var j: usize = 0;
        if (newarr.flags.order == x.flags.order) {
            errdefer cleanup(ReturnType2(op, X, Y), allocator, newarr.data[0..j]);

            const opinfo = @typeInfo(@TypeOf(op));
            for (0..newarr.size) |i| {
                if (opinfo.@"fn".params.len == 2) {
                    newarr.data[i] = op(x.data[i], y);
                } else if (opinfo.@"fn".params.len == 3) {
                    newarr.data[i] = try op(x.data[i], y, .{ .allocator = allocator });
                }

                j += 1;
            }
        } else {
            const iterationOrder: array.IterationOrder = if (newarr.flags.order == .rowMajor) .rightToLeft else .leftToRight;
            const axis: usize = if (newarr.flags.order == .rowMajor) newarr.ndim - 1 else 0;
            errdefer strided.cleanup(ReturnType2(op, X, Y), allocator, &newarr, j, iterationOrder);
            var itero: array.Iterator(ReturnType2(op, X, Y)) = .init(&newarr);
            var iterx = array.Iterator(X).init(&x);
            const opinfo = @typeInfo(@TypeOf(op));
            for (0..newarr.size) |_| {
                if (opinfo.@"fn".params.len == 2) {
                    newarr.data[itero.index] = op(x.data[iterx.index], y);
                } else if (opinfo.@"fn".params.len == 3) {
                    newarr.data[itero.index] = try op(x.data[iterx.index], y, .{ .allocator = allocator });
                }

                j += 1;
                _ = itero.nextAO(axis, iterationOrder);
                _ = iterx.nextAO(axis, iterationOrder);
            }
        }

        return newarr;
    }

    const order: array.Order = options.order orelse if (x.flags.order == y.flags.order) x.flags.order else .rowMajor;
    var xx: Array(X) = undefined;
    var yy: Array(Y) = undefined;
    if (std.mem.eql(usize, x.shape[0..x.ndim], y.shape[0..y.ndim])) {
        if (x.flags.order == y.flags.order and (x.flags.order == order)) {
            // Trivial loop
            var newarr: Array(ReturnType2(op, X, Y)) = try init(allocator, ReturnType2(op, X, Y), x.shape[0..x.ndim], order);
            errdefer newarr.deinit(allocator);
            newarr.flags.writeable = options.writeable;

            var j: usize = 0;
            errdefer cleanup(ReturnType2(op, X, Y), allocator, newarr.data[0..j]);
            const opinfo = @typeInfo(@TypeOf(op));
            for (0..newarr.size) |i| {
                if (opinfo.@"fn".params.len == 2) {
                    newarr.data[i] = op(x.data[i], y.data[i]);
                } else if (opinfo.@"fn".params.len == 3) {
                    newarr.data[i] = try op(x.data[i], y.data[i], .{ .allocator = allocator });
                }

                j += 1;
            }

            return newarr;
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

    var newarr: Array(ReturnType2(op, X, Y)) = try init(allocator, ReturnType2(op, X, Y), xx.shape[0..xx.ndim], order);
    errdefer newarr.deinit(allocator);
    newarr.flags.writeable = options.writeable;

    var j: usize = 0;
    const iterationOrder: array.IterationOrder = if (newarr.flags.order == .rowMajor) .rightToLeft else .leftToRight;
    const axis: usize = if (newarr.flags.order == .rowMajor) newarr.ndim - 1 else 0;
    errdefer strided.cleanup(ReturnType2(op, X, Y), allocator, &newarr, j, iterationOrder);
    var itero: array.Iterator(ReturnType2(op, X, Y)) = .init(&newarr);
    var iterx = array.Iterator(X).init(&xx);
    var itery = array.Iterator(Y).init(&yy);
    const opinfo = @typeInfo(@TypeOf(op));
    for (0..newarr.size) |_| {
        if (opinfo.@"fn".params.len == 2) {
            newarr.data[itero.index] = op(x.data[iterx.index], y.data[itery.index]);
        } else if (opinfo.@"fn".params.len == 3) {
            newarr.data[itero.index] = try op(x.data[iterx.index], y.data[itery.index], .{ .allocator = allocator });
        }

        j += 1;
        _ = itero.nextAO(axis, iterationOrder);
        _ = iterx.nextAO(axis, iterationOrder);
        _ = itery.nextAO(axis, iterationOrder);
    }

    return newarr;
}

pub fn apply2_(
    comptime O: type,
    o: anytype,
    comptime Y: type,
    y: anytype,
    comptime op_: anytype,
    allocator: ?std.mem.Allocator,
) !void {
    if (comptime !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        const opinfo = @typeInfo(@TypeOf(op_));
        for (0..o.size) |i| {
            if (opinfo.@"fn".params.len == 2) {
                op_(&o.data[i], y);
            } else if (opinfo.@"fn".params.len == 3) {
                try op_(&o.data[i], y, .{ .allocator = allocator });
            }
        }

        return;
    }

    var yy: Array(Y) = undefined;
    if (std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim])) {
        if (o.flags.order == y.flags.order) {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (opinfo.@"fn".params.len == 2) {
                    op_(&o.data[i], y.data[i]);
                } else if (opinfo.@"fn".params.len == 3) {
                    try op_(&o.data[i], y.data[i], .{ .allocator = allocator });
                }
            }

            return;
        } else {
            // Different order, but same shape
            yy = y;
        }
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], y.shape[0..y.ndim] });
        if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
            return array.Error.NotBroadcastable;
        }

        yy = try broadcast(Y, &y, bct.shape[0..bct.ndim]);
    }

    const iterationOrder: array.IterationOrder = if (o.flags.order == .rowMajor) .rightToLeft else .leftToRight;
    const axis: usize = if (o.flags.order == .rowMajor) o.ndim - 1 else 0;
    var itero: array.Iterator(O) = .init(o);
    var itery = array.Iterator(Y).init(&yy);
    const opinfo = @typeInfo(@TypeOf(op_));
    for (0..o.size) |_| {
        if (opinfo.@"fn".params.len == 2) {
            op_(&o.data[itero.index], y.data[itery.index]);
        } else if (opinfo.@"fn".params.len == 3) {
            try op_(&o.data[itero.index], y.data[itery.index], .{ .allocator = allocator });
        }

        _ = itero.nextAO(axis, iterationOrder);
        _ = itery.nextAO(axis, iterationOrder);
    }

    return;
}

pub fn apply2_to(
    comptime O: type,
    o: anytype,
    comptime X: type,
    x: anytype,
    comptime Y: type,
    y: anytype,
    comptime op_to: anytype,
    allocator: ?std.mem.Allocator,
) !void {
    if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x)) and
        !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x)))
    {
        const opinfo = @typeInfo(@TypeOf(op_to));
        if (opinfo.@"fn".params.len == 3) {
            op_to(&o.data[0], x, y);
        } else if (opinfo.@"fn".params.len == 4) {
            try op_to(&o.data[o], x, y, .{ .allocator = allocator });
        }

        for (1..o.size) |i| {
            try ops.set(&o.data[i], o.data[0], .{ .allocator = allocator });
        }

        return;
    } else if (comptime !types.isArray(@TypeOf(x)) and !types.isSlice(@TypeOf(x))) {
        if (std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim]) and
            o.flags.order == y.flags.order)
        {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_to));
            for (0..o.size) |i| {
                if (opinfo.@"fn".params.len == 3) {
                    op_to(&o.data[i], x, y.data[i]);
                } else if (opinfo.@"fn".params.len == 4) {
                    try op_to(&o.data[i], x, y.data[i], .{ .allocator = allocator });
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

        const iterationOrder: array.IterationOrder = if (o.flags.order == .rowMajor) .rightToLeft else .leftToRight;
        const axis: usize = if (o.flags.order == .rowMajor) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        var itery = array.Iterator(Y).init(&yy);
        const opinfo = @typeInfo(@TypeOf(op_to));
        for (0..o.size) |_| {
            if (opinfo.@"fn".params.len == 3) {
                op_to(&o.data[itero.index], x, yy.data[itery.index]);
            } else if (opinfo.@"fn".params.len == 4) {
                try op_to(&o.data[itero.index], x, yy.data[itery.index], .{ .allocator = allocator });
            }

            _ = itero.nextAO(axis, iterationOrder);
            _ = itery.nextAO(axis, iterationOrder);
        }

        return;
    } else if (comptime !types.isArray(@TypeOf(y)) and !types.isSlice(@TypeOf(y))) {
        if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
            o.flags.order == x.flags.order)
        {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_to));
            for (0..o.size) |i| {
                if (opinfo.@"fn".params.len == 3) {
                    op_to(&o.data[i], x.data[i], y);
                } else if (opinfo.@"fn".params.len == 4) {
                    try op_to(&o.data[i], x.data[i], y, .{ .allocator = allocator });
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

        const iterationOrder: array.IterationOrder = if (o.flags.order == .rowMajor) .rightToLeft else .leftToRight;
        const axis: usize = if (o.flags.order == .rowMajor) o.ndim - 1 else 0;
        var itero: array.Iterator(O) = .init(o);
        var iterx = array.Iterator(X).init(&xx);
        const opinfo = @typeInfo(@TypeOf(op_to));
        for (0..o.size) |_| {
            if (opinfo.@"fn".params.len == 3) {
                op_to(&o.data[itero.index], xx.data[iterx.index], y);
            } else if (opinfo.@"fn".params.len == 4) {
                try op_to(&o.data[itero.index], xx.data[iterx.index], y, .{ .allocator = allocator });
            }

            _ = itero.nextAO(axis, iterationOrder);
            _ = iterx.nextAO(axis, iterationOrder);
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
            const opinfo = @typeInfo(@TypeOf(op_to));
            for (0..o.size) |i| {
                if (opinfo.@"fn".params.len == 3) {
                    op_to(&o.data[i], x.data[i], y.data[i]);
                } else if (opinfo.@"fn".params.len == 4) {
                    try op_to(&o.data[i], x.data[i], y.data[i], .{ .allocator = allocator });
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

    const iterationOrder: array.IterationOrder = if (o.flags.order == .rowMajor) .rightToLeft else .leftToRight;
    const axis: usize = if (o.flags.order == .rowMajor) o.ndim - 1 else 0;
    var itero: array.Iterator(O) = .init(o);
    var iterx = array.Iterator(X).init(&xx);
    var itery = array.Iterator(Y).init(&yy);
    const opinfo = @typeInfo(@TypeOf(op_to));
    for (0..o.size) |_| {
        if (opinfo.@"fn".params.len == 3) {
            op_to(&o.data[itero.index], xx.data[iterx.index], yy.data[itery.index]);
        } else if (opinfo.@"fn".params.len == 4) {
            try op_to(&o.data[itero.index], xx.data[iterx.index], yy.data[itery.index], .{ .allocator = allocator });
        }

        _ = itero.nextAO(axis, iterationOrder);
        _ = iterx.nextAO(axis, iterationOrder);
        _ = itery.nextAO(axis, iterationOrder);
    }

    return;
}
