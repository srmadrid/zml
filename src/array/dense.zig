const std = @import("std");

const types = @import("../types.zig");
const scast = types.scast;
const cast = types.cast;
const needsAllocator = types.needsAllocator;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");
const Array = array.Array;
const Range = array.Range;

pub inline fn index(ndim: usize, strides: [array.maxDimensions]usize, position: []const usize) usize {
    var idx: usize = 0;
    for (0..ndim) |i| {
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

            strides[idx] = size;
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
    const value_casted: T = try cast(T, value, .{ .allocator = allocator });
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
    const positive_step: bool = try ops.gt(step, 0, .{});
    if (try ops.eq(step, 0, .{}) or
        (try ops.lt(stop, start, .{}) and positive_step) or
        (try ops.gt(stop, start, .{}) and !positive_step))
    {
        return array.Error.InvalidRange;
    }

    const a: if (needsAllocator(T)) std.mem.Allocator else void = comptime if (needsAllocator(T)) allocator else {};
    var start_casted: T = try cast(T, start, .{ .allocator = a, .copy = true });
    errdefer ops.deinit(&start_casted, .{ .allocator = a });
    var stop_casted: T = try cast(T, stop, .{ .allocator = a, .copy = true });
    errdefer ops.deinit(&stop_casted, .{ .allocator = a });
    var step_casted: T = try cast(T, step, .{ .allocator = a, .copy = true });
    errdefer ops.deinit(&step_casted, .{ .allocator = a });

    var diff: T = if (positive_step)
        try ops.sub(stop_casted, start_casted, .{ .allocator = a })
    else
        try ops.sub(start_casted, stop_casted, .{ .allocator = a });
    errdefer ops.deinit(&diff, .{ .allocator = a });

    std.debug.print("diff: {d}\n", .{diff});

    var len_T: T = try ops.div(diff, step_casted, .{ .allocator = a });
    try ops.abs_(&len_T, len_T, .{ .allocator = a });
    try ops.ceil_(&len_T, len_T, .{ .allocator = a });
    const len: usize = try cast(usize, len_T, .{ .allocator = a });

    if (len == 0) {
        ops.deinit(&start_casted, .{ .allocator = a });
        ops.deinit(&stop_casted, .{ .allocator = a });
        ops.deinit(&step_casted, .{ .allocator = a });
        ops.deinit(&diff, .{ .allocator = a });
        ops.deinit(&len_T, .{ .allocator = a });

        return array.Error.InvalidRange;
    }

    var arr: Array(T) = try init(allocator, T, &.{len}, .rowMajor);
    arr.flags.writeable = writeable;
    switch (len) {
        1 => {
            ops.deinit(&stop_casted, .{ .allocator = a });
            ops.deinit(&step_casted, .{ .allocator = a });
            ops.deinit(&diff, .{ .allocator = a });
            ops.deinit(&len_T, .{ .allocator = a });

            arr.data[0] = start_casted;

            return arr;
        },
        2 => {
            ops.deinit(&stop_casted, .{ .allocator = a });
            ops.deinit(&diff, .{ .allocator = a });
            ops.deinit(&len_T, .{ .allocator = a });

            arr.data[0] = start_casted;
            try ops.add_(&step_casted, arr.data[0], step_casted, .{ .allocator = a });
            arr.data[1] = step_casted;

            return arr;
        },
        3 => {
            ops.deinit(&diff, .{ .allocator = a });
            ops.deinit(&len_T, .{ .allocator = a });

            arr.data[0] = start_casted;
            try ops.add_(&stop_casted, arr.data[0], step_casted, .{ .allocator = a });
            arr.data[1] = stop_casted;
            try ops.add_(&step_casted, arr.data[1], step_casted, .{ .allocator = a });
            arr.data[2] = step_casted;

            return arr;
        },
        4 => {
            ops.deinit(&len_T, .{ .allocator = a });

            arr.data[0] = start_casted;
            try ops.add_(&diff, arr.data[0], step_casted, .{ .allocator = a });
            arr.data[1] = diff;
            try ops.add_(&stop_casted, arr.data[1], step_casted, .{ .allocator = a });
            arr.data[2] = stop_casted;
            try ops.add_(&step_casted, arr.data[2], step_casted, .{ .allocator = a });
            arr.data[3] = step_casted;

            return arr;
        },
        else => {
            arr.data[0] = start_casted;
            try ops.add_(&len_T, arr.data[0], step_casted, .{ .allocator = a });
            arr.data[1] = len_T;
            try ops.add_(&diff, arr.data[1], step_casted, .{ .allocator = a });
            arr.data[2] = diff;
            try ops.add_(&stop_casted, arr.data[2], step_casted, .{ .allocator = a });
            arr.data[3] = stop_casted;
        },
    }

    for (4..len - 1) |i| {
        arr.data[i] = try ops.add(arr.data[i - 1], step_casted, .{ .allocator = a });
    }

    try ops.add_(&step_casted, arr.data[len - 2], step_casted, .{ .allocator = a });
    arr.data[len - 1] = step_casted;

    return arr;
}

pub inline fn set(comptime T: type, arr: *Array(T), position: []const usize, value: anytype) !void {
    try checkPosition(arr.ndim, arr.shape, position);

    arr.data[index(arr.ndim, arr.metadata.dense.strides, position)] = try cast(T, value, .{ .allocator = arr.flags.writeable });
}

pub inline fn get(comptime T: type, arr: *const Array(T), position: []const usize) !*T {
    try checkPosition(arr.ndim, arr.shape, position);

    return &arr.data[index(arr.ndim, arr.metadata.dense.strides, position)];
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
