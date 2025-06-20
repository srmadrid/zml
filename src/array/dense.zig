const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
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

inline fn cleanup(
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

    arr.data[index(arr.ndim, arr.metadata.dense.strides, position)] = value;
}

pub inline fn get(comptime T: type, arr: Array(T), position: []const usize) !*T {
    try checkPosition(arr.ndim, arr.shape, position);

    return &arr.data[index(arr.ndim, arr.metadata.dense.strides, position)];
}

pub inline fn reshape(comptime T: type, arr: *const Array(T), shape: []const usize) !Array(T) {
    var new_size: usize = 1;
    var new_shape: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    var new_strides: [array.maxDimensions]isize = .{0} ** array.maxDimensions;
    if (shape.len > 0) {
        for (0..shape.len) |i| {
            const idx: usize = if (arr.order == .rowMajor) shape.len - i - 1 else i;

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
