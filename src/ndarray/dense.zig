const std = @import("std");

const types = @import("../types.zig");
const scast = types.scast;
const cast = types.cast;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const ndarray = @import("../ndarray.zig");
const NDArray = ndarray.NDArray;
const Range = ndarray.Range;

pub inline fn index(ndim: usize, strides: [ndarray.maxDimensions]usize, position: []const usize) usize {
    var idx: usize = 0;
    for (0..ndim) |i| {
        idx += position[i] * strides[i];
    }

    return idx;
}

pub inline fn checkPosition(ndim: usize, shape: [ndarray.maxDimensions]usize, position: []const usize) !void {
    if (position.len > ndim) {
        return ndarray.Error.DimensionMismatch;
    }

    for (0..position.len) |i| {
        if (position[i] >= shape[i]) {
            return ndarray.Error.PositionOutOfBounds;
        }
    }
}

pub inline fn init(
    allocator: std.mem.Allocator,
    comptime T: type,
    shape: []const usize,
    order: ndarray.Order,
) !NDArray(T) {
    var size: usize = 1;
    var shapes: [ndarray.maxDimensions]usize = .{0} ** ndarray.maxDimensions;
    var strides: [ndarray.maxDimensions]usize = .{0} ** ndarray.maxDimensions;
    if (shape.len > 0) {
        for (0..shape.len) |i| {
            const idx: usize = if (order == .rowMajor) shape.len - i - 1 else i;

            strides[idx] = size;
            size *= shape[idx];

            shapes[i] = shape[i];
        }
    }

    return NDArray(T){
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
    order: ndarray.Order,
) !NDArray(T) {
    var array = try init(allocator, T, shape, order);
    for (0..array.size) |i| {
        array.data[i] = cast(T, value, .{ .allocator = allocator });
    }

    return array;
}

pub inline fn arange(
    allocator: std.mem.Allocator,
    comptime T: type,
    start: T,
    stop: T,
    step: T,
) !NDArray(T) {
    _ = allocator; // Unused in this function, but required for the signature
    if (step == 0 or
        (ops.lt(stop, start) and ops.gt(step, 0)) or
        (ops.gt(stop, start) and ops.lt(step, 0)))
    {
        return ndarray.Error.InvalidRange;
    }
}

pub inline fn set(comptime T: type, array: *NDArray(T), position: []const usize, value: T) !void {
    try checkPosition(array.ndim, array.shape, position);

    array.data[index(array.ndim, array.metadata.dense.strides, position)] = value;
}

pub inline fn get(comptime T: type, array: *const NDArray(T), position: []const usize) !*T {
    try checkPosition(array.ndim, array.shape, position);

    return &array.data[index(array.ndim, array.metadata.dense.strides, position)];
}

pub inline fn slice(comptime T: type, array: *const NDArray(T), ranges: []const Range) !NDArray(T) {
    var ndim: usize = array.ndim;
    var size: usize = 1;
    var shape: [ndarray.maxDimensions]usize = .{0} ** ndarray.maxDimensions;
    var strides: [ndarray.maxDimensions]isize = .{0} ** ndarray.maxDimensions;
    var offset: usize = if (array.flags.storage == .dense) 0 else array.metadata.strided.offset;

    var i: usize = 0;
    var j: usize = 0;
    while (i < array.ndim) {
        const stride: isize = scast(isize, array.metadata.dense.strides[i]);

        if (i >= ranges.len) {
            shape[j] = array.shape[i];
            strides[j] = stride;
            size *= array.shape[i];
            j += 1;
            i += 1;
            continue;
        }

        var range: Range = ranges[i];
        if (range.start != int.max(usize) and range.start == range.stop) {
            return ndarray.Error.InvalidRange;
        } else if (range.step > 0) {
            if (range.start != int.max(usize) and range.start >= array.shape[i] or
                (range.stop != int.max(usize) and range.stop > array.shape[i]))
                return ndarray.Error.RangeOutOfBounds;
        } else if (range.step < 0) {
            if ((range.stop != int.max(usize) and range.stop >= array.shape[i]) or
                (range.start != int.max(usize) and range.start > array.shape[i]))
                return ndarray.Error.RangeOutOfBounds;
        }

        var len_adjustment: usize = 0;
        if (range.step > 0) {
            if (range.start == int.max(usize)) {
                range.start = 0;
            }

            if (range.stop == int.max(usize)) {
                range.stop = array.shape[i];
            }
        } else if (range.step < 0) {
            if (range.start == int.max(usize)) {
                range.start = array.shape[i] - 1;
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

    return NDArray(T){
        .data = array.data,
        .ndim = ndim,
        .shape = shape,
        .size = size,
        .base = if (array.flags.ownsData) array else array.base,
        .flags = .{
            .order = array.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .storage = .strided,
            .ownsData = false,
            .writeable = array.flags.writeable,
        },
        .metadata = .{ .strided = .{
            .strides = strides,
            .offset = offset,
        } },
    };
}
