const std = @import("std");

const types = @import("../types.zig");
const scast = types.scast;
const cast = types.cast;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const ndarray = @import("../ndarray.zig");
const NDArray = ndarray.NDArray;
const Range = ndarray.Range;

pub inline fn index(ndim: usize, strides: [ndarray.maxDimensions]isize, offset: usize, position: []const usize) usize {
    var idx: usize = offset;
    for (0..ndim) |i| {
        const stride: isize = strides[i];
        if (stride < 0) {
            idx -= position[i] * scast(usize, int.abs(stride));
        } else {
            idx += position[i] * scast(usize, stride);
        }
    }

    return idx;
}

pub inline fn checkPosition(ndim: usize, shape: [ndarray.maxDimensions]usize, position: []const usize) !void {
    if (position.len != ndim) {
        return ndarray.Error.DimensionMismatch;
    }

    for (0..position.len) |i| {
        if (position[i] >= shape[i]) {
            return ndarray.Error.PositionOutOfBounds;
        }
    }
}

pub inline fn set(comptime T: type, array: *NDArray(T), position: []const usize, value: T) !void {
    try checkPosition(array.ndim, array.shape, position);

    array.data[index(array.ndim, array.metadata.strided.strides, array.metadata.strided.offset, position)] = value;
}

pub fn get(comptime T: type, array: *const NDArray(T), position: []const usize) !*T {
    try checkPosition(array.ndim, array.shape, position);

    return &array.data[index(array.ndim, array.metadata.strided.strides, array.metadata.strided.offset, position)];
}

pub inline fn slice(comptime T: type, array: *const NDArray(T), ranges: []const Range) !NDArray(T) {
    var ndim: usize = array.ndim;
    var size: usize = 1;
    var shapes: [ndarray.maxDimensions]usize = .{0} ** ndarray.maxDimensions;
    var strides: [ndarray.maxDimensions]isize = .{0} ** ndarray.maxDimensions;
    var offset: usize = array.metadata.strided.offset;

    var i: usize = 0;
    var j: usize = 0;
    while (i < array.ndim) {
        const stride: isize = array.metadata.strided.strides[i];

        if (i >= ranges.len) {
            shapes[j] = array.shape[i];
            strides[j] = stride;
            size *= array.shape[i];
            j += 1;
            i += 1;
            continue;
        }

        var range: Range = ranges[i];
        if (range.step > 0 and range.stop != int.max(usize)) {
            if (range.start >= array.shape[i] or range.stop > array.shape[i]) {
                return ndarray.Error.RangeOutOfBounds;
            }
        } else if (range.step < 0 and range.start != int.max(usize)) {
            if (range.stop >= array.shape[i] or range.start > array.shape[i]) {
                return ndarray.Error.RangeOutOfBounds;
            }
        }

        var len_adjustment: usize = 0;
        if (range.step > 0 and range.stop == int.max(usize)) {
            range.stop = array.shape[i];
        } else if (range.step < 0 and range.start == int.max(usize)) {
            range.start = array.shape[i] - 1;
            len_adjustment = 1;
        }

        const len: usize = range.len() + len_adjustment;
        if (len == 1) {
            ndim -= 1;
            if (stride < 0) {
                offset -= range.start * scast(usize, int.abs(stride));
            } else {
                offset += range.start * scast(usize, stride);
            }
        } else {
            shapes[j] = len;
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
        .shape = shapes,
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
