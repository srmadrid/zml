const std = @import("std");

const types = @import("../types.zig");
const Scalar = types.Scalar;
const ReturnType1 = types.ReturnType1;
const scast = types.scast;
const cast = types.cast;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");
const Array = array.Array;
const Range = array.Range;

const dense = @import("dense.zig");

pub inline fn index(ndim: usize, strides: [array.maxDimensions]isize, offset: usize, position: []const usize) usize {
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

pub inline fn checkPosition(ndim: usize, shape: [array.maxDimensions]usize, position: []const usize) !void {
    if (position.len != ndim) {
        return array.Error.DimensionMismatch;
    }

    for (0..position.len) |i| {
        if (position[i] >= shape[i]) {
            return array.Error.PositionOutOfBounds;
        }
    }
}

pub inline fn set(comptime T: type, arr: *Array(T), position: []const usize, value: T) !void {
    try checkPosition(arr.ndim, arr.shape, position);

    arr.data[index(arr.ndim, arr.metadata.strided.strides, arr.metadata.strided.offset, position)] = value;
}

pub fn get(comptime T: type, arr: Array(T), position: []const usize) !*T {
    try checkPosition(arr.ndim, arr.shape, position);

    return &arr.data[index(arr.ndim, arr.metadata.strided.strides, arr.metadata.strided.offset, position)];
}

pub inline fn slice(comptime T: type, arr: *const Array(T), ranges: []const Range) !Array(T) {
    var ndim: usize = arr.ndim;
    var size: usize = 1;
    var shape: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    var strides: [array.maxDimensions]isize = .{0} ** array.maxDimensions;
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

pub inline fn reshape(comptime T: type, arr: *const Array(T), shape: []const usize) !Array(T) {
    var new_size: usize = 1;
    var new_shape: [array.maxDimensions]usize = .{0} ** array.maxDimensions;
    var new_strides: [array.maxDimensions]isize = .{0} ** array.maxDimensions;
    if (shape.len > 0) {
        for (0..shape.len) |i| {
            const idx: usize = if (arr.order == .rowMajor) shape.len - i - 1 else i;

            new_strides[idx] = new_size;
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
            .storage = .strided,
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
    var newarr: Array(ReturnType1(op, T)) = try dense.init(allocator, ReturnType1(op, T), arr.shape[0..arr.ndim], arr.flags.order);
    errdefer newarr.deinit(allocator);
    newarr.flags.writeable = writeable;

    var j: usize = 0;
    errdefer dense.cleanup(ReturnType1(op, T), allocator, newarr.data[0..j]);
    var iter = array.Iterator(T).init(&arr);

    const opinfo = @typeInfo(@TypeOf(op));
    for (0..newarr.size) |i| {
        if (opinfo.@"fn".params.len == 1) {
            newarr.data[i] = op(arr.data[iter.index]);
        } else if (opinfo.@"fn".params.len == 2) {
            newarr.data[i] = try op(arr.data[iter.index], .{ .allocator = allocator });
        }

        j += 1;
        _ = iter.next();
    }

    return newarr;
}

pub inline fn apply1_(
    comptime T: type,
    arr: *Array(T),
    comptime op_: anytype,
    allocator: ?std.mem.Allocator,
) !void {
    var iter = array.Iterator(T).init(arr);

    const opinfo = @typeInfo(@TypeOf(op_));
    for (0..arr.size) |_| {
        if (opinfo.@"fn".params.len == 1) {
            op_(&arr.data[iter.index]);
        } else if (opinfo.@"fn".params.len == 2) {
            try op_(&arr.data[iter.index], .{ .allocator = allocator });
        }

        _ = iter.next();
    }

    return;
}

// For apply2 use recursive loop function
