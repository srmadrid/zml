//! Storage scheme:
//!
//! Full `n`-by-`n` storage only accessing the upper or lower triangular part of
//! the matrix.

const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");

const array = @import("../../array.zig");
const Dense = @import("../../array.zig").Dense;

const dense = @import("../dense.zig");

pub inline fn init(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
    order: array.Order,
    upper: bool,
) !Dense(T) {
    if (shape[0] != shape[1] or shape[0] == 0)
        return array.Error.InvalidShape;

    const size: u32 = shape[0] * shape[1];
    return Dense(T){
        .data = try allocator.alloc(T, size),
        .ndim = 2,
        .shape = .{ shape[0], shape[1] } ++ .{0} ** (array.max_dimensions - 2),
        .size = size,
        .strides = if (order == .col_major)
            .{ 1, shape[1] } ++ .{0} ** (array.max_dimensions - 2)
        else
            .{ shape[0], 1 } ++ .{0} ** (array.max_dimensions - 2),
        .base = null,
        .flags = .{
            .order = order,
            .owns_data = true,
        },
        .kind = .{ .symmetric = .{ .upper = upper } },
    };
}

pub inline fn full(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
    value: anytype,
    order: array.Order,
    upper: bool,
    ctx: anytype,
) !Dense(T) {
    _ = ctx;
    var arr: Dense(T) = try init(T, allocator, shape, order, upper);
    errdefer arr.deinit(allocator);

    if (comptime !types.isArbitraryPrecision(T)) {
        const casted_value: T = types.scast(T, value);
        const n: u32 = shape[0];

        if (order == .col_major) {
            if (upper) {
                var j: u32 = 0;
                while (j < n) : (j += 1) {
                    var i: u32 = 0;
                    while (i <= j) : (i += 1) {
                        arr.data[i + j * n] = casted_value;
                    }
                }
            } else {
                var j: u32 = 0;
                while (j < n) : (j += 1) {
                    var i: u32 = j;
                    while (i < n) : (i += 1) {
                        arr.data[i + j * n] = casted_value;
                    }
                }
            }
        } else {
            if (upper) {
                var i: u32 = 0;
                while (i < n) : (i += 1) {
                    var j: u32 = i;
                    while (j < n) : (j += 1) {
                        arr.data[i * n + j] = casted_value;
                    }
                }
            } else {
                var i: u32 = 0;
                while (i < n) : (i += 1) {
                    var j: u32 = 0;
                    while (j <= i) : (j += 1) {
                        arr.data[i * n + j] = casted_value;
                    }
                }
            }
        }
    } else {
        @compileError("Arbitrary precision types not implemented yet");
    }

    return arr;
}

inline fn index(comptime T: type, arr: *const Dense(T), position: []const u32) u32 {
    var i: u32 = position[0];
    var j: u32 = position[1];
    // Clamp indices to the valid range for symmetric matrices
    if (arr.kind.symmetric.upper) {
        if (i > j) {
            const temp: u32 = i;
            i = j;
            j = temp;
        }
    } else {
        if (i < j) {
            const temp: u32 = i;
            i = j;
            j = temp;
        }
    }

    return i * arr.strides[0] + j * arr.strides[1];
}

inline fn _index(comptime T: type, arr: *const Dense(T), position: []const u32) u32 {
    // Does not check if swapping is needed.
    return position[0] * arr.strides[0] + position[1] * arr.strides[1];
}

inline fn checkPosition(comptime T: type, arr: *const Dense(T), position: []const u32) !void {
    if (position.len != arr.ndim) // 2
        return array.Error.DimensionMismatch;

    if (position[0] >= arr.shape[0] or position[1] >= arr.shape[1])
        return array.Error.PositionOutOfBounds;
}

pub inline fn set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) !void {
    try checkPosition(T, arr, position);

    arr.data[index(T, arr, position)] = value;
}

pub inline fn _set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) void {
    // Assumes position is valid, i.e., within matrix bounds, and in the correct
    // triangular part.
    arr.data[_index(T, arr, position)] = value;
}

pub inline fn get(comptime T: type, arr: *const Dense(T), position: []const u32) !T {
    try checkPosition(T, arr, position);

    return arr.data[index(T, arr, position)];
}

pub inline fn _get(comptime T: type, arr: *const Dense(T), position: []const u32) T {
    // Assumes position is valid, i.e., within matrix bounds, and in the correct
    // triangular part.
    return arr.data[_index(T, arr, position)];
}

pub inline fn transpose(
    comptime T: type,
    arr: *Dense(T),
    axes: []const u32,
) !Dense(T) {
    if (axes[0] != 1 and axes[1] != 0) {
        return array.Error.InvalidAxes; // axes must be the standard permutation [1, 0]
    }

    return Dense(T){
        .data = arr.data,
        .ndim = arr.ndim,
        .shape = arr.shape,
        .strides = .{ arr.strides[1], arr.strides[0] } ++ .{0} ** (array.max_dimensions - 2),
        .size = arr.size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Underlying order remains the same.
            .owns_data = false,
        },
        .kind = .{
            .symmetric = .{
                .upper = !arr.kind.symmetric.upper, // Transpose flips the upper/lower part.
            },
        },
    };
}
