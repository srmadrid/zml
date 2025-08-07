//! Storage scheme:
//!
//! Three 1-d arrays of lengths `n - 1`, `n`, and `n - 1` are used to store the
//! subdiagonal, diagonal, and superdiagonal elements of an `n x n` tridiagonal
//! matrix. They are stored in a single 1-d array of length `3 * n - 2` in that
//! order.

const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");

const array = @import("../../array.zig");
const Dense = @import("../../array.zig").Dense;

const dense = @import("../dense.zig");

pub inline fn init(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
) !Dense(T) {
    if (shape[0] != shape[1] or shape[0] == 0)
        return array.Error.InvalidShape;

    const size: u32 = shape[0] * shape[1];
    return Dense(T){
        .data = try allocator.alloc(T, 3 * shape[0] - 2),
        .ndim = 2,
        .shape = .{ shape[0], shape[1] } ++ .{0} ** (array.max_dimensions - 2),
        .size = size,
        .strides = .{0} ** array.max_dimensions,
        .base = null,
        .flags = .{
            .owns_data = true,
        },
        .kind = .{ .tridiagonal = .{
            .superdiagonal_offset = shape[0] + (shape[0] - 1),
        } },
    };
}

pub inline fn full(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
    value: anytype,
    ctx: anytype,
) !Dense(T) {
    _ = ctx;
    var arr: Dense(T) = try init(T, allocator, shape);
    errdefer arr.deinit(allocator);

    if (comptime !types.isArbitraryPrecision(T)) {
        const casted_value: T = types.scast(T, value);

        for (0..arr.size) |i| {
            arr.data[i] = casted_value;
        }
    } else {
        @compileError("Arbitrary precision types not implemented yet");
    }

    return arr;
}

inline fn index(comptime T: type, arr: *const Dense(T), position: []const u32) u32 {
    const diff: i32 = types.scast(i32, position[1]) - types.scast(i32, position[0]);
    return switch (diff) {
        -1 => position[0] + arr.shape[0] + (arr.shape[0] - 1) - arr.kind.tridiagonal.superdiagonal_offset - 1, // Subdiagonal
        0 => position[0] + (arr.shape[0] - 1), // Diagonal
        1 => position[1] + arr.kind.tridiagonal.superdiagonal_offset - 1, // Superdiagonal
        else => 0, // Trash value, should not happen
    };
}

inline fn checkPosition(comptime T: type, arr: *const Dense(T), position: []const u32) !void {
    if (position.len != arr.ndim) // 2
        return array.Error.DimensionMismatch;

    if (position[0] >= arr.shape[0] or position[1] >= arr.shape[0])
        return array.Error.PositionOutOfBounds;
}

inline fn checkBounds(position: []const u32) !void {
    const diff: i32 = types.scast(i32, position[1]) - types.scast(i32, position[0]);
    if (diff < -1 or diff > 1)
        return array.Error.PositionOutOfBounds;
}

pub inline fn set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) !void {
    try checkPosition(T, arr, position);
    try checkBounds(position);

    arr.data[index(T, arr, position)] = value;
}

pub inline fn _set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) void {
    // Assumes position is valid, i.e., within matrix bounds, and in banded
    // range.
    arr.data[index(T, arr, position)] = value;
}

pub inline fn get(comptime T: type, arr: *const Dense(T), position: []const u32) !T {
    try checkPosition(T, arr, position);

    const diff: i32 = types.scast(i32, position[1]) - types.scast(i32, position[0]);
    if (diff < -1 or diff > 1) {
        return constants.zero(T, .{}) catch unreachable;
    } else {
        return arr.data[index(T, arr, position)];
    }
}

pub inline fn _get(comptime T: type, arr: *const Dense(T), position: []const u32) T {
    // Assumes position is valid, i.e., within matrix bounds, and in banded
    // range.
    return arr.data[index(T, arr, position)];
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
        .strides = arr.strides,
        .size = arr.size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Underlying order remains the same.
            .owns_data = false,
        },
        .kind = .{
            .tridiagonal = .{
                .superdiagonal_offset = arr.shape[0] + (arr.shape[0] - 1) - arr.kind.tridiagonal.superdiagonal_offset,
            },
        },
    };
}
