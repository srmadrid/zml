//! Storage scheme:
//!
//! 1-d array of length `n` holding the diagonal elements.

const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const array = @import("../../array.zig");
const Dense = @import("../../array.zig").Dense;

const dense = @import("../dense.zig");
const strided = @import("../strided.zig");
const Strided = strided.Strided;

pub inline fn init(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
) !Dense(T) {
    const size: u32 = shape[0] * shape[1];
    return Dense(T){
        .data = try allocator.alloc(T, int.min(shape[0], shape[1])),
        .ndim = 2,
        .shape = .{ shape[0], shape[1] } ++ .{0} ** (array.max_dimensions - 2),
        .size = size,
        .strides = .{0} ** array.max_dimensions,
        .base = null,
        .flags = .{
            .owns_data = true,
        },
        .kind = .{ .diagonal = .{} },
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

inline fn checkPosition(comptime T: type, arr: *const Dense(T), position: []const u32) !void {
    if (position.len != arr.ndim) // 2
        return array.Error.DimensionMismatch;

    if (position[0] >= arr.shape[0] or position[1] >= arr.shape[1])
        return array.Error.PositionOutOfBounds;
}

inline fn checkBounds(position: []const u32) !void {
    if (position[0] != position[1])
        return array.Error.PositionOutOfBounds;
}

pub inline fn set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) !void {
    try checkPosition(T, arr, position);
    try checkBounds(position);

    arr.data[position[0]] = value;
}

pub inline fn _set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) void {
    // Assumes position is valid, i.e., within matrix bounds, and on the
    // diagonal.
    arr.data[position[0]] = value;
}

pub inline fn get(comptime T: type, arr: *const Dense(T), position: []const u32) !T {
    try checkPosition(T, arr, position);

    if (position[0] == position[1]) {
        return arr.data[position[0]];
    } else {
        return constants.zero(T, .{}) catch unreachable;
    }
}

pub inline fn _get(comptime T: type, arr: *const Dense(T), position: []const u32) T {
    // Assumes position is valid, i.e., within matrix bounds, and on the
    // diagonal.
    return arr.data[position[0]];
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
        .shape = .{ arr.shape[1], arr.shape[0] } ++ .{0} ** (array.max_dimensions - 2),
        .strides = .{0} ** array.max_dimensions,
        .size = arr.size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Underlying order remains the same.
            .owns_data = false,
        },
        .kind = .diagonal,
    };
}
