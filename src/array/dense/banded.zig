//! Storage scheme:
//!
//! A banded matrix with `lower` subdiagonals and `upper` superdiagonals is
//! stored in a rectangular array of size `n × (lower + upper + 1)`, where `n`
//! is the matrix dimension. The storage layout depends on the matrix order
//! (row-major vs column-major):
//! - For row-major: Each row of storage corresponds to a matrix row, each
//! column to a diagonal band.
//! - For col-major: Each column of storage corresponds to a matrix column, each
//! row to a diagonal band.
//!
//! Example 1: 5×5 matrix with lower=2, upper=1 (bandwidth = 4), order=row-major
//!
//! ```zig
//! Original matrix:        Banded storage:
//! [d₀ u₀  0  0  0]       [*  *  d₀ u₀]  ← row 0
//! [l₀ d₁ u₁  0  0]       [*  l₀ d₁ u₁]  ← row 1
//! [l₁ l₂ d₂ u₂  0]   →   [l₁ l₂ d₂ u₂]  ← row 2
//! [ 0 l₃ l₄ d₃ u₃]       [l₃ l₄ d₃ u₃]  ← row 3
//! [ 0  0 l₅ l₆ d₄]       [l₅ l₆ d₄  *]  ← row 4
//! ```
//!
//! Storage mapping: matrix element (i,j) → storage position (i, lower + j - i)
//!
//! Example 2: 5×5 matrix with lower=2, upper=1 (bandwidth = 4), order=col-major
//!
//! ```zig
//! Original matrix:        Banded storage:
//! [d₀ u₀  0  0  0]       [*    l₀   l₁    0    0]
//! [l₀ d₁ u₁  0  0]       [*    d₁   l₂   l₃    0]
//! [l₁ l₂ d₂ u₂  0]   →   [d₀   u₁   d₂   l₄   l₅]
//! [ 0 l₃ l₄ d₃ u₃]       [u₀   u₂   u₃   d₃   l₆]
//! [ 0  0 l₅ l₆ d₄]       [ 0    0    0   u₃   d₄]
//!                          ↑    ↑    ↑    ↑    ↑
//!                        col0 col1 col2 col3 col4
//! ```
//!
//! Storage mapping: matrix element (i,j) → storage position (upper + i - j, j)
//!
//! Where `*` represents unused storage locations. The storage uses
//! `n × (lower + upper + 1)` = `5 × 4` = 20 elements total.
//!
//! row-major column mapping:
//! - Column 0: 2nd subdiagonal   (j - i = -2)
//! - Column 1: 1st subdiagonal   (j - i = -1)
//! - Column 2: Main diagonal     (j - i =  0)
//! - Column 3: 1st superdiagonal (j - i = +1)
//!
//! col-major row mapping:
//! - Row 0: 1st superdiagonal (i - j = -1)
//! - Row 1: Main diagonal     (i - j =  0)
//! - Row 2: 1st subdiagonal   (i - j = +1)
//! - Row 3: 2nd subdiagonal   (i - j = +2)

const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const array = @import("../../array.zig");
const Dense = @import("../../array.zig").Dense;

const dense = @import("../dense.zig");

pub inline fn init(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
    order: array.Order,
    lower: u16,
    upper: u16,
) !Dense(T) {
    if (lower >= shape[0] or upper >= shape[1])
        return array.Error.InvalidFlags;

    const size: u32 = shape[0] * shape[1];
    return Dense(T){
        .data = try allocator.alloc(T, (lower + upper + 1) * (if (order == .col_major) shape[1] else shape[0])),
        .ndim = 2,
        .shape = .{ shape[0], shape[1] } ++ .{0} ** (array.max_dimensions - 2),
        .size = size,
        .strides = if (order == .col_major)
            .{ 1, lower + upper + 1 } ++ .{0} ** (array.max_dimensions - 2)
        else
            .{ upper + lower + 1, 1 } ++ .{0} ** (array.max_dimensions - 2),
        .base = null,
        .flags = .{
            .order = order,
            .owns_data = true,
        },
        .kind = .{ .banded = .{
            .lower = lower,
            .upper = upper,
        } },
    };
}

pub inline fn full(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
    value: anytype,
    order: array.Order,
    lower: u32,
    upper: bool,
    ctx: anytype,
) !Dense(T) {
    _ = ctx;
    var arr: Dense(T) = try init(T, allocator, shape, order, upper);
    errdefer arr.deinit(allocator);

    if (comptime !types.isArbitraryPrecision(T)) {
        const casted_value: T = types.scast(T, value);
        const m: u32 = shape[0];
        const n: u32 = shape[1];

        if (order == .col_major) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                var i: u32 = if (j < upper) 0 else j - upper;
                while (i <= int.min(m - 1, j + lower)) : (i += 1) {
                    arr.data[(upper + i - j) + j * (lower + upper + 1)] = casted_value;
                }
            }
        } else {
            var i: u32 = 0;
            while (i < m) : (i += 1) {
                var j: u32 = if (i < lower) 0 else i - lower;
                while (j <= int.min(n - 1, i + upper)) : (j += 1) {
                    arr.data[i * (lower + upper + 1) + (lower + j - i)] = casted_value;
                }
            }
        }
    } else {
        @compileError("Arbitrary precision types not implemented yet");
    }

    return arr;
}

inline fn index(comptime T: type, arr: *const Dense(T), position: []const u32) u32 {
    if (arr.flags.order == .col_major) {
        return (arr.kind.banded.upper + position[0] - position[1]) * arr.strides[0] + position[1] * arr.strides[1];
    } else {
        return position[0] * arr.strides[0] + (arr.kind.banded.lower + position[1] - position[0]) * arr.strides[1];
    }
}

inline fn checkPosition(comptime T: type, arr: *const Dense(T), position: []const u32) !void {
    if (position.len != arr.ndim) // 2
        return array.Error.DimensionMismatch;

    if (position[0] >= arr.shape[0] or position[1] >= arr.shape[1])
        return array.Error.PositionOutOfBounds;
}

inline fn checkBounds(comptime T: type, arr: *const Dense(T), position: []const u32) !void {
    if (position[1] + arr.kind.banded.lower < position[0] or
        position[1] > position[0] + arr.kind.banded.upper)
        return array.Error.PositionOutOfBounds;
}

pub inline fn set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) !void {
    try checkPosition(T, arr, position);
    try checkBounds(T, arr, position);

    arr.data[index(T, arr, position)] = value;
}

pub inline fn _set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) void {
    // Assumes position is valid, i.e., within matrix bounds, and in banded
    // range.
    arr.data[index(T, arr, position)] = value;
}

pub inline fn get(comptime T: type, arr: *const Dense(T), position: []const u32) !T {
    try checkPosition(T, arr, position);

    if (position[1] + arr.kind.banded.lower < position[0] or
        position[1] > position[0] + arr.kind.banded.upper)
    {
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
        .shape = .{ arr.shape[1], arr.shape[0] } ++ .{0} ** (array.max_dimensions - 2),
        .strides = .{ arr.strides[1], arr.strides[0] } ++ .{0} ** (array.max_dimensions - 2),
        .size = arr.size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = if (arr.flags.order == .col_major) .row_major else .col_major, // Must switch order
            .owns_data = false,
        },
        .kind = .{
            .banded = .{
                .lower = arr.kind.banded.upper,
                .upper = arr.kind.banded.lower,
            },
        },
    };
}
