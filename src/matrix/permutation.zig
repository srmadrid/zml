//! Storage scheme:
//!
//! A 1-d array of size `n` stores the permutation of 0..n-1. If `direction` is
//! `.forward`, the element at index `i` indicates the column index of the 1 in
//! row `i`, i.e., if `data[i] = j`, then the element at row `i` and column `j`
//! is 1, and all other elements in row `i` are 0. If `direction` is
//! `.backward`, the same applies but for columns, i.e., if `data[j] = i`,
//! then the element at row `i` and column `j` is 1, and all other elements in
//! column `j` are 0.

const std = @import("std");

const types = @import("../types.zig");
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const ReturnType2 = types.ReturnType2;
const Order = types.Order;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../array.zig");

pub const Direction = enum {
    forward,
    backward,
};

pub fn Permutation(T: type) type {
    if (!types.isNumeric(T))
        @compileError("Permutation requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]u32,
        size: u32,
        direction: Direction = .forward,
        flags: Flags = .{},

        pub const empty: Permutation(T) = .{
            .data = &.{},
            .size = 0,
            .direction = .forward,
            .flags = .{ .owns_data = false },
        };

        pub fn tp() type {
            return T;
        }

        pub fn init(
            allocator: std.mem.Allocator,
            size: u32,
        ) !Permutation(T) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(u32, size)).ptr,
                .size = size,
                .direction = .forward,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn eye(
            allocator: std.mem.Allocator,
            size: u32,
        ) !Permutation(T) {
            const mat: Permutation(T) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            while (i < size) : (i += 1) {
                mat.data[i] = i;
            }

            return mat;
        }

        pub fn deinit(self: *Permutation(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0..self.size]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Permutation(T), row: u32, col: u32) !T {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            if (self.direction == .forward) {
                if (self.data[row] == col) {
                    return constants.one(T, .{}) catch unreachable;
                } else {
                    return constants.zero(T, .{}) catch unreachable;
                }
            } else {
                if (self.data[col] == row) {
                    return constants.one(T, .{}) catch unreachable;
                } else {
                    return constants.zero(T, .{}) catch unreachable;
                }
            }
        }

        pub inline fn at(self: *const Permutation(T), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and
            // in banded range.
            if (self.direction == .forward) {
                if (self.data[row] == col) {
                    return constants.one(T, .{}) catch unreachable;
                } else {
                    return constants.zero(T, .{}) catch unreachable;
                }
            } else {
                if (self.data[col] == row) {
                    return constants.one(T, .{}) catch unreachable;
                } else {
                    return constants.zero(T, .{}) catch unreachable;
                }
            }
        }

        // pub fn set(self: *Permutation(T), row: u32, col: u32, value: u32) !void {
        //     if (row >= self.size or col >= self.size)
        //         return matrix.Error.PositionOutOfBounds;

        //     if (value != 0 and value != 1)
        //         return matrix.Error.BreaksStructure;
        // }

        // pub inline fn put(self: *Permutation(T), row: u32, col: u32, value: u32) void {
        //     // Unchecked version of set. Assumes row and col are valid and
        //     // in banded range.
        //     if (value == 1) {
        //         self.data[row] = col;
        //     }
        // }

        pub fn copyToGeneralDenseMatrix(
            self: Permutation(T),
            allocator: std.mem.Allocator,
            comptime order: Order,
            ctx: anytype,
        ) !matrix.general.Dense(T, order) {
            var result: matrix.general.Dense(T, order) = try .init(allocator, self.size, self.size);
            errdefer result.deinit(allocator);

            if (comptime order == .col_major) {
                if (self.direction == .forward) {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i + j * result.ld] = if (self.data[i] == j)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[i + j * result.ld] = if (self.data[j] == i)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                if (self.direction == .forward) {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[i * result.ld + j] = if (self.data[i] == j)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i * result.ld + j] = if (self.data[j] == i)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            }

            return result;
        }

        pub fn copyToDenseArray(
            self: *const Permutation(T),
            allocator: std.mem.Allocator,
            comptime order: Order,
            ctx: anytype,
        ) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.size, self.size });
            errdefer result.deinit(allocator);

            if (comptime order == .col_major) {
                if (self.direction == .forward) {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i + j * result.strides[0]] = if (self.data[i] == j)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[i + j * result.strides[0]] = if (self.data[j] == i)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                if (self.direction == .forward) {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[i * result.strides[0] + j] = if (self.data[i] == j)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i * result.strides[0] + j] = if (self.data[j] == i)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            }

            return result;
        }

        pub fn transpose(self: Permutation(T)) Permutation(T) {
            return .{
                .data = self.data,
                .size = self.size,
                .direction = if (self.direction == .forward) .backward else .forward,
                .flags = self.flags,
            };
        }

        // pub fn submatrix(
        //     self: *const Permutation(T),
        //     start: u32,
        //     end: u32,
        // ) !? {
        //     if (start >= self.size or end > self.size or start >= end)
        //         return matrix.Error.InvalidRange;

        //     const sub_size = end - start;

        //     return .{
        //         .data = self.data,
        //         .size = sub_size,
        //         .osize = self.osize,
        //         .offset = self.offset + start,
        //         .sdoffset = self.sdoffset,
        //         .flags = .{
        //             .owns_data = false,
        //         },
        //     };
        // }
    };
}
