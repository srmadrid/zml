//! Storage scheme:
//!
//! A 1-d array of size `n` stores the permutation of 0..n-1. The element at
//! index `i` indicates the column index of the 1 in row `i`, i.e., if
//! `data[i] = j`, then the element at row `i` and column `j` is 1, and all
//! other elements in row `i` are 0.

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
const General = matrix.General;
const Flags = matrix.Flags;

const array = @import("../array.zig");
const Dense = array.Dense;

pub fn Permutation(T: type) type {
    if (!types.isNumeric(T))
        @compileError("Permutation requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]u32,
        size: u32,
        flags: Flags = .{},

        pub const empty: Permutation(T) = .{
            .data = &.{},
            .size = 0,
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

            if (self.data[row] == col) {
                return constants.one(T, .{}) catch unreachable;
            } else {
                return constants.zero(T, .{}) catch unreachable;
            }
        }

        pub inline fn at(self: *const Permutation(T), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and
            // in banded range.
            if (self.data[row] == col) {
                return constants.one(T, .{}) catch unreachable;
            } else {
                return constants.zero(T, .{}) catch unreachable;
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

        pub fn toGeneral(self: Permutation(T), allocator: std.mem.Allocator, comptime order: Order, ctx: anytype) !General(T, order) {
            var result: General(T, order) = try .init(allocator, self.size, self.size);
            errdefer result.deinit(allocator);

            if (comptime order == .col_major) {
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
                        result.data[i * result.ld + j] = if (self.data[i] == j)
                            constants.one(T, ctx) catch unreachable
                        else
                            constants.zero(T, ctx) catch unreachable;
                    }
                }
            }

            return result;
        }

        pub fn toDenseArray(self: *const Permutation(T), allocator: std.mem.Allocator, comptime order: Order, ctx: anytype) !Dense(T, order) {
            var result: Dense(T, order) = try .init(allocator, &.{ self.size, self.size });
            errdefer result.deinit(allocator);

            if (comptime order == .col_major) {
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
                        result.data[i * result.strides[0] + j] = if (self.data[i] == j)
                            constants.one(T, ctx) catch unreachable
                        else
                            constants.zero(T, ctx) catch unreachable;
                    }
                }
            }

            return result;
        }

        pub fn transpose(self: Permutation(T), allocator: std.mem.Allocator) !Permutation(T) {
            var result: Permutation(T) = try .init(allocator, self.size);
            errdefer result.deinit(allocator);

            var i: u32 = 0;
            while (i < self.size) : (i += 1) {
                result.data[self.data[i]] = i;
            }

            return result;
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

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));
    const R: type = EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y));

    if (comptime !types.isPermutationMatrix(@TypeOf(x))) {
        var result: R = try .init(allocator, y.size, y.size);
        errdefer result.deinit(allocator);

        var j: u32 = 0;
        while (j < y.size) : (j += 1) {
            var i: u32 = 0;
            while (i < y.size) : (i += 1) {
                if (comptime @typeInfo(@TypeOf(op)).@"fn".params.len == 2) {
                    result.data[i + j * result.ld] = op(x, y.at(i, j));
                } else if (comptime @typeInfo(@TypeOf(op)).@"fn".params.len == 3) {
                    result.data[i + j * result.ld] = try op(x, y.at(i, j), ctx);
                }
            }
        }

        return result;
    } else if (comptime !types.isPermutationMatrix(@TypeOf(y))) {
        var result: R = try .init(allocator, x.size, x.size);
        errdefer result.deinit(allocator);

        var j: u32 = 0;
        while (j < x.size) : (j += 1) {
            var i: u32 = 0;
            while (i < x.size) : (i += 1) {
                if (comptime @typeInfo(@TypeOf(op)).@"fn".params.len == 2) {
                    result.data[i + j * result.ld] = op(x.at(i, j), y);
                } else if (comptime @typeInfo(@TypeOf(op)).@"fn".params.len == 3) {
                    result.data[i + j * result.ld] = try op(x.at(i, j), y, ctx);
                }
            }
        }

        return result;
    }

    if (x.size != y.size)
        return matrix.Error.DimensionMismatch;

    var result: R = try .init(allocator, x.size, x.size);
    errdefer result.deinit(allocator);

    var j: u32 = 0;
    while (j < x.size) : (j += 1) {
        var i: u32 = 0;
        while (i < x.size) : (i += 1) {
            if (comptime @typeInfo(@TypeOf(op)).@"fn".params.len == 2) {
                result.data[i + j * result.ld] = op(x.at(i, j), y.at(i, j));
            } else if (comptime @typeInfo(@TypeOf(op)).@"fn".params.len == 3) {
                result.data[i + j * result.ld] = try op(x.at(i, j), y.at(i, j), ctx);
            }
        }
    }

    return result;
}
