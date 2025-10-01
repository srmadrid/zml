//! Storage scheme:
//!
//! Three 1-d arrays of lengths `n - 1`, `n`, and `n - 1` are used to store the
//! subdiagonal, diagonal, and superdiagonal elements of an `n x n` tridiagonal
//! matrix. They are stored in a single 1-d array of length `3 * n - 2` in that
//! order.

const std = @import("std");

const types = @import("../../types.zig");
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const ReturnType2 = types.ReturnType2;
const Order = types.Order;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

pub fn Tridiagonal(T: type) type {
    if (!types.isNumeric(T))
        @compileError("Tridiagonal requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        osize: u32,
        offset: u32,
        sdoffset: u32,
        flags: Flags = .{},

        pub const empty: Tridiagonal(T) = .{
            .data = &.{},
            .size = 0,
            .osize = 0,
            .offset = 0,
            .sdoffset = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            size: u32,
        ) !Tridiagonal(T) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(T, 3 * size - 2)).ptr,
                .size = size,
                .osize = size,
                .offset = 0,
                .sdoffset = (size - 1) + size,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            size: u32,
            value: anytype,
            ctx: anytype,
        ) !Tridiagonal(T) {
            const mat: Tridiagonal(T) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                var i: u32 = 0;
                while (i < 3 * size - 2) : (i += 1) {
                    mat.data[i] = value_casted;
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn eye(
            allocator: std.mem.Allocator,
            size: u32,
            ctx: anytype,
        ) !Tridiagonal(T) {
            const mat: Tridiagonal(T) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                var i: u32 = 0;
                while (i < size - 1) : (i += 1) {
                    mat.data[i] = constants.zero(T, .{}) catch unreachable;
                    mat.data[mat.sdoffset + i] = constants.zero(T, .{}) catch unreachable;
                }

                i = 0;
                while (i < size) : (i += 1) {
                    mat.data[(mat.size - 1) + i] = constants.one(T, .{}) catch unreachable;
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn deinit(self: *Tridiagonal(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0..(3 * self.size - 2)]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Tridiagonal(T), row: u32, col: u32) !T {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            const diff: i32 = types.scast(i32, col) - types.scast(i32, row);
            if (diff < -1 or diff > 1)
                return constants.zero(T, .{}) catch unreachable;

            const idx: u32 = self.offset + switch (diff) {
                -1 => row + self.osize + (self.osize - 1) - self.sdoffset - 1, // Subdiagonal
                0 => row + (self.osize - 1), // Diagonal
                1 => col + self.sdoffset - 1, // Superdiagonal
                else => 0, // Trash value, should not happen
            };

            return self.data[idx];
        }

        pub inline fn at(self: *const Tridiagonal(T), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and
            // in banded range.
            const diff: i32 = types.scast(i32, col) - types.scast(i32, row);
            const idx: u32 = self.offset + switch (diff) {
                -1 => row + self.osize + (self.osize - 1) - self.sdoffset - 1, // Subdiagonal
                0 => row + (self.osize - 1), // Diagonal
                1 => col + self.sdoffset - 1, // Superdiagonal
                else => 0, // Trash value, should not happen
            };

            return self.data[idx];
        }

        pub fn set(self: *Tridiagonal(T), row: u32, col: u32, value: T) !void {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            const diff: i32 = types.scast(i32, col) - types.scast(i32, row);
            if (diff < -1 or diff > 1)
                return matrix.Error.PositionOutOfBounds;

            const idx: u32 = self.offset + switch (diff) {
                -1 => row + self.osize + (self.osize - 1) - self.sdoffset - 1, // Subdiagonal
                0 => row + (self.osize - 1), // Diagonal
                1 => col + self.sdoffset - 1, // Superdiagonal
                else => 0, // Trash value, should not happen
            };

            self.data[idx] = value;
        }

        pub inline fn put(self: *Tridiagonal(T), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid and
            // in banded range.
            const diff: i32 = types.scast(i32, col) - types.scast(i32, row);
            const idx: u32 = self.offset + switch (diff) {
                -1 => row + self.osize + (self.osize - 1) - self.sdoffset - 1, // Subdiagonal
                0 => row + (self.osize - 1), // Diagonal
                1 => col + self.sdoffset - 1, // Superdiagonal
                else => 0, // Trash value, should not happen
            };

            self.data[idx] = value;
        }

        pub fn toGeneralDenseMatrix(self: Tridiagonal(T), allocator: std.mem.Allocator, comptime order: Order, ctx: anytype) !matrix.dense.General(T, order) {
            var result: matrix.dense.General(T, order) = try .init(allocator, self.size, self.size);
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            if (i == j) { // Diagonal
                                result.data[j + j * result.ld] = self.data[self.offset + j + (self.osize - 1)];
                            } else if (i == j + 1) { // Subdiagonal
                                result.data[i + j * result.ld] = self.data[self.offset + i + self.osize + (self.osize - 1) - self.sdoffset - 1];
                            } else if (i + 1 == j) { // Superdiagonal
                                result.data[i + j * result.ld] = self.data[self.offset + j + self.sdoffset - 1];
                            } else {
                                result.data[i + j * result.ld] = constants.zero(T, .{}) catch unreachable;
                            }
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            if (i == j) { // Diagonal
                                result.data[i * result.ld + i] = self.data[self.offset + i + (self.osize - 1)];
                            } else if (i == j + 1) { // Subdiagonal
                                result.data[i * result.ld + j] = self.data[self.offset + i + self.osize + (self.osize - 1) - self.sdoffset - 1];
                            } else if (i + 1 == j) { // Superdiagonal
                                result.data[i * result.ld + j] = self.data[self.offset + j + self.sdoffset - 1];
                            } else {
                                result.data[i * result.ld + j] = constants.zero(T, .{}) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn toDenseArray(self: *const Tridiagonal(T), allocator: std.mem.Allocator, comptime order: Order, ctx: anytype) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.size, self.size });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            if (i == j) { // Diagonal
                                result.data[j + j * result.strides[1]] = self.data[self.offset + j + (self.osize - 1)];
                            } else if (i == j + 1) { // Subdiagonal
                                result.data[i + j * result.strides[1]] = self.data[self.offset + i + self.osize + (self.osize - 1) - self.sdoffset - 1];
                            } else if (i + 1 == j) { // Superdiagonal
                                result.data[i + j * result.strides[1]] = self.data[self.offset + j + self.sdoffset - 1];
                            } else {
                                result.data[i + j * result.strides[1]] = constants.zero(T, .{}) catch unreachable;
                            }
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            if (i == j) { // Diagonal
                                result.data[i * result.strides[0] + i] = self.data[self.offset + i + (self.osize - 1)];
                            } else if (i == j + 1) { // Subdiagonal
                                result.data[i * result.strides[0] + j] = self.data[self.offset + i + self.osize + (self.osize - 1) - self.sdoffset - 1];
                            } else if (i + 1 == j) { // Superdiagonal
                                result.data[i * result.strides[0] + j] = self.data[self.offset + j + self.sdoffset - 1];
                            } else {
                                result.data[i * result.strides[0] + j] = constants.zero(T, .{}) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn transpose(self: Tridiagonal(T)) Tridiagonal(T) {
            return .{
                .data = self.data,
                .size = self.size,
                .osize = self.osize,
                .offset = self.offset,
                .sdoffset = self.size + (self.size - 1) - self.sdoffset,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Tridiagonal(T),
            start: u32,
            end: u32,
        ) !Tridiagonal(T) {
            if (start >= self.size or end > self.size or start >= end)
                return matrix.Error.InvalidRange;

            const sub_size = end - start;

            return .{
                .data = self.data,
                .size = sub_size,
                .osize = self.osize,
                .offset = self.offset + start,
                .sdoffset = self.sdoffset,
                .flags = .{
                    .owns_data = false,
                },
            };
        }
    };
}
