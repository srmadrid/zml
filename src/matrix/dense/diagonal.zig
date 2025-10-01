//! Storage scheme:
//!
//! 1-d array of length `n` holding the diagonal elements.

const std = @import("std");

const types = @import("../../types.zig");
const EnsureMatrix = types.EnsureMatrix;
const Numeric = types.Numeric;
const ReturnType2 = types.ReturnType2;
const Coerce = types.Coerce;
const Order = types.Order;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

pub fn Diagonal(T: type) type {
    if (!types.isNumeric(T))
        @compileError("Diagonal requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        flags: Flags = .{},

        pub const empty: Diagonal(T) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
        ) !Diagonal(T) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return Diagonal(T){
                .data = (try allocator.alloc(T, int.min(rows, cols))).ptr,
                .rows = rows,
                .cols = cols,
                .flags = .{
                    .owns_data = true,
                },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            value: anytype,
            ctx: anytype,
        ) !Diagonal(T) {
            const mat: Diagonal(T) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                var i: u32 = 0;
                while (i < int.min(rows, cols)) : (i += 1) {
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
        ) !Diagonal(T) {
            const mat: Diagonal(T) = try .init(allocator, size, size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                var i: u32 = 0;
                while (i < size) : (i += 1) {
                    mat.data[i] = constants.one(T, ctx) catch unreachable;
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn deinit(self: *Diagonal(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0..(int.min(self.rows, self.cols))]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Diagonal(T), row: u32, col: u32) !T {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (row != col)
                return constants.zero(T, .{}) catch unreachable;

            return self.data[row];
        }

        pub inline fn at(self: *const Diagonal(T), row: u32, col: u32) T {
            _ = col;
            // Unchecked version of get. Assumes row and col are equal and
            // within bounds.
            return self.data[row];
        }

        pub fn set(self: *Diagonal(T), row: u32, col: u32, value: T) !void {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (row != col)
                return matrix.Error.PositionOutOfBounds;

            self.data[row] = value;
        }

        pub inline fn put(self: *Diagonal(T), row: u32, col: u32, value: T) void {
            _ = col;
            // Unchecked version of set. Assumes row and col are equal and
            // within bounds.
            self.data[row] = value;
        }

        pub fn toGeneralDenseMatrix(self: Diagonal(T), allocator: std.mem.Allocator, comptime order: Order, ctx: anytype) !matrix.dense.General(T, order) {
            var result: matrix.dense.General(T, order) = try .init(allocator, self.rows, self.cols);
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, self.rows)) : (i += 1) {
                            result.data[i + j * result.ld] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (j < int.min(self.rows, self.cols)) {
                            result.data[j * result.ld + j] = self.data[j];
                        }

                        i = j + 1;
                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * result.ld] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, self.cols)) : (j += 1) {
                            result.data[i * result.ld + j] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (i < int.min(self.rows, self.cols)) {
                            result.data[i * result.ld + i] = self.data[i];
                        }

                        j = i + 1;
                        while (j < self.cols) : (j += 1) {
                            result.data[i * result.ld + j] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn toDenseArray(self: *const Diagonal(T), allocator: std.mem.Allocator, comptime order: Order, ctx: anytype) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.rows, self.cols });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, self.rows)) : (i += 1) {
                            result.data[i + j * result.ld] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (j < int.min(self.rows, self.cols)) {
                            result.data[j * result.ld + j] = self.data[j];
                        }

                        i = j + 1;
                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * result.ld] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, self.cols)) : (j += 1) {
                            result.data[i * result.ld + j] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (i < int.min(self.rows, self.cols)) {
                            result.data[i * result.ld + i] = self.data[i];
                        }

                        j = i + 1;
                        while (j < self.cols) : (j += 1) {
                            result.data[i * result.ld + j] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn transpose(self: Diagonal(T)) Diagonal(T) {
            return Diagonal(T){
                .data = self.data,
                .rows = self.cols,
                .cols = self.rows,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Diagonal(T),
            start: u32,
            row_end: u32,
            col_end: u32,
        ) !Diagonal(T) {
            if (start >= int.min(self.rows, self.cols) or
                row_end > self.rows or col_end > self.cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - start;
            const sub_cols = col_end - start;

            return Diagonal(T){
                .data = self.data + start,
                .rows = sub_rows,
                .cols = sub_cols,
                .flags = .{
                    .owns_data = false,
                },
            };
        }
    };
}
