//! Storage scheme:
//!
//! 1-d array of length `n` holding the diagonal elements.

const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const General = matrix.General;
const Flags = matrix.Flags;

const array = @import("../array.zig");
const Dense = array.Dense;

pub fn Diagonal(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("Diagonal requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        flags: Flags = .{},

        pub const empty: Diagonal(T) = .{
            .data = &.{},
            .size = 0,
            .flags = .{ .order = .col_major, .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            size: u32,
        ) !Diagonal(T) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return Diagonal(T){
                .data = (try allocator.alloc(T, size)).ptr,
                .size = size,
                .flags = .{
                    .order = .col_major,
                    .owns_data = true,
                },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            size: u32,
            value: anytype,
            ctx: anytype,
        ) !Diagonal(T) {
            const mat: Diagonal(T) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                var i: u32 = 0;
                while (i < size) : (i += 1) {
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
            const mat: Diagonal(T) = try .init(allocator, size);
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
                allocator.?.free(self.data[0..self.size]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Diagonal(T), row: u32, col: u32) !T {
            if (row >= self.size or col >= self.size)
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
            if (row >= self.size or col >= self.size)
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

        pub fn toGeneral(self: Diagonal(T), allocator: std.mem.Allocator, ctx: anytype) !General(T) {
            var result: General(T) = try .init(allocator, self.size, self.size, .{ .order = self.flags.order });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                var j: u32 = 0;
                while (j < self.size) : (j += 1) {
                    var i: u32 = 0;
                    while (i < j) : (i += 1) {
                        result.data[i + j * self.size] = constants.zero(T, ctx) catch unreachable;
                    }

                    result.data[j * self.size + j] = self.data[j];

                    i = j + 1;
                    while (i < self.size) : (i += 1) {
                        result.data[i + j * self.size] = constants.zero(T, ctx) catch unreachable;
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn toDense(self: *const Diagonal(T), allocator: std.mem.Allocator, ctx: anytype) !Dense(T) {
            var result: Dense(T) = try .init(allocator, &.{ self.size, self.size }, .{ .order = self.flags.order });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                var j: u32 = 0;
                while (j < self.size) : (j += 1) {
                    var i: u32 = 0;
                    while (i < j) : (i += 1) {
                        result.data[i + j * self.size] = constants.zero(T, ctx) catch unreachable;
                    }

                    result.data[j * self.size + j] = self.data[j];

                    i = j + 1;
                    while (i < self.size) : (i += 1) {
                        result.data[i + j * self.size] = constants.zero(T, ctx) catch unreachable;
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
                .size = self.size,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Diagonal(T),
            start: u32,
            end: u32,
        ) !Diagonal(T) {
            if (start >= self.size or end > self.size or start >= end)
                return matrix.Error.InvalidRange;

            const sub_size = end - start;

            return Diagonal(T){
                .data = self.data + start,
                .size = sub_size,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }
    };
}
