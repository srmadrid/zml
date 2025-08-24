//! Storage scheme:
//!
//! 1-d array of length `n` holding the diagonal elements.

const std = @import("std");

const types = @import("../types.zig");
const EnsureMatrix = types.EnsureMatrix;
const Numeric = types.Numeric;
const ReturnType2 = types.ReturnType2;
const Coerce = types.Coerce;
const Order = types.Order;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const General = matrix.General;
const Flags = matrix.Flags;

const array = @import("../array.zig");
const Dense = array.Dense;

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

        pub fn toGeneral(self: Diagonal(T), allocator: std.mem.Allocator, comptime order: Order, ctx: anytype) !General(T, order) {
            var result: General(T, order) = try .init(allocator, self.rows, self.cols);
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < j and i < self.rows) : (i += 1) {
                            result.data[i + j * self.size] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (j < int.min(self.rows, self.cols)) {
                            result.data[j * self.size + j] = self.data[j];
                        }

                        i = j + 1;
                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * self.size] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < i and j < self.cols) : (j += 1) {
                            result.data[i * self.size + j] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (i < int.min(self.rows, self.cols)) {
                            result.data[i * self.size + i] = self.data[i];
                        }

                        j = i + 1;
                        while (j < self.cols) : (j += 1) {
                            result.data[i * self.size + j] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn toDenseArray(self: *const Diagonal(T), allocator: std.mem.Allocator, comptime order: Order, ctx: anytype) !Dense(T, order) {
            var result: Dense(T, order) = try .init(allocator, &.{ self.rows, self.cols });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < j and i < self.rows) : (i += 1) {
                            result.data[i + j * self.size] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (j < int.min(self.rows, self.cols)) {
                            result.data[j * self.size + j] = self.data[j];
                        }

                        i = j + 1;
                        while (i < self.size) : (i += 1) {
                            result.data[i + j * self.size] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < i and j < self.cols) : (j += 1) {
                            result.data[i * self.size + j] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (i < int.min(self.rows, self.cols)) {
                            result.data[i * self.size + i] = self.data[i];
                        }

                        j = i + 1;
                        while (j < self.cols) : (j += 1) {
                            result.data[i * self.size + j] = constants.zero(T, ctx) catch unreachable;
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

    if (comptime !types.isDiagonalMatrix(@TypeOf(x))) {
        var result: R = try .init(allocator, y.rows, y.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < int.min(result.rows, result.cols)) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x, y.data[i]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x, y.data[i], ctx);
            }
        }

        return result;
    } else if (comptime !types.isDiagonalMatrix(@TypeOf(y))) {
        var result: R = try .init(allocator, x.rows, x.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < int.min(result.rows, result.cols)) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x.data[i], y);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x.data[i], y, ctx);
            }
        }

        return result;
    }

    if (x.rows != y.rows or x.cols != y.cols)
        return matrix.Error.DimensionMismatch;

    var result: R = try .init(allocator, x.rows, x.cols);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    var i: u32 = 0;
    while (i < int.min(result.rows, result.cols)) : (i += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            result.data[i] = op(x.data[i], y.data[i]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            result.data[i] = try op(x.data[i], y.data[i], ctx);
        }
    }

    return result;
}
