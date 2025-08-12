//! Storage scheme:
//!
//! Full `n`-by-`n` storage only accessing the upper or lower triangular part of
//! the matrix.

const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const Flags = matrix.Flags;

pub fn Hermitian(comptime T: type) type {
    if (!types.isNumeric(T) and !types.isComplex(T))
        @compileError("Hermitian requires a complex numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        strides: [2]u32,
        uplo: Uplo,
        flags: Flags = .{},

        pub const empty: Hermitian(T) = .{
            .data = &.{},
            .size = 0,
            .strides = .{ 0, 0 },
            .uplo = .upper,
            .flags = .{ .order = .col_major, .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            size: u32,
            uplo: Uplo,
            order: Order,
        ) !Hermitian(T) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return Hermitian(T){
                .data = (try allocator.alloc(T, size * size)).ptr,
                .size = size,
                .strides = if (order == .col_major) .{ 1, size } else .{ size, 1 },
                .uplo = uplo,
                .flags = .{ .order = order, .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            size: u32,
            uplo: Uplo,
            value: anytype,
            order: Order,
            ctx: anytype,
        ) !Hermitian(T) {
            const mat: Hermitian(T) = try .init(allocator, size, uplo, order);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                if (order == .col_major) {
                    if (uplo == .upper) {
                        var j: u32 = 0;
                        while (j < size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                mat.data[i + j * size] = value_casted;
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < size) : (j += 1) {
                            var i: u32 = j;
                            while (i < size) : (i += 1) {
                                mat.data[i + j * size] = value_casted;
                            }
                        }
                    }
                } else {
                    if (uplo == .upper) {
                        var i: u32 = 0;
                        while (i < size) : (i += 1) {
                            var j: u32 = i;
                            while (j < size) : (j += 1) {
                                mat.data[i * size + j] = value_casted;
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                mat.data[i * size + j] = value_casted;
                            }
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn eye(
            allocator: std.mem.Allocator,
            size: u32,
            uplo: Uplo,
            order: Order,
            ctx: anytype,
        ) !Hermitian(T) {
            const mat: Hermitian(T) = try .init(allocator, size, uplo, order);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (order == .col_major) {
                    if (uplo == .upper) {
                        var j: u32 = 0;
                        while (j < size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                            }

                            mat.data[j + j * size] = constants.one(T, ctx) catch unreachable;
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < size) : (j += 1) {
                            mat.data[j + j * size] = constants.one(T, ctx) catch unreachable;

                            var i: u32 = j + 1;
                            while (i < size) : (i += 1) {
                                mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    }
                } else {
                    if (uplo == .upper) {
                        var i: u32 = 0;
                        while (i < size) : (i += 1) {
                            mat.data[i * size + i] = constants.one(T, ctx) catch unreachable;

                            var j: u32 = i + 1;
                            while (j < size) : (j += 1) {
                                mat.data[i * size + j] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                mat.data[i * size + j] = constants.zero(T, ctx) catch unreachable;
                            }

                            mat.data[i * size + i] = constants.one(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn deinit(self: *Hermitian(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. self.size * self.size]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Hermitian(T), row: u32, col: u32) !T {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = row;
            var j: u32 = col;
            var noconj: bool = true;
            if (self.uplo == .upper) {
                if (i > j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                    noconj = false;
                }
            } else {
                if (i < j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                    noconj = false;
                }
            }

            return if (noconj)
                self.data[i * self.strides[0] + j * self.strides[1]]
            else
                ops.conjugate(self.data[i * self.strides[0] + j * self.strides[1]], .{}) catch unreachable;
        }

        pub inline fn at(self: *const Hermitian(T), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and on
            // the correct triangular part.
            return self.data[row * self.strides[0] + col * self.strides[1]];
        }

        pub fn set(self: *Hermitian(T), row: u32, col: u32, value: T) !void {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            if (row == col and value.im != 0)
                return matrix.Error.BreaksStructure;

            var i: u32 = row;
            var j: u32 = col;
            var noconj: bool = true;
            if (self.uplo == .upper) {
                if (i > j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                    noconj = false;
                }
            } else {
                if (i < j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                    noconj = false;
                }
            }

            self.data[i * self.strides[0] + j * self.strides[1]] = if (noconj) {
                value;
            } else {
                ops.conjugate(value, .{}) catch unreachable;
            };
        }

        pub inline fn put(self: *Hermitian(T), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid and on
            // the correct triangular part.
            self.data[row * self.strides[0] + col * self.strides[1]] = value;
        }

        pub fn transpose(self: Hermitian(T)) Hermitian(T) {
            return Hermitian(T){
                .data = self.data,
                .size = self.size,
                .strides = .{ self.strides[1], self.strides[0] },
                .uplo = self.uplo.invert(),
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Hermitian(T),
            start: u32,
            end: u32,
        ) !Hermitian(T) {
            if (start >= self.size or end > self.size or start >= end)
                return matrix.Error.InvalidRange;

            const sub_size = end - start;

            return Hermitian(T){
                .data = self.data + (start * self.strides[0] + start * self.strides[1]),
                .size = sub_size,
                .strides = self.strides,
                .uplo = self.uplo,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }
    };
}
