//! Storage scheme:
//!
//! Full `n`-by-`n` storage only accessing the upper or lower triangular part of
//! the matrix.

const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const General = matrix.General;
const Flags = matrix.Flags;

pub fn Symmetric(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("Symmetric requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        strides: [2]u32,
        uplo: Uplo,
        flags: Flags = .{},

        pub const empty: Symmetric(T) = .{
            .data = &.{},
            .size = 0,
            .strides = .{ 0, 0 },
            .uplo = .upper,
            .flags = .{ .order = .col_major, .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            size: u32,
            opts: struct {
                uplo: Uplo = .upper,
                order: Order = .col_major,
            },
        ) !Symmetric(T) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return Symmetric(T){
                .data = (try allocator.alloc(T, size * size)).ptr,
                .size = size,
                .strides = if (opts.order == .col_major) .{ 1, size } else .{ size, 1 },
                .uplo = opts.uplo,
                .flags = .{ .order = opts.order, .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            size: u32,
            value: anytype,
            opts: struct {
                uplo: Uplo = .upper,
                order: Order = .col_major,
            },
            ctx: anytype,
        ) !Symmetric(T) {
            const mat: Symmetric(T) = try .init(allocator, size, opts.uplo, .{ .order = opts.order });
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                if (opts.order == .col_major) {
                    if (opts.uplo == .upper) {
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
                    if (opts.uplo == .upper) {
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
            opts: struct {
                uplo: Uplo = .upper,
                order: Order = .col_major,
            },
            ctx: anytype,
        ) !Symmetric(T) {
            const mat: Symmetric(T) = try .init(allocator, size, opts.uplo, .{ .order = opts.order });
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (opts.order == .col_major) {
                    if (opts.uplo == .upper) {
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
                    if (opts.uplo == .upper) {
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

        pub fn deinit(self: *Symmetric(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. self.size * self.size]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Symmetric(T), row: u32, col: u32) !T {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = row;
            var j: u32 = col;
            if (self.uplo == .upper) {
                if (i > j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                }
            } else {
                if (i < j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                }
            }

            return self.data[i * self.strides[0] + j * self.strides[1]];
        }

        pub inline fn at(self: *const Symmetric(T), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and on
            // the correct triangular part.
            return self.data[row * self.strides[0] + col * self.strides[1]];
        }

        pub fn set(self: *Symmetric(T), row: u32, col: u32, value: T) !void {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = row;
            var j: u32 = col;
            if (self.uplo == .upper) {
                if (i > j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                }
            } else {
                if (i < j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                }
            }

            self.data[i * self.strides[0] + j * self.strides[1]] = value;
        }

        pub inline fn put(self: *Symmetric(T), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid and on
            // the correct triangular part.
            self.data[row * self.strides[0] + col * self.strides[1]] = value;
        }

        pub fn toGeneral(self: Symmetric(T), allocator: std.mem.Allocator, ctx: anytype) !General(T) {
            var result: General(T) = try .init(allocator, self.size, self.size, .{ .order = self.flags.order });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (self.flags.order == .col_major) {
                    if (self.uplo == .upper) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.strides[1]];
                                result.data[j + i * result.strides[1]] = self.data[i + j * self.strides[1]];
                            }

                            result.data[j + j * result.strides[1]] = self.data[j + j * self.strides[1]];
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[j + j * result.strides[1]] = self.data[j + j * self.strides[1]];

                            var i: u32 = j + 1;
                            while (i < self.size) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.strides[1]];
                                result.data[j + i * result.strides[1]] = self.data[i + j * self.strides[1]];
                            }
                        }
                    }
                } else {
                    if (self.uplo == .upper) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i * result.strides[0] + i] = self.data[i * self.strides[0] + i];

                            var j: u32 = i + 1;
                            while (j < self.size) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + j];
                                result.data[j * result.strides[0] + i] = self.data[i * self.strides[0] + j];
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + j];
                                result.data[j * result.strides[0] + i] = self.data[i * self.strides[0] + j];
                            }

                            result.data[i * result.strides[0] + i] = self.data[i * self.strides[0] + i];
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn transpose(self: Symmetric(T)) Symmetric(T) {
            return Symmetric(T){
                .data = self.data,
                .size = self.size,
                .strides = .{ self.strides[1], self.strides[0] },
                .uplo = self.uplo.invert(),
                .flags = .{
                    .order = self.flags.order.invert(),
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Symmetric(T),
            start: u32,
            end: u32,
        ) !Symmetric(T) {
            if (start >= self.size or end > self.size or start >= end)
                return matrix.Error.InvalidRange;

            const sub_size = end - start;

            return Symmetric(T){
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
