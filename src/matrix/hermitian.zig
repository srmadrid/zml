//! Storage scheme:
//!
//! Full `n`-by-`n` storage only accessing the upper or lower triangular part of
//! the matrix.

const std = @import("std");

const types = @import("../types.zig");
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const Coerce = types.Coerce;
const EnsureMatrix = types.EnsureMatrix;
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const General = matrix.General;
const Flags = matrix.Flags;

const array = @import("../array.zig");
const Dense = array.Dense;

pub fn Hermitian(comptime T: type) type {
    if (!types.isNumeric(T) or !types.isComplex(T))
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
            opts: struct {
                uplo: Uplo = .upper,
                order: Order = .col_major,
            },
        ) !Hermitian(T) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return Hermitian(T){
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
        ) !Hermitian(T) {
            const mat: Hermitian(T) = try .init(allocator, size, .{ .uplo = opts.uplo, .order = opts.order });
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
        ) !Hermitian(T) {
            const mat: Hermitian(T) = try .init(allocator, size, .{ .uplo = opts.uplo, .order = opts.order });
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

        pub fn toGeneral(self: Hermitian(T), allocator: std.mem.Allocator, ctx: anytype) !General(T) {
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
                                result.data[j + i * result.strides[1]] = ops.conjugate(self.data[i + j * self.strides[1]], ctx) catch unreachable;
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
                                result.data[j + i * result.strides[1]] = ops.conjugate(self.data[i + j * self.strides[1]], ctx) catch unreachable;
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
                                result.data[j * result.strides[0] + i] = ops.conjugate(self.data[i * self.strides[0] + j], ctx) catch unreachable;
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + j];
                                result.data[j * result.strides[0] + i] = ops.conjugate(self.data[i * self.strides[0] + j], ctx) catch unreachable;
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

        pub fn toDenseArray(self: *const Hermitian(T), allocator: std.mem.Allocator, ctx: anytype) !Dense(T) {
            var result: Dense(T) = try .init(allocator, &.{ self.size, self.size }, .{ .order = self.flags.order });
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
                                result.data[j + i * result.strides[1]] = ops.conjugate(self.data[i + j * self.strides[1]], ctx) catch unreachable;
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
                                result.data[j + i * result.strides[1]] = ops.conjugate(self.data[i + j * self.strides[1]], ctx) catch unreachable;
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
                                result.data[j * result.strides[0] + i] = ops.conjugate(self.data[i * self.strides[0] + j], ctx) catch unreachable;
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + j];
                                result.data[j * result.strides[0] + i] = ops.conjugate(self.data[i * self.strides[0] + j], ctx) catch unreachable;
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

        pub fn transpose(self: Hermitian(T)) Hermitian(T) {
            return Hermitian(T){
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

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    opts: if (types.isHermitianMatrix(Coerce(@TypeOf(x), @TypeOf(y))))
        struct {
            uplo: ?Uplo = null,
            order: ?Order = null,
        }
    else
        struct {
            order: ?Order = null,
        },
    ctx: anytype,
) !EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));

    if (comptime !types.isHermitianMatrix(@TypeOf(x))) {
        const ruplo: Uplo = if (comptime types.isHermitianMatrix(Coerce(@TypeOf(x), @TypeOf(y))))
            opts.uplo orelse y.uplo
        else
            y.uplo;
        var result: Coerce(@TypeOf(x), @TypeOf(y)) = if (comptime types.isHermitianMatrix(Coerce(@TypeOf(x), @TypeOf(y))))
            try .init( // Hermitian
                allocator,
                y.size,
                .{
                    .uplo = ruplo,
                    .order = opts.order orelse y.flags.order,
                },
            )
        else
            try .init( // General
                allocator,
                y.size,
                y.size,
                .{ .order = opts.order orelse y.flags.order },
            );
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (result.flags.order == .col_major) {
            if (y.flags.order == .col_major) {
                if (ruplo == .upper) {
                    if (y.uplo == .upper) { // cu cu
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(x, y.data[i + j * y.strides[1]]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(x, y.data[i + j * y.strides[1]], ctx);
                                }

                                if (comptime types.isComplex(X)) { // Result is a general matrix
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x, y.data[j + j * y.strides[1]]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x, y.data[j + j * y.strides[1]], ctx);
                            }
                        }
                    } else { // cu cl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(x, y.data[j + i * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(x, y.data[j + i * y.strides[1]], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x, y.data[j + j * y.strides[1]]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x, y.data[j + j * y.strides[1]], ctx);
                            }
                        }
                    }
                } else {
                    if (y.uplo == .upper) { // cl cu
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x, y.data[j + j * y.strides[1]]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x, y.data[j + j * y.strides[1]], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(x, y.data[j + i * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(x, y.data[j + i * y.strides[1]], ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl cl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x, y.data[j + j * y.strides[1]]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x, y.data[j + j * y.strides[1]], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(x, y.data[i + j * y.strides[1]]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(x, y.data[i + j * y.strides[1]], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (ruplo == .upper) {
                    if (y.uplo == .upper) { // cu ru
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(x, y.data[i * y.strides[0] + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(x, y.data[i * y.strides[0] + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.strides[0] + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.strides[0] + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x, y.data[j * y.strides[0] + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x, y.data[j * y.strides[0] + j], ctx);
                            }
                        }
                    } else { // cu rl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(x, y.data[j * y.strides[0] + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(x, y.data[j * y.strides[0] + i], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x, y.data[j * y.strides[0] + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x, y.data[j * y.strides[0] + j], ctx);
                            }
                        }
                    }
                } else {
                    if (y.uplo == .upper) { // cl ru
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x, y.data[j * y.strides[0] + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x, y.data[j * y.strides[0] + j], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(x, y.data[j * y.strides[0] + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(x, y.data[j * y.strides[0] + i], ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl rl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x, y.data[j * y.strides[0] + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x, y.data[j * y.strides[0] + j], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(x, y.data[i * y.strides[0] + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(x, y.data[i * y.strides[0] + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.strides[0] + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.strides[0] + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (y.flags.order == .col_major) {
                if (ruplo == .upper) {
                    if (y.uplo == .upper) { // ru cu
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x, y.data[i + i * y.strides[1]]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x, y.data[i + i * y.strides[1]], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(x, y.data[i + j * y.strides[1]]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(x, y.data[i + j * y.strides[1]], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else { // ru cl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x, y.data[i + i * y.strides[1]]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x, y.data[i + i * y.strides[1]], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(x, y.data[j + i * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(x, y.data[j + i * y.strides[1]], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (y.uplo == .upper) { // rl cu
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(x, y.data[j + i * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(x, y.data[j + i * y.strides[1]], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x, y.data[i + i * y.strides[1]]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x, y.data[i + i * y.strides[1]], ctx);
                            }
                        }
                    } else { // rl cl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(x, y.data[i + j * y.strides[1]]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(x, y.data[i + j * y.strides[1]], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x, y.data[i + i * y.strides[1]]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x, y.data[i + i * y.strides[1]], ctx);
                            }
                        }
                    }
                }
            } else {
                if (ruplo == .upper) {
                    if (y.uplo == .upper) { // ru ru
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x, y.data[i * y.strides[0] + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x, y.data[i * y.strides[0] + i], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(x, y.data[i * y.strides[0] + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(x, y.data[i * y.strides[0] + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.strides[0] + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.strides[0] + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else { // ru rl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x, y.data[i * y.strides[0] + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x, y.data[i * y.strides[0] + i], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(x, y.data[j * y.strides[0] + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(x, y.data[j * y.strides[0] + i], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (y.uplo == .upper) { // rl ru
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(x, y.data[j * y.strides[0] + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(x, y.data[j * y.strides[0] + i], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x, y.data[i * y.strides[0] + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x, y.data[i * y.strides[0] + i], ctx);
                            }
                        }
                    } else { // rl rl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(x, y.data[i * y.strides[0] + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(x, y.data[i * y.strides[0] + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.strides[0] + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.strides[0] + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x, y.data[i * y.strides[0] + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x, y.data[i * y.strides[0] + i], ctx);
                            }
                        }
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isHermitianMatrix(@TypeOf(y))) {
        const ruplo: Uplo = if (comptime types.isHermitianMatrix(Coerce(@TypeOf(x), @TypeOf(y))))
            opts.uplo orelse x.uplo
        else
            x.uplo;
        var result: Coerce(@TypeOf(x), @TypeOf(y)) = if (comptime types.isHermitianMatrix(Coerce(@TypeOf(x), @TypeOf(y))))
            try .init( // Hermitian
                allocator,
                x.size,
                .{
                    .uplo = ruplo,
                    .order = opts.order orelse x.flags.order,
                },
            )
        else
            try .init( // General
                allocator,
                x.size,
                x.size,
                .{ .order = opts.order orelse x.flags.order },
            );
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (result.flags.order == .col_major) {
            if (x.flags.order == .col_major) {
                if (ruplo == .upper) {
                    if (x.uplo == .upper) { // cu cu
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y, ctx);
                                }

                                if (comptime types.isComplex(X)) { // Result is a general matrix
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(
                                            ops.conjugate(x.data[i + j * x.strides[1]], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(
                                            ops.conjugate(x.data[i + j * x.strides[1]], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x.data[j + j * x.strides[1]], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x.data[j + j * x.strides[1]], y, ctx);
                            }
                        }
                    } else { // cu cl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(
                                        ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(
                                        ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(x.data[j + i * x.strides[1]], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(x.data[j + i * x.strides[1]], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x.data[j + j * x.strides[1]], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x.data[j + j * x.strides[1]], y, ctx);
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) { // cl cu
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x.data[j + j * x.strides[1]], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x.data[j + j * x.strides[1]], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(
                                        ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(
                                        ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(x.data[j + i * x.strides[1]], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(x.data[j + i * x.strides[1]], y, ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl cl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x.data[j + j * x.strides[1]], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x.data[j + j * x.strides[1]], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(
                                            ops.conjugate(x.data[i + j * x.strides[1]], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(
                                            ops.conjugate(x.data[i + j * x.strides[1]], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (ruplo == .upper) {
                    if (x.uplo == .upper) { // cu ru
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(
                                            ops.conjugate(x.data[i * x.strides[0] + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(
                                            ops.conjugate(x.data[i * x.strides[0] + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x.data[j * x.strides[0] + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x.data[j * x.strides[0] + j], y, ctx);
                            }
                        }
                    } else { // cu rl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(
                                        ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(
                                        ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(x.data[j * x.strides[0] + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(x.data[j * x.strides[0] + i], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x.data[j * x.strides[0] + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x.data[j * x.strides[0] + j], y, ctx);
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) { // cl ru
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x.data[j * x.strides[0] + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x.data[j * x.strides[0] + j], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(
                                        ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(
                                        ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(x.data[j * x.strides[0] + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(x.data[j * x.strides[0] + i], y, ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl rl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.strides[1]] = op(x.data[j * x.strides[0] + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.strides[1]] = try op(x.data[j * x.strides[0] + j], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = op(
                                            ops.conjugate(x.data[i * x.strides[0] + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = try op(
                                            ops.conjugate(x.data[i * x.strides[0] + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (x.flags.order == .col_major) {
                if (ruplo == .upper) {
                    if (x.uplo == .upper) { // ru cu
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x.data[i + i * x.strides[1]], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x.data[i + i * x.strides[1]], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(
                                            ops.conjugate(x.data[i + j * x.strides[1]], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(
                                            ops.conjugate(x.data[i + j * x.strides[1]], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else { // ru cl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x.data[i + i * x.strides[1]], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x.data[i + i * x.strides[1]], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(
                                        ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(
                                        ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(x.data[j + i * x.strides[1]], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(x.data[j + i * x.strides[1]], y, ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) { // rl cu
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(
                                        ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(
                                        ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(x.data[j + i * x.strides[1]], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(x.data[j + i * x.strides[1]], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x.data[i + i * x.strides[1]], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x.data[i + i * x.strides[1]], y, ctx);
                            }
                        }
                    } else { // rl cl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(
                                            ops.conjugate(x.data[i + j * x.strides[1]], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(
                                            ops.conjugate(x.data[i + j * x.strides[1]], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x.data[i + i * x.strides[1]], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x.data[i + i * x.strides[1]], y, ctx);
                            }
                        }
                    }
                }
            } else {
                if (ruplo == .upper) {
                    if (x.uplo == .upper) { // ru ru
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x.data[i * x.strides[0] + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x.data[i * x.strides[0] + i], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(x.data[i * x.strides[0] + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(x.data[i * x.strides[0] + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(
                                            ops.conjugate(x.data[i * x.strides[0] + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(
                                            ops.conjugate(x.data[i * x.strides[0] + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else { // ru rl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x.data[i * x.strides[0] + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x.data[i * x.strides[0] + i], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(
                                        ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(
                                        ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(x.data[j * x.strides[0] + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(x.data[j * x.strides[0] + i], y, ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) { // rl ru
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(
                                        ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(
                                        ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(x.data[j * x.strides[0] + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(x.data[j * x.strides[0] + i], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x.data[i * x.strides[0] + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x.data[i * x.strides[0] + i], y, ctx);
                            }
                        }
                    } else { // rl rl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.strides[0] + j] = op(x.data[i * x.strides[0] + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.strides[0] + j] = try op(x.data[i * x.strides[0] + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = op(
                                            ops.conjugate(x.data[i * x.strides[0] + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = try op(
                                            ops.conjugate(x.data[i * x.strides[0] + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.strides[0] + i] = op(x.data[i * x.strides[0] + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.strides[0] + i] = try op(x.data[i * x.strides[0] + i], y, ctx);
                            }
                        }
                    }
                }
            }
        }

        return result;
    }

    if (x.size != y.size)
        return matrix.Error.DimensionMismatch;

    var result: Hermitian(ReturnType2(op, X, Y)) = try .init(allocator, x.size, .{
        .uplo = opts.uplo orelse x.uplo.resolve2(y.uplo),
        .order = opts.order orelse x.flags.order.resolve2(y.flags.order),
    });
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (result.flags.order == .col_major) {
        if (x.flags.order == .col_major) {
            if (y.flags.order == .col_major) {
                if (result.uplo == .upper) {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // cu cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx);
                                    }
                                }
                            }
                        } else { // cu cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // cu cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cu cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // cl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // cl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cl cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (result.uplo == .upper) {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // cu cu ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx);
                                    }
                                }
                            }
                        } else { // cu cu rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // cu cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cu cl rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // cl cu ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl cu rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // cl cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cl cl rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (y.flags.order == .col_major) {
                if (result.uplo == .upper) {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // cu ru cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx);
                                    }
                                }
                            }
                        } else { // cu ru cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // cu rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cu rl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // cl ru cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl ru cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // cl rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cl rl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (result.uplo == .upper) {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // cu ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx);
                                    }
                                }
                            }
                        } else { // cu ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // cu rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cu rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // cl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.strides[1]] = ops.conjugate(try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // cl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cl rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        if (x.flags.order == .col_major) {
            if (y.flags.order == .col_major) {
                if (result.uplo == .upper) {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // ru cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx);
                                    }
                                }
                            }
                        } else { // ru cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // ru cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // ru cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // rl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // rl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // rl cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (result.uplo == .upper) {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // ru cu ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx);
                                    }
                                }
                            }
                        } else { // ru cu rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // ru cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // ru cl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // rl cu ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl cu rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            ops.conjugate(x.data[j + i * x.strides[1]], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // rl cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            x.data[i + j * x.strides[1]],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // rl cl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (y.flags.order == .col_major) {
                if (result.uplo == .upper) {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // ru ru cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx);
                                    }
                                }
                            }
                        } else { // ru ru cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // ru rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // ru rl cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // rl ru cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl ru cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i + j * y.strides[1]],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // rl rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j + i * y.strides[1]], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // rl rl cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (result.uplo == .upper) {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // ru ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx);
                                    }
                                }
                            }
                        } else { // ru ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // ru rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // ru rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (x.uplo == .upper) {
                        if (y.uplo == .upper) { // rl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.strides[0] + i] = ops.conjugate(try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            ops.conjugate(x.data[j * x.strides[0] + i], .{}) catch unreachable,
                                            y.data[i * y.strides[0] + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (y.uplo == .upper) { // rl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(
                                            x.data[i * x.strides[0] + j],
                                            ops.conjugate(y.data[j * y.strides[0] + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // rl rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.strides[0] + j] = op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.strides[0] + j] = try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return result;
}
