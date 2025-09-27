//! Storage scheme:
//!
//! Full `n`-by-`n` storage only accessing the upper or lower triangular part of
//! the matrix.

const std = @import("std");

const types = @import("../../types.zig");
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

const linalg = @import("../../linalg.zig");

pub fn Symmetric(T: type, uplo: Uplo, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("Symmetric requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        ld: u32, // leading dimension
        flags: Flags = .{},

        pub const empty: Symmetric(T, uplo, order) = .{
            .data = &.{},
            .size = 0,
            .ld = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            size: u32,
        ) !Symmetric(T, uplo, order) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(T, size * size)).ptr,
                .size = size,
                .ld = size,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            size: u32,
            value: anytype,
            ctx: anytype,
        ) !Symmetric(T, uplo, order) {
            const mat: Symmetric(T, uplo, order) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                mat.data[i + j * size] = value_casted;
                            }
                        }
                    } else { // cl
                        var j: u32 = 0;
                        while (j < size) : (j += 1) {
                            var i: u32 = j;
                            while (i < size) : (i += 1) {
                                mat.data[i + j * size] = value_casted;
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < size) : (i += 1) {
                            var j: u32 = i;
                            while (j < size) : (j += 1) {
                                mat.data[i * size + j] = value_casted;
                            }
                        }
                    } else { // rl
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
            ctx: anytype,
        ) !Symmetric(T, uplo, order) {
            const mat: Symmetric(T, uplo, order) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                            }

                            mat.data[j + j * size] = constants.one(T, ctx) catch unreachable;
                        }
                    } else { // cl
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
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < size) : (i += 1) {
                            mat.data[i * size + i] = constants.one(T, ctx) catch unreachable;

                            var j: u32 = i + 1;
                            while (j < size) : (j += 1) {
                                mat.data[i * size + j] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    } else { // rl
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

        pub fn deinit(self: *Symmetric(T, uplo, order), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. self.size * self.size]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Symmetric(T, uplo, order), row: u32, col: u32) !T {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = row;
            var j: u32 = col;
            if (comptime uplo == .upper) {
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

            return if (comptime order == .col_major)
                self.data[i + j * self.ld]
            else
                self.data[i * self.ld + j];
        }

        pub inline fn at(self: *const Symmetric(T, uplo, order), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and on
            // the correct triangular part.
            return if (comptime order == .col_major)
                self.data[row + col * self.ld]
            else
                self.data[row * self.ld + col];
        }

        pub fn set(self: *Symmetric(T, uplo, order), row: u32, col: u32, value: T) !void {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = row;
            var j: u32 = col;
            if (comptime uplo == .upper) {
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

            if (comptime order == .col_major) {
                self.data[i + j * self.ld] = value;
            } else {
                self.data[i * self.ld + j] = value;
            }
        }

        pub inline fn put(self: *Symmetric(T, uplo, order), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid and on
            // the correct triangular part.
            if (comptime order == .col_major) {
                self.data[row + col * self.ld] = value;
            } else {
                self.data[row * self.ld + col] = value;
            }
        }

        pub fn copy(self: *const Symmetric(T, uplo, order), allocator: std.mem.Allocator, ctx: anytype) !Symmetric(T, uplo, order) {
            var mat: Symmetric(T, uplo, order) = try .init(allocator, self.size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            try linalg.blas.copy(
                                types.scast(i32, j + 1),
                                self.data + (0 + j * self.ld),
                                1,
                                mat.data + (0 + j * mat.ld),
                                1,
                                ctx,
                            );
                        }
                    } else { // cl
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            try linalg.blas.copy(
                                types.scast(i32, self.size - j),
                                self.data + (j + j * self.ld),
                                1,
                                mat.data + (j + j * mat.ld),
                                1,
                                ctx,
                            );
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            try linalg.blas.copy(
                                types.scast(i32, self.size - i),
                                self.data + (i * self.ld + i),
                                1,
                                mat.data + (i * mat.ld + i),
                                1,
                                ctx,
                            );
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            try linalg.blas.copy(
                                types.scast(i32, i + 1),
                                self.data + (i * self.ld + 0),
                                1,
                                mat.data + (i * mat.ld + 0),
                                1,
                                ctx,
                            );
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn copyInverseUplo(
            self: *const Symmetric(T, uplo, order),
            allocator: std.mem.Allocator,
            ctx: anytype,
        ) !Symmetric(T, uplo.invert(), order) {
            var mat: Symmetric(T, uplo.invert(), order) = try .init(allocator, self.size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu -> cl
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            try linalg.blas.copy(
                                types.scast(i32, j + 1),
                                self.data + (0 + j * self.ld),
                                1,
                                mat.data + j,
                                types.scast(i32, mat.ld),
                                ctx,
                            );
                        }
                    } else { // cl -> cu
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            try linalg.blas.copy(
                                types.scast(i32, self.size - j),
                                self.data + (j + j * self.ld),
                                1,
                                mat.data + (j + j * self.ld),
                                types.scast(i32, mat.ld),
                                ctx,
                            );
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru -> rl
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            try linalg.blas.copy(
                                types.scast(i32, self.size - i),
                                self.data + (i * self.ld + i),
                                1,
                                mat.data + (i * mat.ld + i),
                                types.scast(i32, mat.ld),
                                ctx,
                            );
                        }
                    } else { // rl -> ru
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            try linalg.blas.copy(
                                types.scast(i32, i + 1),
                                self.data + (i * self.ld + 0),
                                1,
                                mat.data + i,
                                types.scast(i32, mat.ld),
                                ctx,
                            );
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn toGeneralDenseMatrix(self: Symmetric(T, uplo, order), allocator: std.mem.Allocator, ctx: anytype) !matrix.dense.General(T, order) {
            var result: matrix.dense.General(T, order) = try .init(allocator, self.size, self.size);
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                result.data[i + j * result.ld] = self.data[i + j * self.ld];
                                result.data[j + i * result.ld] = self.data[i + j * self.ld];
                            }

                            result.data[j + j * result.ld] = self.data[j + j * self.ld];
                        }
                    } else { // cl
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[j + j * result.ld] = self.data[j + j * self.ld];

                            var i: u32 = j + 1;
                            while (i < self.size) : (i += 1) {
                                result.data[i + j * result.ld] = self.data[i + j * self.ld];
                                result.data[j + i * result.ld] = self.data[i + j * self.ld];
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i * result.ld + i] = self.data[i * self.ld + i];

                            var j: u32 = i + 1;
                            while (j < self.size) : (j += 1) {
                                result.data[i * result.ld + j] = self.data[i * self.ld + j];
                                result.data[j * result.ld + i] = self.data[i * self.ld + j];
                            }
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.ld + j] = self.data[i * self.ld + j];
                                result.data[j * result.ld + i] = self.data[i * self.ld + j];
                            }

                            result.data[i * result.ld + i] = self.data[i * self.ld + i];
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn toDenseArray(self: *const Symmetric(T, uplo, order), allocator: std.mem.Allocator, ctx: anytype) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.size, self.size });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                                result.data[j + i * result.strides[1]] = self.data[i + j * self.ld];
                            }

                            result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];
                        }
                    } else { // cl
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];

                            var i: u32 = j + 1;
                            while (i < self.size) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                                result.data[j + i * result.strides[1]] = self.data[i + j * self.ld];
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];

                            var j: u32 = i + 1;
                            while (j < self.size) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                                result.data[j * result.strides[0] + i] = self.data[i * self.ld + j];
                            }
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                                result.data[j * result.strides[0] + i] = self.data[i * self.ld + j];
                            }

                            result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn transpose(self: Symmetric(T, uplo, order)) Symmetric(T, uplo.invert(), order.invert()) {
            return .{
                .data = self.data,
                .size = self.size,
                .ld = self.ld,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Symmetric(T, uplo, order),
            start: u32,
            end: u32,
        ) !Symmetric(T, uplo, order) {
            if (start >= self.size or end > self.size or start >= end)
                return matrix.Error.InvalidRange;

            const sub_size: u32 = end - start;

            return .{
                .data = self.data + (start + start * self.ld),
                .size = sub_size,
                .ld = self.ld,
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

    if (comptime !types.isSymmetricMatrix(@TypeOf(x))) {
        var result: R = try .init(allocator, y.size);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[j + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[j + i * y.ld], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[j + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[j + i * y.ld], ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[j * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[j * y.ld + i], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[j * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[j * y.ld + i], ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[j + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[j + i * y.ld], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[j + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[j + i * y.ld], ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[j * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[j * y.ld + i], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[j * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[j * y.ld + i], ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isSymmetricMatrix(@TypeOf(y))) {
        var result: R = try .init(allocator, x.size);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y, ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y, ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y, ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y, ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                                }
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

    var result: R = try .init(allocator, x.size);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cu cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cu cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cl cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu cu ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cu cu rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cu cl rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl cu ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cl cu rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        } else { // cl cl rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu ru cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cu ru cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cu rl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl ru cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cl ru cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cl rl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cu ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cu rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        } else { // cl rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // ru cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // ru cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // rl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // rl cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru cu ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // ru cu rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // ru cl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl cu ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // rl cu rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        } else { // rl cl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru ru cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // ru ru cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // ru rl cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl ru cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // rl ru cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // rl rl cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // ru ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // ru rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // rl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        } else { // rl rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
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
