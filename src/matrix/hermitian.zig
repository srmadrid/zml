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

const linalg = @import("../linalg.zig");

pub fn Hermitian(T: type, uplo: Uplo, order: Order) type {
    if (!types.isNumeric(T) or !types.isComplex(T))
        @compileError("Hermitian requires a complex numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        ld: u32, // leading dimension
        flags: Flags = .{},

        pub const empty: Hermitian(T, uplo, order) = .{
            .data = &.{},
            .size = 0,
            .ld = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            size: u32,
        ) !Hermitian(T, uplo, order) {
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
        ) !Hermitian(T, uplo, order) {
            const mat: Hermitian(T, uplo, order) = try .init(allocator, size);
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
        ) !Hermitian(T, uplo, order) {
            const mat: Hermitian(T, uplo, order) = try .init(allocator, size);
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

        pub fn deinit(self: *Hermitian(T, uplo, order), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. self.size * self.size]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Hermitian(T, uplo, order), row: u32, col: u32) !T {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = row;
            var j: u32 = col;
            var noconj: bool = true;
            if (comptime uplo == .upper) {
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

            if (noconj)
                return if (comptime order == .col_major)
                    self.data[i + j * self.ld]
                else
                    self.data[i * self.ld + j]
            else
                return if (comptime order == .col_major)
                    ops.conjugate(self.data[i + j * self.ld], .{}) catch unreachable
                else
                    ops.conjugate(self.data[i * self.ld + j], .{}) catch unreachable;
        }

        pub inline fn at(self: *const Hermitian(T, uplo, order), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and on
            // the correct triangular part.
            return if (comptime order == .col_major)
                self.data[row + col * self.ld]
            else
                self.data[row * self.ld + col];
        }

        pub fn set(self: *Hermitian(T, uplo, order), row: u32, col: u32, value: T) !void {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            if (row == col and value.im != 0)
                return matrix.Error.BreaksStructure;

            var i: u32 = row;
            var j: u32 = col;
            var noconj: bool = true;
            if (comptime uplo == .upper) {
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

            if (comptime order == .col_major) {
                self.data[i + j * self.ld] = if (noconj) {
                    value;
                } else {
                    ops.conjugate(value, .{}) catch unreachable;
                };
            } else {
                self.data[i * self.ld + j] = if (noconj) {
                    value;
                } else {
                    ops.conjugate(value, .{}) catch unreachable;
                };
            }
        }

        pub inline fn put(self: *Hermitian(T, uplo, order), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid and on
            // the correct triangular part.
            if (comptime order == .col_major) {
                self.data[row + col * self.ld] = value;
            } else {
                self.data[row * self.ld + col] = value;
            }
        }

        pub fn copy(self: *const Hermitian(T, uplo, order), allocator: std.mem.Allocator, ctx: anytype) !Hermitian(T, uplo, order) {
            var mat: Hermitian(T, uplo, order) = try .init(allocator, self.size);
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
            self: *const Hermitian(T, uplo, order),
            allocator: std.mem.Allocator,
            ctx: anytype,
        ) !Hermitian(T, uplo.invert(), order) {
            var mat: Hermitian(T, uplo.invert(), order) = try .init(allocator, self.size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu -> cl
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                try ops.set(
                                    &mat.data[j + i * mat.ld],
                                    ops.conjugate(self.data[i + j * self.ld], ctx) catch unreachable,
                                    ctx,
                                );
                            }

                            ops.set(
                                &mat.data[j + j * mat.ld],
                                self.data[j + j * self.ld].re,
                                ctx,
                            ) catch unreachable;
                        }
                    } else { // cl -> cu
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            ops.set(
                                &mat.data[j + j * mat.ld],
                                self.data[j + j * self.ld].re,
                                ctx,
                            ) catch unreachable;

                            var i: u32 = j + 1;
                            while (i < self.size) : (i += 1) {
                                try ops.set(
                                    &mat.data[j + i * mat.ld],
                                    ops.conjugate(self.data[i + j * self.ld], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru -> rl
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            ops.set(
                                &mat.data[i * mat.ld + i],
                                self.data[i * self.ld + i].re,
                                ctx,
                            ) catch unreachable;

                            var j: u32 = i + 1;
                            while (j < self.size) : (j += 1) {
                                try ops.set(
                                    &mat.data[j * mat.ld + i],
                                    ops.conjugate(self.data[i * self.ld + j], ctx) catch unreachable,
                                    ctx,
                                );
                            }
                        }
                    } else { // rl -> ru
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                try ops.set(
                                    &mat.data[j * mat.ld + i],
                                    ops.conjugate(self.data[i * self.ld + j], ctx) catch unreachable,
                                    ctx,
                                );
                            }

                            ops.set(
                                &mat.data[i * mat.ld + i],
                                self.data[i * self.ld + i].re,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn toGeneral(self: Hermitian(T, uplo, order), allocator: std.mem.Allocator, ctx: anytype) !General(T, order) {
            var result: General(T, order) = try .init(allocator, self.size, self.size);
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
                                result.data[j + i * result.ld] = ops.conjugate(self.data[i + j * self.ld], ctx) catch unreachable;
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
                                result.data[j + i * result.ld] = ops.conjugate(self.data[i + j * self.ld], ctx) catch unreachable;
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
                                result.data[j * result.ld + i] = ops.conjugate(self.data[i * self.ld + j], ctx) catch unreachable;
                            }
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.ld + j] = self.data[i * self.ld + j];
                                result.data[j * result.ld + i] = ops.conjugate(self.data[i * self.ld + j], ctx) catch unreachable;
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

        pub fn toDenseArray(self: *const Hermitian(T, uplo, order), allocator: std.mem.Allocator, ctx: anytype) !Dense(T, order) {
            var result: Dense(T, order) = try .init(allocator, &.{ self.size, self.size });
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
                                result.data[j + i * result.strides[1]] = ops.conjugate(self.data[i + j * self.ld], ctx) catch unreachable;
                            }

                            result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];

                            var i: u32 = j + 1;
                            while (i < self.size) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                                result.data[j + i * result.strides[1]] = ops.conjugate(self.data[i + j * self.ld], ctx) catch unreachable;
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
                                result.data[j * result.strides[0] + i] = ops.conjugate(self.data[i * self.ld + j], ctx) catch unreachable;
                            }
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                                result.data[j * result.strides[0] + i] = ops.conjugate(self.data[i * self.ld + j], ctx) catch unreachable;
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

        pub fn transpose(self: Hermitian(T, uplo, order)) Hermitian(T, uplo.invert(), order.invert()) {
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
            self: *const Hermitian(T, uplo, order),
            start: u32,
            end: u32,
        ) !Hermitian(T, uplo, order) {
            if (start >= self.size or end > self.size or start >= end)
                return matrix.Error.InvalidRange;

            const sub_size = end - start;

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

    if (comptime !types.isHermitianMatrix(@TypeOf(x))) {
        const ruplo: Uplo = comptime types.uploOf(@TypeOf(y));
        var result: R = if (comptime types.isHermitianMatrix(R))
            try .init( // Hermitian
                allocator,
                y.size,
            )
        else
            try .init( // General
                allocator,
                y.size,
                y.size,
            );
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cu
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                                }

                                if (comptime types.isComplex(X)) { // Result is a general matrix
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                            }
                        }
                    } else { // cu cl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x, y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x, y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cu
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x, y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x, y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl cl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu ru
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                            }
                        }
                    } else { // cu rl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x, y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x, y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl ru
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x, y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x, y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl rl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
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
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cu
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
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
                                result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x, y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x, y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cu
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x, y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x, y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                            }
                        }
                    } else { // rl cl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru ru
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
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
                                result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x, y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x, y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl ru
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x, y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x, y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                            }
                        }
                    } else { // rl rl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                            }
                        }
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isHermitianMatrix(@TypeOf(y))) {
        const ruplo: Uplo = comptime types.uploOf(@TypeOf(x));
        var result: R = if (comptime types.isHermitianMatrix(R))
            try .init( // Hermitian
                allocator,
                x.size,
            )
        else
            try .init( // General
                allocator,
                x.size,
                x.size,
            );
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cu cu
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                                }

                                if (comptime types.isComplex(X)) { // Result is a general matrix
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                            }
                        }
                    } else { // cu cl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[j + i * x.ld], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[j + i * x.ld], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cl cu
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[j + i * x.ld], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[j + i * x.ld], y, ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl cl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
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
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cu ru
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                            }
                        }
                    } else { // cu rl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[j * x.ld + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[j * x.ld + i], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cl ru
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[j * x.ld + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[j * x.ld + i], y, ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl rl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
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
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // ru cu
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
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
                                result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[j + i * x.ld], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[j + i * x.ld], y, ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // rl cu
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[j + i * x.ld], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[j + i * x.ld], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                            }
                        }
                    } else { // rl cl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // ru ru
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
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
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[j * x.ld + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[j * x.ld + i], y, ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // rl ru
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[j * x.ld + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[j * x.ld + i], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                            }
                        }
                    } else { // rl rl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
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
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cu cu
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
                                        result.data[i + j * result.ld] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
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
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
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
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cu ru
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
                                        result.data[i + j * result.ld] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
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
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cu ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl cu rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
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
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu ru cu
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
                                        result.data[i + j * result.ld] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
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
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl ru cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl ru cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
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
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu ru ru
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
                                        result.data[i + j * result.ld] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
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
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
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
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cu cu
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
                                        result.data[i * result.ld + j] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
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
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
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
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cu ru
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
                                        result.data[i * result.ld + j] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
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
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cu ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl cu rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
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
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru ru cu
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
                                        result.data[i * result.ld + j] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
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
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl ru cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl ru cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
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
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru ru ru
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
                                        result.data[i * result.ld + j] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
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
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
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
