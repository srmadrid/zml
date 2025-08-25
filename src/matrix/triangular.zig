//! Storage scheme:
//!
//! Full `n`-by-`n` storage only accessing the upper or lower triangular part of
//! the matrix, with an option for unit triangular matrices.

const std = @import("std");

const types = @import("../types.zig");
const Numeric = types.Numeric;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Order = types.Order;
const Uplo = types.Uplo;
const Diag = types.Diag;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const General = matrix.General;
const Flags = matrix.Flags;

const array = @import("../array.zig");
const Dense = array.Dense;

pub fn Triangular(T: type, uplo: Uplo, diag: Diag, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("Triangular requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        ld: u32, // leading dimension
        flags: Flags = .{},

        pub const empty: Triangular(T, uplo, diag, order) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .ld = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
        ) !Triangular(T, uplo, diag, order) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(T, rows * cols)).ptr,
                .rows = rows,
                .cols = cols,
                .ld = if (order == .col_major) rows else cols,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            value: anytype,
            ctx: anytype,
        ) !Triangular(T, uplo, diag, order) {
            const mat: Triangular(T, uplo, diag, order) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) {
                        if (comptime diag == .unit) { // cuu
                            var j: u32 = 0;
                            while (j < cols) : (j += 1) {
                                var i: u32 = 0;
                                while (i < int.min(j, rows)) : (i += 1) {
                                    mat.data[i + j * rows] = value_casted;
                                }
                            }
                        } else { // cun
                            var j: u32 = 0;
                            while (j < cols) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= int.min(j, rows - 1)) : (i += 1) {
                                    mat.data[i + j * rows] = value_casted;
                                }
                            }
                        }
                    } else {
                        if (comptime diag == .unit) { // clu
                            var j: u32 = 0;
                            while (j < cols) : (j += 1) {
                                var i: u32 = j + 1;
                                while (i < rows) : (i += 1) {
                                    mat.data[i + j * rows] = value_casted;
                                }
                            }
                        } else { // cln
                            var j: u32 = 0;
                            while (j < cols) : (j += 1) {
                                var i: u32 = j;
                                while (i < rows) : (i += 1) {
                                    mat.data[i + j * rows] = value_casted;
                                }
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) {
                        if (comptime diag == .unit) { // ruu
                            var i: u32 = 0;
                            while (i < rows) : (i += 1) {
                                var j: u32 = i + 1;
                                while (j < cols) : (j += 1) {
                                    mat.data[i * cols + j] = value_casted;
                                }
                            }
                        } else { // run
                            var i: u32 = 0;
                            while (i < rows) : (i += 1) {
                                var j: u32 = i;
                                while (j < cols) : (j += 1) {
                                    mat.data[i * cols + j] = value_casted;
                                }
                            }
                        }
                    } else {
                        if (comptime diag == .unit) { // rlu
                            var i: u32 = 0;
                            while (i < rows) : (i += 1) {
                                var j: u32 = 0;
                                while (j < int.min(i, cols)) : (j += 1) {
                                    mat.data[i * cols + j] = value_casted;
                                }
                            }
                        } else { // rln
                            var i: u32 = 0;
                            while (i < rows) : (i += 1) {
                                var j: u32 = 0;
                                while (j <= int.min(i, cols - 1)) : (j += 1) {
                                    mat.data[i * cols + j] = value_casted;
                                }
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
        ) !Triangular(T, uplo, diag, order) {
            const mat: Triangular(T, uplo, diag, order) = try .init(allocator, size, size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) {
                        if (comptime diag == .unit) { // cuu
                            var j: u32 = 0;
                            while (j < size) : (j += 1) {
                                var i: u32 = 0;
                                while (i < j) : (i += 1) {
                                    mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                                }
                            }
                        } else { // cun
                            var j: u32 = 0;
                            while (j < size) : (j += 1) {
                                var i: u32 = 0;
                                while (i < j) : (i += 1) {
                                    mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                                }

                                mat.data[j + j * size] = constants.one(T, ctx) catch unreachable;
                            }
                        }
                    } else {
                        if (comptime diag == .unit) { // clu
                            var j: u32 = 0;
                            while (j < size) : (j += 1) {
                                var i: u32 = j + 1;
                                while (i < size) : (i += 1) {
                                    mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                                }
                            }
                        } else { // cln
                            var j: u32 = 0;
                            while (j < size) : (j += 1) {
                                mat.data[j + j * size] = constants.one(T, ctx) catch unreachable;

                                var i: u32 = j + 1;
                                while (i < size) : (i += 1) {
                                    mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) {
                        if (comptime diag == .unit) { // ruu
                            var i: u32 = 0;
                            while (i < size) : (i += 1) {
                                var j: u32 = i + 1;
                                while (j < size) : (j += 1) {
                                    mat.data[i * size + j] = constants.zero(T, ctx) catch unreachable;
                                }
                            }
                        } else { // run
                            var i: u32 = 0;
                            while (i < size) : (i += 1) {
                                mat.data[i * size + i] = constants.one(T, ctx) catch unreachable;

                                var j: u32 = i + 1;
                                while (j < size) : (j += 1) {
                                    mat.data[i * size + j] = constants.zero(T, ctx) catch unreachable;
                                }
                            }
                        }
                    } else {
                        if (comptime diag == .unit) { // rlu
                            var i: u32 = 0;
                            while (i < size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < i) : (j += 1) {
                                    mat.data[i * size + j] = constants.zero(T, ctx) catch unreachable;
                                }
                            }
                        } else { // rln
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
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn deinit(self: *Triangular(T, uplo, diag, order), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Triangular(T, uplo, diag, order), row: u32, col: u32) !T {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (row > col)
                    return constants.zero(T, .{}) catch unreachable;
            } else {
                if (row < col)
                    return constants.zero(T, .{}) catch unreachable;
            }

            if (comptime diag == .unit) {
                if (row == col)
                    return constants.one(T, .{}) catch unreachable;
            }

            return if (comptime order == .col_major)
                self.data[row + col * self.ld]
            else
                self.data[row * self.ld + col];
        }

        pub inline fn at(self: *const Triangular(T, uplo, diag, order), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and on
            // the correct triangular part, and outside the diagonal if diag
            // triangular.
            return if (comptime order == .col_major)
                self.data[row + col * self.ld]
            else
                self.data[row * self.ld + col];
        }

        pub fn set(self: *Triangular(T, uplo, diag, order), row: u32, col: u32, value: T) !void {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (row > col)
                    return matrix.Error.PositionOutOfBounds;
            } else {
                if (row < col)
                    return matrix.Error.PositionOutOfBounds;
            }

            if (comptime diag == .unit) {
                if (row == col)
                    return matrix.Error.PositionOutOfBounds;
            }

            if (comptime order == .col_major) {
                self.data[row + col * self.ld] = value;
            } else {
                self.data[row * self.ld + col] = value;
            }
        }

        pub inline fn put(self: *Triangular(T, uplo, diag, order), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid and on
            // the correct triangular part, and outside the diagonal if diag
            // triangular.
            if (comptime order == .col_major) {
                self.data[row + col * self.ld] = value;
            } else {
                self.data[row * self.ld + col] = value;
            }
        }

        pub fn toGeneral(self: Triangular(T, uplo, diag, order), allocator: std.mem.Allocator, ctx: anytype) !General(T, order) {
            var result: General(T, order) = try .init(allocator, self.rows, self.cols);
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, self.rows)) : (i += 1) {
                                result.data[i + j * result.ld] = self.data[i + j * self.ld];
                            }

                            if (j < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[j + j * result.ld] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[j + j * result.ld] = self.data[j + j * self.ld];
                                }
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.ld] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    } else { // cl
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, self.rows)) : (i += 1) {
                                result.data[i + j * result.ld] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (j < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[j + j * result.ld] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[j + j * result.ld] = self.data[j + j * self.ld];
                                }
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.ld] = self.data[i + j * self.ld];
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, self.cols)) : (j += 1) {
                                result.data[i * result.ld + j] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (i < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[i * result.ld + i] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[i * result.ld + i] = self.data[i * self.ld + i];
                                }
                            }

                            j = i + 1;
                            while (j < self.cols) : (j += 1) {
                                result.data[i * result.ld + j] = self.data[i * self.ld + j];
                            }
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, self.cols)) : (j += 1) {
                                result.data[i * result.ld + j] = self.data[i * self.ld + j];
                            }

                            if (i < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[i * result.ld + i] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[i * result.ld + i] = self.data[i * self.ld + i];
                                }
                            }

                            j = i + 1;
                            while (j < self.cols) : (j += 1) {
                                result.data[i * result.ld + j] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn toDenseArray(self: *const Triangular(T, uplo, diag, order), allocator: std.mem.Allocator, ctx: anytype) !Dense(T, order) {
            var result: Dense(T, order) = try .init(allocator, &.{ self.rows, self.cols });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, self.rows)) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                            }

                            if (j < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[j + j * result.strides[1]] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];
                                }
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    } else { // cl
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, self.rows)) : (i += 1) {
                                result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (j < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[j + j * result.strides[1]] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];
                                }
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, self.cols)) : (j += 1) {
                                result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (i < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[i * result.strides[0] + i] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];
                                }
                            }

                            j = i + 1;
                            while (j < self.cols) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                            }
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, self.cols)) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                            }

                            if (i < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[i * result.strides[0] + i] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];
                                }
                            }

                            j = i + 1;
                            while (j < self.cols) : (j += 1) {
                                result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn transpose(self: Triangular(T, uplo, diag, order)) Triangular(T, uplo.invert(), diag, order.invert()) {
            return .{
                .data = self.data,
                .rows = self.cols,
                .cols = self.rows,
                .ld = self.ld,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Triangular(T, uplo, diag, order),
            start: u32,
            row_end: u32,
            col_end: u32,
        ) !Triangular(T, uplo, diag, order) {
            if (start >= int.min(self.rows, self.cols) or
                row_end > self.rows or col_end > self.cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - start;
            const sub_cols = col_end - start;

            return .{
                .data = self.data + (start + start * self.ld),
                .rows = sub_rows,
                .cols = sub_cols,
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
    const R: type = EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y))));

    if (comptime !types.isTriangularMatrix(@TypeOf(x))) {
        var result: R = try .init(allocator, y.rows, y.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu
                    var j: u32 = 0;
                    while (j < y.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, y.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                            }
                        }

                        if (j < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                                }
                            }
                        }
                    }
                } else { // c cl
                    var j: u32 = 0;
                    while (j < y.cols) : (j += 1) {
                        if (j < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < y.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c ru
                    var j: u32 = 0;
                    while (j < y.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, y.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                            }
                        }

                        if (j < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                                }
                            }
                        }
                    }
                } else { // c rl
                    var j: u32 = 0;
                    while (j < y.cols) : (j += 1) {
                        if (j < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < y.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r cu
                    var i: u32 = 0;
                    while (i < y.rows) : (i += 1) {
                        if (i < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < y.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                            }
                        }
                    }
                } else { // r cl
                    var i: u32 = 0;
                    while (i < y.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, y.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                            }
                        }

                        if (i < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru
                    var i: u32 = 0;
                    while (i < y.rows) : (i += 1) {
                        if (i < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < y.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                            }
                        }
                    }
                } else { // r rl
                    var i: u32 = 0;
                    while (i < y.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, y.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                            }
                        }

                        if (i < int.min(y.rows, y.cols)) {
                            if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x, constants.one(Y, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x, constants.one(Y, ctx) catch unreachable, ctx);
                                }
                            } else {
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
        }

        return result;
    } else if (comptime !types.isTriangularMatrix(@TypeOf(y))) {
        var result: R = try .init(allocator, x.rows, x.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // c cu
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, x.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                            }
                        }

                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                } else { // c cl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // c ru
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, x.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                            }
                        }

                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                                }
                            }
                        }
                    }
                } else { // c rl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // r cu
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                            }
                        }
                    }
                } else { // r cl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, x.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                            }
                        }

                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // r ru
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                            }
                        }
                    }
                } else { // r rl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, x.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                            }
                        }

                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y, ctx);
                                }
                            } else {
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
        }

        return result;
    }

    if (x.rows != y.rows or x.cols != y.cols)
        return matrix.Error.DimensionMismatch;

    var result: R = try .init(allocator, x.rows, x.cols);
    errdefer result.deinit(allocator);

    // Two cases:
    // - Result is triangular
    // - Result is general
    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.isTriangularMatrix(@TypeOf(result))) { // uploOf(x) == uploOf(y)
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cu cu cu
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, x.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                            }
                        }

                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else { // cl cl cl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                                    }
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cu cu ru
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, x.rows)) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                            }
                        }

                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else { // cl cl rl
                    var j: u32 = 0;
                    while (j < x.cols) : (j += 1) {
                        if (j < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                                    }
                                }
                            }
                        }

                        var i: u32 = j + 1;
                        while (i < x.rows) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // ru ru cu
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                                    }
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                            }
                        }
                    }
                } else { // rl rl cl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, x.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                            }
                        }

                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) { // ru ru ru
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                                    }
                                }
                            }
                        }

                        var j: u32 = i + 1;
                        while (j < x.cols) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                            }
                        }
                    }
                } else { // rl rl rl
                    var i: u32 = 0;
                    while (i < x.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, x.cols)) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                            }
                        }

                        if (i < int.min(x.rows, x.cols)) {
                            if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                                    }
                                }
                            } else {
                                if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                    }
                                } else {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else { // uploOf(x) != uploOf(y)
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu cu
                        unreachable; // Handled above
                    } else { // c cu cl
                        var j: u32 = 0;
                        while (j < x.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, x.rows)) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }

                            if (j < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                                        }
                                    }
                                }
                            }

                            i = j + 1;
                            while (i < x.rows) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cl cu
                        var j: u32 = 0;
                        while (j < x.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, x.rows)) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                                }
                            }

                            if (j < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j + j * y.ld], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j + j * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j + j * y.ld], ctx);
                                        }
                                    }
                                }
                            }

                            i = j + 1;
                            while (i < x.rows) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }
                        }
                    } else { // c cl cl
                        unreachable; // Handled above
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cu ru
                        unreachable; // Handled above
                    } else { // c cu rl
                        var j: u32 = 0;
                        while (j < x.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, x.rows)) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }

                            if (j < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                                        }
                                    }
                                }
                            }

                            i = j + 1;
                            while (i < x.rows) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // c cl ru
                        var j: u32 = 0;
                        while (j < x.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, x.rows)) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                                }
                            }

                            if (j < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(constants.one(X, ctx) catch unreachable, y.data[j * y.ld + j], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y.data[j * y.ld + j]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y.data[j * y.ld + j], ctx);
                                        }
                                    }
                                }
                            }

                            i = j + 1;
                            while (i < x.rows) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }
                        }
                    } else { // c cl rl
                        unreachable; // Handled above
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru cu
                        unreachable; // Handled above
                    } else { // r ru cl
                        var i: u32 = 0;
                        while (i < x.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, x.cols)) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                                }
                            }

                            if (i < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                                        }
                                    }
                                }
                            }

                            j = i + 1;
                            while (j < x.cols) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r rl cu
                        var i: u32 = 0;
                        while (i < x.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, x.cols)) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }

                            if (i < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i + i * y.ld], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i + i * y.ld]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i + i * y.ld], ctx);
                                        }
                                    }
                                }
                            }

                            j = i + 1;
                            while (j < x.cols) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    } else { // r rl cl
                        unreachable; // Handled above
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r ru ru
                        unreachable; // Handled above
                    } else { // r ru rl
                        var i: u32 = 0;
                        while (i < x.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, x.cols)) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                                }
                            }

                            if (i < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                                        }
                                    }
                                }
                            }

                            j = i + 1;
                            while (j < x.cols) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // r rl ru
                        var i: u32 = 0;
                        while (i < x.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, x.cols)) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], constants.zero(Y, ctx) catch unreachable, ctx);
                                }
                            }

                            if (i < int.min(x.rows, x.cols)) {
                                if (comptime types.diagOf(@TypeOf(x)) == .unit) {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(constants.one(X, ctx) catch unreachable, y.data[i * y.ld + i], ctx);
                                        }
                                    }
                                } else {
                                    if (comptime types.diagOf(@TypeOf(y)) == .unit) {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], constants.one(Y, ctx) catch unreachable, ctx);
                                        }
                                    } else {
                                        if (comptime opinfo.@"fn".params.len == 2) {
                                            result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y.data[i * y.ld + i]);
                                        } else if (comptime opinfo.@"fn".params.len == 3) {
                                            result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y.data[i * y.ld + i], ctx);
                                        }
                                    }
                                }
                            }

                            j = i + 1;
                            while (j < x.cols) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(constants.zero(X, ctx) catch unreachable, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    } else { // r rl rl
                        unreachable; // Handled above
                    }
                }
            }
        }
    }

    return result;
}
