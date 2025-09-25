const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const Uplo = types.Uplo;
const Diag = types.Diag;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../array.zig");
const Dense = array.Dense;

const linalg = @import("../linalg.zig");

pub fn General(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("General requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        ld: u32, // leading dimension
        flags: Flags = .{},

        pub const empty: General(T, order) = .{
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
        ) !General(T, order) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(T, rows * cols)).ptr,
                .rows = rows,
                .cols = cols,
                .ld = if (comptime order == .col_major) rows else cols,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            value: anytype,
            ctx: anytype,
        ) !General(T, order) {
            const mat: General(T, order) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                var i: u32 = 0;
                while (i < rows * cols) : (i += 1) {
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
        ) !General(T, order) {
            const mat: General(T, order) = try .init(allocator, size, size);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < size) : (j += 1) {
                        mat.data[j + j * size] = constants.one(T, ctx) catch unreachable;

                        var i: u32 = j + 1;
                        while (i < size) : (i += 1) {
                            mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                            mat.data[j + i * size] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < size) : (i += 1) {
                        mat.data[i * size + i] = constants.one(T, ctx) catch unreachable;

                        var j: u32 = i + 1;
                        while (j < size) : (j += 1) {
                            mat.data[i * size + j] = constants.zero(T, ctx) catch unreachable;
                            mat.data[j * size + i] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn deinit(self: *General(T, order), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const General(T, order), row: u32, col: u32) !T {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            return if (comptime order == .col_major)
                self.data[row + col * self.ld]
            else
                self.data[row * self.ld + col];
        }

        pub inline fn at(self: *const General(T, order), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid.
            return if (comptime order == .col_major)
                self.data[row + col * self.ld]
            else
                self.data[row * self.ld + col];
        }

        pub fn set(self: *General(T, order), row: u32, col: u32, value: T) !void {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime order == .col_major) {
                self.data[row + col * self.ld] = value;
            } else {
                self.data[row * self.ld + col] = value;
            }
        }

        pub inline fn put(self: *General(T, order), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid.
            if (comptime order == .col_major) {
                self.data[row + col * self.ld] = value;
            } else {
                self.data[row * self.ld + col] = value;
            }
        }

        pub fn copy(self: *const General(T, order), allocator: std.mem.Allocator, ctx: anytype) !General(T, order) {
            var mat: General(T, order) = try .init(allocator, self.rows, self.cols);
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.cols) : (j += 1) {
                        try linalg.blas.copy(
                            types.scast(i32, self.rows),
                            self.data + j * self.ld,
                            1,
                            mat.data + j * mat.ld,
                            1,
                            ctx,
                        );
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        try linalg.blas.copy(
                            types.scast(i32, self.cols),
                            self.data + i * self.ld,
                            1,
                            mat.data + i * mat.ld,
                            1,
                            ctx,
                        );
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn toDenseArray(self: *const General(T, order), allocator: std.mem.Allocator, ctx: anytype) !Dense(T, order) {
            var result: Dense(T, order) = try .init(allocator, &.{ self.rows, self.cols });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn transpose(self: General(T, order)) General(T, order.invert()) {
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
            self: *const General(T, order),
            row_start: u32,
            row_end: u32,
            col_start: u32,
            col_end: u32,
        ) !General(T, order) {
            if (row_start >= self.rows or col_start >= self.cols or
                row_end > self.rows or col_end > self.cols or
                row_start >= row_end or col_start >= col_end)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - row_start;
            const sub_cols = col_end - col_start;

            return .{
                .data = self.data + if (comptime order == .col_major)
                    row_start + col_start * self.ld
                else
                    row_start * self.ld + col_start,
                .rows = sub_rows,
                .cols = sub_cols,
                .ld = self.ld,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        // pub fn row(self: *const General(T, order), r: u32) !vector.Vector(T) {
        //     if (r >= self.rows)
        //         return matrix.Error.PositionOutOfBounds;

        //     return array.Vector(T){
        //         .data = if (comptime order == .col_major)
        //             self.data + row
        //         else
        //             self.data + row * self.ld,
        //         .len = self.cols,
        //         .stride = if (comptime order == .col_major) self.ld else 1,
        //     };
        // }

        pub fn asSymmetric(self: *const General(T, order), comptime uplo: Uplo) !matrix.Symmetric(T, uplo, order) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            return .{
                .data = self.data,
                .size = self.rows,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn asHermitian(self: *const General(T, order), comptime uplo: Uplo) !matrix.Hermitian(T, uplo, order) {
            comptime if (!types.isComplex(T))
                @compileError("Hermitian matrices require a complex type, got " ++ @typeName(T));

            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            return .{
                .data = self.data,
                .size = self.rows,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn asTriangular(self: *const General(T, order), comptime uplo: Uplo, comptime diag: Diag) !matrix.Triangular(T, uplo, diag, order) {
            return .{
                .data = self.data,
                .rows = self.rows,
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn asDenseArray(self: *const General(T, order)) Dense(T, order) {
            return Dense(T){
                .data = self.data,
                .ndim = 2,
                .shape = .{ self.rows, self.cols } ++ .{0} ** (array.max_dim - 2),
                .strides = if (comptime order == .col_major)
                    .{ 1, self.ld } ++ .{0} ** (array.max_dim - 2)
                else
                    .{ self.ld, 1 } ++ .{0} ** (array.max_dim - 2),
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
    ctx: anytype,
) !EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));
    const R: type = EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y));

    if (comptime !types.isGeneralMatrix(@TypeOf(x))) {
        var result: R = try .init(allocator, y.rows, y.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // c c
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                        }
                    }
                }
            } else { // c r
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // r c
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                        }
                    }
                }
            } else { // r r
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                        }
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isGeneralMatrix(@TypeOf(y))) {
        var result: R = try .init(allocator, x.rows, x.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) { // c c
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                        }
                    }
                }
            } else { // c r
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) { // r c
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                        }
                    }
                }
            } else { // r r
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
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

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // c c c
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                        }
                    }
                }
            } else { // c c r
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // c r c
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                        }
                    }
                }
            } else { // c r r
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                        }
                    }
                }
            }
        }
    } else {
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // r c c
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                        }
                    }
                }
            } else { // r c r
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // r r c
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                        }
                    }
                }
            } else { // r r r
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
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

    return result;
}
