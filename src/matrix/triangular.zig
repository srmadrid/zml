//! Storage scheme:
//!
//! Full `n`-by-`n` storage only accessing the upper or lower triangular part of
//! the matrix, with an option for unit triangular matrices.

const std = @import("std");

const types = @import("../types.zig");
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

pub fn Triangular(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("Triangular requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        strides: [2]u32,
        uplo: Uplo,
        diag: Diag,
        flags: Flags = .{},

        pub const empty: Triangular(T) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .strides = .{ 0, 0 },
            .uplo = .upper,
            .diag = .non_unit,
            .flags = .{ .order = .col_major, .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            opts: struct {
                uplo: Uplo = .upper,
                diag: Diag = .non_unit,
                order: Order = .col_major,
            },
        ) !Triangular(T) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return Triangular(T){
                .data = (try allocator.alloc(T, rows * cols)).ptr,
                .rows = rows,
                .cols = cols,
                .strides = if (opts.order == .col_major) .{ 1, rows } else .{ cols, 1 },
                .uplo = opts.uplo,
                .diag = opts.diag,
                .flags = .{ .order = opts.order, .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            value: anytype,
            opts: struct {
                uplo: Uplo = .upper,
                diag: Diag = .non_unit,
                order: Order = .col_major,
            },
            ctx: anytype,
        ) !Triangular(T) {
            const mat: Triangular(T) = try .init(allocator, rows, cols, .{ .uplo = opts.uplo, .diag = opts.diag, .order = opts.order });
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                if (opts.order == .col_major) {
                    if (opts.uplo == .upper) {
                        if (opts.diag == .unit) {
                            var j: u32 = 0;
                            while (j < cols) : (j += 1) {
                                var i: u32 = 0;
                                while (i < j) : (i += 1) {
                                    mat.data[i + j * rows] = value_casted;
                                }
                            }
                        } else {
                            var j: u32 = 0;
                            while (j < cols) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    mat.data[i + j * rows] = value_casted;
                                }
                            }
                        }
                    } else {
                        if (opts.diag == .unit) {
                            var j: u32 = 0;
                            while (j < cols) : (j += 1) {
                                var i: u32 = j + 1;
                                while (i < rows) : (i += 1) {
                                    mat.data[i + j * rows] = value_casted;
                                }
                            }
                        } else {
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
                    if (opts.uplo == .upper) {
                        if (opts.diag == .unit) {
                            var i: u32 = 0;
                            while (i < rows) : (i += 1) {
                                var j: u32 = i + 1;
                                while (j < cols) : (j += 1) {
                                    mat.data[i * cols + j] = value_casted;
                                }
                            }
                        } else {
                            var i: u32 = 0;
                            while (i < rows) : (i += 1) {
                                var j: u32 = i;
                                while (j < cols) : (j += 1) {
                                    mat.data[i * cols + j] = value_casted;
                                }
                            }
                        }
                    } else {
                        if (opts.diag == .unit) {
                            var i: u32 = 0;
                            while (i < rows) : (i += 1) {
                                var j: u32 = 0;
                                while (j < i) : (j += 1) {
                                    mat.data[i * cols + j] = value_casted;
                                }
                            }
                        } else {
                            var i: u32 = 0;
                            while (i < rows) : (i += 1) {
                                var j: u32 = 0;
                                while (j <= i) : (j += 1) {
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
            opts: struct {
                uplo: Uplo = .upper,
                diag: Diag = .non_unit,
                order: Order = .col_major,
            },
            ctx: anytype,
        ) !Triangular(T) {
            const mat: Triangular(T) = try .init(allocator, size, size, .{ .uplo = opts.uplo, .diag = opts.diag, .order = opts.order });
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (opts.order == .col_major) {
                    if (opts.uplo == .upper) {
                        if (opts.diag == .unit) {
                            var j: u32 = 0;
                            while (j < size) : (j += 1) {
                                var i: u32 = 0;
                                while (i < j) : (i += 1) {
                                    mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                                }
                            }
                        } else {
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
                        if (opts.diag == .unit) {
                            var j: u32 = 0;
                            while (j < size) : (j += 1) {
                                var i: u32 = j + 1;
                                while (i < size) : (i += 1) {
                                    mat.data[i + j * size] = constants.zero(T, ctx) catch unreachable;
                                }
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
                    }
                } else {
                    if (opts.uplo == .upper) {
                        if (opts.diag == .unit) {
                            var i: u32 = 0;
                            while (i < size) : (i += 1) {
                                var j: u32 = i + 1;
                                while (j < size) : (j += 1) {
                                    mat.data[i * size + j] = constants.zero(T, ctx) catch unreachable;
                                }
                            }
                        } else {
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
                        if (opts.diag == .unit) {
                            var i: u32 = 0;
                            while (i < size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < i) : (j += 1) {
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
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn deinit(self: *Triangular(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Triangular(T), row: u32, col: u32) !T {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (self.uplo == .upper) {
                if (row > col)
                    return constants.zero(T, .{}) catch unreachable;
            } else {
                if (row < col)
                    return constants.zero(T, .{}) catch unreachable;
            }

            if (self.diag == .unit) {
                if (row == col)
                    return constants.one(T, .{}) catch unreachable;
            }

            return self.data[row * self.strides[0] + col * self.strides[1]];
        }

        pub inline fn at(self: *const Triangular(T), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid and on
            // the correct triangular part, and outside the diagonal if diag
            // triangular.
            return self.data[row * self.strides[0] + col * self.strides[1]];
        }

        pub fn set(self: *Triangular(T), row: u32, col: u32, value: T) !void {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (self.uplo == .upper) {
                if (row > col)
                    return matrix.Error.PositionOutOfBounds;
            } else {
                if (row < col)
                    return matrix.Error.PositionOutOfBounds;
            }

            if (self.diag == .unit) {
                if (row == col)
                    return matrix.Error.PositionOutOfBounds;
            }

            self.data[row * self.strides[0] + col * self.strides[1]] = value;
        }

        pub inline fn put(self: *Triangular(T), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid and on
            // the correct triangular part, and outside the diagonal if diag
            // triangular.
            self.data[row * self.strides[0] + col * self.strides[1]] = value;
        }

        pub fn toGeneral(self: Triangular(T), allocator: std.mem.Allocator, ctx: anytype) !General(T) {
            var result: General(T) = try .init(allocator, self.rows, self.cols, .{ .order = self.flags.order });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (self.flags.order == .col_major) {
                    if (self.uplo == .upper) {
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.strides[1]];
                            }

                            if (self.diag == .unit) {
                                result.data[j + j * result.strides[1]] = constants.one(T, ctx) catch unreachable;
                            } else {
                                result.data[j + j * result.strides[1]] = self.data[j + j * self.strides[1]];
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (self.diag == .unit) {
                                result.data[j + j * result.strides[1]] = constants.one(T, ctx) catch unreachable;
                            } else {
                                result.data[j + j * result.strides[1]] = self.data[j + j * self.strides[1]];
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.strides[1]];
                            }
                        }
                    }
                } else {
                    if (self.uplo == .upper) {
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (self.diag == .unit) {
                                result.data[i * result.strides[0] + i] = constants.one(T, ctx) catch unreachable;
                            } else {
                                result.data[i * result.strides[0] + i] = self.data[i * self.strides[0] + i];
                            }

                            j = i + 1;
                            while (j < self.cols) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + j];
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + j];
                            }

                            if (self.diag == .unit) {
                                result.data[i * result.strides[0] + i] = constants.one(T, ctx) catch unreachable;
                            } else {
                                result.data[i * result.strides[0] + i] = self.data[i * self.strides[0] + i];
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

        pub fn toDense(self: *const Triangular(T), allocator: std.mem.Allocator, ctx: anytype) !Dense(T) {
            var result: Dense(T) = try .init(allocator, &.{ self.rows, self.cols }, .{ .order = self.flags.order });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (self.flags.order == .col_major) {
                    if (self.uplo == .upper) {
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.strides[1]];
                            }

                            if (self.diag == .unit) {
                                result.data[j + j * result.strides[1]] = constants.one(T, ctx) catch unreachable;
                            } else {
                                result.data[j + j * result.strides[1]] = self.data[j + j * self.strides[1]];
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (self.diag == .unit) {
                                result.data[j + j * result.strides[1]] = constants.one(T, ctx) catch unreachable;
                            } else {
                                result.data[j + j * result.strides[1]] = self.data[j + j * self.strides[1]];
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.strides[1]];
                            }
                        }
                    }
                } else {
                    if (self.uplo == .upper) {
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (self.diag == .unit) {
                                result.data[i * result.strides[0] + i] = constants.one(T, ctx) catch unreachable;
                            } else {
                                result.data[i * result.strides[0] + i] = self.data[i * self.strides[0] + i];
                            }

                            j = i + 1;
                            while (j < self.cols) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + j];
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + j];
                            }

                            if (self.diag == .unit) {
                                result.data[i * result.strides[0] + i] = constants.one(T, ctx) catch unreachable;
                            } else {
                                result.data[i * result.strides[0] + i] = self.data[i * self.strides[0] + i];
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

        pub fn transpose(self: Triangular(T)) Triangular(T) {
            return Triangular(T){
                .data = self.data,
                .rows = self.cols,
                .cols = self.rows,
                .strides = .{ self.strides[1], self.strides[0] },
                .uplo = self.uplo.invert(),
                .diag = self.diag,
                .flags = .{
                    .order = self.flags.order.invert(),
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Triangular(T),
            start: u32,
            row_end: u32,
            col_end: u32,
        ) !Triangular(T) {
            if (start >= int.min(self.rows, self.cols) or
                row_end > self.rows or col_end > self.cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - start;
            const sub_cols = col_end - start;

            return Triangular(T){
                .data = self.data + (start * self.strides[0] + start * self.strides[1]),
                .rows = sub_rows,
                .cols = sub_cols,
                .strides = self.strides,
                .uplo = self.uplo,
                .diag = self.diag,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }
    };
}
