//! Storage scheme:
//!
//! A banded matrix with `lower` subdiagonals and `upper` superdiagonals is
//! stored in a rectangular array of size `n × (lower + upper + 1)`, where `n`
//! is the matrix dimension. The storage layout depends on the matrix order
//! (row-major vs column-major):
//! - For row-major: Each row of storage corresponds to a matrix row, each
//! column to a diagonal band.
//! - For col-major: Each column of storage corresponds to a matrix column, each
//! row to a diagonal band.
//!
//! Example 1: 5×5 matrix with lower=2, upper=1 (bandwidth = 4), order=row-major
//!
//! ```zig
//! Original matrix:        Banded storage:
//! [d₀ u₀  0  0  0]       [*  *  d₀ u₀]  ← row 0
//! [l₀ d₁ u₁  0  0]       [*  l₀ d₁ u₁]  ← row 1
//! [l₁ l₂ d₂ u₂  0]   →   [l₁ l₂ d₂ u₂]  ← row 2
//! [ 0 l₃ l₄ d₃ u₃]       [l₃ l₄ d₃ u₃]  ← row 3
//! [ 0  0 l₅ l₆ d₄]       [l₅ l₆ d₄  *]  ← row 4
//! ```
//!
//! Storage mapping: matrix element (i,j) → storage position (i, lower + j - i)
//!
//! Example 2: 5×5 matrix with lower=2, upper=1 (bandwidth = 4), order=col-major
//!
//! ```zig
//! Original matrix:        Banded storage:
//! [d₀ u₀  0  0  0]       [*    l₀   l₁    0    0]
//! [l₀ d₁ u₁  0  0]       [*    d₁   l₂   l₃    0]
//! [l₁ l₂ d₂ u₂  0]   →   [d₀   u₁   d₂   l₄   l₅]
//! [ 0 l₃ l₄ d₃ u₃]       [u₀   u₂   u₃   d₃   l₆]
//! [ 0  0 l₅ l₆ d₄]       [ 0    0    0   u₃   d₄]
//!                          ↑    ↑    ↑    ↑    ↑
//!                        col0 col1 col2 col3 col4
//! ```
//!
//! Storage mapping: matrix element (i,j) → storage position (upper + i - j, j)
//!
//! Where `*` represents unused storage locations. The storage uses
//! `n × (lower + upper + 1)` = `5 × 4` = 20 elements total.
//!
//! row-major column mapping:
//! - Column 0: 2nd subdiagonal   (j - i = -2)
//! - Column 1: 1st subdiagonal   (j - i = -1)
//! - Column 2: Main diagonal     (j - i =  0)
//! - Column 3: 1st superdiagonal (j - i = +1)
//!
//! col-major row mapping:
//! - Row 0: 1st superdiagonal (i - j = -1)
//! - Row 1: Main diagonal     (i - j =  0)
//! - Row 2: 1st subdiagonal   (i - j = +1)
//! - Row 3: 2nd subdiagonal   (i - j = +2)

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

pub fn Banded(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("Banded requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        strides: [2]u32,
        lower: u32,
        upper: u32,
        flags: Flags = .{},

        pub const empty: Banded(T) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .strides = .{ 0, 0 },
            .lower = 0,
            .upper = 0,
            .flags = .{ .order = .col_major, .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            lower: u32,
            upper: u32,
            opts: struct {
                order: Order = .col_major,
            },
        ) !Banded(T) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return Banded(T){
                .data = (try allocator.alloc(T, (lower + upper + 1) * (if (opts.order == .col_major) cols else rows))).ptr,
                .rows = rows,
                .cols = cols,
                .strides = if (opts.order == .col_major) .{ 1, lower + upper + 1 } else .{ lower + upper + 1, 1 },
                .lower = lower,
                .upper = upper,
                .flags = .{ .order = opts.order, .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            lower: u32,
            upper: u32,
            value: anytype,
            opts: struct {
                order: Order = .col_major,
            },
            ctx: anytype,
        ) !Banded(T) {
            const mat: Banded(T) = try .init(allocator, rows, cols, lower, upper, .{ .order = opts.order });
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                if (opts.order == .col_major) {
                    var j: u32 = 0;
                    while (j < cols) : (j += 1) {
                        var i: u32 = if (j < upper) 0 else j - upper;
                        while (i <= int.min(rows - 1, j + lower)) : (i += 1) {
                            mat.data[(upper + i - j) + j * (lower + upper + 1)] = value_casted;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < rows) : (i += 1) {
                        var j: u32 = if (i < lower) 0 else i - lower;
                        while (j <= int.min(cols - 1, i + upper)) : (j += 1) {
                            mat.data[i * (lower + upper + 1) + (lower + j - i)] = value_casted;
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
            lower: u32,
            upper: u32,
            opts: struct {
                order: Order = .col_major,
            },
            ctx: anytype,
        ) !Banded(T) {
            const mat: Banded(T) = try .init(allocator, size, size, lower, upper, .{ .order = opts.order });
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (opts.order == .col_major) {
                    var j: u32 = 0;
                    while (j < size) : (j += 1) {
                        var i: u32 = if (j < upper) 0 else j - upper;
                        while (i <= int.min(size - 1, j + lower)) : (i += 1) {
                            if (i == j) {
                                mat.data[(upper + i - j) + j * (lower + upper + 1)] = constants.one(T, ctx) catch unreachable;
                            } else {
                                mat.data[(upper + i - j) + j * (lower + upper + 1)] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < size) : (i += 1) {
                        var j: u32 = if (i < lower) 0 else i - lower;
                        while (j <= int.min(size - 1, i + upper)) : (j += 1) {
                            if (i == j) {
                                mat.data[i * (lower + upper + 1) + (lower + j - i)] = constants.one(T, ctx) catch unreachable;
                            } else {
                                mat.data[i * (lower + upper + 1) + (lower + j - i)] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return mat;
        }

        pub fn deinit(self: *Banded(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. (self.lower + self.upper + 1) * (if (self.flags.order == .col_major) self.cols else self.rows)]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Banded(T), row: u32, col: u32) !T {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (col + self.lower < row or col > row + self.upper) {
                return constants.zero(T, .{}) catch unreachable;
            } else {
                if (self.flags.order == .col_major) {
                    return self.data[(self.upper + row - col) * self.strides[0] + col * self.strides[1]];
                } else {
                    return self.data[row * self.strides[0] + (self.lower + col - row) * self.strides[1]];
                }
            }
        }

        pub inline fn at(self: *const Banded(T), row: u32, col: u32) T {
            // Unchecked version of get. Assumes position is valid, i.e., within
            // matrix bounds, and in banded range.
            if (self.flags.order == .col_major) {
                return self.data[(self.upper + row - col) * self.strides[0] + col * self.strides[1]];
            } else {
                return self.data[row * self.strides[0] + (self.lower + col - row) * self.strides[1]];
            }
        }

        pub fn set(self: *Banded(T), row: u32, col: u32, value: T) !void {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (col + self.lower < row or col > row + self.upper) {
                return matrix.Error.PositionOutOfBounds;
            } else {
                if (self.flags.order == .col_major) {
                    self.data[(self.upper + row - col) * self.strides[0] + col * self.strides[1]] = value;
                } else {
                    self.data[row * self.strides[0] + (self.lower + col - row) * self.strides[1]] = value;
                }
            }
        }

        pub inline fn put(self: *Banded(T), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes position is valid, i.e., within
            // matrix bounds, and in banded range.
            if (self.flags.order == .col_major) {
                self.data[(self.upper + row - col) * self.strides[0] + col * self.strides[1]] = value;
            } else {
                self.data[row * self.strides[0] + (self.lower + col - row) * self.strides[1]] = value;
            }
        }

        pub fn toGeneral(self: Banded(T), allocator: std.mem.Allocator, ctx: anytype) !General(T) {
            var result: General(T) = try .init(allocator, self.rows, self.cols, .{ .order = self.flags.order });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (self.flags.order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < if (j < self.upper) 0 else j - self.upper) : (i += 1) {
                            result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                        }

                        while (i <= int.min(self.rows - 1, j + self.lower)) : (i += 1) {
                            result.data[i + j * result.strides[1]] = self.data[(self.upper + i - j) + j * self.strides[1]];
                        }

                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < if (i < self.lower) 0 else i - self.lower) : (j += 1) {
                            result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                        }

                        while (j <= int.min(self.cols - 1, i + self.upper)) : (j += 1) {
                            result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + (self.lower + j - i)];
                        }

                        while (j < self.cols) : (j += 1) {
                            result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn toDense(self: *const Banded(T), allocator: std.mem.Allocator, ctx: anytype) !Dense(T) {
            var result: Dense(T) = try .init(allocator, &.{ self.rows, self.cols }, .{ .order = self.flags.order });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (self.flags.order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < if (j < self.upper) 0 else j - self.upper) : (i += 1) {
                            result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                        }

                        while (i <= int.min(self.rows - 1, j + self.lower)) : (i += 1) {
                            result.data[i + j * result.strides[1]] = self.data[(self.upper + i - j) + j * self.strides[1]];
                        }

                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < if (i < self.lower) 0 else i - self.lower) : (j += 1) {
                            result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                        }

                        while (j <= int.min(self.cols - 1, i + self.upper)) : (j += 1) {
                            result.data[i * result.strides[0] + j] = self.data[i * self.strides[0] + (self.lower + j - i)];
                        }

                        while (j < self.cols) : (j += 1) {
                            result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn transpose(self: Banded(T)) Banded(T) {
            return Banded(T){
                .data = self.data,
                .rows = self.cols,
                .cols = self.rows,
                .strides = .{ self.strides[1], self.strides[0] },
                .lower = self.upper,
                .upper = self.lower,
                .flags = .{
                    .order = self.flags.order.invert(),
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const Banded(T),
            start: u32,
            row_end: u32,
            col_end: u32,
        ) !Banded(T) {
            std.debug.print("Banded submatrix does not work!\n", .{});
            if (start >= int.min(self.rows, self.cols) or
                row_end > self.rows or col_end > self.cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - start;
            const sub_cols = col_end - start;

            return Banded(T){
                .data = self.data + start * self.strides[0] + start * self.strides[1] -
                    (if (self.flags.order == .col_major) self.upper else self.lower),
                .rows = sub_rows,
                .cols = sub_cols,
                .strides = self.strides,
                .lower = int.min(self.lower, sub_rows - 1),
                .upper = int.min(self.upper, sub_cols - 1),
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }
    };
}
