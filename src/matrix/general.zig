const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const Order = types.Order;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../array.zig");
const Dense = array.Dense;

pub fn General(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("General requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        strides: [2]u32,
        flags: Flags = .{},

        pub const empty: General(T) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .strides = .{ 0, 0 },
            .flags = .{ .order = .col_major, .owns_data = false },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            opts: struct {
                order: Order = .col_major,
            },
        ) !General(T) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return General(T){
                .data = (try allocator.alloc(T, rows * cols)).ptr,
                .rows = rows,
                .cols = cols,
                .strides = if (opts.order == .col_major) .{ 1, rows } else .{ cols, 1 },
                .flags = .{ .order = opts.order, .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            value: anytype,
            opts: struct {
                order: Order = .col_major,
            },
            ctx: anytype,
        ) !General(T) {
            const mat: General(T) = try .init(allocator, rows, cols, .{ .order = opts.order });
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
            opts: struct {
                order: Order = .col_major,
            },
            ctx: anytype,
        ) !General(T) {
            const mat: General(T) = try .init(allocator, size, size, .{ .order = opts.order });
            errdefer mat.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (opts.order == .col_major) {
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

        pub fn deinit(self: *General(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const General(T), row: u32, col: u32) !T {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            return self.data[row * self.strides[0] + col * self.strides[1]];
        }

        pub inline fn at(self: *const General(T), row: u32, col: u32) T {
            // Unchecked version of get. Assumes row and col are valid.
            return self.data[row * self.strides[0] + col * self.strides[1]];
        }

        pub fn set(self: *General(T), row: u32, col: u32, value: T) !void {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            self.data[row * self.strides[0] + col * self.strides[1]] = value;
        }

        pub inline fn put(self: *General(T), row: u32, col: u32, value: T) void {
            // Unchecked version of set. Assumes row and col are valid.
            self.data[row * self.strides[0] + col * self.strides[1]] = value;
        }

        pub fn toDenseArray(self: *const General(T), allocator: std.mem.Allocator, ctx: anytype) !Dense(T) {
            var result: Dense(T) = try .init(allocator, &.{ self.rows, self.cols }, .{ .order = self.flags.order });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                var i: u32 = 0;
                while (i < self.rows * self.cols) : (i += 1) {
                    result.data[i] = self.data[i];
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub fn transpose(self: General(T)) General(T) {
            return General(T){
                .data = self.data,
                .rows = self.cols,
                .cols = self.rows,
                .strides = .{ self.strides[1], self.strides[0] },
                .flags = .{
                    .order = self.flags.order.invert(),
                    .owns_data = false,
                },
            };
        }

        pub fn submatrix(
            self: *const General(T),
            row_start: u32,
            row_end: u32,
            col_start: u32,
            col_end: u32,
        ) !General(T) {
            if (row_start >= self.rows or col_start >= self.cols or
                row_end > self.rows or col_end > self.cols or
                row_start >= row_end or col_start >= col_end)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - row_start;
            const sub_cols = col_end - col_start;

            return General(T){
                .data = self.data + (row_start * self.strides[0] + col_start * self.strides[1]),
                .rows = sub_rows,
                .cols = sub_cols,
                .strides = self.strides,
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }

        pub fn asDenseArray(self: *const General(T)) Dense(T) {
            return Dense(T){
                .data = self.data,
                .ndim = 2,
                .shape = .{ self.rows, self.cols } ++ .{0} ** (array.max_dim - 2),
                .strides = .{ self.strides[0], self.strides[1] } ++ .{0} ** (array.max_dim - 2),
                .base = null,
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
    opts: struct {
        order: ?Order = null,
    },
    ctx: anytype,
) !General(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));

    if (comptime !types.isGeneralMatrix(@TypeOf(x))) {
        var result: General(ReturnType2(op, X, Y)) = try .init(allocator, y.rows, y.cols, .{ .order = opts.order orelse y.flags.order });
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (result.flags.order == .col_major) {
            if (y.flags.order == .col_major) {
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.strides[1]] = op(x, y.data[i + j * y.strides[1]]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.strides[1]] = try op(x, y.data[i + j * y.strides[1]], ctx);
                        }
                    }
                }
            } else {
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.strides[1]] = op(x, y.data[i * y.strides[0] + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.strides[1]] = try op(x, y.data[i * y.strides[0] + j], ctx);
                        }
                    }
                }
            }
        } else {
            if (y.flags.order == .col_major) {
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.strides[0] + j] = op(x, y.data[i + j * y.strides[1]]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.strides[0] + j] = try op(x, y.data[i + j * y.strides[1]], ctx);
                        }
                    }
                }
            } else {
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.strides[0] + j] = op(x, y.data[i * y.strides[0] + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.strides[0] + j] = try op(x, y.data[i * y.strides[0] + j], ctx);
                        }
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isGeneralMatrix(@TypeOf(y))) {
        var result: General(ReturnType2(op, X, Y)) = try .init(allocator, x.rows, x.cols, .{ .order = opts.order orelse x.flags.order });
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (result.flags.order == .col_major) {
            if (x.flags.order == .col_major) {
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y, ctx);
                        }
                    }
                }
            } else {
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y, ctx);
                        }
                    }
                }
            }
        } else {
            if (x.flags.order == .col_major) {
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y, ctx);
                        }
                    }
                }
            } else {
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.strides[0] + j] = op(x.data[i * x.strides[0] + j], y);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.strides[0] + j] = try op(x.data[i * x.strides[0] + j], y, ctx);
                        }
                    }
                }
            }
        }

        return result;
    }

    if (x.rows != y.rows or x.cols != y.cols)
        return matrix.Error.DimensionMismatch;

    var result: General(ReturnType2(op, X, Y)) = try .init(allocator, x.rows, x.cols, .{ .order = opts.order orelse x.flags.order.resolve2(y.flags.order) });
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (result.flags.order == .col_major) {
        if (x.flags.order == .col_major) {
            if (y.flags.order == .col_major) {
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx);
                        }
                    }
                }
            } else {
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.strides[1]] = op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.strides[1]] = try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx);
                        }
                    }
                }
            }
        } else {
            if (y.flags.order == .col_major) {
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx);
                        }
                    }
                }
            } else {
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i + j * result.strides[1]] = op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i + j * result.strides[1]] = try op(x.data[i * x.strides[0] + j], y.data[i * y.strides[0] + j], ctx);
                        }
                    }
                }
            }
        }
    } else {
        if (x.flags.order == .col_major) {
            if (y.flags.order == .col_major) {
                var j: u32 = 0;
                while (j < result.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < result.rows) : (i += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y.data[i + j * y.strides[1]], ctx);
                        }
                    }
                }
            } else {
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.strides[0] + j] = op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.strides[0] + j] = try op(x.data[i + j * x.strides[1]], y.data[i * y.strides[0] + j], ctx);
                        }
                    }
                }
            }
        } else {
            if (y.flags.order == .col_major) {
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
                        if (comptime opinfo.@"fn".params.len == 2) {
                            result.data[i * result.strides[0] + j] = op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]]);
                        } else if (comptime opinfo.@"fn".params.len == 3) {
                            result.data[i * result.strides[0] + j] = try op(x.data[i * x.strides[0] + j], y.data[i + j * y.strides[1]], ctx);
                        }
                    }
                }
            } else {
                var i: u32 = 0;
                while (i < result.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < result.cols) : (j += 1) {
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

    return result;
}
