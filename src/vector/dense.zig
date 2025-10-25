const std = @import("std");

const types = @import("../types.zig");
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const ops = @import("../ops.zig");
const linalg = @import("../linalg.zig");

const vector = @import("../vector.zig");
const Flags = vector.Flags;

pub fn Dense(T: type) type {
    if (!types.isNumeric(T))
        @compileError("vector.Dense requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        len: u32,
        inc: i32,
        flags: Flags = .{},

        pub const empty = Dense(T){
            .data = &.{},
            .len = 0,
            .inc = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(allocator: std.mem.Allocator, len: u32) !Dense(T) {
            if (len == 0)
                return vector.Error.ZeroLength;

            return .{
                .data = (try allocator.alloc(T, len)).ptr,
                .len = len,
                .inc = 1,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            len: u32,
            value: anytype,
            ctx: anytype,
        ) !Dense(T) {
            var vec: Dense(T) = try .init(allocator, len);
            errdefer vec.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                const value_casted: T = types.scast(T, value);

                var i: u32 = 0;
                while (i < len) : (i += 1) {
                    vec.data[i] = value_casted;
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return vec;
        }

        pub fn deinit(self: *Dense(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.len]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Dense(T), index: u32) !T {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            return if (self.inc > 0)
                self.data[index * types.scast(u32, self.inc)]
            else
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))];
        }

        pub inline fn at(self: *const Dense(T), index: u32) T {
            // Unchecked version of get. Assumes index is valid.
            return if (self.inc > 0)
                self.data[index * types.scast(u32, self.inc)]
            else
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))];
        }

        pub fn set(self: *Dense(T), index: u32, value: T) !void {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            if (self.inc > 0) {
                self.data[index * types.scast(u32, self.inc)] = value;
            } else {
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))] = value;
            }
        }

        pub inline fn put(self: *Dense(T), index: u32, value: T) void {
            // Unchecked version of set. Assumes index is valid.
            if (self.inc > 0) {
                self.data[index * types.scast(u32, self.inc)] = value;
            } else {
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))] = value;
            }
        }

        pub fn asDiagonal(self: *const Dense(T), rows: u32, cols: u32) !Dense(T) {
            if (rows == 0 or cols == 0 or
                (rows > self.len and cols > self.len))
                return vector.Error.ZeroDimension;

            if (self.inc != 1)
                return vector.Error.NonContiguousData;

            return .{
                .data = self.data,
                .rows = rows,
                .cols = cols,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn cleanup(self: *Dense(T), ctx: anytype) void {
            return _cleanup(self, self.len, ctx);
        }

        pub fn _cleanup(self: *Dense(T), num_elems: u32, ctx: anytype) void {
            if (comptime types.isArbitraryPrecision(T)) {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                if (self.inc == 1) {
                    var i: u32 = 0;
                    while (i < num_elems) : (i += 1) {
                        ops.deinit(self.data[i], ctx);
                    }
                } else {
                    var is: i32 = if (self.inc < 0) (-types.scast(i32, self.len) + 1) * self.inc else 0;
                    var i: u32 = 0;
                    while (i < num_elems) : (i += 1) {
                        ops.deinit(self.data[types.scast(u32, is)], ctx);

                        is += self.inc;
                    }
                }
            } else {
                comptime types.validateContext(@TypeOf(ctx), .{});
            }
        }
    };
}

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !Dense(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = ReturnType2(op, Numeric(X), Numeric(Y));

    if (comptime !types.isDenseVector(@TypeOf(x))) {
        var result: Dense(R) = try .init(allocator, y.len);
        errdefer result.deinit(allocator);

        var i: u32 = 0;

        errdefer result._cleanup(i, ctx);

        const opinfo = @typeInfo(@TypeOf(op));
        if (y.inc == 1) {
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x, y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x, y.data[i], ctx);
                }
            }
        } else {
            var iy: i32 = if (y.inc < 0) (-types.scast(i32, y.len) + 1) * y.inc else 0;
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x, y.data[types.scast(u32, iy)]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x, y.data[types.scast(u32, iy)], ctx);
                }

                iy += y.inc;
            }
        }

        return result;
    } else if (comptime !types.isDenseVector(@TypeOf(y))) {
        var result: Dense(R) = try .init(allocator, x.len);
        errdefer result.deinit(allocator);

        var i: u32 = 0;

        errdefer result._cleanup(i, ctx);

        const opinfo = @typeInfo(@TypeOf(op));
        if (x.inc == 1) {
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], y, ctx);
                }
            }
        } else {
            var ix: i32 = if (x.inc < 0) (-types.scast(i32, x.len) + 1) * x.inc else 0;
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[types.scast(u32, ix)], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[types.scast(u32, ix)], y, ctx);
                }

                ix += x.inc;
            }
        }

        return result;
    }

    if (x.len != y.len)
        return vector.Error.DimensionMismatch;

    var result: Dense(R) = try .init(allocator, x.len);
    errdefer result.deinit(allocator);

    var i: u32 = 0;

    errdefer result._cleanup(i, ctx);

    const opinfo = @typeInfo(@TypeOf(op));
    if (x.inc == 1 and y.inc == 1) {
        while (i < result.len) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x.data[i], y.data[i]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x.data[i], y.data[i], ctx);
            }
        }
    } else {
        var ix: i32 = if (x.inc < 0) (-types.scast(i32, x.len) + 1) * x.inc else 0;
        var iy: i32 = if (y.inc < 0) (-types.scast(i32, y.len) + 1) * y.inc else 0;
        while (i < result.len) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x.data[types.scast(u32, ix)], y.data[types.scast(u32, iy)]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x.data[types.scast(u32, ix)], y.data[types.scast(u32, iy)], ctx);
            }

            ix += x.inc;
            iy += y.inc;
        }
    }

    return result;
}
