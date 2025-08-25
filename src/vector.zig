const std = @import("std");

const types = @import("types.zig");
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const int = @import("int.zig");

const matrix = @import("matrix.zig");
const Diagonal = matrix.Diagonal;

pub const Flags = packed struct {
    owns_data: bool = true,
};

pub fn Vector(T: type) type {
    if (!types.isNumeric(T))
        @compileError("Vector requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        len: u32,
        inc: i32,
        flags: Flags = .{},

        pub const empty = Vector(T){
            .data = &.{},
            .len = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(allocator: std.mem.Allocator, len: u32) !Vector(T) {
            if (len == 0)
                return Error.ZeroLength;

            return .{
                .data = (try allocator.alloc(T, len)).ptr,
                .len = len,
                .inc = 1,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn deinit(self: *Vector(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.len]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Vector(T), index: u32) !T {
            if (index >= self.len)
                return Error.PositionOutOfBounds;

            return if (self.inc > 0)
                self.data[index * types.scast(u32, self.inc)]
            else
                self.data[(-self.len + 1) * types.scast(u32, int.abs(self.inc)) - index * types.scast(u32, int.abs(self.inc))];
        }

        pub inline fn at(self: *const Vector(T), index: u32) T {
            // Unchecked version of get. Assumes index is valid.
            return if (self.inc > 0)
                self.data[index * types.scast(u32, self.inc)]
            else
                self.data[(-self.len + 1) * types.scast(u32, int.abs(self.inc)) - index * types.scast(u32, int.abs(self.inc))];
        }

        pub fn set(self: *Vector(T), index: u32, value: T) !void {
            if (index >= self.len)
                return Error.PositionOutOfBounds;

            if (self.inc > 0) {
                self.data[index * types.scast(u32, self.inc)] = value;
            } else {
                self.data[(-self.len + 1) * types.scast(u32, int.abs(self.inc)) - index * types.scast(u32, int.abs(self.inc))] = value;
            }
        }

        pub inline fn put(self: *Vector(T), index: u32, value: T) void {
            // Unchecked version of set. Assumes index is valid.
            if (self.inc > 0) {
                self.data[index * types.scast(u32, self.inc)] = value;
            } else {
                self.data[(-self.len + 1) * types.scast(u32, int.abs(self.inc)) - index * types.scast(u32, int.abs(self.inc))] = value;
            }
        }

        pub fn asDiagonal(self: *const Vector(T), rows: u32, cols: u32) !Diagonal(T) {
            if (rows == 0 or cols == 0 or
                (rows > self.len and cols > self.len))
                return Error.ZeroDimension;

            if (self.inc != 1)
                return Error.NonContiguousData;

            return .{
                .data = self.data,
                .rows = rows,
                .cols = cols,
                .flags = .{ .owns_data = false },
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
) !Vector(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));
    const R: type = Vector(ReturnType2(op, X, Y));

    comptime if (!types.isVector(X) and !types.isVector(Y))
        @compileError("apply2: at least one of x or y must be a vector, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op)).@"fn".params.len != 3))
        @compileError("apply2: op must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op)));

    if (comptime !types.isVector(@TypeOf(x))) {
        var result: R = try .init(allocator, y.len);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (y.inc == 1) {
            var i: u32 = 0;
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x, y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x, y.data[i], ctx);
                }
            }
        } else {
            var iy: i32 = if (y.inc < 0) (-types.scast(i32, y.len) + 1) * y.inc else 0;
            var i: u32 = 0;
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
    } else if (comptime !types.isVector(@TypeOf(y))) {
        var result: R = try .init(allocator, x.len);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (x.inc == 1) {
            var i: u32 = 0;
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], y, ctx);
                }
            }
        } else {
            var ix: i32 = if (x.inc < 0) (-types.scast(i32, x.len) + 1) * x.inc else 0;
            var i: u32 = 0;
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
        return Error.DimensionMismatch;

    var result: R = try .init(allocator, x.len);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (x.inc == 1 and y.inc == 1) {
        var i: u32 = 0;
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
        var i: u32 = 0;
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

pub const Error = error{
    ZeroLength,
    PositionOutOfBounds,
    DimensionMismatch,
    NonContiguousData,
};
