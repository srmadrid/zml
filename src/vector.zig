const std = @import("std");

const types = @import("types.zig");
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const Coerce = types.Coerce;
const MulCoerce = types.MulCoerce;
const int = @import("int.zig");
const ops = @import("ops.zig");
const linalg = @import("linalg.zig");

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
            .inc = 0,
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
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))];
        }

        pub inline fn at(self: *const Vector(T), index: u32) T {
            // Unchecked version of get. Assumes index is valid.
            return if (self.inc > 0)
                self.data[index * types.scast(u32, self.inc)]
            else
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))];
        }

        pub fn set(self: *Vector(T), index: u32, value: T) !void {
            if (index >= self.len)
                return Error.PositionOutOfBounds;

            if (self.inc > 0) {
                self.data[index * types.scast(u32, self.inc)] = value;
            } else {
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))] = value;
            }
        }

        pub inline fn put(self: *Vector(T), index: u32, value: T) void {
            // Unchecked version of set. Assumes index is valid.
            if (self.inc > 0) {
                self.data[index * types.scast(u32, self.inc)] = value;
            } else {
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))] = value;
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
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Vector(ReturnType2(op, Numeric(X), Numeric(Y)));

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

pub inline fn add(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (!types.isVector(@TypeOf(x)) or !types.isVector(@TypeOf(y)))
        @compileError("Both arguments to add must be vector types");

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            types.validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2(
        allocator,
        x,
        y,
        ops.add,
        ctx,
    );
}

pub inline fn sub(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (!types.isVector(@TypeOf(x)) or !types.isVector(@TypeOf(y)))
        @compileError("Both arguments to sub must be vector types");

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            types.validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2(
        allocator,
        x,
        y,
        ops.sub,
        ctx,
    );
}

pub inline fn mul(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !MulCoerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(Numeric(X), Numeric(Y));

    comptime if (!types.isVector(X) and !types.isVector(Y))
        @compileError("At least one of the arguments must be a vector type");

    if (comptime types.isVector(X) and types.isVector(Y)) { // vector * vector
        comptime if (types.isArbitraryPrecision(C)) {
            @compileError("Arbitrary precision types not implemented yet");
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        };

        return linalg.dot(x, y, ctx);
    } else {
        comptime if (types.isArbitraryPrecision(C)) { // scalar * vector  or  vector * scalar
            @compileError("Arbitrary precision types not implemented yet");
        } else {
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            }
        };

        return apply2(
            allocator,
            x,
            y,
            ops.mul,
            ctx,
        );
    }
}

pub inline fn div(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(Numeric(X), Numeric(Y));

    comptime if (!types.isVector(X) and types.isVector(Y))
        @compileError("First argument must be a vector type and second argument must be a scalar type");

    comptime if (types.isArbitraryPrecision(C)) {
        @compileError("Arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.div,
        ctx,
    );
}

pub const Error = error{
    ZeroLength,
    PositionOutOfBounds,
    DimensionMismatch,
    NonContiguousData,
};
