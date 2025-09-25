const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const MulCoerce = types.MulCoerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const EnsureVector = types.EnsureVector;
const ReturnType2 = types.ReturnType2;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const matrix = @import("../array.zig");

const dense = @import("dense.zig");
const sparse = @import("sparse.zig");

const desp = @import("ops/desp.zig");
const spde = @import("ops/spde.zig");

const linalg = @import("../linalg.zig");

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureVector(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isVector(X) and !types.isVector(Y))
        @compileError("apply2: at least one of x or y must be a vector, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op)).@"fn".params.len != 3))
        @compileError("apply2: op must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op)));

    if (comptime !types.isVector(X)) {
        switch (comptime types.vectorType(Y)) {
            .dense => return dense.apply2(allocator, x, y, op, ctx),
            .sparse => return sparse.apply2(allocator, x, y, op, ctx),
            .numeric => unreachable,
        }
    } else if (comptime !types.isVector(Y)) {
        switch (comptime types.vectorType(X)) {
            .dense => return dense.apply2(allocator, x, y, op, ctx),
            .sparse => return sparse.apply2(allocator, x, y, op, ctx),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.vectorType(X)) {
            .dense => switch (comptime types.vectorType(Y)) {
                .dense => return dense.apply2(allocator, x, y, op, ctx),
                .sparse => return desp.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .sparse => switch (comptime types.vectorType(Y)) {
                .dense => return spde.apply2(allocator, x, y, op, ctx),
                .sparse => return sparse.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .numeric => unreachable,
        }
    }
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
