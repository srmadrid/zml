const std = @import("std");

const types = @import("../../types.zig");
const EnsureVector = types.EnsureVector;
const Coerce = types.Coerce;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;

const dense = @import("../dense.zig");
const sparse = @import("../sparse.zig");

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
                .sparse => return @import("apply2/desp.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .sparse => switch (comptime types.vectorType(Y)) {
                .dense => return @import("apply2/spde.zig").apply2(allocator, x, y, op, ctx),
                .sparse => return sparse.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .numeric => unreachable,
        }
    }
}
