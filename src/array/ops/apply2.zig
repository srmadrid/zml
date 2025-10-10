const std = @import("std");

const types = @import("../../types.zig");
const EnsureArray = types.EnsureArray;
const Coerce = types.Coerce;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;

const dense = @import("../dense.zig");
const strided = @import("../strided.zig");
// const sparse = @import("../sparse.zig");

///
pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y))
        @compileError("apply2: at least one of x or y must be an array, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op)).@"fn".params.len != 3))
        @compileError("apply2: op must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op)));

    if (comptime !types.isArray(X)) {
        switch (comptime types.arrayType(Y)) {
            .dense => return dense.apply2(allocator, x, y, op, ctx),
            .strided => return strided.apply2(allocator, x, y, op, ctx),
            .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else if (comptime !types.isArray(Y)) {
        switch (comptime types.arrayType(X)) {
            .dense => return dense.apply2(allocator, x, y, op, ctx),
            .strided => return strided.apply2(allocator, x, y, op, ctx),
            .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.arrayType(X)) {
            .dense => switch (comptime types.arrayType(Y)) {
                .dense => return dense.apply2(allocator, x, y, op, ctx),
                .strided => return strided.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(Y)) {
                .dense => return strided.apply2(allocator, x, y, op, ctx),
                .strided => return strided.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    }
}
