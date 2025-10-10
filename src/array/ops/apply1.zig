const std = @import("std");

const types = @import("../../types.zig");
const EnsureArray = types.EnsureArray;
const ReturnType1 = types.ReturnType1;
const Numeric = types.Numeric;

const dense = @import("../dense.zig");
const strided = @import("../strided.zig");

///
pub fn apply1(
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureArray(@TypeOf(x), ReturnType1(op, Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X))
        @compileError("apply1: x must be an array, got " ++ @typeName(X));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 1 and @typeInfo(@TypeOf(op)).@"fn".params.len != 2))
        @compileError("apply1: op must be a function of one argument, or a function of two arguments with the second argument being a context, got " ++ @typeName(@TypeOf(op)));

    switch (comptime types.arrayType(@TypeOf(x))) {
        .dense => return dense.apply1(allocator, x, op, ctx),
        .strided => return strided.apply1(allocator, x, op, ctx),
        .sparse => @compileError("apply1 not implemented for sparse arrays yet"),
        .numeric => unreachable,
    }
}
