const std = @import("std");

// Maybe in place versions (_) also require allocator if the coerced type of the inputs is arbitrary precision

const types = @import("../types.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const EnsureMatrix = types.EnsureMatrix;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const matrix = @import("../array.zig");

const general = @import("general.zig");
const symmetric = @import("symmetric.zig");
const hermitian = @import("hermitian.zig");
const triangular = @import("triangular.zig");
const diagonal = @import("diagonal.zig");
const banded = @import("banded.zig");
const tridiagonal = @import("tridiagonal.zig");
const sparse = @import("sparse.zig");

pub fn apply1(
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    opts: struct {
        order: ?types.Order = null,
        // eventually will add axis, axes and keepdims opts, to apply2 as well
    },
    ctx: anytype,
) !EnsureMatrix(@TypeOf(x), ReturnType1(op, Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isMatrix(X))
        @compileError("apply1: x must be a matrix, got " ++ @typeName(X));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 1 and @typeInfo(@TypeOf(op)).@"fn".params.len != 2))
        @compileError("apply1: op must be a function of one argument, or a function of two arguments with the second argument being a context, got " ++ @typeName(@TypeOf(op)));

    switch (comptime types.matrixType(@TypeOf(x))) {
        .general => return general.apply1(allocator, x, op, .{ .order = opts.order }, ctx),
        .sparse => @compileError("apply1: sparse matrices are not implemented yet"),
        else => @compileError("apply1 is only defined for general and sparse matrices, got " ++ @typeName(@TypeOf(x))),
    }
}
