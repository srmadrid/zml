const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const MulCoerce = types.MulCoerce;
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
const permutation = @import("permutation.zig");
const sparse = @import("sparse.zig");

const gesy = @import("ops/gesy.zig");
const gehe = @import("ops/gehe.zig");
const getr = @import("ops/getr.zig");
const gedi = @import("ops/gedi.zig");
const geba = @import("ops/geba.zig");
const gegt = @import("ops/gegt.zig");
const gepe = @import("ops/gepe.zig");
const syge = @import("ops/syge.zig");
const syhe = @import("ops/syhe.zig");
const sytr = @import("ops/sytr.zig");
const sydi = @import("ops/sydi.zig");
const syba = @import("ops/syba.zig");
const sygt = @import("ops/sygt.zig");
const sype = @import("ops/sype.zig");
const hege = @import("ops/hege.zig");
const hesy = @import("ops/hesy.zig");
const hetr = @import("ops/hetr.zig");
const hedi = @import("ops/hedi.zig");
const heba = @import("ops/heba.zig");
const hegt = @import("ops/hegt.zig");
const hepe = @import("ops/hepe.zig");
const trge = @import("ops/trge.zig");
const trsy = @import("ops/trsy.zig");
const trhe = @import("ops/trhe.zig");
const trdi = @import("ops/trdi.zig");
const trba = @import("ops/trba.zig");
const trgt = @import("ops/trgt.zig");
const trpe = @import("ops/trpe.zig");
const dige = @import("ops/dige.zig");
const disy = @import("ops/disy.zig");
const dihe = @import("ops/dihe.zig");
const ditr = @import("ops/ditr.zig");
const diba = @import("ops/diba.zig");
const digt = @import("ops/digt.zig");
const dipe = @import("ops/dipe.zig");
const bage = @import("ops/bage.zig");
const basy = @import("ops/basy.zig");
const bahe = @import("ops/bahe.zig");
const batr = @import("ops/batr.zig");
const badi = @import("ops/badi.zig");
const bagt = @import("ops/bagt.zig");
const bape = @import("ops/bape.zig");
const gtge = @import("ops/gtge.zig");
const gtsy = @import("ops/gtsy.zig");
const gthe = @import("ops/gthe.zig");
const gttr = @import("ops/gttr.zig");
const gtdi = @import("ops/gtdi.zig");
const gtba = @import("ops/gtba.zig");
const gtpe = @import("ops/gtpe.zig");
const pege = @import("ops/pege.zig");
const pesy = @import("ops/pesy.zig");
const pehe = @import("ops/pehe.zig");
const petr = @import("ops/petr.zig");
const pedi = @import("ops/pedi.zig");
const peba = @import("ops/peba.zig");
const pegt = @import("ops/pegt.zig");

const linalg = @import("../linalg.zig");

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isMatrix(X) and !types.isMatrix(Y))
        @compileError("apply2: at least one of x or y must be a matrix, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op)).@"fn".params.len != 3))
        @compileError("apply2: op must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op)));

    if (comptime !types.isMatrix(X)) {
        switch (comptime types.matrixType(Y)) {
            .general => return general.apply2(allocator, x, y, op, ctx),
            .symmetric => return symmetric.apply2(allocator, x, y, op, ctx),
            .hermitian => return hermitian.apply2(allocator, x, y, op, ctx),
            .triangular => return triangular.apply2(allocator, x, y, op, ctx),
            .diagonal => return diagonal.apply2(allocator, x, y, op, ctx),
            .banded => return banded.apply2(allocator, x, y, op, ctx),
            .tridiagonal => return tridiagonal.apply2(allocator, x, y, op, ctx),
            .permutation => return permutation.apply2(allocator, x, y, op, ctx),
            .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .numeric => unreachable,
        }
    } else if (comptime !types.isMatrix(Y)) {
        switch (comptime types.matrixType(X)) {
            .general => return general.apply2(allocator, x, y, op, ctx),
            .symmetric => return symmetric.apply2(allocator, x, y, op, ctx),
            .hermitian => return hermitian.apply2(allocator, x, y, op, ctx),
            .triangular => return triangular.apply2(allocator, x, y, op, ctx),
            .diagonal => return diagonal.apply2(allocator, x, y, op, ctx),
            .banded => return banded.apply2(allocator, x, y, op, ctx),
            .tridiagonal => return tridiagonal.apply2(allocator, x, y, op, ctx),
            .permutation => return permutation.apply2(allocator, x, y, op, ctx),
            .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.matrixType(X)) {
            .general => switch (comptime types.matrixType(Y)) {
                .general => return general.apply2(allocator, x, y, op, ctx),
                .symmetric => return gesy.apply2(allocator, x, y, op, ctx),
                .hermitian => return gehe.apply2(allocator, x, y, op, ctx),
                .triangular => return getr.apply2(allocator, x, y, op, ctx),
                .diagonal => return gedi.apply2(allocator, x, y, op, ctx),
                .banded => return geba.apply2(allocator, x, y, op, ctx),
                .tridiagonal => return gegt.apply2(allocator, x, y, op, ctx),
                .permutation => return gepe.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .symmetric => switch (comptime types.matrixType(Y)) {
                .general => return syge.apply2(allocator, x, y, op, ctx),
                .symmetric => return symmetric.apply2(allocator, x, y, op, ctx),
                .hermitian => return syhe.apply2(allocator, x, y, op, ctx),
                .triangular => return sytr.apply2(allocator, x, y, op, ctx),
                .diagonal => return sydi.apply2(allocator, x, y, op, ctx),
                .banded => return syba.apply2(allocator, x, y, op, ctx),
                .tridiagonal => return sygt.apply2(allocator, x, y, op, ctx),
                .permutation => return sype.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .hermitian => switch (comptime types.matrixType(Y)) {
                .general => return hege.apply2(allocator, x, y, op, ctx),
                .symmetric => return hesy.apply2(allocator, x, y, op, ctx),
                .hermitian => return hermitian.apply2(allocator, x, y, op, ctx),
                .triangular => return hetr.apply2(allocator, x, y, op, ctx),
                .diagonal => return hedi.apply2(allocator, x, y, op, ctx),
                .banded => return heba.apply2(allocator, x, y, op, ctx),
                .tridiagonal => return hegt.apply2(allocator, x, y, op, ctx),
                .permutation => return hepe.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .triangular => switch (comptime types.matrixType(Y)) {
                .general => return trge.apply2(allocator, x, y, op, ctx),
                .symmetric => return trsy.apply2(allocator, x, y, op, ctx),
                .hermitian => return trhe.apply2(allocator, x, y, op, ctx),
                .triangular => return triangular.apply2(allocator, x, y, op, ctx),
                .diagonal => return trdi.apply2(allocator, x, y, op, ctx),
                .banded => return trba.apply2(allocator, x, y, op, ctx),
                .tridiagonal => return trgt.apply2(allocator, x, y, op, ctx),
                .permutation => return trpe.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .diagonal => switch (comptime types.matrixType(Y)) {
                .general => return dige.apply2(allocator, x, y, op, ctx),
                .symmetric => return disy.apply2(allocator, x, y, op, ctx),
                .hermitian => return dihe.apply2(allocator, x, y, op, ctx),
                .triangular => return ditr.apply2(allocator, x, y, op, ctx),
                .diagonal => return diagonal.apply2(allocator, x, y, op, ctx),
                .banded => return diba.apply2(allocator, x, y, op, ctx),
                .tridiagonal => return digt.apply2(allocator, x, y, op, ctx),
                .permutation => return dipe.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .banded => switch (comptime types.matrixType(Y)) {
                .general => return bage.apply2(allocator, x, y, op, ctx),
                .symmetric => return basy.apply2(allocator, x, y, op, ctx),
                .hermitian => return bahe.apply2(allocator, x, y, op, ctx),
                .triangular => return batr.apply2(allocator, x, y, op, ctx),
                .diagonal => return badi.apply2(allocator, x, y, op, ctx),
                .banded => return banded.apply2(allocator, x, y, op, ctx),
                .tridiagonal => return bagt.apply2(allocator, x, y, op, ctx),
                .permutation => return bape.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .tridiagonal => switch (comptime types.matrixType(Y)) {
                .general => return gtge.apply2(allocator, x, y, op, ctx),
                .symmetric => return gtsy.apply2(allocator, x, y, op, ctx),
                .hermitian => return gthe.apply2(allocator, x, y, op, ctx),
                .triangular => return gttr.apply2(allocator, x, y, op, ctx),
                .diagonal => return gtdi.apply2(allocator, x, y, op, ctx),
                .banded => return gtba.apply2(allocator, x, y, op, ctx),
                .tridiagonal => return tridiagonal.apply2(allocator, x, y, op, ctx),
                .permutation => return gtpe.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .permutation => switch (comptime types.matrixType(Y)) {
                .general => return pege.apply2(allocator, x, y, op, ctx),
                .symmetric => return pesy.apply2(allocator, x, y, op, ctx),
                .hermitian => return pehe.apply2(allocator, x, y, op, ctx),
                .triangular => return petr.apply2(allocator, x, y, op, ctx),
                .diagonal => return pedi.apply2(allocator, x, y, op, ctx),
                .banded => return peba.apply2(allocator, x, y, op, ctx),
                .tridiagonal => return pegt.apply2(allocator, x, y, op, ctx),
                .permutation => return permutation.apply2(allocator, x, y, op, ctx),
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
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

    comptime if (!types.isMatrix(@TypeOf(x)) or !types.isMatrix(@TypeOf(y)))
        @compileError("Both arguments to add must be matrix types");

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

    comptime if (!types.isMatrix(@TypeOf(x)) or !types.isMatrix(@TypeOf(y)))
        @compileError("Both arguments to sub must be matrix types");

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

    comptime if (!types.isMatrix(X) and !types.isMatrix(Y))
        @compileError("At least one of the arguments must be a matrix type");

    if (comptime (types.isMatrix(X) and types.isMatrix(Y)) or
        types.isVector(X) or types.isVector(Y))
    { // matrix * matrix  or  vector * matrix  or  matrix * vector
        comptime if (types.isArbitraryPrecision(C)) {
            @compileError("Arbitrary precision types not implemented yet");
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        };

        @compileError("Matrix multiplication not implemented yet");
        // return linalg.matmul(...);
    } else {
        comptime if (types.isArbitraryPrecision(C)) { // scalar * matrix  or  matrix * scalar
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

    comptime if (!types.isMatrix(X) and types.isMatrix(Y))
        @compileError("First argument must be a matrix type and second argument must be a scalar type");

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
