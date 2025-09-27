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

const matrix = @import("../matrix.zig");
const dgeneral = @import("dense/general.zig");
const dsymmetric = @import("dense/symmetric.zig");
const dhermitian = @import("dense/hermitian.zig");
const dtriangular = @import("dense/triangular.zig");
const ddiagonal = @import("dense/diagonal.zig");
const dbanded = @import("dense/banded.zig");
const dtridiagonal = @import("dense/tridiagonal.zig");
const sgeneral = @import("sparse/general.zig");
const ssymmetric = @import("sparse/symmetric.zig");
const shermitian = @import("sparse/hermitian.zig");
const striangular = @import("sparse/triangular.zig");
const sbanded = @import("sparse/banded.zig");
const sblockgeneral = @import("sparse/block/general.zig");
const sblocksymmetric = @import("sparse/block/symmetric.zig");
const sblockhermitian = @import("sparse/block/hermitian.zig");
const spermutation = @import("sparse/permutation.zig");

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
            .dense_general => return dgeneral.apply2(allocator, x, y, op, ctx),
            .dense_symmetric => return dsymmetric.apply2(allocator, x, y, op, ctx),
            .dense_hermitian => return dhermitian.apply2(allocator, x, y, op, ctx),
            .dense_triangular => return dtriangular.apply2(allocator, x, y, op, ctx),
            .dense_diagonal => return ddiagonal.apply2(allocator, x, y, op, ctx),
            .dense_banded => return dbanded.apply2(allocator, x, y, op, ctx),
            .dense_tridiagonal => return dtridiagonal.apply2(allocator, x, y, op, ctx),
            .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_permutation => return spermutation.apply2(allocator, x, y, op, ctx),
            .numeric => unreachable,
        }
    } else if (comptime !types.isMatrix(Y)) {
        switch (comptime types.matrixType(X)) {
            .dense_general => return dgeneral.apply2(allocator, x, y, op, ctx),
            .dense_symmetric => return dsymmetric.apply2(allocator, x, y, op, ctx),
            .dense_hermitian => return dhermitian.apply2(allocator, x, y, op, ctx),
            .dense_triangular => return dtriangular.apply2(allocator, x, y, op, ctx),
            .dense_diagonal => return ddiagonal.apply2(allocator, x, y, op, ctx),
            .dense_banded => return dbanded.apply2(allocator, x, y, op, ctx),
            .dense_tridiagonal => return dtridiagonal.apply2(allocator, x, y, op, ctx),
            .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_permutation => return spermutation.apply2(allocator, x, y, op, ctx),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.matrixType(X)) {
            .dense_general => switch (comptime types.matrixType(Y)) {
                .dense_general => return dgeneral.apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return gesy.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return gehe.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return getr.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return gedi.apply2(allocator, x, y, op, ctx),
                .dense_banded => return geba.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return gegt.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return gepe.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_symmetric => switch (comptime types.matrixType(Y)) {
                .dense_general => return syge.apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return dsymmetric.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return syhe.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return sytr.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return sydi.apply2(allocator, x, y, op, ctx),
                .dense_banded => return syba.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return sygt.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return sype.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_hermitian => switch (comptime types.matrixType(Y)) {
                .dense_general => return hege.apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return hesy.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return dhermitian.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return hetr.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return hedi.apply2(allocator, x, y, op, ctx),
                .dense_banded => return heba.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return hegt.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return hepe.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_triangular => switch (comptime types.matrixType(Y)) {
                .dense_general => return trge.apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return trsy.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return trhe.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return dtriangular.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return trdi.apply2(allocator, x, y, op, ctx),
                .dense_banded => return trba.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return trgt.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return trpe.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_diagonal => switch (comptime types.matrixType(Y)) {
                .dense_general => return dige.apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return disy.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return dihe.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return ditr.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return ddiagonal.apply2(allocator, x, y, op, ctx),
                .dense_banded => return diba.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return digt.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return dipe.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_banded => switch (comptime types.matrixType(Y)) {
                .dense_general => return bage.apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return basy.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return bahe.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return batr.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return badi.apply2(allocator, x, y, op, ctx),
                .dense_banded => return dbanded.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return bagt.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return bape.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_tridiagonal => switch (comptime types.matrixType(Y)) {
                .dense_general => return gtge.apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return gtsy.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return gthe.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return gttr.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return gtdi.apply2(allocator, x, y, op, ctx),
                .dense_banded => return gtba.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return dtridiagonal.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return gtpe.apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .sparse_permutation => switch (comptime types.matrixType(Y)) {
                .dense_general => return pege.apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return pesy.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return pehe.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return petr.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return pedi.apply2(allocator, x, y, op, ctx),
                .dense_banded => return peba.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return pegt.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return spermutation.apply2(allocator, x, y, op, ctx),
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

        return linalg.matmul(allocator, x, y, ctx);
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
