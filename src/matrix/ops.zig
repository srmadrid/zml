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
        @compileError("apply2: op must be a function of two arguments, or a function of three arguments with tdhe third argument being a context, got " ++ @typeName(@TypeOf(op)));

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
                .dense_symmetric => return @import("ops/dgedsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dgedhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dgedtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dgeddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dgedba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dgedgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return @import("ops/dgespe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_symmetric => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dsydge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return dsymmetric.apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dsydhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dsydtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dsyddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dsydba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dsydgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return @import("ops/dsyspe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_hermitian => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dhedge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dhedsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return dhermitian.apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dhedtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dheddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dhedba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dhedgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return @import("ops/dhespe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_triangular => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dtrdge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dtrdsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dtrdhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return dtriangular.apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dtrddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dtrdba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dtrdgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return @import("ops/dtrspe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_diagonal => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/ddidge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/ddidsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/ddidhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/ddidtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return ddiagonal.apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/ddidba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/ddidgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return @import("ops/ddispe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_banded => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dbadge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dbadsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dbadhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dbadtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dbaddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return dbanded.apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dbadgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return @import("ops/dbaspe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_tridiagonal => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dgtdge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dgtdsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dgtdhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dgtdtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dgtddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dgtdba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return dtridiagonal.apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_permutation => return @import("ops/dgtspe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .sparse_permutation => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/spedge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/spedsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/spedhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/spedtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/speddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/spedba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/spedgt.zig").apply2(allocator, x, y, op, ctx),
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
            else => @compileError("apply2 not implemented for matrix type " ++ @typeName(X) ++ " yet"),
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

pub inline fn ddiv(
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
