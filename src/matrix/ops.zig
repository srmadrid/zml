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
    const R: type = ReturnType2(op, Numeric(X), Numeric(Y));

    comptime if (!types.isMatrix(X) and !types.isMatrix(Y))
        @compileError("apply2: at least one of x or y must be a matrix, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op)).@"fn".params.len != 3))
        @compileError("apply2: op must be a function of two arguments, or a function of three arguments with tdhe third argument being a context, got " ++ @typeName(@TypeOf(op)));

    if (comptime !types.isMatrix(X)) {
        switch (comptime types.matrixType(Y)) {
            .dense_general => return @import("ops/dge.zig").apply2(allocator, x, y, op, ctx),
            .dense_symmetric => return @import("ops/dsy.zig").apply2(allocator, x, y, op, ctx),
            .dense_hermitian => return @import("ops/dhe.zig").apply2(allocator, x, y, op, ctx),
            .dense_triangular => return @import("ops/dtr.zig").apply2(allocator, x, y, op, ctx),
            .dense_diagonal => return @import("ops/ddi.zig").apply2(allocator, x, y, op, ctx),
            .dense_banded => return @import("ops/dba.zig").apply2(allocator, x, y, op, ctx),
            .dense_tridiagonal => return @import("ops/dgt.zig").apply2(allocator, x, y, op, ctx),
            .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
            .permutation => @compileError("apply2 not implemented for permutation matrices yet"),
            .numeric => unreachable,
        }
    } else if (comptime !types.isMatrix(Y)) {
        switch (comptime types.matrixType(X)) {
            .dense_general => return @import("ops/dge.zig").apply2(allocator, x, y, op, ctx),
            .dense_symmetric => return @import("ops/dsy.zig").apply2(allocator, x, y, op, ctx),
            .dense_hermitian => return @import("ops/dhe.zig").apply2(allocator, x, y, op, ctx),
            .dense_triangular => return @import("ops/dtr.zig").apply2(allocator, x, y, op, ctx),
            .dense_diagonal => return @import("ops/ddi.zig").apply2(allocator, x, y, op, ctx),
            .dense_banded => return @import("ops/dba.zig").apply2(allocator, x, y, op, ctx),
            .dense_tridiagonal => return @import("ops/dgt.zig").apply2(allocator, x, y, op, ctx),
            .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
            .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
            .permutation => @compileError("apply2 not implemented for permutation matrices yet"),
            .numeric => unreachable,
        }
    } else {
        const m: u32 = if (comptime types.isSquareMatrix(X)) x.size else x.rows;
        const n: u32 = if (comptime types.isSquareMatrix(X)) x.size else x.cols;

        if (m != (if (comptime types.isSquareMatrix(Y)) y.size else y.rows))
            return linalg.Error.DimensionMismatch;

        if (n != (if (comptime types.isSquareMatrix(Y)) y.size else y.cols))
            return linalg.Error.DimensionMismatch;

        switch (comptime types.matrixType(X)) {
            .dense_general => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dgedsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dgedhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dgedtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dgeddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dgedba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dgedgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_symmetric => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_hermitian => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_triangular => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_banded => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_symmetric => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_hermitian => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .permutation => return @import("ops/dgepe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_symmetric => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dsydge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dsydhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dsydtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dsyddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dsydba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dsydgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_symmetric => {
                    var result: matrix.dense.Symmetric(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDSY(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_hermitian => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.dense.Hermitian(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .sparse_triangular => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_banded => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_symmetric => {
                    var result: matrix.dense.Symmetric(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDSY(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_hermitian => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.dense.Hermitian(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .permutation => return @import("ops/dsype.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_hermitian => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dhedge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dhedsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dhedtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dheddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dhedba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dhedgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_symmetric => {
                    if (comptime !types.isComplex(Numeric(Y))) {
                        var result: matrix.dense.Hermitian(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .sparse_hermitian => {
                    var result: matrix.dense.Hermitian(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDHE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_triangular => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_banded => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_symmetric => {
                    if (comptime !types.isComplex(Numeric(Y))) {
                        var result: matrix.dense.Hermitian(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .sparse_block_hermitian => {
                    var result: matrix.dense.Hermitian(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDHE(&result, x, y, op, ctx);

                    return result;
                },
                .permutation => return @import("ops/dhepe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_triangular => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dtrdge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dtrdsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dtrdhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dtrddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => {
                    var result: matrix.dense.Banded(R, types.orderOf(X)) = try .init(
                        allocator,
                        m,
                        n,
                        if (comptime types.uploOf(X) == .upper) y.lower else m - 1,
                        if (comptime types.uploOf(X) == .upper) n - 1 else y.upper,
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowDBA(&result, x, y, op, ctx);

                    return result;
                },
                .dense_tridiagonal => return @import("ops/dtrdgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_symmetric => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_hermitian => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_triangular => {
                    if (comptime types.uploOf(X) == types.uploOf(Y)) {
                        var result: matrix.dense.Triangular(R, types.uploOf(X), .non_unit, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDTR(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_symmetric => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_hermitian => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .permutation => return @import("ops/dtrpe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_diagonal => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/ddidge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/ddidsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/ddidhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/ddidtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/ddi.zig").apply2(allocator, x, y, op, ctx),
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
                .permutation => return @import("ops/ddipe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_banded => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dbadge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dbadsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dbadhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => {
                    var result: matrix.dense.Banded(R, types.orderOf(X)) = try .init(
                        allocator,
                        m,
                        n,
                        if (comptime types.uploOf(Y) == .upper) x.lower else m - 1,
                        if (comptime types.uploOf(Y) == .upper) n - 1 else x.upper,
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowDBA(&result, x, y, op, ctx);

                    return result;
                },
                .dense_diagonal => return @import("ops/dbaddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dbadgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_symmetric => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_hermitian => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_triangular => {
                    var result: matrix.dense.Banded(R, types.orderOf(X)) = try .init(
                        allocator,
                        m,
                        n,
                        if (comptime types.uploOf(Y) == .upper) x.lower else m - 1,
                        if (comptime types.uploOf(Y) == .upper) n - 1 else x.upper,
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowDBA(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_symmetric => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .sparse_block_hermitian => {
                    var result: matrix.dense.General(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .permutation => return @import("ops/dbape.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .dense_tridiagonal => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/dgtdge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/dgtdsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/dgtdhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/dgtdtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/dgtddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/dgtdba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/dgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => {
                    var result: matrix.sparse.Builder(R, types.orderOf(X)) = try .init(allocator, m, n, 3 * m - 2 + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .sparse_symmetric => {
                    var result: matrix.sparse.Builder(R, types.orderOf(X)) = try .init(allocator, m, n, 3 * m - 2 + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => return @import("ops/dgtpe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .permutation => switch (comptime types.matrixType(Y)) {
                .dense_general => return @import("ops/pedge.zig").apply2(allocator, x, y, op, ctx),
                .dense_symmetric => return @import("ops/pedsy.zig").apply2(allocator, x, y, op, ctx),
                .dense_hermitian => return @import("ops/pedhe.zig").apply2(allocator, x, y, op, ctx),
                .dense_triangular => return @import("ops/pedtr.zig").apply2(allocator, x, y, op, ctx),
                .dense_diagonal => return @import("ops/peddi.zig").apply2(allocator, x, y, op, ctx),
                .dense_banded => return @import("ops/pedba.zig").apply2(allocator, x, y, op, ctx),
                .dense_tridiagonal => return @import("ops/pedgt.zig").apply2(allocator, x, y, op, ctx),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for two permutation matrices yet"),
                .numeric => unreachable,
            },
            .sparse_general => switch (comptime types.matrixType(Y)) {
                .dense_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse_symmetric => switch (comptime types.matrixType(Y)) {
                .dense_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse_hermitian => switch (comptime types.matrixType(Y)) {
                .dense_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse_triangular => switch (comptime types.matrixType(Y)) {
                .dense_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse_banded => switch (comptime types.matrixType(Y)) {
                .dense_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse_block_general => switch (comptime types.matrixType(Y)) {
                .dense_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse_block_symmetric => switch (comptime types.matrixType(Y)) {
                .dense_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse_block_hermitian => switch (comptime types.matrixType(Y)) {
                .dense_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .dense_tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_triangular => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_general => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_symmetric => @compileError("apply2 not implemented for sparse matrices yet"),
                .sparse_block_hermitian => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .numeric => unreachable,
        }
    }
}

fn defaultSlowDGE(result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    const m: u32 = result.rows;
    const n: u32 = result.cols;
    const opinfo = @typeInfo(@TypeOf(op));

    var j: u32 = 0;
    while (j < n) : (j += 1) {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                try result.set(i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
            } else if (comptime opinfo.@"fn".params.len == 3) {
                try result.set(i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
            }
        }
    }
}

fn defaultSlowDSY(result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    const n: u32 = result.size;
    const opinfo = @typeInfo(@TypeOf(op));

    if (comptime types.uploOf(types.Child(@TypeOf(result))) == .upper) {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var i: u32 = 0;
            while (i <= j) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    try result.set(i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try result.set(i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
                }
            }
        }
    } else {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var i: u32 = j;
            while (i < n) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    try result.set(i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try result.set(i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
                }
            }
        }
    }
}

fn defaultSlowDHE(result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    const n: u32 = result.size;
    const opinfo = @typeInfo(@TypeOf(op));

    if (comptime types.uploOf(types.Child(@TypeOf(result))) == .upper) {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var i: u32 = 0;
            while (i <= j) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    try result.set(i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try result.set(i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
                }
            }
        }
    } else {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var i: u32 = j;
            while (i < n) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    try result.set(i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try result.set(i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
                }
            }
        }
    }
}

fn defaultSlowDTR(result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    // Always non-unit triangular
    const m: u32 = result.rows;
    const n: u32 = result.cols;
    const opinfo = @typeInfo(@TypeOf(op));

    if (comptime types.uploOf(types.Child(@TypeOf(x))) == .upper) {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var i: u32 = 0;
            while (i <= int.min(j, m - 1)) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    try result.set(i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try result.set(i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
                }
            }
        }
    } else {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var i: u32 = j;
            while (i < m) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    try result.set(i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try result.set(i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
                }
            }
        }
    }
}

fn defaultSlowDBA(result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    const m: u32 = result.rows;
    const n: u32 = result.cols;
    const kl: u32 = result.lower;
    const ku: u32 = result.upper;
    const opinfo = @typeInfo(@TypeOf(op));

    var j: u32 = 0;
    while (j < n) : (j += 1) {
        var i: u32 = if (j < ku) 0 else j - ku;
        while (i <= int.min(m + 1, j + kl)) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                try result.set(i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
            } else if (comptime opinfo.@"fn".params.len == 3) {
                try result.set(i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
            }
        }
    }
}

fn defaultSlowSGE(comptime order: types.Order, allocator: std.mem.Allocator, result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    const m: u32 = result.rows;
    const n: u32 = result.cols;
    const opinfo = @typeInfo(@TypeOf(op));

    if (comptime order == .col_major) {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var i: u32 = 0;
            while (i < m) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    try result.set(allocator, i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try result.set(allocator, i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
                }
            }
        }
    } else {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    try result.set(allocator, i, j, op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable));
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try result.set(allocator, i, j, try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx));
                }
            }
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
