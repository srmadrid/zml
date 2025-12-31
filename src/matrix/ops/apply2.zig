const std = @import("std");

const types = @import("../../types.zig");
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

const matops = @import("../ops.zig");

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
            .general_dense => return @import("apply2/dge.zig").apply2(allocator, x, y, op, ctx),
            .symmetric_dense => return @import("apply2/dsy.zig").apply2(allocator, x, y, op, ctx),
            .hermitian_dense => return @import("apply2/dhe.zig").apply2(allocator, x, y, op, ctx),
            .triangular_dense => return @import("apply2/dtr.zig").apply2(allocator, x, y, op, ctx),
            .diagonal => return @import("apply2/ddi.zig").apply2(allocator, x, y, op, ctx),
            .banded => return @import("apply2/dba.zig").apply2(allocator, x, y, op, ctx),
            .tridiagonal => return @import("apply2/dgt.zig").apply2(allocator, x, y, op, ctx),
            .general_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .symmetric_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .hermitian_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .triangular_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .general_block => @compileError("apply2 not implemented for sparse matrices yet"),
            .symmetric_block => @compileError("apply2 not implemented for sparse matrices yet"),
            .hermitian_block => @compileError("apply2 not implemented for sparse matrices yet"),
            .permutation => @compileError("apply2 not implemented for permutation matrices yet"),
            .numeric => unreachable,
        }
    } else if (comptime !types.isMatrix(Y)) {
        switch (comptime types.matrixType(X)) {
            .general_dense => return @import("apply2/dge.zig").apply2(allocator, x, y, op, ctx),
            .symmetric_dense => return @import("apply2/dsy.zig").apply2(allocator, x, y, op, ctx),
            .hermitian_dense => return @import("apply2/dhe.zig").apply2(allocator, x, y, op, ctx),
            .triangular_dense => return @import("apply2/dtr.zig").apply2(allocator, x, y, op, ctx),
            .diagonal => return @import("apply2/ddi.zig").apply2(allocator, x, y, op, ctx),
            .banded => return @import("apply2/dba.zig").apply2(allocator, x, y, op, ctx),
            .tridiagonal => return @import("apply2/dgt.zig").apply2(allocator, x, y, op, ctx),
            .general_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .symmetric_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .hermitian_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .triangular_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .general_block => @compileError("apply2 not implemented for sparse matrices yet"),
            .symmetric_block => @compileError("apply2 not implemented for sparse matrices yet"),
            .hermitian_block => @compileError("apply2 not implemented for sparse matrices yet"),
            .permutation => @compileError("apply2 not implemented for permutation matrices yet"),
            .numeric => unreachable,
        }
    } else {
        const m: u32 = if (comptime types.isSquareMatrix(X)) x.size else x.rows;
        const n: u32 = if (comptime types.isSquareMatrix(X)) x.size else x.cols;

        if (m != (if (comptime types.isSquareMatrix(Y)) y.size else y.rows))
            return matrix.Error.DimensionMismatch;

        if (n != (if (comptime types.isSquareMatrix(Y)) y.size else y.cols))
            return matrix.Error.DimensionMismatch;

        switch (comptime types.matrixType(X)) {
            .general_dense => switch (comptime types.matrixType(Y)) {
                .general_dense => return @import("apply2/dge.zig").apply2(allocator, x, y, op, ctx),
                .symmetric_dense => return @import("apply2/dgedsy.zig").apply2(allocator, x, y, op, ctx),
                .hermitian_dense => return @import("apply2/dgedhe.zig").apply2(allocator, x, y, op, ctx),
                .triangular_dense => return @import("apply2/dgedtr.zig").apply2(allocator, x, y, op, ctx),
                .general_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .triangular_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .general_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .diagonal => return @import("apply2/dgeddi.zig").apply2(allocator, x, y, op, ctx),
                .banded => return @import("apply2/dgedba.zig").apply2(allocator, x, y, op, ctx),
                .tridiagonal => return @import("apply2/dgedgt.zig").apply2(allocator, x, y, op, ctx),
                .permutation => return @import("apply2/dgepe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .symmetric_dense => switch (comptime types.matrixType(Y)) {
                .general_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_dense => return @import("apply2/dsy.zig").apply2(allocator, x, y, op, ctx),
                .hermitian_dense => return @import("apply2/dsydhe.zig").apply2(allocator, x, y, op, ctx),
                .triangular_dense => return @import("apply2/dsydtr.zig").apply2(allocator, x, y, op, ctx),
                .general_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_sparse => {
                    var result: matrix.symmetric.Dense(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDSY(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_sparse => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.hermitian.Dense(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .triangular_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .general_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_block => {
                    var result: matrix.symmetric.Dense(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDSY(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_block => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.hermitian.Dense(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .diagonal => return @import("apply2/dsyddi.zig").apply2(allocator, x, y, op, ctx),
                .banded => return @import("apply2/dsydba.zig").apply2(allocator, x, y, op, ctx),
                .tridiagonal => return @import("apply2/dsydgt.zig").apply2(allocator, x, y, op, ctx),
                .permutation => return @import("apply2/dsype.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .hermitian_dense => switch (comptime types.matrixType(Y)) {
                .general_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_dense => return @import("apply2/dhedsy.zig").apply2(allocator, x, y, op, ctx),
                .hermitian_dense => return @import("apply2/dhe.zig").apply2(allocator, x, y, op, ctx),
                .triangular_dense => return @import("apply2/dhedtr.zig").apply2(allocator, x, y, op, ctx),
                .general_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_sparse => {
                    if (comptime !types.isComplex(Numeric(Y))) {
                        var result: matrix.hermitian.Dense(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .hermitian_sparse => {
                    var result: matrix.hermitian.Dense(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDHE(&result, x, y, op, ctx);

                    return result;
                },
                .triangular_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .general_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_block => {
                    if (comptime !types.isComplex(Numeric(Y))) {
                        var result: matrix.hermitian.Dense(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .hermitian_block => {
                    var result: matrix.hermitian.Dense(R, types.uploOf(X), types.orderOf(X)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDHE(&result, x, y, op, ctx);

                    return result;
                },
                .diagonal => return @import("apply2/dheddi.zig").apply2(allocator, x, y, op, ctx),
                .banded => return @import("apply2/dhedba.zig").apply2(allocator, x, y, op, ctx),
                .tridiagonal => return @import("apply2/dhedgt.zig").apply2(allocator, x, y, op, ctx),
                .permutation => return @import("apply2/dhepe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .triangular_dense => switch (comptime types.matrixType(Y)) {
                .general_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .triangular_dense => return @import("apply2/dtr.zig").apply2(allocator, x, y, op, ctx),
                .general_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .triangular_sparse => {
                    if (comptime types.uploOf(X) == types.uploOf(Y)) {
                        var result: matrix.triangular.Dense(R, types.uploOf(X), .non_unit, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDTR(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .general_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .diagonal => return @import("apply2/dtrddi.zig").apply2(allocator, x, y, op, ctx),
                .banded => {
                    var result: matrix.Banded(R, types.orderOf(X)) = try .init(
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
                .tridiagonal => return @import("apply2/dtrdgt.zig").apply2(allocator, x, y, op, ctx),
                .permutation => return @import("apply2/dtrpe.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .general_sparse => switch (comptime types.matrixType(Y)) {
                .general_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .triangular_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .general_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + 2 * y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .hermitian_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + 2 * y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .triangular_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + y.nnz + (if (comptime types.diagOf(Y) == .unit) int.min(m, n) else 0));
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .general_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + 2 * y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .hermitian_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + 2 * y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .diagonal => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + int.min(m, n));
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .banded => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .tridiagonal => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + 3 * n - 2);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return result;
                },
                .permutation => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + n);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .numeric => unreachable,
            },
            .symmetric_sparse => switch (comptime types.matrixType(Y)) {
                .general_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_dense => {
                    var result: matrix.symmetric.Dense(R, types.uploOf(Y), types.orderOf(Y)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDSY(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_dense => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.hermitian.Dense(R, types.uploOf(Y), types.orderOf(Y)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .triangular_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .general_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 2 * x.nnz + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSSY(types.uploOf(X), types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compileSymmetric(allocator, types.uploOf(X));
                },
                .hermitian_sparse => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + y.nnz);
                        errdefer result.deinit(allocator);

                        try defaultSlowSHE(types.uploOf(X), types.orderOf(X), allocator, &result, x, y, op, ctx);

                        return try result.compileHermitian(allocator, types.uploOf(X));
                    } else {
                        var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 2 * x.nnz + 2 * y.nnz);
                        errdefer result.deinit(allocator);

                        try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                        return try result.compile(allocator);
                    }
                },
                .triangular_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 2 * x.nnz + y.nnz + (if (comptime types.diagOf(Y) == .unit) int.min(m, n) else 0));
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .general_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 2 * x.nnz + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSSY(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compileSymmetric(allocator, types.uploOf(X));
                },
                .hermitian_block => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + y.nnzb * y.bsize * y.bsize);
                        errdefer result.deinit(allocator);

                        try defaultSlowSHE(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                        return try result.compileHermitian(allocator, types.uploOf(Y));
                    } else {
                        var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 2 * x.nnz + 2 * y.nnzb * y.bsize * y.bsize);
                        errdefer result.deinit(allocator);

                        try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                        return try result.compile(allocator);
                    }
                },
                .diagonal => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + int.min(m, n));
                    errdefer result.deinit(allocator);

                    try defaultSlowSSY(types.uploOf(X), types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compileSymmetric(allocator, types.uploOf(X));
                },
                .banded => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .tridiagonal => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 2 * x.nnz + 3 * n - 2);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .permutation => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + n);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .numeric => unreachable,
            },
            .hermitian_sparse => switch (comptime types.matrixType(Y)) {
                .general_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_dense => {
                    if (comptime !types.isComplex(Numeric(Y))) {
                        var result: matrix.hermitian.Dense(R, types.uploOf(Y), types.orderOf(Y)) = try .init(allocator, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDHE(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .hermitian_dense => {
                    var result: matrix.hermitian.Dense(R, types.uploOf(Y), types.orderOf(Y)) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDHE(&result, x, y, op, ctx);

                    return result;
                },
                .triangular_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .general_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 2 * x.nnz + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_sparse => {
                    if (comptime !types.isComplex(Numeric(Y))) {
                        var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + y.nnz);
                        errdefer result.deinit(allocator);

                        try defaultSlowSSY(types.uploOf(X), types.orderOf(X), allocator, &result, x, y, op, ctx);

                        return try result.compileSymmetric(allocator, types.uploOf(X));
                    } else {
                        var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 2 * x.nnz + 2 * y.nnz);
                        errdefer result.deinit(allocator);

                        try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                        return try result.compile(allocator);
                    }
                },
                .hermitian_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSHE(types.uploOf(X), types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compileHermitian(allocator, types.uploOf(X));
                },
                .triangular_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 2 * x.nnz + y.nnz + (if (comptime types.diagOf(Y) == .unit) int.min(m, n) else 0));
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .general_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 2 * x.nnz + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_block => {
                    if (comptime !types.isComplex(Numeric(Y))) {
                        var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + y.nnzb * y.bsize * y.bsize);
                        errdefer result.deinit(allocator);

                        try defaultSlowSHE(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                        return try result.compileHermitian(allocator, types.uploOf(Y));
                    } else {
                        var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 2 * x.nnz + 2 * y.nnzb * y.bsize * y.bsize);
                        errdefer result.deinit(allocator);

                        try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                        return try result.compile(allocator);
                    }
                },
                .hermitian_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSHE(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compileHermitian(allocator, types.uploOf(Y));
                },
                .diagonal => {
                    if (comptime !types.isComplex(Numeric(Y))) {
                        var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + int.min(m, n));
                        errdefer result.deinit(allocator);

                        try defaultSlowSHE(types.uploOf(X), types.orderOf(X), allocator, &result, x, y, op, ctx);

                        return try result.compileHermitian(allocator, types.uploOf(X));
                    } else {
                        var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 2 * x.nnz + int.min(m, n));
                        errdefer result.deinit(allocator);

                        try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                        return try result.compile(allocator);
                    }
                },
                .banded => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .tridiagonal => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 2 * x.nnz + 3 * n - 2);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .permutation => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, x.nnz + n);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .numeric => unreachable,
            },
            .triangular_sparse => switch (comptime types.matrixType(Y)) {
                .general_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_dense => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .triangular_dense => {
                    if (comptime types.uploOf(X) == types.uploOf(Y)) {
                        var result: matrix.triangular.Dense(R, types.uploOf(Y), .non_unit, types.orderOf(Y)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDTR(&result, x, y, op, ctx);

                        return result;
                    } else {
                        var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowDGE(&result, x, y, op, ctx);

                        return result;
                    }
                },
                .general_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + (if (comptime types.diagOf(X) == .unit) int.min(m, n) else 0) + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + (if (comptime types.diagOf(X) == .unit) int.min(m, n) else 0) + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .hermitian_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + (if (comptime types.diagOf(X) == .unit) int.min(m, n) else 0) + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .triangular_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .general_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + (if (comptime types.diagOf(X) == .unit) int.min(m, n) else 0) + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + (if (comptime types.diagOf(X) == .unit) int.min(m, n) else 0) + 2 * y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .hermitian_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, x.nnz + (if (comptime types.diagOf(X) == .unit) int.min(m, n) else 0) + 2 * y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .general_block => switch (comptime types.matrixType(Y)) {
                .general_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .triangular_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .general_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .triangular_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .general_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .symmetric_block => switch (comptime types.matrixType(Y)) {
                .general_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .triangular_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .general_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .triangular_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .general_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .hermitian_block => switch (comptime types.matrixType(Y)) {
                .general_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .triangular_dense => @compileError("apply2 not implemented for sparse matrices yet"),
                .general_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .triangular_sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .general_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .symmetric_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .hermitian_block => @compileError("apply2 not implemented for sparse matrices yet"),
                .diagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .banded => @compileError("apply2 not implemented for sparse matrices yet"),
                .tridiagonal => @compileError("apply2 not implemented for sparse matrices yet"),
                .permutation => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .diagonal => switch (comptime types.matrixType(Y)) {
                .general_dense => return @import("apply2/ddidge.zig").apply2(allocator, x, y, op, ctx),
                .symmetric_dense => return @import("apply2/ddidsy.zig").apply2(allocator, x, y, op, ctx),
                .hermitian_dense => return @import("apply2/ddidhe.zig").apply2(allocator, x, y, op, ctx),
                .triangular_dense => return @import("apply2/ddidtr.zig").apply2(allocator, x, y, op, ctx),
                .general_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, int.min(m, n) + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return result.compile(allocator);
                },
                .symmetric_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, int.min(m, n) + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSSY(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return result.compileSymmetric(allocator, types.uploOf(Y));
                },
                .hermitian_sparse => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, int.min(m, n) + y.nnz);
                        errdefer result.deinit(allocator);

                        try defaultSlowSHE(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                        return result.compileHermitian(allocator, types.uploOf(Y));
                    } else {
                        var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, int.min(m, n) + 2 * y.nnz);
                        errdefer result.deinit(allocator);

                        try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                        return result.compile(allocator);
                    }
                },
                .triangular_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, (if (comptime types.diagOf(Y) == .unit) 0 else int.min(m, n)) + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSTR(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return result.compileTriangular(allocator, types.uploOf(Y), .non_unit);
                },
                .general_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, int.min(m, n) + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return result.compile(allocator);
                },
                .symmetric_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, int.min(m, n) + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSSY(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return result.compileSymmetric(allocator, types.uploOf(Y));
                },
                .hermitian_block => {
                    if (comptime !types.isComplex(Numeric(X))) {
                        var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, int.min(m, n) + y.nnzb * y.bsize * y.bsize);
                        errdefer result.deinit(allocator);

                        try defaultSlowSHE(types.uploOf(Y), types.orderOf(Y), allocator, &result, x, y, op, ctx);

                        return result.compileHermitian(allocator, types.uploOf(Y));
                    } else {
                        var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, int.min(m, n) + 2 * y.nnzb * y.bsize * y.bsize);
                        errdefer result.deinit(allocator);

                        try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                        return result.compile(allocator);
                    }
                },
                .diagonal => return @import("apply2/ddi.zig").apply2(allocator, x, y, op, ctx),
                .banded => return @import("apply2/ddidba.zig").apply2(allocator, x, y, op, ctx),
                .tridiagonal => return @import("apply2/ddidgt.zig").apply2(allocator, x, y, op, ctx),
                .permutation => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, int.min(m, n) + n);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .numeric => unreachable,
            },
            .banded => switch (comptime types.matrixType(Y)) {
                .general_dense => return @import("apply2/dbadge.zig").apply2(allocator, x, y, op, ctx),
                .symmetric_dense => return @import("apply2/dbadsy.zig").apply2(allocator, x, y, op, ctx),
                .hermitian_dense => return @import("apply2/dbadhe.zig").apply2(allocator, x, y, op, ctx),
                .triangular_dense => {
                    var result: matrix.Banded(R, types.orderOf(Y)) = try .init(
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
                .general_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_sparse => {
                    var result: matrix.general.Dense(R, types.orderOf(Y)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .triangular_sparse => {
                    var result: matrix.Banded(R, types.orderOf(Y)) = try .init(
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
                .general_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .symmetric_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .hermitian_block => {
                    var result: matrix.general.Dense(R, types.orderOf(X)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowDGE(&result, x, y, op, ctx);

                    return result;
                },
                .diagonal => return @import("apply2/dbaddi.zig").apply2(allocator, x, y, op, ctx),
                .banded => return @import("apply2/dba.zig").apply2(allocator, x, y, op, ctx),
                .tridiagonal => return @import("apply2/dbadgt.zig").apply2(allocator, x, y, op, ctx),
                .permutation => return @import("apply2/dbape.zig").apply2(allocator, x, y, op, ctx),
                .numeric => unreachable,
            },
            .tridiagonal => switch (comptime types.matrixType(Y)) {
                .general_dense => return @import("apply2/dgtdge.zig").apply2(allocator, x, y, op, ctx),
                .symmetric_dense => return @import("apply2/dgtdsy.zig").apply2(allocator, x, y, op, ctx),
                .hermitian_dense => return @import("apply2/dgtdhe.zig").apply2(allocator, x, y, op, ctx),
                .triangular_dense => return @import("apply2/dgtdtr.zig").apply2(allocator, x, y, op, ctx),
                .general_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 3 * m - 2 + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 3 * m - 2 + 2 * y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .hermitian_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 3 * m - 2 + 2 * y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .triangular_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, (if (comptime types.diagOf(Y) == .unit) 2 else 3) * m - 2 + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .general_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 3 * m - 2 + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 3 * m - 2 + 2 * y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .hermitian_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, 3 * m - 2 + 2 * y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .diagonal => return @import("apply2/dgtddi.zig").apply2(allocator, x, y, op, ctx),
                .banded => return @import("apply2/dgtdba.zig").apply2(allocator, x, y, op, ctx),
                .tridiagonal => return @import("apply2/dgt.zig").apply2(allocator, x, y, op, ctx),
                .permutation => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 3 * m - 2 + n);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .numeric => unreachable,
            },
            .permutation => switch (comptime types.matrixType(Y)) {
                .general_dense => return @import("apply2/pedge.zig").apply2(allocator, x, y, op, ctx),
                .symmetric_dense => return @import("apply2/pedsy.zig").apply2(allocator, x, y, op, ctx),
                .hermitian_dense => return @import("apply2/pedhe.zig").apply2(allocator, x, y, op, ctx),
                .triangular_dense => return @import("apply2/pedtr.zig").apply2(allocator, x, y, op, ctx),
                .general_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, n + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, n + 2 * y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .hermitian_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, n + 2 * y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .triangular_sparse => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, n + y.nnz);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .general_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, n + y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .symmetric_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, n + 2 * y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .hermitian_block => {
                    var result: matrix.builder.Sparse(R, types.orderOf(Y)) = try .init(allocator, m, n, n + 2 * y.nnzb * y.bsize * y.bsize);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(Y), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .diagonal => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, n + int.min(m, n));
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .banded => return @import("apply2/pedba.zig").apply2(allocator, x, y, op, ctx),
                .tridiagonal => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, n + 3 * n - 2);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
                .permutation => {
                    var result: matrix.builder.Sparse(R, types.orderOf(X)) = try .init(allocator, m, n, 2 * n);
                    errdefer result.deinit(allocator);

                    try defaultSlowSGE(types.orderOf(X), allocator, &result, x, y, op, ctx);

                    return try result.compile(allocator);
                },
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

    if (comptime types.uploOf(types.Child(@TypeOf(result))) == .upper) {
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
                    const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                    if (try ops.ne(rij, 0, ctx)) {
                        try result.set(allocator, i, j, rij);
                    }
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                    if (try ops.ne(rij, 0, ctx)) {
                        try result.set(allocator, i, j, rij);
                    }
                }
            }
        }
    } else {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                    if (try ops.ne(rij, 0, ctx)) {
                        try result.set(allocator, i, j, rij);
                    }
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                    if (try ops.ne(rij, 0, ctx)) {
                        try result.set(allocator, i, j, rij);
                    }
                }
            }
        }
    }
}

fn defaultSlowSSY(comptime uplo: types.Uplo, comptime order: types.Order, allocator: std.mem.Allocator, result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    const n: u32 = result.rows;
    const opinfo = @typeInfo(@TypeOf(op));

    if (comptime uplo == .upper) {
        if (comptime order == .col_major) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                var i: u32 = 0;
                while (i <= j) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        } else {
            var i: u32 = 0;
            while (i < n) : (i += 1) {
                var j: u32 = i;
                while (j < n) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        }
    } else {
        if (comptime order == .col_major) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                var i: u32 = j;
                while (i < n) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        } else {
            var i: u32 = 0;
            while (i < n) : (i += 1) {
                var j: u32 = 0;
                while (j <= i) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        }
    }
}

fn defaultSlowSHE(comptime uplo: types.Uplo, comptime order: types.Order, allocator: std.mem.Allocator, result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    const n: u32 = result.rows;
    const opinfo = @typeInfo(@TypeOf(op));

    if (comptime uplo == .upper) {
        if (comptime order == .col_major) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                var i: u32 = 0;
                while (i <= j) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        } else {
            var i: u32 = 0;
            while (i < n) : (i += 1) {
                var j: u32 = i;
                while (j < n) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        }
    } else {
        if (comptime order == .col_major) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                var i: u32 = j;
                while (i < n) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        } else {
            var i: u32 = 0;
            while (i < n) : (i += 1) {
                var j: u32 = 0;
                while (j <= i) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        }
    }
}

fn defaultSlowSTR(comptime uplo: types.Uplo, comptime order: types.Order, allocator: std.mem.Allocator, result: anytype, x: anytype, y: anytype, comptime op: anytype, ctx: anytype) !void {
    // Always non-unit triangular
    const m: u32 = result.rows;
    const n: u32 = result.cols;
    const opinfo = @typeInfo(@TypeOf(op));

    if (comptime uplo == .upper) {
        if (comptime order == .col_major) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                var i: u32 = 0;
                while (i <= int.min(j, m - 1)) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        } else {
            var i: u32 = 0;
            while (i < m) : (i += 1) {
                var j: u32 = i;
                while (j < n) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        }
    } else {
        if (comptime order == .col_major) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                var i: u32 = j;
                while (i < m) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        } else {
            var i: u32 = 0;
            while (i < m) : (i += 1) {
                var j: u32 = 0;
                while (j <= int.min(i, n - 1)) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        const rij = op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        const rij = try op(x.get(i, j) catch unreachable, y.get(i, j) catch unreachable, ctx);

                        if (try ops.ne(rij, 0, ctx)) {
                            try result.set(allocator, i, j, rij);
                        }
                    }
                }
            }
        }
    }
}
