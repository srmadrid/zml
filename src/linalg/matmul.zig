const std = @import("std");

const types = @import("../types.zig");
const Numeric = types.Numeric;
const Coerce = types.Coerce;
const MulCoerce = types.MulCoerce;
const int = @import("../int.zig");
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");

const linalg = @import("../linalg.zig");

const blas = @import("blas.zig");
const lapack = @import("lapack.zig");

pub inline fn matmul(allocator: std.mem.Allocator, a: anytype, b: anytype, ctx: anytype) !MulCoerce(@TypeOf(a), @TypeOf(b)) {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = Coerce(Numeric(A), Numeric(B));

    comptime if (!((types.isMatrix(A) and types.isMatrix(B)) or
        (types.isMatrix(A) and types.isVector(B)) or
        (types.isVector(A) and types.isMatrix(B))))
        @compileError("matmul: at least one argument must be a matrix, the other must be a matrix or vector, got " ++ @typeName(A) ++ " and " ++ @typeName(B));

    comptime if (types.isArbitraryPrecision(Numeric(A)) or types.isArbitraryPrecision(Numeric(B))) {
        // When implemented, expand if
        @compileError("zml.linalg.matmul not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime !types.isMatrix(A)) { // vector * matrix
        const m: u32 = if (comptime types.isSquareMatrix(B)) b.size else b.rows;
        const n: u32 = if (comptime types.isSquareMatrix(B)) b.size else b.cols;

        if (a.len != m)
            return linalg.Error.DimensionMismatch;

        switch (comptime types.vectorType(A)) {
            .dense => switch (comptime types.matrixType(B)) {
                .dense_general => return @import("matmul/dvdge.zig").vm(allocator, a, b, ctx), // dense vector * dense general matrix
                .dense_symmetric => return @import("matmul/dvdsy.zig").vm(allocator, a, b, ctx), // dense vector * dense symmetric matrix
                .dense_hermitian => return @import("matmul/dvdhe.zig").vm(allocator, a, b, ctx), // dense vector * dense hermitian matrix
                .dense_triangular => return @import("matmul/dvdtr.zig").vm(allocator, a, b, ctx), // dense vector * dense triangular matrix
                .dense_diagonal => return @import("matmul/dvddi.zig").vm(allocator, a, b, ctx), // dense vector * dense diagonal matrix
                .dense_banded => return @import("matmul/dvdba.zig").vm(allocator, a, b, ctx), // dense vector * dense banded matrix
                .dense_tridiagonal => { // dense vector * dense tridiagonal matrix
                    var result: vector.Dense(C) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    // lapack.lagtm
                    try defaultSlowVM(&result, a, b, ctx);

                    return result;
                },
                .permutation => return @import("matmul/dvpe.zig").vm(allocator, a, b, ctx), // dense vector * permutation matrix
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .sparse => switch (comptime types.matrixType(B)) {
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .numeric => unreachable,
        }
    } else if (comptime !types.isMatrix(B)) { // matrix * vector
        const m: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.rows;
        const n: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.cols;

        if (b.len != n)
            return linalg.Error.DimensionMismatch;

        switch (comptime types.matrixType(A)) {
            .dense_general => switch (comptime types.vectorType(B)) {
                .dense => return @import("matmul/dgedv.zig").mv(allocator, a, b, ctx), // dense general matrix * dense vector
                .sparse => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_symmetric => switch (comptime types.vectorType(B)) {
                .dense => return @import("matmul/dsydv.zig").mv(allocator, a, b, ctx), // dense symmetric matrix * dense vector
                .sparse => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_hermitian => switch (comptime types.vectorType(B)) {
                .dense => return @import("matmul/dhedv.zig").mv(allocator, a, b, ctx), // dense hermitian matrix * dense vector
                .sparse => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_triangular => switch (comptime types.vectorType(B)) {
                .dense => return @import("matmul/dtrdv.zig").mv(allocator, a, b, ctx), // dense triangular matrix * dense vector
                .sparse => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_diagonal => switch (comptime types.vectorType(B)) {
                .dense => return @import("matmul/ddidv.zig").mv(allocator, a, b, ctx), // dense diagonal matrix * dense vector
                .sparse => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_banded => switch (comptime types.vectorType(B)) {
                .dense => return @import("matmul/dbadv.zig").mv(allocator, a, b, ctx), // dense banded matrix * dense vector
                .sparse => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_tridiagonal => switch (comptime types.vectorType(B)) {
                .dense => { // dense tridiagonal matrix * dense vector
                    var result: vector.Dense(C) = try .init(allocator, m);
                    errdefer result.deinit(allocator);

                    // lapack.lagtm
                    try defaultSlowMV(&result, a, b, ctx);

                    return result;
                },
                .sparse => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .permutation => switch (comptime types.vectorType(B)) {
                .dense => return @import("matmul/pedv.zig").mv(allocator, a, b, ctx), // permutation matrix * dense vector
                .sparse => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
            .numeric => unreachable,
        }
    } else {
        const m: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.rows;
        const k: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.cols;
        const n: u32 = if (comptime types.isSquareMatrix(B)) b.size else b.cols;

        if (k != (if (comptime types.isSquareMatrix(B)) b.size else b.rows))
            return linalg.Error.DimensionMismatch;

        switch (comptime types.matrixType(A)) {
            .dense_general => switch (comptime types.matrixType(B)) {
                .dense_general => return @import("matmul/dgedge.zig").mm(allocator, a, b, ctx), // dense general matrix * dense general matrix
                .dense_symmetric => return @import("matmul/dgedsy.zig").mm(allocator, a, b, ctx), // dense general matrix * dense symmetric matrix
                .dense_hermitian => return @import("matmul/dgedhe.zig").mm(allocator, a, b, ctx), // dense general matrix * dense hermitian matrix
                .dense_triangular => { // dense general matrix * dense triangular matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .full(allocator, m, n, 0, ctx);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_diagonal => return @import("matmul/dgeddi.zig").mm(allocator, a, b, ctx), // dense general matrix * dense diagonal matrix
                .dense_banded => { // dense general matrix * dense banded matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_tridiagonal => { // dense general matrix * dense tridiagonal matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.lagtm
                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // dense general matrix * permutation matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.lapmt: convert first to inverse permutation and then 1 based indexing
                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_symmetric => switch (comptime types.matrixType(B)) {
                .dense_general => return @import("matmul/dsydge.zig").mm(allocator, a, b, ctx), // dense symmetric matrix * dense general matrix
                .dense_symmetric => { // dense symmetric matrix * dense symmetric matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_hermitian => { // dense symmetric matrix * dense hermitian matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_triangular => { // dense symmetric matrix * dense triangular matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_diagonal => { // dense symmetric matrix * dense diagonal matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_banded => { // dense symmetric matrix * dense banded matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_tridiagonal => { // dense symmetric matrix * dense tridiagonal matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // dense symmetric matrix * permutation matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_hermitian => switch (comptime types.matrixType(B)) {
                .dense_general => return @import("matmul/dhedge.zig").mm(allocator, a, b, ctx), // dense hermitian matrix * dense general matrix
                .dense_symmetric => { // dense hermitian matrix * dense symmetric matrix
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_hermitian => { // dense hermitian matrix * dense hermitian matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_triangular => { // dense hermitian matrix * dense triangular matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_diagonal => { // dense hermitian matrix * dense diagonal matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_banded => { // dense hermitian matrix * dense banded matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_tridiagonal => { // dense hermitian matrix * dense tridiagonal matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // dense hermitian matrix * permutation matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_triangular => switch (comptime types.matrixType(B)) {
                .dense_general => { // dense triangular matrix * dense general matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // Change to blas calls
                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_symmetric => { // dense triangular matrix * dense symmetric matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_hermitian => { // dense triangular matrix * dense hermitian matrix
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_triangular => { // dense triangular matrix * dense triangular matrix
                    if (comptime types.uploOf(A) == types.uploOf(B)) {
                        var result: matrix.dense.Triangular(
                            C,
                            types.uploOf(A),
                            if (types.diagOf(A) == types.diagOf(B))
                                types.diagOf(A)
                            else
                                .non_unit,
                            types.orderOf(A),
                        ) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowTMM(&result, a, b, ctx);

                        return result;
                    } else {
                        var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowMM(&result, a, b, ctx);

                        return result;
                    }
                },
                .dense_diagonal => { // dense triangular matrix * dense diagonal matrix
                    var result: matrix.dense.Triangular(C, types.uploOf(A), .non_unit, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowTMM(&result, a, b, ctx);

                    return result;
                },
                .dense_banded => { // dense triangular matrix * dense banded matrix
                    var result: matrix.dense.Banded(C, types.orderOf(A)) = try .init(
                        allocator,
                        m,
                        n,
                        // lC ≤ min{m − 1, lA ​+ lB​}, uC ​≤ min{n − 1, uA ​+ uB​}.​
                        if (comptime types.uploOf(A) == .upper) int.min(m - 1, b.lower) else m - 1,
                        if (comptime types.uploOf(A) == .upper) int.min(n - 1, k - 1 + b.upper) else int.min(n - 1, b.upper),
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .dense_tridiagonal => { // dense triangular matrix * dense tridiagonal matrix
                    var result: matrix.dense.Banded(C, types.orderOf(A)) = try .init(
                        allocator,
                        m,
                        n,
                        // lC ≤ min{m − 1, lA ​+ lB​}, uC ​≤ min{n − 1, uA ​+ uB​}.​
                        if (comptime types.uploOf(A) == .upper) 1 else m - 1,
                        if (comptime types.uploOf(A) == .upper) n - 1 else 1,
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // dense triangular matrix * permutation matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_diagonal => switch (comptime types.matrixType(B)) {
                .dense_general => return @import("matmul/ddidge.zig").mm(allocator, a, b, ctx), // dense diagonal matrix * dense general matrix
                .dense_symmetric => { // dense diagonal matrix * dense symmetric matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_hermitian => { // dense diagonal matrix * dense hermitian matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_triangular => { // dense diagonal matrix * dense triangular matrix
                    var result: matrix.dense.Triangular(C, types.uploOf(B), .non_unit, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowTMM(&result, a, b, ctx);

                    return result;
                },
                .dense_diagonal => return @import("matmul/ddiddi.zig").mm(allocator, a, b, ctx), // dense diagonal matrix * dense diagonal matrix
                .dense_banded => { // dense diagonal matrix * dense banded matrix
                    var result: matrix.dense.Banded(C, types.orderOf(B)) = try .init(allocator, m, n, b.lower, b.upper);
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .dense_tridiagonal => { // dense diagonal matrix * dense tridiagonal matrix
                    var result: matrix.dense.Banded(C, types.orderOf(B)) = try .init(allocator, m, n, 1, 1);
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // dense diagonal matrix * permutation matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_banded => switch (comptime types.matrixType(B)) {
                .dense_general => { // dense banded matrix * dense general matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_symmetric => { // dense banded matrix * dense symmetric matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_hermitian => { // dense banded matrix * dense hermitian matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_triangular => { // densed banded matrix * dense triangular matrix
                    var result: matrix.dense.Banded(C, types.orderOf(A)) = try .init(
                        allocator,
                        m,
                        n,
                        // lC ≤ min{m − 1, lA ​+ lB​}, uC ​≤ min{n − 1, uA ​+ uB​}.​
                        if (comptime types.uploOf(B) == .upper) int.min(m - 1, a.lower) else int.min(m - 1, a.lower + k - 1),
                        if (comptime types.uploOf(B) == .upper) n - 1 else int.min(n - 1, a.upper),
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .dense_diagonal => { // dense banded matrix * dense diagonal matrix
                    var result: matrix.dense.Banded(C, types.orderOf(A)) = try .init(allocator, m, n, a.lower, a.upper);
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .dense_banded => { // dense banded matrix * dense banded matrix
                    var result: matrix.dense.Banded(C, types.orderOf(A)) = try .init(
                        allocator,
                        m,
                        n,
                        // lC ≤ min{m − 1, lA ​+ lB​}, uC ​≤ min{n − 1, uA ​+ uB​}.​
                        int.min(m - 1, a.lower + b.lower),
                        int.min(n - 1, a.upper + b.upper),
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .dense_tridiagonal => { // dense banded matrix * dense tridiagonal matrix
                    var result: matrix.dense.Banded(C, types.orderOf(A)) = try .init(
                        allocator,
                        m,
                        n,
                        // lC ≤ min{m − 1, lA ​+ lB​}, uC ​≤ min{n − 1, uA ​+ uB​}.​
                        int.min(m - 1, a.lower + 1),
                        int.min(n - 1, a.upper + 1),
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // dense banded matrix * permutation matrix
                    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .dense_tridiagonal => switch (comptime types.matrixType(B)) {
                .dense_general => { // dense tridiagonal matrix * dense general matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.lagtm
                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_symmetric => { // dense tridiagonal matrix * dense symmetric matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_hermitian => { // dense tridiagonal matrix * dense hermitian matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .dense_triangular => { // dense tridiagonal matrix * dense triangular matrix
                    var result: matrix.dense.Banded(C, types.orderOf(B)) = try .init(
                        allocator,
                        m,
                        n,
                        // lC ≤ min{m − 1, lA ​+ lB​}, uC ​≤ min{n − 1, uA ​+ uB​}.​
                        if (comptime types.uploOf(B) == .upper) 1 else m - 1,
                        if (comptime types.uploOf(B) == .upper) n - 1 else 1,
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .dense_diagonal => { // dense tridiagonal matrix * dense diagonal matrix
                    var result: matrix.dense.Banded(C, types.orderOf(B)) = try .init(allocator, m, n, 1, 1);
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .dense_banded => { // dense tridiagonal matrix * dense banded matrix
                    var result: matrix.dense.Banded(C, types.orderOf(B)) = try .init(
                        allocator,
                        m,
                        n,
                        // lC ≤ min{m − 1, lA ​+ lB​}, uC ​≤ min{n − 1, uA ​+ uB​}.​
                        int.min(m - 1, 1 + b.lower),
                        int.min(n - 1, 1 + b.upper),
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .dense_tridiagonal => { // dense tridiagonal matrix * dense tridiagonal matrix
                    var result: matrix.dense.Banded(C, types.orderOf(B)) = try .init(
                        allocator,
                        n,
                        n,
                        // lC ≤ min{m − 1, lA ​+ lB​}, uC ​≤ min{n − 1, uA ​+ uB​}.​
                        int.min(n - 1, 2),
                        int.min(n - 1, 2),
                    );
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // dense tridiagonal matrix * permutation matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            .permutation => switch (comptime types.matrixType(B)) {
                .dense_general => { // permutation matrix * dense general matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.laswp
                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .dense_symmetric => { // permutation matrix * dense symmetric matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .dense_hermitian => { // permutation matrix * dense hermitian matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .dense_triangular => { // permutation matrix * dense triangular matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .dense_diagonal => { // permutation matrix * dense diagonal matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .dense_banded => { // permutation matrix * dense banded matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .dense_tridiagonal => { // permutation matrix * dense tridiagonal matrix
                    var result: matrix.dense.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // permutation matrix * permutation matrix
                    var result: matrix.Permutation(C) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    var i: u32 = 0;
                    while (i < n) : (i += 1) {
                        result.data[i] = b.data[a.data[i]];
                    }

                    return result;
                },
                else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
                .numeric => unreachable,
            },
            else => @compileError("matmul not implemented for " ++ @typeName(A) ++ " and " ++ @typeName(B)),
            .numeric => unreachable,
        }
    }
}

fn defaultSlowVM(result: anytype, x: anytype, a: anytype, ctx: anytype) !void {
    const A: type = @TypeOf(a);
    const C: type = Numeric(types.Child(@TypeOf(result)));

    const m: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.cols;

    var i: u32 = 0;
    while (i < n) : (i += 1) {
        try result.set(i, try constants.zero(C, ctx));

        var j: u32 = 0;
        while (j < m) : (j += 1) {
            try result.set(
                i,
                try ops.add(
                    result.get(i) catch unreachable,
                    try ops.mul(
                        a.get(j, i) catch unreachable,
                        x.get(j) catch unreachable,
                        ctx,
                    ),
                    ctx,
                ),
            );
        }
    }
}

fn defaultSlowMV(result: anytype, a: anytype, x: anytype, ctx: anytype) !void {
    const A: type = @TypeOf(a);
    const C: type = Numeric(types.Child(@TypeOf(result)));

    const m: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.cols;

    var i: u32 = 0;
    while (i < m) : (i += 1) {
        try result.set(i, try constants.zero(C, ctx));

        var j: u32 = 0;
        while (j < n) : (j += 1) {
            try result.set(
                i,
                try ops.add(
                    result.get(i) catch unreachable,
                    try ops.mul(
                        a.get(i, j) catch unreachable,
                        x.get(j) catch unreachable,
                        ctx,
                    ),
                    ctx,
                ),
            );
        }
    }
}

fn defaultSlowMM(result: anytype, a: anytype, b: anytype, ctx: anytype) !void {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = Numeric(types.Child(@TypeOf(result)));

    const m: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.rows;
    const k: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.cols;
    const n: u32 = if (comptime types.isSquareMatrix(B)) b.size else b.cols;

    var i: u32 = 0;
    while (i < m) : (i += 1) {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            try result.set(i, j, try constants.zero(C, ctx));

            var kk: u32 = 0;
            while (kk < k) : (kk += 1) {
                try result.set(
                    i,
                    j,
                    try ops.add(
                        result.get(i, j) catch unreachable,
                        try ops.mul(
                            a.get(i, kk) catch unreachable,
                            b.get(kk, j) catch unreachable,
                            ctx,
                        ),
                        ctx,
                    ),
                );
            }
        }
    }
}

fn defaultSlowTMM(result: anytype, a: anytype, b: anytype, ctx: anytype) !void {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = Numeric(types.Child(@TypeOf(result)));

    const m: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.rows;
    const k: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.cols;
    const n: u32 = if (comptime types.isSquareMatrix(B)) b.size else b.cols;

    if (comptime types.uploOf(types.Child(@TypeOf(result))) == .upper) {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = undefined;
            if (comptime types.diagOf(types.Child(@TypeOf(result))) == .unit) {
                j = i + 1;
            } else {
                j = i;
            }

            while (j < n) : (j += 1) {
                try result.set(i, j, try constants.zero(C, ctx));

                var kk: u32 = 0;
                while (kk < k) : (kk += 1) {
                    try result.set(
                        i,
                        j,
                        try ops.add(
                            result.get(i, j) catch unreachable,
                            try ops.mul(
                                a.get(i, kk) catch unreachable,
                                b.get(kk, j) catch unreachable,
                                ctx,
                            ),
                            ctx,
                        ),
                    );
                }
            }
        }
    } else {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = 0;
            while (j < if (comptime types.diagOf(types.Child(@TypeOf(result))) == .unit) i else i + 1) : (j += 1) {
                try result.set(i, j, try constants.zero(C, ctx));

                var kk: u32 = 0;
                while (kk < k) : (kk += 1) {
                    try result.set(
                        i,
                        j,
                        try ops.add(
                            result.get(i, j) catch unreachable,
                            try ops.mul(
                                a.get(i, kk) catch unreachable,
                                b.get(kk, j) catch unreachable,
                                ctx,
                            ),
                            ctx,
                        ),
                    );
                }
            }
        }
    }
}

fn defaultSlowBMM(result: anytype, a: anytype, b: anytype, ctx: anytype) !void {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = Numeric(types.Child(@TypeOf(result)));

    const m: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.rows;
    const k: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.cols;
    const n: u32 = if (comptime types.isSquareMatrix(B)) b.size else b.cols;

    var i: u32 = 0;
    while (i < m) : (i += 1) {
        var j: u32 = if (i < result.lower) 0 else i - result.lower;
        while (j <= int.min(n - 1, i + result.upper)) : (j += 1) {
            try result.set(i, j, try constants.zero(C, ctx));

            var kk: u32 = 0;
            while (kk < k) : (kk += 1) {
                try result.set(
                    i,
                    j,
                    try ops.add(
                        result.get(i, j) catch unreachable,
                        try ops.mul(
                            a.get(i, kk) catch unreachable,
                            b.get(kk, j) catch unreachable,
                            ctx,
                        ),
                        ctx,
                    ),
                );
            }
        }
    }
}

fn defaultSlowMP(result: anytype, a: anytype, b: anytype, ctx: anytype) !void {
    // Since P is on the right and we follow LAPACK convention (row permutation
    // array), we have to get the inverse of the permutation to correctly
    // permute the columns of A.
    const A: type = @TypeOf(a);

    const m: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSquareMatrix(A)) a.size else a.cols;

    if (b.direction == .forward) {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var k: u32 = 0;
            while (k < n) : (k += 1) {
                if (b.data[k] == j) {
                    break;
                }
            }

            var i: u32 = 0;
            while (i < m) : (i += 1) {
                try ops.set(
                    &result.data[
                        if (comptime types.orderOf(types.Child(@TypeOf(result))) == .col_major)
                            i + j * result.ld
                        else
                            i * result.ld + j
                    ],
                    a.get(i, k) catch unreachable,
                    ctx,
                );
            }
        }
    } else {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            var i: u32 = 0;
            while (i < m) : (i += 1) {
                try ops.set(
                    &result.data[
                        if (comptime types.orderOf(types.Child(@TypeOf(result))) == .col_major)
                            i + j * result.ld
                        else
                            i * result.ld + j
                    ],
                    a.get(i, b.data[j]) catch unreachable,
                    ctx,
                );
            }
        }
    }
}

fn defaultSlowPM(result: anytype, a: anytype, b: anytype, ctx: anytype) !void {
    const B: type = @TypeOf(b);

    const m: u32 = if (comptime types.isSquareMatrix(B)) b.size else b.rows;
    const n: u32 = if (comptime types.isSquareMatrix(B)) b.size else b.cols;

    if (a.direction == .forward) {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var j: u32 = 0;
            while (j < n) : (j += 1) {
                try ops.set(
                    &result.data[
                        if (comptime types.orderOf(types.Child(@TypeOf(result))) == .col_major)
                            i + j * result.ld
                        else
                            i * result.ld + j
                    ],
                    b.get(a.data[i], j) catch unreachable,
                    ctx,
                );
            }
        }
    } else {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            var k: u32 = 0;
            while (k < m) : (k += 1) {
                if (a.data[k] == i) {
                    break;
                }
            }

            var j: u32 = 0;
            while (j < n) : (j += 1) {
                try ops.set(
                    &result.data[
                        if (comptime types.orderOf(types.Child(@TypeOf(result))) == .col_major)
                            i + j * result.ld
                        else
                            i * result.ld + j
                    ],
                    b.get(k, j) catch unreachable,
                    ctx,
                );
            }
        }
    }
}
