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
        const m: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B) or types.isPermutationMatrix(B))
            b.size
        else
            b.rows;
        const n: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B) or types.isPermutationMatrix(B))
            b.size
        else
            b.cols;

        if (a.len != m)
            return linalg.Error.DimensionMismatch;

        switch (comptime types.matrixType(B)) {
            .general => { // vector * general
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                try blas.gemv(
                    types.orderOf(B),
                    .trans,
                    types.scast(i32, m),
                    types.scast(i32, n),
                    1,
                    b.data,
                    types.scast(i32, b.ld),
                    a.data,
                    a.inc,
                    0,
                    result.data,
                    result.inc,
                    ctx,
                );

                return result;
            },
            .symmetric => { // vector * symmetric
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                try blas.symv(
                    types.orderOf(B).invert(),
                    types.uploOf(B).invert(),
                    types.scast(i32, n),
                    1,
                    b.data,
                    types.scast(i32, b.ld),
                    a.data,
                    a.inc,
                    0,
                    result.data,
                    result.inc,
                    ctx,
                );

                return result;
            },
            .hermitian => { // vector * hermitian
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                try blas.hemv(
                    types.orderOf(B).invert(),
                    types.uploOf(B).invert(),
                    types.scast(i32, n),
                    1,
                    b.data,
                    types.scast(i32, b.ld),
                    a.data,
                    a.inc,
                    0,
                    result.data,
                    result.inc,
                    ctx,
                );

                return result;
            },
            .triangular => { // vector * triangular
                var result: vector.Vector(C) = try .full(allocator, n, 0, ctx);
                errdefer result.deinit(allocator);

                if (m == n) {
                    try blas.copy(
                        types.scast(i32, n),
                        a.data,
                        a.inc,
                        result.data,
                        result.inc,
                        ctx,
                    );

                    try blas.trmv(
                        types.orderOf(B),
                        types.uploOf(B),
                        .trans,
                        types.diagOf(B),
                        types.scast(i32, n),
                        b.data,
                        types.scast(i32, b.ld),
                        result.data,
                        result.inc,
                        ctx,
                    );
                } else {
                    const min_dim: u32 = int.min(m, n);

                    try blas.copy(
                        types.scast(i32, min_dim),
                        a.data,
                        a.inc,
                        result.data,
                        result.inc,
                        ctx,
                    );

                    try blas.trmv(
                        types.orderOf(B),
                        types.uploOf(B),
                        .trans,
                        types.diagOf(B),
                        types.scast(i32, min_dim),
                        b.data,
                        types.scast(i32, b.ld),
                        result.data,
                        result.inc,
                        ctx,
                    );

                    if (comptime types.uploOf(B) == .upper) {
                        if (n > min_dim) {
                            try blas.gemv(
                                types.orderOf(B),
                                .trans,
                                types.scast(i32, m),
                                types.scast(i32, n - min_dim),
                                1,
                                b.data +
                                    if (comptime types.orderOf(B) == .col_major)
                                        min_dim * b.ld
                                    else
                                        min_dim,
                                types.scast(i32, b.ld),
                                a.data,
                                a.inc,
                                0,
                                result.data + min_dim * types.scast(u32, result.inc),
                                result.inc,
                                ctx,
                            );
                        }
                    } else {
                        if (m > min_dim) {
                            try blas.gemv(
                                types.orderOf(B),
                                .trans,
                                types.scast(i32, m - min_dim),
                                types.scast(i32, n),
                                1,
                                b.data +
                                    if (comptime types.orderOf(B) == .col_major)
                                        min_dim
                                    else
                                        min_dim * b.ld,
                                types.scast(i32, b.ld),
                                a.data +
                                    if (a.inc > 0)
                                        min_dim * types.scast(u32, a.inc)
                                    else
                                        0,
                                a.inc,
                                1,
                                result.data,
                                result.inc,
                                ctx,
                            );
                        } else if (n > min_dim) {
                            try blas.scal(
                                types.scast(i32, n - min_dim),
                                0,
                                result.data + min_dim * types.scast(u32, result.inc),
                                result.inc,
                                ctx,
                            );
                        }
                    }
                }

                return result;
            },
            .diagonal => { // vector * diagonal
                var result: vector.Vector(C) = try .full(allocator, n, 0, ctx);
                errdefer result.deinit(allocator);

                try blas.copy(
                    types.scast(i32, m),
                    a.data,
                    a.inc,
                    result.data,
                    result.inc,
                    ctx,
                );

                var i: u32 = 0;
                while (i < int.min(m, n)) : (i += 1) {
                    try ops.mul_( // result[i] *= b[i, i]
                        &result.data[i],
                        result.data[i],
                        b.data[i],
                        ctx,
                    );
                }

                while (i < n) : (i += 1) {
                    result.data[i] = try constants.zero(C, ctx);
                }

                return result;
            },
            .banded => { // vector * banded
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                try blas.gbmv(
                    types.orderOf(B),
                    .trans,
                    types.scast(i32, m),
                    types.scast(i32, n),
                    types.scast(i32, b.lower),
                    types.scast(i32, b.upper),
                    1,
                    b.data,
                    types.scast(i32, b.ld),
                    a.data,
                    a.inc,
                    0,
                    result.data,
                    result.inc,
                    ctx,
                );

                return result;
            },
            .tridiagonal => { // vector * tridiagonal
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                // lapack.lagtm
                try defaultSlowVM(&result, a, b, ctx);

                return result;
            },
            .permutation => { // vector * permutation
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                var i: u32 = 0;
                while (i < n) : (i += 1) {
                    try ops.set(
                        &result.data[b.data[i]],
                        a.get(i) catch unreachable,
                        ctx,
                    );
                }

                return result;
            },
            .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .numeric => unreachable,
        }
    } else if (comptime !types.isMatrix(B)) { // matrix * vector
        const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A) or types.isPermutationMatrix(A))
            a.size
        else
            a.rows;
        const n: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A) or types.isPermutationMatrix(A))
            a.size
        else
            a.cols;

        if (b.len != n)
            return linalg.Error.DimensionMismatch;

        switch (comptime types.matrixType(A)) {
            .general => { // general * vector
                var result: vector.Vector(C) = try .init(allocator, m);
                errdefer result.deinit(allocator);

                try blas.gemv(
                    types.orderOf(A),
                    .no_trans,
                    types.scast(i32, m),
                    types.scast(i32, n),
                    1,
                    a.data,
                    types.scast(i32, a.ld),
                    b.data,
                    b.inc,
                    0,
                    result.data,
                    result.inc,
                    ctx,
                );

                return result;
            },
            .symmetric => { // vector * symmetric
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                try blas.symv(
                    types.orderOf(A),
                    types.uploOf(A),
                    types.scast(i32, n),
                    1,
                    a.data,
                    types.scast(i32, a.ld),
                    b.data,
                    b.inc,
                    0,
                    result.data,
                    result.inc,
                    ctx,
                );

                return result;
            },
            .hermitian => { // vector * hermitian
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                try blas.hemv(
                    types.orderOf(A),
                    types.uploOf(A),
                    types.scast(i32, n),
                    1,
                    a.data,
                    types.scast(i32, a.ld),
                    b.data,
                    b.inc,
                    0,
                    result.data,
                    result.inc,
                    ctx,
                );

                return result;
            },
            .triangular => { // triangular * vector
                var result: vector.Vector(C) = try .full(allocator, m, 0, ctx);
                errdefer result.deinit(allocator);

                if (m == n) {
                    try blas.copy(
                        types.scast(i32, n),
                        b.data,
                        b.inc,
                        result.data,
                        result.inc,
                        ctx,
                    );

                    try blas.trmv(
                        types.orderOf(A),
                        types.uploOf(A),
                        .no_trans,
                        types.diagOf(A),
                        types.scast(i32, n),
                        a.data,
                        types.scast(i32, a.ld),
                        result.data,
                        result.inc,
                        ctx,
                    );
                } else {
                    const min_dim: u32 = int.min(m, n);

                    try blas.copy(
                        types.scast(i32, min_dim),
                        b.data,
                        b.inc,
                        result.data,
                        result.inc,
                        ctx,
                    );

                    try blas.trmv(
                        types.orderOf(A),
                        types.uploOf(A),
                        .no_trans,
                        types.diagOf(A),
                        types.scast(i32, min_dim),
                        a.data,
                        types.scast(i32, a.ld),
                        result.data,
                        result.inc,
                        ctx,
                    );

                    if (comptime types.uploOf(A) == .upper) {
                        if (n > min_dim) { // extra rows (filled)
                            try blas.gemv(
                                types.orderOf(A),
                                .no_trans,
                                types.scast(i32, m),
                                types.scast(i32, n - min_dim),
                                1,
                                a.data +
                                    if (comptime types.orderOf(A) == .col_major)
                                        min_dim * a.ld
                                    else
                                        min_dim,
                                types.scast(i32, a.ld),
                                b.data +
                                    if (b.inc > 0)
                                        min_dim * types.scast(u32, b.inc)
                                    else
                                        0,
                                b.inc,
                                1,
                                result.data,
                                result.inc,
                                ctx,
                            );
                        } else if (m > min_dim) { // extra rows (empty)
                            try blas.scal(
                                types.scast(i32, m - min_dim),
                                0,
                                result.data + min_dim * types.scast(u32, result.inc),
                                result.inc,
                                ctx,
                            );
                        }
                    } else {
                        if (m > min_dim) { // extra rows (filled)
                            try blas.gemv(
                                types.orderOf(A),
                                .no_trans,
                                types.scast(i32, m - min_dim),
                                types.scast(i32, n),
                                1,
                                a.data +
                                    if (comptime types.orderOf(A) == .col_major)
                                        min_dim
                                    else
                                        min_dim * a.ld,
                                types.scast(i32, a.ld),
                                b.data,
                                b.inc,
                                1,
                                result.data + min_dim * types.scast(u32, result.inc),
                                result.inc,
                                ctx,
                            );
                        }
                    }
                }

                return result;
            },
            .diagonal => { // diagonal * vector
                var result: vector.Vector(C) = try .full(allocator, m, 0, ctx);
                errdefer result.deinit(allocator);

                try blas.copy(
                    types.scast(i32, m),
                    b.data,
                    b.inc,
                    result.data,
                    result.inc,
                    ctx,
                );

                var i: u32 = 0;
                while (i < int.min(m, n)) : (i += 1) {
                    try ops.mul_( // result[i] *= a[i, i]
                        &result.data[i],
                        result.data[i],
                        a.data[i],
                        ctx,
                    );
                }

                while (i < m) : (i += 1) {
                    result.data[i] = try constants.zero(C, ctx);
                }

                return result;
            },
            .banded => { // banded * vector
                var result: vector.Vector(C) = try .init(allocator, m);
                errdefer result.deinit(allocator);

                try blas.gbmv(
                    types.orderOf(A),
                    .no_trans,
                    types.scast(i32, m),
                    types.scast(i32, n),
                    types.scast(i32, a.lower),
                    types.scast(i32, a.upper),
                    1,
                    a.data,
                    types.scast(i32, a.ld),
                    b.data,
                    b.inc,
                    0,
                    result.data,
                    result.inc,
                    ctx,
                );

                return result;
            },
            .tridiagonal => { // tridiagonal * vector
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                // lapack.lagtm
                try defaultSlowMV(&result, a, b, ctx);

                return result;
            },
            .permutation => { // permutation * vector
                var result: vector.Vector(C) = try .init(allocator, n);
                errdefer result.deinit(allocator);

                var i: u32 = 0;
                while (i < n) : (i += 1) {
                    try ops.set(
                        &result.data[i],
                        b.get(a.data[i]) catch unreachable,
                        ctx,
                    );
                }

                return result;
            },
            .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .numeric => unreachable,
        }
    } else {
        const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A) or types.isPermutationMatrix(A))
            a.size
        else
            a.rows;
        const k: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A) or types.isPermutationMatrix(A))
            a.size
        else
            a.cols;
        const n: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B) or types.isPermutationMatrix(B))
            b.size
        else
            b.cols;

        if (k != (if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B) or types.isPermutationMatrix(B))
            b.size
        else
            b.rows))
            return linalg.Error.DimensionMismatch;

        switch (comptime types.matrixType(A)) {
            .general => switch (comptime types.matrixType(B)) {
                .general => { // general * general
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try blas.gemm(
                        types.orderOf(A),
                        .no_trans,
                        if (comptime types.orderOf(A) == types.orderOf(B)) .no_trans else .trans,
                        types.scast(i32, m),
                        types.scast(i32, n),
                        types.scast(i32, k),
                        1,
                        a.data,
                        types.scast(i32, a.ld),
                        b.data,
                        types.scast(i32, b.ld),
                        0,
                        result.data,
                        types.scast(i32, result.ld),
                        ctx,
                    );

                    return result;
                },
                .symmetric => { // general * symmetric
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try blas.symm(
                        types.orderOf(A),
                        .right,
                        if (comptime types.orderOf(A) == types.orderOf(B)) types.uploOf(B) else types.uploOf(B).invert(),
                        types.scast(i32, m),
                        types.scast(i32, n),
                        1,
                        b.data,
                        types.scast(i32, b.ld),
                        a.data,
                        types.scast(i32, a.ld),
                        0,
                        result.data,
                        types.scast(i32, result.ld),
                        ctx,
                    );

                    return result;
                },
                .hermitian => { // general * hermitian
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try blas.hemm(
                        types.orderOf(A),
                        .right,
                        if (comptime types.orderOf(A) == types.orderOf(B)) types.uploOf(B) else types.uploOf(B).invert(),
                        types.scast(i32, m),
                        types.scast(i32, n),
                        1,
                        b.data,
                        types.scast(i32, b.ld),
                        a.data,
                        types.scast(i32, a.ld),
                        0,
                        result.data,
                        types.scast(i32, result.ld),
                        ctx,
                    );

                    return result;
                },
                .triangular => { // general * triangular
                    var result: matrix.General(C, types.orderOf(A)) = try .full(allocator, m, n, 0, ctx);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .diagonal => { // general * diagonal
                    var result: matrix.General(C, types.orderOf(A)) = try .full(allocator, m, n, 0, ctx);
                    errdefer result.deinit(allocator);

                    // lapack.lascl2

                    var j: u32 = 0;
                    while (j < int.min(k, n)) : (j += 1) {
                        try blas.axpy(
                            types.scast(i32, m),
                            b.data[j],
                            a.data + if (comptime types.orderOf(A) == .col_major) j * a.ld else j,
                            if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, a.ld),
                            result.data + if (comptime types.orderOf(A) == .col_major) j * result.ld else j,
                            if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, result.ld),
                            ctx,
                        );
                    }

                    while (j < n) : (j += 1) {
                        var i: u32 = 0;
                        while (i < m) : (i += 1) {
                            if (comptime types.orderOf(A) == .col_major) {
                                result.data[i + j * result.ld] = try constants.zero(C, ctx);
                            } else {
                                result.data[i * result.ld + j] = try constants.zero(C, ctx);
                            }
                        }
                    }

                    return result;
                },
                .banded => { // general * banded
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .tridiagonal => { // general * tridiagonal
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.lagtm
                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // general * permutation
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.lapmt: convert first to inverse permutation and then 1 based indexing
                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .symmetric => switch (comptime types.matrixType(B)) {
                .general => {
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    if (comptime types.orderOf(A) == types.orderOf(B)) {
                        try blas.symm(
                            types.orderOf(B),
                            .left,
                            types.uploOf(A),
                            types.scast(i32, m),
                            types.scast(i32, n),
                            1,
                            a.data,
                            types.scast(i32, a.ld),
                            b.data,
                            types.scast(i32, b.ld),
                            0,
                            result.data,
                            types.scast(i32, result.ld),
                            ctx,
                        );

                        return result;
                    } else {
                        var j: u32 = 0;
                        while (j < n) : (j += 1) {
                            try blas.symv(
                                types.orderOf(A),
                                types.uploOf(A),
                                types.scast(i32, m),
                                1,
                                a.data,
                                types.scast(i32, a.ld),
                                b.data +
                                    if (comptime types.orderOf(B) == .col_major)
                                        j * b.ld
                                    else
                                        j,
                                if (comptime types.orderOf(B) == .col_major) 1 else types.scast(i32, b.ld),
                                0,
                                result.data +
                                    if (comptime types.orderOf(A) == .col_major)
                                        j * result.ld
                                    else
                                        j,
                                if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, result.ld),
                                ctx,
                            );
                        }
                    }

                    return result;
                },
                .symmetric => { // symmetric * symmetric
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .hermitian => { // symmetric * hermitian
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .triangular => { // symmetric * triangular
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .diagonal => { // symmetric * diagonal
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .banded => { // symmetric * banded
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .tridiagonal => { // symmetric * tridiagonal
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // symmetric * permutation
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .hermitian => switch (comptime types.matrixType(B)) {
                .general => { // hermitian * general
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    if (comptime types.orderOf(A) == types.orderOf(B)) {
                        try blas.hemm(
                            types.orderOf(B),
                            .left,
                            types.uploOf(A),
                            types.scast(i32, m),
                            types.scast(i32, n),
                            1,
                            a.data,
                            types.scast(i32, a.ld),
                            b.data,
                            types.scast(i32, b.ld),
                            0,
                            result.data,
                            types.scast(i32, result.ld),
                            ctx,
                        );

                        return result;
                    } else {
                        var j: u32 = 0;
                        while (j < n) : (j += 1) {
                            try blas.hemv(
                                types.orderOf(A),
                                types.uploOf(A),
                                types.scast(i32, m),
                                1,
                                a.data,
                                types.scast(i32, a.ld),
                                b.data +
                                    if (comptime types.orderOf(B) == .col_major)
                                        j * b.ld
                                    else
                                        j,
                                if (comptime types.orderOf(B) == .col_major) 1 else types.scast(i32, b.ld),
                                0,
                                result.data +
                                    if (comptime types.orderOf(A) == .col_major)
                                        j * result.ld
                                    else
                                        j,
                                if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, result.ld),
                                ctx,
                            );
                        }
                    }

                    return result;
                },
                .symmetric => { // hermitian * symmetric
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .hermitian => { // hermitian * hermitian
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .triangular => { // hermitian * triangular
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .diagonal => { // hermitian * diagonal
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .banded => { // hermitian * banded
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .tridiagonal => { // hermitian * tridiagonal
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // hermitian * permutation
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .triangular => switch (comptime types.matrixType(B)) {
                .general => { // triangular * general
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // Change to blas calls
                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .symmetric => { // triangular * symmetric
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .hermitian => { // triangular * hermitian
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .triangular => { // triangular * triangular
                    if (comptime types.uploOf(A) == types.uploOf(B)) {
                        var result: matrix.Triangular(
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
                        var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                        errdefer result.deinit(allocator);

                        try defaultSlowMM(&result, a, b, ctx);

                        return result;
                    }
                },
                .diagonal => { // triangular * diagonal
                    var result: matrix.Triangular(C, types.uploOf(A), .non_unit, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowTMM(&result, a, b, ctx);

                    return result;
                },
                .banded => { // triangular * banded
                    var result: matrix.Banded(C, types.orderOf(A)) = try .init(
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
                .tridiagonal => { // triangular * tridiagonal
                    var result: matrix.Banded(C, types.orderOf(A)) = try .init(
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
                .permutation => { // triangular * permutation
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .diagonal => switch (comptime types.matrixType(B)) {
                .general => { // diagonal * general
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.lascl2

                    var i: u32 = 0;
                    while (i < int.min(m, k)) : (i += 1) {
                        try blas.axpy(
                            types.scast(i32, n),
                            a.data[i],
                            b.data + if (comptime types.orderOf(B) == .col_major) i * b.ld else i,
                            if (comptime types.orderOf(B) == .col_major) 1 else types.scast(i32, b.ld),
                            result.data + if (comptime types.orderOf(A) == .col_major) i * result.ld else i,
                            if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, result.ld),
                            ctx,
                        );
                    }

                    while (i < m) : (i += 1) {
                        var j: u32 = 0;
                        while (j < n) : (j += 1) {
                            if (comptime types.orderOf(A) == .col_major) {
                                result.data[i + j * result.ld] = try constants.zero(C, ctx);
                            } else {
                                result.data[i * result.ld + j] = try constants.zero(C, ctx);
                            }
                        }
                    }

                    return result;
                },
                .symmetric => { // diagonal * symmetric
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .hermitian => { // diagonal * hermitian
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .triangular => { // diagonal * triangular
                    var result: matrix.Triangular(C, types.uploOf(B), .non_unit, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowTMM(&result, a, b, ctx);

                    return result;
                },
                .diagonal => { // diagonal * diagonal
                    var result: matrix.Diagonal(C) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    var i: u32 = 0;
                    while (i < int.min(a.rows, b.cols)) : (i += 1) {
                        result.data[i] = try ops.mul(
                            a.data[i],
                            b.data[i],
                            ctx,
                        );
                    }

                    return result;
                },
                .banded => { // diagonal * banded
                    var result: matrix.Banded(C, types.orderOf(B)) = try .init(allocator, m, n, b.lower, b.upper);
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .tridiagonal => { // diagonal * tridiagonal
                    var result: matrix.Banded(C, types.orderOf(B)) = try .init(allocator, m, n, 1, 1);
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // diagonal * permutation
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .banded => switch (comptime types.matrixType(B)) {
                .general => { // banded * general
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .symmetric => { // banded * symmetric
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .hermitian => { // banded * hermitian
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .triangular => { // banded * triangular
                    var result: matrix.Banded(C, types.orderOf(A)) = try .init(
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
                .diagonal => { // banded * diagonal
                    var result: matrix.Banded(C, types.orderOf(A)) = try .init(allocator, m, n, a.lower, a.upper);
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .banded => { // banded * banded
                    var result: matrix.Banded(C, types.orderOf(A)) = try .init(
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
                .tridiagonal => { // banded * tridiagonal
                    var result: matrix.Banded(C, types.orderOf(A)) = try .init(
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
                .permutation => { // banded * permutation
                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .tridiagonal => switch (comptime types.matrixType(B)) {
                .general => { // tridiagonal * general
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.lagtm
                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .symmetric => { // tridiagonal * symmetric
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .hermitian => { // tridiagonal * hermitian
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .triangular => { // tridiagonal * triangular
                    var result: matrix.Banded(C, types.orderOf(B)) = try .init(
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
                .diagonal => { // tridiagonal * diagonal
                    var result: matrix.Banded(C, types.orderOf(B)) = try .init(allocator, m, n, 1, 1);
                    errdefer result.deinit(allocator);

                    try defaultSlowBMM(&result, a, b, ctx);

                    return result;
                },
                .banded => { // tridiagonal * banded
                    var result: matrix.Banded(C, types.orderOf(B)) = try .init(
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
                .tridiagonal => { // tridiagonal * tridiagonal
                    var result: matrix.Banded(C, types.orderOf(B)) = try .init(
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
                .permutation => { // tridiagonal * permutation
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowMP(&result, a, b, ctx);

                    return result;
                },
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .permutation => switch (comptime types.matrixType(B)) {
                .general => { // permutation * general
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    // lapack.laswp
                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .symmetric => { // permutation * symmetric
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .hermitian => { // permutation * hermitian
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .triangular => { // permutation * triangular
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .diagonal => { // permutation * diagonal
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .banded => { // permutation * banded
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, m, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .tridiagonal => { // permutation * tridiagonal
                    var result: matrix.General(C, types.orderOf(B)) = try .init(allocator, n, n);
                    errdefer result.deinit(allocator);

                    try defaultSlowPM(&result, a, b, ctx);

                    return result;
                },
                .permutation => { // permutation * permutation
                    var result: matrix.Permutation(C) = try .init(allocator, n);
                    errdefer result.deinit(allocator);

                    var i: u32 = 0;
                    while (i < n) : (i += 1) {
                        result.data[i] = b.data[a.data[i]];
                    }

                    return result;
                },
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
            .numeric => unreachable,
        }
    }
}

fn defaultSlowVM(result: anytype, x: anytype, a: anytype, ctx: anytype) !void {
    const A: type = @TypeOf(a);
    const C: type = Numeric(types.Child(@TypeOf(result)));

    const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.cols;

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

    const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.cols;

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

    const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A) or types.isPermutationMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B) or types.isPermutationMatrix(B)) b.size else b.cols;
    const k: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A) or types.isPermutationMatrix(A)) a.size else a.cols;

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

    const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B)) b.size else b.cols;
    const k: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.cols;

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

    const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B)) b.size else b.cols;
    const k: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.cols;

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

    const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.cols;

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

    const m: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B)) b.size else b.rows;
    const n: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B)) b.size else b.cols;

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
