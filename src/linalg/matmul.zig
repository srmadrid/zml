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

pub inline fn matmul(allocator: std.mem.Allocator, a: anytype, b: anytype, ctx: anytype) !MulCoerce(@TypeOf(a), @TypeOf(b)) {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = Coerce(Numeric(A), Numeric(B));

    comptime if (!((types.isMatrix(A) and types.isMatrix(B)) or
        (types.isMatrix(A) and types.isVector(B)) or
        (types.isVector(A) and types.isMatrix(B))))
        @compileError("matmul: at least one argument must be a matrix, the other must be a matrix or vector, got " ++ @typeName(A) ++ " and " ++ @typeName(B));

    if (comptime !types.isMatrix(A)) { // vector * matrix
        switch (comptime types.matrixType(B)) {
            .general => { // vector * general
                if (a.len != b.rows)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, b.cols);
                errdefer result.deinit(allocator);

                try blas.gemv(
                    types.orderOf(B),
                    .trans,
                    types.scast(i32, b.rows),
                    types.scast(i32, b.cols),
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
                if (a.len != b.size)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, b.size);
                errdefer result.deinit(allocator);

                try blas.symv(
                    types.orderOf(B).invert(),
                    types.uploOf(B).invert(),
                    types.scast(i32, b.size),
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
                if (a.len != b.size)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, b.size);
                errdefer result.deinit(allocator);

                try blas.hemv(
                    types.orderOf(B).invert(),
                    types.uploOf(B).invert(),
                    types.scast(i32, b.size),
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
                if (a.len != b.rows)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, b.cols);
                errdefer result.deinit(allocator);

                if (b.rows == b.cols) {
                    try blas.copy(
                        types.scast(i32, a.len),
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
                        types.scast(i32, b.rows),
                        b.data,
                        types.scast(i32, b.ld),
                        result.data,
                        result.inc,
                        ctx,
                    );
                } else {
                    const min_dim: u32 = int.min(b.rows, b.cols);

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
                        if (b.cols > min_dim) {
                            try blas.gemv(
                                types.orderOf(B),
                                .trans,
                                types.scast(i32, b.rows),
                                types.scast(i32, b.cols - min_dim),
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
                        if (b.rows > min_dim) {
                            try blas.gemv(
                                types.orderOf(B),
                                .trans,
                                types.scast(i32, b.rows - min_dim),
                                types.scast(i32, b.cols),
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
                        }
                        if (b.cols > min_dim) {
                            try blas.scal(
                                types.scast(i32, b.cols - min_dim),
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
                if (a.len != b.rows)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, b.cols);
                errdefer result.deinit(allocator);

                try blas.copy(
                    types.scast(i32, a.len),
                    a.data,
                    a.inc,
                    result.data,
                    result.inc,
                    ctx,
                );

                var i: u32 = 0;
                while (i < int.min(b.rows, b.cols)) : (i += 1) {
                    try ops.mul_( // result[i] *= b[i, i]
                        &result.data[i],
                        result.data[i],
                        b.data[i],
                        ctx,
                    );
                }

                while (i < b.cols) : (i += 1) {
                    result.data[i] = try constants.zero(C, ctx);
                }

                return result;
            },
            .banded => { // vector * banded
                if (a.len != b.rows)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, b.cols);
                errdefer result.deinit(allocator);

                try blas.gbmv(
                    types.orderOf(B),
                    .trans,
                    types.scast(i32, b.rows),
                    types.scast(i32, b.cols),
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
                if (a.len != b.size)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, b.size);
                errdefer result.deinit(allocator);

                // lapack.lagtm
                try defaultSlowVM(&result, a, b, ctx);

                return result;
            },
            .permutation => { // vector * permutation
                if (a.len != b.size)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, b.size);
                errdefer result.deinit(allocator);

                var i: u32 = 0;
                while (i < b.size) : (i += 1) {
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
        switch (comptime types.matrixType(A)) {
            .general => { // general * vector
                if (a.cols != b.len)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, a.rows);
                errdefer result.deinit(allocator);

                try blas.gemv(
                    types.orderOf(A),
                    .no_trans,
                    types.scast(i32, a.rows),
                    types.scast(i32, a.cols),
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
                if (a.size != b.len)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, a.size);
                errdefer result.deinit(allocator);

                try blas.symv(
                    types.orderOf(A),
                    types.uploOf(A),
                    types.scast(i32, a.size),
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
                if (a.size != b.len)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, a.size);
                errdefer result.deinit(allocator);

                try blas.hemv(
                    types.orderOf(A),
                    types.uploOf(A),
                    types.scast(i32, a.size),
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
                if (a.cols != b.len)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, a.rows);
                errdefer result.deinit(allocator);

                if (a.rows == a.cols) {
                    try blas.copy(
                        types.scast(i32, b.len),
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
                        types.scast(i32, a.rows),
                        a.data,
                        types.scast(i32, a.ld),
                        result.data,
                        result.inc,
                        ctx,
                    );
                } else {
                    const min_dim: u32 = int.min(a.rows, a.cols);

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
                        if (a.cols > min_dim) { // extra rows (filled)
                            try blas.gemv(
                                types.orderOf(A),
                                .no_trans,
                                types.scast(i32, a.rows),
                                types.scast(i32, a.cols - min_dim),
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
                        } else { // extra columns (zeroed)
                            try blas.scal(
                                types.scast(i32, a.rows - min_dim),
                                0,
                                result.data + min_dim * types.scast(u32, result.inc),
                                result.inc,
                                ctx,
                            );
                        }
                    } else {
                        if (a.rows > min_dim) { // extra rows (filled)
                            try blas.gemv(
                                types.orderOf(A),
                                .no_trans,
                                types.scast(i32, a.rows - min_dim),
                                types.scast(i32, a.cols),
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
                if (a.cols != b.len)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, a.rows);
                errdefer result.deinit(allocator);

                try blas.copy(
                    types.scast(i32, b.len),
                    b.data,
                    b.inc,
                    result.data,
                    result.inc,
                    ctx,
                );

                var i: u32 = 0;
                while (i < int.min(a.rows, a.cols)) : (i += 1) {
                    try ops.mul_( // result[i] *= a[i, i]
                        &result.data[i],
                        result.data[i],
                        a.data[i],
                        ctx,
                    );
                }

                while (i < a.rows) : (i += 1) {
                    result.data[i] = try constants.zero(C, ctx);
                }

                return result;
            },
            .banded => { // banded * vector
                if (a.cols != b.len)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, a.rows);
                errdefer result.deinit(allocator);

                try blas.gbmv(
                    types.orderOf(A),
                    .no_trans,
                    types.scast(i32, a.rows),
                    types.scast(i32, a.cols),
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
                if (a.size != b.len)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, a.size);
                errdefer result.deinit(allocator);

                // lapack.lagtm
                try defaultSlowMV(&result, a, b, ctx);

                return result;
            },
            .permutation => { // permutation * vector
                if (a.size != b.len)
                    return linalg.Error.DimensionMismatch;

                var result: vector.Vector(C) = try .init(allocator, a.size);
                errdefer result.deinit(allocator);

                var i: u32 = 0;
                while (i < a.size) : (i += 1) {
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
        switch (comptime types.matrixType(A)) {
            .general => switch (comptime types.matrixType(B)) {
                .general => { // general * general
                    if (a.cols != b.rows)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.rows, b.cols);
                    errdefer result.deinit(allocator);

                    try blas.gemm(
                        types.orderOf(A),
                        .no_trans,
                        if (comptime types.orderOf(A) == types.orderOf(B)) .no_trans else .trans,
                        types.scast(i32, a.rows),
                        types.scast(i32, b.cols),
                        types.scast(i32, a.cols),
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
                    if (a.cols != b.size)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.rows, b.size);
                    errdefer result.deinit(allocator);

                    try blas.symm(
                        types.orderOf(A),
                        .right,
                        if (comptime types.orderOf(A) == types.orderOf(B)) types.uploOf(B) else types.uploOf(B).invert(),
                        types.scast(i32, a.rows),
                        types.scast(i32, b.size),
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
                    if (a.cols != b.size)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.rows, b.size);
                    errdefer result.deinit(allocator);

                    try blas.hemm(
                        types.orderOf(A),
                        .right,
                        if (comptime types.orderOf(A) == types.orderOf(B)) types.uploOf(B) else types.uploOf(B).invert(),
                        types.scast(i32, a.rows),
                        types.scast(i32, b.size),
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
                    if (a.cols != b.rows)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.rows, b.cols);
                    errdefer result.deinit(allocator);

                    if (b.rows == b.cols) {
                        if (comptime types.orderOf(A) == .col_major) {
                            var j: u32 = 0;
                            while (j < b.cols) : (j += 1) {
                                try blas.copy(
                                    types.scast(i32, a.rows),
                                    a.data + j * a.ld,
                                    1,
                                    result.data + j * result.ld,
                                    1,
                                    ctx,
                                );
                            }
                        } else {
                            var i: u32 = 0;
                            while (i < a.rows) : (i += 1) {
                                try blas.copy(
                                    types.scast(i32, b.cols),
                                    a.data + i * a.ld,
                                    1,
                                    result.data + i * result.ld,
                                    1,
                                    ctx,
                                );
                            }
                        }

                        try blas.trmm(
                            types.orderOf(A),
                            .right,
                            if (comptime types.orderOf(A) == types.orderOf(B)) types.uploOf(B) else types.uploOf(B).invert(),
                            if (comptime types.orderOf(A) == types.orderOf(B)) .no_trans else .trans,
                            types.diagOf(B),
                            types.scast(i32, a.rows),
                            types.scast(i32, b.cols),
                            1,
                            b.data,
                            types.scast(i32, b.ld),
                            result.data,
                            types.scast(i32, result.ld),
                            ctx,
                        );
                    } else {
                        const min_dim: u32 = int.min(b.rows, b.cols);

                        if (comptime types.orderOf(A) == .col_major) {
                            var j: u32 = 0;
                            while (j < min_dim) : (j += 1) {
                                try blas.copy(
                                    types.scast(i32, a.rows),
                                    a.data + j * a.ld,
                                    1,
                                    result.data + j * result.ld,
                                    1,
                                    ctx,
                                );
                            }
                        } else {
                            var i: u32 = 0;
                            while (i < a.rows) : (i += 1) {
                                try blas.copy(
                                    types.scast(i32, min_dim),
                                    a.data + i * a.ld,
                                    types.scast(i32, a.ld),
                                    result.data + i * result.ld,
                                    types.scast(i32, result.ld),
                                    ctx,
                                );
                            }
                        }

                        try blas.trmm(
                            types.orderOf(A),
                            .right,
                            if (comptime types.orderOf(A) == types.orderOf(B)) types.uploOf(B) else types.uploOf(B).invert(),
                            if (comptime types.orderOf(A) == types.orderOf(B)) .no_trans else .trans,
                            types.diagOf(B),
                            types.scast(i32, a.rows),
                            types.scast(i32, min_dim),
                            1,
                            b.data,
                            types.scast(i32, b.ld),
                            result.data,
                            types.scast(i32, result.ld),
                            ctx,
                        );

                        if (comptime types.uploOf(B) == .upper) {
                            if (b.cols > min_dim) {
                                try blas.gemm(
                                    types.orderOf(A),
                                    .no_trans,
                                    if (comptime types.orderOf(A) == types.orderOf(B)) .no_trans else .trans,
                                    types.scast(i32, a.rows),
                                    types.scast(i32, b.cols - min_dim),
                                    types.scast(i32, b.rows),
                                    1,
                                    a.data,
                                    types.scast(i32, a.ld),
                                    b.data + if (comptime types.orderOf(B) == .col_major)
                                        min_dim * b.ld
                                    else
                                        min_dim,
                                    types.scast(i32, b.ld),
                                    0,
                                    result.data + if (comptime types.orderOf(A) == .col_major)
                                        min_dim * result.ld
                                    else
                                        min_dim,
                                    types.scast(i32, result.ld),
                                    ctx,
                                );
                            }
                        } else {
                            if (b.rows > min_dim) {
                                try blas.gemm(
                                    types.orderOf(A),
                                    .no_trans,
                                    if (comptime types.orderOf(A) == types.orderOf(B)) .no_trans else .trans,
                                    types.scast(i32, a.rows),
                                    types.scast(i32, b.cols),
                                    types.scast(i32, b.rows - min_dim),
                                    1,
                                    a.data + if (comptime types.orderOf(A) == .col_major)
                                        min_dim
                                    else
                                        min_dim * a.ld,
                                    types.scast(i32, a.ld),
                                    b.data + if (comptime types.orderOf(B) == .col_major)
                                        min_dim
                                    else
                                        min_dim * b.ld,
                                    types.scast(i32, b.ld),
                                    1,
                                    result.data,
                                    types.scast(i32, result.ld),
                                    ctx,
                                );
                            }
                            if (b.cols > min_dim) {
                                var j: u32 = min_dim;
                                while (j < b.cols) : (j += 1) {
                                    try blas.scal(
                                        types.scast(i32, a.rows),
                                        0,
                                        result.data + if (comptime types.orderOf(A) == .col_major)
                                            j * result.ld
                                        else
                                            j,
                                        if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, result.ld),
                                        ctx,
                                    );
                                }
                            }
                        }
                    }

                    return result;
                },
                .diagonal => { // general * diagonal
                    if (a.cols != b.rows)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .full(allocator, a.rows, b.cols, 0, ctx);
                    errdefer result.deinit(allocator);

                    // lapack.lascl2

                    var j: u32 = 0;
                    while (j < int.min(a.cols, b.cols)) : (j += 1) {
                        try blas.axpy(
                            types.scast(i32, a.rows),
                            b.data[j],
                            a.data + if (comptime types.orderOf(A) == .col_major) j * a.ld else j,
                            if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, a.ld),
                            result.data + if (comptime types.orderOf(A) == .col_major) j * result.ld else j,
                            if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, result.ld),
                            ctx,
                        );
                    }

                    while (j < b.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < a.rows) : (i += 1) {
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
                    if (a.cols != b.rows)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.rows, b.cols);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .tridiagonal => { // general * tridiagonal
                    if (a.cols != b.size)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.rows, b.cols);
                    errdefer result.deinit(allocator);

                    // lapack.lagtm
                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .permutation => return linalg.Error.NotImplemented, // lapack.lapmt
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .symmetric => switch (comptime types.matrixType(B)) {
                .general => {
                    if (a.size != b.rows)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.size, b.cols);
                    errdefer result.deinit(allocator);

                    if (comptime types.orderOf(A) == types.orderOf(B)) {
                        try blas.symm(
                            types.orderOf(B),
                            .left,
                            types.uploOf(A),
                            types.scast(i32, a.size),
                            types.scast(i32, b.cols),
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
                        while (j < b.cols) : (j += 1) {
                            try blas.symv(
                                types.orderOf(A),
                                types.uploOf(A),
                                types.scast(i32, a.size),
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
                    if (a.size != b.size)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.size, a.size);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .hermitian => { // symmetric * hermitian
                    if (a.size != b.size)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.size, a.size);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .triangular => { // symmetric * triangular
                    if (a.size != b.rows)
                        return linalg.Error.DimensionMismatch;

                    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, a.size, b.cols);
                    errdefer result.deinit(allocator);

                    try defaultSlowMM(&result, a, b, ctx);

                    return result;
                },
                .diagonal => return linalg.Error.NotImplemented,
                .banded => return linalg.Error.NotImplemented,
                .tridiagonal => return linalg.Error.NotImplemented,
                .permutation => return linalg.Error.NotImplemented,
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .hermitian => switch (comptime types.matrixType(B)) {
                .general => return linalg.Error.NotImplemented,
                .symmetric => return linalg.Error.NotImplemented,
                .hermitian => return linalg.Error.NotImplemented,
                .triangular => return linalg.Error.NotImplemented,
                .diagonal => return linalg.Error.NotImplemented,
                .banded => return linalg.Error.NotImplemented,
                .tridiagonal => return linalg.Error.NotImplemented,
                .permutation => return linalg.Error.NotImplemented,
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .triangular => switch (comptime types.matrixType(B)) {
                .general => return linalg.Error.NotImplemented,
                .symmetric => return linalg.Error.NotImplemented,
                .hermitian => return linalg.Error.NotImplemented,
                .triangular => return linalg.Error.NotImplemented,
                .diagonal => return linalg.Error.NotImplemented,
                .banded => return linalg.Error.NotImplemented,
                .tridiagonal => return linalg.Error.NotImplemented,
                .permutation => return linalg.Error.NotImplemented,
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .diagonal => switch (comptime types.matrixType(B)) {
                .general => return linalg.Error.NotImplemented,
                .symmetric => return linalg.Error.NotImplemented,
                .hermitian => return linalg.Error.NotImplemented,
                .triangular => return linalg.Error.NotImplemented,
                .diagonal => return linalg.Error.NotImplemented,
                .banded => return linalg.Error.NotImplemented,
                .tridiagonal => return linalg.Error.NotImplemented,
                .permutation => return linalg.Error.NotImplemented,
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .banded => switch (comptime types.matrixType(B)) {
                .general => return linalg.Error.NotImplemented,
                .symmetric => return linalg.Error.NotImplemented,
                .hermitian => return linalg.Error.NotImplemented,
                .triangular => return linalg.Error.NotImplemented,
                .diagonal => return linalg.Error.NotImplemented,
                .banded => return linalg.Error.NotImplemented,
                .tridiagonal => return linalg.Error.NotImplemented,
                .permutation => return linalg.Error.NotImplemented,
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .tridiagonal => switch (comptime types.matrixType(B)) {
                .general => return linalg.Error.NotImplemented,
                .symmetric => return linalg.Error.NotImplemented,
                .hermitian => return linalg.Error.NotImplemented,
                .triangular => return linalg.Error.NotImplemented,
                .diagonal => return linalg.Error.NotImplemented,
                .banded => return linalg.Error.NotImplemented,
                .tridiagonal => return linalg.Error.NotImplemented,
                .permutation => return linalg.Error.NotImplemented,
                .sparse => @compileError("apply2 not implemented for sparse matrices yet"),
                .numeric => unreachable,
            },
            .permutation => switch (comptime types.matrixType(B)) {
                .general => return linalg.Error.NotImplemented,
                .symmetric => return linalg.Error.NotImplemented,
                .hermitian => return linalg.Error.NotImplemented,
                .triangular => return linalg.Error.NotImplemented,
                .diagonal => return linalg.Error.NotImplemented,
                .banded => return linalg.Error.NotImplemented,
                .tridiagonal => return linalg.Error.NotImplemented,
                .permutation => return linalg.Error.NotImplemented,
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

    if (x.len != n)
        return linalg.Error.DimensionMismatch;

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

    const m: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.rows;
    const n: u32 = if (comptime types.isSymmetricMatrix(B) or types.isHermitianMatrix(B) or types.isTridiagonalMatrix(B)) b.size else b.cols;
    const k: u32 = if (comptime types.isSymmetricMatrix(A) or types.isHermitianMatrix(A) or types.isTridiagonalMatrix(A)) a.size else a.cols;

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
