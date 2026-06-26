const std = @import("std");

const meta = @import("../meta.zig");

const int = @import("../int.zig");

const numeric = @import("../numeric.zig");
const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");

const linalg = @import("../linalg.zig");

const vecutils = @import("../vector/utils.zig");
const matutils = @import("../matrix/utils.zig");

pub fn Matmul(comptime X: type, comptime Y: type) type {
    comptime if ((!meta.isVector(X) and !meta.isMatrix(X)) or (!meta.isVector(Y) and !meta.isMatrix(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)))
        @compileError("zsl.linalg.Matmul: at least one of X or Y must be a matrix type, the other must be a matrix or a vector type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    comptime if (meta.isBuilderMatrix(X) or meta.isBuilderMatrix(Y))
        @compileError("zsl.linalg.Matmul: builder matrix types are not allowed, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    const R = numeric.Fma(meta.Numeric(X), meta.Numeric(Y), numeric.Mul(meta.Numeric(X), meta.Numeric(Y)));

    comptime if ((meta.isStaticVector(X) or meta.isStaticMatrix(X)) and (meta.isStaticVector(Y) or meta.isStaticMatrix(Y))) {
        if ((if (meta.isVector(X)) X.len else X.cols) != (if (meta.isVector(Y)) Y.len else Y.rows))
            @compileError("zsl.linalg.Matmul: static types X and Y must have compatible inner dimensions (x cols == y rows), got\n\tX = " ++
                @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");
    };

    switch (comptime meta.domain(X)) {
        .matrix => switch (comptime meta.domain(Y)) {
            .matrix => {
                const x_static = comptime meta.matrixStorage(X) == .static;
                const y_static = comptime meta.matrixStorage(Y) == .static;
                const x_square = comptime meta.isSquareMatrix(X);
                const y_square = comptime meta.isSquareMatrix(Y);
                const o_static = comptime (x_static or (x_square and y_static)) and (y_static or (y_square and x_static));
                const o_rows = comptime if (x_static) X.rows else if (y_static) Y.rows else 0;
                const o_cols = comptime if (y_static) Y.cols else if (x_static) X.cols else 0;

                switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense => return if (comptime o_static)
                        matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(X).?)
                    else
                        matrix.general.Dense(R, meta.layoutOf(X).?),
                    .general_sparse, .symmetric_sparse, .hermitian_sparse => return if (comptime o_static)
                        matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(X).?)
                    else if (comptime meta.matrixStorage(Y) == .sparse or meta.matrixKind(Y) == .diagonal or meta.matrixKind(Y) == .permutation)
                        matrix.general.Sparse(R, meta.layoutOf(X).?)
                    else
                        matrix.general.Dense(R, meta.layoutOf(X).?),
                    .triangular_static, .triangular_dense => switch (comptime meta.matrixKind(Y)) {
                        .triangular => switch (comptime meta.uploOf(X).? == meta.uploOf(Y).?) {
                            true => return if (comptime o_static)
                                matrix.triangular.Static(o_rows, o_cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?)
                            else
                                matrix.triangular.Dense(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                            false => return if (comptime o_static)
                                matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(X).?)
                            else
                                matrix.general.Dense(R, meta.layoutOf(X).?),
                        },
                        .diagonal => return if (comptime o_static)
                            matrix.triangular.Static(o_rows, o_cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?)
                        else
                            matrix.triangular.Dense(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                        else => return if (comptime o_static)
                            matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(X).?)
                        else
                            matrix.general.Dense(R, meta.layoutOf(X).?),
                    },
                    .triangular_sparse => switch (comptime meta.matrixKind(Y)) {
                        .triangular => switch (comptime meta.uploOf(X).? == meta.uploOf(Y).?) {
                            true => return if (comptime o_static)
                                matrix.triangular.Static(o_rows, o_cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?)
                            else if (comptime meta.matrixStorage(Y) == .sparse)
                                matrix.triangular.Sparse(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?)
                            else
                                matrix.triangular.Dense(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                            false => return if (comptime o_static)
                                matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(X).?)
                            else if (comptime meta.matrixStorage(Y) == .sparse)
                                matrix.general.Sparse(R, meta.layoutOf(X).?)
                            else
                                matrix.general.Dense(R, meta.layoutOf(X).?),
                        },
                        .diagonal => return if (comptime o_static)
                            matrix.triangular.Static(o_rows, o_cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?)
                        else
                            matrix.triangular.Sparse(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                        else => return if (comptime o_static)
                            matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(X).?)
                        else if (comptime meta.matrixStorage(Y) == .sparse or meta.matrixKind(Y) == .diagonal or meta.matrixKind(Y) == .permutation)
                            matrix.general.Sparse(R, meta.layoutOf(X).?)
                        else
                            matrix.general.Dense(R, meta.layoutOf(X).?),
                    },
                    .diagonal_static => switch (comptime meta.matrixKind(Y)) {
                        .diagonal => return if (comptime o_static)
                            matrix.diagonal.Static(o_rows, o_cols, R)
                        else
                            matrix.diagonal.Sparse(R),
                        .triangular => return if (comptime o_static)
                            matrix.triangular.Static(o_rows, o_cols, R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?)
                        else if (comptime meta.matrixStorage(Y) == .sparse)
                            matrix.triangular.Sparse(R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?)
                        else
                            matrix.triangular.Dense(R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
                        else => return if (comptime o_static)
                            matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(Y) orelse matrix.Layout.default)
                        else if (comptime meta.matrixStorage(Y) == .sparse)
                            matrix.general.Sparse(R, meta.layoutOf(Y) orelse matrix.Layout.default)
                        else
                            matrix.general.Dense(R, meta.layoutOf(Y) orelse matrix.Layout.default),
                    },
                    .diagonal_sparse => switch (comptime meta.matrixKind(Y)) {
                        .diagonal => return if (comptime o_static)
                            matrix.diagonal.Static(o_rows, o_cols, R)
                        else
                            matrix.diagonal.Sparse(R),
                        .triangular => return if (comptime o_static)
                            matrix.triangular.Static(o_rows, o_cols, R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?)
                        else if (comptime meta.matrixStorage(Y) == .sparse)
                            matrix.triangular.Sparse(R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?)
                        else
                            matrix.triangular.Dense(R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
                        else => return if (comptime o_static)
                            matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(Y) orelse matrix.Layout.default)
                        else if (comptime meta.matrixStorage(Y) == .sparse or meta.matrixKind(Y) == .permutation)
                            matrix.general.Sparse(R, meta.layoutOf(Y) orelse matrix.Layout.default)
                        else
                            matrix.general.Dense(R, meta.layoutOf(Y) orelse matrix.Layout.default),
                    },
                    .permutation_static => switch (comptime meta.matrixKind(Y)) {
                        .permutation => return if (comptime o_static)
                            matrix.permutation.Static(o_rows, R, X.direction)
                        else
                            matrix.permutation.Sparse(R, X.direction),
                        else => return if (comptime o_static)
                            matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(Y) orelse matrix.Layout.default)
                        else if (comptime meta.matrixStorage(Y) == .sparse)
                            matrix.general.Sparse(R, meta.layoutOf(Y) orelse matrix.Layout.default)
                        else
                            matrix.general.Dense(R, meta.layoutOf(Y) orelse matrix.Layout.default),
                    },
                    .permutation_sparse => switch (comptime meta.matrixKind(Y)) {
                        .permutation => return if (comptime o_static)
                            matrix.permutation.Static(o_rows, R, X.direction)
                        else
                            matrix.permutation.Sparse(R, X.direction),
                        else => return if (comptime o_static)
                            matrix.general.Static(o_rows, o_cols, R, meta.layoutOf(Y) orelse matrix.Layout.default)
                        else if (comptime meta.matrixStorage(Y) == .sparse or meta.matrixKind(Y) == .diagonal)
                            matrix.general.Sparse(R, meta.layoutOf(Y) orelse matrix.Layout.default)
                        else
                            matrix.general.Dense(R, meta.layoutOf(Y) orelse matrix.Layout.default),
                    },
                    else => unreachable,
                }
            },
            .vector => {
                const x_static = comptime meta.matrixStorage(X) == .static;
                const y_static = comptime meta.isStaticVector(Y);
                const x_square = comptime meta.isSquareMatrix(X);
                const o_static = comptime x_static or (x_square and y_static);
                const o_len = comptime if (x_static) X.rows else if (y_static) Y.len else 0;

                return if (comptime o_static)
                    vector.Static(o_len, R)
                else if (comptime meta.matrixStorage(X) == .sparse and meta.vectorType(Y) == .sparse)
                    vector.Sparse(R)
                else
                    vector.Dense(R);
            },
            else => unreachable,
        },
        .vector => switch (comptime meta.domain(Y)) {
            .matrix => {
                const x_static = comptime meta.isStaticVector(X);
                const y_static = comptime meta.matrixStorage(Y) == .static;
                const y_square = comptime meta.isSquareMatrix(Y);
                const o_static = comptime y_static or (y_square and x_static);
                const o_len = comptime if (y_static) Y.cols else if (x_static) X.len else 0;

                return if (comptime o_static)
                    vector.Static(o_len, R)
                else if (comptime meta.vectorType(X) == .sparse and meta.matrixStorage(Y) == .sparse)
                    vector.Sparse(R)
                else
                    vector.Dense(R);
            },
            else => unreachable,
        },
        else => unreachable,
    }
}

/// Performs matrix multiplication between two matrices, or between a matrix and
/// a vector.
///
/// For matrix outputs, the result inherits its memory layout from the inputs,
/// i.e., if the input layouts mismatch, the left operand (`x`) strictly
/// dictates the output layout, unless it provides no layout information. For
/// more control over layouts, use `linalg.matmulInto`.
///
/// This function is intended for when the result's dimension is known at
/// compile time. For two static inputs, dimension checks are performed at
/// compile time, for any other combination, dimension checks are performed at
/// runtime throught `std.debug.assert`.
///
/// ## Signature
/// ```zig
/// linalg.matmul(x: X, y: Y) linalg.Matmul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `linalg.Matmul(@TypeOf(x), @TypeOf(y))`: The result of the operation.
pub fn matmul(x: anytype, y: anytype) linalg.Matmul(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = linalg.Matmul(X, Y);

    if (comptime meta.isDenseVector(R) or meta.isSparseVector(R) or
        meta.isDenseMatrix(R) or meta.isSparseMatrix(R))
        @compileError("zsl.linalg.matmul: the result cannot be a heap-allocated type, i.e., both inputs must be static, or any vector and a static matrix, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tresult: " ++ @typeName(R) ++ "\nFor these inputs use zsl.linalg.matmulAlloc instead.");

    const x_cols = if (comptime meta.isVector(X))
        (if (comptime meta.isStaticVector(X)) X.len else x.len)
    else
        (if (comptime meta.isStaticMatrix(X)) X.cols else x.cols);
    const y_rows = if (comptime meta.isVector(Y))
        (if (comptime meta.isStaticVector(Y)) Y.len else y.len)
    else
        (if (comptime meta.isStaticMatrix(Y)) Y.rows else y.rows);

    if (comptime !((meta.isVector(X) and meta.isStaticVector(X)) or meta.isStaticMatrix(X)) or
        !((meta.isVector(Y) and meta.isStaticVector(Y)) or meta.isStaticMatrix(Y)))
        std.debug.assert(x_cols == y_rows);

    var result = R.init;

    linalg.matmulIntoUnchecked(&result, x, y);

    return result;
}

/// Performs matrix multiplication between two matrices, or between a matrix and
/// a vector, dynamically allocating memory for the result.
///
/// For matrix outputs, the result inherits its memory layout from the inputs,
/// i.e., if the input layouts mismatch, the left operand (`x`) strictly
/// dictates the output layout, unless it provides no layout information. For
/// more control over layouts, use `linalg.matmulInto`.
///
/// This function is intended for when the result's dimension is known at
/// runtime.
///
/// For two sparse inputs, the operation is only applied to the indices where at
/// least one of the matrices has a non-zero element.
///
/// ## Signature
/// ```zig
/// linalg.matmulAlloc(allocator: std.mem.Allocator, x: X, y: Y) !linalg.Matmul(X, Y)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `linalg.Matmul(@TypeOf(x), @TypeOf(y))`: The result of the operation.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `linalg.Error.DimensionMismatch`: If the inputs do not have compatible
///   dimensions.
pub fn matmulAlloc(allocator: std.mem.Allocator, x: anytype, y: anytype) !linalg.Matmul(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = linalg.Matmul(X, Y);

    const x_rows = if (comptime meta.isVector(X))
        1
    else
        (if (comptime meta.isStaticMatrix(X)) X.rows else x.rows);
    const x_cols = if (comptime meta.isVector(X))
        (if (comptime meta.isStaticVector(X)) X.len else x.len)
    else
        (if (comptime meta.isStaticMatrix(X)) X.cols else x.cols);
    const y_rows = if (comptime meta.isVector(Y))
        (if (comptime meta.isStaticVector(Y)) Y.len else y.len)
    else
        (if (comptime meta.isStaticMatrix(Y)) Y.rows else y.rows);
    const y_cols = if (comptime meta.isVector(Y))
        1
    else
        (if (comptime meta.isStaticMatrix(Y)) Y.cols else y.cols);

    if (comptime !((meta.isVector(X) and meta.isStaticVector(X)) or meta.isStaticMatrix(X)) or
        !((meta.isVector(Y) and meta.isStaticVector(Y)) or meta.isStaticMatrix(Y)))
    {
        if (x_cols != y_rows)
            return matrix.Error.DimensionMismatch;
    }

    var result = switch (comptime meta.domain(R)) {
        .matrix => switch (comptime meta.matrixStorage(R)) {
            .static => R.init,
            .dense => switch (comptime meta.matrixKind((R))) {
                .general, .triangular => try R.init(allocator, x_rows, y_cols),
                else => unreachable,
            },
            .sparse => switch (comptime meta.matrixKind(R)) {
                .general, .triangular => blk: {
                    const work = try allocator.alloc(usize, int.max(x_rows, y_cols));
                    defer allocator.free(work);

                    break :blk try R.init(allocator, x_rows, y_cols, matmulNNZ(R, x_rows, y_cols, x_cols, x, y, work));
                },
                .diagonal => try R.init(allocator, x_rows, y_cols),
                .permutation => try R.init(allocator, x_rows),
                else => unreachable,
            },
            .numeric => unreachable,
        },
        .vector => switch (comptime meta.vectorType(R)) {
            .static => R.init,
            .dense => try R.init(allocator, if (comptime meta.isVector(X)) y_cols else x_rows),
            .sparse => blk: {
                const work = try allocator.alloc(usize, int.max(x_rows, y_cols));
                defer allocator.free(work);

                break :blk try R.init(allocator, if (comptime meta.isVector(X)) y_cols else x_rows, matmulNNZ(R, x_rows, y_cols, x_cols, x, y, work));
            },
            else => unreachable,
        },
        else => unreachable,
    };

    linalg.matmulIntoUnchecked(&result, x, y);

    return result;
}

/// Performs matrix multiplication between two matrices, or between a matrix and
/// a vector, storing the result in an output matrix or vector.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// not permitted and might yield incorrect results.
///
/// For three static inputs, or a static output, a static matrix and any vector,
/// the function cannot return an error.
///
/// For sparse outputs, the operation is only applied to the indices where at
/// least one of the inputs has a non-zero element.
///
/// ## Signature
/// ```zig
/// linalg.matmulInto(*O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left input operand.
/// * `y` (`anytype`): The right input operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.Error.DimensionMismatch`: If the matrices do not have the same
///   dimensions.
pub fn matmulInto(o: anytype, x: anytype, y: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or (!meta.isVector(meta.Child(O)) and !meta.isMatrix(meta.Child(O))) or
        (!meta.isVector(X) and !meta.isMatrix(X)) or (!meta.isVector(Y) and !meta.isMatrix(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)))
        @compileError("zsl.linalg.matmulInto: o must be a mutable one-itme pointer to a matrix or vector, at least one of x or y must be a matrix, the other must be a matrix or a vector, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    O = meta.Child(O);

    const o_rows = if (comptime meta.isVector(O))
        (if (comptime meta.isVector(X)) 1 else (if (comptime meta.isStaticVector(O)) O.len else o.len))
    else
        (if (comptime meta.isStaticMatrix(O)) O.rows else o.rows);
    const o_cols = if (comptime meta.isVector(O))
        (if (comptime meta.isVector(Y)) 1 else (if (comptime meta.isStaticVector(O)) O.len else o.len))
    else
        (if (comptime meta.isStaticMatrix(O)) O.cols else o.cols);
    const x_rows = if (comptime meta.isVector(X))
        1
    else
        (if (comptime meta.isStaticMatrix(X)) X.rows else x.rows);
    const x_cols = if (comptime meta.isVector(X))
        (if (comptime meta.isStaticVector(X)) X.len else x.len)
    else
        (if (comptime meta.isStaticMatrix(X)) X.cols else x.cols);
    const y_rows = if (comptime meta.isVector(Y))
        (if (comptime meta.isStaticVector(Y)) Y.len else y.len)
    else
        (if (comptime meta.isStaticMatrix(Y)) Y.rows else y.rows);
    const y_cols = if (comptime meta.isVector(Y))
        1
    else
        (if (comptime meta.isStaticMatrix(Y)) Y.cols else y.cols);

    if (comptime ((meta.isVector(X) and meta.isStaticVector(X)) or meta.isStaticMatrix(X)) and
        ((meta.isVector(Y) and meta.isStaticVector(Y)) or meta.isStaticMatrix(Y)))
    {
        if (comptime x_cols != y_rows)
            @compileError("zsl.linalg.matmulInto: inner dimensions mismatch (x cols != y rows), got\n\tx: " ++
                @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");
    } else {
        if (x_cols != y_rows)
            return linalg.Error.DimensionMismatch;
    }

    if (comptime ((meta.isVector(O) and (meta.isVector(X) or meta.isStaticVector(O))) or meta.isStaticMatrix(O)) and
        (meta.isVector(X) or meta.isStaticMatrix(X)))
    {
        if (comptime o_rows != x_rows)
            @compileError("zsl.linalg.matmulInto: output rows mismatch (o rows != x rows), got\n\to: *" ++
                @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");
    } else {
        if (o_rows != x_rows)
            return linalg.Error.DimensionMismatch;
    }

    if (comptime ((meta.isVector(O) and (meta.isVector(Y) or meta.isStaticVector(O))) or meta.isStaticMatrix(O)) and
        (meta.isVector(Y) or meta.isStaticMatrix(Y)))
    {
        if (comptime o_cols != y_cols)
            @compileError("zsl.linalg.matmulInto: output cols mismatch (o cols != y cols), got\n\to: *" ++
                @typeName(O) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");
    } else {
        if (o_cols != y_cols)
            return linalg.Error.DimensionMismatch;
    }

    if (comptime meta.isSparseMatrix(O) and !meta.isDiagonalMatrix(O) and !meta.isPermutationMatrix(O)) {
        if (comptime (meta.isStaticMatrix(X) and !meta.isDiagonalMatrix(X) and !meta.isPermutationMatrix(X)) or meta.isDenseMatrix(X) or
            (meta.isStaticMatrix(Y) and !meta.isDiagonalMatrix(Y) and !meta.isPermutationMatrix(Y)) or meta.isDenseMatrix(Y))
            @compileError("zsl.linalg.matmulInto: o cannot point to a sparse matrix if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

        const nnz = matmulNNZ(O, o_rows, o_cols, x_cols, x, y, null);

        if (comptime meta.isBuilderMatrix(O)) {
            if (o._dlen < nnz or o._rlen < nnz or o._clen < nnz)
                return matrix.Error.InsufficientSpace;
        } else {
            if (o._dlen < nnz or o._ilen < nnz)
                return matrix.Error.InsufficientSpace;
        }
    }

    if (comptime meta.isSparseVector(O)) {
        if (comptime meta.isStaticVector(X) or meta.isDenseVector(X) or (meta.isStaticMatrix(X) and !meta.isDiagonalMatrix(X) and !meta.isPermutationMatrix(X)) or meta.isDenseMatrix(X) or
            meta.isStaticVector(Y) or meta.isDenseVector(Y) or (meta.isStaticMatrix(Y) and !meta.isDiagonalMatrix(Y) and !meta.isPermutationMatrix(Y)) or meta.isDenseMatrix(Y))
            @compileError("zsl.linalg.matmulInto: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

        const nnz = matmulNNZ(O, o_rows, o_cols, x, y);

        if (o._dlen < nnz or o._ilen < nnz)
            return vector.Error.InsufficientSpace;
    }

    return matmulIntoUnchecked(o, x, y);
}

/// Performs matrix multiplication between two matrices, or between a matrix and
/// a vector, storing the result in an output matrix or vector, without
/// performing dimension checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// not permitted and might yield incorrect results.
///
/// For sparse outputs, the operation is only applied to the indices where at
/// least one of the inputs has a non-zero element.
///
/// ## Signature
/// ```zig
/// linalg.matmulIntoUnchecked(*O, x: X, y: Y) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left input operand.
/// * `y` (`anytype`): The right input operand.
///
/// ## Returns
/// `void`
pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or (!meta.isVector(meta.Child(O)) and !meta.isMatrix(meta.Child(O))) or
        (!meta.isVector(X) and !meta.isMatrix(X)) or (!meta.isVector(Y) and !meta.isMatrix(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)))
        @compileError("zsl.linalg.matmulInto: o must be a mutable one-itme pointer to a matrix or vector, at least one of x or y must be a matrix, the other must be a matrix or a vector, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    O = meta.Child(O);

    switch (comptime meta.domain(O)) {
        .matrix => switch (comptime meta.matrixType(O)) {
            .general_static => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgensta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .general_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .general_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matgenspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .symmetric_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymsta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .symmetric_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .symmetric_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matsymspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .hermitian_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mathersta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .hermitian_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .hermitian_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matherspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrista_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattriden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_mattrispa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiasta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matdiaspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matpersta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgensta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgensta_matperspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .general_dense => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgensta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .general_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .general_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matgenspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .symmetric_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymsta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .symmetric_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .symmetric_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matsymspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .hermitian_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mathersta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .hermitian_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .hermitian_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matherspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrista_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattriden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_mattrispa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiasta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matdiaspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matpersta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenden_matperspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .general_sparse => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense, .triangular_static, .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .general_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matgenspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matgenspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matgenspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matgenspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matgenspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matgenspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matgenspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matgenspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .symmetric_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matsymspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matsymspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matsymspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matsymspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matsymspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matsymspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matsymspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matsymspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .hermitian_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matherspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matherspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matherspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matherspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matherspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matherspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matherspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matherspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_mattrispa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_mattrispa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_mattrispa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_mattrispa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_mattrispa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_mattrispa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_mattrispa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_mattrispa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiasta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiasta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiasta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiasta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiasta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiasta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiaspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiaspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiaspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiaspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiaspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matdiaspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matpersta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matpersta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matpersta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matpersta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matpersta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matpersta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matpersta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matpersta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matperspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matperspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matperspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matperspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matperspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matperspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matperspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matgenspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matgenspa_matperspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .symmetric_static, .symmetric_dense, .symmetric_sparse => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: symmetric output is mathematically unsafe for general multiplication; use a specialized kernel or a general output\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .hermitian_static, .hermitian_dense, .hermitian_sparse => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: hermitian output is mathematically unsafe for general multiplication; use a specialized kernel or a general output\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => switch (comptime meta.matrixType(Y)) {
                        .builder_sparse => unreachable,
                        else => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    },
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .triangular_static => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .permutation_static, .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrista_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrista_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrista_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrista_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrista_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattriden_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattriden_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattriden_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattriden_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattriden_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrispa_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrispa_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrispa_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrispa_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_mattrispa_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiasta_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiasta_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiasta_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiaspa_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiaspa_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiaspa_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrista_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrista_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .triangular_dense => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .permutation_static, .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrista_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrista_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrista_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrista_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrista_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattriden_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattriden_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattriden_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattriden_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattriden_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrispa_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrispa_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrispa_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrispa_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_mattrispa_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiasta_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiasta_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiasta_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiaspa_mattrista.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_dense => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiaspa_mattriden.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiaspa_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattriden_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattriden_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .triangular_sparse => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .permutation_static, .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_static, .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static, .triangular_dense, .triangular_sparse, .diagonal_static, .diagonal_sparse => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static, .triangular_dense => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_mattrispa_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_mattrispa_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_mattrispa_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static, .triangular_dense => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_matdiasta_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .triangular_static, .triangular_dense => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                            .triangular_sparse => {
                                comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_matdiaspa_mattrispa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_static => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .diagonal_sparse => {
                                comptime if (meta.diagOf(O) == .unit)
                                    @compileError("zsl.linalg.matmulInto: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                                        @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                                return @import("matmul/mattrispa_slow.zig").matmulIntoUnchecked(o, x, y); // return @import("matmul/mattrispa_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y);
                            },
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: triangular output requires both inputs to be triangular or diagonal\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .diagonal_static => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: diagonal output requires both inputs to be diagonal matrices\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .diagonal_static => return @import("matmul/matdiasta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matdiasta_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matdiasta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matdiasta_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: diagonal output requires both inputs to be diagonal matrices\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .diagonal_static => return @import("matmul/matdiasta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matdiasta_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matdiasta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matdiasta_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: diagonal output requires both inputs to be diagonal matrices\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .diagonal_sparse => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: diagonal output requires both inputs to be diagonal matrices\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .diagonal_static => return @import("matmul/matdiaspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matdiaspa_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matdiaspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matdiaspa_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: diagonal output requires both inputs to be diagonal matrices\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .diagonal_static => return @import("matmul/matdiaspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matdiaspa_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matdiaspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matdiaspa_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: diagonal output requires both inputs to be diagonal matrices\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .permutation_static => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .diagonal_static, .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: permutation output requires permutation × permutation inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .permutation_static => return @import("matmul/matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: permutation output requires permutation × permutation inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .permutation_static => return @import("matmul/matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: permutation output requires permutation × permutation inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .permutation_sparse => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .diagonal_static, .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: permutation output requires permutation × permutation inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .permutation_static => return @import("matmul/matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: permutation output requires permutation × permutation inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .permutation_static => return @import("matmul/matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: permutation output requires permutation × permutation inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            .builder_sparse => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense, .triangular_static, .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .general_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matgenspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matgenspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matgenspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matgenspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matgenspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matgenspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matgenspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matgenspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .symmetric_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matsymspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matsymspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matsymspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matsymspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matsymspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matsymspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matsymspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matsymspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .hermitian_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matherspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matherspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matherspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matherspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matherspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matherspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matherspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matherspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_mattrispa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_mattrispa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_mattrispa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_mattrispa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_mattrispa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_mattrispa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_mattrispa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_mattrispa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiasta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiasta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiasta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiasta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiasta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiasta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiasta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiasta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiaspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiaspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiaspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiaspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiaspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiaspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiaspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matdiaspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matpersta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matpersta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matpersta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matpersta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matpersta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matpersta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matpersta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matpersta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matperspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matperspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matperspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matperspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matperspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matperspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matperspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/matbuispa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/matbuispa_matperspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.domain(Y)) {
                    .matrix => @compileError("zsl.linalg.matmulInto: matrix output requires two matrix inputs (shape mismatch)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    .vector => @compileError("zsl.linalg.matmulInto: vector × vector outer products are not supported by matmulInto; use linalg.outer instead\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                    else => unreachable,
                },
                else => unreachable,
            },
            else => unreachable,
        },
        .vector => switch (comptime meta.vectorType(O)) {
            .static => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgensta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgensta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgensta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .general_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgenden_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgenden_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgenden_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .general_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgenspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgenspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matgenspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .symmetric_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymsta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymsta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymsta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .symmetric_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymden_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymden_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymden_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .symmetric_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matsymspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .hermitian_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mathersta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mathersta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mathersta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .hermitian_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matherden_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matherden_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matherden_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .hermitian_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matherspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matherspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matherspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .triangular_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattrista_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattrista_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattrista_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattriden_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattriden_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattriden_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattrispa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattrispa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_mattrispa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matdiasta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matdiasta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matdiasta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matdiaspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matdiaspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matdiaspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matpersta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matpersta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matpersta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matperspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matperspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_matperspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.vectorType(X)) {
                    .static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecsta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/vecsta_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecsta_vecspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                else => unreachable,
            },
            .dense => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgensta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgensta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgensta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .general_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgenden_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgenden_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgenden_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .general_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgenspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgenspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matgenspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .symmetric_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymsta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymsta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymsta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .symmetric_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymden_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymden_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymden_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .symmetric_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matsymspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .hermitian_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mathersta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mathersta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mathersta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .hermitian_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matherden_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matherden_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matherden_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .hermitian_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matherspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matherspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matherspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .triangular_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattrista_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattrista_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattrista_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattriden_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattriden_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattriden_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattrispa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattrispa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_mattrispa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matdiasta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matdiasta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matdiasta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matdiaspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matdiaspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matdiaspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matpersta_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matpersta_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matpersta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matperspa_vecsta.zig").matmulIntoUnchecked(o, x, y),
                            .dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matperspa_vecden.zig").matmulIntoUnchecked(o, x, y),
                            .sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_matperspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.vectorType(X)) {
                    .static => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecsta_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .dense => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecden_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matgensta.zig").matmulIntoUnchecked(o, x, y),
                            .general_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matgenden.zig").matmulIntoUnchecked(o, x, y),
                            .general_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matsymsta.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matsymden.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_mathersta.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matherden.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_mattrista.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_dense => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_mattriden.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/vecden_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecden_vecspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            else => unreachable,
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                else => unreachable,
            },
            .sparse => switch (comptime meta.domain(X)) {
                .matrix => switch (comptime meta.matrixType(X)) {
                    .general_static, .general_dense, .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense, .triangular_static, .triangular_dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .general_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_matgenspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        else => unreachable,
                    },
                    .symmetric_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_matsymspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        else => unreachable,
                    },
                    .hermitian_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_matherspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        else => unreachable,
                    },
                    .triangular_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_mattrispa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        else => unreachable,
                    },
                    .diagonal_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_matdiasta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        else => unreachable,
                    },
                    .diagonal_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_matdiaspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        else => unreachable,
                    },
                    .permutation_static => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_matpersta_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        else => unreachable,
                    },
                    .permutation_sparse => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => switch (comptime meta.vectorType(Y)) {
                            .sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_matperspa_vecspa.zig").matmulIntoUnchecked(o, x, y),
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        else => unreachable,
                    },
                    else => unreachable,
                },
                .vector => switch (comptime meta.vectorType(X)) {
                    .static, .dense => switch (comptime meta.domain(Y)) {
                        .matrix => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        .vector => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    .sparse => switch (comptime meta.domain(Y)) {
                        .matrix => switch (comptime meta.matrixType(Y)) {
                            .general_sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_vecspa_matgenspa.zig").matmulIntoUnchecked(o, x, y),
                            .symmetric_sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_vecspa_matsymspa.zig").matmulIntoUnchecked(o, x, y),
                            .hermitian_sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_vecspa_matherspa.zig").matmulIntoUnchecked(o, x, y),
                            .triangular_sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_vecspa_mattrispa.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_static => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_vecspa_matdiasta.zig").matmulIntoUnchecked(o, x, y),
                            .diagonal_sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_vecspa_matdiaspa.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_static => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_vecspa_matpersta.zig").matmulIntoUnchecked(o, x, y),
                            .permutation_sparse => return @import("matmul/vecspa_slow.zig").matmulIntoUnchecked(o, x, y), // return @import("matmul/vecspa_vecspa_matperspa.zig").matmulIntoUnchecked(o, x, y),
                            .builder_sparse => unreachable,
                            else => @compileError("zsl.linalg.matmulInto: sparse output requires sparse, diagonal, or permutation inputs (no dense or static operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        },
                        .vector => @compileError("zsl.linalg.matmulInto: vector output requires exactly one matrix and one vector as inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                        else => unreachable,
                    },
                    else => unreachable,
                },
                else => unreachable,
            },
            else => unreachable,
        },
        else => unreachable,
    }
}

/// Extracts the non-zero indices of a specific mathematical line (row or
/// column), utilizing fast contiguous reads when the natural storage layout
/// permits.
fn yieldLine(comptime O: type, m: anytype, line: usize, dim_len: usize, comptime want_row: bool, comptime callback: anytype, ctx: anytype) void {
    const M = @TypeOf(m);

    if (comptime meta.isVector(M)) {
        if (dim_len > 1 and line == 0) {
            for (0..m.nnz) |p| {
                callback(ctx, m.idx[p]);
            }
            return;
        }

        for (0..dim_len) |other| {
            const r = if (want_row) line else other;
            const c = if (want_row) other else line;
            if (vecutils.hasNZ(m, r, c)) callback(ctx, other);
        }

        return;
    }

    if (comptime !matutils.isCompressed(M)) {
        for (0..dim_len) |other| {
            const r = if (want_row) line else other;
            const c = if (want_row) other else line;
            if (matutils.hasNZ(m, r, c)) callback(ctx, other);
        }
        return;
    }

    if (comptime matutils.isCompressed(M)) {
        const is_row_major = comptime meta.layoutOf(M) == .row_major;
        const mirrored = comptime matutils.needsMirror(O, M);
        const m_diag = comptime matutils.hasImplicitDiag(M);

        if (is_row_major == want_row) {
            const p_start = m.ptr[line];
            const p_end = m.ptr[line + 1];
            for (p_start..p_end) |p| {
                callback(ctx, m.idx[p]);
            }
            if (m_diag and line < dim_len) {
                callback(ctx, line);
            }
            if (!mirrored) return;
        }

        for (0..dim_len) |other| {
            const r = if (want_row) line else other;
            const c = if (want_row) other else line;
            if (matutils.hasNZ(m, r, c)) callback(ctx, other);
        }

        return;
    }
}

/// Resolves the exact non-zero capacity required for the product `x * y`.
/// `work` of length `max(o_rows, o_cols)` for a more efficient algorithm.
pub fn matmulNNZ(comptime O: type, o_rows: usize, o_cols: usize, inner_dim: usize, x: anytype, y: anytype, work: ?[]usize) usize {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    if (work) |w| {
        var total: usize = 0;
        @memset(w, std.math.maxInt(usize));

        const x_row_major = if (comptime meta.isVector(X)) true else if (comptime matutils.isCompressed(X)) meta.layoutOf(X) == .row_major else true;
        const y_row_major = if (comptime meta.isVector(Y)) false else if (comptime matutils.isCompressed(Y)) meta.layoutOf(Y) == .row_major else true;

        const use_row_major = if (x_row_major == y_row_major)
            x_row_major
        else if (comptime meta.isVector(X))
            y_row_major
        else if (comptime meta.isVector(Y))
            x_row_major
        else
            x_row_major;

        if (use_row_major) {
            const RowCtx = struct {
                y: Y,
                w: []usize,
                total: *usize,
                i: usize,
                o_cols: usize,
            };

            for (0..o_rows) |i| {
                var ctx = RowCtx{ .y = y, .w = w, .total = &total, .i = i, .o_cols = o_cols };

                const YieldX = struct {
                    fn cb(cctx: *RowCtx, k: usize) void {
                        const YieldY = struct {
                            fn cb(cctx2: *RowCtx, j: usize) void {
                                if (cctx2.w[j] != cctx2.i) {
                                    cctx2.w[j] = cctx2.i;
                                    cctx2.total.* += 1;
                                }
                            }
                        };

                        yieldLine(O, cctx.y, k, cctx.o_cols, true, YieldY.cb, cctx);
                    }
                };

                yieldLine(O, x, i, inner_dim, true, YieldX.cb, &ctx);
            }
        } else {
            const ColCtx = struct {
                x: X,
                w: []usize,
                total: *usize,
                j: usize,
                o_rows: usize,
            };

            for (0..o_cols) |j| {
                var ctx = ColCtx{ .x = x, .w = w, .total = &total, .j = j, .o_rows = o_rows };

                const YieldY = struct {
                    fn cb(cctx: *ColCtx, k: usize) void {
                        const YieldX = struct {
                            fn cb(cctx2: *ColCtx, i: usize) void {
                                if (cctx2.w[i] != cctx2.j) {
                                    cctx2.w[i] = cctx2.j;
                                    cctx2.total.* += 1;
                                }
                            }
                        };

                        yieldLine(O, cctx.x, k, cctx.o_rows, false, YieldX.cb, cctx);
                    }
                };

                yieldLine(O, y, j, inner_dim, false, YieldY.cb, &ctx);
            }
        }

        return total;
    } else {
        var total: usize = 0;
        for (0..o_rows) |i| {
            for (0..o_cols) |j| {
                for (0..inner_dim) |k| {
                    const x_nz = if (comptime meta.isVector(X)) vecutils.hasNZ(x, i, k) else matutils.hasNZ(x, i, k);
                    const y_nz = if (comptime meta.isVector(Y)) vecutils.hasNZ(y, k, j) else matutils.hasNZ(y, k, j);
                    if (x_nz and y_nz) {
                        total += 1;
                        break;
                    }
                }
            }
        }

        return total;
    }
}
