const std = @import("std");

const meta = @import("../../meta.zig");

const int = @import("../../int.zig");

const numeric = @import("../../numeric.zig");
const matrix = @import("../../matrix.zig");

const matops = @import("../ops.zig");
const matutils = @import("../utils.zig");

pub fn Apply2(comptime X: type, comptime Y: type, comptime op: anytype) type {
    const Op = @TypeOf(op);
    const opinfo = @typeInfo(Op);

    comptime if ((!meta.isMatrix(X) and !meta.isNumeric(X)) or (!meta.isMatrix(Y) and !meta.isNumeric(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 2)
        @compileError("zsl.matrix.Apply2: at least one of X or Y must be a matrix type, the other must be a matrix or a numeric type, and op must be a function of two arguments, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tOp = " ++ @typeName(Op) ++ "\n");

    comptime if (meta.isBuilderMatrix(X) or meta.isBuilderMatrix(Y))
        @compileError("zsl.matrix.Apply2: builder matrix types are not allowed, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tOp = " ++ @typeName(Op) ++ "\n");

    comptime var R = meta.ReturnTypeFromInputs(op, &.{ meta.Numeric(X), meta.Numeric(Y) });
    const rinfo = @typeInfo(R);
    if (rinfo == .error_union)
        R = rinfo.error_union.payload;

    comptime if (!meta.isNumeric(R))
        @compileError("zsl.matrix.Apply2: calling op with arguments of types " ++ @typeName(meta.Numeric(X)) ++ " and " ++ @typeName(meta.Numeric(Y)) ++ " must return a numeric, got\n\tR = " ++ @typeName(R) ++ "\n");

    comptime if (op != numeric.add and op != numeric.sub and op != numeric.mul and op != numeric.div)
        @compileError("zsl.matrix.Apply2: op must be zsl.numeric.add, zsl.numeric.sub, zsl.numeric.mul or zsl.numeric.div, got\n\top: " ++ @typeName(Op) ++ "\n");

    comptime if (meta.isStaticMatrix(X) and meta.isStaticMatrix(Y) and (X.rows != Y.rows or X.cols != Y.cols))
        @compileError("zsl.matrix.Apply2: static matrix types X and Y must have equal dimensions, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n");

    switch (comptime meta.matrixType(X)) {
        .general_static => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(X).?),
        .general_dense => switch (comptime meta.matrixStorage(Y)) {
            .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
            else => return matrix.general.Dense(R, meta.layoutOf(X).?),
        },
        .general_sparse => switch (comptime meta.matrixStorage(Y)) {
            .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
            .dense => return matrix.general.Dense(R, meta.layoutOf(X).?),
            else => return matrix.general.Sparse(R, meta.layoutOf(X).?),
        },
        .symmetric_static => switch (comptime meta.matrixKind(Y)) {
            .symmetric, .diagonal, .numeric => return matrix.symmetric.Static(X.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
            .hermitian => switch (comptime meta.isComplex(meta.Numeric(X))) {
                true => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(X).?),
                false => return matrix.hermitian.Static(X.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
            },
            else => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(X).?),
        },
        .symmetric_dense => switch (comptime meta.matrixKind(Y)) {
            .symmetric, .diagonal, .numeric => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.symmetric.Static(Y.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
                else => return matrix.symmetric.Dense(R, meta.uploOf(X).?, meta.layoutOf(X).?),
            },
            .hermitian => switch (comptime meta.isComplex(meta.Numeric(X))) {
                true => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                    else => return matrix.general.Dense(R, meta.layoutOf(X).?),
                },
                false => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.hermitian.Static(Y.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
                    else => return matrix.hermitian.Dense(R, meta.uploOf(X).?, meta.layoutOf(X).?),
                },
            },
            else => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                else => return matrix.general.Dense(R, meta.layoutOf(X).?),
            },
        },
        .symmetric_sparse => switch (comptime meta.matrixKind(Y)) {
            .symmetric, .diagonal, .numeric => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.symmetric.Static(Y.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
                .dense => return matrix.symmetric.Dense(R, meta.uploOf(X).?, meta.layoutOf(X).?),
                else => return matrix.symmetric.Sparse(R, meta.uploOf(X).?, meta.layoutOf(X).?),
            },
            .hermitian => switch (comptime meta.isComplex(meta.Numeric(X))) {
                true => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                    .dense => return matrix.general.Dense(R, meta.layoutOf(X).?),
                    else => return matrix.general.Sparse(R, meta.layoutOf(X).?),
                },
                false => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.hermitian.Static(Y.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
                    .dense => return matrix.hermitian.Dense(R, meta.uploOf(X).?, meta.layoutOf(X).?),
                    else => return matrix.hermitian.Sparse(R, meta.uploOf(X).?, meta.layoutOf(X).?),
                },
            },
            else => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                .dense => return matrix.general.Dense(R, meta.layoutOf(X).?),
                else => return matrix.general.Sparse(R, meta.layoutOf(X).?),
            },
        },
        .hermitian_static => switch (comptime meta.matrixKind(Y)) {
            .symmetric, .diagonal, .numeric => switch (comptime meta.isComplex(meta.Numeric(Y))) {
                true => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(X).?),
                false => return matrix.hermitian.Static(X.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
            },
            .hermitian => return matrix.hermitian.Static(X.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
            else => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(X).?),
        },
        .hermitian_dense => switch (comptime meta.matrixKind(Y)) {
            .symmetric, .diagonal, .numeric => switch (comptime meta.isComplex(meta.Numeric(Y))) {
                true => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                    else => return matrix.general.Dense(R, meta.layoutOf(X).?),
                },
                false => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.hermitian.Static(Y.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
                    else => return matrix.hermitian.Dense(R, meta.uploOf(X).?, meta.layoutOf(X).?),
                },
            },
            .hermitian => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.hermitian.Static(Y.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
                else => return matrix.hermitian.Dense(R, meta.uploOf(X).?, meta.layoutOf(X).?),
            },
            else => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                else => return matrix.general.Dense(R, meta.layoutOf(X).?),
            },
        },
        .hermitian_sparse => switch (comptime meta.matrixKind(Y)) {
            .symmetric, .diagonal, .numeric => switch (comptime meta.isComplex(meta.Numeric(Y))) {
                true => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                    .dense => return matrix.general.Dense(R, meta.layoutOf(X).?),
                    else => return matrix.general.Sparse(R, meta.layoutOf(X).?),
                },
                false => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.hermitian.Static(Y.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
                    .dense => return matrix.hermitian.Dense(R, meta.uploOf(X).?, meta.layoutOf(X).?),
                    else => return matrix.hermitian.Sparse(R, meta.uploOf(X).?, meta.layoutOf(X).?),
                },
            },
            .hermitian => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.hermitian.Static(Y.rows, R, meta.uploOf(X).?, meta.layoutOf(X).?),
                .dense => return matrix.hermitian.Dense(R, meta.uploOf(X).?, meta.layoutOf(X).?),
                else => return matrix.hermitian.Sparse(R, meta.uploOf(X).?, meta.layoutOf(X).?),
            },
            else => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                .dense => return matrix.general.Dense(R, meta.layoutOf(X).?),
                else => return matrix.general.Sparse(R, meta.layoutOf(X).?),
            },
        },
        .triangular_static => switch (comptime meta.matrixKind(Y)) {
            .triangular => switch (comptime meta.uploOf(X).? == meta.uploOf(Y).?) {
                true => return matrix.triangular.Static(X.rows, X.cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                false => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(X).?),
            },
            .diagonal, .numeric => return matrix.triangular.Static(X.rows, X.cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
            else => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(X).?),
        },
        .triangular_dense => switch (comptime meta.matrixKind(Y)) {
            .triangular => switch (comptime meta.uploOf(X).? == meta.uploOf(Y).?) {
                true => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.triangular.Static(Y.rows, Y.cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                    else => return matrix.triangular.Dense(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                },
                false => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                    else => return matrix.general.Dense(R, meta.layoutOf(X).?),
                },
            },
            .diagonal, .numeric => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.triangular.Static(Y.rows, Y.cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                else => return matrix.triangular.Dense(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
            },
            else => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                else => return matrix.general.Dense(R, meta.layoutOf(X).?),
            },
        },
        .triangular_sparse => switch (comptime meta.matrixKind(Y)) {
            .triangular => switch (comptime meta.uploOf(X).? == meta.uploOf(Y).?) {
                true => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.triangular.Static(Y.rows, Y.cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                    .dense => return matrix.triangular.Dense(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                    else => return matrix.triangular.Sparse(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                },
                false => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                    .dense => return matrix.general.Dense(R, meta.layoutOf(X).?),
                    else => return matrix.general.Sparse(R, meta.layoutOf(X).?),
                },
            },
            .diagonal, .numeric => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.triangular.Static(Y.rows, Y.cols, R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                .dense => return matrix.triangular.Dense(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
                else => return matrix.triangular.Sparse(R, meta.uploOf(X).?, .non_unit, meta.layoutOf(X).?),
            },
            else => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(X).?),
                .dense => return matrix.general.Dense(R, meta.layoutOf(X).?),
                else => return matrix.general.Sparse(R, meta.layoutOf(X).?),
            },
        },
        .diagonal_static => switch (comptime meta.matrixKind(Y)) {
            .symmetric => return matrix.symmetric.Static(X.rows, R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
            .hermitian => switch (comptime meta.isComplex(meta.Numeric(X))) {
                true => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(Y).?),
                false => return matrix.hermitian.Static(X.rows, R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
            },
            .triangular => return matrix.triangular.Static(X.rows, X.cols, R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
            .diagonal, .numeric => return matrix.diagonal.Static(X.rows, X.cols, R),
            else => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(Y) orelse matrix.Layout.default),
        },
        .diagonal_sparse => switch (comptime meta.matrixKind(Y)) {
            .symmetric => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.symmetric.Static(Y.rows, R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                .dense => return matrix.symmetric.Dense(R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                else => return matrix.symmetric.Sparse(R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
            },
            .hermitian => switch (comptime meta.isComplex(meta.Numeric(X))) {
                true => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(Y).?),
                    .dense => return matrix.general.Dense(R, meta.layoutOf(Y).?),
                    else => return matrix.general.Sparse(R, meta.layoutOf(Y).?),
                },
                false => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.hermitian.Static(Y.rows, R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                    .dense => return matrix.hermitian.Dense(R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                    else => return matrix.hermitian.Sparse(R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                },
            },
            .triangular => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.triangular.Static(Y.rows, Y.cols, R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
                .dense => return matrix.triangular.Dense(R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
                else => return matrix.triangular.Sparse(R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
            },
            .diagonal, .numeric => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.diagonal.Static(Y.rows, Y.cols, R),
                else => return matrix.diagonal.Sparse(R),
            },
            else => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(Y) orelse matrix.Layout.default),
                .dense => return matrix.general.Dense(R, meta.layoutOf(Y) orelse matrix.Layout.default),
                else => return matrix.general.Sparse(R, meta.layoutOf(Y) orelse matrix.Layout.default),
            },
        },
        .permutation_static => return matrix.general.Static(X.rows, X.cols, R, meta.layoutOf(Y) orelse matrix.Layout.default),
        .permutation_sparse => switch (comptime meta.matrixStorage(Y)) {
            .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(Y) orelse matrix.Layout.default),
            .dense => return matrix.general.Dense(R, meta.layoutOf(Y).?),
            else => return matrix.general.Sparse(R, meta.layoutOf(Y) orelse matrix.Layout.default),
        },
        .builder_sparse => unreachable,
        .numeric => switch (comptime meta.matrixKind(Y)) {
            .general, .permutation => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(Y) orelse matrix.Layout.default),
                .dense => return matrix.general.Dense(R, meta.layoutOf(Y).?),
                else => return matrix.general.Sparse(R, meta.layoutOf(Y) orelse matrix.Layout.default),
            },
            .symmetric => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.symmetric.Static(Y.rows, R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                .dense => return matrix.symmetric.Dense(R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                else => return matrix.symmetric.Sparse(R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
            },
            .hermitian => switch (comptime meta.isComplex(X)) {
                true => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.general.Static(Y.rows, Y.cols, R, meta.layoutOf(Y).?),
                    .dense => return matrix.general.Dense(R, meta.layoutOf(Y).?),
                    else => return matrix.general.Sparse(R, meta.layoutOf(Y).?),
                },
                false => switch (comptime meta.matrixStorage(Y)) {
                    .static => return matrix.hermitian.Static(Y.rows, R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                    .dense => return matrix.hermitian.Dense(R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                    else => return matrix.hermitian.Sparse(R, meta.uploOf(Y).?, meta.layoutOf(Y).?),
                },
            },
            .triangular => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.triangular.Static(Y.rows, Y.cols, R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
                .dense => return matrix.triangular.Dense(R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
                else => return matrix.triangular.Sparse(R, meta.uploOf(Y).?, .non_unit, meta.layoutOf(Y).?),
            },
            .diagonal => switch (comptime meta.matrixStorage(Y)) {
                .static => return matrix.diagonal.Static(Y.rows, Y.cols, R),
                else => return matrix.diagonal.Sparse(R),
            },
            else => unreachable,
        },
    }
}

/// Applies a binary operation elementwise between two matrices, or between a
/// matrix and a numeric.
///
/// The result inherits its memory layout from the inputs, i.e., if the input
/// layouts mismatch, the left operand (`x`) strictly dictates the output
/// layout, unless it provides no layout information. For more control over
/// layouts, use `matrix.apply2Into`.
///
/// This function is intended for when the result matrix's dimension is known at
/// compile time, i.e., at least one of the inputs is a static matrix. For two
/// static matrices, or a static matrix and a numeric, dimension checks are
/// performed at compile time, for any other combination, dimension checks are
/// performed at runtime throught `std.debug.assert`.
///
/// ## Signature
/// ```zig
/// matrix.apply2(x: X, y: Y, op: Op) !matrix.Apply2(X, Y, op)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
/// * `op` (`comptime anytype`): A binary numeric function to apply elementwise
///   to `x` and `y`.
///
/// ## Returns
/// `matrix.Apply2(@TypeOf(x), @TypeOf(y), op)`: The result of the operation.
///
/// ## Errors
/// * `matrix.Error.DimensionMismatch`: If the two matrices do not have the same
///   dimensions. Can only happen if both operands are matrices.
pub fn apply2(x: anytype, y: anytype, comptime op: anytype) !matops.Apply2(@TypeOf(x), @TypeOf(y), op) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Op: type = @TypeOf(op);
    const R: type = matops.Apply2(X, Y, op);

    if (comptime meta.isDenseMatrix(R) or meta.isSparseMatrix(R))
        @compileError("zsl.matrix.apply2: the result cannot be a heap-allocated matrix type, i.e., at least one input must be a static matrix, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n\tresult: " ++ @typeName(R) ++ "\nFor these inputs use zsl.matrix.apply2Alloc instead.");

    const x_rows_optional: ?usize = if (comptime meta.isMatrix(X)) (if (comptime meta.isStaticMatrix(X)) X.rows else x.rows) else null;
    const x_cols_optional: ?usize = if (comptime meta.isMatrix(X)) (if (comptime meta.isStaticMatrix(X)) X.cols else x.cols) else null;
    const y_rows = if (comptime meta.isMatrix(Y)) (if (comptime meta.isStaticMatrix(Y)) Y.rows else y.rows) else x_rows_optional.?;
    const y_cols = if (comptime meta.isMatrix(Y)) (if (comptime meta.isStaticMatrix(Y)) Y.cols else y.cols) else x_cols_optional.?;
    const x_rows = x_rows_optional orelse y_rows;
    const x_cols = x_cols_optional orelse y_cols;

    if (comptime !(meta.isStaticMatrix(X) or meta.isNumeric(X)) or !(meta.isStaticMatrix(Y) or meta.isNumeric(Y))) {
        if (x_rows != y_rows or x_cols != y_cols)
            return matrix.Error.DimensionMismatch;
    }

    var result = R.init;

    matops.apply2IntoUnchecked(
        &result,
        x,
        y,
        switch (comptime op) {
            numeric.add => numeric.addInto,
            numeric.sub => numeric.subInto,
            numeric.mul => numeric.mulInto,
            numeric.div => numeric.divInto,
            else => unreachable,
        },
    );

    return result;
}

/// Applies a binary operation elementwise between two matrices, or between a
/// matrix and a numeric, without performing any dimension checks.
///
/// The result inherits its memory layout from the inputs, i.e., if the input
/// layouts mismatch, the left operand (`x`) strictly dictates the output
/// layout, unless it provides no layout information. For more control over
/// layouts, use `matrix.apply2Into`.
///
/// This function is intended for when the result matrix's dimension is known at
/// compile time, i.e., at least one of the inputs is a static matrix.
///
/// ## Signature
/// ```zig
/// matrix.apply2Unchecked(x: X, y: Y, op: Op) matrix.Apply2(X, Y, op)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
/// * `op` (`comptime anytype`): A binary numeric function to apply elementwise
///   to `x` and `y`.
///
/// ## Returns
/// `matrix.Apply2(@TypeOf(x), @TypeOf(y), op)`: The result of the operation.
pub fn apply2Unchecked(x: anytype, y: anytype, comptime op: anytype) matops.Apply2(@TypeOf(x), @TypeOf(y), op) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Op: type = @TypeOf(op);
    const R: type = matops.Apply2(X, Y, op);

    if (comptime meta.isDenseMatrix(R) or meta.isSparseMatrix(R))
        @compileError("zsl.matrix.apply2: the result cannot be a heap-allocated matrix type, i.e., at least one input must be a static matrix, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top: " ++ @typeName(Op) ++ "\n\tresult: " ++ @typeName(R) ++ "\nFor these inputs use zsl.matrix.apply2Alloc instead.");

    var result = R.init;

    matops.apply2IntoUnchecked(
        &result,
        x,
        y,
        switch (comptime op) {
            numeric.add => numeric.addInto,
            numeric.sub => numeric.subInto,
            numeric.mul => numeric.mulInto,
            numeric.div => numeric.divInto,
            else => unreachable,
        },
    );

    return result;
}

/// Applies a binary operation elementwise between two matrices, or between a
/// matrix and a numeric, dynamically allocating memory for the result.
///
/// The result inherits its memory layout from the inputs, i.e., if the input
/// layouts mismatch, the left operand (`x`) strictly dictates the output
/// layout, unless it provides no layout information. For more control over
/// layouts, use `matrix.apply2Into`.
///
/// This function is intended for when the result matrix's dimension is known
/// at runtime, i.e., none of the inputs is a static matrix.
///
/// For two sparse matrices, or a sparse matrix and a numeric, the operation is
/// only applied to the indices where at least one of the matrices has a
/// non-zero element.
///
/// ## Signature
/// ```zig
/// matrix.apply2Alloc(allocator: std.mem.Allocator, x: X, y: Y, op: Op) !matrix.Apply2(X, Y, op)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
/// * `op` (`comptime anytype`): A binary numeric function to apply elementwise
///   to `x` and `y`.
///
/// ## Returns
/// `matrix.Apply2(@TypeOf(x), @TypeOf(y), op)`: The result of the operation.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `matrix.Error.DimensionMismatch`: If the two matrices do not have the same
///   dimensions. Can only happen if both operands are matrices.
pub fn apply2Alloc(allocator: std.mem.Allocator, x: anytype, y: anytype, comptime op: anytype) !matops.Apply2(@TypeOf(x), @TypeOf(y), op) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = matops.Apply2(X, Y, op);

    const x_rows_optional: ?usize = if (comptime meta.isMatrix(X)) (if (comptime meta.isStaticMatrix(X)) X.rows else x.rows) else null;
    const x_cols_optional: ?usize = if (comptime meta.isMatrix(X)) (if (comptime meta.isStaticMatrix(X)) X.cols else x.cols) else null;
    const y_rows = if (comptime meta.isMatrix(Y)) (if (comptime meta.isStaticMatrix(Y)) Y.rows else y.rows) else x_rows_optional.?;
    const y_cols = if (comptime meta.isMatrix(Y)) (if (comptime meta.isStaticMatrix(Y)) Y.cols else y.cols) else x_cols_optional.?;
    const x_rows = x_rows_optional orelse y_rows;
    const x_cols = x_cols_optional orelse y_cols;

    if (comptime !(meta.isStaticMatrix(X) or meta.isNumeric(X)) or !(meta.isStaticMatrix(Y) or meta.isNumeric(Y))) {
        if (x_rows != y_rows or x_cols != y_cols)
            return matrix.Error.DimensionMismatch;
    }

    var result = switch (comptime meta.matrixStorage(R)) {
        .static => R.init,
        .dense => switch (comptime meta.matrixKind((R))) {
            .general, .triangular => try R.init(allocator, x_rows, x_cols),
            .symmetric, .hermitian => try R.init(allocator, x_rows),
            else => unreachable,
        },
        .sparse => switch (comptime meta.matrixKind(R)) {
            .general, .triangular => try R.init(allocator, x_rows, x_cols, apply2NNZ(R, x_rows, x_cols, x, y)),
            .symmetric, .hermitian => try R.init(allocator, x_rows, apply2NNZ(R, x_rows, x_cols, x, y)),
            .diagonal => try R.init(allocator, x_rows, x_cols),
            else => unreachable,
        },
        .numeric => unreachable,
    };

    matops.apply2IntoUnchecked(
        &result,
        x,
        y,
        switch (comptime op) {
            numeric.add => numeric.addInto,
            numeric.sub => numeric.subInto,
            numeric.mul => numeric.mulInto,
            numeric.div => numeric.divInto,
            else => unreachable,
        },
    );

    return result;
}

/// Applies a binary into operation elementwise between an output and two input
/// matrices, or between an output matrix, an input matrix and an input numeric.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted for static or dense outputs, or sparse outputs when both inputs
/// are not matrices. Any other form of memory overlap might yield incorrect
/// results.
///
/// For three static matrices, or a static output matrix, an input static matrix
/// and an input numeric, the function cannot return an error unless opInto can.
///
/// For two input sparse matrices, or an input sparse matrix and an input
/// numeric, the operation is only applied to the indices where at least one of
/// the matrices has a non-zero element.
///
/// ## Signature
/// ```zig
/// matrix.apply2Into(*O, x: X, y: Y, opInto: OpInto) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left input operand.
/// * `y` (`anytype`): The right input operand.
/// * `opInto` (`comptime anytype`): An ito binary numeric function to apply
///   elementwise to `o`, `x` and `y`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `matrix.Error.DimensionMismatch`: If the matrices do not have the same
///   dimensions.
pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const OpInto: type = @TypeOf(opInto);
    const opinfo = @typeInfo(OpInto);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isMatrix(meta.Child(O)) or
        (!meta.isMatrix(X) and !meta.isNumeric(X)) or (!meta.isMatrix(Y) and !meta.isNumeric(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 3)
        @compileError("zsl.matrix.apply2Into: o must be a mutable one-itme pointer to a matrix, at least one of x or y must be a matrix, the other must be a matrix or a numeric, and opInto must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    comptime if (meta.isBuilderMatrix(X) or meta.isBuilderMatrix(Y))
        @compileError("zsl.matrix.apply2Into: builder matrix types are not allowed as inputs, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (opInto != numeric.addInto and opInto != numeric.subInto and opInto != numeric.mulInto and opInto != numeric.divInto)
        @compileError("zsl.matrix.apply2IntoUnchecked: opInto must be zsl.numeric.addInto, zsl.numeric.subInto, zsl.numeric.mulInto or zsl.numeric.divInto, got\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    O = meta.Child(O);

    const o_rows = if (comptime meta.isStaticMatrix(O)) O.rows else o.rows;
    const o_cols = if (comptime meta.isStaticMatrix(O)) O.cols else o.cols;
    const x_rows = if (comptime meta.isMatrix(X)) (if (comptime meta.isStaticMatrix(X)) X.rows else x.rows) else o_rows;
    const x_cols = if (comptime meta.isMatrix(X)) (if (comptime meta.isStaticMatrix(X)) X.cols else x.cols) else o_cols;
    const y_rows = if (comptime meta.isMatrix(Y)) (if (comptime meta.isStaticMatrix(Y)) Y.rows else y.rows) else o_rows;
    const y_cols = if (comptime meta.isMatrix(Y)) (if (comptime meta.isStaticMatrix(Y)) Y.cols else y.cols) else o_cols;

    if (comptime meta.isStaticMatrix(O) and
        (meta.isStaticMatrix(X) or meta.isNumeric(X)) and
        (meta.isStaticMatrix(Y) or meta.isNumeric(Y)))
    {
        if (comptime o_rows != x_rows or o_rows != y_rows or o_cols != x_cols or o_cols != y_cols)
            if (comptime meta.isMatrix(X) and meta.isMatrix(Y))
                @compileError("zsl.matrix.apply2: static matrices o, x and y must have equal compile time dimensions, got\n\to: *" ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n")
            else if (comptime meta.isMatrix(X))
                @compileError("zsl.matrix.apply2: static matrices o and x must have equal compile time dimensions, got\n\to: *" ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n")
            else
                @compileError("zsl.matrix.apply2: static matrices o and y must have equal compile time dimensions, got\n\to: *" ++
                    @typeName(O) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");
    } else {
        if (o_rows != x_rows or o_rows != y_rows or o_cols != x_cols or o_cols != y_cols)
            return matrix.Error.DimensionMismatch;
    }

    if (comptime meta.isSparseMatrix(O) and !meta.isDiagonalMatrix(O) and !meta.isPermutationMatrix(O)) {
        if (comptime (meta.isStaticMatrix(X) and !meta.isDiagonalMatrix(X) and !meta.isPermutationMatrix(X)) or meta.isDenseMatrix(X) or
            (meta.isStaticMatrix(Y) and !meta.isDiagonalMatrix(Y) and !meta.isPermutationMatrix(Y)) or meta.isDenseMatrix(Y))
            @compileError("zsl.matrix.apply2Into: o cannot point to a sparse matrix if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

        const nnz = apply2NNZ(O, o_rows, o_cols, x, y);

        if (comptime meta.isBuilderMatrix(O)) {
            if (o._dlen < nnz or o._rlen < nnz or o._clen < nnz)
                return matrix.Error.InsufficientSpace;
        } else {
            if (o._dlen < nnz or o._ilen < nnz)
                return matrix.Error.InsufficientSpace;
        }
    }

    return apply2IntoUnchecked(o, x, y, opInto);
}

/// Applies a binary into operation elementwise between an output and two input
/// matrices, or between an output matrix, an input matrix and an input numeric,
/// without performing any dimension checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted for static or dense outputs, or sparse outputs when both inputs
/// are not matrices. Any other form of memory overlap might yield incorrect
/// results.
///
/// For two input sparse matrices, or an input sparse matrix and an input
/// numeric, the operation is only applied to the indices where at least one of
/// the matrices has a non-zero element.
///
/// ## Signature
/// ```zig
/// matrix.apply2IntoUnchecked(*O, x: X, y: Y, opInto: OpInto) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left input operand.
/// * `y` (`anytype`): The right input operand.
/// * `opInto` (`comptime anytype`): An ito binary numeric function to apply
///   elementwise to `o`, `x` and `y`.
///
/// ## Returns
/// `void`
pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Op: type = @TypeOf(opInto);
    const opinfo = @typeInfo(Op);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isMatrix(meta.Child(O)) or
        (!meta.isMatrix(X) and !meta.isNumeric(X)) or (!meta.isMatrix(Y) and !meta.isNumeric(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 3)
        @compileError("zsl.matrix.apply2Into: o must be a mutable one-itme pointer to a matrix, at least one of x or y must be a matrix, the other must be a matrix or a numeric, and opInto must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(Op) ++ "\n");

    comptime if (meta.isBuilderMatrix(X) or meta.isBuilderMatrix(Y))
        @compileError("zsl.matrix.apply2Into: builder matrix types are not allowed as inputs, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    O = meta.Child(O);

    switch (comptime meta.matrixType(O)) {
        .general_static => switch (comptime meta.matrixType(X)) {
            .general_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgensta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .general_dense => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .general_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matgenspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .symmetric_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymsta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .hermitian_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mathersta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .hermitian_dense => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matherspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .triangular_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrista_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .triangular_dense => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattriden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_mattrispa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .permutation_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matpersta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .permutation_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_matperspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .builder_sparse => unreachable,
            .numeric => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgensta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgensta_num_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => unreachable,
            },
        },
        .general_dense => switch (comptime meta.matrixType(X)) {
            .general_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgensta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .general_dense => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_matgenden_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_matgenden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_matgenden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_matgenden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_matgenden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenden_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_matgenden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .general_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matgenspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .symmetric_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymsta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_matsymden_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_matsymden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_matsymden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_matsymden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_matsymden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymden_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_matsymden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .hermitian_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mathersta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .hermitian_dense => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_matherden_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_matherden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_matherden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_matherden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_matherden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherden_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_matherden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matherspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .triangular_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrista_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .triangular_dense => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_mattriden_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_mattriden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_mattriden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_mattriden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_mattriden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattriden_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_mattriden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_mattrispa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_matdiaspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_matdiaspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_matdiaspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_matdiaspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matdiaspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .permutation_static => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matpersta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .permutation_sparse => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_matperspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
            },
            .builder_sparse => unreachable,
            .numeric => switch (comptime meta.matrixType(Y)) {
                .general_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_matgensta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_dense => return @import("apply2/matgenden_num_matgenden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .general_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matgenden_num_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matgenden_num_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_dense => return @import("apply2/matgenden_num_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenden_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenden_num_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => unreachable,
            },
        },
        .general_sparse => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense, .triangular_static, .triangular_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .general_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matgenspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matherspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_mattrispa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .permutation_static => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matpersta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .permutation_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_matperspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .builder_sparse => unreachable,
            .numeric => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_num_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_num_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_num_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_num_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matgenspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matgenspa_num_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
        },
        .symmetric_static => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .symmetric_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymsta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymsta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymsta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymsta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymsta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymsta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiasta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiasta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiaspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiaspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_num_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_num_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymsta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymsta_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .symmetric_dense => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .symmetric_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymsta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymsta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymsta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymsta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymsta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymsta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiasta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiasta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiaspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiaspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_num_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_dense => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_num_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymden_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .symmetric_sparse => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .symmetric_static, .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense, .symmetric_sparse, .diagonal_static, .diagonal_sparse, .numeric => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .builder_sparse => unreachable,
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matsymspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matsymspa_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: symmetric output requires symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .hermitian_static => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => switch (comptime meta.matrixType(Y)) {
                .builder_sparse => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymsta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymden_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_mathersta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_mathersta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_mathersta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_mathersta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_mathersta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_mathersta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_mathersta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_mathersta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_mathersta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_matherden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_matherden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_matherden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherden_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_matherspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_matherspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/mathersta_matherspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matherspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_num_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_num_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_num_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_num_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_num_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mathersta_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mathersta_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .hermitian_dense => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .symmetric_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymsta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymden_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_mathersta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_mathersta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_mathersta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_mathersta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_mathersta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_mathersta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_mathersta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_mathersta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_mathersta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherden_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherden_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherden_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_matherden_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_matherden_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_matherden_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherden_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_matherspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_dense => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_matherspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherden_matherspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matherspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .symmetric_static => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_num_matsymsta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_dense => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_num_matsymden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_static => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_num_mathersta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_num_matherden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_num_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherden_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .hermitian_sparse => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense => switch (comptime meta.matrixType(Y)) {
                .general_static, .general_dense, .general_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .builder_sparse => unreachable,
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matsymspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matherspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matherspa_matherspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matherspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matherspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matherspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiasta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiaspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .symmetric_sparse => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_num_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2Into: hermitian output requires symmetric, diagonal and numeric inputs to have a real element type\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/matherspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/matherspa_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: hermitian output requires hermitian, symmetric, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .triangular_static => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .triangular_static => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrista_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrista_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrista_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrista_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrista_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrista_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .triangular_dense => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattriden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattriden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattriden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattriden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattriden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattriden_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrispa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrispa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrispa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrispa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrispa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_mattrispa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiasta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiasta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiasta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiaspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiaspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiaspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_num_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_num_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_num_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrista_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrista_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .triangular_dense => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .triangular_static => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrista_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrista_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrista_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrista_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrista_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrista_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .triangular_dense => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattriden_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattriden_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattriden_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattriden_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattriden_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattriden_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrispa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrispa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrispa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrispa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrispa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_mattrispa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiasta_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiasta_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiasta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiaspa_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiaspa_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiaspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .triangular_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_num_mattrista.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_num_mattriden.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_num_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattriden_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattriden_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .triangular_sparse => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .triangular_static, .triangular_dense => switch (comptime meta.matrixType(Y)) {
                .triangular_static, .triangular_dense, .triangular_sparse, .diagonal_static, .diagonal_sparse, .numeric => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .builder_sparse => unreachable,
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .triangular_static, .triangular_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_mattrispa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_mattrispa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_mattrispa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_mattrispa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .triangular_static, .triangular_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_matdiasta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .triangular_static, .triangular_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_matdiaspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse => unreachable,
                .numeric => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .triangular_static, .triangular_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_num_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_static => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .diagonal_sparse => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2Into: triangular operands must share uplo, and output must not be unit-diagonal\n\to: *" ++
                            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

                    return @import("apply2/mattrispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto); // return @import("apply2/mattrispa_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto);
                },
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: triangular output requires triangular, diagonal, or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .diagonal_static => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: diagonal output requires diagonal or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .diagonal_static => return @import("apply2/matdiasta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiasta_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matdiasta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiasta_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matdiasta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiasta_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: diagonal output requires diagonal or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .diagonal_static => return @import("apply2/matdiasta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiasta_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matdiasta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiasta_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matdiasta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiasta_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: diagonal output requires diagonal or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .diagonal_static => return @import("apply2/matdiasta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiasta_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matdiasta_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiasta_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: diagonal output requires diagonal or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .diagonal_sparse => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .general_sparse, .symmetric_static, .symmetric_dense, .symmetric_sparse, .hermitian_static, .hermitian_dense, .hermitian_sparse, .triangular_static, .triangular_dense, .triangular_sparse, .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: diagonal output requires diagonal or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .diagonal_static => return @import("apply2/matdiaspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiaspa_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matdiaspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiaspa_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matdiaspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiaspa_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: diagonal output requires diagonal or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .diagonal_static => return @import("apply2/matdiaspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiaspa_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matdiaspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiaspa_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matdiaspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiaspa_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: diagonal output requires diagonal or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .diagonal_static => return @import("apply2/matdiaspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiaspa_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matdiaspa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matdiaspa_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: diagonal output requires diagonal or numeric inputs\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .permutation_static, .permutation_sparse => @compileError("zsl.matrix.apply2Into: permutation matrices cannot be used as output\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
        .builder_sparse => switch (comptime meta.matrixType(X)) {
            .general_static, .general_dense, .symmetric_static, .symmetric_dense, .hermitian_static, .hermitian_dense, .triangular_static, .triangular_dense => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            .general_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matgenspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matsymspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matherspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_mattrispa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_static => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiasta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .diagonal_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matdiaspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .permutation_static => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matpersta_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .permutation_sparse => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_matperspa_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .general_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_num_matgenspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .symmetric_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_num_matsymspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .hermitian_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_num_matherspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .triangular_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_num_mattrispa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_num_matdiasta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .diagonal_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_num_matdiaspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_static => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_num_matpersta.zig").apply2IntoUnchecked(o, x, y, opInto),
                .permutation_sparse => return @import("apply2/matbuispa_slow.zig").apply2IntoUnchecked(o, x, y, opInto), // return @import("apply2/matbuispa_num_matperspa.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse, .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2Into: sparse output requires sparse or numeric inputs (no dense operand)\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y)),
            },
            else => unreachable,
        },
        .numeric => unreachable,
    }
}

/// Derives the non-zero capacity of a trivial matrix (diagonal or permutation)
/// bounded by the target output dimensions.
fn trivialCount(t: anytype, o_rows: usize, o_cols: usize) usize {
    return switch (comptime meta.matrixType(@TypeOf(t))) {
        .diagonal_static, .diagonal_sparse => int.min(o_rows, o_cols),
        .permutation_static, .permutation_sparse => o_rows,
        else => unreachable,
    };
}

/// Computes the total materialized non-zero count of a compressed matrix `m`,
/// expanding off-diagonal entries if `O` requires it.
fn ownCount(comptime O: type, o_rows: usize, o_cols: usize, m: anytype) usize {
    const M = @TypeOf(m);

    if (comptime !matutils.needsMirror(O, M)) {
        var total: usize = m.nnz;
        if (matutils.hasImplicitDiag(M)) total += int.min(o_rows, o_cols);
        return total;
    }

    const row_major = comptime meta.layoutOf(M) == .row_major;
    const lines: usize = if (row_major) o_rows else o_cols;
    var total: usize = m.nnz;
    var line: usize = 0;
    while (line < lines) : (line += 1) {
        var p: usize = m.ptr[line];
        const end: usize = m.ptr[line + 1];
        while (p < end) : (p += 1) {
            const other = m.idx[p];
            const row = if (row_major) line else other;
            const col = if (row_major) other else line;
            if (row != col) total += 1;
        }
    }

    return total;
}

/// Calculates the effective non-zero contribution of a single matrix operand
/// when projected into the target output structure `O`.
fn singleOperandCount(comptime O: type, o_rows: usize, o_cols: usize, m: anytype) usize {
    if (comptime matutils.isCompressed(@TypeOf(m)))
        return ownCount(O, o_rows, o_cols, m);

    return trivialCount(m, o_rows, o_cols);
}

/// A single-line coordinate iterator that yields implicit unit-diagonal entries
/// during sequential reads.
const LineCursor = struct {
    idx: [*]usize,
    p: usize,
    end: usize,
    diag: ?usize,

    fn peek(self: *const LineCursor) ?usize {
        if (self.diag) |d| {
            if (self.p >= self.end or d <= self.idx[self.p]) return d;
            return self.idx[self.p];
        }

        return if (self.p < self.end) self.idx[self.p] else null;
    }

    fn advance(self: *LineCursor) void {
        if (self.diag) |d| {
            if (self.p >= self.end or d <= self.idx[self.p]) {
                self.diag = null;
                return;
            }
        }

        self.p += 1;
    }

    fn remaining(self: *const LineCursor) usize {
        const explicit = self.end - self.p;

        return if (self.diag != null) explicit + 1 else explicit;
    }
};

/// Performs a linear O(nnz_x + nnz_y) coordinate union scan across two
/// structurally compatible compressed matrices.
fn mergeCount(o_rows: usize, o_cols: usize, x: anytype, y: anytype) usize {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    const row_major = comptime meta.layoutOf(X) == .row_major;
    const lines: usize = if (row_major) o_rows else o_cols;
    const other_dim: usize = if (row_major) o_cols else o_rows;

    var total: usize = 0;
    var line: usize = 0;
    while (line < lines) : (line += 1) {
        var cx = LineCursor{
            .idx = x.idx,
            .p = x.ptr[line],
            .end = x.ptr[line + 1],
            .diag = if (comptime matutils.hasImplicitDiag(X)) (if (line < other_dim) line else null) else null,
        };

        var cy = LineCursor{
            .idx = y.idx,
            .p = y.ptr[line],
            .end = y.ptr[line + 1],
            .diag = if (comptime matutils.hasImplicitDiag(Y)) (if (line < other_dim) line else null) else null,
        };

        while (true) {
            const ix = cx.peek() orelse break;
            const iy = cy.peek() orelse break;

            total += 1;
            if (ix == iy) {
                cx.advance();
                cy.advance();
            } else if (ix < iy) {
                cx.advance();
            } else {
                cy.advance();
            }
        }

        total += cx.remaining() + cy.remaining();
    }

    return total;
}

/// Computes the non-zero union of a compressed matrix and a trivial matrix in a
/// single pass using O(1) point probes.
fn compressedTrivialCount(comptime O: type, o_rows: usize, o_cols: usize, c: anytype, t: anytype) usize {
    const C = @TypeOf(c);

    const row_major = comptime meta.layoutOf(C) == .row_major;
    const lines: usize = if (row_major) o_rows else o_cols;
    const other_dim: usize = if (row_major) o_cols else o_rows;
    const c_diag = comptime matutils.hasImplicitDiag(C);

    var own: usize = c.nnz;
    if (c_diag) own += int.min(o_rows, o_cols);
    var inter: usize = 0;

    var line: usize = 0;
    while (line < lines) : (line += 1) {
        var p: usize = c.ptr[line];
        const end: usize = c.ptr[line + 1];
        while (p < end) : (p += 1) {
            const other = c.idx[p];
            const row = if (row_major) line else other;
            const col = if (row_major) other else line;
            if (matutils.hasNZ(t, row, col)) inter += 1;
            if ((comptime matutils.needsMirror(O, C)) and row != col) {
                own += 1;
                if (matutils.hasNZ(t, col, row)) inter += 1;
            }
        }

        if (c_diag and line < other_dim and matutils.hasNZ(t, line, line)) inter += 1;
    }

    return own + trivialCount(t, o_rows, o_cols) - inter;
}

/// Counts the overlapping non-zero coordinates between a trivial matrix `small`
/// and any arbitrary matrix `big`.
fn trivialIntersect(o_rows: usize, o_cols: usize, small: anytype, big: anytype) usize {
    const S = @TypeOf(small);
    var inter: usize = 0;
    switch (comptime meta.matrixType(S)) {
        .diagonal_static, .diagonal_sparse => {
            const n = int.min(o_rows, o_cols);
            var i: usize = 0;
            while (i < n) : (i += 1) {
                if (matutils.hasNZ(big, i, i)) inter += 1;
            }
        },
        .permutation_static, .permutation_sparse => {
            var i: usize = 0;
            while (i < o_rows) : (i += 1) {
                const row = if (S.direction == .forward) i else small.idx[i];
                const col = if (S.direction == .forward) small.idx[i] else i;
                if (matutils.hasNZ(big, row, col)) inter += 1;
            }
        },
        else => unreachable,
    }
    return inter;
}

/// Computes the non-zero union between two trivial (diagonal or permutation)
/// matrices.
fn trivialTrivialCount(o_rows: usize, o_cols: usize, a: anytype, b: anytype) usize {
    const ca = trivialCount(a, o_rows, o_cols);
    const cb = trivialCount(b, o_rows, o_cols);
    const inter = if (ca <= cb) trivialIntersect(o_rows, o_cols, a, b) else trivialIntersect(o_rows, o_cols, b, a);
    return ca + cb - inter;
}

/// Iterates over every valid mathematical coordinate in `m`, automatically
/// materializing symmetric mirrors and unit diagonals for `O`.
fn visitCompressedPositions(comptime O: type, o_rows: usize, o_cols: usize, m: anytype, comptime callback: anytype, ctx: anytype) void {
    const M = @TypeOf(m);

    const row_major = comptime meta.layoutOf(M) == .row_major;
    const lines: usize = if (row_major) o_rows else o_cols;
    const other_dim: usize = if (row_major) o_cols else o_rows;
    const mirrored = comptime matutils.needsMirror(O, M);
    const m_diag = comptime matutils.hasImplicitDiag(M);

    var line: usize = 0;
    while (line < lines) : (line += 1) {
        var p: usize = m.ptr[line];
        const end: usize = m.ptr[line + 1];
        while (p < end) : (p += 1) {
            const other = m.idx[p];
            const row = if (row_major) line else other;
            const col = if (row_major) other else line;
            callback(ctx, row, col);
            if (mirrored and row != col) callback(ctx, col, row);
        }

        if (m_diag and line < other_dim) callback(ctx, line, line);
    }
}

/// Counts overlapping coordinates between two compressed matrices by sweeping
/// `small` and binary-searching `big`.
fn intersectCompressed(comptime O: type, o_rows: usize, o_cols: usize, small: anytype, big: anytype) usize {
    const Ctx = struct { big: @TypeOf(big), n: usize };
    var ctx = Ctx{ .big = big, .n = 0 };
    visitCompressedPositions(O, o_rows, o_cols, small, struct {
        fn cb(ctxp: *Ctx, row: usize, col: usize) void {
            if (matutils.hasNZ(ctxp.big, row, col)) ctxp.n += 1;
        }
    }.cb, &ctx);
    return ctx.n;
}

/// Fallback non-zero union counter for incompatible sparse layouts via the
/// inclusion-exclusion principle: |X| + |Y| - |X ∩ Y|.
fn inclusionExclusionCount(comptime O: type, o_rows: usize, o_cols: usize, x: anytype, y: anytype) usize {
    const cx = ownCount(O, o_rows, o_cols, x);
    const cy = ownCount(O, o_rows, o_cols, y);
    const inter = if (cx <= cy) intersectCompressed(O, o_rows, o_cols, x, y) else intersectCompressed(O, o_rows, o_cols, y, x);
    return cx + cy - inter;
}

/// Resolves the exact non-zero capacity required for the elementwise operation
/// `x op y`.
fn apply2NNZ(comptime O: type, o_rows: usize, o_cols: usize, x: anytype, y: anytype) usize {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.isNumeric(X)) return singleOperandCount(O, o_rows, o_cols, y);
    if (comptime meta.isNumeric(Y)) return singleOperandCount(O, o_rows, o_cols, x);

    const x_compressed = comptime matutils.isCompressed(X);
    const y_compressed = comptime matutils.isCompressed(Y);

    if (comptime x_compressed and y_compressed) {
        if (comptime meta.layoutOf(X) == meta.layoutOf(Y) and !matutils.needsMirror(O, X) and !matutils.needsMirror(O, Y))
            return mergeCount(o_rows, o_cols, x, y);

        return inclusionExclusionCount(O, o_rows, o_cols, x, y);
    }

    if (comptime x_compressed) return compressedTrivialCount(O, o_rows, o_cols, x, y);
    if (comptime y_compressed) return compressedTrivialCount(O, o_rows, o_cols, y, x);

    return trivialTrivialCount(o_rows, o_cols, x, y);
}
