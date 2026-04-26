const std = @import("std");

const types = @import("../types.zig");

const int = @import("../int.zig");
const dyadic = @import("../dyadic.zig");
const Dyadic = dyadic.Dyadic;
const complex = @import("../complex.zig");
const Complex = complex.Complex;

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");
const expression = @import("../expression.zig");
const Expression = expression.Expression;

/// Coerces the input types to the smallest type that can represent both types.
///
/// This function takes two types `X` and `Y` and returns the smallest type that
/// can represent both types without loss of information.
///
/// For two ints, if they have different signedness, the result is a signed int.
/// The bit-width of the result is either the larger of the two bit-widths (if
/// the signed type is larger) or the larger of the two bit-widths plus one (if
/// the unsigned type is larger). If both ints are "standard" (see
/// `types.standard_integer_types`), the result is the next larger standard
/// type that can hold both values.
///
/// For two matrices, the coerced type uses the order of the denser operand.
/// Density is ranked (most to least) as:
///   `general.Dense`, `symmetric.Dense`/`hermitian.Dense`, `triangular.Dense`,
///   `general.Sparse`, `symmetric.Sparse`/`hermitian.Sparse`,
///   `triangular.Sparse`.
/// If both operands fall in the same rank but have different orders, the result
/// uses the left operand’s order. `Diagonal`, `Tridiagonal`, and `Permutation`
/// do not contribute order information. If the denser operand is one of these
/// (or both are), the result inherits the other operand’s order; if neither
/// provides an order and the resulting type requires one, the default is
/// `.col_major`.
///
/// If either `X` or `Y` is a custom type, it must implement the required
/// `Coerce` method. The expected signature and behavior of `Coerce` are as
/// follows:
/// * `fn Coerce(comptime X: type, comptime Y: type) type`: This function should
///   return the coerced type that can represent both `X` and `Y`.
///
/// ## Arguments
/// * `X` (`comptime type`): The first type to coerce. Must supported type.
/// * `Y` (`comptime type`): The second type to coerce. Must supported type.
/// * `Ctx` (`comptime type`): The type of a context parameter.
///
/// ## Returns
/// `type`: The coerced type that can represent both `X` and `Y`.
pub fn Coerce(comptime X: type, comptime Y: type) type {
    if (comptime X == Y and
        !types.isTriangularDenseMatrix(X) and
        !types.isPermutationMatrix(X))
        return X;

    if (comptime types.isCustomNumeric(X)) {
        if (comptime types.isCustomNumeric(Y)) {
            if (comptime types.hasMethod(X, "Coerce", fn (type, type) type, &.{}))
                return X.Coerce(X, Y);

            if (comptime types.hasMethod(Y, "Coerce", fn (type, type) type, &.{}))
                return Y.Coerce(X, Y);

            @compileError("zml.types.Coerce: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Coerce(type, type) type`");
        }

        if (comptime types.hasMethod(X, "Coerce", fn (type, type) type, &.{}))
            return X.Coerce(X, Y);

        @compileError("zml.types.Coerce: " ++ @typeName(X) ++ " must implement `fn Coerce(type, type) type`");
    } else if (comptime types.isCustomNumeric(Y)) {
        if (comptime types.hasMethod(Y, "Coerce", fn (type, type) type, &.{}))
            return Y.Coerce(X, Y);

        @compileError("zml.types.Coerce: " ++ @typeName(Y) ++ " must implement `fn Coerce(type, type) type`");
    }

    switch (comptime types.domain(X)) {
        .numeric => switch (comptime types.domain(Y)) {
            .numeric => {}, // Dealt with later
            .vector => switch (comptime types.vectorType(Y)) {
                .dense => return vector.Dense(Coerce(X, types.Numeric(Y))), // numeric + dense vector
                .sparse => return vector.Sparse(Coerce(X, types.Numeric(Y))), // numeric + sparse vector
                .custom => unreachable,
                .numeric => unreachable,
            },
            .matrix => switch (comptime types.matrixType(Y)) {
                .general_dense => return matrix.general.Dense(Coerce(X, types.Numeric(Y)), types.layoutOf(Y)), // numeric + general dense matrix
                .general_sparse => return matrix.general.Sparse(Coerce(X, types.Numeric(Y)), types.layoutOf(Y)), // numeric + general sparse matrix
                .symmetric_dense => return matrix.symmetric.Dense(Coerce(X, types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)), // numeric + symmetric dense matrix
                .symmetric_sparse => return matrix.symmetric.Sparse(Coerce(X, types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)), // numeric + symmetric sparse matrix
                .hermitian_dense => {
                    if (comptime types.isComplex(X))
                        return matrix.general.Dense(Coerce(X, types.Numeric(Y)), types.layoutOf(Y)) // numeric (complex) + hermitian dense matrix
                    else
                        return matrix.hermitian.Dense(Coerce(X, types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)); // numeric (real) + hermitian dense matrix

                },
                .hermitian_sparse => {
                    if (comptime types.isComplex(X))
                        return matrix.general.Sparse(Coerce(X, types.Numeric(Y)), types.layoutOf(Y)) // numeric (complex) + hermitian sparse matrix
                    else
                        return matrix.hermitian.Sparse(Coerce(X, types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)); // numeric (real) + hermitian sparse matrix

                },
                .triangular_dense => return matrix.triangular.Dense(Coerce(X, types.Numeric(Y)), types.uploOf(Y), .non_unit, types.layoutOf(Y)), // numeric + triangular dense matrix
                .triangular_sparse => return matrix.triangular.Sparse(Coerce(X, types.Numeric(Y)), types.uploOf(Y), .non_unit, types.layoutOf(Y)), // numeric + triangular sparse matrix
                .diagonal => return matrix.Diagonal(Coerce(X, types.Numeric(Y))), // numeric + diagonal matrix
                .permutation => return matrix.general.Sparse(Coerce(X, types.Numeric(Y)), types.layoutOf(Y)), // numeric + permutation matrix
                .custom => unreachable,
                .numeric => unreachable,
            },
            .array => switch (comptime types.arrayType(Y)) {
                .dense => return array.Dense(Coerce(X, types.Numeric(Y)), types.layoutOf(Y)), // numeric + dense array
                .strided => return array.Dense(Coerce(X, types.Numeric(Y))), // numeric + strided array
                .sparse => return array.Sparse(Coerce(X, types.Numeric(Y)), types.layoutOf(Y)), // numeric + sparse array
                .custom => unreachable,
                .numeric => unreachable,
            },
            .expression => Expression, // numeric + expression
        },
        .vector => switch (comptime types.vectorType(X)) {
            .dense => switch (comptime types.domain(Y)) {
                .numeric => return vector.Dense(Coerce(types.Numeric(X), Y)), // dense vector + numeric
                .vector => switch (comptime types.vectorType(Y)) {
                    .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector + dense vector
                    .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector + sparse vector
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .matrix => @compileError("Cannot coerce vector and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense vector + matrix
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense vector + array
                .expression => Expression, // dense vector + expression
            },
            .sparse => switch (comptime types.domain(Y)) {
                .numeric => return vector.Sparse(Coerce(types.Numeric(X), Y)), // sparse vector + numeric
                .vector => switch (comptime types.vectorType(Y)) {
                    .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector + dense vector
                    .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector + sparse vector
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .matrix => @compileError("Cannot coerce vector and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector + matrix
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector + array
                .expression => Expression, // sparse vector + expression
            },
            .custom => unreachable,
            .numeric => unreachable,
        },
        .matrix => switch (comptime types.matrixType(X)) {
            .general_dense => switch (comptime types.domain(Y)) {
                .numeric => return matrix.general.Dense(Coerce(types.Numeric(X), Y), types.layoutOf(X)), // general dense matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general dense matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + general dense matrix
                    .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + triangular sparse matrix
                    .diagonal => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + diagonal matrix
                    .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general dense matrix + array
                .expression => Expression, // general dense matrix + expression
            },
            .general_sparse => switch (comptime types.domain(Y)) {
                .numeric => return matrix.general.Sparse(Coerce(types.Numeric(X), Y), types.layoutOf(X)), // general sparse matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general sparse matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // general sparse matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // general sparse matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // general sparse matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // general sparse matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix + triangular sparse matrix
                    .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general sparse matrix + array
                .expression => Expression, // general sparse matrix + expression
            },
            .symmetric_dense => switch (comptime types.domain(Y)) {
                .numeric => return matrix.symmetric.Dense(Coerce(types.Numeric(X), Y), types.uploOf(X), types.layoutOf(X)), // symmetric dense matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric dense matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric dense matrix + general dense matrix
                    .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix + general sparse matrix
                    .symmetric_dense => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)), // symmetric dense matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)), // symmetric dense matrix + symmetric sparse matrix
                    .hermitian_dense => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)) // symmetric dense matrix (complex) + hermitian dense matrix
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)); // symmetric dense matrix (real) + hermitian dense matrix
                    },
                    .hermitian_sparse => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)) // symmetric dense matrix (complex) + hermitian sparse matrix
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)); // symmetric dense matrix (real) + hermitian sparse matrix
                    },
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix + triangular sparse matrix
                    .diagonal => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)), // symmetric dense matrix + diagonal matrix
                    .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric dense matrix + array
                .expression => Expression, // symmetric dense matrix + expression
            },
            .symmetric_sparse => switch (comptime types.domain(Y)) {
                .numeric => return matrix.symmetric.Sparse(Coerce(types.Numeric(X), Y), types.uploOf(X), types.layoutOf(X)), // symmetric sparse matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric sparse matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric sparse matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric sparse matrix + general sparse matrix
                    .symmetric_dense => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)), // symmetric sparse matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.symmetric.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)), // symmetric sparse matrix + symmetric sparse matrix
                    .hermitian_dense => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)) // symmetric sparse matrix (complex) + hermitian dense matrix
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)); // symmetric sparse matrix (real) + hermitian dense matrix

                    },
                    .hermitian_sparse => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)) // symmetric sparse matrix (complex) + hermitian sparse matrix
                        else
                            return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)); // symmetric sparse matrix (real) + hermitian sparse matrix
                    },
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric sparse matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric sparse matrix + triangular sparse matrix
                    .diagonal => return matrix.symmetric.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)), // symmetric sparse matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric sparse matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric sparse matrix + array
                .expression => Expression, // symmetric sparse matrix + expression
            },
            .hermitian_dense => switch (comptime types.domain(Y)) {
                .numeric => {
                    if (comptime types.isComplex(Y))
                        return matrix.general.Dense(Coerce(types.Numeric(X), Y), types.layoutOf(X)) // hermitian dense matrix + numeric (complex)
                    else
                        return matrix.hermitian.Dense(Coerce(types.Numeric(X), Y), types.uploOf(X), types.layoutOf(X)); // hermitian dense matrix + numeric (real)
                },
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian dense matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian dense matrix + general dense matrix
                    .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix + general sparse matrix
                    .symmetric_dense => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)) // hermitian dense matrix + symmetric dense matrix (complex)
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)); // hermitian dense matrix + symmetric dense matrix (real)
                    },
                    .symmetric_sparse => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)) // hermitian dense matrix + symmetric sparse matrix (complex)
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)); // hermitian dense matrix + symmetric sparse matrix (real)
                    },
                    .hermitian_dense => return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)), // hermitian dense matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)), // hermitian dense matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix + triangular sparse matrix
                    .diagonal => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)) // hermitian dense matrix + diagonal matrix (complex)
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)); // hermitian dense matrix + diagonal matrix (real)
                    },
                    .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian dense matrix + array
                .expression => Expression, // hermitian dense matrix + expression
            },
            .hermitian_sparse => switch (comptime types.domain(Y)) {
                .numeric => {
                    if (comptime types.isComplex(Y))
                        return matrix.general.Sparse(Coerce(types.Numeric(X), Y), types.layoutOf(X)) // hermitian sparse matrix + numeric (complex)
                    else
                        return matrix.hermitian.Sparse(Coerce(types.Numeric(X), Y), types.uploOf(X), types.layoutOf(X)); // hermitian sparse matrix + numeric (real)
                },
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian sparse matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian sparse matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian sparse matrix + general sparse matrix
                    .symmetric_dense => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)) // hermitian sparse matrix + symmetric dense matrix (complex)
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)); // hermitian sparse matrix + symmetric dense matrix (real)
                    },
                    .symmetric_sparse => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)) // hermitian sparse matrix + symmetric sparse matrix (complex)
                        else
                            return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)); // hermitian sparse matrix + symmetric sparse matrix (real)
                    },
                    .hermitian_dense => return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)), // hermitian sparse matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)), // hermitian sparse matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian sparse matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(X)), types.layoutOf(X)), // hermitian sparse matrix + triangular sparse matrix
                    .diagonal => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)) // hermitian sparse matrix + diagonal matrix (complex)
                        else
                            return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.layoutOf(X)); // hermitian sparse matrix + diagonal matrix (real)
                    },
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian sparse matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian sparse matrix + array
                .expression => Expression, // hermitian sparse matrix + expression
            },
            .triangular_dense => switch (comptime types.domain(Y)) {
                .numeric => return matrix.triangular.Dense(Coerce(types.Numeric(X), Y), types.uploOf(X), .non_unit, types.layoutOf(X)), // triangular dense matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular dense matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular dense matrix + general dense matrix
                    .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular dense matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular dense matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular dense matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular dense matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular dense matrix + hermitian sparse matrix
                    .triangular_dense => {
                        if (comptime types.uploOf(X) == types.uploOf(Y))
                            return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)) // triangular dense matrix + triangular dense matrix (same uplo)
                        else
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)); // triangular dense matrix + triangular dense matrix (different uplo)
                    },
                    .triangular_sparse => {
                        if (comptime types.uploOf(X) == types.uploOf(Y))
                            return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)) // triangular dense matrix + triangular sparse matrix (same uplo)
                        else
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)); // triangular dense matrix + triangular sparse matrix (different uplo)
                    },
                    .diagonal => return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)), // triangular dense matrix + diagonal matrix
                    .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular dense matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular dense matrix + array
                .expression => Expression, // triangular dense matrix + expression
            },
            .triangular_sparse => switch (comptime types.domain(Y)) {
                .numeric => return matrix.triangular.Sparse(Coerce(types.Numeric(X), Y), types.uploOf(X), .non_unit, types.layoutOf(X)), // triangular sparse matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular sparse matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix + hermitian sparse matrix
                    .triangular_dense => {
                        if (comptime types.uploOf(X) == types.uploOf(Y))
                            return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(Y)) // triangular sparse matrix + triangular dense matrix (same uplo)
                        else
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)); // triangular sparse matrix + triangular dense matrix (different uplo)
                    },
                    .triangular_sparse => {
                        if (comptime types.uploOf(X) == types.uploOf(Y))
                            return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)) // triangular sparse matrix + triangular sparse matrix (same uplo)
                        else
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)); // triangular sparse matrix + triangular sparse matrix (different uplo)
                    },
                    .diagonal => return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)), // triangular sparse matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular sparse matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular sparse matrix + array
                .expression => Expression, // triangular sparse matrix + expression
            },
            .diagonal => switch (comptime types.domain(Y)) {
                .numeric => return matrix.Diagonal(Coerce(types.Numeric(X), Y)), // diagonal matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // diagonal matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // diagonal matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // diagonal matrix + general sparse matrix
                    .symmetric_dense => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)), // diagonal matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.symmetric.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)), // diagonal matrix + symmetric sparse matrix
                    .hermitian_dense => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)) // diagonal matrix (complex) + hermitian dense matrix
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)); // diagonal matrix (real) + hermitian dense matrix
                    },
                    .hermitian_sparse => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)) // diagonal matrix (complex) + hermitian sparse matrix
                        else
                            return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.layoutOf(Y)); // diagonal matrix (real) + hermitian sparse matrix
                    },
                    .triangular_dense => return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), .non_unit, types.layoutOf(Y)), // diagonal matrix + triangular dense matrix
                    .triangular_sparse => return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), .non_unit, types.layoutOf(Y)), // diagonal matrix + triangular sparse matrix
                    .diagonal => return matrix.Diagonal(Coerce(types.Numeric(X), types.Numeric(Y))), // diagonal matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // diagonal matrix + array
                .expression => Expression, // diagonal matrix + expression
            },
            .permutation => switch (comptime types.domain(Y)) {
                .numeric => return matrix.general.Sparse(Coerce(types.Numeric(X), Y), types.layoutOf(X)), // permutation matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // permutation matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // permutation matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // permutation matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // permutation matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // permutation matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // permutation matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // permutation matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // permutation matrix + triangular sparse matrix
                    .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix + permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix + array
                .expression => Expression, // permutation matrix + expression
            },
            .custom => unreachable,
            .numeric => unreachable,
        },
        .array => switch (comptime types.arrayType(X)) {
            .dense => switch (comptime types.domain(Y)) {
                .numeric => return array.Dense(Coerce(types.Numeric(X), Y), types.layoutOf(X)), // dense + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense + matrix
                .array => switch (comptime types.arrayType(Y)) {
                    .dense => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // dense + dense
                    .strided => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // dense + strided
                    .sparse => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // dense + sparse
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .expression => Expression, // dense + expression
            },
            .strided => switch (comptime types.domain(Y)) {
                .numeric => return array.Dense(Coerce(types.Numeric(X), Y), types.layoutOf(X)), // strided + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // strided + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // strided + matrix
                .array => switch (comptime types.arrayType(Y)) {
                    .dense => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // strided + dense
                    .strided => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // strided + strided
                    .sparse => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // strided + sparse
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .expression => Expression, // strided + expression
            },
            .sparse => switch (comptime types.domain(Y)) {
                .numeric => return array.Sparse(Coerce(types.Numeric(X), Y), types.layoutOf(X)), // sparse + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse + matrix
                .array => switch (comptime types.arrayType(Y)) {
                    .dense => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // sparse + dense
                    .strided => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // sparse + strided
                    .sparse => return array.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // sparse + sparse
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .expression => Expression, // sparse + expression
            },
            .custom => unreachable,
            .numeric => unreachable,
        },
        .expression => Expression, // expression + anything
    }
}

/// Coerces the input types to the smallest type that can represent the result
/// of their multiplication.
///
/// For scalar or array types, this is equivalent to `Coerce`. For matrices and
/// vectors, this function takes into account the rules of linear algebra to
/// determine the appropriate resulting type.
///
/// If either `X` or `Y` is a custom type, it must implement the required
/// `MulCoerce` method. The expected signature and behavior of `MulCoerce` are
/// as follows:
/// * `fn Coerce(comptime X: type, comptime Y: type) type`: This function should
///   return the coerced type that can represent the result of multiplying `X`
///   and `Y`.
///
/// ## Arguments
/// * `X` (`comptime type`): The first type to coerce. Must be a supported type.
/// * `Y` (`comptime type`): The second type to coerce. Must be a supported
///   type.
///
/// ## Returns
/// `type`: The coerced type that can represent the result of multiplying `X`
/// and `Y`.
pub fn MulCoerce(comptime X: type, comptime Y: type) type {
    if (comptime types.isCustomNumeric(X)) {
        if (comptime types.isCustomNumeric(Y)) {
            if (comptime types.hasMethod(X, "MulCoerce", fn (type, type) type, &.{}))
                return X.MulCoerce(X, Y);

            if (comptime types.hasMethod(Y, "MulCoerce", fn (type, type) type, &.{}))
                return Y.MulCoerce(X, Y);

            @compileError("zml.types.MulCoerce: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn MulCoerce(type, type) type`");
        }

        if (comptime types.hasMethod(X, "MulCoerce", fn (type, type) type, &.{}))
            return X.MulCoerce(X, Y);

        @compileError("zml.types.MulCoerce: " ++ @typeName(X) ++ " must implement `fn MulCoerce(type, type) type`");
    } else if (comptime types.isCustomNumeric(Y)) {
        if (comptime types.hasMethod(Y, "MulCoerce", fn (type, type) type, &.{}))
            return Y.MulCoerce(X, Y);

        @compileError("zml.types.MulCoerce: " ++ @typeName(Y) ++ " must implement `fn MulCoerce(type, type) type`");
    }

    switch (comptime types.domain(X)) {
        .numeric => return Coerce(X, Y), // Same as Coerce
        .vector => switch (comptime types.vectorType(X)) {
            .dense => switch (comptime types.domain(Y)) {
                .numeric => return Coerce(X, Y), // Same as Coerce
                .vector => switch (comptime types.vectorType(Y)) {
                    .dense => return Coerce(types.Numeric(X), types.Numeric(Y)), // dense vector * dense vector
                    .sparse => return Coerce(types.Numeric(X), types.Numeric(Y)), // dense vector * sparse vector
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * general dense matrix
                    .general_sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * general sparse matrix
                    .symmetric_dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * symmetric dense matrix
                    .symmetric_sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * symmetric sparse matrix
                    .hermitian_dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * hermitian dense matrix
                    .hermitian_sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * hermitian sparse matrix
                    .triangular_dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * triangular dense matrix
                    .triangular_sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * triangular sparse matrix
                    .diagonal => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * diagonal matrix
                    .permutation => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // dense vector * permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense vector * array
                .expression => Expression, // dense vector * expression
            },
            .sparse => switch (comptime types.domain(Y)) {
                .numeric => return Coerce(X, Y), // Same as Coerce
                .vector => switch (comptime types.vectorType(Y)) {
                    .dense => return Coerce(types.Numeric(X), types.Numeric(Y)), // sparse vector * dense vector
                    .sparse => return Coerce(types.Numeric(X), types.Numeric(Y)), // sparse vector * sparse vector
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * general dense matrix
                    .general_sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * general sparse matrix
                    .symmetric_dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * symmetric dense matrix
                    .symmetric_sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * symmetric sparse matrix
                    .hermitian_dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * hermitian dense matrix
                    .hermitian_sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * hermitian sparse matrix
                    .triangular_dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * triangular dense matrix
                    .triangular_sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * triangular sparse matrix
                    .diagonal => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * diagonal matrix
                    .permutation => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // sparse vector * permutation matrix
                    .custom => unreachable,
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector * array
                .expression => Expression, // sparse vector * expression
            },
            .custom => unreachable,
            .numeric => unreachable,
        },
        .matrix => {
            switch (comptime types.matrixType(X)) {
                .general_dense => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // general dense matrix * dense vector
                        .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // general dense matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * general dense matrix
                        .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * diagonal matrix
                        .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general dense matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general dense matrix * array
                    .expression => Expression, // general dense matrix * expression
                },
                .general_sparse => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // general sparse matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // general sparse matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // general sparse matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // general sparse matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // general sparse matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // general sparse matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // general sparse matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general sparse matrix * array
                    .expression => Expression, // general sparse matrix * expression
                },
                .symmetric_dense => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // symmetric dense matrix * dense vector
                        .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // symmetric dense matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric dense matrix * general dense matrix
                        .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * diagonal matrix
                        .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric dense matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric dense matrix * array
                    .expression => Expression, // symmetric dense matrix * expression
                },
                .symmetric_sparse => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // symmetric sparse matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // symmetric sparse matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric sparse matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric sparse matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric sparse matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric sparse matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric sparse matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric sparse matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // symmetric sparse matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric sparse matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric sparse matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // symmetric sparse matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric sparse matrix * array
                    .expression => Expression, // symmetric sparse matrix * expression
                },
                .hermitian_dense => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // hermitian dense matrix * dense vector
                        .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // hermitian dense matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian dense matrix * general dense matrix
                        .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * diagonal matrix
                        .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian dense matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian dense matrix * array
                    .expression => Expression, // hermitian dense matrix * expression
                },
                .hermitian_sparse => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // hermitian sparse matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // hermitian sparse matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian sparse matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian sparse matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian sparse matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian sparse matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian sparse matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian sparse matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // hermitian sparse matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian sparse matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian sparse matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // hermitian sparse matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian sparse matrix * array
                    .expression => Expression, // hermitian sparse matrix * expression
                },
                .triangular_dense => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // triangular dense matrix * dense vector
                        .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // triangular dense matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular dense matrix * general dense matrix
                        .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular dense matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular dense matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular dense matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular dense matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular dense matrix * hermitian sparse matrix
                        .triangular_dense => {
                            if (comptime types.uploOf(X) == types.uploOf(Y)) {
                                if (comptime types.diagOf(X) == types.diagOf(Y))
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.diagOf(X), types.layoutOf(X)) // triangular dense matrix * triangular dense matrix (same uplo and diag)
                                else
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)); // triangular dense matrix * triangular dense matrix (same uplo, different diag)
                            } else {
                                return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)); // triangular dense matrix * triangular dense matrix (different uplo)
                            }
                        },
                        .triangular_sparse => {
                            if (comptime types.uploOf(X) == types.uploOf(Y)) {
                                if (comptime types.diagOf(X) == types.diagOf(Y))
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.diagOf(X), types.layoutOf(X)) // triangular dense matrix * triangular sparse matrix (same uplo and diag)
                                else
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)); // triangular dense matrix * triangular sparse matrix (same uplo, different diag)
                            } else {
                                return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)); // triangular dense matrix * triangular sparse matrix (different uplo)
                            }
                        },
                        .diagonal => return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)), // triangular dense matrix * diagonal matrix
                        .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular dense matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular dense matrix * array
                    .expression => Expression, // triangular dense matrix * expression
                },
                .triangular_sparse => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // triangular sparse matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // triangular sparse matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)), // triangular sparse matrix * hermitian sparse matrix
                        .triangular_dense => {
                            if (comptime types.uploOf(X) == types.uploOf(Y)) {
                                if (comptime types.diagOf(X) == types.diagOf(Y))
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.diagOf(X), types.layoutOf(Y)) // triangular sparse matrix * triangular dense matrix (same uplo and diag)
                                else
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(Y)); // triangular sparse matrix * triangular dense matrix (same uplo, different diag)
                            } else {
                                return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(Y)); // triangular sparse matrix * triangular dense matrix (different uplo)
                            }
                        },
                        .triangular_sparse => {
                            if (comptime types.uploOf(X) == types.uploOf(Y)) {
                                if (comptime types.diagOf(X) == types.diagOf(Y))
                                    return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.diagOf(X), types.layoutOf(X)) // triangular sparse matrix * triangular sparse matrix (same uplo and diag)
                                else
                                    return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)); // triangular sparse matrix * triangular sparse matrix (same uplo, different diag)
                            } else {
                                return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)); // triangular sparse matrix * triangular sparse matrix (different uplo)
                            }
                        },
                        .diagonal => return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.layoutOf(X)), // triangular sparse matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // triangular sparse matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular sparse matrix * array
                    .expression => Expression, // triangular sparse matrix * expression
                },
                .diagonal => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // diagonal matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // diagonal matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), .non_unit, types.layoutOf(X)), // diagonal matrix * triangular dense matrix
                        .triangular_sparse => return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), .non_unit, types.layoutOf(X)), // diagonal matrix * triangular sparse matrix
                        .diagonal => return matrix.Diagonal(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // diagonal matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // diagonal matrix * array
                    .expression => Expression, // diagonal matrix * expression
                },
                .permutation => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // permutation matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // permutation matrix * sparse vector
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * diagonal matrix
                        .permutation => return matrix.Permutation(Coerce(types.Numeric(X), types.Numeric(Y)), types.layoutOf(X)), // permutation matrix * permutation matrix
                        .custom => unreachable,
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix * array
                    .expression => Expression, // permutation matrix * expression
                },
                .custom => unreachable,
                .numeric => unreachable,
            }
        },
        .array => {
            switch (comptime types.domain(X)) {
                .numeric => return Coerce(X, Y), // array * numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // array * vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // array * matrix
                .array => return Coerce(X, Y), // array * array
                .expression => Expression,
            }
        },
        .expression => Expression,
    }

    return Coerce(X, Y);
}

/// Checks if if `F` can be coerced to `T` without loss of information. This is
/// a more flexible version of `Coerce`, as it does not require `T` to be to the
/// smallest type that can represent both types. The only requirement is that
/// `T` can represent all values of the first two types.
///
/// ## Arguments
/// * `F` (`comptime type`): The type to check if it can be coerced. Must be a
///   supported numeric type.
/// * `T` (`comptime type`): The target type. Must be a supported numeric type.
///
/// ## Returns
/// `bool`: `true` if `F` can be coerced to `T` without loss of information,
/// `false` otherwise.
pub fn canCoerce(comptime F: type, comptime T: type) bool {
    comptime if (!types.isNumeric(F))
        @compileError("zml.types.canCoerce: " ++ @typeName(F) ++ " is not a supported numeric type");

    comptime if (!types.isNumeric(T))
        @compileError("zml.types.canCoerce: " ++ @typeName(T) ++ " is not a supported numeric type");

    return Coerce(F, T) == T;
}

/// Coerces the second type to a vector, matrix, array or expression type based
/// on the first type.
///
/// This function is useful for ensuring that the second type is always a
/// vector, matrix or an array when the first is.
///
/// ## Arguments
/// * `X` (`comptime type`): The type to check. Must be a supported type.
/// * `Y` (`comptime type`): The type to coerce. Must be a numeric type.
///
/// Returns
/// -------
/// `type`: The coerced type.
pub fn EnsureDomain(comptime X: type, comptime Y: type) type {
    if (comptime !types.isNumeric(Y))
        @compileError("zml.types.EnsureDomain: " ++ @typeName(Y) ++ " is not a supported numeric type");

    switch (comptime types.domain(X)) {
        .numeric => return Y,
        .vector => return EnsureVector(X, Y),
        .matrix => return EnsureMatrix(X, Y),
        .array => return EnsureArray(X, Y),
        .expression => return Expression,
    }
}

pub fn EnsureMatrix(comptime X: type, comptime Y: type) type {
    if (comptime !types.isNumeric(Y))
        @compileError("zml.types.EnsureMatrix: " ++ @typeName(Y) ++ " is not a supported numeric type");

    switch (comptime types.domain(X)) {
        .numeric => return Y,
        .vector => return Y,
        .matrix => switch (types.matrixType(X)) {
            .general_dense => return matrix.general.Dense(Y, types.layoutOf(X)),
            .general_sparse => return matrix.general.Sparse(Y, types.layoutOf(X)),
            .symmetric_dense => return matrix.symmetric.Dense(Y, types.uploOf(X), types.layoutOf(X)),
            .symmetric_sparse => return matrix.symmetric.Sparse(Y, types.uploOf(X), types.layoutOf(X)),
            .hermitian_dense => return matrix.hermitian.Dense(Y, types.uploOf(X), types.layoutOf(X)),
            .hermitian_sparse => return matrix.hermitian.Sparse(Y, types.uploOf(X), types.layoutOf(X)),
            .triangular_dense => return matrix.triangular.Dense(Y, types.uploOf(X), types.diagOf(X), types.layoutOf(X)),
            .triangular_sparse => return matrix.triangular.Sparse(Y, types.uploOf(X), types.diagOf(X), types.layoutOf(X)),
            .diagonal => return matrix.Diagonal(Y),
            .permutation => return matrix.Permutation(Y),
            .custom => {
                if (comptime types.hasMethod(types.matrixType(X), "EnsureMatrix", fn (type, type) type, &.{}))
                    return types.matrixType(X).EnsureMatrix(X, Y);

                @compileError("zml.types.EnsureMatrix: " ++ @typeName(types.matrixType(X)) ++ " must implement `fn EnsureMatrix(type, type) type`");
            },
            .numeric => unreachable,
        },
        .array => return Y,
        .expression => return Y,
    }
}

pub fn EnsureArray(comptime X: type, comptime Y: type) type {
    if (comptime !types.isNumeric(Y))
        @compileError("zml.types.EnsureArray: " ++ @typeName(Y) ++ " is not a supported numeric type");

    switch (comptime types.domain(X)) {
        .numeric => return Y,
        .vector => return Y,
        .matrix => return Y,
        .array => switch (types.arrayType(X)) {
            .dense => return array.Dense(Y, types.layoutOf(X)),
            .strided => return array.Dense(Y, types.layoutOf(X)),
            .sparse => return array.Sparse(Y, types.layoutOf(X)),
            .custom => {
                if (comptime types.hasMethod(types.arrayType(X), "EnsureArray", fn (type, type) type, &.{}))
                    return types.arrayType(X).EnsureArray(X, Y);

                @compileError("zml.types.EnsureArray: " ++ @typeName(types.arrayType(X)) ++ " must implement `fn EnsureArray(type, type) type`");
            },
            .numeric => unreachable,
        },
        .expression => return Y,
    }
}
