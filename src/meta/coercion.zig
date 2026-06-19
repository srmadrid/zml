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
