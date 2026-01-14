const std = @import("std");

const types = @import("../types.zig");

const int = @import("../int.zig");
const dyadic = @import("../dyadic.zig");
const Dyadic = dyadic.Dyadic;
const cfloat = @import("../cfloat.zig");
const Cfloat = cfloat.Cfloat;
const integer = @import("../integer.zig");
const Integer = integer.Integer;
const rational = @import("../rational.zig");
const Rational = rational.Rational;
const real = @import("../real.zig");
const Real = real.Real;
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
/// `types.standard_integer_types`), the result is the next larger standard type
/// if needed.
///
/// For two matrices, the coerced type will use the order of the denser operand.
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
/// Parameters
/// ----------
/// comptime X (`type`): The first type to coerce. Must be a supported numeric
/// type, a vector, a matrix, or an array.
///
/// comptime Y (`type`): The second type to coerce. Must be a supported numeric
/// type, a vector, a matrix, or an array.
///
/// Returns
/// -------
/// `type`: The coerced type that can represent both `X` and `Y`.
pub fn Coerce(comptime X: type, comptime Y: type) type {
    if (comptime X == Y and !types.isTriangularDenseMatrix(X) and !types.isPermutationMatrix(X))
        return X;

    switch (comptime types.domain(X)) {
        .numeric => switch (comptime types.domain(Y)) {
            .numeric => {},
            .vector => switch (comptime types.vectorType(Y)) {
                .dense => return vector.Dense(Coerce(X, types.Numeric(Y))), // numeric + dense vector
                .sparse => return vector.Sparse(Coerce(X, types.Numeric(Y))), // numeric + sparse vector
                .numeric => unreachable,
            },
            .matrix => switch (comptime types.matrixType(Y)) {
                .general_dense => return matrix.general.Dense(Coerce(X, types.Numeric(Y)), types.orderOf(Y)), // numeric + general dense matrix
                .general_sparse => return matrix.general.Sparse(Coerce(X, types.Numeric(Y)), types.orderOf(Y)), // numeric + general sparse matrix
                .symmetric_dense => return matrix.symmetric.Dense(Coerce(X, types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)), // numeric + symmetric dense matrix
                .symmetric_sparse => return matrix.symmetric.Sparse(Coerce(X, types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)), // numeric + symmetric sparse matrix
                .hermitian_dense => {
                    if (comptime types.isComplex(X))
                        return matrix.general.Dense(Coerce(X, types.Numeric(Y)), types.orderOf(Y)) // numeric (complex) + hermitian dense matrix
                    else
                        return matrix.hermitian.Dense(Coerce(X, types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)); // numeric (real) + hermitian dense matrix

                },
                .hermitian_sparse => {
                    if (comptime types.isComplex(X))
                        return matrix.general.Sparse(Coerce(X, types.Numeric(Y)), types.orderOf(Y)) // numeric (complex) + hermitian sparse matrix
                    else
                        return matrix.hermitian.Sparse(Coerce(X, types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)); // numeric (real) + hermitian sparse matrix

                },
                .triangular_dense => return matrix.triangular.Dense(Coerce(X, types.Numeric(Y)), types.uploOf(Y), .non_unit, types.orderOf(Y)), // numeric + triangular dense matrix
                .triangular_sparse => return matrix.triangular.Sparse(Coerce(X, types.Numeric(Y)), types.uploOf(Y), .non_unit, types.orderOf(Y)), // numeric + triangular sparse matrix
                .diagonal => return matrix.Diagonal(Coerce(X, types.Numeric(Y))), // numeric + diagonal matrix
                .permutation => return matrix.general.Sparse(Coerce(X, types.Numeric(Y)), types.orderOf(Y)), // numeric + permutation matrix
                .numeric => unreachable,
            },
            .array => switch (comptime types.arrayType(Y)) {
                .dense => return array.Dense(Coerce(X, types.Numeric(Y)), types.orderOf(Y)), // numeric + dense array
                .strided => return array.Dense(Coerce(X, types.Numeric(Y))), // numeric + strided array
                .sparse => return array.Sparse(Coerce(X, types.Numeric(Y)), types.orderOf(Y)), // numeric + sparse array
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
                    .numeric => unreachable,
                },
                .matrix => @compileError("Cannot coerce vector and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector + matrix
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector + array
                .expression => Expression, // sparse vector + expression
            },
            .numeric => unreachable,
        },
        .matrix => switch (comptime types.matrixType(X)) {
            .general_dense => switch (comptime types.domain(Y)) {
                .numeric => return matrix.general.Dense(Coerce(types.Numeric(X), Y), types.orderOf(X)), // general dense matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general dense matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + general dense matrix
                    .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + triangular sparse matrix
                    .diagonal => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + diagonal matrix
                    .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general dense matrix + array
                .expression => Expression, // general dense matrix + expression
            },
            .general_sparse => switch (comptime types.domain(Y)) {
                .numeric => return matrix.general.Sparse(Coerce(types.Numeric(X), Y), types.orderOf(X)), // general sparse matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general sparse matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // general sparse matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // general sparse matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // general sparse matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // general sparse matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix + triangular sparse matrix
                    .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general sparse matrix + array
                .expression => Expression, // general sparse matrix + expression
            },
            .symmetric_dense => switch (comptime types.domain(Y)) {
                .numeric => return matrix.symmetric.Dense(Coerce(types.Numeric(X), Y), types.uploOf(X), types.orderOf(X)), // symmetric dense matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric dense matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric dense matrix + general dense matrix
                    .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix + general sparse matrix
                    .symmetric_dense => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)), // symmetric dense matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)), // symmetric dense matrix + symmetric sparse matrix
                    .hermitian_dense => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)) // symmetric dense matrix (complex) + hermitian dense matrix
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)); // symmetric dense matrix (real) + hermitian dense matrix
                    },
                    .hermitian_sparse => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)) // symmetric dense matrix (complex) + hermitian sparse matrix
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)); // symmetric dense matrix (real) + hermitian sparse matrix
                    },
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix + triangular sparse matrix
                    .diagonal => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)), // symmetric dense matrix + diagonal matrix
                    .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric dense matrix + array
                .expression => Expression, // symmetric dense matrix + expression
            },
            .symmetric_sparse => switch (comptime types.domain(Y)) {
                .numeric => return matrix.symmetric.Sparse(Coerce(types.Numeric(X), Y), types.uploOf(X), types.orderOf(X)), // symmetric sparse matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric sparse matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric sparse matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric sparse matrix + general sparse matrix
                    .symmetric_dense => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)), // symmetric sparse matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.symmetric.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)), // symmetric sparse matrix + symmetric sparse matrix
                    .hermitian_dense => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)) // symmetric sparse matrix (complex) + hermitian dense matrix
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)); // symmetric sparse matrix (real) + hermitian dense matrix

                    },
                    .hermitian_sparse => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)) // symmetric sparse matrix (complex) + hermitian sparse matrix
                        else
                            return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)); // symmetric sparse matrix (real) + hermitian sparse matrix
                    },
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric sparse matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric sparse matrix + triangular sparse matrix
                    .diagonal => return matrix.symmetric.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)), // symmetric sparse matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric sparse matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric sparse matrix + array
                .expression => Expression, // symmetric sparse matrix + expression
            },
            .hermitian_dense => switch (comptime types.domain(Y)) {
                .numeric => {
                    if (comptime types.isComplex(Y))
                        return matrix.general.Dense(Coerce(types.Numeric(X), Y), types.orderOf(X)) // hermitian dense matrix + numeric (complex)
                    else
                        return matrix.hermitian.Dense(Coerce(types.Numeric(X), Y), types.uploOf(X), types.orderOf(X)); // hermitian dense matrix + numeric (real)
                },
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian dense matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian dense matrix + general dense matrix
                    .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix + general sparse matrix
                    .symmetric_dense => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)) // hermitian dense matrix + symmetric dense matrix (complex)
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)); // hermitian dense matrix + symmetric dense matrix (real)
                    },
                    .symmetric_sparse => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)) // hermitian dense matrix + symmetric sparse matrix (complex)
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)); // hermitian dense matrix + symmetric sparse matrix (real)
                    },
                    .hermitian_dense => return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)), // hermitian dense matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)), // hermitian dense matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix + triangular sparse matrix
                    .diagonal => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)) // hermitian dense matrix + diagonal matrix (complex)
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)); // hermitian dense matrix + diagonal matrix (real)
                    },
                    .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian dense matrix + array
                .expression => Expression, // hermitian dense matrix + expression
            },
            .hermitian_sparse => switch (comptime types.domain(Y)) {
                .numeric => {
                    if (comptime types.isComplex(Y))
                        return matrix.general.Sparse(Coerce(types.Numeric(X), Y), types.orderOf(X)) // hermitian sparse matrix + numeric (complex)
                    else
                        return matrix.hermitian.Sparse(Coerce(types.Numeric(X), Y), types.uploOf(X), types.orderOf(X)); // hermitian sparse matrix + numeric (real)
                },
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian sparse matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian sparse matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian sparse matrix + general sparse matrix
                    .symmetric_dense => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)) // hermitian sparse matrix + symmetric dense matrix (complex)
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)); // hermitian sparse matrix + symmetric dense matrix (real)
                    },
                    .symmetric_sparse => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)) // hermitian sparse matrix + symmetric sparse matrix (complex)
                        else
                            return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)); // hermitian sparse matrix + symmetric sparse matrix (real)
                    },
                    .hermitian_dense => return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)), // hermitian sparse matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)), // hermitian sparse matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian sparse matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(X)), types.orderOf(X)), // hermitian sparse matrix + triangular sparse matrix
                    .diagonal => {
                        if (comptime types.isComplex(types.Numeric(Y)))
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)) // hermitian sparse matrix + diagonal matrix (complex)
                        else
                            return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.orderOf(X)); // hermitian sparse matrix + diagonal matrix (real)
                    },
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian sparse matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian sparse matrix + array
                .expression => Expression, // hermitian sparse matrix + expression
            },
            .triangular_dense => switch (comptime types.domain(Y)) {
                .numeric => return matrix.triangular.Dense(Coerce(types.Numeric(X), Y), types.uploOf(X), .non_unit, types.orderOf(X)), // triangular dense matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular dense matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular dense matrix + general dense matrix
                    .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular dense matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular dense matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular dense matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular dense matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular dense matrix + hermitian sparse matrix
                    .triangular_dense => {
                        if (comptime types.uploOf(X) == types.uploOf(Y))
                            return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)) // triangular dense matrix + triangular dense matrix (same uplo)
                        else
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)); // triangular dense matrix + triangular dense matrix (different uplo)
                    },
                    .triangular_sparse => {
                        if (comptime types.uploOf(X) == types.uploOf(Y))
                            return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)) // triangular dense matrix + triangular sparse matrix (same uplo)
                        else
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)); // triangular dense matrix + triangular sparse matrix (different uplo)
                    },
                    .diagonal => return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)), // triangular dense matrix + diagonal matrix
                    .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular dense matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular dense matrix + array
                .expression => Expression, // triangular dense matrix + expression
            },
            .triangular_sparse => switch (comptime types.domain(Y)) {
                .numeric => return matrix.triangular.Sparse(Coerce(types.Numeric(X), Y), types.uploOf(X), .non_unit, types.orderOf(X)), // triangular sparse matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular sparse matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix + hermitian sparse matrix
                    .triangular_dense => {
                        if (comptime types.uploOf(X) == types.uploOf(Y))
                            return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(Y)) // triangular sparse matrix + triangular dense matrix (same uplo)
                        else
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)); // triangular sparse matrix + triangular dense matrix (different uplo)
                    },
                    .triangular_sparse => {
                        if (comptime types.uploOf(X) == types.uploOf(Y))
                            return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)) // triangular sparse matrix + triangular sparse matrix (same uplo)
                        else
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)); // triangular sparse matrix + triangular sparse matrix (different uplo)
                    },
                    .diagonal => return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)), // triangular sparse matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular sparse matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular sparse matrix + array
                .expression => Expression, // triangular sparse matrix + expression
            },
            .diagonal => switch (comptime types.domain(Y)) {
                .numeric => return matrix.Diagonal(Coerce(types.Numeric(X), Y)), // diagonal matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // diagonal matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // diagonal matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // diagonal matrix + general sparse matrix
                    .symmetric_dense => return matrix.symmetric.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)), // diagonal matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.symmetric.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)), // diagonal matrix + symmetric sparse matrix
                    .hermitian_dense => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)) // diagonal matrix (complex) + hermitian dense matrix
                        else
                            return matrix.hermitian.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)); // diagonal matrix (real) + hermitian dense matrix
                    },
                    .hermitian_sparse => {
                        if (comptime types.isComplex(types.Numeric(X)))
                            return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)) // diagonal matrix (complex) + hermitian sparse matrix
                        else
                            return matrix.hermitian.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), types.orderOf(Y)); // diagonal matrix (real) + hermitian sparse matrix
                    },
                    .triangular_dense => return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), .non_unit, types.orderOf(Y)), // diagonal matrix + triangular dense matrix
                    .triangular_sparse => return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), .non_unit, types.orderOf(Y)), // diagonal matrix + triangular sparse matrix
                    .diagonal => return matrix.Diagonal(Coerce(types.Numeric(X), types.Numeric(Y))), // diagonal matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // diagonal matrix + array
                .expression => Expression, // diagonal matrix + expression
            },
            .permutation => switch (comptime types.domain(Y)) {
                .numeric => return matrix.general.Sparse(Coerce(types.Numeric(X), Y), types.orderOf(X)), // permutation matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix + vector
                .matrix => switch (comptime types.matrixType(Y)) {
                    .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // permutation matrix + general dense matrix
                    .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // permutation matrix + general sparse matrix
                    .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // permutation matrix + symmetric dense matrix
                    .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // permutation matrix + symmetric sparse matrix
                    .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // permutation matrix + hermitian dense matrix
                    .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // permutation matrix + hermitian sparse matrix
                    .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // permutation matrix + triangular dense matrix
                    .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // permutation matrix + triangular sparse matrix
                    .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix + diagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix + array
                .expression => Expression, // permutation matrix + expression
            },
            .numeric => unreachable,
        },
        .array => switch (comptime types.arrayType(X)) {
            .dense => switch (comptime types.domain(Y)) {
                .numeric => return array.Dense(Coerce(types.Numeric(X), Y), types.orderOf(X)), // dense + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense + matrix
                .array => switch (comptime types.arrayType(Y)) {
                    .dense => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // dense + dense
                    .strided => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // dense + strided
                    .sparse => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // dense + sparse
                    .numeric => unreachable,
                },
                .expression => Expression, // dense + expression
            },
            .strided => switch (comptime types.domain(Y)) {
                .numeric => return array.Dense(Coerce(types.Numeric(X), Y), types.orderOf(X)), // strided + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // strided + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // strided + matrix
                .array => switch (comptime types.arrayType(Y)) {
                    .dense => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // strided + dense
                    .strided => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // strided + strided
                    .sparse => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // strided + sparse
                    .numeric => unreachable,
                },
                .expression => Expression, // strided + expression
            },
            .sparse => switch (comptime types.domain(Y)) {
                .numeric => return array.Sparse(Coerce(types.Numeric(X), Y), types.orderOf(X)), // sparse + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse + matrix
                .array => switch (comptime types.arrayType(Y)) {
                    .dense => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // sparse + dense
                    .strided => return array.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // sparse + strided
                    .sparse => return array.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // sparse + sparse
                    .numeric => unreachable,
                },
                .expression => Expression, // sparse + expression
            },
            .numeric => unreachable,
        },
        .expression => Expression, // expression + anything
    }

    // Else two numeric types
    const xnumeric = types.numericType(X);
    const ynumeric = types.numericType(Y);

    switch (xnumeric) { // Manage comptime ints and floats everywhere we don't do yet
        .bool => switch (ynumeric) {
            .bool => return bool,
            else => return Y,
        },
        .int => switch (ynumeric) {
            .bool => return Coerce(Y, X),
            .int => {
                if (X == comptime_int)
                    return Y;

                if (Y == comptime_int)
                    return X;

                comptime var xinfo = @typeInfo(X);
                comptime var yinfo = @typeInfo(Y);

                if (xinfo.int.signedness == .unsigned) {
                    if (yinfo.int.signedness == .unsigned) {
                        // Both unsigned
                        if (xinfo.int.bits > yinfo.int.bits)
                            return X
                        else
                            return Y;
                    } else {
                        // X unsigned, Y signed
                        if (xinfo.int.bits > yinfo.int.bits) {
                            // Unsigned is larger than signed
                            if (std.mem.indexOfScalar(type, &types.standard_integer_types, X) != null and
                                std.mem.indexOfScalar(type, &types.standard_integer_types, Y) != null)
                            {
                                // Both are standard integers, need to double
                                // bits, unless already at max, then only add 1
                                if (xinfo.int.bits == 128) {
                                    yinfo.int.bits = xinfo.int.bits + 1; // 129
                                    return @Type(yinfo);
                                } else {
                                    yinfo.int.bits = xinfo.int.bits * 2; // Another standard size
                                    return @Type(yinfo);
                                }
                            } else {
                                // One of the types is not a standard integer,
                                // only need to increase max bits by 1
                                yinfo.int.bits = xinfo.int.bits + 1;
                                return @Type(yinfo);
                            }
                        } else if (xinfo.int.bits == yinfo.int.bits) {
                            // Equal size, need to add 1 bit to signed
                            yinfo.int.bits += 1;
                            return @Type(yinfo);
                        } else {
                            // Signed is larger than unsigned
                            return Y;
                        }
                    }
                } else {
                    if (yinfo.int.signedness == .unsigned) {
                        // X signed, Y unsigned
                        if (yinfo.int.bits > xinfo.int.bits) {
                            // Unsigned is larger than signed
                            if (std.mem.indexOfScalar(type, &types.standard_integer_types, X) != null and
                                std.mem.indexOfScalar(type, &types.standard_integer_types, Y) != null)
                            {
                                // Both are standard integers, need to double
                                // bits, unless already at max, then only add 1
                                if (yinfo.int.bits == 128) {
                                    xinfo.int.bits = yinfo.int.bits + 1; // 129
                                    return @Type(xinfo);
                                } else {
                                    xinfo.int.bits = yinfo.int.bits * 2; // Another standard size
                                    return @Type(xinfo);
                                }
                            } else {
                                // One of the types is not a standard integer,
                                // only need to increase max bits by 1
                                xinfo.int.bits = yinfo.int.bits + 1;
                                return @Type(xinfo);
                            }
                        } else if (xinfo.int.bits == yinfo.int.bits) {
                            // Equal size, need to add 1 bit to signed
                            xinfo.int.bits += 1;
                            return @Type(xinfo);
                        } else {
                            // Signed is larger than unsigned
                            return X;
                        }
                    } else {
                        // Both signed
                        if (xinfo.int.bits > yinfo.int.bits)
                            return X
                        else
                            return Y;
                    }
                }
            },
            .float => {
                if (Y == comptime_float)
                    return EnsureFloat(X);

                return Coerce(EnsureFloat(X), Y);
            },
            .dyadic => {
                const xinfo = @typeInfo(X);

                return Dyadic(
                    int.max(xinfo.int.bits, Y.mantissaBits()),
                    Y.exponentBits(),
                );
            },
            .cfloat => {
                return Cfloat(Coerce(X, types.Scalar(Y)));
            },
            else => return Y,
        },
        .float => {
            switch (ynumeric) {
                .bool => return Coerce(Y, X),
                .int => return Coerce(Y, X),
                .float => {
                    if (X == comptime_float)
                        return Y;

                    if (Y == comptime_float)
                        return X;

                    const xinfo = @typeInfo(X);
                    const yinfo = @typeInfo(Y);

                    if (xinfo.float.bits > yinfo.float.bits)
                        return X
                    else
                        return Y;
                },
                .dyadic => {
                    const x_mantissa_bits = std.math.floatMantissaBits(X) + 1;
                    const x_exponent_bits = std.math.floatExponentBits(X);

                    return Dyadic(
                        int.max(x_mantissa_bits, Y.mantissaBits()),
                        int.max(x_exponent_bits, Y.exponentBits()),
                    );
                },
                .cfloat => {
                    return Cfloat(Coerce(X, types.Scalar(Y)));
                },
                else => return Y,
            }
        },
        .dyadic => switch (ynumeric) {
            .bool => return Coerce(Y, X),
            .int => return Coerce(Y, X),
            .float => return Coerce(Y, X),
            .dyadic => {
                return Dyadic(
                    int.max(X.mantissaBits(), Y.mantissaBits()),
                    int.max(X.exponentBits(), Y.exponentBits()),
                );
            },
            .cfloat => {
                return Cfloat(Coerce(X, types.Scalar(Y)));
            },
            else => return Y,
        },
        .cfloat => {
            switch (ynumeric) {
                .bool => return Coerce(Y, X),
                .int => return Coerce(Y, X),
                .float => return Coerce(Y, X),
                .dyadic => return Coerce(Y, X),
                .cfloat => {
                    return cfloat.Cfloat(Coerce(types.Scalar(X), types.Scalar(Y)));
                },
                .integer => return Complex(Rational),
                .rational => return Complex(Rational),
                .real => return Complex(Real),
                .complex => return Y,
            }
        },
        .integer => switch (ynumeric) {
            .bool => return Coerce(Y, X),
            .int => return Coerce(Y, X),
            .float => return Coerce(Y, X),
            .cfloat => return Coerce(Y, X),
            .integer => return Integer,
            else => return Y,
        },
        .rational => switch (ynumeric) {
            .bool => return Coerce(Y, X),
            .int => return Coerce(Y, X),
            .float => return Coerce(Y, X),
            .cfloat => return Coerce(Y, X),
            .integer => return Coerce(Y, X),
            .rational => return X,
            else => return Y,
        },
        .real => switch (ynumeric) {
            .bool => return Coerce(Y, X),
            .int => return Coerce(Y, X),
            .float => return Coerce(Y, X),
            .cfloat => return Coerce(Y, X),
            .integer => return Coerce(Y, X),
            .rational => return Coerce(Y, X),
            .real => return Real,
            .complex => return Complex(Real),
        },
        .complex => switch (ynumeric) {
            .bool => return Coerce(Y, X),
            .int => return Coerce(Y, X),
            .float => return Coerce(Y, X),
            .cfloat => return Coerce(Y, X),
            .integer => return Coerce(Y, X),
            .rational => return Coerce(Y, X),
            .real => return Coerce(Y, X),
            .complex => return Complex(Coerce(types.Scalar(X), types.Scalar(Y))),
        },
    }
}

/// Coerces the input types to the smallest type that can represent the result
/// of their multiplication.
///
/// For scalar or array types, this is equivalent to `Coerce`. For matrices and
/// vectors, this function takes into account the rules of linear algebra to
/// determine the appropriate resulting type.
///
/// Parameters
/// ----------
/// comptime X (`type`): The first type to coerce. Must be a supported numeric
/// type, a vector, a matrix, or an array.
///
/// comptime Y (`type`): The second type to coerce. Must be a supported numeric
/// type, a vector, a matrix, or an array.
///
/// Returns
/// -------
/// `type`: The coerced type that can represent the result of multiplying `X`
/// and `Y`.
pub fn MulCoerce(comptime X: type, comptime Y: type) type {
    switch (comptime types.domain(X)) {
        .numeric => return Coerce(X, Y), // Same as Coerce
        .vector => switch (comptime types.vectorType(X)) {
            .dense => switch (comptime types.domain(Y)) {
                .numeric => return Coerce(X, Y), // Same as Coerce
                .vector => switch (comptime types.vectorType(Y)) {
                    .dense => return Coerce(types.Numeric(X), types.Numeric(Y)), // dense vector * dense vector
                    .sparse => return Coerce(types.Numeric(X), types.Numeric(Y)), // dense vector * sparse vector
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
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector * array
                .expression => Expression, // sparse vector * expression
            },
            .numeric => unreachable,
        },
        .matrix => {
            switch (comptime types.matrixType(X)) {
                .general_dense => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.domain(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // general dense matrix * dense vector
                        .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // general dense matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * general dense matrix
                        .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * diagonal matrix
                        .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general dense matrix * permutation matrix
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general dense matrix * array
                    .expression => Expression, // general dense matrix * expression
                },
                .general_sparse => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.domain(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // general sparse matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // general sparse matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // general sparse matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // general sparse matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // general sparse matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // general sparse matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // general sparse matrix * permutation matrix
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general sparse matrix * array
                    .expression => Expression, // general sparse matrix * expression
                },
                .symmetric_dense => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.domain(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // symmetric dense matrix * dense vector
                        .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // symmetric dense matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric dense matrix * general dense matrix
                        .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * diagonal matrix
                        .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric dense matrix * permutation matrix
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric dense matrix * array
                    .expression => Expression, // symmetric dense matrix * expression
                },
                .symmetric_sparse => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.domain(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // symmetric sparse matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // symmetric sparse matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric sparse matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric sparse matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric sparse matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric sparse matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric sparse matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric sparse matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // symmetric sparse matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric sparse matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric sparse matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // symmetric sparse matrix * permutation matrix
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric sparse matrix * array
                    .expression => Expression, // symmetric sparse matrix * expression
                },
                .hermitian_dense => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.domain(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // hermitian dense matrix * dense vector
                        .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // hermitian dense matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian dense matrix * general dense matrix
                        .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * diagonal matrix
                        .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian dense matrix * permutation matrix
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian dense matrix * array
                    .expression => Expression, // hermitian dense matrix * expression
                },
                .hermitian_sparse => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.domain(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // hermitian sparse matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // hermitian sparse matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian sparse matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian sparse matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian sparse matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian sparse matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian sparse matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian sparse matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // hermitian sparse matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian sparse matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian sparse matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // hermitian sparse matrix * permutation matrix
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // hermitian sparse matrix * array
                    .expression => Expression, // hermitian sparse matrix * expression
                },
                .triangular_dense => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.domain(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // triangular dense matrix * dense vector
                        .sparse => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // triangular dense matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular dense matrix * general dense matrix
                        .general_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular dense matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular dense matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular dense matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular dense matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular dense matrix * hermitian sparse matrix
                        .triangular_dense => {
                            if (comptime types.uploOf(X) == types.uploOf(Y)) {
                                if (comptime types.diagOf(X) == types.diagOf(Y))
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.diagOf(X), types.orderOf(X)) // triangular dense matrix * triangular dense matrix (same uplo and diag)
                                else
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)); // triangular dense matrix * triangular dense matrix (same uplo, different diag)
                            } else {
                                return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)); // triangular dense matrix * triangular dense matrix (different uplo)
                            }
                        },
                        .triangular_sparse => {
                            if (comptime types.uploOf(X) == types.uploOf(Y)) {
                                if (comptime types.diagOf(X) == types.diagOf(Y))
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.diagOf(X), types.orderOf(X)) // triangular dense matrix * triangular sparse matrix (same uplo and diag)
                                else
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)); // triangular dense matrix * triangular sparse matrix (same uplo, different diag)
                            } else {
                                return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)); // triangular dense matrix * triangular sparse matrix (different uplo)
                            }
                        },
                        .diagonal => return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)), // triangular dense matrix * diagonal matrix
                        .permutation => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular dense matrix * permutation matrix
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // triangular dense matrix * array
                    .expression => Expression, // triangular dense matrix * expression
                },
                .triangular_sparse => switch (comptime types.domain(Y)) {
                    .numeric => return Coerce(X, Y), // Same as Coerce
                    .vector => switch (comptime types.domain(Y)) {
                        .dense => return vector.Dense(Coerce(types.Numeric(X), types.Numeric(Y))), // triangular sparse matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(types.Numeric(X), types.Numeric(Y))), // triangular sparse matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)), // triangular sparse matrix * hermitian sparse matrix
                        .triangular_dense => {
                            if (comptime types.uploOf(X) == types.uploOf(Y)) {
                                if (comptime types.diagOf(X) == types.diagOf(Y))
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.diagOf(X), types.orderOf(Y)) // triangular sparse matrix * triangular dense matrix (same uplo and diag)
                                else
                                    return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(Y)); // triangular sparse matrix * triangular dense matrix (same uplo, different diag)
                            } else {
                                return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(Y)); // triangular sparse matrix * triangular dense matrix (different uplo)
                            }
                        },
                        .triangular_sparse => {
                            if (comptime types.uploOf(X) == types.uploOf(Y)) {
                                if (comptime types.diagOf(X) == types.diagOf(Y))
                                    return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), types.diagOf(X), types.orderOf(X)) // triangular sparse matrix * triangular sparse matrix (same uplo and diag)
                                else
                                    return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)); // triangular sparse matrix * triangular sparse matrix (same uplo, different diag)
                            } else {
                                return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)); // triangular sparse matrix * triangular sparse matrix (different uplo)
                            }
                        },
                        .diagonal => return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(X), .non_unit, types.orderOf(X)), // triangular sparse matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // triangular sparse matrix * permutation matrix
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
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.triangular.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), .non_unit, types.orderOf(X)), // diagonal matrix * triangular dense matrix
                        .triangular_sparse => return matrix.triangular.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.uploOf(Y), .non_unit, types.orderOf(X)), // diagonal matrix * triangular sparse matrix
                        .diagonal => return matrix.Diagonal(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // diagonal matrix * permutation matrix
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
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime types.matrixType(Y)) {
                        .general_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * general dense matrix
                        .general_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * general sparse matrix
                        .symmetric_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * symmetric dense matrix
                        .symmetric_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * symmetric sparse matrix
                        .hermitian_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * hermitian dense matrix
                        .hermitian_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * hermitian sparse matrix
                        .triangular_dense => return matrix.general.Dense(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * triangular dense matrix
                        .triangular_sparse => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * triangular sparse matrix
                        .diagonal => return matrix.general.Sparse(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * diagonal matrix
                        .permutation => return matrix.Permutation(Coerce(types.Numeric(X), types.Numeric(Y)), types.orderOf(X)), // permutation matrix * permutation matrix
                        .numeric => unreachable,
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix * array
                    .expression => Expression, // permutation matrix * expression
                },
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
/// Parameters
/// ----------
/// comptime F (`type`): The type to check if it can be coerced. Must be a
/// supported numeric type.
///
/// comptime T (`type`): The target type. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if `F` can be coerced to `T` without loss of information,
/// `false` otherwise.
pub fn canCoerce(comptime F: type, comptime T: type) bool {
    comptime if (!types.isNumeric(F) or !types.isNumeric(T))
        @compileError("canCoerce can only be used with numeric types, got F: " ++ @typeName(F) ++ " and T: " ++ @typeName(T));

    return Coerce(F, T) == T;
}

/// Coerces the input type to the smallest non-integral numeric type that can
/// represent all values of the input type.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to coerce. Must be a supported numeric type, a
/// vector, a matrix, or an array.
///
/// Returns
/// -------
/// `type`: The coerced type.
pub fn EnsureFloat(comptime T: type) type {
    if (types.isExpression(T)) {
        return Expression;
    } else if (types.isArray(T)) {
        switch (types.arrayType(T)) {
            .dense => return array.Dense(EnsureFloat(types.Numeric(T)), types.orderOf(T)),
            .strided => return array.Strided(EnsureFloat(types.Numeric(T)), types.orderOf(T)),
            .sparse => return array.Sparse(EnsureFloat(types.Numeric(T)), types.orderOf(T)),
            .numeric => unreachable,
        }
    } else if (types.isMatrix(T)) {
        switch (types.matrixType(T)) {
            .general_dense => return matrix.general.Dense(EnsureFloat(types.Numeric(T)), types.orderOf(T)),
            .general_sparse => return matrix.general.Sparse(EnsureFloat(types.Numeric(T)), types.orderOf(T)),
            .symmetric_dense => return matrix.symmetric.Dense(EnsureFloat(types.Numeric(T)), types.uploOf(T), types.orderOf(T)),
            .symmetric_sparse => return matrix.symmetric.Sparse(EnsureFloat(types.Numeric(T)), types.uploOf(T), types.orderOf(T)),
            .hermitian_dense => return matrix.hermitian.Dense(EnsureFloat(types.Numeric(T)), types.uploOf(T), types.orderOf(T)),
            .hermitian_sparse => return matrix.hermitian.Sparse(EnsureFloat(types.Numeric(T)), types.uploOf(T), types.orderOf(T)),
            .triangular_dense => return matrix.triangular.Dense(EnsureFloat(types.Numeric(T)), types.uploOf(T), types.diagOf(T), types.orderOf(T)),
            .triangular_sparse => return matrix.triangular.Sparse(EnsureFloat(types.Numeric(T)), types.uploOf(T), types.diagOf(T), types.orderOf(T)),
            .diagonal => return matrix.Diagonal(EnsureFloat(types.Numeric(T))),
            .permutation => return matrix.Permutation(EnsureFloat(types.Numeric(T))),
            .numeric => unreachable,
        }
    } else if (types.isVector(T)) {
        switch (types.vectorType(T)) {
            .dense => return vector.Dense(EnsureFloat(types.Numeric(T))),
            .sparse => return vector.Sparse(EnsureFloat(types.Numeric(T))),
            .numeric => unreachable,
        }
    }

    switch (types.numericType(T)) {
        .bool => return f32, // Smallest float type
        .int => {
            const tinfo = @typeInfo(T);
            if (tinfo.int.bits <= 11)
                return f16;

            if (tinfo.int.bits <= 24)
                return f32;

            if (tinfo.int.bits <= 64) // Lossy past 53, but to not explode width
                return f64;

            return f128;
        },
        .float => return T,
        .dyadic => return T,
        .cfloat => return T,
        .integer => return Rational,
        .rational => return T,
        .real => return T,
        .complex => return T,
    }
}

/// Coerces the second type to a vector, matrix or array type based on the first
/// type.
///
/// This function is useful for ensuring that the second type is always a
/// vector, matrix or an array when the first is.
///
/// Parameters
/// ----------
/// comptime X (`type`): The type to check. Must be a supported numeric type, a
/// vector, a matrix, or an array.
///
/// comptime Y (`type`): The type to coerce. Must be a supported numeric type.
///
/// Returns
/// -------
/// `type`: The coerced type.
pub fn EnsureDomain(comptime X: type, comptime Y: type) type {
    if (!types.isNumeric(Y))
        @compileError("EnsureDomain requires Y to be a numeric type, got: " ++ @typeName(Y));

    switch (comptime types.domain(X)) {
        .numeric => return Y,
        .vector => return EnsureVector(X, Y),
        .matrix => return EnsureMatrix(X, Y),
        .array => return EnsureArray(X, Y),
        .expression => return Expression,
    }
}

pub fn EnsureVector(comptime X: type, comptime Y: type) type {
    if (!types.isNumeric(Y))
        @compileError("EnsureVector requires Y to be a numeric type, got: " ++ @typeName(Y));

    switch (comptime types.domain(X)) {
        .numeric => return Y,
        .vector => switch (types.vectorType(X)) {
            .dense => return vector.Dense(Y),
            .sparse => return vector.Sparse(Y),
            .numeric => unreachable,
        },
        .matrix => return Y,
        .array => return Y,
        .expression => return Y,
    }
}

pub fn EnsureMatrix(comptime X: type, comptime Y: type) type {
    if (!types.isNumeric(Y))
        @compileError("EnsureMatrix requires Y to be a numeric type, got: " ++ @typeName(Y));

    switch (comptime types.domain(X)) {
        .numeric => return Y,
        .vector => return Y,
        .matrix => switch (types.matrixType(X)) {
            .general_dense => return matrix.general.Dense(Y, types.orderOf(X)),
            .general_sparse => return matrix.general.Sparse(Y, types.orderOf(X)),
            .symmetric_dense => return matrix.symmetric.Dense(Y, types.uploOf(X), types.orderOf(X)),
            .symmetric_sparse => return matrix.symmetric.Sparse(Y, types.uploOf(X), types.orderOf(X)),
            .hermitian_dense => return matrix.hermitian.Dense(Y, types.uploOf(X), types.orderOf(X)),
            .hermitian_sparse => return matrix.hermitian.Sparse(Y, types.uploOf(X), types.orderOf(X)),
            .triangular_dense => return matrix.triangular.Dense(Y, types.uploOf(X), types.diagOf(X), types.orderOf(X)),
            .triangular_sparse => return matrix.triangular.Sparse(Y, types.uploOf(X), types.diagOf(X), types.orderOf(X)),
            .diagonal => return matrix.Diagonal(Y),
            .permutation => return matrix.Permutation(Y),
            .numeric => unreachable,
        },
        .array => return Y,
        .expression => return Y,
    }
}

/// Coerces the second type to an array type based on the first type.
///
/// This function is useful for ensuring that the second type is always an array
/// when the first is. If the first type is a strided array, the second type
/// will be coerced to a dense array type.
///
/// Parameters
/// ----------
/// comptime X (`type`): The type to check. Must be a supported numeric type or
/// an array.
///
/// comptime Y (`type`): The type to coerce. Must be a supported numeric type.
///
/// Returns
/// -------
/// `type`: The coerced type.
pub fn EnsureArray(comptime X: type, comptime Y: type) type {
    if (!types.isNumeric(Y))
        @compileError("EnsureArray requires Y to be a numeric type, got: " ++ @typeName(Y));

    switch (comptime types.domain(X)) {
        .numeric => return Y,
        .vector => return Y,
        .matrix => return Y,
        .array => switch (types.arrayType(X)) {
            .dense => return array.Dense(Y, types.orderOf(X)),
            .strided => return array.Dense(Y, types.orderOf(X)),
            .sparse => return array.Sparse(Y, types.orderOf(X)),
            .numeric => unreachable,
        },
        .expression => return Y,
    }
}
