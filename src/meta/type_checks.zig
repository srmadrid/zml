//! Namespace for type checking functions.

const std = @import("std");

const meta = @import("../meta.zig");

/// Checks if the input type is a supported type (numeric, vector, matrix,
/// array, polynomial).
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a supported type, `false` otherwise.
pub fn isSupportedType(comptime T: type) bool {
    return isNumeric(T) or
        isVector(T) or
        isMatrix(T) or
        isArray(T) or
        isPoly(T);
}

/// Checks if the input type is a one-item pointer.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a one-item pointer, `false` otherwise.
pub fn isPointer(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .pointer => |info| blk: {
            if (info.size != .one and
                info.size != .c) break :blk false;

            break :blk true;
        },
        else => false,
    };
}

/// Checks if the input type is a many-item pointer.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a many-item pointer, `false` otherwise.
pub fn isManyItemPointer(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .pointer => |info| blk: {
            if (info.size != .many and
                info.size != .c) break :blk false;

            break :blk true;
        },
        else => false,
    };
}

/// Checks if the input type is a constant pointer. Works for one-item pointers,
/// many-item pointers, and slices.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a constant one-item pointer, `false`
/// otherwise.
pub fn isConstPointer(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .pointer => |info| if (info.is_const)
            true
        else
            false,
        else => false,
    };
}

/// Checks if the input type is a slice.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a slice, `false` otherwise.
pub fn isSlice(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .pointer => |info| blk: {
            if (info.size != .slice)
                break :blk false;

            break :blk true;
        },
        else => false,
    };
}

/// Checks if the input type is a simd vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a simd vector, `false` otherwise.
pub fn isSimdVector(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .vector => true,
        else => false,
    };
}

/// Checks if the input type `T` is a numeric type.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a numeric type, `false` otherwise.
pub fn isNumeric(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .bool => true,
        .int, .comptime_int => true,
        .float, .comptime_float => true,
        .@"struct", .@"union" => meta.hasBoolDecl(T, "is_numeric"),
        else => false,
    };
}

/// Checks if the input type is a custom numeric type.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a custom numeric type, `false` otherwise.
pub fn isCustomNumeric(comptime T: type) bool {
    return isNumeric(T) and meta.numericType(T) == .custom;
}

/// Checks if the input type is an instance of a vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a vector, `false` otherwise.
pub fn isVector(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_vector"),
        else => false,
    };
}

/// Checks if the input type is an instance of a static vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a static vector, `false` otherwise.
pub fn isStaticVector(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_vector") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a dense vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a dense vector, `false` otherwise.
pub fn isDenseVector(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_vector") and meta.hasBoolDecl(T, "is_dense"),
        else => false,
    };
}

/// Checks if the input type is an instance of a sparse vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a sparse vector, `false` otherwise.
pub fn isSparseVector(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_vector") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a matrix, `false` otherwise.
pub fn isMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix"),
        else => false,
    };
}

/// Checks if the input type is a square matrix, i.e., a matrix type that is
/// always square (symmetric or hermitian).
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a square matrix, `false` otherwise.
pub fn isSquareMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_symmetric") or meta.hasBoolDecl(T, "is_hermitian") or meta.hasBoolDecl(T, "is_permutation"),
        else => false,
    };
}

/// Checks if the input type is an instance of a general matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a general matrix, `false` otherwise.
pub fn isGeneralMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_general"),
        else => false,
    };
}

/// Checks if the input type is an instance of a general static matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a general static matrix, `false` otherwise.
pub fn isGeneralStaticMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_general") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a general dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a general dense matrix, `false` otherwise.
pub fn isGeneralDenseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_general") and meta.hasBoolDecl(T, "is_dense"),
        else => false,
    };
}

/// Checks if the input type is an instance of a general sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a general sparse matrix, `false` otherwise.
pub fn isGeneralSparseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_general") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a symmetric matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a symmetric matrix, `false` otherwise.
pub fn isSymmetricMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_symmetric"),
        else => false,
    };
}

/// Checks if the input type is an instance of a symmetric static matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a symmetric static matrix, `false` otherwise.
pub fn isSymmetricStaticMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_symmetric") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a symmetric dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a symmetric dense matrix, `false` otherwise.
pub fn isSymmetricDenseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_symmetric") and meta.hasBoolDecl(T, "is_dense"),
        else => false,
    };
}

/// Checks if the input type is an instance of a symmetric sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a symmetric sparse matrix, `false` otherwise.
pub fn isSymmetricSparseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_symmetric") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a hermitian matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a hermitian matrix, `false` otherwise.
pub fn isHermitianMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_hermitian"),
        else => false,
    };
}

/// Checks if the input type is an instance of a hermitian static matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a hermitian static matrix, `false` otherwise.
pub fn isHermitianStaticMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_hermitian") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a hermitian dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a hermitian dense matrix, `false` otherwise.
pub fn isHermitianDenseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_hermitian") and meta.hasBoolDecl(T, "is_dense"),
        else => false,
    };
}

/// Checks if the input type is an instance of a hermitian sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a hermitian sparse matrix, `false` otherwise.
pub fn isHermitianSparseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_hermitian") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a triangular matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a triangular matrix, `false` otherwise.
pub fn isTriangularMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_triangular"),
        else => false,
    };
}

/// Checks if the input type is an instance of a triangular static matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a triangular static matrix, `false` otherwise.
pub fn isTriangularStaticMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_triangular") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a triangular dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a triangular dense matrix, `false` otherwise.
pub fn isTriangularDenseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_triangular") and meta.hasBoolDecl(T, "is_dense"),
        else => false,
    };
}

/// Checks if the input type is an instance of a triangular sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a triangular sparse matrix, `false` otherwise.
pub fn isTriangularSparseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_triangular") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a diagonal matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a diagonal matrix, `false` otherwise.
pub fn isDiagonalMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_diagonal"),
        else => false,
    };
}

/// Checks if the input type is an instance of a diagonal static matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a diagonal static matrix, `false` otherwise.
pub fn isDiagonalStaticMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_diagonal") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a diagonal sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a diagonal sparse matrix, `false` otherwise.
pub fn isDiagonalSparseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_diagonal") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a permutation matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a permutation matrix, `false` otherwise.
pub fn isPermutationMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_permutation"),
        else => false,
    };
}

/// Checks if the input type is an instance of a permutation static matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a permutation static matrix, `false` otherwise.
pub fn isPermutationStaticMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_permutation") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a permutation sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a permutation sparse matrix, `false` otherwise.
pub fn isPermutationSparseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_permutation") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a builder matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a builder matrix, `false` otherwise.
pub fn isBuilderMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_builder"),
        else => false,
    };
}

/// Checks if the input type is an instance of a builder sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a builder sparse matrix, `false` otherwise.
pub fn isBuilderSparseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_builder") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is a static matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a static matrix, `false` otherwise.
pub fn isStaticMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is a dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a dense matrix, `false` otherwise.
pub fn isDenseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_dense"),
        else => false,
    };
}

/// Checks if the input type is a sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a sparse matrix, `false` otherwise.
pub fn isSparseMatrix(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_matrix") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of an array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is an array, `false` otherwise.
pub fn isArray(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_array"),
        else => false,
    };
}

/// Checks if the input type is an instance of a static array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a static array, `false` otherwise.
pub fn isStaticArray(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_array") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a dense array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a dense array, `false` otherwise.
pub fn isDenseArray(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_array") and meta.hasBoolDecl(T, "is_dense"),
        else => false,
    };
}

/// Checks if the input type is an instance of a sparse array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a sparse array, `false` otherwise.
pub fn isSparseArray(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_array") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a builder sparse array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a builder sparse array, `false` otherwise.
pub fn isBuilderSparseArray(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_array") and meta.hasBoolDecl(T, "is_builder") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input type is an instance of a polynomial.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a polynomial, `false` otherwise.
pub fn isPoly(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_poly"),
        else => false,
    };
}

/// Checks if the input type is an instance of a static polynomial.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a static polynomial, `false` otherwise.
pub fn isStaticPoly(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_poly") and meta.hasBoolDecl(T, "is_static"),
        else => false,
    };
}

/// Checks if the input type is an instance of a dense polynomial.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a dense polynomial, `false` otherwise.
pub fn isDensePoly(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_poly") and meta.hasBoolDecl(T, "is_dense"),
        else => false,
    };
}

/// Checks if the input type is an instance of a sparse polynomial.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a sparse polynomial, `false` otherwise.
pub fn isSparsePoly(comptime T: type) bool {
    return switch (comptime @typeInfo(T)) {
        .@"struct" => meta.hasBoolDecl(T, "is_poly") and meta.hasBoolDecl(T, "is_sparse"),
        else => false,
    };
}

/// Checks if the input numeric type is integral.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check. Must be a supported numeric
///   type.
///
/// ## Returns
/// `bool`: `true` if the type is integral, `false` otherwise.
pub fn isIntegral(comptime N: type) bool {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.meta.isIntegral: " ++ @typeName(N) ++ " is not a supported numeric type");

    return switch (comptime meta.numericType(N)) {
        .bool => true,
        .int => true,
        .float => false,
        else => meta.hasBoolDecl(N, "is_integral"),
    };
}

/// Checks if the input numeric type is non-integral.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check. Must be a supported numeric
///   type.
///
/// ## Returns
/// `bool`: `true` if the type is non-integral, `false` otherwise.
pub fn isNonIntegral(comptime N: type) bool {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.meta.isNonIntegral: " ++ @typeName(N) ++ " is not a supported numeric type");

    return !isIntegral(N);
}

/// Checks if the input numeric type is real.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check. Must be a supported numeric
///   type.
///
/// ## Returns
/// `bool`: `true` if the type is real, `false` otherwise.
pub fn isReal(comptime N: type) bool {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.meta.isReal: " ++ @typeName(N) ++ " is not a supported numeric type");

    return !isComplex(N);
}

/// Checks if the input numeric type is complex.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check. Must be a supported numeric
///   type.
///
/// ## Returns
/// `bool`: `true` if the type is complex, `false` otherwise.
pub fn isComplex(comptime N: type) bool {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.meta.isComplex: " ++ @typeName(N) ++ " is not a supported numeric type");

    return switch (comptime meta.numericType(N)) {
        .bool => false,
        .int => false,
        .float => false,
        else => meta.hasBoolDecl(N, "is_complex"),
    };
}

/// Checks if the input numeric type is signed.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check. Must be a supported numeric
///   type.
///
/// ## Returns
/// `bool`: `true` if the type is signed, `false` otherwise.
pub fn isSigned(comptime N: type) bool {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.meta.isSigned: " ++ @typeName(N) ++ " is not a supported numeric type");

    return !isUnsigned(N);
}

/// Checks if the input numeric type is unsigned.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check. Must be a supported numeric
///   type.
///
/// ## Returns
/// `bool`: `true` if the type is unsigned, `false` otherwise.
pub fn isUnsigned(comptime N: type) bool {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.meta.isUnsigned: " ++ @typeName(N) ++ " is not a supported numeric type");

    return switch (comptime meta.numericType(N)) {
        .bool => true,
        .int => blk: {
            switch (comptime @typeInfo(N)) {
                .int => |info| break :blk info.signedness == .unsigned,
                .comptime_int => break :blk false,
                else => unreachable,
            }
        },
        .float => false,
        else => meta.hasBoolDecl(N, "is_unsigned"),
    };
}
