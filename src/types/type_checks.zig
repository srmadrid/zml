//! Namespace for type checking functions.

const std = @import("std");

const types = @import("../types.zig");

/// Checks if the input type is a supported type (numeric, vector, matrix,
/// array, or expression).
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
        isExpression(T);
}

/// Checks if the input type is a one-item pointer.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a one-item pointer, `false` otherwise.
pub fn isPointer(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .one and
                info.size != .c) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input type is a many-item pointer.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a many-item pointer, `false` otherwise.
pub fn isManyPointer(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .many and
                info.size != .c) return false;

            return true;
        },
        else => return false,
    }
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
    switch (comptime @typeInfo(T)) {
        .pointer => |info| {
            if (info.is_const) {
                return true;
            } else {
                return false;
            }
        },
        else => return false,
    }
}

/// Checks if the input type is a slice.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a slice, `false` otherwise.
pub fn isSlice(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .slice) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input type is a simd vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a simd vector, `false` otherwise.
pub fn isSimdVector(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .vector => return true,
        else => return false,
    }
}

/// Checks if the input type `T` is a numeric type.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a numeric type, `false` otherwise.
pub fn isNumeric(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .bool => return true,
        .int, .comptime_int => return true,
        .float, .comptime_float => return true,
        .@"struct" => {
            if (comptime @hasDecl(T, "is_dyadic"))
                return true;

            if (comptime (@hasDecl(T, "is_cfloat")) or
                T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or
                T == std.math.Complex(f80) or T == std.math.Complex(f128) or T == std.math.Complex(comptime_float))
                return true;

            if (comptime @hasDecl(T, "is_integer"))
                return true;

            if (comptime @hasDecl(T, "is_rational"))
                return true;

            if (comptime @hasDecl(T, "is_real"))
                return true;

            if (comptime @hasDecl(T, "is_complex"))
                return true;

            if (comptime @hasDecl(T, "is_numeric"))
                return true;

            return false;
        },
        else => return false,
    }
}

/// Checks if the input type is an instance of a vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a vector, `false` otherwise.
pub fn isVector(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a dense vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a dense vector, `false` otherwise.
pub fn isDenseVector(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a sparse vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a sparse vector, `false` otherwise.
pub fn isSparseVector(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a custom vector.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a custom vector, `false` otherwise.
pub fn isCustomVector(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector") and !@hasDecl(T, "is_dense") and !@hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a matrix, `false` otherwise.
pub fn isMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix"),
        else => return false,
    }
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
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_symmetric") or @hasDecl(T, "is_hermitian"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a general matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a general matrix, `false` otherwise.
pub fn isGeneralMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a general dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a general dense matrix, `false` otherwise.
pub fn isGeneralDenseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a general sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a general sparse matrix, `false` otherwise.
pub fn isGeneralSparseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a symmetric matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a symmetric matrix, `false` otherwise.
pub fn isSymmetricMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a symmetric dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a symmetric dense matrix, `false` otherwise.
pub fn isSymmetricDenseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a symmetric sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a symmetric sparse matrix, `false` otherwise.
pub fn isSymmetricSparseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a hermitian matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a hermitian matrix, `false` otherwise.
pub fn isHermitianMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a hermitian dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a hermitian dense matrix, `false` otherwise.
pub fn isHermitianDenseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a hermitian sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a hermitian sparse matrix, `false` otherwise.
pub fn isHermitianSparseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a triangular matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a triangular matrix, `false` otherwise.
pub fn isTriangularMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_triangular"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a triangular dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a triangular dense matrix, `false` otherwise.
pub fn isTriangularDenseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_triangular") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a triangular sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a triangular sparse matrix, `false` otherwise.
pub fn isTriangularSparseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_triangular") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a diagonal matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a diagonal matrix, `false` otherwise.
pub fn isDiagonalMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_diagonal"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a permutation matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a permutation matrix, `false` otherwise.
pub fn isPermutationMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_permutation"),
        else => return false,
    }
}

/// Checks if the input type is a dense matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a dense matrix, `false` otherwise.
pub fn isDenseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is a sparse matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a sparse matrix, `false` otherwise.
pub fn isSparseMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a custom matrix.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a custom matrix, `false` otherwise.
pub fn isCustomMatrix(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and !@hasDecl(T, "is_dense") and !@hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of an array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is an array, `false` otherwise.
pub fn isArray(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a dense array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a dense array, `false` otherwise.
pub fn isDenseArray(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a strided array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a strided array, `false` otherwise.
pub fn isStridedArray(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and @hasDecl(T, "is_strided"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a sparse array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a sparse array, `false` otherwise.
pub fn isSparseArray(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a custom array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `bool`: `true` if the type is a custom array, `false` otherwise.
pub fn isCustomArray(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and !@hasDecl(T, "is_dense") and !@hasDecl(T, "is_strided") and !@hasDecl(T, "is_sparse"),
        else => return false,
    }
}

pub fn isExpression(comptime T: type) bool {
    _ = T;
    return false;
    // return T == Expression;
}

/// Checks if the input numeric type requires allocation.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check. Must be a supported numeric
///   type.
///
/// ## Returns
/// `bool`: `true` if the type requires allocation, `false` otherwise.
pub fn isAllocated(comptime N: type) bool {
    if (comptime !types.isNumeric(N))
        @compileError("zml.types.isAllocated: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (comptime @typeInfo(N)) {
        .@"struct" => return @hasDecl(N, "is_allocated") and N.is_allocated,
        else => return false,
    }
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
    if (comptime !types.isNumeric(N))
        @compileError("zml.types.isIntegral: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (comptime types.numericType(N)) {
        .bool => return true,
        .int => return true,
        else => return @hasDecl(N, "is_integral"),
    }
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
    if (comptime !types.isNumeric(N))
        @compileError("zml.types.isNonIntegral: " ++ @typeName(N) ++ " is not a supported numeric type");

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
pub fn isRealType(comptime N: type) bool {
    if (comptime !types.isNumeric(N))
        @compileError("zml.types.isRealType: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (comptime types.numericType(N)) {
        .bool => return true,
        .int => return true,
        .float => return true,
        else => return @hasDecl(N, "is_real_type"),
    }
}

/// Checks if the input numeric type is complex.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check. Must be a supported numeric
///   type.
///
/// ## Returns
/// `bool`: `true` if the type is complex, `false` otherwise.
pub fn isComplexType(comptime N: type) bool {
    if (comptime !types.isNumeric(N))
        @compileError("zml.types.isComplexType: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (comptime types.numericType(N)) {
        .bool => return false,
        .int => return false,
        .float => return false,
        else => return @hasDecl(N, "is_complex_type"),
    }
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
    if (comptime !types.isNumeric(N))
        @compileError("zml.types.isSigned: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (comptime types.numericType(N)) {
        .int => {
            switch (comptime @typeInfo(N)) {
                .int => |info| return info.signedness == .signed,
                .comptime_int => return true,
                else => unreachable,
            }
        },
        .float => return true,
        else => return @hasDecl(N, "is_signed"),
    }
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
    if (comptime !types.isNumeric(N))
        @compileError("zml.types.isUnsigned: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (comptime types.numericType(N)) {
        .bool => return true,
        .int => {
            switch (comptime @typeInfo(N)) {
                .int => |info| return info.signedness == .unsigned,
                .comptime_int => return false,
                else => unreachable,
            }
        },
        .float => return false,
        else => return @hasDecl(N, "is_unsigned"),
    }
}
