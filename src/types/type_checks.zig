//! Namespace for type checking functions.

const types = @import("../types.zig");

/// Checks if the input type is a supported type (numeric, vector, matrix,
/// array, or expression).
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a supported type, `false` otherwise.
pub fn isSupportedType(comptime T: type) bool {
    if (comptime types.isNumeric(T))
        return true;

    if (comptime isVector(T))
        return true;

    if (comptime isMatrix(T))
        return true;

    if (comptime isArray(T))
        return true;

    if (comptime isExpression(T))
        return true;

    return false;
}

/// Checks if the input type is a one-item pointer.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a one-item pointer, `false` otherwise.
pub fn isPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
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
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a many-item pointer, `false` otherwise.
pub fn isManyPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
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
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a constant one-item pointer, `false`
/// otherwise.
pub fn isConstPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
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
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a slice, `false` otherwise.
pub fn isSlice(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .slice) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input type is a simd vector.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a simd vector, `false` otherwise.
pub fn isSimdVector(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .vector => return true,
        else => return false,
    }
}

/// Checks if the input type is an instance of a vector.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a vector, `false` otherwise.
pub fn isVector(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a dense vector.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a dense vector, `false` otherwise.
pub fn isDenseVector(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a sparse vector.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a sparse vector, `false` otherwise.
pub fn isSparseVector(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a matrix, `false` otherwise.
pub fn isMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix"),
        else => return false,
    }
}

/// Checks if the input type is a square matrix, i.e., a matrix type that is
/// always square (symmetric or hermitian).
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a square matrix, `false` otherwise.
pub fn isSquareMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_symmetric") or @hasDecl(T, "is_hermitian"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a general matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a general matrix, `false` otherwise.
pub fn isGeneralMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a general dense matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a general dense matrix, `false` otherwise.
pub fn isGeneralDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a general sparse matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a general sparse matrix, `false` otherwise.
pub fn isGeneralSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a symmetric matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a symmetric matrix, `false` otherwise.
pub fn isSymmetricMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a symmetric dense matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a symmetric dense matrix, `false` otherwise.
pub fn isSymmetricDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a symmetric sparse matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a symmetric sparse matrix, `false` otherwise.
pub fn isSymmetricSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a hermitian matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a hermitian matrix, `false` otherwise.
pub fn isHermitianMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a hermitian dense matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a hermitian dense matrix, `false` otherwise.
pub fn isHermitianDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a hermitian sparse matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a hermitian sparse matrix, `false` otherwise.
pub fn isHermitianSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a triangular matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a triangular matrix, `false` otherwise.
pub fn isTriangularMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_triangular"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a triangular dense matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a triangular dense matrix, `false` otherwise.
pub fn isTriangularDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_triangular") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a triangular sparse matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a triangular sparse matrix, `false` otherwise.
pub fn isTriangularSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_triangular") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a diagonal matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a diagonal matrix, `false` otherwise.
pub fn isDiagonalMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_diagonal"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a permutation matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a permutation matrix, `false` otherwise.
pub fn isPermutationMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_permutation"),
        else => return false,
    }
}

/// Checks if the input type is a dense matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a dense matrix, `false` otherwise.
pub fn isDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is a sparse matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a sparse matrix, `false` otherwise.
pub fn isSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of an array.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is an array, `false` otherwise.
pub fn isArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a dense array.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a dense array, `false` otherwise.
pub fn isDenseArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a strided array.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a strided array, `false` otherwise.
pub fn isStridedArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and @hasDecl(T, "is_strided"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a sparse array.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a sparse array, `false` otherwise.
pub fn isSparseArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

pub fn isExpression(comptime T: type) bool {
    _ = T;
    return false;
    // return T == Expression;
}

/// Checks if the input numeric type is of fixed precision.
///
/// This function checks if the input numeric type is of fixed precision, which
/// includes types labeled as `bool`, `int`, `float`, `dyadic` and `cfloat`.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is of fixed precision, `false` otherwise.
pub fn isFixedPrecision(comptime N: type) bool {
    if (!types.isNumeric(N))
        @compileError("zml.types.isFixedPrecision: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (types.numericType(N)) {
        .bool => return true,
        .int => return true,
        .float => return true,
        .dyadic => return true,
        .cfloat => return true,
        else => return false,
    }
}

/// Checks if the input numeric type is of arbitrary precision.
///
/// This function checks if the input numeric type is of arbitrary precision,
/// which includes types labeled as `integer`, `rational`, `real` and `complex`.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is of arbitrary precision, `false` otherwise.
pub fn isArbitraryPrecision(comptime N: type) bool {
    if (!types.isNumeric(N))
        @compileError("zml.types.isArbitraryPrecision: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (types.numericType(N)) {
        .integer => return true,
        .rational => return true,
        .real => return true,
        .complex => return true,
        else => return false,
    }
}

/// Checks if the input numeric type is integral.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is integral, `false` otherwise.
pub fn isIntegral(comptime N: type) bool {
    if (!types.isNumeric(N))
        @compileError("zml.types.isIntegral: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (types.numericType(N)) {
        .bool => return true,
        .int => return true,
        .integer => return true,
        else => return false,
    }
}

/// Checks if the input numeric type is non-integral.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is non-integral, `false` otherwise.
pub fn isNonIntegral(comptime N: type) bool {
    if (!types.isNumeric(N))
        @compileError("zml.types.isNonIntegral: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (types.numericType(N)) {
        .float => return true,
        .dyadic => return true,
        .cfloat => return true,
        .rational => return true,
        .real => return true,
        .complex => return true,
        else => return false,
    }
}

/// Checks if the input numeric type is real.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is real, `false` otherwise.
pub fn isReal(comptime N: type) bool {
    if (!types.isNumeric(N))
        @compileError("zml.types.isReal: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (types.numericType(N)) {
        .bool => return true,
        .int => return true,
        .float => return true,
        .dyadic => return true,
        .integer => return true,
        .rational => return true,
        .real => return true,
        else => return false,
    }
}

/// Checks if the input numeric type is complex.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is complex, `false` otherwise.
pub fn isComplex(comptime N: type) bool {
    if (!types.isNumeric(N))
        @compileError("zml.types.isComplex: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (types.numericType(N)) {
        .cfloat => return true,
        .complex => return true,
        else => return false,
    }
}

/// Checks if the input numeric type is signed.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is signed, `false` otherwise.
pub fn isSigned(comptime N: type) bool {
    if (!types.isNumeric(N))
        @compileError("zml.types.isSigned: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (types.numericType(N)) {
        .int => {
            switch (@typeInfo(N)) {
                .int => |info| return info.signedness == .signed,
                .comptime_int => return true,
                else => unreachable,
            }
        },
        .float => return true,
        .dyadic => return true,
        .cfloat => return true,
        .integer => return true,
        .rational => return true,
        .real => return true,
        .complex => return true,
        else => return false,
    }
}

/// Checks if the input numeric type is unsigned.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is unsigned, `false` otherwise.
pub fn isUnsigned(comptime N: type) bool {
    if (!types.isNumeric(N))
        @compileError("zml.types.isUnsigned: " ++ @typeName(N) ++ " is not a supported numeric type");

    switch (types.numericType(N)) {
        .bool => return true,
        .int => {
            switch (@typeInfo(N)) {
                .int => |info| return info.signedness == .unsigned,
                .comptime_int => return false,
                else => unreachable,
            }
        },
        else => return false,
    }
}
