//! Namespace for type checking functions.

const types = @import("../types.zig");

/// Checks if the input type is a supported type (numeric, vector, matrix,
/// array, or expression).
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
pub fn isDenseVector(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

pub fn isSparseVector(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_vector") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a matrix, i.e., an instance of a
/// `matrix.General`, `matrix.Symmetric`, `matrix.Hermitian`,
/// `matrix.Triangular`, `matrix.Diagonal`, `matrix.Banded`,
/// `matrix.Tridiagonal`, or `matrix.sparse` matrix.
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

pub fn isSquareMatrix(comptime T: type) bool {
    return isSymmetricDenseMatrix(T) or isHermitianDenseMatrix(T) or
        isGeneralTridiagonalMatrix(T) or isSymmetricTridiagonalMatrix(T) or
        isHermitianTridiagonalMatrix(T) or isSymmetricSparseMatrix(T) or
        isHermitianSparseMatrix(T) or isSymmetricBlockMatrix(T) or
        isHermitianBlockMatrix(T) or isPermutationMatrix(T);
}

/// Checks if the input type is an instance of a `matrix.general.Dense`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.general.Dense`, `false` otherwise.
pub fn isGeneralDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.symmetric.Dense`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.symmetric.Dense`, `false` otherwise.
pub fn isSymmetricDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.hermitian.Dense`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.hermitian.Dense`, `false` otherwise.
pub fn isHermitianDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.triangular.Dense`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.triangular.Dense`, `false` otherwise.
pub fn isTriangularDenseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_triangular") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.general.Banded`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.general.Banded`, `false` otherwise.
pub fn isGeneralBandedMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_banded"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.symmetric.Banded`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.symmetric.Banded`, `false` otherwise.
pub fn isSymmetricBandedMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_banded"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.hermitian.Banded`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.hermitian.Banded`, `false` otherwise.
pub fn isHermitianBandedMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_banded"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.general.Tridiagonal`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.general.Tridiagonal`, `false` otherwise.
pub fn isGeneralTridiagonalMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_tridiagonal"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.symmetric.Tridiagonal`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.symmetric.Tridiagonal`, `false` otherwise.
pub fn isSymmetricTridiagonalMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_tridiagonal"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.hermitian.Tridiagonal`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.hermitian.Tridiagonal`, `false` otherwise.
pub fn isHermitianTridiagonalMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_tridiagonal"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.general.Sparse`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.general.Sparse`, `false` otherwise.
pub fn isGeneralSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.symmetric.Sparse`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.symmetric.Sparse`, `false` otherwise.
pub fn isSymmetricSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.hermitian.Sparse`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.hermitian.Sparse`, `false` otherwise.
pub fn isHermitianSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.triangular.Sparse`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.triangular.Sparse`, `false` otherwise.
pub fn isTriangularSparseMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_triangular") and @hasDecl(T, "is_sparse"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.general.Block`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.general.Block`, `false` otherwise.
pub fn isGeneralBlockMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_general") and @hasDecl(T, "is_block"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.symmetric.Block`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.symmetric.Block`, `false` otherwise.
pub fn isSymmetricBlockMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_symmetric") and @hasDecl(T, "is_block"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.hermitian.Block`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.hermitian.Block`, `false` otherwise.
pub fn isHermitianBlockMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_hermitian") and @hasDecl(T, "is_block"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.Diagonal`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.Diagonal`, `false` otherwise.
pub fn isDiagonalMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_diagonal"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.Permutation`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.Permutation`, `false` otherwise.
pub fn isPermutationMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => @hasDecl(T, "is_matrix") and @hasDecl(T, "is_permutation"),
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

/// Checks if the input type is a banded matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a banded matrix, `false` otherwise.
pub fn isBandedMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_banded"),
        else => return false,
    }
}

/// Checks if the input type is a tridiagonal matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a tridiagonal matrix, `false` otherwise.
pub fn isTridiagonalMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_tridiagonal"),
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

/// Checks if the input type is a block matrix.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a block matrix, `false` otherwise.
pub fn isBlockMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_matrix") and @hasDecl(T, "is_block"),
        else => return false,
    }
}

/// Checks if the input type is an instance of an array, i.e., an instance of a
/// `array.Dense`, `array.Strided`, or `array.Sparse` array.
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

/// Checks if the input type is an instance of a `array.Dense`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `array.Dense`, `false` otherwise.
pub fn isDenseArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and @hasDecl(T, "is_dense"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `array.Strided`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `array.Strided`, `false` otherwise.
pub fn isStridedArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_array") and @hasDecl(T, "is_strided"),
        else => return false,
    }
}

/// Checks if the input type is an instance of a `array.Sparse`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `array.Sparse`, `false` otherwise.
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
/// This function checks if the input type is a fixed precision numeric type,
/// which includes types labeled as `bool`, `int`, `float`, and `cfloat`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is of fixed precision, `false` otherwise.
pub fn isFixedPrecision(comptime T: type) bool {
    switch (types.numericType(T)) {
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
/// This function checks if the input type is an arbitrary precision numeric
/// type, which includes types labeled as `integer`, `rational`, `real`,
/// `complex`, and `expression`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is of arbitrary precision, `false` otherwise.
pub fn isArbitraryPrecision(comptime T: type) bool {
    switch (types.numericType(T)) {
        .integer => return true,
        .rational => return true,
        .real => return true,
        .complex => return true,
        else => return false,
    }
}

/// Checks if the input type is integral.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is integral, `false` otherwise.
pub fn isIntegral(comptime T: type) bool {
    switch (types.numericType(T)) {
        .bool => return true,
        .int => return true,
        .integer => return true,
        else => return false,
    }
}

/// Checks if the input type is non-integral.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is non-integral, `false` otherwise.
pub fn isNonIntegral(comptime T: type) bool {
    switch (types.numericType(T)) {
        .float => return true,
        .dyadic => return true,
        .cfloat => return true,
        .rational => return true,
        .real => return true,
        .complex => return true,
        else => return false,
    }
}

/// Checks if the input type is real.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is real, `false` otherwise.
pub fn isReal(comptime T: type) bool {
    switch (types.numericType(T)) {
        .float => return true,
        .dyadic => return true,
        .integer => return true,
        .rational => return true,
        .real => return true,
        else => return false,
    }
}

/// Checks if the input type is complex.
///
/// This function checks if the input type is a complex numeric type, which
/// includes types labeled as `cfloat`, `complex`, and `expression`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is complex, `false` otherwise.
pub fn isComplex(comptime T: type) bool {
    switch (types.numericType(T)) {
        .cfloat => return true,
        .complex => return true,
        else => return false,
    }
}

/// Checks if the input type is signed.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is signed, `false` otherwise.
pub fn isSigned(comptime T: type) bool {
    switch (types.numericType(T)) {
        .int => {
            switch (@typeInfo(T)) {
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

/// Checks if the input type is unsigned.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is unsigned, `false` otherwise.
pub fn isUnsigned(comptime T: type) bool {
    switch (types.numericType(T)) {
        .bool => return true,
        .int => {
            switch (@typeInfo(T)) {
                .int => |info| return info.signedness == .unsigned,
                .comptime_int => return false,
                else => unreachable,
            }
        },
        else => return false,
    }
}
