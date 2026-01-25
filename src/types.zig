//! Namespace for type definitions and utilities.

const std = @import("std");

const constants = @import("constants.zig");

pub const default_uint = u32;
pub const default_int = i32;
pub const default_float = f64;
pub const default_layout = Layout.col_major;
pub const default_uplo = Uplo.upper;
pub const default_diag = Diag.non_unit;

pub const standard_integer_types: [10]type = .{
    u8,   u16,
    u32,  u64,
    u128, i8,
    i16,  i32,
    i64,  i128,
};

const dyadic = @import("dyadic.zig");
const Dyadic = dyadic.Dyadic;
const cfloat = @import("cfloat.zig");
const Cfloat = cfloat.Cfloat;
const cf16 = @import("cfloat.zig").cf16;
const cf32 = @import("cfloat.zig").cf32;
const cf64 = @import("cfloat.zig").cf64;
const cf80 = @import("cfloat.zig").cf80;
const cf128 = @import("cfloat.zig").cf128;
const comptime_cfloat = @import("cfloat.zig").comptime_cfloat;
const integer = @import("integer.zig");
const Integer = integer.Integer;
const rational = @import("rational.zig");
const Rational = rational.Rational;
const real = @import("real.zig");
const Real = real.Real;
const complex = @import("complex.zig");
const Complex = complex.Complex;

const vector = @import("vector.zig");
const matrix = @import("matrix.zig");
const array = @import("array.zig");
const Expression = @import("expression.zig").Expression;

pub const Cmp = enum(u2) {
    gt,
    eq,
    lt,

    pub inline fn invert(self: Cmp) Cmp {
        return switch (self) {
            .gt => .lt,
            .eq => .eq,
            .lt => .gt,
        };
    }
};

pub const Layout = enum(u1) {
    row_major,
    col_major,

    pub inline fn toCUInt(self: Layout) c_uint {
        return switch (self) {
            .row_major => 101,
            .col_major => 102,
        };
    }

    pub inline fn toCInt(self: Layout) c_int {
        return switch (self) {
            .row_major => 101,
            .col_major => 102,
        };
    }

    pub inline fn toIterationOrder(self: Layout) IterationOrder {
        return switch (self) {
            .row_major => .right_to_left,
            .col_major => .left_to_right,
        };
    }

    pub inline fn invert(self: Layout) Layout {
        return switch (self) {
            .row_major => .col_major,
            .col_major => .row_major,
        };
    }

    pub fn resolve3(self: Layout, other1: Layout, other2: Layout) Layout {
        if (self == other1 and self == other2)
            return self;

        if (self == other1 or self == other2)
            return self;

        if (other1 == other2)
            return other1;

        return default_layout;
    }
};

pub const Uplo = enum(u1) {
    upper,
    lower,

    pub inline fn toCUInt(self: Uplo) c_uint {
        return switch (self) {
            .upper => 121,
            .lower => 122,
        };
    }

    pub inline fn toCInt(self: Uplo) c_int {
        return switch (self) {
            .upper => 121,
            .lower => 122,
        };
    }

    pub inline fn toChar(self: Uplo) u8 {
        return switch (self) {
            .upper => 'U',
            .lower => 'L',
        };
    }

    pub inline fn invert(self: Uplo) Uplo {
        return switch (self) {
            .upper => .lower,
            .lower => .upper,
        };
    }
};

pub const Diag = enum(u1) {
    non_unit,
    unit,

    pub inline fn toCUInt(self: Diag) c_uint {
        return switch (self) {
            .non_unit => 131,
            .unit => 132,
        };
    }

    pub inline fn toCInt(self: Diag) c_int {
        return switch (self) {
            .non_unit => 131,
            .unit => 132,
        };
    }

    pub inline fn toChar(self: Diag) u8 {
        return switch (self) {
            .non_unit => 'N',
            .unit => 'U',
        };
    }

    pub inline fn invert(self: Diag) Diag {
        return switch (self) {
            .non_unit => .unit,
            .unit => .non_unit,
        };
    }
};

pub const IterationOrder = enum {
    left_to_right,
    right_to_left,
};

/// `types.NumericType` is an enum that represents the different numeric types
/// supported by the library. It is used to categorize types based on their
/// properties and capabilities, such as whether they are integers, floats,
/// complex numbers, etc.
///
/// This enum is used in various places in the library to determine how to
/// handle different types of numeric data. It allows for type checking and
/// coercion between different numeric types, ensuring that operations are
/// performed correctly and efficiently.
///
/// ## Values
/// * `bool`: Represents the boolean type (`bool`).
/// * `int`: Represents integer types:
///   * `usize`, `u8`, `u16`, `u32`, `u64`, `u128`
///   * `isize`, `i8`, `i16`, `i32`, `i64`, `i128`
///   * `uX`, `iX` (where X is any bit size)
///   * `comptime_int`
/// * `float`: Represents floating*point types:
///   * `f16`, `f32`, `f64`, `f80`, `f128`
///   * `comptime_float`
/// * `dyadic`: Represents dyadic rational types (`Dyadic`).
/// * `cfloat`: Represents complex floating*point types:
///   * `cf16`, `cf32`, `cf64`, `cf80`, `cf128`
///   * `Cfloat(Dyadic(...))`
///   * `comptime_cfloat`
/// * `integer`: Represents the arbitrary precision integer type (`Integer`).
/// * `rational`: Represents the arbitrary precision rational type (`Rational`).
/// * `real`: Represents the arbitrary precision real type (`Real`).
/// * `complex`: Represents complex arbitrary precision types:
///   * `Complex(Rational)`
///   * `Complex(Real)`
pub const NumericType = enum {
    bool,
    int,
    float,
    dyadic,
    cfloat,
    integer,
    rational,
    real,
    complex,
    custom,

    pub fn lt(self: NumericType, other: NumericType) bool {
        return @intFromEnum(self) < @intFromEnum(other);
    }

    pub fn le(self: NumericType, other: NumericType) bool {
        return @intFromEnum(self) <= @intFromEnum(other);
    }

    pub fn gt(self: NumericType, other: NumericType) bool {
        return @intFromEnum(self) > @intFromEnum(other);
    }

    pub fn ge(self: NumericType, other: NumericType) bool {
        return @intFromEnum(self) >= @intFromEnum(other);
    }
};

pub const VectorType = enum {
    dense,
    sparse,
    custom,
    numeric, // Fallback for numeric types that are not vectors
};

pub const MatrixType = enum {
    general_dense,
    general_sparse,
    symmetric_dense,
    symmetric_sparse,
    hermitian_dense,
    hermitian_sparse,
    triangular_dense,
    triangular_sparse,
    diagonal,
    permutation,
    custom,
    numeric, // Fallback for numeric types that are not matrices
};

pub const MatrixKind = enum {
    general,
    symmetric,
    hermitian,
    triangular,
    diagonal,
    permutation,
};

pub const MatrixStorage = enum {
    dense,
    sparse,
};

pub const ArrayType = enum {
    dense,
    strided,
    sparse,
    custom,
    numeric, // Fallback for numeric types that are not arrays
};

pub const Domain = enum {
    numeric,
    vector,
    matrix,
    array,
    expression,
};

/// A useless allocator that does nothing, always signalling allocation failure.
pub const useless_allocator: std.mem.Allocator = .{
    .ptr = undefined,
    .vtable = &vtable,
};

const vtable: std.mem.Allocator.VTable = .{
    .alloc = alloc,
    .resize = resize,
    .remap = remap,
    .free = free,
};

fn alloc(context: *anyopaque, len: usize, alignment: std.mem.Alignment, ra: usize) ?[*]u8 {
    _ = context;
    _ = len;
    _ = alignment;
    _ = ra;

    return null;
}

fn resize(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ra: usize) bool {
    _ = context;
    _ = memory;
    _ = alignment;
    _ = new_len;
    _ = ra;

    return false;
}

fn remap(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ra: usize) ?[*]u8 {
    _ = context;
    _ = memory;
    _ = alignment;
    _ = new_len;
    _ = ra;

    return null;
}

fn free(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, ra: usize) void {
    _ = context;
    _ = memory;
    _ = alignment;
    _ = ra;

    return;
}

pub fn empty(comptime T: type) T {
    if (comptime !isSupportedType(T))
        @compileError("zml.types.empty: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .numeric => switch (comptime numericType(T)) {
            .bool => return false,
            .int => return 0,
            .float => return 0.0,
            .dyadic => return .zero,
            .cfloat => return .{ .re = 0.0, .im = 0.0 },
            .integer => return .empty,
            .rational => return .empty,
            .real => return .empty,
            .complex => return .empty,
            .custom => {
                if (comptime !@hasDecl(T, "empty"))
                    @compileError("zml.types.empty: custom numeric type " ++ @typeName(T) ++ " must have an `empty` declaration");

                return .empty;
            },
        },
        .vector => return .empty,
        .matrix => return .empty,
        .array => return .empty,
        .expression => return .empty,
    }
}

/// Checks the the input type `N` and returns the corresponding
/// `types.NumericType`.
///
/// Checks that the input type is a supported numeric type and returns the
/// corresponding `types.NumericType` enum value. If the type is not supported,
/// it will raise a compile error.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check.
///
/// ## Returns
/// `types.NumericType`: The corresponding `types.NumericType` enum value.
pub inline fn numericType(comptime N: type) NumericType {
    // Without inline, functions calling this fail miserably. I have no idea why.
    switch (comptime @typeInfo(N)) {
        .bool => return .bool,
        .int, .comptime_int => return .int,
        .float, .comptime_float => return .float,
        else => {
            if (comptime @hasDecl(N, "is_dyadic"))
                return .dyadic;

            if (comptime (@hasDecl(N, "is_cfloat")) or
                N == std.math.Complex(f16) or N == std.math.Complex(f32) or N == std.math.Complex(f64) or
                N == std.math.Complex(f80) or N == std.math.Complex(f128) or N == std.math.Complex(comptime_float))
                return .cfloat;

            if (comptime @hasDecl(N, "is_integer"))
                return .integer;

            if (comptime @hasDecl(N, "is_rational"))
                return .rational;

            if (comptime @hasDecl(N, "is_real"))
                return .real;

            if (comptime @hasDecl(N, "is_complex"))
                return .complex;

            if (comptime @hasDecl(N, "is_numeric"))
                return .custom;

            @compileError("zml.types.numericType: " ++ @typeName(N) ++ " is not a supported numeric type");
        },
    }
}

/// Determines the vector type of the input type `V`.
///
/// ## Arguments
/// * `V` (`comptime type`): The type to check.
///
/// ## Returns
/// `types.VectorType`: The corresponding `types.VectorType` enum value.
pub inline fn vectorType(comptime V: type) VectorType {
    if (comptime isDenseVector(V))
        return .dense;

    if (comptime isSparseVector(V))
        return .sparse;

    if (comptime isCustomVector(V))
        return .custom;

    return .numeric; // Fallback for numeric types that are not vectors
}

/// Determines the matrix type of the input type `M`.
///
/// ## Arguments
/// * `M` (`comptime type`): The type to check.
///
/// ## Returns
/// `types.MatrixType`: The corresponding `types.MatrixType` enum value.
pub inline fn matrixType(comptime M: type) MatrixType {
    if (comptime isGeneralDenseMatrix(M))
        return .general_dense;

    if (comptime isSymmetricDenseMatrix(M))
        return .symmetric_dense;

    if (comptime isHermitianDenseMatrix(M))
        return .hermitian_dense;

    if (comptime isTriangularDenseMatrix(M))
        return .triangular_dense;

    if (comptime isGeneralSparseMatrix(M))
        return .general_sparse;

    if (comptime isSymmetricSparseMatrix(M))
        return .symmetric_sparse;

    if (comptime isHermitianSparseMatrix(M))
        return .hermitian_sparse;

    if (comptime isTriangularSparseMatrix(M))
        return .triangular_sparse;

    if (comptime isDiagonalMatrix(M))
        return .diagonal;

    if (comptime isPermutationMatrix(M))
        return .permutation;

    if (comptime isCustomMatrix(M))
        return .custom;

    return .numeric; // Fallback for numeric types that are not matrices
}

/// Determines the array type of the input type `A`.
///
/// ## Arguments
/// * `A` (`comptime type`): The type to check.
///
/// ## Returns
/// `types.ArrayType`: The corresponding `types.ArrayType` enum value.
pub inline fn arrayType(comptime A: type) ArrayType {
    if (comptime isDenseArray(A))
        return .dense;

    if (comptime isStridedArray(A))
        return .strided;

    if (comptime isSparseArray(A))
        return .sparse;

    if (comptime isCustomArray(A))
        return .custom;

    return .numeric; // Fallback for numeric types that are not arrays
}

/// Determines the domain of the input type `T`.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `types.Domain`: The corresponding `types.Domain` enum value.
pub inline fn domain(comptime T: type) Domain {
    if (comptime isNumeric(T))
        return .numeric;

    if (comptime isVector(T))
        return .vector;

    if (comptime isMatrix(T))
        return .matrix;

    if (comptime isArray(T))
        return .array;

    if (comptime isExpression(T))
        return .expression;

    @compileError("zml.types.domain: " ++ @typeName(T) ++ " does not belong to any supported domain");
}

const type_checks = @import("types/type_checks.zig");
pub const isSupportedType = type_checks.isSupportedType;
pub const isPointer = type_checks.isPointer;
pub const isManyPointer = type_checks.isManyPointer;
pub const isConstPointer = type_checks.isConstPointer;
pub const isSlice = type_checks.isSlice;
pub const isSimdVector = type_checks.isSimdVector;
pub const isNumeric = type_checks.isNumeric;
pub const isVector = type_checks.isVector;
pub const isDenseVector = type_checks.isDenseVector;
pub const isSparseVector = type_checks.isSparseVector;
pub const isCustomVector = type_checks.isCustomVector;
pub const isMatrix = type_checks.isMatrix;
pub const isSquareMatrix = type_checks.isSquareMatrix;
pub const isGeneralDenseMatrix = type_checks.isGeneralDenseMatrix;
pub const isSymmetricDenseMatrix = type_checks.isSymmetricDenseMatrix;
pub const isHermitianDenseMatrix = type_checks.isHermitianDenseMatrix;
pub const isTriangularDenseMatrix = type_checks.isTriangularDenseMatrix;
pub const isGeneralSparseMatrix = type_checks.isGeneralSparseMatrix;
pub const isSymmetricSparseMatrix = type_checks.isSymmetricSparseMatrix;
pub const isHermitianSparseMatrix = type_checks.isHermitianSparseMatrix;
pub const isTriangularSparseMatrix = type_checks.isTriangularSparseMatrix;
pub const isDiagonalMatrix = type_checks.isDiagonalMatrix;
pub const isPermutationMatrix = type_checks.isPermutationMatrix;
pub const isGeneralMatrix = type_checks.isGeneralMatrix;
pub const isSymmetricMatrix = type_checks.isSymmetricMatrix;
pub const isHermitianMatrix = type_checks.isHermitianMatrix;
pub const isTriangularMatrix = type_checks.isTriangularMatrix;
pub const isDenseMatrix = type_checks.isDenseMatrix;
pub const isSparseMatrix = type_checks.isSparseMatrix;
pub const isCustomMatrix = type_checks.isCustomMatrix;
pub const isArray = type_checks.isArray;
pub const isDenseArray = type_checks.isDenseArray;
pub const isStridedArray = type_checks.isStridedArray;
pub const isSparseArray = type_checks.isSparseArray;
pub const isCustomArray = type_checks.isCustomArray;
pub const isExpression = type_checks.isExpression;
pub const isAllocated = type_checks.isAllocated;
pub const isIntegral = type_checks.isIntegral;
pub const isNonIntegral = type_checks.isNonIntegral;
pub const isRealType = type_checks.isRealType;
pub const isComplexType = type_checks.isComplexType;
pub const isSigned = type_checks.isSigned;
pub const isUnsigned = type_checks.isUnsigned;

const coercion = @import("types/coercion.zig");
pub const Coerce = coercion.Coerce;
pub const MulCoerce = coercion.MulCoerce;
pub const canCoerce = coercion.canCoerce;
pub const EnsureDomain = coercion.EnsureDomain;
pub const EnsureVector = coercion.EnsureVector;
pub const EnsureMatrix = coercion.EnsureMatrix;
pub const EnsureArray = coercion.EnsureArray;
pub const EnsureFloat = coercion.EnsureFloat;

const casting = @import("types/casting.zig");
pub const scast = casting.scast;
pub const cast = casting.cast;

/// Returns the input type as is, without any modifications.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to return.
///
/// ## Returns
/// `type`: The input type `T`.
pub fn Identity(comptime T: type) type {
    return T;
}

/// Returns the scalar type of a given numeric type, vector, matrix or array.
///
/// This function returns the scalar type of a given numeric type, vector,
/// matrix, or array. If the input type is a vector, a matrix or an array, it
/// returns the element type (equivalent to `types.Numeric`). If the input type
/// is a numeric type, it returns the type itself, unless it is a complex type,
/// in which case it returns the scalar type of the complex type.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to get the scalar type of. Must be a
/// supported numeric type, vector, matrix, or array.
///
/// ## Returns
/// `type`: The scalar type of the input type.
pub fn Scalar(comptime T: type) type {
    if (comptime !isSupportedType(T))
        @compileError("zml.types.Scalar: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .numeric => switch (comptime numericType(T)) {
            .bool => return T,
            .int => return T,
            .float => return T,
            .dyadic => return T,
            .cfloat => switch (T) {
                std.math.Complex(f16) => return f16,
                std.math.Complex(f32) => return f32,
                std.math.Complex(f64) => return f64,
                std.math.Complex(f80) => return f80,
                std.math.Complex(f128) => return f128,
                std.math.Complex(comptime_float) => return comptime_float,
                else => return T.Scalar,
            },
            .integer => return Integer,
            .rational => return Rational,
            .real => return Real,
            .complex => return T.Scalar,
            .custom => {
                if (comptime !@hasDecl(T, "Scalar"))
                    @compileError("zml.types.Scalar: custom numeric type " ++ @typeName(T) ++ " must have a `Scalar` declaration");

                return T.Scalar;
            },
        },
        .vector => return Numeric(T),
        .matrix => return Numeric(T),
        .array => return Numeric(T),
        .expression => return Expression,
    }
}

/// Returns the underlying numeric type of a given numeric type, vector, matrix
/// or array.
///
/// This function returns the underlying numeric type of a given numeric type,
/// vector, matrix or array. If the input type is a vector, matrix or array, it
/// returns the element type. If the input type is a numeric type, it returns
/// the type itself.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to get the numeric type of. Must be a
/// supported numeric type, vector, matrix, or array.
///
/// ## Returns
/// `type`: The underlying numeric type of the input type.
pub fn Numeric(comptime T: type) type {
    if (comptime !isSupportedType(T))
        @compileError("zml.types.Numeric: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .numeric => return T,
        .vector => switch (comptime vectorType(T)) {
            .dense => return T.Numeric,
            .sparse => return T.Numeric,
            .custom => {
                if (comptime !@hasDecl(T, "Numeric"))
                    @compileError("zml.types.Numeric: custom vector type " ++ @typeName(T) ++ " must have a `Numeric` declaration");

                return T.Numeric;
            },
            .numeric => return T,
        },
        .matrix => switch (comptime matrixType(T)) {
            .general_dense => return T.Numeric,
            .general_sparse => return T.Numeric,
            .symmetric_dense => return T.Numeric,
            .symmetric_sparse => return T.Numeric,
            .hermitian_dense => return T.Numeric,
            .hermitian_sparse => return T.Numeric,
            .triangular_dense => return T.Numeric,
            .triangular_sparse => return T.Numeric,
            .diagonal => return T.Numeric,
            .permutation => return T.Numeric,
            .custom => {
                if (comptime !@hasDecl(T, "Numeric"))
                    @compileError("zml.types.Numeric: custom matrix type " ++ @typeName(T) ++ " must have a `Numeric` declaration");

                return T.Numeric;
            },
            .numeric => return T,
        },
        .array => switch (comptime arrayType(T)) {
            .dense => return T.Numeric,
            .strided => return T.Numeric,
            .sparse => return T.Numeric,
            .custom => {
                if (comptime !@hasDecl(T, "Numeric"))
                    @compileError("zml.types.Numeric: custom array type " ++ @typeName(T) ++ " must have a `Numeric` declaration");

                return T.Numeric;
            },
            .numeric => return T,
        },
        .expression => return Expression,
    }
}

pub fn layoutOf(comptime T: type) Layout {
    if (comptime !isSupportedType(T))
        @compileError("zml.types.layoutOf: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .matrix => switch (comptime matrixType(T)) {
            .general_dense => return T.storage_layout,
            .general_sparse => return T.storage_layout,
            .symmetric_dense => return T.storage_layout,
            .symmetric_sparse => return T.storage_layout,
            .hermitian_dense => return T.storage_layout,
            .hermitian_sparse => return T.storage_layout,
            .triangular_dense => return T.storage_layout,
            .triangular_sparse => return T.storage_layout,
            .diagonal => return T.storage_layout,
            .permutation => return T.storage_layout,
            .custom => {
                if (comptime !@hasDecl(T, "storage_layout"))
                    @compileError("zml.types.layoutOf: custom matrix type " ++ @typeName(T) ++ " must have a `storage_layout` declaration");

                return T.storage_layout;
            },
            .numeric => unreachable,
        },
        .array => switch (comptime arrayType(T)) {
            .dense => return T.storage_layout,
            .strided => return T.storage_layout,
            .sparse => return T.storage_layout,
            .custom => {
                if (comptime !@hasDecl(T, "storage_layout"))
                    @compileError("zml.types.layoutOf: custom array type " ++ @typeName(T) ++ " must have a `storage_layout` declaration");

                return T.storage_layout;
            },
            .numeric => unreachable,
        },
        else => @compileError("zml.types.layoutOf: T must be a matrix or array type, got " ++ @typeName(T)),
    }
}

pub fn uploOf(comptime T: type) Uplo {
    if (comptime !isSupportedType(T))
        @compileError("zml.types.uploOf: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .matrix => switch (comptime matrixType(T)) {
            .general_dense => return T.storage_uplo,
            .general_sparse => return T.storage_uplo,
            .symmetric_dense => return T.storage_uplo,
            .symmetric_sparse => return T.storage_uplo,
            .hermitian_dense => return T.storage_uplo,
            .hermitian_sparse => return T.storage_uplo,
            .triangular_dense => return T.storage_uplo,
            .triangular_sparse => return T.storage_uplo,
            .diagonal => return default_uplo,
            .permutation => return default_uplo,
            .custom => {
                if (comptime !@hasDecl(T, "storage_uplo"))
                    @compileError("zml.types.uploOf: custom matrix type " ++ @typeName(T) ++ " must have a `storage_uplo` declaration");

                return T.storage_uplo;
            },
            .numeric => unreachable,
        },
        else => @compileError("zml.types.uploOf: T must be a matrix type, got " ++ @typeName(T)),
    }
}

pub fn diagOf(comptime T: type) Diag {
    if (comptime !isSupportedType(T))
        @compileError("zml.types.diagOf: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .matrix => switch (comptime matrixType(T)) {
            .general_dense => return T.storage_diag,
            .general_sparse => return T.storage_diag,
            .symmetric_dense => return T.storage_diag,
            .symmetric_sparse => return T.storage_diag,
            .hermitian_dense => return T.storage_diag,
            .hermitian_sparse => return T.storage_diag,
            .triangular_dense => return T.storage_diag,
            .triangular_sparse => return T.storage_diag,
            .diagonal => return default_diag,
            .permutation => return default_diag,
            .custom => {
                if (comptime !@hasDecl(T, "storage_diag"))
                    @compileError("zml.types.diagOf: custom matrix type " ++ @typeName(T) ++ " must have a `storage_diag` declaration");

                return T.storage_diag;
            },
            .numeric => unreachable,
        },
        else => @compileError("zml.types.diagOf: T must be a matrix type, got " ++ @typeName(T)),
    }
}

/// Returns the pointer child type of a given pointer type, or the type itself
/// if the input type is not a pointer.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to get the child type of. Must be a
///   pointer type.
///
/// ## Returns
/// `type`: The child type of the input pointer type.
pub fn Child(comptime T: type) type {
    switch (comptime @typeInfo(T)) {
        .pointer => |info| {
            return info.child;
        },
        .vector => |info| {
            return info.child;
        },
        else => return T,
    }
}

const context_checks = @import("types/context_checks.zig");
pub const validateContext = context_checks.validateContext;
pub const partialValidateContext = context_checks.partialValidateContext;
pub const ctxHasField = context_checks.ctxHasField;
pub const getFieldOrDefault = context_checks.getFieldOrDefault;
pub const MixStructFields = context_checks.MixStructFields;
pub const mixStructFields = context_checks.mixStructFields;
pub const StripStructFields = context_checks.StripStructFields;
pub const stripStructFields = context_checks.stripStructFields;
pub const RenameStructFields = context_checks.RenameStructFields;
pub const renameStructFields = context_checks.renameStructFields;
pub const KeepStructFields = context_checks.KeepStructFields;
pub const keepStructFields = context_checks.keepStructFields;

pub fn ReturnTypeFromInputs(
    comptime func: anytype,
    comptime input_types: []const type,
) type {
    comptime var inputs: std.meta.Tuple(input_types) = undefined;

    inline for (input_types, 0..) |input_type, i| {
        inputs[i] = if (input_type == std.mem.Allocator)
            useless_allocator
        else
            empty(input_type);
    }

    switch (comptime input_types.len) {
        0 => return @TypeOf(func()),
        1 => return @TypeOf(func(inputs[0])),
        2 => return @TypeOf(func(inputs[0], inputs[1])),
        3 => return @TypeOf(func(inputs[0], inputs[1], inputs[2])),
        4 => return @TypeOf(func(inputs[0], inputs[1], inputs[2], inputs[3])),
        5 => return @TypeOf(func(inputs[0], inputs[1], inputs[2], inputs[3], inputs[4])),
        else => @compileError("zml.types.ReturnTypeFromInputs: functions with more than 5 parameters are not supported"),
    }
}

/// Checks if the type `T` has a method with the given name and type. `anytype`
/// parameters are counted as matching any type.
///
/// For  allocated numeric types, the method is expected to have an allocator as
/// first parameter and return an error union, even though they are not included
/// in `method_type`. For example, for a `method_type` of
/// `fn (a: i32, b: i32) i32`, the actual method signature should be
/// `fn (allocator: std.mem.Allocator, a: i32, b: i32) !i32`.
///
/// If the method has no parameters (i.e., is a constant method), the allocator
/// parameter must be optional.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
/// * `method_name` (`comptime []const u8`): The name of the method to check.
/// * `method_type` (`comptime type`): The expected type of the method,
///   including parameters and return type. Must not include the allocator
///   parameter for allocated types, and the return type must not be an error
///   union.
/// * `input_types` (`comptime []const type`): The types of the inputs to the
///   method, used to determine the return type for methods with inferred return
///   types.
///
/// ## Returns
/// `bool`: `true` if the method exists and has the correct type, `false`
/// otherwise.
pub fn hasMethod(
    comptime T: type,
    comptime method_name: []const u8,
    comptime method_type: type,
    comptime input_types: []const type,
) bool {
    if (comptime !isSupportedType(T))
        @compileError("zml.types.hasMethod: " ++ @typeName(T) ++ " is not a supported type");

    if (comptime !@hasDecl(T, method_name))
        return false;

    // Test that the method has the correct type
    const info_spec = @typeInfo(method_type);
    if (comptime info_spec != .@"fn")
        @compileError("zml.types.hasMethod: method_type must be a function type");

    const info_method = @typeInfo(@TypeOf(@field(T, method_name)));
    if (comptime info_method != .@"fn")
        return false;

    const spec_params = info_spec.@"fn".params;
    const spec_return = info_spec.@"fn".return_type.?;
    comptime var method_params = info_method.@"fn".params;
    comptime var method_return = if (info_method.@"fn".return_type) |r|
        r
    else
        ReturnTypeFromInputs(@field(T, method_name), input_types);
    const is_constant_method = spec_params.len == 0;

    if (comptime std.mem.eql(u8, method_name, "deinit")) {
        // Special case for deinit methods: they must have allocator as the second
        // parameter and return void
        if (comptime method_params.len != spec_params.len + 1)
            return false;

        if (method_params[0].type.? != spec_params[0].type.?)
            return false;

        if (method_params[1].type.? != std.mem.Allocator)
            return false;

        if (method_return != void)
            return false;

        return true;
    }

    switch (comptime domain(T)) {
        .numeric => {
            if (comptime isAllocated(T)) {
                // Allocated numeric types have an allocator as first parameter
                if (method_params.len != spec_params.len + 1)
                    return false;

                if (is_constant_method) {
                    // For constant methods, the allocator is optional
                    if (method_params[0].type.? != ?std.mem.Allocator)
                        return false;
                } else {
                    if (method_params[0].type.? != std.mem.Allocator)
                        return false;
                }

                // Remove the allocator parameter for comparison
                method_params = method_params[1..];
            } else {
                if (method_params.len != spec_params.len)
                    return false;
            }

            // Check parameter types
            inline for (spec_params, method_params) |spec_param, method_param| {
                if (comptime spec_param.is_generic or method_param.is_generic)
                    continue;

                if (comptime spec_param.type.? != method_param.type.?)
                    return false;
            }

            // Check return type
            if (comptime isAllocated(T)) {
                // Allocated numeric types return error unions
                const return_info = @typeInfo(method_return);
                if (return_info != .error_union)
                    return false;

                method_return = return_info.error_union.payload;
            }

            if (comptime spec_return != method_return)
                return false;

            return true;
        },
        else => @compileError("zml.types.hasMethod: only implemented for numeric types so far"),
    }
}
