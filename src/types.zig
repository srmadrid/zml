//! Namespace for type definitions and utilities.

const std = @import("std");

const constants = @import("constants.zig");

pub const default_uint = u32;
pub const default_int = i32;
pub const default_float = f64;

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

pub const Order = enum(u1) {
    row_major,
    col_major,

    pub inline fn toCUInt(self: Order) c_uint {
        return switch (self) {
            .row_major => 101,
            .col_major => 102,
        };
    }

    pub inline fn toCInt(self: Order) c_int {
        return switch (self) {
            .row_major => 101,
            .col_major => 102,
        };
    }

    pub inline fn toIterationOrder(self: Order) IterationOrder {
        return switch (self) {
            .row_major => .right_to_left,
            .col_major => .left_to_right,
        };
    }

    pub inline fn invert(self: Order) Order {
        return switch (self) {
            .row_major => .col_major,
            .col_major => .row_major,
        };
    }

    pub fn resolve3(self: Order, other1: Order, other2: Order) Order {
        if (self == other1 and self == other2)
            return self;

        if (self == other1 or self == other2)
            return self;

        if (other1 == other2)
            return other1;

        return .col_major; // default order
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

/// `NumericType` is an enum that represents the different numeric types
/// supported by the library. It is used to categorize types based on their
/// properties and capabilities, such as whether they are integers, floats,
/// complex numbers, etc.
///
/// This enum is used in various places in the library to determine how to
/// handle different types of numeric data. It allows for type checking and
/// coercion between different numeric types, ensuring that operations are
/// performed correctly and efficiently.
///
/// Values
/// ------
/// - `bool`: Represents the boolean type (`bool`).
/// - `int`: Represents integer types:
///   - `usize`, `u8`, `u16`, `u32`, `u64`, `u128`
///   - `isize`, `i8`, `i16`, `i32`, `i64`, `i128`
///   - `comptime_int`
/// - `float`: Represents floating-point types:
///   - `f16`, `f32`, `f64`, `f80`, `f128`
///   - `comptime_float`
/// - `dyadic`: Represents dyadic rational types (`Dyadic`).
/// - `cfloat`: Represents complex floating-point types:
///   - `Cfloat(Dyadic(...))`
///   - `cf16`, `cf32`, `cf64`, `cf80`, `cf128`
///   - `comptime_cfloat`
/// - `integer`: Represents arbitrary precision integer type (`Integer`).
/// - `rational`: Represents arbitrary precision rational type (`Rational`).
/// - `real`: Represents arbitrary precision real type (`Real`).
/// - `complex`: Represents complex arbitrary precision types:
///   - `Complex(Rational)`
///   - `Complex(Real)`
/// - `expression`: Represents symbolic expressions (Expression).
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

    pub fn match(comptime self: NumericType) fn (type) bool {
        return switch (self) {
            .bool => {
                const tmp = struct {
                    fn isBool(comptime T: type) bool {
                        switch (@typeInfo(T)) {
                            .bool => return true,
                            else => return false,
                        }
                    }
                };
                return tmp.isBool;
            },
            .int => {
                const tmp = struct {
                    fn isInt(comptime T: type) bool {
                        switch (@typeInfo(T)) {
                            .int, .comptime_int => return true,
                            else => return false,
                        }
                    }
                };
                return tmp.isInt;
            },
            .float => {
                const tmp = struct {
                    fn isFloat(comptime T: type) bool {
                        switch (@typeInfo(T)) {
                            .float, .comptime_float => return true,
                            else => return false,
                        }
                    }
                };
                return tmp.isFloat;
            },
            .dyadic => {
                const tmp = struct {
                    fn isDyadic(comptime T: type) bool {
                        return @hasDecl(T, "is_dyadic");
                    }
                };
                return tmp.isDyadic;
            },
            .cfloat => {
                const tmp = struct {
                    fn isCfloat(comptime T: type) bool {
                        if ((@hasDecl(T, "is_cfloat")) or
                            T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or
                            T == std.math.Complex(f80) or T == std.math.Complex(f128) or T == std.math.Complex(comptime_float))
                            return true;

                        return false;
                    }
                };
                return tmp.isCfloat;
            },
            .integer => {
                const tmp = struct {
                    fn isInteger(comptime T: type) bool {
                        return T == Integer;
                    }
                };
                return tmp.isInteger;
            },
            .rational => {
                const tmp = struct {
                    fn isRational(comptime T: type) bool {
                        return T == Rational;
                    }
                };
                return tmp.isRational;
            },
            .real => {
                const tmp = struct {
                    fn isReal(comptime T: type) bool {
                        return T == Real;
                    }
                };
                return tmp.isReal;
            },
            .complex => {
                const tmp = struct {
                    fn isComplex(comptime T: type) bool {
                        if (@hasDecl(T, "is_complex"))
                            return true;

                        return false;
                    }
                };
                return tmp.isComplex;
            },
        };
    }
};

const supported_numeric_types: [33]type = .{
    bool,    u8,
    u16,     u32,
    u64,     u128,
    usize,   c_uint,
    i8,      i16,
    i32,     i64,
    i128,    isize,
    c_int,   comptime_int,
    f16,     f32,
    f64,     f80,
    f128,    comptime_float,
    cf16,    cf32,
    cf64,    cf80,
    cf128,   comptime_cfloat,
    Integer, Rational,
    Real,    Complex(Rational),
    Complex(Real),
    // Expression,
};

const supported_complex_types: [8]type = .{
    cf16,  cf32,
    cf64,  cf80,
    cf128, comptime_cfloat,
    Complex(Rational), Complex(Real), // Expression,
};

pub const VectorType = enum {
    dense,
    sparse,
    numeric, // Fallback for numeric types that are not vectors

    pub fn match(comptime self: VectorType) fn (type) bool {
        return switch (self) {
            .dense => isDenseVector,
            .sparse => isSparseVector,
            .numeric => unreachable,
        };
    }

    pub fn toString(self: VectorType) []const u8 {
        return switch (self) {
            .dense => "vector.Dense",
            .sparse => "vector.Sparse",
            .numeric => unreachable,
        };
    }
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
    numeric, // Fallback for numeric types that are not matrices

    pub fn match(comptime self: MatrixType) fn (type) bool {
        return switch (self) {
            .general_dense => isGeneralDenseMatrix,
            .general_sparse => isGeneralSparseMatrix,
            .symmetric_dense => isSymmetricDenseMatrix,
            .symmetric_sparse => isSymmetricSparseMatrix,
            .hermitian_dense => isHermitianDenseMatrix,
            .hermitian_sparse => isHermitianSparseMatrix,
            .triangular_dense => isTriangularDenseMatrix,
            .triangular_sparse => isTriangularSparseMatrix,
            .diagonal => isDiagonalMatrix,
            .permutation => isPermutationMatrix,
            .numeric => unreachable,
        };
    }

    pub fn toString(self: MatrixType) []const u8 {
        return switch (self) {
            .general_dense => "matrix.general.Dense",
            .general_sparse => "matrix.general.Sparse",
            .symmetric_dense => "matrix.symmetric.Dense",
            .symmetric_sparse => "matrix.symmetric.Sparse",
            .hermitian_dense => "matrix.hermitian.Dense",
            .hermitian_sparse => "matrix.hermitian.Sparse",
            .triangular_dense => "matrix.triangular.Dense",
            .triangular_sparse => "matrix.triangular.Sparse",
            .diagonal => "matrix.Diagonal",
            .permutation => "matrix.Permutation",
            .numeric => unreachable,
        };
    }
};

pub const MatrixKind = enum {
    general,
    symmetric,
    hermitian,
    triangular,
    diagonal,
    permutation,

    pub fn match(comptime self: MatrixKind) fn (type) bool {
        return switch (self) {
            .general => isGeneralMatrix,
            .symmetric => isSymmetricMatrix,
            .hermitian => isHermitianMatrix,
            .triangular => isTriangularMatrix,
            .diagonal => isDiagonalMatrix,
            .permutation => isPermutationMatrix,
        };
    }

    pub fn toString(self: MatrixKind) []const u8 {
        return switch (self) {
            .general => "general matrix",
            .symmetric => "symmetric matrix",
            .hermitian => "hermitian matrix",
            .triangular => "triangular matrix",
            .diagonal => "diagonal matrix",
            .permutation => "permutation matrix",
        };
    }
};

pub const MatrixStorage = enum {
    dense,
    sparse,

    pub fn match(comptime self: MatrixStorage) fn (type) bool {
        return switch (self) {
            .dense => isDenseMatrix,
            .sparse => isSparseMatrix,
        };
    }

    pub fn toString(self: MatrixStorage) []const u8 {
        return switch (self) {
            .dense => "dense matrix",
            .sparse => "sparse matrix",
        };
    }
};

pub const ArrayType = enum {
    dense,
    strided,
    sparse,
    numeric, // Fallback for numeric types that are not arrays

    pub fn match(comptime self: ArrayType) fn (type) bool {
        return switch (self) {
            .dense => isDenseArray,
            .strided => isStridedArray,
            .sparse => isSparseArray,
            .numeric => unreachable,
        };
    }

    pub fn toString(self: ArrayType) []const u8 {
        return switch (self) {
            .dense => "array.Dense",
            .strided => "array.Strided",
            .sparse => "array.Sparse",
            .numeric => unreachable,
        };
    }
};

pub const Domain = enum {
    numeric,
    vector,
    matrix,
    array,
    expression,

    pub fn match(comptime self: Domain) fn (type) bool {
        return switch (self) {
            .numeric => isNumeric,
            .vector => isVector,
            .matrix => isMatrix,
            .array => isArray,
            .expression => isExpression,
        };
    }

    pub fn toString(self: Domain) []const u8 {
        return switch (self) {
            .numeric => "numeric",
            .vector => "vector",
            .matrix => "matrix",
            .array => "array",
            .expression => "expression",
        };
    }
};

/// A useless allocator that does nothing, always signalling allocation failure.
pub const useless_allocator: std.mem.Allocator = .{
    .ptr = undefined,
    .vtable = &vtable,
};

pub const vtable: std.mem.Allocator.VTable = .{
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

/// Checks the the input type `N` and returns the corresponding `NumericType`.
///
/// Checks that the input type is a supported numeric type and returns the
/// corresponding `NumericType` enum value. If the type is not supported, it
/// will raise a compile error.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check.
///
/// Returns
/// -------
/// `NumericType`: The corresponding `NumericType` enum value.
pub inline fn numericType(comptime N: type) NumericType {
    // Without inline functions calling this fail miserably. I have no idea why.

    @setEvalBranchQuota(10000000);

    switch (@typeInfo(N)) {
        .bool => return .bool,
        .int, .comptime_int => return .int,
        .float, .comptime_float => return .float,
        else => {
            if (@hasDecl(N, "is_dyadic"))
                return .dyadic;

            if ((@hasDecl(N, "is_cfloat")) or
                N == std.math.Complex(f16) or N == std.math.Complex(f32) or N == std.math.Complex(f64) or
                N == std.math.Complex(f80) or N == std.math.Complex(f128) or N == std.math.Complex(comptime_float))
                return .cfloat;

            if (N == Integer)
                return .integer;

            if (N == Rational)
                return .rational;

            if (N == Real)
                return .real;

            if (@hasDecl(N, "is_complex"))
                return .complex;

            @compileError("Unsupported numeric type: " ++ @typeName(N));
        },
    }
}

/// Checks if the input type `N` is a supported numeric type, but does not
/// return the corresponding `NumericType` or raise a compile error if not.
///
/// Parameters
/// ----------
/// comptime N (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a supported numeric type, `false` otherwise.
pub fn isNumeric(comptime N: type) bool {
    @setEvalBranchQuota(10000000);

    switch (@typeInfo(N)) {
        .bool => return true,
        .int, .comptime_int => return true,
        .float, .comptime_float => return true,
        else => {
            if (@hasDecl(N, "is_dyadic"))
                return true;

            if ((@hasDecl(N, "is_cfloat")) or
                N == std.math.Complex(f16) or N == std.math.Complex(f32) or N == std.math.Complex(f64) or
                N == std.math.Complex(f80) or N == std.math.Complex(f128) or N == std.math.Complex(comptime_float))
                return true;

            if (N == Integer)
                return true;

            if (N == Rational)
                return true;

            if (N == Real)
                return true;

            if (@hasDecl(N, "is_complex"))
                return true;

            return false;
        },
    }
}

pub inline fn vectorType(comptime V: type) VectorType {
    @setEvalBranchQuota(10000);

    if (isDenseVector(V))
        return .dense;

    if (isSparseVector(V))
        return .sparse;

    return .numeric; // Fallback for numeric types that are not vectors
}

pub inline fn matrixType(comptime M: type) MatrixType {
    @setEvalBranchQuota(10000);

    if (isGeneralDenseMatrix(M))
        return .general_dense;

    if (isSymmetricDenseMatrix(M))
        return .symmetric_dense;

    if (isHermitianDenseMatrix(M))
        return .hermitian_dense;

    if (isTriangularDenseMatrix(M))
        return .triangular_dense;

    if (isGeneralSparseMatrix(M))
        return .general_sparse;

    if (isSymmetricSparseMatrix(M))
        return .symmetric_sparse;

    if (isHermitianSparseMatrix(M))
        return .hermitian_sparse;

    if (isTriangularSparseMatrix(M))
        return .triangular_sparse;

    if (isDiagonalMatrix(M))
        return .diagonal;

    if (isPermutationMatrix(M))
        return .permutation;

    return .numeric; // Fallback for numeric types that are not matrices
}

pub inline fn arrayType(comptime A: type) ArrayType {
    @setEvalBranchQuota(10000);

    if (isDenseArray(A))
        return .dense;

    if (isStridedArray(A))
        return .strided;

    if (isSparseArray(A))
        return .sparse;

    return .numeric; // Fallback for numeric types that are not arrays
}

pub inline fn domain(comptime T: type) Domain {
    @setEvalBranchQuota(10000);

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

    @compileError(@typeName(T) ++ " does not belong to any supported domain");
}

const type_checks = @import("types/type_checks.zig");
pub const isSupportedType = type_checks.isSupportedType;
pub const isPointer = type_checks.isPointer;
pub const isManyPointer = type_checks.isManyPointer;
pub const isConstPointer = type_checks.isConstPointer;
pub const isSlice = type_checks.isSlice;
pub const isSimdVector = type_checks.isSimdVector;
pub const isVector = type_checks.isVector;
pub const isDenseVector = type_checks.isDenseVector;
pub const isSparseVector = type_checks.isSparseVector;
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
pub const isArray = type_checks.isArray;
pub const isDenseArray = type_checks.isDenseArray;
pub const isStridedArray = type_checks.isStridedArray;
pub const isSparseArray = type_checks.isSparseArray;
pub const isExpression = type_checks.isExpression;
pub const isFixedPrecision = type_checks.isFixedPrecision;
pub const isArbitraryPrecision = type_checks.isArbitraryPrecision;
pub const isIntegral = type_checks.isIntegral;
pub const isNonIntegral = type_checks.isNonIntegral;
pub const isReal = type_checks.isReal;
pub const isComplex = type_checks.isComplex;
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

const casting = @import("types/casting.zig");
pub const scast = casting.scast;
pub const cast = casting.cast;

/// Returns the scalar type of a given numeric type, vector, matrix or array.
///
/// This function returns the scalar type of a given numeric type, vector,
/// matrix, or array. If the input type is a vector, a matrix or an array, it
/// returns the element type (equivalent to `Numeric`). If the input type is a
/// numeric type, it returns the type itself, unless it is a complex type, in
/// which case it returns the scalar type of the complex type.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to get the scalar type of. Must be a supported
/// numeric type, vector, matrix, or array.
///
/// Returns
/// -------
/// `type`: The scalar type of the input type.
pub fn Scalar(comptime T: type) type {
    if (isArray(T) or isMatrix(T) or isVector(T)) {
        return Numeric(T);
    }

    switch (comptime numericType(T)) {
        .bool => return T,
        .int => return T,
        .float => return T,
        .cfloat => switch (T) {
            cf16 => return f16,
            cf32 => return f32,
            cf64 => return f64,
            cf80 => return f80,
            cf128 => return f128,
            comptime_cfloat => return comptime_float,
            std.math.Complex(f16) => return f16,
            std.math.Complex(f32) => return f32,
            std.math.Complex(f64) => return f64,
            std.math.Complex(f80) => return f80,
            std.math.Complex(f128) => return f128,
            std.math.Complex(comptime_float) => return comptime_float,
            else => unreachable,
        },
        .integer => return Integer,
        .rational => return Rational,
        .real => return Real,
        .complex => switch (T) {
            Complex(Rational) => return Rational,
            Complex(Real) => return Real,
            else => unreachable,
        },
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
/// Parameters
/// ----------
/// comptime T (`type`): The type to get the numeric type of. Must be a
/// supported numeric type, vector, matrix, or array.
///
/// Returns
/// -------
/// `type`: The underlying numeric type of the input type.
pub fn Numeric(comptime T: type) type {
    if (isExpression(T))
        return Expression;

    if (isArray(T) or isMatrix(T) or isVector(T))
        return T.Numeric;

    if (!isNumeric(T)) {
        @compileError("Expected a numeric type, but got " ++ @typeName(T));
    }

    return T;
}

pub fn orderOf(comptime T: type) Order {
    if (comptime isArray(T)) {
        if (T == array.Dense(Numeric(T), .row_major) or
            T == array.Strided(Numeric(T), .row_major) or
            T == array.Sparse(Numeric(T), .row_major))
            return .row_major
        else
            return .col_major;
    } else if (comptime isMatrix(T)) {
        if (comptime isComplex(Numeric(T))) {
            if (T == matrix.general.Dense(Numeric(T), .row_major) or
                T == matrix.general.Sparse(Numeric(T), .row_major) or
                T == matrix.symmetric.Dense(Numeric(T), .upper, .row_major) or
                T == matrix.symmetric.Dense(Numeric(T), .lower, .row_major) or
                T == matrix.symmetric.Sparse(Numeric(T), .upper, .row_major) or
                T == matrix.symmetric.Sparse(Numeric(T), .lower, .row_major) or
                T == matrix.hermitian.Dense(Numeric(T), .upper, .row_major) or
                T == matrix.hermitian.Dense(Numeric(T), .lower, .row_major) or
                T == matrix.hermitian.Sparse(Numeric(T), .upper, .row_major) or
                T == matrix.hermitian.Sparse(Numeric(T), .lower, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .upper, .non_unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .upper, .unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .non_unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .upper, .non_unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .upper, .unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .non_unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .unit, .row_major))
                return .row_major
            else
                return .col_major;
        } else {
            if (T == matrix.general.Dense(Numeric(T), .row_major) or
                T == matrix.general.Sparse(Numeric(T), .row_major) or
                T == matrix.symmetric.Dense(Numeric(T), .upper, .row_major) or
                T == matrix.symmetric.Dense(Numeric(T), .lower, .row_major) or
                T == matrix.symmetric.Sparse(Numeric(T), .upper, .row_major) or
                T == matrix.symmetric.Sparse(Numeric(T), .lower, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .upper, .non_unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .upper, .unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .non_unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .upper, .non_unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .upper, .unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .non_unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .unit, .row_major))
                return .row_major
            else
                return .col_major;
        }
    } else {
        @compileError("Use `orderOf` only with matrix or array types.");
    }
}

pub fn uploOf(comptime T: type) Uplo {
    if (comptime isMatrix(T)) {
        if (comptime isComplex(Numeric(T))) {
            if (T == matrix.symmetric.Dense(Numeric(T), .lower, orderOf(T)) or
                T == matrix.symmetric.Sparse(Numeric(T), .lower, orderOf(T)) or
                T == matrix.hermitian.Dense(Numeric(T), .lower, orderOf(T)) or
                T == matrix.hermitian.Sparse(Numeric(T), .lower, orderOf(T)) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .unit, orderOf(T)) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .non_unit, orderOf(T)) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .unit, orderOf(T)) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .non_unit, orderOf(T)))
                return .lower
            else
                return .upper;
        } else {
            if (T == matrix.symmetric.Dense(Numeric(T), .lower, orderOf(T)) or
                T == matrix.symmetric.Sparse(Numeric(T), .lower, orderOf(T)) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .unit, orderOf(T)) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .non_unit, orderOf(T)) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .unit, orderOf(T)) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .non_unit, orderOf(T)))
                return .lower
            else
                return .upper;
        }
    } else {
        @compileError("Use `uploOf` only with matrix types.");
    }
}

pub fn diagOf(comptime T: type) Diag {
    if (comptime isMatrix(T)) {
        if (T == matrix.triangular.Dense(Numeric(T), uploOf(T), .unit, orderOf(T)) or
            T == matrix.triangular.Sparse(Numeric(T), uploOf(T), .unit, orderOf(T)))
            return .unit
        else
            return .non_unit;
    } else {
        @compileError("Use `diagOf` only with matrix types.");
    }
}

/// Returns the pointer child type of a given pointer type, or the type itself
/// if the input type is not a pointer.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to get the child type of. Must be a pointer
/// type.
///
/// Returns
/// -------
/// `type`: The child type of the input pointer type.
pub fn Child(comptime T: type) type {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            return info.child;
        },
        .vector => |info| {
            return info.child;
        },
        else => return T,
    }
}

/// Returns the return type of a function with one parameter when called with a
/// value of type `X`.
///
/// This function is useful when the return type of the function depends on
/// the type of the parameter passed to it. The function may have a second
/// parameter that is a context struct.
pub fn ReturnType1(comptime op: anytype, comptime X: type) type {
    const opinfo = @typeInfo(@TypeOf(op));

    comptime if (opinfo.@"fn".params.len != 1 and opinfo.@"fn".params.len != 2)
        @compileError("ReturnType1: op must be a function with one parameter, or two if the second is an context struct");

    const val: X = if (isArbitraryPrecision(X))
        .empty
    else if (isComplex(X))
        .{ .re = 1, .im = 1 }
    else
        1;

    const result_type: type = if (opinfo.@"fn".params.len == 1)
        @TypeOf(op(val))
    else // We must pass an allocator, although it is not used
        if (isArbitraryPrecision(X))
            @TypeOf(op(val, .{ .allocator = useless_allocator }))
        else
            @TypeOf(op(val, .{}));

    const resinfo = @typeInfo(result_type);
    switch (resinfo) {
        .error_union => return resinfo.error_union.payload,
        else => return result_type,
    }
}

/// Returns the return type of a function with two parameters when called with
/// values of types `X` and `Y`.
///
/// This function is useful when the return type of the function depends on
/// the types of the parameters passed to it. The function may have a third
/// parameter that is a context struct.
pub fn ReturnType2(comptime op: anytype, comptime X: type, comptime Y: type) type {
    const opinfo = @typeInfo(@TypeOf(op));

    comptime if (opinfo.@"fn".params.len != 2 and opinfo.@"fn".params.len != 3)
        @compileError("ReturnType2: op must be a function with two parameters, or three if the third is an context struct");

    const val1: X = if (isArbitraryPrecision(X))
        .empty
    else if (isComplex(X))
        .{ .re = 1, .im = 1 }
    else
        1;
    const val2: Y = if (isArbitraryPrecision(Y))
        .empty
    else if (isComplex(Y))
        .{ .re = 1, .im = 1 }
    else
        1;

    const result_type: type = if (opinfo.@"fn".params.len == 2)
        @TypeOf(op(val1, val2))
    else // We must pass an allocator, although it is not used
        if (isArbitraryPrecision(X) or isArbitraryPrecision(Y))
            @TypeOf(op(val1, val2, .{ .allocator = useless_allocator }))
        else
            @TypeOf(op(val1, val2, .{}));

    const resinfo = @typeInfo(result_type);
    switch (resinfo) {
        .error_union => return resinfo.error_union.payload,
        else => return result_type,
    }
}

const context_checks = @import("types/context_checks.zig");
pub const validateContext = context_checks.validateContext;
pub const partialValidateContext = context_checks.partialValidateContext;
pub const ctxHasField = context_checks.ctxHasField;
pub const getFieldOrDefault = context_checks.getFieldOrDefault;
pub const MixStructs = context_checks.MixStructs;
pub const mixStructs = context_checks.mixStructs;
pub const StripStruct = context_checks.StripStruct;
pub const stripStruct = context_checks.stripStruct;
pub const RenameStructFields = context_checks.RenameStructFields;
pub const renameStructFields = context_checks.renameStructFields;
pub const KeepStructFields = context_checks.KeepStructFields;
pub const keepStructFields = context_checks.keepStructFields;
