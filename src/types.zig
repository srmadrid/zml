//! Namespace for type definitions and utilities.

const std = @import("std");

const constants = @import("constants.zig");

pub const default_uint = u32;
pub const default_int = i32;
pub const default_float = f64;

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
};

pub const MatrixType = enum {
    dense_general,
    dense_symmetric,
    dense_hermitian,
    dense_triangular,
    banded_general,
    banded_symmetric,
    banded_hermitian,
    tridiagonal_general,
    tridiagonal_symmetric,
    tridiagonal_hermitian,
    sparse_general,
    sparse_symmetric,
    sparse_hermitian,
    sparse_triangular,
    block_general,
    block_symmetric,
    block_hermitian,
    diagonal,
    permutation,
    numeric, // Fallback for numeric types that are not matrices

    pub fn match(comptime self: MatrixType) fn (type) bool {
        return switch (self) {
            .dense_general => isGeneralDenseMatrix,
            .dense_symmetric => isSymmetricDenseMatrix,
            .dense_hermitian => isHermitianDenseMatrix,
            .dense_triangular => isTriangularDenseMatrix,
            .banded_general => isGeneralBandedMatrix,
            .banded_symmetric => isSymmetricBandedMatrix,
            .banded_hermitian => isHermitianBandedMatrix,
            .tridiagonal_general => isGeneralTridiagonalMatrix,
            .tridiagonal_symmetric => isSymmetricTridiagonalMatrix,
            .tridiagonal_hermitian => isHermitianTridiagonalMatrix,
            .sparse_general => isGeneralSparseMatrix,
            .sparse_symmetric => isSymmetricSparseMatrix,
            .sparse_hermitian => isHermitianSparseMatrix,
            .sparse_triangular => isTriangularSparseMatrix,
            .block_general => isGeneralBlockMatrix,
            .block_symmetric => isSymmetricBlockMatrix,
            .block_hermitian => isHermitianBlockMatrix,
            .diagonal => isDiagonalMatrix,
            .permutation => isPermutationMatrix,
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
};

pub const MatrixStorage = enum {
    dense,
    banded,
    tridiagonal,
    sparse,
    block,

    pub fn match(comptime self: MatrixStorage) fn (type) bool {
        return switch (self) {
            .dense => isDenseMatrix,
            .banded => isBandedMatrix,
            .tridiagonal => isTridiagonalMatrix,
            .sparse => isSparseMatrix,
            .block => isBlockMatrix,
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
        return .dense_general;

    if (isSymmetricDenseMatrix(M))
        return .dense_symmetric;

    if (isHermitianDenseMatrix(M))
        return .dense_hermitian;

    if (isTriangularDenseMatrix(M))
        return .dense_triangular;

    if (isGeneralBandedMatrix(M))
        return .banded_general;

    if (isSymmetricBandedMatrix(M))
        return .banded_symmetric;

    if (isHermitianBandedMatrix(M))
        return .banded_hermitian;

    if (isGeneralTridiagonalMatrix(M))
        return .tridiagonal_general;

    if (isSymmetricTridiagonalMatrix(M))
        return .tridiagonal_symmetric;

    if (isHermitianTridiagonalMatrix(M))
        return .tridiagonal_hermitian;

    if (isGeneralSparseMatrix(M))
        return .sparse_general;

    if (isSymmetricSparseMatrix(M))
        return .sparse_symmetric;

    if (isHermitianSparseMatrix(M))
        return .sparse_hermitian;

    if (isTriangularSparseMatrix(M))
        return .sparse_triangular;

    if (isGeneralBlockMatrix(M))
        return .block_general;

    if (isSymmetricBlockMatrix(M))
        return .block_symmetric;

    if (isHermitianBlockMatrix(M))
        return .block_hermitian;

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

pub inline fn domainType(comptime T: type) Domain {
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

    @compileError("Unsupported type for domainType: " ++ @typeName(T));
}

/// Checks if the input type is a supported type (numeric, vector, matrix,
/// array, or expression).
pub fn isSupportedType(comptime T: type) bool {
    if (comptime isNumeric(T))
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
    switch (numericType(T)) {
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
    switch (numericType(T)) {
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
    switch (numericType(T)) {
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
    switch (numericType(T)) {
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
    switch (numericType(T)) {
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
    switch (numericType(T)) {
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
    switch (numericType(T)) {
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
    switch (numericType(T)) {
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

const CheckedDomain = struct {
    first_level: FirstLevel, // The first level type (domain, numeric, vector, matrix_kind, matrix_storage, matrix, array)
    element_type: ?[]const u8, // The element type expression, if any (e.g., "float" in "matrix.dense(float)")

    const FirstLevel = union(enum) {
        domain: Domain,
        numeric: []const u8, // The numeric type expression as a string since we allow more complex expressions here
        vector: VectorType,
        matrix_kind: MatrixKind,
        matrix_storage: MatrixStorage,
        matrix: MatrixType,
        array: ArrayType,

        pub fn fromAny(comptime s: anytype) FirstLevel {
            return switch (@TypeOf(s)) {
                Domain => .{ .domain = s },
                NumericType => .{ .numeric = s },
                VectorType => .{ .vector = s },
                MatrixKind => .{ .matrix_kind = s },
                MatrixStorage => .{ .matrix_storage = s },
                MatrixType => .{ .matrix = s },
                ArrayType => .{ .array = s },
                else => unreachable,
            };
        }
    };
};

fn mixSubdomains(comptime base: Domain, kind: anytype, storage: anytype) switch (base) {
    .matrix => MatrixType,
    else => unreachable,
} {
    switch (base) {
        .matrix => {
            const K = @TypeOf(kind);
            const S = @TypeOf(storage);

            if (K != MatrixKind or S != MatrixStorage) {
                @compileError("For base domain 'matrix', kind and storage must be of types MatrixKind and MatrixStorage, respectively.");
            }

            switch (kind) {
                .general => switch (storage) {
                    .dense => return .dense_general,
                    .banded => return .banded_general,
                    .tridiagonal => return .tridiagonal_general,
                    .sparse => return .sparse_general,
                    .block => return .block_general,
                },
                .symmetric => switch (storage) {
                    .dense => return .dense_symmetric,
                    .banded => return .banded_symmetric,
                    .tridiagonal => return .tridiagonal_symmetric,
                    .sparse => return .sparse_symmetric,
                    .block => return .block_symmetric,
                },
                .hermitian => switch (storage) {
                    .dense => return .dense_hermitian,
                    .banded => return .banded_hermitian,
                    .tridiagonal => return .tridiagonal_hermitian,
                    .sparse => return .sparse_hermitian,
                    .block => return .block_hermitian,
                },
                .triangular => switch (storage) {
                    .dense => return .dense_triangular,
                    .banded => unreachable,
                    .tridiagonal => unreachable,
                    .sparse => return .sparse_triangular,
                    .block => unreachable,
                },
                .diagonal => return .diagonal,
                .permutation => return .permutation,
            }
        },
        else => unreachable,
    }
}

const CategoryOperator = enum {
    /// bool, int, integer
    integral,
    /// float, dyadic, cfloat, rational, real, complex
    nonintegral,
    /// bool, int, float, dyadic, cfloat
    fixed,
    /// integer, rational, real, complex
    arbitrary,
    /// bool, int, float, dyadic, integer, rational, real
    real,
    /// cfloat, complex
    complex,
    /// signed int, float, dyadic, cfloat, integer, rational, real, complex
    signed,
    /// bool, unsigned int
    unsigned,

    pub fn match(self: CategoryOperator) fn (type) bool {
        return switch (self) {
            .integral => return isIntegral,
            .nonintegral => return isNonIntegral,
            .fixed => return isFixedPrecision,
            .arbitrary => return isArbitraryPrecision,
            .real => return isReal,
            .complex => return isComplex,
            .signed => return isSigned,
            .unsigned => return isUnsigned,
        };
    }
};

/// Given a type signature (as a string) and a type, checks whether the type
/// matches the signature, triggering a compile-time error if not.
///
/// Parameters
/// ----------
/// comptime type_signature (`[]const u8`): The type signature to check against.
/// The type signature must start with:
///
/// - "": pass by value
/// - "*": mutable one-item pointer
/// - "* const": constant one-item pointer
/// - "[*]": mutable many-item pointer
/// - "[*] const": constant many-item pointer
/// - "[]": mutable slice
/// - "[] const": constant slice
///
/// Then, it must specify the type:
///
/// - "numeric": any supported numeric type. For specific numeric types, use the
///   names defined in the `NumericType` enum.
/// - "vector": any vector type. For specific vector types, use "vector.<storage>",
///   where `<storage>` is one of:
///
///     - "dense": dense vector
///     - "sparse": sparse vector
///
/// - "matrix": any matrix type. For specific matrix types, use "matrix.<kind>.<storage>",
///   where `<kind>` is one of:
///
///     - "general": general matrix
///     - "symmetric": symmetric matrix
///     - "hermitian": hermitian matrix
///     - "triangular": triangular matrix
///     - "diagonal": diagonal matrix (<storage> is ignored)
///     - "permutation": permutation matrix (<storage> is ignored)
///
///   and `<storage>` is one of:
///
///     - "dense": dense matrix
///     - "banded": banded matrix
///     - "tridiagonal": tridiagonal matrix
///     - "sparse": sparse matrix
///     - "block": block sparse matrix
///   Setting only the kind (e.g., "matrix.symmetric") or only the storage
///   (e.g., "matrix.sparse") is also accepted, allowing for more general checks,
///   otherwise the order must be `<kind>.<storage>`.
/// - "array": any array type. For specific array types, use "array.<storage>",
///   where `<storage>` is one of:
///     - "dense": dense array
///     - "strided": strided array
///     - "sparse": sparse array
/// - "expression": just use `Expression` as the type, why are you even checking this?
/// - "any": any type is accepted
///
/// In any container type (vector, matrix, array, etc.), after the type name and
/// without spaces, you can add a numeric type expression in angle brackets to
/// specify the element type. For example, "matrix.dense(float)" specifies a
/// dense matrix with float elements. If no element type is specified, any
/// element type is accepted. Inside the angle brackets, the same rules as for
/// numeric types apply, including the more granular constraints described below.
/// For more granular control over numeric types, constraint operators can be
/// used. Any operator <op> must be placed immediately before the type it
/// applies to, and with no spaces, like <op><type>. We define three kinds of
/// operators, principal operators (which can only be used once or twice at the
/// beginning of the numeric type expression), additional operators (which can
/// be used multiple times and in any order after the principal operators), and
/// final operators (which must be at the end of the expression and apply to
/// everything before them). The principal operators are:
///
/// - Ordering operators (only for numeric types, following the order in the
///   `NumericType` enum):
///     - "<": only for numeric types, accepts types lower than the specified type.
///     - "<=": only for numeric types, accepts types lower than or equal to the
///     specified type.
///     - ">": only for numeric types, accepts types higher than the specified type.
///     - ">=": only for numeric types, accepts types higher than or equal to the
///     specified type.
/// - "=": accepts only the specified type. It is equivalent to not using any operator
///
/// but may speed up compilation in some cases.
/// If two ordering operators are used, the first must be a lower bound ("<" or "<=")
/// and the second an upper bound (">" or ">="). For example, "<=dyadic >=int" accepts
/// any numeric type between int and float, inclusive (in this case, ints, floats and
/// dyadics). The additional operators can also be used for different domains. However,
/// if used in this way and a numeric type is specified as a base type (not inside a container),
/// the numeric part must be at the end of the signature. For example,
/// "matrix.sparse(@real) <op> matrix.block(@fixed @real) <op> >=int <=dyadic" is valid,
/// but ">=int <=dyadic <op> matrix.sparse(@real) <op> matrix.block(@fixed @real)" is not.
/// The additional operators are:
///
/// - "!": exclusion operator, for excluding types. For example, "numeric !float !complex"
/// accepts any numeric type except float and complex types. The principal part
/// must be "numeric" or a range when using exclusion operators on numeric types, or
/// "any" or a general container type (vector, matrix, array) when excluding container types,
/// for example, "matrix !matrix.sparse" accepts any matrix type except sparse matrices, or
/// "any !numeric" accepts any type except numeric types.
/// - "|": inclusion operator, for including only specific types. The principal
/// part must be a specific type or a range when using inclusion operators.
/// For instance, ">=integer |int" accepts anything from integer upwards, but
/// also int (which is below integer). For container types this works similarly,
/// e.g., "matrix |vector" accepts any matrix or vector type.
///
/// The final operators are all category operators, and are prefixed with "@":
///
/// - "@integral": accepts only integral types (bool, int, and integer).
/// - "@nonintegral": accepts only non-integral types (rational, real, float, dyadic, cfloat, complex, and expression).
/// - "@fixed": accepts only fixed precision types (bool, int, float, dyadic, and cfloat).
/// - "@arbitrary": accepts only arbitrary precision types (integer, rational,
/// real, and complex).
/// - "@real": accepts only real types (bool, int, integer, rational, real, float).
/// - "@complex": accepts only complex types (cfloat, complex).
/// - "@signed": accepts only signed types (int, integer, rational, real, float, cfloat,
/// complex).
/// - "@unsigned": accepts only unsigned types (bool and uint).
///
/// Final operators can be combined. For example, "numeric @fixed @real" accepts
/// any fixed precision real type (bool, int, float, dyadic). They also do not need
/// the principal part, so "@real" alone is valid and accepts any real numeric type.
/// However, ain important constraint is that final operators apply to everything else
/// within the same numeric type expression. For example, "@real !float |complex"
/// still excludes complex types, because the category operator is applied last.
/// The type signature is case-insensitive and ignores repeated spaces. Below are
/// examples of valid type signatures. We divide by the nature of the examples:
///
/// - Basic signatures (pointer prefixes):
///     - "numeric": any numeric type
///     - "* numeric": mutable one-item pointer to any numeric type
///     - "* const numeric": constant one-item pointer to any numeric type
///     - "[] numeric": mutable slice of any numeric type
/// - Basic numeric expressions:
///     - "=float": only float, equivalent to just "float"
///     - ">=int": int or any type higher than int
///     - ">=int <=dyadic": any type between int and dyadic, inclusive
///     - "numeric !float !complex": any numeric type except float and complex types
///     - ">=integer |int": anything from integer upwards, but also int
///     - "numeric @real": any real numeric type (not cfloat or complex)
///     - ">=int @signed": any signed type equal to or higher than int
///     - "numeric @fixed @real": any fixed precision real type (bool, int, float, dyadic)
/// - Basic container types:
///     - "vector": any vector type holding any element type
///     - "matrix.symmetric.": any symmetric matrix type holding any element type
///     - "matrix.general.dense": a `matrix.general.Dense` holding any element type
/// - Container types with numeric expressions:
///     - "vector(float)": any vector type holding float elements
///     - "matrix(@complex)": any matrix type holding complex element types
/// - Complex combinations:
///     - "[] matrix.sparse(@real)": mutable slice of sparse matrices holding real element types
///     - "* const matrix.symmetric(rational |integer)": constant one-item pointer to symmetric matrices
///     holding rational or integer element types
///     - "matrix.general.dense<numeric !float @signed>": a `matrix.general.Dense`
///     holding any signed numeric element type except float
///     - "matrix |vector": any matrix or vector type holding any element type
///
/// comptime T (`type`): The type to check against the signature.
///
/// comptime fn_name (`[]const u8`): The name of the function from which this
/// check is being called. Used for error messages.
///
/// comptime param_name (`[]const u8`): The name of the parameter being checked.
/// Used for error messages.
///
/// Returns
/// -------
/// `void`
pub fn checkParameterType(
    comptime type_signature: []const u8,
    comptime T: type,
    comptime fn_name: []const u8,
    comptime param_name: []const u8,
) void {
    comptime var VT: type = T;

    // Preprocess type signature: lowercase and trim leading spaces
    comptime var signature_array: [type_signature.len + 1]u8 = undefined;
    for (type_signature, 0..) |c, i|
        signature_array[i] = std.ascii.toLower(c);
    signature_array[type_signature.len] = 0;
    comptime var signature: []const u8 = signature_array[0..signature_array.len];

    // Pointer or slice prefix matching
    signature = std.mem.trimStart(u8, signature, " ");

    comptime var pointer: ?std.builtin.Type.Pointer.Size = null;
    comptime var constant: bool = false;
    switch (signature[0]) {
        '*' => {
            pointer = .one;

            if (signature.len >= 7 and
                std.mem.eql(u8, signature[1..7], " const"))
            {
                // constant one-item pointer
                if (!isPointer(VT)) // any non-const pointer can be coerced to const
                    @compileError(
                        std.fmt.comptimePrint(
                            "{s} requires parameter '{s}' to be a one-item pointer, but got '{s}'.",
                            .{ fn_name, param_name, @typeName(T) },
                        ),
                    );

                constant = true;
                VT = Child(VT);
                signature = signature[7..];
            } else {
                // mutable one-item pointer
                if (!isPointer(VT) or isConstPointer(VT))
                    @compileError(
                        std.fmt.comptimePrint(
                            "{s} requires parameter '{s}' to be a mutable one-item pointer, but got '{s}'.",
                            .{ fn_name, param_name, @typeName(T) },
                        ),
                    );

                VT = Child(VT);
                signature = signature[1..];
            }
        },
        '[' => {
            switch (signature[1]) {
                '*' => {
                    pointer = .many;

                    if (signature.len >= 8 and
                        std.mem.eql(u8, signature[2..8], "] const"))
                    {
                        // constant many-item pointer
                        if (!isManyPointer(VT)) // any non-const pointer can be coerced to const
                            @compileError(
                                std.fmt.comptimePrint(
                                    "{s} requires parameter '{s}' to be a many-item pointer, but got '{s}'.",
                                    .{ fn_name, param_name, @typeName(T) },
                                ),
                            );

                        VT = Child(VT);
                        signature = signature[8..];
                    } else {
                        // mutable many-item pointer
                        if (!isManyPointer(VT) or isConstPointer(VT))
                            @compileError(
                                std.fmt.comptimePrint(
                                    "{s} requires parameter '{s}' to be a mutable many-item pointer, but got '{s}'.",
                                    .{ fn_name, param_name, @typeName(T) },
                                ),
                            );

                        VT = Child(VT);
                        signature = signature[2..];
                    }
                },
                ']' => {
                    pointer = .slice;

                    if (signature.len >= 7 and
                        std.mem.eql(u8, signature[2..7], " const"))
                    {
                        // constant slice
                        if (!isSlice(VT))
                            @compileError(
                                std.fmt.comptimePrint(
                                    "{s} requires parameter '{s}' to be a slice, but got '{s}'.",
                                    .{ fn_name, param_name, @typeName(T) },
                                ),
                            );

                        VT = Child(VT);
                        signature = signature[7..];
                    } else {
                        // mutable slice
                        if (!isSlice(VT) or isConstPointer(VT))
                            @compileError(
                                std.fmt.comptimePrint(
                                    "{s} requires parameter '{s}' to be a mutable slice, but got '{s}'.",
                                    .{ fn_name, param_name, @typeName(T) },
                                ),
                            );

                        VT = Child(VT);
                        signature = signature[2..];
                    }
                },
                else => @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized pointer type in type signature '{s}' in parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                ),
            }
        },
        else => {},
    }

    signature = std.mem.trimStart(u8, signature, " ");

    // Actual type checking
    comptime var matched_any_positive = false; // Must match at least one positive type ("" or "|")
    comptime var matched_any_negative = false; // Must not match any negative type ("!")
    comptime var finished = false;
    comptime var checked = 0; // Number of principal checks done
    comptime var domains_checked: [32]?CheckedDomain = .{null} ** 32; // For compilation error messages
    while (!finished) : (checked += 1) {
        signature = std.mem.trimStart(u8, signature, " ");

        if (signature.len == 0 or signature[0] == 0) {
            finished = true;
            break;
        }

        // Check for "|" and "!" operators
        comptime var positive_operator = true;
        if (signature.len >= 1 and signature[0] == '|') {
            if (checked == 0) {
                @compileError(
                    std.fmt.comptimePrint(
                        "The '|' operator cannot be used as the first operator in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );
            }

            positive_operator = true;
            signature = signature[1..];
        } else if (signature.len >= 1 and signature[0] == '!') {
            if (checked == 0) {
                @compileError(
                    std.fmt.comptimePrint(
                        "The '!' operator cannot be used as the first operator in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );
            }

            positive_operator = false;
            signature = signature[1..];
        }

        // Check for "any". Only valid if nothing else has been checked yet
        if (std.mem.eql(u8, signature[0..3], "any")) {
            if (checked == 0) {
                // "any" matches everything supported by the library, ignore rest of signature
                if (!isSupportedType(VT))
                    @compileError(
                        std.fmt.comptimePrint(
                            "{s} requires parameter '{s}' to be any supported type, but got '{s}'.",
                            .{ fn_name, param_name, @typeName(T) },
                        ),
                    )
                else
                    return; // Matched anything, ignore rest of signature
            } else {
                @compileError(
                    std.fmt.comptimePrint(
                        "The 'any' type can only be used as the first principal type, in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );
            }
        }

        comptime var matched_domain: Domain = .numeric;

        // Check for domain types
        inline for (std.meta.fields(Domain)) |domain| {
            if (domain.name.len <= signature.len and
                std.mem.eql(u8, signature[0..domain.name.len], domain.name))
            {
                matched_domain = @enumFromInt(domain.value);

                if (matched_domain == .numeric) {
                    // Numeric types are handled later
                    continue;
                }

                // Consume matched part
                domains_checked[checked] = .{ .first_level = .{ .domain = matched_domain }, .element_type = null };
                signature = signature[domain.name.len..];

                break;
            }
        }

        if (matched_domain == .numeric) {
            // If numeric domain, either we matched numeric or nothing yet
            const numeric_constraint = checkNumericConstraints(
                type_signature,
                signature,
                T,
                fn_name,
                param_name,
                &domains_checked[checked],
            );

            if (positive_operator)
                matched_any_positive = matched_any_positive or numeric_constraint
            else
                matched_any_negative = matched_any_negative or numeric_constraint;

            // Numeric domain must be the only one or the last one
            finished = true;
        } else {
            // Need to validate domain match

            // Domain refinement
            switch (matched_domain) {
                .vector => {
                    // Vector type refinement, can only have one level
                    const one_level_check = oneLevelCheck(
                        type_signature,
                        &signature,
                        VT,
                        fn_name,
                        param_name,
                        .vector,
                        VectorType,
                        &domains_checked[checked],
                    );

                    if (positive_operator)
                        matched_any_positive = matched_any_positive or one_level_check
                    else
                        matched_any_negative = matched_any_negative or one_level_check;
                },
                .matrix => {
                    // Matrix type refinement
                    const two_level_check = twoLevelCheck(
                        type_signature,
                        &signature,
                        VT,
                        fn_name,
                        param_name,
                        .matrix,
                        MatrixKind,
                        MatrixStorage,
                        &domains_checked[checked],
                    );

                    if (positive_operator)
                        matched_any_positive = matched_any_positive or two_level_check
                    else
                        matched_any_negative = matched_any_negative or two_level_check;
                },
                .array => {
                    // Array type refinement, can only have one level
                    const one_level_check = oneLevelCheck(
                        type_signature,
                        &signature,
                        VT,
                        fn_name,
                        param_name,
                        .array,
                        ArrayType,
                        &domains_checked[checked],
                    );

                    if (positive_operator)
                        matched_any_positive = matched_any_positive or one_level_check
                    else
                        matched_any_negative = matched_any_negative or one_level_check;
                },
                else => unreachable,
            }
        }
    }

    if (!matched_any_positive or matched_any_negative) {
        // Build error message
        // Domains with | are separated by " or "
        // Domains with ! are like any except ...
        // Numerics are put with the string expression directly

        // For now, just print the type signature
        const error_message = std.fmt.comptimePrint(
            "{s} requires parameter '{s}' to be of type '{s}', but got '{s}'.",
            .{ fn_name, param_name, type_signature, @typeName(T) },
        );

        @compileError(error_message);
    }
}

fn oneLevelCheck(
    comptime type_signature: []const u8,
    comptime signature: *[]const u8,
    comptime T: type,
    comptime fn_name: []const u8,
    comptime param_name: []const u8,
    comptime domain: Domain,
    comptime DomainType: type,
    comptime checked_domain: *?CheckedDomain,
) bool {
    comptime var type_base_name: [32:0]u8 = .{0} ** 32;
    @memcpy(type_base_name[0..domain.toString().len], domain.toString());
    comptime var filled_length = domain.toString().len;

    comptime var matched_subtype: ?DomainType = null;
    if (signature.*[0] == '.') {
        signature.* = signature.*[1..]; // Consume "."
        inline for (std.meta.fields(DomainType)) |subtype| {
            if (subtype.name.len <= signature.len and
                std.mem.eql(u8, signature.*[0..subtype.name.len], subtype.name))
            {
                // Consume matched part
                signature.* = signature.*[subtype.name.len..];
                matched_subtype = @enumFromInt(subtype.value);

                break;
            }
        }

        if (matched_subtype == null)
            @compileError(
                std.fmt.comptimePrint(
                    "Unrecognized {s} type in type signature '{s}' for parameter '{s}' of function '{s}'.",
                    .{ std.mem.span(@as([*:0]const u8, @ptrCast(&type_base_name))), type_signature, param_name, fn_name },
                ),
            );

        type_base_name[filled_length] = '.';
        const tag_name = @tagName(matched_subtype.?);
        type_base_name[filled_length + 1] = std.ascii.toUpper(tag_name[0]);
        filled_length += 2;
        @memcpy(type_base_name[filled_length .. filled_length + tag_name.len - 1], tag_name[1..]);
        filled_length += tag_name.len - 1;
    }

    // Check for element type specification
    comptime var element_type_specified: bool = false;
    comptime var correct_element_type: bool = true; // Default to true for no specification
    comptime var element_type_signature: []const u8 = "";
    if (signature.len > 0 and signature.*[0] == '(') {
        const closing_paren_index = std.mem.indexOfScalar(u8, signature.*, ')');
        if (closing_paren_index == null) {
            @compileError(
                std.fmt.comptimePrint(
                    "Unmatched '(' in type signature '{s}' for parameter '{s}' of function '{s}'.",
                    .{ type_signature, param_name, fn_name },
                ),
            );
        }

        element_type_specified = true;
        element_type_signature = signature.*[1..closing_paren_index.?];

        // Check element type
        correct_element_type = checkNumericConstraints(
            type_signature,
            element_type_signature,
            Numeric(T),
            fn_name,
            param_name,
            checked_domain,
        );

        // Consume element type part, including parentheses
        signature.* = signature.*[closing_paren_index.? + 1 ..];

        checked_domain.*.?.element_type = element_type_signature;
    }

    // Validate type match
    comptime var correct_subdomain_type: bool = false;
    if (matched_subtype) |subtype| {
        checked_domain.*.?.first_level = CheckedDomain.FirstLevel.fromAny(subtype);
        correct_subdomain_type = subtype.match()(T);
    } else {
        checked_domain.*.?.first_level = .{ .domain = domain };
        correct_subdomain_type = domain.match()(T);
    }

    return correct_subdomain_type and correct_element_type;
}

fn twoLevelCheck(
    comptime type_signature: []const u8,
    comptime signature: *[]const u8,
    comptime T: type,
    comptime fn_name: []const u8,
    comptime param_name: []const u8,
    comptime domain: Domain,
    comptime DomainType1: type, // Can be skipped if only one level is present
    comptime DomainType2: type,
    comptime checked_domain: *?CheckedDomain,
) bool {
    comptime var type_base_name: [32:0]u8 = .{0} ** 32;
    @memcpy(type_base_name[0..domain.toString().len], domain.toString());
    comptime var filled_length = domain.toString().len;

    comptime var matched_subtype1: ?DomainType1 = null;
    comptime var matched_subtype2: ?DomainType2 = null;

    // First subtype
    if (signature.*[0] == '.') {
        signature.* = signature.*[1..]; // Consume "."

        // Test first DomainType1
        inline for (std.meta.fields(DomainType1)) |subtype1| {
            if (subtype1.name.len <= signature.len and
                std.mem.eql(u8, signature.*[0..subtype1.name.len], subtype1.name))
            {
                // Consume matched part
                signature.* = signature.*[subtype1.name.len..];
                matched_subtype1 = @enumFromInt(subtype1.value);

                break;
            }
        }

        if (matched_subtype1 == null) {
            // Test second DomainType2 if first failed
            inline for (std.meta.fields(DomainType2)) |subtype2| {
                if (subtype2.name.len <= signature.len and
                    std.mem.eql(u8, signature.*[0..subtype2.name.len], subtype2.name))
                {
                    // Consume matched part
                    signature.* = signature.*[subtype2.name.len..];
                    matched_subtype2 = @enumFromInt(subtype2.value);

                    break;
                }
            }

            if (matched_subtype2 == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized {s} type in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ std.mem.span(@as([*:0]const u8, @ptrCast(&type_base_name))), type_signature, param_name, fn_name },
                    ),
                );

            type_base_name[filled_length] = '.';
            const tag_name2 = @tagName(matched_subtype2.?);
            filled_length += 1;
            @memcpy(type_base_name[filled_length .. filled_length + tag_name2.len], tag_name2);
            filled_length += tag_name2.len;
        } else {
            type_base_name[filled_length] = '.';
            const tag_name1 = @tagName(matched_subtype1.?);
            filled_length += 1;
            @memcpy(type_base_name[filled_length .. filled_length + tag_name1.len], tag_name1);
            filled_length += tag_name1.len;
        }
    }

    // Second subtype
    if (matched_subtype2 == null) {
        if (signature.len > 0 and signature.*[0] == '.') {
            signature.* = signature.*[1..]; // Consume "."

            // Second is always DomainType2
            inline for (std.meta.fields(DomainType2)) |subtype2| {
                if (subtype2.name.len <= signature.len and
                    std.mem.eql(u8, signature.*[0..subtype2.name.len], subtype2.name))
                {
                    // Consume matched part
                    signature.* = signature.*[subtype2.name.len..];
                    matched_subtype2 = @enumFromInt(subtype2.value);

                    break;
                }
            }

            if (matched_subtype2 == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized {s} type in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ std.mem.span(@as([*:0]const u8, @ptrCast(&type_base_name))), type_signature, param_name, fn_name },
                    ),
                );

            type_base_name[filled_length] = '.';
            const tag_name2 = @tagName(matched_subtype2.?);
            type_base_name[filled_length + 1] = std.ascii.toUpper(tag_name2[0]);
            filled_length += 2;
            @memcpy(type_base_name[filled_length .. filled_length + tag_name2.len - 1], tag_name2[1..]);
            filled_length += tag_name2.len - 1;
        }
    } else if (signature.len > 0 and signature.*[0] == '.') {
        @compileError(
            std.fmt.comptimePrint(
                "Found unexpected second {s} type in type signature '{s}' for parameter '{s}' of function '{s}'.",
                .{ std.mem.span(@as([*:0]const u8, @ptrCast(&type_base_name))), type_signature, param_name, fn_name },
            ),
        );
    }

    // Check for element type specification
    comptime var element_type_specified: bool = false;
    comptime var correct_element_type: bool = true; // Default to true for no specification
    comptime var element_type_signature: []const u8 = "";
    if (signature.len > 0 and signature.*[0] == '(') {
        const closing_paren_index = std.mem.indexOfScalar(u8, signature.*, ')');
        if (closing_paren_index == null) {
            @compileError(
                std.fmt.comptimePrint(
                    "Unmatched '(' in type signature '{s}' for parameter '{s}' of function '{s}'.",
                    .{ type_signature, param_name, fn_name },
                ),
            );
        }

        element_type_specified = true;
        element_type_signature = signature.*[1..closing_paren_index.?];

        // Check element type
        correct_element_type = checkNumericConstraints(
            type_signature,
            element_type_signature,
            Numeric(T),
            fn_name,
            param_name,
            checked_domain,
        );

        // Consume element type part, including parentheses
        signature.* = signature.*[closing_paren_index.? + 1 ..];
    }

    // Validate type match
    comptime var correct_subdomain_type: bool = false;
    if (matched_subtype1 != null and matched_subtype2 != null) {
        // Both active
        checked_domain.*.?.first_level = CheckedDomain.FirstLevel.fromAny(mixSubdomains(domain, matched_subtype1.?, matched_subtype2.?));
        correct_subdomain_type = matched_subtype1.?.match()(T) and matched_subtype2.?.match()(T);
    } else if (matched_subtype1) |subtype1| {
        // Only matched_subtype1 active
        checked_domain.*.?.first_level = CheckedDomain.FirstLevel.fromAny(subtype1);
        correct_subdomain_type = subtype1.match()(T);
    } else if (matched_subtype2) |subtype2| {
        // Only matched_subtype2 active
        checked_domain.*.?.first_level = CheckedDomain.FirstLevel.fromAny(subtype2);
        correct_subdomain_type = subtype2.match()(T);
    } else {
        checked_domain.*.?.first_level = .{ .domain = domain };
        correct_subdomain_type = domain.match()(T);
    }

    return correct_subdomain_type and correct_element_type;
}

/// Given a numeric type signature (as a string) and a numeric type, checks
/// whether the type matches the signature, returning `true` if it does,
/// `false` otherwise. Does not trigger compile-time errors.
///
/// Parameters
/// ----------
/// comptime signature (`[]const u8`): The numeric type signature to check against.
/// See the documentation of `checkParameterType` for the format of type signatures.
///
/// comptime T (`type`): The numeric type to check against the signature.
/// Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type matches the signature, `false` otherwise.
fn checkNumericConstraints(
    comptime type_signature: []const u8,
    comptime numeric_signature: []const u8,
    comptime T: type,
    comptime fn_name: []const u8,
    comptime param_name: []const u8,
    comptime checked_domain: *?CheckedDomain, // Is null if called from numeric domain, else only element_type is null
) bool {
    // First check T is numeric
    const is_numeric = isNumeric(T);

    _ = checked_domain; // Currently unused
    comptime var signature: []const u8 = std.mem.trimStart(u8, numeric_signature, " ");
    const actual_type: ?NumericType = if (is_numeric) numericType(T) else null;

    // Remove trailing 0 (implicit null terminator in string literals)
    signature = std.mem.trimEnd(u8, signature, &.{0});

    // If signature is only one word "numeric", accept any numeric type
    if (7 <= signature.len and
        std.mem.eql(u8, signature, "numeric"))
        return true and is_numeric;

    // Test principal operators
    // "=" means the same as no operator, but can only have one numeric type category, or all of them ("numeric")
    if (std.mem.startsWith(u8, signature, "=")) {
        signature = signature[1..];

        if (7 == signature.len and
            std.mem.eql(u8, signature, "numeric"))
            return true and is_numeric;

        // Exact match
        comptime var matched_type: ?NumericType = null;
        inline for (std.meta.fields(NumericType)) |num_type| {
            if (num_type.name.len <= signature.len and
                std.mem.eql(u8, signature[0..num_type.name.len], num_type.name))
            {
                // Consume matched part
                signature = signature[num_type.name.len..];
                matched_type = @enumFromInt(num_type.value);

                break;
            }
        }

        if (matched_type == null)
            return false and is_numeric;

        return is_numeric and matched_type.? == actual_type.?;
    }

    // Other principal operators: "<", "<=", ">", ">="
    comptime var lower_bound: ?NumericType = null;
    comptime var lower_inclusive: bool = false;
    comptime var upper_bound: ?NumericType = null;
    comptime var upper_inclusive: bool = false;
    comptime var finish_principal = false;
    while (!finish_principal) {
        signature = std.mem.trimStart(u8, signature, " ");

        if (std.mem.indexOfScalar(u8, signature, '<') == null and
            std.mem.indexOfScalar(u8, signature, '>') == null)
        {
            // Must be at the beginning, so no more principal operators
            finish_principal = true;
            break;
        }

        // If not <= or <, must be >= or >
        const is_lower = std.mem.startsWith(u8, signature, ">");
        signature = signature[1..]; // Consume "<" or ">"
        const is_inclusive = std.mem.startsWith(u8, signature, "=");
        signature = if (is_inclusive) signature[1..] else signature;

        // Match numeric type
        comptime var matched_type: ?NumericType = null;
        inline for (std.meta.fields(NumericType)) |num_type| {
            if (num_type.name.len <= signature.len and
                std.mem.eql(u8, signature[0..num_type.name.len], num_type.name))
            {
                // Consume matched part
                signature = signature[num_type.name.len..];
                matched_type = @enumFromInt(num_type.value);

                break;
            }
        }

        if (matched_type == null)
            return false and is_numeric;

        if (is_lower) {
            if (lower_bound != null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Multiple lower bound operators in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );

            lower_bound = matched_type;
            lower_inclusive = is_inclusive;
        } else {
            upper_bound = matched_type;
            upper_inclusive = is_inclusive;
        }
    }

    // Now check if actual_type matches the bounds, if they were specified, or the principal type
    comptime var matched_all_bounds = true;
    comptime var checked = 0;
    if (is_numeric) { // Ensure actual_type is not null
        if (lower_bound != null and upper_bound != null) {
            // Both bounds present
            if (lower_inclusive)
                matched_all_bounds = matched_all_bounds and actual_type.?.ge(lower_bound.?)
            else
                matched_all_bounds = matched_all_bounds and actual_type.?.gt(lower_bound.?);

            if (upper_inclusive)
                matched_all_bounds = matched_all_bounds and actual_type.?.le(upper_bound.?)
            else
                matched_all_bounds = matched_all_bounds and actual_type.?.lt(upper_bound.?);

            checked += 2;
        } else if (lower_bound != null) {
            // Only lower bound present
            if (lower_inclusive)
                matched_all_bounds = matched_all_bounds and actual_type.?.ge(lower_bound.?)
            else
                matched_all_bounds = matched_all_bounds and actual_type.?.gt(lower_bound.?);

            checked += 1;
        } else if (upper_bound != null) {
            // Only upper bound present
            if (upper_inclusive)
                matched_all_bounds = matched_all_bounds and actual_type.?.le(upper_bound.?)
            else
                matched_all_bounds = matched_all_bounds and actual_type.?.lt(upper_bound.?);

            checked += 1;
        } // No bounds present -> either we have a standalone principal, or a category as the principal
    }

    // Now check remaining signature for positive/negative types and categories ("@")
    comptime var matched_any_positive = false; // Must match at least one positive type ("" or "|")
    comptime var matched_any_negative = false; // Must not match any negative type ("!")
    comptime var matched_all_categories = true; // All categories must match
    comptime var finished = false;
    while (!finished) : (checked += 1) {
        signature = std.mem.trimStart(u8, signature, " ");
        if (signature.len == 0 or signature[0] == 0) {
            // No more types to check
            finished = true;
            break;
        }

        // Check for "|" and "!" operators. Invalid if checked == 0
        if (signature[0] == '|' or signature[0] == '!') {
            comptime var positive_operator = true;
            if (signature[0] == '|') {
                if (checked == 0) {
                    @compileError(
                        std.fmt.comptimePrint(
                            "Operator '|' cannot be used before any principal numeric type or bound in type signature '{s}' for parameter '{s}' of function '{s}'.",
                            .{ type_signature, param_name, fn_name },
                        ),
                    );
                }

                positive_operator = true;
                signature = signature[1..];
            } else if (signature[0] == '!') {
                if (checked == 0) {
                    @compileError(
                        std.fmt.comptimePrint(
                            "Operator '!' cannot be used before any principal numeric type or bound in type signature '{s}' for parameter '{s}' of function '{s}'.",
                            .{ type_signature, param_name, fn_name },
                        ),
                    );
                }

                positive_operator = false;
                signature = signature[1..];
            }

            // Update checked count
            checked += 1;

            // Match numeric type
            comptime var matched_type: ?NumericType = null;
            inline for (std.meta.fields(NumericType)) |num_type| {
                if (num_type.name.len <= signature.len and
                    std.mem.eql(u8, signature[0..num_type.name.len], num_type.name))
                {
                    // Consume matched part
                    signature = signature[num_type.name.len..];
                    matched_type = @enumFromInt(num_type.value);
                    break;
                }
            }

            if (matched_type == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized numeric type in numeric type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );

            if (positive_operator) {
                if (is_numeric and matched_type.? == actual_type.?)
                    matched_any_positive = true;
            } else {
                if (is_numeric and matched_type.? == actual_type.?)
                    matched_any_negative = true;
            }

            continue;
        } else if (signature[0] == '@') {
            // Check for "@" category operator. Valid for any checked value
            signature = signature[1..]; // Consume "@"

            // Match numeric category
            comptime var matched_category: ?CategoryOperator = null;
            inline for (std.meta.fields(CategoryOperator)) |num_cat| {
                if (num_cat.name.len <= signature.len and
                    std.mem.eql(u8, signature[0..num_cat.name.len], num_cat.name))
                {
                    // Consume matched part
                    signature = signature[num_cat.name.len..];
                    matched_category = @enumFromInt(num_cat.value);

                    break;
                }
            }

            if (matched_category == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized numeric category in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );

            // Check category match
            const category_match = matched_category.?.match()(T);
            matched_all_categories = matched_all_categories and category_match;
            matched_any_positive = matched_any_positive or category_match;

            continue;
        }

        // Check for no operator. Only valid if checked == 0
        if (checked == 0) {
            // Match numeric type
            comptime var matched_type: ?NumericType = null;
            inline for (std.meta.fields(NumericType)) |num_type| {
                if (num_type.name.len <= signature.len and
                    std.mem.eql(u8, signature[0..num_type.name.len], num_type.name))
                {
                    // Consume matched part
                    signature = signature[num_type.name.len..];
                    matched_type = @enumFromInt(num_type.value);

                    break;
                }
            }

            if (matched_type == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized numeric type in numeric type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );

            if (is_numeric and matched_type.? == actual_type.?)
                matched_any_positive = true;
        } else {
            @compileError(
                std.fmt.comptimePrint(
                    "Missing operator in numeric type signature '{s}' for parameter '{s}' of function '{s}'.",
                    .{ type_signature, param_name, fn_name },
                ),
            );
        }
    }

    return matched_all_bounds and matched_any_positive and !matched_any_negative and matched_all_categories and is_numeric;
}

/// Checks if the input types match the corresponding type signatures,
/// triggering compile-time errors if not. Also checks the relationships
/// between the types according to the provided type rules.
///
/// Parameters
/// ----------
/// comptime type_signatures (`[][]const u8`): The type signatures to check against.
/// See the documentation of `checkParameterType` for the format of type signatures.
/// If only one type signature is provided, it is used for all types.
///
/// comptime type_rules (`[][]const u8`): The type rules to check between the types.
// pub fn checkParameterTypes(
//     comptime type_signatures: [][]const u8,
//     comptime type_rules: [][]const u8,
//     comptime Ts: []type,
//     comptime fn_name: []const u8,
//     comptime param_names: [][]const u8,
// ) void {}

/// Coerces the input types to the smallest type that can represent both types.
///
/// This function takes two types `X` and `Y` and returns the smallest type that
/// can represent both types without loss of information.
///
/// For two matrices, the coerced type will use the order of the denser operand.
/// Density is ranked (most to least) as:
///   `general.Dense`, `symmetric.Dense`/`hermitian.Dense`, `triangular.Dense`,
///   `Banded`, `general.Block`, `symmetric.Block`/`hermitian.Block`,
///   `general.Sparse`, `symmetric.Sparse`/`hermitian.Sparse`, `triangular.Sparse`.
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
    if (comptime X == Y and !isTriangularDenseMatrix(X) and !isPermutationMatrix(X))
        return X;

    switch (comptime domainType(X)) {
        .numeric => switch (comptime domainType(Y)) {
            .numeric => {},
            .vector => switch (comptime vectorType(Y)) {
                .dense => return vector.Dense(Coerce(X, Numeric(Y))), // numeric + dense vector
                .sparse => return vector.Sparse(Coerce(X, Numeric(Y))), // numeric + sparse vector
                .numeric => unreachable,
            },
            .matrix => switch (comptime matrixType(Y)) {
                .dense_general => return matrix.general.Dense(Coerce(X, Numeric(Y)), orderOf(Y)), // numeric +  dense general matrix
                .dense_symmetric => return matrix.symmetric.Dense(Coerce(X, Numeric(Y)), uploOf(Y), orderOf(Y)), // numeric + dense symmetric matrix
                .dense_hermitian => {
                    if (comptime isComplex(X)) {
                        return matrix.general.Dense(Coerce(X, Numeric(Y)), orderOf(Y)); // numeric (complex) + dense hermitian matrix
                    } else {
                        return matrix.hermitian.Dense(Coerce(X, Numeric(Y)), uploOf(Y), orderOf(Y)); // numeric (real) + dense hermitian matrix
                    }
                },
                .dense_triangular => return matrix.triangular.Dense(Coerce(X, Numeric(Y)), uploOf(Y), .non_unit, orderOf(Y)), // numeric + dense triangular matrix
                .sparse_general => return matrix.general.Sparse(Coerce(X, Numeric(Y)), orderOf(Y)), // numeric + sparse general matrix
                .sparse_symmetric => return matrix.symmetric.Sparse(Coerce(X, Numeric(Y)), uploOf(Y), orderOf(Y)), // numeric + sparse symmetric matrix
                .sparse_hermitian => {
                    if (comptime isComplex(X)) {
                        return matrix.general.Sparse(Coerce(X, Numeric(Y)), orderOf(Y)); // numeric (complex) + sparse hermitian matrix
                    } else {
                        return matrix.hermitian.Sparse(Coerce(X, Numeric(Y)), uploOf(Y), orderOf(Y)); // numeric (real) + sparse hermitian matrix
                    }
                },
                .sparse_triangular => return matrix.triangular.Sparse(Coerce(X, Numeric(Y)), uploOf(Y), .non_unit, orderOf(Y)), // numeric + sparse triangular matrix
                .block_general => return matrix.general.Block(Coerce(X, Numeric(Y)), borderOf(Y), orderOf(Y)), // numeric + sparse block general matrix
                .block_symmetric => return matrix.symmetric.Block(Coerce(X, Numeric(Y)), uploOf(Y), borderOf(Y), orderOf(Y)), // numeric + sparse block symmetric matrix
                .block_hermitian => {
                    if (comptime isComplex(X)) {
                        return matrix.general.Block(Coerce(X, Numeric(Y)), borderOf(Y), orderOf(Y)); // numeric (complex) + sparse block hermitian matrix
                    } else {
                        return matrix.hermitian.Block(Coerce(X, Numeric(Y)), uploOf(Y), borderOf(Y), orderOf(Y)); // numeric (real) + sparse block hermitian matrix
                    }
                },
                .diagonal => return matrix.Diagonal(Coerce(X, Numeric(Y))), // numeric + diagonal matrix
                .banded => return matrix.Banded(Coerce(X, Numeric(Y)), orderOf(Y)), // numeric + banded matrix
                .tridiagonal => return matrix.Tridiagonal(Coerce(X, Numeric(Y))), // numeric + tridiagonal matrix
                .permutation => return matrix.general.Sparse(Coerce(X, Numeric(Y)), orderOf(Y)), // numeric + permutation matrix
                .numeric => unreachable,
            },
            .array => switch (comptime arrayType(Y)) {
                .dense => return array.Dense(Coerce(X, Numeric(Y)), orderOf(Y)), // numeric + dense array
                .strided => return array.Dense(Coerce(X, Numeric(Y)), orderOf(Y)), // numeric + strided array
                .sparse => return array.Sparse(Coerce(X, Numeric(Y)), orderOf(Y)), // numeric + sparse array
                .numeric => unreachable,
            },
            .expression => Expression, // numeric + expression
        },
        .vector => switch (comptime vectorType(X)) {
            .dense => switch (comptime domainType(Y)) {
                .numeric => return vector.Dense(Coerce(Numeric(X), Y)), // dense vector + numeric
                .vector => switch (comptime vectorType(Y)) {
                    .dense => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense vector + dense vector
                    .sparse => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense vector + sparse vector
                    .numeric => unreachable,
                },
                .matrix => @compileError("Cannot coerce vector and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense vector + matrix
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense vector + array
                .expression => Expression, // dense vector + expression
            },
            .sparse => switch (comptime domainType(Y)) {
                .numeric => return vector.Sparse(Coerce(Numeric(X), Y)), // sparse vector + numeric
                .vector => switch (comptime vectorType(Y)) {
                    .dense => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse vector + dense vector
                    .sparse => return vector.Sparse(Coerce(Numeric(X), Numeric(Y))), // sparse vector + sparse vector
                    .numeric => unreachable,
                },
                .matrix => @compileError("Cannot coerce vector and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector + matrix
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector + array
                .expression => Expression, // sparse vector + expression
            },
            .numeric => unreachable,
        },
        .matrix => switch (comptime matrixType(X)) {
            .dense_general => switch (comptime domainType(Y)) {
                .numeric => return matrix.general.Dense(Coerce(Numeric(X), Y), orderOf(X)), // dense general matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense general matrix + vector
                .matrix => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense general matrix + matrix
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense general matrix + array
                .expression => Expression, // dense general matrix + expression
            },
            .dense_symmetric => switch (comptime domainType(Y)) {
                .numeric => return matrix.symmetric.Dense(Coerce(Numeric(X), Y), uploOf(X), orderOf(X)), // dense symmetric matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense symmetric matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // dense symmetric matrix + dense general matrix
                    .dense_symmetric => return matrix.symmetric.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // dense symmetric matrix + dense symmetric matrix
                    .dense_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense symmetric matrix (complex) + dense hermitian matrix
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // dense symmetric matrix (real) + dense hermitian matrix
                        }
                    },
                    .sparse_symmetric => return matrix.symmetric.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // dense symmetric matrix + sparse symmetric matrix
                    .sparse_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense symmetric matrix (complex) + sparse hermitian matrix
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // dense symmetric matrix (real) + sparse hermitian matrix
                        }
                    },
                    .block_symmetric => return matrix.symmetric.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // dense symmetric matrix + sparse block symmetric matrix
                    .block_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense symmetric matrix (complex) + sparse block hermitian matrix
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // dense symmetric matrix (real) + sparse block hermitian matrix
                        }
                    },
                    .diagonal => return matrix.symmetric.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // dense symmetric matrix + diagonal matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense symmetric matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense symmetric matrix + array
                .expression => Expression, // dense symmetric matrix + expression
            },
            .dense_hermitian => switch (comptime domainType(Y)) {
                .numeric => {
                    if (comptime isComplex(Y)) {
                        return matrix.general.Dense(Coerce(Numeric(X), Y), orderOf(X)); // dense hermitian matrix + numeric (complex)
                    } else {
                        return matrix.hermitian.Dense(Coerce(Numeric(X), Y), uploOf(X), orderOf(X)); // dense hermitian matrix + numeric (real)
                    }
                },
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense hermitian matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // dense hermitian matrix + dense general matrix
                    .dense_symmetric => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense hermitian matrix + dense symmetric matrix (complex)
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // dense hermitian matrix + dense symmetric matrix (real)
                        }
                    },
                    .dense_hermitian => return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // dense hermitian matrix + dense hermitian matrix
                    .sparse_symmetric => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense hermitian matrix + sparse symmetric matrix (complex)
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // dense hermitian matrix + sparse symmetric matrix (real)
                        }
                    },
                    .sparse_hermitian => return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // dense hermitian matrix + sparse hermitian matrix
                    .block_symmetric => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense hermitian matrix + sparse block symmetric matrix (complex)
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // dense hermitian matrix + sparse block symmetric matrix (real)
                        }
                    },
                    .block_hermitian => return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // dense hermitian matrix + sparse block hermitian matrix
                    .diagonal => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense hermitian matrix + diagonal matrix (complex)
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // dense hermitian matrix + diagonal matrix (real)
                        }
                    },
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense hermitian matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense hermitian matrix + array
                .expression => Expression, // dense hermitian matrix + expression
            },
            .dense_triangular => switch (comptime domainType(Y)) {
                .numeric => return matrix.triangular.Dense(Coerce(Numeric(X), Y), uploOf(X), .non_unit, orderOf(X)), // dense triangular matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense triangular matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // dense triangular matrix + dense general matrix
                    .dense_symmetric => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // dense triangular matrix + dense symmetric matrix
                    .dense_hermitian => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // dense triangular matrix + dense hermitian matrix
                    .dense_triangular => {
                        if (comptime uploOf(X) == uploOf(Y)) {
                            return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)); // dense triangular matrix + dense triangular matrix (same uplo)
                        } else {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense triangular matrix + dense triangular matrix (different uplo)
                        }
                    },
                    .sparse_triangular => {
                        if (comptime uploOf(X) == uploOf(Y)) {
                            return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)); // dense triangular matrix + sparse triangular matrix (same uplo)
                        } else {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense triangular matrix + sparse triangular matrix (different uplo)
                        }
                    },
                    .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense triangular matrix + banded matrix
                    .diagonal => return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)), // dense triangular matrix + diagonal matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense triangular matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense triangular matrix + array
                .expression => Expression, // dense triangular matrix + expression
            },
            .sparse_general => switch (comptime domainType(Y)) {
                .numeric => return matrix.general.Sparse(Coerce(Numeric(X), Y), orderOf(X)), // sparse general matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse general matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse general matrix + dense general matrix
                    .dense_symmetric => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse general matrix + dense symmetric matrix
                    .dense_hermitian => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse general matrix + dense hermitian matrix
                    .dense_triangular => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse general matrix + dense triangular matrix
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix + sparse general matrix
                    .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix + sparse symmetric matrix
                    .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix + sparse hermitian matrix
                    .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse general matrix + sparse block general matrix
                    .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse general matrix + sparse block symmetric matrix
                    .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse general matrix + sparse block hermitian matrix
                    .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix + diagonal matrix
                    .banded => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse general matrix + banded matrix
                    .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix + tridiagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse general matrix + array
                .expression => Expression, // sparse general matrix + expression
            },
            .sparse_symmetric => switch (comptime domainType(Y)) {
                .numeric => return matrix.symmetric.Sparse(Coerce(Numeric(X), Y), uploOf(X), orderOf(X)), // sparse symmetric matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse symmetric matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse symmetric matrix + dense general matrix
                    .dense_symmetric => return matrix.symmetric.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // sparse symmetric matrix + dense symmetric matrix
                    .dense_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // sparse symmetric matrix (complex) + dense hermitian matrix
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // sparse symmetric matrix (real) + dense hermitian matrix
                        }
                    },
                    .dense_triangular => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse symmetric matrix + dense triangular matrix
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse symmetric matrix + sparse general matrix
                    .sparse_symmetric => return matrix.symmetric.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // sparse symmetric matrix + sparse symmetric matrix
                    .sparse_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse symmetric matrix (complex) + sparse hermitian matrix
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // sparse symmetric matrix (real) + sparse hermitian matrix
                        }
                    },
                    .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse symmetric matrix + sparse block general matrix
                    .block_symmetric => return matrix.symmetric.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // sparse symmetric matrix + sparse block symmetric matrix
                    .block_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // sparse symmetric matrix (complex) + sparse block hermitian matrix
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // sparse symmetric matrix (real) + sparse block hermitian matrix
                        }
                    },
                    .diagonal => return matrix.symmetric.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // sparse symmetric matrix + diagonal matrix
                    .banded => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse symmetric matrix + banded matrix
                    .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse symmetric matrix + tridiagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix + permutation matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse symmetric matrix + array
                .expression => Expression, // sparse symmetric matrix + expression
            },
            .sparse_hermitian => switch (comptime domainType(Y)) {
                .numeric => {
                    if (comptime isComplex(Y)) {
                        return matrix.general.Sparse(Coerce(Numeric(X), Y), orderOf(X)); // sparse hermitian matrix + numeric (complex)
                    } else {
                        return matrix.hermitian.Sparse(Coerce(Numeric(X), Y), uploOf(X), orderOf(X)); // sparse hermitian matrix + numeric (real)
                    }
                },
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse hermitian matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse hermitian matrix + dense general matrix
                    .dense_symmetric => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // sparse hermitian matrix + dense symmetric matrix (complex)
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // sparse hermitian matrix + dense symmetric matrix (real)
                        }
                    },
                    .dense_hermitian => return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // sparse hermitian matrix + dense hermitian matrix
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse hermitian matrix + sparse general matrix
                    .sparse_symmetric => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse hermitian matrix + sparse symmetric matrix (complex)
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // sparse hermitian matrix + sparse symmetric matrix (real)
                        }
                    },
                    .sparse_hermitian => return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // sparse hermitian matrix + sparse hermitian matrix
                    .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse hermitian matrix + sparse block general matrix
                    .block_symmetric => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // sparse hermitian matrix (complex) + sparse block symmetric matrix
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // sparse hermitian matrix (real) + sparse block symmetric matrix
                        }
                    },
                    .block_hermitian => return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // sparse hermitian matrix + sparse block hermitian matrix
                    .diagonal => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse hermitian matrix (complex) + diagonal matrix
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // sparse hermitian matrix (real) + diagonal matrix
                        }
                    },
                    .banded => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse hermitian matrix + banded matrix
                    .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix + tridiagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix + permutation matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse hermitian matrix + array
                .expression => Expression, // sparse hermitian matrix + expression
            },
            .sparse_triangular => switch (comptime domainType(Y)) {
                .numeric => return matrix.triangular.Sparse(Coerce(Numeric(X), Y), uploOf(X), .non_unit, orderOf(X)), // sparse triangular matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse triangular matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + dense general matrix
                    .dense_symmetric => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + dense symmetric matrix
                    .dense_hermitian => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + dense hermitian matrix
                    .dense_triangular => {
                        if (comptime uploOf(X) == uploOf(Y)) {
                            return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), .non_unit, orderOf(Y)); // sparse triangular matrix + dense triangular matrix (same uplo)
                        } else {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // sparse triangular matrix + dense triangular matrix (different uplo)
                        }
                    },
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + sparse general matrix
                    .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + sparse symmetric matrix
                    .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + sparse hermitian matrix
                    .sparse_triangular => {
                        if (comptime uploOf(X) == uploOf(Y)) {
                            return matrix.triangular.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)); // sparse triangular matrix + sparse triangular matrix (same uplo)
                        } else {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse triangular matrix + sparse triangular matrix (different uplo)
                        }
                    },
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + sparse block general matrix
                    .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + sparse block symmetric matrix
                    .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + sparse block hermitian matrix
                    .diagonal => return matrix.triangular.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)), // sparse triangular matrix + diagonal matrix
                    .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + banded matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + permutation matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse triangular matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse triangular matrix + array
                .expression => Expression, // sparse triangular matrix + expression
            },
            .block_general => switch (comptime domainType(Y)) {
                .numeric => return matrix.general.Block(Coerce(Numeric(X), Y), borderOf(X), orderOf(X)), // sparse block matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block matrix + dense general matrix
                    .dense_symmetric => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block matrix + dense symmetric matrix
                    .dense_hermitian => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block matrix + dense hermitian matrix
                    .dense_triangular => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block matrix + dense triangular matrix
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + sparse general matrix
                    .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + sparse symmetric matrix
                    .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + sparse hermitian matrix
                    .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + sparse block general matrix
                    .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + sparse block symmetric matrix
                    .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + sparse block hermitian matrix
                    .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + diagonal matrix
                    .banded => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block matrix + banded matrix
                    .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + tridiagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + permutation matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block matrix + array
                .expression => Expression, // sparse block matrix + expression
            },
            .block_symmetric => switch (comptime domainType(Y)) {
                .numeric => return matrix.symmetric.Block(Coerce(Numeric(X), Y), borderOf(X), orderOf(X)), // sparse block symmetric matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block symmetric matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block symmetric matrix + dense general matrix
                    .dense_symmetric => return matrix.symmetric.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // sparse block symmetric matrix + dense symmetric matrix
                    .dense_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // sparse block symmetric matrix (complex) + dense hermitian matrix
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // sparse block symmetric matrix (real) + dense hermitian matrix
                        }
                    },
                    .dense_triangular => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block symmetric matrix + dense triangular matrix
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix + sparse general matrix
                    .sparse_symmetric => return matrix.symmetric.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // sparse block symmetric matrix + sparse symmetric matrix
                    .sparse_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse block symmetric matrix (complex) + sparse hermitian matrix
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // sparse block symmetric matrix (real) + sparse hermitian matrix
                        }
                    },
                    .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block symmetric matrix + sparse block general matrix
                    .block_symmetric => return matrix.symmetric.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // sparse block symmetric matrix + sparse block symmetric matrix
                    .block_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse block symmetric matrix (complex) + sparse block hermitian matrix
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // sparse block symmetric matrix (real) + sparse block hermitian matrix
                        }
                    },
                    .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix + diagonal matrix
                    .banded => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block symmetric matrix + banded matrix
                    .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix + tridiagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix + permutation matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block symmetric matrix + array
                .expression => Expression, // sparse block symmetric matrix + expression
            },
            .block_hermitian => switch (comptime domainType(Y)) {
                .numeric => {
                    if (comptime isComplex(Y)) {
                        return matrix.general.Block(Coerce(Numeric(X), Y), borderOf(X), orderOf(X)); // sparse block hermitian matrix + numeric (complex)
                    } else {
                        return matrix.hermitian.Block(Coerce(Numeric(X), Y), uploOf(X), borderOf(X), orderOf(X)); // sparse block hermitian matrix + numeric (real)
                    }
                },
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block hermitian matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block hermitian matrix + dense general matrix
                    .dense_symmetric => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // sparse block hermitian matrix + dense symmetric matrix (complex)
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // sparse block hermitian matrix + dense symmetric matrix (real)
                        }
                    },
                    .dense_hermitian => return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // sparse block hermitian matrix + dense hermitian matrix
                    .dense_triangular => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block hermitian matrix + dense triangular matrix
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix + sparse general matrix
                    .sparse_symmetric => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse block hermitian matrix + sparse symmetric matrix (complex)
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // sparse block hermitian matrix + sparse symmetric matrix (real)
                        }
                    },
                    .sparse_hermitian => return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // sparse block hermitian matrix + sparse hermitian matrix
                    .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block hermitian matrix + sparse block general matrix
                    .block_symmetric => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse block hermitian matrix (complex) + sparse block symmetric matrix
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // sparse block hermitian matrix (real) + sparse block symmetric matrix
                        }
                    },
                    .block_hermitian => return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)), // sparse block hermitian matrix + sparse block hermitian matrix
                    .diagonal => {
                        if (comptime isComplex(Numeric(Y))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse block hermitian matrix (complex) + diagonal matrix
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), orderOf(X)); // sparse block hermitian matrix (real) + diagonal matrix
                        }
                    },
                    .banded => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // sparse block hermitian matrix + banded matrix
                    .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix + tridiagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix + permutation matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block hermitian matrix + array
                .expression => Expression, // sparse block hermitian matrix + expression
            },
            .diagonal => switch (comptime domainType(Y)) {
                .numeric => return matrix.Diagonal(Coerce(Numeric(X), Y)), // diagonal matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // diagonal matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix + dense general matrix
                    .dense_symmetric => return matrix.symmetric.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // diagonal matrix + dense symmetric matrix
                    .dense_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // diagonal matrix (complex) + dense hermitian matrix
                        } else {
                            return matrix.hermitian.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // diagonal matrix (real) + dense hermitian matrix
                        }
                    },
                    .dense_triangular => return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), .non_unit, orderOf(Y)), // diagonal matrix + dense triangular matrix
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix + sparse general matrix
                    .sparse_symmetric => return matrix.symmetric.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // diagonal matrix + sparse symmetric matrix
                    .sparse_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // diagonal matrix (complex) + sparse hermitian matrix
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // diagonal matrix (real) + sparse hermitian matrix
                        }
                    },
                    .sparse_triangular => return matrix.triangular.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), .non_unit, orderOf(Y)), // diagonal matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix + sparse block general matrix
                    .block_symmetric => return matrix.symmetric.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)), // diagonal matrix + sparse block symmetric matrix
                    .block_hermitian => {
                        if (comptime isComplex(Numeric(X))) {
                            return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // diagonal matrix (complex) + sparse block hermitian matrix
                        } else {
                            return matrix.hermitian.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), orderOf(Y)); // diagonal matrix (real) + sparse block hermitian matrix
                        }
                    },
                    .diagonal => return matrix.Diagonal(Coerce(Numeric(X), Numeric(Y))), // diagonal matrix + diagonal matrix
                    .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix + banded matrix
                    .tridiagonal => return matrix.Tridiagonal(Coerce(Numeric(X), Numeric(Y))), // diagonal matrix + tridiagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // diagonal matrix + permutation matrix
                    .numeric => unreachable,
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // diagonal matrix + array
                .expression => Expression, // diagonal matrix + expression
            },
            .banded => switch (comptime domainType(Y)) {
                .numeric => return matrix.Banded(Coerce(Numeric(X), Y), orderOf(X)), // banded matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // banded matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_general => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // banded matrix + dense general matrix
                    .dense_symmetric => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // banded matrix + dense symmetric matrix
                    .dense_hermitian => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // banded matrix + dense hermitian matrix
                    .dense_triangular => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // banded matrix + dense triangular matrix
                    .sparse_triangular => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix + sparse triangular matrix
                    .diagonal => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix + diagonal matrix
                    .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix + banded matrix
                    .tridiagonal => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix + tridiagonal matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // banded matrix + array
                .expression => Expression, // banded matrix + expression
            },
            .tridiagonal => switch (comptime domainType(Y)) {
                .numeric => return matrix.Tridiagonal(Coerce(Numeric(X), Y)), // tridiagonal matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // tridiagonal matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .dense_triangular => return matrix.Banded(Coerce(Numeric(X), Numeric(Y))), // tridiagonal matrix + dense triangular matrix
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix + sparse general matrix
                    .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix + sparse symmetric matrix
                    .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix + sparse hermitian matrix
                    .sparse_triangular => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix + sparse block general matrix
                    .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix + sparse block symmetric matrix
                    .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix + sparse block hermitian matrix
                    .diagonal => return matrix.Tridiagonal(Coerce(Numeric(X), Numeric(Y))), // tridiagonal matrix + diagonal matrix
                    .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix + banded matrix
                    .tridiagonal => return matrix.Tridiagonal(Coerce(Numeric(X), Numeric(Y))), // tridiagonal matrix + tridiagonal matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // tridiagonal matrix + permutation matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // tridiagonal + array
                .expression => Expression, // tridiagonal + expression
            },
            .permutation => switch (comptime domainType(Y)) {
                .numeric => return matrix.general.Sparse(Coerce(Numeric(X), Y), orderOf(X)), // permutation matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix + vector
                .matrix => switch (comptime matrixType(Y)) {
                    .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + sparse general matrix
                    .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + sparse symmetric matrix
                    .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + sparse hermitian matrix
                    .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + sparse triangular matrix
                    .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + sparse block general matrix
                    .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + sparse block symmetric matrix
                    .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + sparse block hermitian matrix
                    .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + permutation matrix
                    else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // permutation matrix + rest of matrices
                },
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix + array
                .expression => Expression, // permutation matrix + expression
            },
            .numeric => unreachable,
        },
        .array => switch (comptime arrayType(X)) {
            .dense => switch (comptime domainType(Y)) {
                .numeric => return array.Dense(Coerce(Numeric(X), Y), orderOf(X)), // dense + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense + matrix
                .array => return array.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense + array
                .expression => Expression, // dense + expression
            },
            .strided => switch (comptime domainType(Y)) {
                .numeric => return array.Dense(Coerce(Numeric(X), Y), orderOf(X)), // strided + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // strided + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // strided + matrix
                .array => return array.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // strided + array
                .expression => Expression, // strided + expression
            },
            .sparse => switch (comptime domainType(Y)) {
                .numeric => return array.Sparse(Coerce(Numeric(X), Y), orderOf(X)), // sparse + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse + matrix
                .array => switch (comptime arrayType(Y)) {
                    .sparse => return array.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse + sparse
                    else => return array.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse + rest of arrays
                },
                .expression => Expression, // sparse + expression
            },
            .numeric => unreachable,
        },
        .expression => Expression, // expression + anything
    }

    // Else two numeric types
    const xnumeric = numericType(X);
    const ynumeric = numericType(Y);

    switch (xnumeric) {
        .bool => switch (ynumeric) {
            .bool => return bool,
            else => return Y,
        },
        .int => switch (ynumeric) {
            .bool => return X,
            .int => {
                if (X == comptime_int) {
                    return Y;
                } else if (Y == comptime_int) {
                    return X;
                }

                comptime var xinfo = @typeInfo(X);
                comptime var yinfo = @typeInfo(Y);

                if (xinfo.int.signedness == .unsigned) {
                    if (yinfo.int.signedness == .unsigned) {
                        if (xinfo.int.bits > yinfo.int.bits) {
                            return X;
                        } else {
                            return Y;
                        }
                    } else {
                        if (xinfo.int.bits > yinfo.int.bits) {
                            if (xinfo.int.bits == 128) {
                                @compileError("Cannot coerce " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " as " ++ @typeName(X) ++ " already has the maximum amount of bits.");
                            } else {
                                yinfo.int.bits = xinfo.int.bits * 2;
                                return @Type(yinfo);
                            }
                        } else if (xinfo.int.bits == yinfo.int.bits) {
                            yinfo.int.bits *= 2;
                            return @Type(yinfo);
                        } else {
                            return Y;
                        }
                    }
                } else {
                    if (yinfo.int.signedness == .unsigned) {
                        if (yinfo.int.bits > xinfo.int.bits) {
                            if (yinfo.int.bits == 128) {
                                @compileError("Cannot coerce " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " as " ++ @typeName(Y) ++ " already has the maximum amount of bits.");
                            } else {
                                xinfo.int.bits = yinfo.int.bits * 2;
                                return @Type(xinfo);
                            }
                        } else if (xinfo.int.bits == yinfo.int.bits) {
                            xinfo.int.bits *= 2;
                            return @Type(xinfo);
                        } else {
                            return X;
                        }
                    } else {
                        if (xinfo.int.bits > yinfo.int.bits) {
                            return X;
                        } else {
                            return Y;
                        }
                    }
                }
            },
            .float => {
                if (Y == comptime_float) {
                    return EnsureFloat(X);
                }

                return Y;
            },
            // .cfloat -> internal float management
            else => return Y,
        },
        .float => {
            switch (ynumeric) {
                .bool, .int => {
                    if (X == comptime_float) {
                        return EnsureFloat(Y);
                    }

                    return X;
                },
                .float => {
                    if (X == comptime_float) {
                        return Y;
                    } else if (Y == comptime_float) {
                        return X;
                    }

                    const xinfo = @typeInfo(X);
                    const yinfo = @typeInfo(Y);

                    if (xinfo.float.bits > yinfo.float.bits) {
                        return X;
                    } else {
                        return Y;
                    }
                },
                .cfloat => {
                    const xinfo = @typeInfo(X);
                    const yinfo = @typeInfo(Scalar(Y));

                    if (xinfo.float.bits > yinfo.float.bits) {
                        return cfloat.Cfloat(X);
                    } else {
                        return Y;
                    }
                },
                else => return Y,
            }
        },
        .cfloat => {
            switch (ynumeric) {
                .bool => return X,
                .int => return X,
                .float => {
                    const xinfo = @typeInfo(Scalar(X));
                    const yinfo = @typeInfo(Y);

                    if (xinfo.float.bits > yinfo.float.bits) {
                        return X;
                    } else {
                        return cfloat.Cfloat(Y);
                    }
                },
                .cfloat => {
                    const x: X = .{ .re = 0, .im = 0 };
                    const y: Y = .{ .re = 0, .im = 0 };
                    return cfloat.Cfloat(Coerce(@TypeOf(@field(x, "re")), @TypeOf(@field(y, "re"))));
                },
                .integer => return Complex(Rational),
                .rational => return Complex(Rational),
                .real => return Complex(Real),
                .complex => return Y,
            }
        },
        .integer => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return Complex(Rational),
            .integer => return X,
            else => return Y,
        },
        .rational => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return Complex(Rational),
            .integer => return X,
            .rational => return X,
            .complex => return Y,
            else => return Y,
        },
        .real => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return Complex(Real),
            .integer => return X,
            .rational => return X,
            .real => return X,
            .complex => {
                if (Y == Complex(Rational)) {
                    return Complex(Real);
                } else {
                    return Y;
                }
            },
            else => return Y,
        },
        .complex => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return X,
            .integer => return X,
            .rational => return X,
            .real => return X,
            .complex => {
                const x: X = .empty;
                const y: Y = .empty;
                return Complex(Coerce(@TypeOf(@field(x, "re")), @TypeOf(@field(y, "re"))));
            },
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
    switch (comptime domainType(X)) {
        .numeric => {}, // Same as Coerce
        .vector => switch (comptime vectorType(X)) {
            .dense => switch (comptime domainType(Y)) {
                .numeric => {}, // Same as Coerce
                .vector => Coerce(Numeric(X), Numeric(Y)), // dense vector * vector
                .matrix => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense vector * matrix
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense vector * array
                .expression => Expression, // dense vector * expression
            },
            .sparse => switch (comptime domainType(Y)) {
                .numeric => {}, // Same as Coerce
                .vector => Coerce(Numeric(X), Numeric(Y)), // sparse vector * vector
                .matrix => switch (comptime matrixType(Y)) {
                    .diagonal => return vector.Sparse(Coerce(Numeric(X), Numeric(Y))), // sparse vector * diagonal matrix
                    else => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse vector * rest of matrices
                },
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector * array
                .expression => Expression, // sparse vector * expression
            },
            .numeric => unreachable,
        },
        .matrix => {
            switch (comptime matrixType(X)) {
                .dense_general => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense general matrix * vector
                    .matrix => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense general matrix * matrix
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense general matrix * array
                    .expression => Expression, // dense general matrix * expression
                },
                .dense_symmetric => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense symmetric matrix * vector
                    .matrix => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense symmetric matrix * matrix
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense symmetric matrix * array
                    .expression => Expression, // dense symmetric matrix * expression
                },
                .dense_hermitian => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense hermitian matrix * vector
                    .matrix => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense hermitian matrix * matrix
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense hermitian matrix * array
                    .expression => Expression, // dense hermitian matrix * expression
                },
                .dense_triangular => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense triangular matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .dense_triangular => {
                            if (comptime uploOf(X) == uploOf(Y)) {
                                if (comptime diagOf(X) == diagOf(Y)) {
                                    return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), diagOf(X), orderOf(X)); // dense triangular matrix * dense triangular matrix (same uplo and diag)
                                } else {
                                    return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)); // dense triangular matrix * dense triangular matrix (same uplo, different diag)
                                }
                            } else {
                                return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // dense triangular matrix * dense triangular matrix (different uplo)
                            }
                        },
                        .sparse_triangular => {
                            if (comptime uploOf(X) == uploOf(Y)) {
                                if (comptime diagOf(X) == diagOf(Y)) {
                                    return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), diagOf(X), orderOf(Y)); // dense triangular matrix * sparse triangular matrix (same uplo and diag)
                                } else {
                                    return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(Y)); // dense triangular matrix * sparse triangular matrix (same uplo, different diag)
                                }
                            } else {
                                return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // dense triangular matrix * sparse triangular matrix (different uplo)
                            }
                        },
                        .diagonal => return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)), // dense triangular matrix * diagonal matrix
                        .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense triangular matrix * banded matrix
                        .tridiagonal => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense triangular matrix * tridiagonal matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // triangular * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense triangular matrix * array
                    .expression => Expression, // dense triangular matrix * expression
                },
                .sparse_general => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse general matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * sparse hermitian matrix
                        .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * sparse triangular matrix
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * sparse block hermitian matrix
                        .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * diagonal matrix
                        .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * tridiagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse general matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse general matrix * array
                    .expression => Expression, // sparse general matrix * expression
                },
                .sparse_symmetric => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse symmetric matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * sparse hermitian matrix
                        .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * sparse triangular matrix
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * sparse block hermitian matrix
                        .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * diagonal matrix
                        .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * tridiagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse symmetric matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse symmetric matrix * array
                    .expression => Expression, // sparse symmetric matrix * expression
                },
                .sparse_hermitian => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse hermitian matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * sparse hermitian matrix
                        .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * sparse triangular matrix
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * sparse block hermitian matrix
                        .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * diagonal matrix
                        .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * tridiagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse hermitian matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse hermitian matrix * array
                    .expression => Expression, // sparse hermitian matrix * expression
                },
                .sparse_triangular => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse triangular matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .dense_triangular => {
                            if (comptime uploOf(X) == uploOf(Y)) {
                                if (comptime diagOf(X) == diagOf(Y)) {
                                    return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), diagOf(X), orderOf(X)); // sparse triangular matrix * dense triangular matrix (same uplo and diag)
                                } else {
                                    return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)); // sparse triangular matrix * dense triangular matrix (same uplo, different diag)
                                }
                            } else {
                                return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)); // sparse triangular matrix * dense triangular matrix (different uplo)
                            }
                        },
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * sparse hermitian matrix
                        .sparse_triangular => {
                            if (comptime uploOf(X) == uploOf(Y)) {
                                if (comptime diagOf(X) == diagOf(Y)) {
                                    return matrix.triangular.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), diagOf(X), orderOf(Y)); // sparse triangular matrix * sparse triangular matrix (same uplo and diag)
                                } else {
                                    return matrix.triangular.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(Y)); // sparse triangular matrix * sparse triangular matrix (same uplo, different diag)
                                }
                            } else {
                                return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)); // sparse triangular matrix * sparse triangular matrix (different uplo)
                            }
                        },
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * sparse block hermitian matrix
                        .diagonal => return matrix.triangular.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(X), .non_unit, orderOf(X)), // sparse triangular matrix * diagonal matrix
                        .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * banded matrix
                        .tridiagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * tridiagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse triangular matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // triangular * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse triangular matrix * array
                    .expression => Expression, // sparse triangular matrix * expression
                },
                .block_general => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse block general matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * sparse hermitian matrix
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * sparse block hermitian matrix
                        .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block general matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block general matrix * array
                    .expression => Expression, // sparse block general matrix * expression
                },
                .block_symmetric => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse block symmetric matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * sparse hermitian matrix
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * sparse block hermitian matrix
                        .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block symmetric matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block symmetric matrix * array
                    .expression => Expression, // sparse block symmetric matrix * expression
                },
                .block_hermitian => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse block hermitian matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * sparse hermitian matrix
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * sparse block hermitian matrix
                        .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * diagonal matrix
                        .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse block hermitian matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse block hermitian matrix * array
                    .expression => Expression, // sparse block hermitian matrix * expression
                },
                .diagonal => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => switch (comptime vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // diagonal matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(Numeric(X), Numeric(Y))), // diagonal matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime matrixType(Y)) {
                        .dense_triangular => return matrix.triangular.Dense(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), .non_unit, orderOf(Y)), // diagonal matrix * dense triangular matrix
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * sparse hermitian matrix
                        .sparse_triangular => return matrix.triangular.Sparse(Coerce(Numeric(X), Numeric(Y)), uploOf(Y), .non_unit, orderOf(Y)), // diagonal matrix * sparse triangular matrix
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * sparse block hermitian matrix
                        .diagonal => return matrix.Diagonal(Coerce(Numeric(X), Numeric(Y))), // diagonal matrix * diagonal matrix
                        .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * banded matrix
                        .tridiagonal => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // diagonal matrix * tridiagonal matrix. Should be tridiagonal, but since the diagonal can be rectangular and tridiagonal must be square, we return banded
                        .permutation => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // diagonal matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // diagonal matrix * array
                    .expression => Expression, // diagonal matrix * expression
                },
                .banded => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // banded matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .dense_triangular => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix * dense triangular matrix
                        .sparse_triangular => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix * sparse triangular matrix
                        .diagonal => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix * diagonal matrix
                        .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix * banded matrix
                        .tridiagonal => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix * tridiagonal matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // banded matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // banded matrix * array
                    .expression => Expression, // banded matrix * expression
                },
                .tridiagonal => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // tridiagonal matrix * vector
                    .matrix => switch (comptime matrixType(Y)) {
                        .dense_triangular => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * dense triangular matrix
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * sparse hermitian matrix
                        .sparse_triangular => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * sparse triangular matrix
                        .block_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * sparse block general matrix
                        .block_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * sparse block symmetric matrix
                        .block_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * sparse block hermitian matrix
                        .diagonal => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // tridiagonal matrix * diagonal matrix
                        .banded => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * banded matrix
                        .tridiagonal => return matrix.Banded(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * tridiagonal matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(Y)), // tridiagonal matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // tridiagonal matrix * array
                    .expression => Expression, // tridiagonal matrix * array
                },
                .permutation => switch (comptime domainType(Y)) {
                    .numeric => return matrix.general.Sparse(Coerce(Numeric(X), Y), orderOf(X)), // permutation matrix * numeric
                    .vector => switch (comptime vectorType(Y)) {
                        .dense => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // permutation matrix * dense vector
                        .sparse => return vector.Sparse(Coerce(Numeric(X), Numeric(Y))), // permutation matrix * sparse vector
                        .numeric => unreachable,
                    },
                    .matrix => switch (comptime matrixType(Y)) {
                        .sparse_general => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // permutation matrix * sparse general matrix
                        .sparse_symmetric => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // permutation matrix * sparse symmetric matrix
                        .sparse_hermitian => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // permutation matrix * sparse hermitian matrix
                        .sparse_block => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // permutation matrix * sparse block matrix
                        .diagonal => return matrix.general.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // permutation matrix * diagonal matrix
                        .permutation => return matrix.Permutation(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // permutation matrix * permutation matrix
                        else => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // permutation matrix * rest of matrices
                    },
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // permutation matrix * array
                    .expression => Expression, // permutation matrix * expression
                },
                .numeric => unreachable,
            }
        },
        .array => {
            switch (comptime domainType(X)) {
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

/// Checks if if `K` can be coerced to `V` without loss of information. This is
/// a more flexible version of `Coerce`, as it does not require `V` to be to the
/// smallest type that can represent both types. The only requirement is that
/// `V` can represent all values of the first two types.
///
/// Unused right now, kept for possible future use.
pub fn canCoerce(comptime K: type, comptime V: type) bool {
    comptime if (isArray(K) or isArray(V) or
        isMatrix(K) or isMatrix(V) or
        isVector(K) or isVector(V))
    {
        @compileError("Cannot use `canCoerce` with arrays, matrices or vectors.");
    };

    const T1: type = K;
    const T2: type = V;

    const t1numeric = numericType(T1);
    const t2numeric = numericType(T2);

    comptime var T3: type = undefined;

    if (T1 == T2) {
        return true;
    }

    switch (t1numeric) {
        .bool => switch (t2numeric) {
            .bool => T3 = bool,
            else => T3 = T2,
        },
        .int => switch (t2numeric) {
            .bool => T3 = T1,
            .int => {
                if (T1 == comptime_int or T2 == comptime_int) {
                    T3 = comptime_int;
                } else {
                    comptime var t1info = @typeInfo(T1);
                    comptime var t2info = @typeInfo(T2);

                    if (t1info.int.signedness == .unsigned) {
                        if (t2info.int.signedness == .unsigned) {
                            if (t1info.int.bits > t2info.int.bits) {
                                T3 = T1;
                            } else {
                                T3 = T2;
                            }
                        } else {
                            if (t1info.int.bits > t2info.int.bits) {
                                if (t1info.int.bits == 128) {
                                    return false;
                                } else {
                                    t2info.int.bits = t1info.int.bits * 2;
                                    T3 = @Type(t2info);
                                }
                            } else if (t1info.int.bits == t2info.int.bits) {
                                t2info.int.bits *= 2;
                                T3 = @Type(t2info);
                            } else {
                                T3 = T2;
                            }
                        }
                    } else {
                        if (t2info.int.signedness == .unsigned) {
                            if (t2info.int.bits > t1info.int.bits) {
                                if (t2info.int.bits == 128) {
                                    return false;
                                } else {
                                    t1info.int.bits = t2info.int.bits * 2;
                                    T3 = @Type(t1info);
                                }
                            } else if (t1info.int.bits == t2info.int.bits) {
                                t1info.int.bits *= 2;
                                T3 = @Type(t1info);
                            } else {
                                T3 = T1;
                            }
                        } else {
                            if (t1info.int.bits > t2info.int.bits) {
                                T3 = T1;
                            } else {
                                T3 = T2;
                            }
                        }
                    }
                }
            },
            else => T3 = T2,
        },
        .float => {
            switch (t2numeric) {
                .bool => T3 = T1,
                .int => T3 = T1,
                .float => {
                    if (T1 == comptime_float) {
                        T3 = T2;
                    } else {
                        const t1info = @typeInfo(T1);
                        const t2info = @typeInfo(T2);

                        if (t1info.float.bits > t2info.float.bits) {
                            T3 = T1;
                        } else {
                            T3 = T2;
                        }
                    }
                },
                .cfloat => {
                    const t1info = @typeInfo(T1);
                    const t2: T2 = .{ .re = 0, .im = 0 };
                    const t2info = @typeInfo(@TypeOf(@field(t2, "re")));

                    if (t1info.float.bits > t2info.float.bits) {
                        T3 = cfloat.Cfloat(T1);
                    } else {
                        T3 = T2;
                    }
                },
                else => T3 = T2,
            }
        },
        .cfloat => {
            switch (t2numeric) {
                .bool => T3 = T1,
                .int => T3 = T1,
                .float => T3 = T1,
                .cfloat => {
                    const t1: T1 = .{ .re = 0, .im = 0 };
                    const t2: T2 = .{ .re = 0, .im = 0 };
                    T3 = cfloat.Cfloat(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                },
                .integer => T3 = Complex(Rational),
                .rational => T3 = Complex(Rational),
                .real => T3 = Complex(Real),
                .complex => T3 = T2,
                else => T3 = T2,
            }
        },
        .integer => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = Complex(Rational),
            .integer => T3 = T1,
            else => T3 = T2,
        },
        .rational => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = Complex(Rational),
            .integer => T3 = T1,
            .rational => T3 = T1,
            .complex => T3 = T2,
            else => T3 = T2,
        },
        .real => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = Complex(Real),
            .integer => T3 = T1,
            .rational => T3 = T1,
            .real => T3 = T1,
            .complex => {
                if (T2 == Complex(Rational)) {
                    T3 = Complex(Real);
                } else {
                    T3 = T2;
                }
            },
            else => T3 = T2,
        },
        .complex => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => {
                if (T1 == Complex(Integer)) {
                    T3 = Complex(Rational);
                } else {
                    T3 = T1;
                }
            },
            .integer => T3 = T1,
            .rational => T3 = T1,
            .real => T3 = T1,
            .complex => {
                const t1: T1 = .empty;
                const t2: T2 = .empty;
                if (@hasField(T1, "allocator") or @hasField(T2, "allocator")) {
                    T3 = Complex(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                } else {
                    T3 = Complex(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                }
            },
            else => unreachable,
        },
    }

    return T2 == T3;
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
    if (isExpression(X)) {
        return Expression;
    } else if (isArray(X)) {
        switch (arrayType(X)) {
            .dense => return array.Dense(Y, orderOf(X)),
            .strided => return array.Dense(Y, orderOf(X)),
            .sparse => return array.Sparse(Y, orderOf(X)),
            .numeric => unreachable,
        }
    } else if (isMatrix(X)) {
        switch (matrixType(X)) {
            .dense_general => return matrix.general.Dense(Y, orderOf(X)),
            .dense_symmetric => return matrix.symmetric.Dense(Y, uploOf(X), orderOf(X)),
            .dense_hermitian => return matrix.hermitian.Dense(Y, uploOf(X), orderOf(X)),
            .dense_triangular => return matrix.triangular.Dense(Y, uploOf(X), diagOf(X), orderOf(X)),
            .sparse_general => return matrix.general.Sparse(Y, orderOf(X)),
            .sparse_symmetric => return matrix.symmetric.Sparse(Y, uploOf(X), orderOf(X)),
            .sparse_hermitian => return matrix.hermitian.Sparse(Y, uploOf(X), orderOf(X)),
            .sparse_triangular => return matrix.triangular.Sparse(Y, uploOf(X), diagOf(X), orderOf(X)),
            .block_general => return matrix.general.Block(Y, borderOf(X), orderOf(X)),
            .block_symmetric => return matrix.symmetric.Block(Y, uploOf(X), borderOf(X), orderOf(X)),
            .block_hermitian => return matrix.hermitian.Block(Y, uploOf(X), borderOf(X), orderOf(X)),
            .diagonal => return matrix.Diagonal(Y),
            .banded => return matrix.Banded(Y, orderOf(X)),
            .tridiagonal => return matrix.Tridiagonal(Y),
            .permutation => return matrix.Permutation(Y, orderOf(X)),
            .numeric => unreachable,
        }
    } else if (isVector(X)) {
        switch (vectorType(X)) {
            .dense => return vector.Dense(Y),
            .sparse => return vector.Sparse(Y),
            .numeric => unreachable,
        }
    } else {
        return Y;
    }
}

pub fn EnsureVector(comptime X: type, comptime Y: type) type {
    if (isVector(X)) {
        switch (vectorType(X)) {
            .dense => return vector.Dense(Y),
            .sparse => return vector.Sparse(Y),
            .numeric => unreachable,
        }
    } else if (isMatrix(X)) {
        return Y;
    } else if (isArray(X)) {
        return Y;
    } else {
        return Y;
    }
}

pub fn EnsureMatrix(comptime X: type, comptime Y: type) type {
    if (isVector(X)) {
        return Y;
    } else if (isMatrix(X)) {
        switch (matrixType(X)) {
            .dense_general => return matrix.general.Dense(Y, orderOf(X)),
            .dense_symmetric => return matrix.symmetric.Dense(Y, uploOf(X), orderOf(X)),
            .dense_hermitian => return matrix.hermitian.Dense(Y, uploOf(X), orderOf(X)),
            .dense_triangular => return matrix.triangular.Dense(Y, uploOf(X), diagOf(X), orderOf(X)),
            .sparse_general => return matrix.general.Sparse(Y, orderOf(X)),
            .sparse_symmetric => return matrix.symmetric.Sparse(Y, uploOf(X), orderOf(X)),
            .sparse_hermitian => return matrix.hermitian.Sparse(Y, uploOf(X), orderOf(X)),
            .sparse_triangular => return matrix.triangular.Sparse(Y, uploOf(X), diagOf(X), orderOf(X)),
            .block_general => return matrix.general.Block(Y, borderOf(X), orderOf(X)),
            .block_symmetric => return matrix.symmetric.Block(Y, uploOf(X), borderOf(X), orderOf(X)),
            .block_hermitian => return matrix.hermitian.Block(Y, uploOf(X), borderOf(X), orderOf(X)),
            .diagonal => return matrix.Diagonal(Y),
            .banded => return matrix.Banded(Y, orderOf(X)),
            .tridiagonal => return matrix.Tridiagonal(Y),
            .permutation => return matrix.Permutation(Y, orderOf(X)),
            .numeric => unreachable,
        }
    } else if (isArray(X)) {
        return Y;
    } else {
        return Y;
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
    if (isVector(X)) {
        return Y;
    } else if (isMatrix(X)) {
        return Y;
    } else if (isArray(X)) {
        switch (arrayType(X)) {
            .dense => return array.Dense(Y, orderOf(X)),
            .strided => return array.Dense(Y, orderOf(X)),
            .sparse => return array.Sparse(Y, orderOf(X)),
            .numeric => unreachable,
        }
    } else {
        return Y;
    }
}

/// Checks if `X` can be safely cast to `Y`, meaning that the cast will not
/// result in a runtime panic.
///
/// This function checks if the type `X` can be safely cast to type `Y` without
/// losing information or causing a runtime panic. It considers various numeric
/// types and their properties, such as signedness and bit width, to determine
/// if the cast is safe.
///
/// For example, casting a signed integer to an unsigned integer is not
/// considered safe, as it may lead to unexpected results if the value is
/// negative.
///
/// Parameters
/// ----------
/// comptime X (`type`): The type to check if it can be cast safely. Must be a
/// supported numeric type.
///
/// comptime Y (`type`): The type to check if `X` can be cast to. Must be a
/// supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the cast is safe, `false` otherwise.
pub fn canCastSafely(comptime X: type, comptime Y: type) bool {
    if (X == Y) {
        return true;
    }

    const xnumeric = numericType(X);
    const ynumeric = numericType(Y);

    switch (xnumeric) {
        .bool => switch (ynumeric) {
            .bool => return true,
            .int => return true,
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
        },
        .int => switch (ynumeric) {
            .bool => return true,
            .int => {
                const linfo = @typeInfo(X);
                const rinfo = @typeInfo(Y);

                if (linfo.int.signedness == .snsigned and rinfo.int.signedness == .unsigned) {
                    // Casting signed to unsigned is not safe
                    return false;
                }

                if (linfo.int.signedness == .unsigned and rinfo.int.signedness == .signed) {
                    if (linfo.int.bits >= rinfo.int.bits) {
                        // Casting unsigned to signed is not safe if the unsigned type has more or equal bits
                        return false;
                    }
                }

                if (linfo.int.bits > rinfo.int.bits) {
                    // Casting to a smaller integer type is not safe
                    return false;
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
        },
        .float => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting float to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
        },
        .cfloat => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting cfloat to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
        },
        .integer => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting integer to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
        },
        .rational => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting rational to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
        },
        .real => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting real to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
        },
        .complex => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting complex to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
        },
    }
}

/// Coerces the input type to a floating point type if it is not already a
/// higher range type.
pub fn EnsureFloat(comptime T: type) type {
    if (isExpression(T)) {
        return Expression;
    } else if (isArray(T)) {
        switch (arrayType(T)) {
            .dense => return array.Dense(EnsureFloat(Numeric(T)), orderOf(T)),
            .strided => return array.Strided(EnsureFloat(Numeric(T)), orderOf(T)),
            .sparse => return array.Sparse(EnsureFloat(Numeric(T)), orderOf(T)),
            .numeric => unreachable,
        }
    } else if (isMatrix(T)) {
        switch (matrixType(T)) {
            .dense_general => return matrix.general.Dense(EnsureFloat(Numeric(T)), orderOf(T)),
            .dense_symmetric => return matrix.symmetric.Dense(EnsureFloat(Numeric(T)), uploOf(T), orderOf(T)),
            .dense_hermitian => return matrix.hermitian.Dense(EnsureFloat(Numeric(T)), uploOf(T), orderOf(T)),
            .dense_triangular => return matrix.triangular.Dense(EnsureFloat(Numeric(T)), uploOf(T), diagOf(T), orderOf(T)),
            .sparse_general => return matrix.general.Sparse(EnsureFloat(Numeric(T)), orderOf(T)),
            .sparse_symmetric => return matrix.symmetric.Sparse(EnsureFloat(Numeric(T)), uploOf(T), orderOf(T)),
            .sparse_hermitian => return matrix.hermitian.Sparse(EnsureFloat(Numeric(T)), uploOf(T), orderOf(T)),
            .sparse_triangular => return matrix.triangular.Sparse(EnsureFloat(Numeric(T)), uploOf(T), diagOf(T), orderOf(T)),
            .block_general => return matrix.general.Block(EnsureFloat(Numeric(T)), borderOf(T), orderOf(T)),
            .block_symmetric => return matrix.symmetric.Block(EnsureFloat(Numeric(T)), borderOf(T), uploOf(T), orderOf(T)),
            .block_hermitian => return matrix.hermitian.Block(EnsureFloat(Numeric(T)), borderOf(T), uploOf(T), orderOf(T)),
            .diagonal => return matrix.Diagonal(EnsureFloat(Numeric(T))),
            .banded => return matrix.Banded(EnsureFloat(Numeric(T)), orderOf(T)),
            .tridiagonal => return matrix.Tridiagonal(EnsureFloat(Numeric(T))),
            .permutation => return matrix.Permutation(EnsureFloat(Numeric(T)), orderOf(T)),
            .numeric => unreachable,
        }
    } else if (isVector(T)) {
        switch (vectorType(T)) {
            .dense => return vector.Dense(EnsureFloat(Numeric(T))),
            .sparse => return vector.Sparse(EnsureFloat(Numeric(T))),
            .numeric => unreachable,
        }
    }

    switch (numericType(T)) {
        .bool => return default_float,
        .int => return default_float,
        .float => return T,
        .cfloat => return T,
        .integer => return Rational,
        .rational => return T,
        .real => return T,
        .complex => return T,
    }
}

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
                T == matrix.symmetric.Dense(Numeric(T), .upper, .row_major) or
                T == matrix.symmetric.Dense(Numeric(T), .lower, .row_major) or
                T == matrix.hermitian.Dense(Numeric(T), .upper, .row_major) or
                T == matrix.hermitian.Dense(Numeric(T), .lower, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .upper, .non_unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .upper, .unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .non_unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .unit, .row_major) or
                T == matrix.general.Sparse(Numeric(T), .row_major) or
                T == matrix.symmetric.Sparse(Numeric(T), .upper, .row_major) or
                T == matrix.symmetric.Sparse(Numeric(T), .lower, .row_major) or
                T == matrix.hermitian.Sparse(Numeric(T), .upper, .row_major) or
                T == matrix.hermitian.Sparse(Numeric(T), .lower, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .upper, .non_unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .upper, .unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .non_unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .unit, .row_major) or
                T == matrix.general.Block(Numeric(T), .col_major, .row_major) or
                T == matrix.general.Block(Numeric(T), .row_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .upper, .col_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .lower, .col_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .upper, .row_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .lower, .row_major, .row_major) or
                T == matrix.hermitian.Block(Numeric(T), .upper, .col_major, .row_major) or
                T == matrix.hermitian.Block(Numeric(T), .lower, .col_major, .row_major) or
                T == matrix.hermitian.Block(Numeric(T), .upper, .row_major, .row_major) or
                T == matrix.hermitian.Block(Numeric(T), .lower, .row_major, .row_major) or
                T == matrix.Banded(Numeric(T), .row_major))
                return .row_major
            else
                return .col_major;
        } else {
            if (T == matrix.general.Dense(Numeric(T), .row_major) or
                T == matrix.symmetric.Dense(Numeric(T), .upper, .row_major) or
                T == matrix.symmetric.Dense(Numeric(T), .lower, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .upper, .non_unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .upper, .unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .non_unit, .row_major) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .unit, .row_major) or
                T == matrix.general.Sparse(Numeric(T), .row_major) or
                T == matrix.symmetric.Sparse(Numeric(T), .upper, .row_major) or
                T == matrix.symmetric.Sparse(Numeric(T), .lower, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .upper, .non_unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .upper, .unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .non_unit, .row_major) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .unit, .row_major) or
                T == matrix.general.Block(Numeric(T), .col_major, .row_major) or
                T == matrix.general.Block(Numeric(T), .row_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .upper, .col_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .lower, .col_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .upper, .row_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .lower, .row_major, .row_major) or
                T == matrix.Banded(Numeric(T), .row_major))
                return .row_major
            else
                return .col_major;
        }
    } else {
        @compileError("Cannot get order of a non-matrix or non-array type. Use `orderOf` only with matrix or array types.");
    }
}

pub fn borderOf(comptime T: type) Order {
    if (comptime isMatrix(T)) {
        if (comptime isComplex(Numeric(T))) {
            if (T == matrix.general.Block(Numeric(T), .col_major, .row_major) or
                T == matrix.general.Block(Numeric(T), .col_major, .col_major) or
                T == matrix.symmetric.Block(Numeric(T), .upper, .col_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .lower, .col_major, .row_major) or
                T == matrix.hermitian.Block(Numeric(T), .upper, .col_major, .row_major) or
                T == matrix.hermitian.Block(Numeric(T), .lower, .col_major, .row_major))
                return .col_major
            else
                return .row_major;
        } else {
            if (T == matrix.general.Block(Numeric(T), .col_major, .row_major) or
                T == matrix.general.Block(Numeric(T), .col_major, .col_major) or
                T == matrix.symmetric.Block(Numeric(T), .upper, .col_major, .row_major) or
                T == matrix.symmetric.Block(Numeric(T), .lower, .col_major, .row_major))
                return .col_major
            else
                return .row_major;
        }
    } else {
        @compileError("Cannot get block order of a non-matrix type. Use `borderOf` only with block matrix types.");
    }
}

pub fn uploOf(comptime T: type) Uplo {
    if (comptime isMatrix(T)) {
        if (comptime isComplex(Numeric(T))) {
            if (T == matrix.symmetric.Dense(Numeric(T), .lower, orderOf(T)) or
                T == matrix.hermitian.Dense(Numeric(T), .lower, orderOf(T)) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .unit, orderOf(T)) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .non_unit, orderOf(T)) or
                T == matrix.symmetric.Sparse(Numeric(T), .lower, orderOf(T)) or
                T == matrix.hermitian.Sparse(Numeric(T), .lower, orderOf(T)) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .unit, orderOf(T)) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .non_unit, orderOf(T)) or
                T == matrix.symmetric.Block(Numeric(T), .lower, borderOf(T), orderOf(T)) or
                T == matrix.hermitian.Block(Numeric(T), .lower, borderOf(T), orderOf(T)))
                return .lower
            else
                return .upper;
        } else {
            if (T == matrix.symmetric.Dense(Numeric(T), .lower, orderOf(T)) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .unit, orderOf(T)) or
                T == matrix.triangular.Dense(Numeric(T), .lower, .non_unit, orderOf(T)) or
                T == matrix.symmetric.Sparse(Numeric(T), .lower, orderOf(T)) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .unit, orderOf(T)) or
                T == matrix.triangular.Sparse(Numeric(T), .lower, .non_unit, orderOf(T)) or
                T == matrix.symmetric.Block(Numeric(T), .lower, borderOf(T), orderOf(T)))
                return .lower
            else
                return .upper;
        }
    } else {
        @compileError("Cannot get uplo of a non-matrix type. Use `uploOf` only with matrix types.");
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
        @compileError("Cannot get diag of a non-matrix type. Use `diagOf` only with matrix types.");
    }
}

/// Casts a value of any numeric type to any fixed precision numeric type.
///
/// It is a more concise version of `cast` that does not require any options and
/// cannot return an error.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to cast to. Must be a fixed precision numeric
/// type.
///
/// value (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`): The value to cast.
///
/// Returns
/// -------
/// `T`: The value casted to the type `T`.
///
/// Notes
/// -----
/// This function does not check if the cast is safe, which could lead to
/// runtime panics if the value cannot be represented in the target type.
pub inline fn scast(
    comptime T: type,
    value: anytype,
) T {
    const I: type = @TypeOf(value);
    const O: type = T;

    comptime if (!isFixedPrecision(O))
        @compileError("Expected a fixed precision type, but got " ++ @typeName(O));

    if (comptime I == O) {
        return value;
    }

    switch (comptime numericType(I)) {
        .bool => switch (comptime numericType(O)) {
            .bool => unreachable,
            .int => return if (value) 1 else 0,
            .float => return if (value) 1.0 else 0.0,
            .cfloat => return .{
                .re = if (value) 1.0 else 0.0,
                .im = 0.0,
            },
            else => unreachable,
        },
        .int => switch (comptime numericType(O)) {
            .bool => return value != 0,
            .int => return @intCast(value),
            .float => return @floatFromInt(value),
            .cfloat => return .{
                .re = @floatFromInt(value),
                .im = 0.0,
            },
            else => unreachable,
        },
        .float => switch (comptime numericType(O)) {
            .bool => return value != 0.0,
            .int => return @intFromFloat(value),
            .float => return @floatCast(value),
            .cfloat => return if (I == Scalar(O)) .{
                .re = value,
                .im = 0.0,
            } else .{
                .re = @floatCast(value),
                .im = 0.0,
            },
            else => unreachable,
        },
        .cfloat => switch (comptime numericType(O)) {
            .bool => return if (value.re != 0 or value.im != 0) true else false,
            .int => return @intFromFloat(value.re),
            .float => return if (Scalar(I) == O) value.re else @floatCast(value.re),
            .cfloat => return .{
                .re = @floatCast(value.re),
                .im = @floatCast(value.im),
            },
            else => unreachable,
        },
        .integer => switch (comptime numericType(O)) {
            .bool => return integer.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .cfloat => return .{
                .re = value.toFloat(Scalar(O)),
                .im = 0.0,
            },
            else => unreachable,
        },
        .rational => switch (comptime numericType(O)) {
            .bool => return rational.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .cfloat => return .{
                .re = value.toFloat(Scalar(O)),
                .im = 0.0,
            },
            else => unreachable,
        },
        .real => @compileError("Not implemented yet."),
        .complex => switch (comptime numericType(O)) {
            .bool => return complex.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .cfloat => return value.toCFloat(O),
            else => unreachable,
        },
    }
}

/// Casts a value of any numeric type to any other numeric type.
///
/// Parameters
/// ----------
/// comptime `T` (`type`): The type to cast to. Must be a supported numeric type.
///
/// `value` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`):The value to cast.
///
/// `ctx` (`struct`): The context for the cast operation.
///
/// Returns
/// -------
/// `T`: The value casted to the type `T`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value.
///
/// Notes
/// -----
/// This function does not check if the cast is safe, which could lead to
/// runtime panics.
pub inline fn cast(
    comptime T: type,
    value: anytype,
    ctx: anytype,
) !T {
    const I: type = @TypeOf(value);
    const O: type = T;

    if (I == O) {
        switch (comptime numericType(O)) {
            .bool, .int, .float, .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                return value;
            },
            .integer => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime validateContext(@TypeOf(ctx), spec);

                if (getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return value.copy(allocator);
                } else {
                    return value;
                }
            },
            .rational => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime validateContext(@TypeOf(ctx), spec);

                if (getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return value.copy(allocator);
                } else {
                    return value;
                }
            },
            .real => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime validateContext(@TypeOf(ctx), spec);

                if (getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return value.copy(allocator);
                } else {
                    return value;
                }
            },
            .complex => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime validateContext(@TypeOf(ctx), spec);

                if (getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return value.copy(allocator);
                } else {
                    return value;
                }
            },
        }
        return value;
    }

    switch (comptime numericType(I)) {
        .bool => switch (comptime numericType(O)) {
            .bool => unreachable,
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                return if (value) 1 else 0;
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                return if (value) 1.0 else 0.0;
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                return .{
                    .re = if (value) 1.0 else 0.0,
                    .im = 0.0,
                };
            },
            .integer => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime validateContext(@TypeOf(ctx), spec);

                if (getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return if (value)
                        constants.one(Integer, .{ .allocator = allocator })
                    else
                        constants.zero(Integer, .{ .allocator = allocator });
                } else {
                    return if (value)
                        constants.one(Integer, .{}) catch unreachable
                    else
                        constants.zero(Integer, .{}) catch unreachable;
                }
            },
            .rational => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime validateContext(@TypeOf(ctx), spec);

                if (getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return if (value)
                        constants.one(Rational, .{ .allocator = allocator })
                    else
                        constants.zero(Rational, .{ .allocator = allocator });
                } else {
                    return if (value)
                        constants.one(Rational, .{}) catch unreachable
                    else
                        constants.zero(Rational, .{}) catch unreachable;
                }
            },
            .real => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime validateContext(@TypeOf(ctx), spec);

                if (getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return if (value)
                        constants.one(Real, .{ .allocator = allocator })
                    else
                        constants.zero(Real, .{ .allocator = allocator });
                } else {
                    return if (value)
                        constants.one(Real, .{}) catch unreachable
                    else
                        constants.zero(Real, .{}) catch unreachable;
                }
            },
            .complex => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime validateContext(@TypeOf(ctx), spec);

                if (getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return if (value)
                        constants.one(O, .{ .allocator = allocator })
                    else
                        constants.zero(O, .{ .allocator = allocator });
                } else {
                    return if (value)
                        constants.one(O, .{}) catch unreachable
                    else
                        constants.zero(O, .{}) catch unreachable;
                }
            },
        },
        .int => switch (comptime numericType(O)) {
            .bool => {
                comptime validateContext(@TypeOf(ctx), .{});

                return value != 0;
            },
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                return @intCast(value);
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                return @floatFromInt(value);
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                return .{
                    .re = @floatFromInt(value),
                    .im = 0.0,
                };
            },
            .integer => @compileError("Not implemented yet: casting from int to integer"),
            .rational => @compileError("Not implemented yet: casting from int to rational"),
            .real => @compileError("Not implemented yet: casting from int to real"),
            .complex => @compileError("Not implemented yet: casting from int to complex"),
        },
        .float => switch (comptime numericType(O)) {
            .bool => {
                comptime validateContext(@TypeOf(ctx), .{});

                return value != 0.0;
            },
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                return @intFromFloat(value);
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                return @floatCast(value);
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                return if (comptime I == Scalar(O)) .{
                    .re = value,
                    .im = 0.0,
                } else .{
                    .re = @floatCast(value),
                    .im = 0.0,
                };
            },
            .integer => @compileError("Not implemented yet: casting from float to integer"),
            .rational => @compileError("Not implemented yet: casting from float to rational"),
            .real => @compileError("Not implemented yet: casting from float to real"),
            .complex => @compileError("Not implemented yet: casting from float to complex"),
        },
        .cfloat => switch (comptime numericType(O)) {
            .bool => {
                comptime validateContext(@TypeOf(ctx), .{});

                return value.re != 0.0 or value.im != 0.0;
            },
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                return @intFromFloat(value.re);
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                return if (comptime Scalar(I) == O)
                    value.re
                else
                    @floatCast(value.re);
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                return .{
                    .re = @floatCast(value.re),
                    .im = @floatCast(value.im),
                };
            },
            .integer => @compileError("Not implemented yet: casting from cfloat to integer"),
            .rational => @compileError("Not implemented yet: casting from cfloat to rational"),
            .real => @compileError("Not implemented yet: casting from cfloat to real"),
            .complex => @compileError("Not implemented yet: casting from cfloat to complex"),
        },
        .integer => switch (comptime numericType(O)) {
            .bool => @compileError("Not implemented yet: casting from integer to bool"),
            .int => @compileError("Not implemented yet: casting from integer to int"),
            .float => @compileError("Not implemented yet: casting from integer to float"),
            .cfloat => @compileError("Not implemented yet: casting from integer to cfloat"),
            .integer => unreachable,
            .rational => @compileError("Not implemented yet: casting from integer to rational"),
            .real => @compileError("Not implemented yet: casting from integer to real"),
            .complex => @compileError("Not implemented yet: casting from integer to complex"),
        },
        .rational => switch (comptime numericType(O)) {
            .bool => @compileError("Not implemented yet: casting from rational to bool"),
            .int => @compileError("Not implemented yet: casting from rational to int"),
            .float => @compileError("Not implemented yet: casting from rational to float"),
            .cfloat => @compileError("Not implemented yet: casting from rational to cfloat"),
            .integer => @compileError("Not implemented yet: casting from rational to integer"),
            .rational => unreachable,
            .real => @compileError("Not implemented yet: casting from rational to real"),
            .complex => @compileError("Not implemented yet: casting from rational to complex"),
        },
        .real => switch (comptime numericType(O)) {
            .bool => @compileError("Not implemented yet: casting from real to bool"),
            .int => @compileError("Not implemented yet: casting from real to int"),
            .float => @compileError("Not implemented yet: casting from real to float"),
            .cfloat => @compileError("Not implemented yet: casting from real to cfloat"),
            .integer => @compileError("Not implemented yet: casting from real to integer"),
            .rational => @compileError("Not implemented yet: casting from real to rational"),
            .real => unreachable,
            .complex => @compileError("Not implemented yet: casting from real to complex"),
        },
        .complex => switch (comptime numericType(O)) {
            .bool => @compileError("Not implemented yet: casting from complex to bool"),
            .int => @compileError("Not implemented yet: casting from complex to int"),
            .float => @compileError("Not implemented yet: casting from complex to float"),
            .cfloat => @compileError("Not implemented yet: casting from complex to cfloat"),
            .integer => @compileError("Not implemented yet: casting from complex to integer"),
            .rational => @compileError("Not implemented yet: casting from complex to rational"),
            .real => @compileError("Not implemented yet: casting from complex to real"),
            .complex => @compileError("Not implemented yet: casting from complex to complex"),
        },
    }
}

/// Returns the pointer child type of a given pointer type.
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
        else => @compileError("Expected a pointer type, but got " ++ @typeName(T)),
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

/// Validates a context struct against a specification struct.
///
/// This function checks that the context struct matches the specification
/// struct, ensuring that all required fields are present and have the correct
/// types. It also checks for unexpected fields in the context struct.
///
/// Parameters
/// ----------
/// ctx (`anytype`): The context struct to validate. Must be a struct type.
/// spec (`anytype`): The specification struct that defines the expected fields
/// and their types. Must have the following structure:
/// ```zig
/// .{
///     .field_name = .{ .type = type, .required = bool },
///     ...
/// }
/// ```
/// where `field_name` is the name of the field, `type` is the expected type of
/// the field, and `required` is a boolean indicating whether the field is
/// required or not. If `required` is `false`, then another field named
/// `default` must be present, specifying the default value for the field.
///
/// Returns
/// -------
/// `void`: If the context struct is valid according to the specification.
///
/// Notes
/// -----
/// If the specification struct has no fields, the context struct is allowed to
/// be `void`.
pub fn validateContext(comptime Ctx: type, comptime spec: anytype) void {
    const ctxinfo = @typeInfo(Ctx);

    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    if (specinfo.@"struct".fields.len == 0 and Ctx == void)
        return; // No context needed, nothing to validate

    if (ctxinfo != .@"struct")
        @compileError("Expected struct for context, got " ++ @typeName(Ctx));

    // Check that all required fields exist and types match
    inline for (specinfo.@"struct".fields) |field| {
        const field_spec = @field(spec, field.name);
        const required = @hasField(@TypeOf(field_spec), "required") and field_spec.required;
        const expected_type = field_spec.type;

        if (@hasField(Ctx, field.name)) {
            const actual_type = @FieldType(Ctx, field.name);
            const type_info = @typeInfo(expected_type);
            const types_match = if (actual_type == @TypeOf(.enum_literal)) blk: { // Special case for enum literals
                break :blk type_info == .@"enum" or (type_info == .optional and @typeInfo(type_info.optional.child) == .@"enum");
            } else if (type_info == .optional)
                actual_type == type_info.optional.child or actual_type == @TypeOf(null)
            else
                actual_type == expected_type;

            if (!types_match)
                @compileError(formatSpecCtxMismatch(Ctx, spec));
        } else if (required) {
            @compileError(formatSpecCtxMismatch(Ctx, spec));
        }
    }

    // Check for unexpected fields in context
    inline for (ctxinfo.@"struct".fields) |ctx_field| {
        if (!@hasField(SpecType, ctx_field.name)) {
            @compileError(formatSpecCtxMismatch(Ctx, spec));
        }
    }
}

fn formatSpecCtxMismatch(
    comptime Ctx: type,
    comptime spec: anytype,
) []const u8 {
    const ctxinfo = @typeInfo(Ctx);

    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    comptime var spec_str: []const u8 = "";
    if (specinfo.@"struct".fields.len == 0) {
        spec_str = ".{}\n";
    } else {
        comptime var max_name_len: usize = 0;
        comptime var max_type_len: usize = 0;
        inline for (specinfo.@"struct".fields) |field| {
            const field_type = @typeName(@field(@field(spec, field.name), "type"));
            if (field.name.len > max_name_len) max_name_len = field.name.len;
            if (field_type.len > max_type_len) max_type_len = field_type.len;
        }

        spec_str = spec_str ++ ".{\n";
        inline for (specinfo.@"struct".fields) |field| {
            const field_type = @typeName(@field(@field(spec, field.name), "type"));
            const required = @hasField(@TypeOf(@field(spec, field.name)), "required") and @field(@field(spec, field.name), "required");
            spec_str = spec_str ++ std.fmt.comptimePrint(
                "    .{s}: {s}",
                .{
                    field.name,
                    field_type,
                },
            );

            for (0..((max_name_len + max_type_len) - field.name.len - field_type.len)) |_| {
                spec_str = spec_str ++ " ";
            }

            spec_str = spec_str ++ std.fmt.comptimePrint(
                "    ({s})\n",
                .{
                    if (required)
                        "required"
                    else
                        std.fmt.comptimePrint(
                            "optional, default = {}",
                            .{@field(@field(spec, field.name), "default")},
                        ),
                },
            );
        }
        spec_str = spec_str ++ "}\n";
    }

    comptime var ctx_str: []const u8 = "";
    if (ctxinfo.@"struct".fields.len == 0) {
        ctx_str = ".{}\n";
    } else {
        ctx_str = ctx_str ++ ".{\n";
        inline for (ctxinfo.@"struct".fields) |field| {
            const field_type = @typeName(@FieldType(Ctx, field.name));
            ctx_str = ctx_str ++ std.fmt.comptimePrint(
                "    .{s}: {s}\n",
                .{
                    field.name,
                    field_type,
                },
            );
        }
        ctx_str = ctx_str ++ "}\n";
    }

    return std.fmt.comptimePrint(
        "Context struct must have the following structure:\n{s}\nGot:\n{s}",
        .{ spec_str, ctx_str },
    );
}

/// Validates that a context struct contains all required fields as defined
/// by a specification struct.
///
/// This function differs from `validateContext` in that it does not raise
/// a compile error if the context struct contains unexpected fields.
///
/// Parameters
/// ----------
/// ctx (`anytype`): The context struct to validate. Must be a struct type.
///
/// spec (`anytype`): The specification struct that defines the expected fields
/// and their types. Must have the following structure:
/// ```zig
/// .{
///     .field_name = .{ .type = type, .required = bool },
///     ...
/// }
/// ```
/// where `field_name` is the name of the field, `type` is the expected type of
/// the field, and `required` is a boolean indicating whether the field is
/// required or not.
///
/// Returns
/// -------
/// `void`: If the context struct is valid according to the specification.
pub fn partialValidateContext(comptime Ctx: type, comptime spec: anytype) void {
    const ctxinfo = @typeInfo(Ctx);

    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    if (specinfo.@"struct".fields.len == 0 and Ctx == void)
        return; // No context needed, nothing to validate

    if (ctxinfo != .@"struct")
        @compileError("Expected struct for context, got " ++ @typeName(Ctx));

    // Check that all required fields exist and types match
    inline for (specinfo.@"struct".fields) |field| {
        const field_spec = @field(spec, field.name);
        const required = @hasField(@TypeOf(field_spec), "required") and field_spec.required;
        const expected_type = field_spec.type;

        if (@hasField(Ctx, field.name)) {
            const actual_type = @FieldType(Ctx, field.name);
            const types_match = if (actual_type == @TypeOf(.enum_literal)) blk: { // Special case for enum literals
                const type_info = @typeInfo(expected_type);
                break :blk type_info == .@"enum" or (type_info == .optional and @typeInfo(type_info.optional.child) == .@"enum");
            } else actual_type == expected_type;

            if (!types_match)
                @compileError(formatCtxRequiredFields(Ctx, spec));
        } else if (required) {
            @compileError(formatCtxRequiredFields(Ctx, spec));
        }
    }
}

fn formatCtxRequiredFields(
    comptime Ctx: type,
    comptime spec: anytype,
) []const u8 {
    const ctxinfo = @typeInfo(Ctx);

    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    comptime var spec_str: []const u8 = "";
    if (specinfo.@"struct".fields.len == 0) {
        spec_str = ".{}\n";
    } else {
        comptime var max_name_len: usize = 0;
        comptime var max_type_len: usize = 0;
        inline for (specinfo.@"struct".fields) |field| {
            const field_type = @typeName(@field(@field(spec, field.name), "type"));
            if (field.name.len > max_name_len) max_name_len = field.name.len;
            if (field_type.len > max_type_len) max_type_len = field_type.len;
        }

        spec_str = spec_str ++ ".{\n";
        inline for (specinfo.@"struct".fields) |field| {
            const field_type = @typeName(@field(@field(spec, field.name), "type"));
            const required = @hasField(@TypeOf(@field(spec, field.name)), "required") and @field(@field(spec, field.name), "required");
            spec_str = spec_str ++ std.fmt.comptimePrint(
                "    .{s}: {s}",
                .{
                    field.name,
                    field_type,
                },
            );

            for (0..((max_name_len + max_type_len) - field.name.len - field_type.len)) |_| {
                spec_str = spec_str ++ " ";
            }

            spec_str = spec_str ++ std.fmt.comptimePrint(
                "    ({s})\n",
                .{
                    if (required) "required" else "optional",
                },
            );
        }
        spec_str = spec_str ++ "}\n";
    }

    comptime var ctx_str: []const u8 = "";
    if (ctxinfo.@"struct".fields.len == 0) {
        ctx_str = ".{}\n";
    } else {
        ctx_str = ctx_str ++ ".{\n";
        inline for (ctxinfo.@"struct".fields) |field| {
            const field_type = @typeName(@FieldType(Ctx, field.name));
            ctx_str = ctx_str ++ std.fmt.comptimePrint(
                "    .{s}: {s}\n",
                .{
                    field.name,
                    field_type,
                },
            );
        }
        ctx_str = ctx_str ++ "}\n";
    }

    return std.fmt.comptimePrint(
        "Context struct must contain at least the following fields:\n{s}\nGot:\n{s}",
        .{ spec_str, ctx_str },
    );
}

pub fn ctxHasField(
    comptime T: type,
    comptime field_name: []const u8,
    comptime FieldType: type,
) bool {
    const has_field: bool = @hasField(T, field_name);

    if (has_field) {
        if (@FieldType(T, field_name) == FieldType or
            (@FieldType(T, field_name) == @TypeOf(.enum_literal) and
                (@typeInfo(FieldType) == .@"enum" or
                    (@typeInfo(FieldType) == .optional and @typeInfo(@typeInfo(FieldType).optional.child) == .@"enum"))) or
            (@FieldType(T, field_name) == @TypeOf(null) and
                @typeInfo(FieldType) == .optional))
            return true;
    }

    return false;
}

pub fn getFieldOrDefault(ctx: anytype, comptime spec: anytype, comptime field_name: []const u8) @field(spec, field_name).type {
    const T = @TypeOf(ctx);
    const FieldType = @field(spec, field_name).type;

    if (@hasField(T, field_name)) {
        const actual_type = @FieldType(T, field_name);
        const type_info = @typeInfo(FieldType);
        const expected_type = FieldType;
        const types_match = if (actual_type == @TypeOf(.enum_literal)) blk: { // Special case for enum literals
            break :blk type_info == .@"enum" or (type_info == .optional and @typeInfo(type_info.optional.child) == .@"enum");
        } else if (type_info == .optional)
            actual_type == type_info.optional.child or actual_type == @TypeOf(null)
        else
            actual_type == expected_type;

        if (!types_match)
            @compileError("Field '" ++ field_name ++ "' has type " ++ @typeName(@FieldType(T, field_name)) ++ ", expected " ++ @typeName(FieldType));

        return @field(ctx, field_name);
    }

    return @field(spec, field_name).default;
}

pub fn MixStructs(comptime S1: type, comptime S2: type) type {
    const info1 = @typeInfo(S1);
    const info2 = @typeInfo(S2);

    if (info1 != .@"struct" or info2 != .@"struct")
        @compileError("MixStructs: both types must be structs");

    for (info1.@"struct".fields) |field| {
        if (@hasField(S2, field.name))
            @compileError("Field '" ++ field.name ++ "' already exists in the second struct");
    }

    return @Type(.{
        .@"struct" = .{
            .layout = .auto,
            .fields = info1.@"struct".fields ++ info2.@"struct".fields,
            .decls = &.{},
            .is_tuple = false,
        },
    });
}

pub fn mixStructs(s1: anytype, s2: anytype) MixStructs(@TypeOf(s1), @TypeOf(s2)) {
    const S1 = @TypeOf(s1);
    const S2 = @TypeOf(s2);

    const info1 = @typeInfo(S1);
    const info2 = @typeInfo(S2);

    if (info1 != .@"struct" or info2 != .@"struct")
        @compileError("mixStructs: both types must be structs");

    var result: MixStructs(S1, S2) = undefined;
    inline for (info1.@"struct".fields) |field| {
        if (@hasField(S2, field.name))
            @compileError("Field '" ++ field.name ++ "' already exists in the second struct");

        @field(result, field.name) = @field(s1, field.name);
    }
    inline for (info2.@"struct".fields) |field| {
        @field(result, field.name) = @field(s2, field.name);
    }

    return result;
}

pub fn StripStruct(comptime S: type, comptime fields_to_remove: []const []const u8) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    // Calculate how many fields will remain
    comptime var remaining_count: usize = 0;
    inline for (info.@"struct".fields) |field| {
        comptime var should_remove: bool = false;
        inline for (fields_to_remove) |field_to_remove| {
            if (std.mem.eql(u8, field.name, field_to_remove)) {
                should_remove = true;
                break;
            }
        }

        if (!should_remove) {
            remaining_count += 1;
        }
    }

    // Create new fields array
    comptime var new_fields: [remaining_count]std.builtin.Type.StructField = undefined;
    comptime var new_field_index: usize = 0;
    inline for (info.@"struct".fields) |field| {
        var should_remove = false;
        inline for (fields_to_remove) |field_to_remove| {
            if (std.mem.eql(u8, field.name, field_to_remove)) {
                should_remove = true;
                break;
            }
        }

        if (!should_remove) {
            new_fields[new_field_index] = field;
            new_field_index += 1;
        }
    }

    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &new_fields,
        .decls = &.{},
        .is_tuple = false,
    } });
}

pub fn stripStruct(s: anytype, comptime fields_to_remove: []const []const u8) StripStruct(@TypeOf(s), fields_to_remove) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    var result: StripStruct(S, fields_to_remove) = undefined;
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        @field(result, field.name) = @field(s, field.name);
    }

    return result;
}

pub fn RenameStructFields(comptime S: type, comptime fields_to_rename: anytype) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    const F: type = @TypeOf(fields_to_rename);
    const finfo = @typeInfo(F);

    if (finfo != .@"struct")
        @compileError("fields_to_rename must be a struct");

    // Check that all new names do not exist in the original struct
    inline for (finfo.@"struct".fields) |field| {
        if (@hasField(S, @field(fields_to_rename, field.name)))
            @compileError("Field '" ++ @field(fields_to_rename, field.name) ++ "' already exists in the struct");
    }

    // Create new fields array
    comptime var new_fields: [info.@"struct".fields.len]std.builtin.Type.StructField = undefined;
    comptime var new_field_index: usize = 0;
    inline for (info.@"struct".fields) |field| {
        comptime var renamed = false;
        inline for (finfo.@"struct".fields) |ffield| {
            if (comptime std.mem.eql(u8, field.name, ffield.name)) {
                new_fields[new_field_index] = .{
                    .name = @as([:0]const u8, @field(fields_to_rename, ffield.name)),
                    .type = field.type,
                    .default_value_ptr = field.default_value_ptr,
                    .is_comptime = field.is_comptime,
                    .alignment = field.alignment,
                };
                renamed = true;
                break;
            }
        }

        if (!renamed) {
            new_fields[new_field_index] = field;
        }
        new_field_index += 1;
    }

    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &new_fields,
        .decls = &.{},
        .is_tuple = false,
    } });
}

/// Renames fields of a struct according to a mapping provided in another struct.
///
/// If a field specified in `fields_to_rename` does not exist in `s`, it is
/// ignored. Fields in `s` that are not specified in `fields_to_rename` retain
/// their original names.
///
/// Parameters
/// ----------
/// `s` (`anytype`): The struct instance whose fields are to be renamed.
/// Must be a struct type.
///
/// `fields_to_rename` (`anytype`): A struct instance that defines the mapping
/// of old field names to new field names. Each field in this struct should
/// have the name of the field to be renamed in `s`, and its value
/// should be a string representing the new name. I.e.:
/// ```zig
/// .{
///     .old_field_name1 = "new_field_name1",
///     .old_field_name2 = "new_field_name2",
///     ...
/// }
/// ```
///
/// Returns
/// -------
/// The new struct instance with fields renamed according to the mapping.
pub fn renameStructFields(s: anytype, comptime fields_to_rename: anytype) RenameStructFields(@TypeOf(s), fields_to_rename) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    const F: type = @TypeOf(fields_to_rename);
    const finfo = @typeInfo(F);

    var result: RenameStructFields(S, fields_to_rename) = undefined;
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        comptime var renamed = false;
        inline for (finfo.@"struct".fields) |ffield| {
            if (comptime std.mem.eql(u8, field.name, @field(fields_to_rename, ffield.name))) {
                @field(result, field.name) = @field(s, ffield.name);
                renamed = true;
                break;
            }
        }

        if (!renamed) {
            @field(result, field.name) = @field(s, field.name);
        }
    }

    return result;
}

pub fn KeepStructFields(comptime S: type, comptime fields_to_keep: []const []const u8) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    // Create new fields array
    comptime var temp_new_fields: [fields_to_keep.len]std.builtin.Type.StructField = undefined;
    comptime var real_num_fields: comptime_int = 0;
    inline for (fields_to_keep) |field_to_keep| {
        comptime var field_info: std.builtin.Type.StructField = undefined;
        inline for (info.@"struct".fields) |field| {
            if (comptime std.mem.eql(u8, field.name, field_to_keep)) {
                field_info = field;
                temp_new_fields[real_num_fields] = field_info;
                real_num_fields += 1;
                break;
            }
        }
    }

    comptime var new_fields: [real_num_fields]std.builtin.Type.StructField = undefined;
    inline for (0..real_num_fields) |i| {
        new_fields[i] = temp_new_fields[i];
    }

    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &new_fields,
        .decls = &.{},
        .is_tuple = false,
    } });
}

/// Creates a new struct type by keeping only the specified fields from the
/// original struct type.
///
/// The fields specified in `fields_to_keep` need not exist in the original
/// struct; any non-existing fields will simply be ignored.
///
/// Parameters
/// ----------
/// `s` (`anytype`): The struct instance from which to keep fields. Must be a
/// struct type.
///
/// `fields_to_keep` (`[]const []const u8`): An array of field names to keep
/// in the new struct type.
///
/// Returns
/// -------
/// The new struct instance containing only the specified fields.
pub fn keepStructFields(s: anytype, comptime fields_to_keep: []const []const u8) KeepStructFields(@TypeOf(s), fields_to_keep) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    var result: KeepStructFields(S, fields_to_keep) = undefined;
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        @field(result, field.name) = @field(s, field.name);
    }

    return result;
}
