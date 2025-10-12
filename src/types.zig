const std = @import("std");

pub const default_uint = usize;
pub const default_int = isize;
pub const default_float = f64;

const cfloat = @import("cfloat.zig");
const cf16 = @import("cfloat.zig").cf16;
const cf32 = @import("cfloat.zig").cf32;
const cf64 = @import("cfloat.zig").cf64;
const cf80 = @import("cfloat.zig").cf80;
const cf128 = @import("cfloat.zig").cf128;
const comptime_complex = @import("cfloat.zig").comptime_complex;
const integer = @import("integer.zig");
const Integer = integer.Integer;
const rational = @import("rational.zig");
const Rational = rational.Rational;
const real = @import("real.zig");
const Real = real.Real;
const complex = @import("complex.zig");
const Complex = complex.Complex;
//pub const Expression = @import("../expression/expression.zig").Expression;

const vector = @import("vector.zig");
const matrix = @import("matrix.zig");
const array = @import("array.zig");

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
/// - `cfloat`: Represents complex floating-point types:
///   - `cf16`, `cf32`, `cf64`, `cf80`, `cf128`
///   - `comptime_complex`
/// - `integer`: Represents arbitrary precision integer type (`Integer`).
/// - `rational`: Represents arbitrary precision rational type (`Rational`).
/// - `real`: Represents arbitrary precision real type (`Real`).
/// - `complex`: Represents complex arbitrary precision types:
///   - `Complex(Integer)`
///   - `Complex(Rational)`
///   - `Complex(Real)`
/// - `expression`: Represents symbolic expressions (Expression).
pub const NumericType = enum {
    bool,
    int,
    float,
    cfloat,
    integer,
    rational,
    real,
    complex,
    expression,
};

const supported_numeric_types: [34]type = .{
    bool,              u8,
    u16,               u32,
    u64,               u128,
    usize,             c_uint,
    i8,                i16,
    i32,               i64,
    i128,              isize,
    c_int,             comptime_int,
    f16,               f32,
    f64,               f80,
    f128,              comptime_float,
    cf16,              cf32,
    cf64,              cf80,
    cf128,             comptime_complex,
    Integer,           Rational,
    Real,              Complex(Integer),
    Complex(Rational),
    Complex(Real),
    // Expression,
};

const supported_complex_types: [9]type = .{
    cf16,             cf32,
    cf64,             cf80,
    cf128,            comptime_complex,
    Complex(Integer), Complex(Rational),
    Complex(Real), // Expression,
};

pub const VectorType = enum {
    dense,
    sparse,
    numeric, // Fallback for numeric types that are not vectors
};

pub const MatrixType = enum {
    dense_general,
    dense_symmetric,
    dense_hermitian,
    dense_triangular,
    sparse_general,
    sparse_symmetric,
    sparse_hermitian,
    sparse_triangular,
    block_general,
    block_symmetric,
    block_hermitian,
    diagonal,
    banded,
    tridiagonal,
    permutation,
    numeric, // Fallback for numeric types that are not matrices
};

pub const ArrayType = enum {
    dense,
    strided,
    sparse,
    numeric, // Fallback for numeric types that are not arrays
};

pub const Domain = enum {
    numeric,
    vector,
    matrix,
    array,
};

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

/// Checks the the input type `T` and returns the corresponding `NumericType`.
///
/// Checks that the input type is a supported numeric type and returns the
/// corresponding `NumericType` enum value. If the type is not supported, it
/// will raise a compile error.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `NumericType`: The corresponding `NumericType` enum value.
pub inline fn numericType(comptime T: type) NumericType {
    // Without inline functions calling this fail miserably. I have no idea why.

    @setEvalBranchQuota(10000000);

    switch (@typeInfo(T)) {
        .bool => return .bool,
        .int, .comptime_int => {
            if (T != u8 and T != u16 and T != u32 and T != u64 and T != u128 and T != usize and T != c_uint and T != i8 and T != i16 and T != i32 and T != i64 and T != i128 and T != isize and T != c_int and T != comptime_int)
                @compileError("Unsupported integer type: " ++ @typeName(T));

            return .int;
        },
        .float, .comptime_float => return .float,
        else => {
            if (T == cf16 or T == cf32 or T == cf64 or T == cf80 or T == cf128 or T == comptime_complex or
                T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or
                T == std.math.Complex(f80) or T == std.math.Complex(f128) or T == std.math.Complex(comptime_float))
            {
                return .cfloat;
            } else if (T == Integer) {
                return .integer;
            } else if (T == Rational) {
                return .rational;
            } else if (T == Real) {
                return .real;
            } else if (T == Complex(Integer) or T == Complex(Rational) or T == Complex(Real)) {
                return .complex;
                //} else if (T == Expression) {
                //    return .expression;
            } else {
                @compileError("Unsupported numeric type: " ++ @typeName(T));
            }
        },
    }
}

pub fn isNumeric(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .bool => return true,
        .int, .comptime_int => {
            if (T != u8 and T != u16 and T != u32 and T != u64 and T != u128 and T != usize and T != c_uint and T != i8 and T != i16 and T != i32 and T != i64 and T != i128 and T != isize and T != c_int and T != comptime_int)
                return false;

            return true;
        },
        .float, .comptime_float => return true,
        else => {
            if (T == cf16 or T == cf32 or T == cf64 or T == cf80 or T == cf128 or T == comptime_complex or
                T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or
                T == std.math.Complex(f80) or T == std.math.Complex(f128) or T == std.math.Complex(comptime_float))
            {
                return true;
            } else if (T == Integer) {
                return true;
            } else if (T == Rational) {
                return true;
            } else if (T == Real) {
                return true;
            } else if (T == Complex(Integer) or T == Complex(Rational) or T == Complex(Real)) {
                return true;
                //} else if (T == Expression) {
                //    return .expression;
            } else {
                return false;
            }
        },
    }
}

pub inline fn vectorType(comptime T: type) VectorType {
    @setEvalBranchQuota(10000);

    if (isDenseVector(T)) {
        return .dense;
    } else if (isSparseVector(T)) {
        return .sparse;
    }

    return .numeric; // Fallback for numeric types that are not vectors
}

pub inline fn matrixType(comptime T: type) MatrixType {
    @setEvalBranchQuota(10000);

    if (isGeneralDenseMatrix(T)) {
        return .dense_general;
    } else if (isSymmetricDenseMatrix(T)) {
        return .dense_symmetric;
    } else if (isHermitianDenseMatrix(T)) {
        return .dense_hermitian;
    } else if (isTriangularDenseMatrix(T)) {
        return .dense_triangular;
    } else if (isGeneralSparseMatrix(T)) {
        return .sparse_general;
    } else if (isSymmetricSparseMatrix(T)) {
        return .sparse_symmetric;
    } else if (isHermitianSparseMatrix(T)) {
        return .sparse_hermitian;
    } else if (isTriangularSparseMatrix(T)) {
        return .sparse_triangular;
    } else if (isGeneralBlockMatrix(T)) {
        return .block_general;
    } else if (isSymmetricBlockMatrix(T)) {
        return .block_symmetric;
    } else if (isHermitianBlockMatrix(T)) {
        return .block_hermitian;
    } else if (isDiagonalMatrix(T)) {
        return .diagonal;
    } else if (isBandedMatrix(T)) {
        return .banded;
    } else if (isTridiagonalMatrix(T)) {
        return .tridiagonal;
    } else if (isPermutationMatrix(T)) {
        return .permutation;
    }

    return .numeric; // Fallback for numeric types that are not matrices
}

pub inline fn arrayType(comptime T: type) ArrayType {
    @setEvalBranchQuota(10000);

    if (isDenseArray(T)) {
        return .dense;
    } else if (isStridedArray(T)) {
        return .strided;
    } else if (isSparseArray(T)) {
        return .sparse;
    }

    return .numeric; // Fallback for numeric types that are not arrays
}

pub inline fn domainType(comptime T: type) Domain {
    @setEvalBranchQuota(10000);

    if (comptime isNumeric(T)) {
        return .numeric;
    } else if (comptime isVector(T)) {
        return .vector;
    } else if (comptime isMatrix(T)) {
        return .matrix;
    } else if (comptime isArray(T)) {
        return .array;
    }

    @compileError("Unsupported type for domainType: " ++ @typeName(T));
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
    @setEvalBranchQuota(10000);
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == vector.Dense(numeric_type) or
                    T == vector.Sparse(numeric_type)) return true;
            }

            return false;
        },
        else => return false,
    }
}

pub fn isDenseVector(comptime T: type) bool {
    @setEvalBranchQuota(10000);
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == vector.Dense(numeric_type)) return true;
            }

            return false;
        },
        else => return false,
    }
}

pub fn isSparseVector(comptime T: type) bool {
    @setEvalBranchQuota(10000);
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == vector.Sparse(numeric_type)) return true;
            }

            return false;
        },
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
    @setEvalBranchQuota(100000);
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                const matrix_types: [if (isComplex(numeric_type)) 61 else 45]type = if (comptime isComplex(numeric_type)) .{
                    matrix.general.Dense(numeric_type, .col_major),
                    matrix.general.Dense(numeric_type, .row_major),
                    matrix.general.Sparse(numeric_type, .col_major),
                    matrix.general.Sparse(numeric_type, .row_major),
                    matrix.general.Block(numeric_type, .col_major, .col_major),
                    matrix.general.Block(numeric_type, .row_major, .col_major),
                    matrix.general.Block(numeric_type, .col_major, .row_major),
                    matrix.general.Block(numeric_type, .row_major, .row_major),
                    matrix.symmetric.Dense(numeric_type, .upper, .col_major),
                    matrix.symmetric.Dense(numeric_type, .lower, .col_major),
                    matrix.symmetric.Dense(numeric_type, .upper, .row_major),
                    matrix.symmetric.Dense(numeric_type, .lower, .row_major),
                    matrix.symmetric.Sparse(numeric_type, .upper, .col_major),
                    matrix.symmetric.Sparse(numeric_type, .lower, .col_major),
                    matrix.symmetric.Sparse(numeric_type, .upper, .row_major),
                    matrix.symmetric.Sparse(numeric_type, .lower, .row_major),
                    matrix.symmetric.Block(numeric_type, .upper, .col_major, .col_major),
                    matrix.symmetric.Block(numeric_type, .lower, .col_major, .col_major),
                    matrix.symmetric.Block(numeric_type, .upper, .row_major, .col_major),
                    matrix.symmetric.Block(numeric_type, .lower, .row_major, .col_major),
                    matrix.symmetric.Block(numeric_type, .upper, .col_major, .row_major),
                    matrix.symmetric.Block(numeric_type, .lower, .col_major, .row_major),
                    matrix.symmetric.Block(numeric_type, .upper, .row_major, .row_major),
                    matrix.symmetric.Block(numeric_type, .lower, .row_major, .row_major),
                    matrix.hermitian.Dense(numeric_type, .upper, .col_major),
                    matrix.hermitian.Dense(numeric_type, .lower, .col_major),
                    matrix.hermitian.Dense(numeric_type, .upper, .row_major),
                    matrix.hermitian.Dense(numeric_type, .lower, .row_major),
                    matrix.hermitian.Sparse(numeric_type, .upper, .col_major),
                    matrix.hermitian.Sparse(numeric_type, .lower, .col_major),
                    matrix.hermitian.Sparse(numeric_type, .upper, .row_major),
                    matrix.hermitian.Sparse(numeric_type, .lower, .row_major),
                    matrix.hermitian.Block(numeric_type, .upper, .col_major, .col_major),
                    matrix.hermitian.Block(numeric_type, .lower, .col_major, .col_major),
                    matrix.hermitian.Block(numeric_type, .upper, .row_major, .col_major),
                    matrix.hermitian.Block(numeric_type, .lower, .row_major, .col_major),
                    matrix.hermitian.Block(numeric_type, .upper, .col_major, .row_major),
                    matrix.hermitian.Block(numeric_type, .lower, .col_major, .row_major),
                    matrix.hermitian.Block(numeric_type, .upper, .row_major, .row_major),
                    matrix.hermitian.Block(numeric_type, .lower, .row_major, .row_major),
                    matrix.triangular.Dense(numeric_type, .upper, .non_unit, .col_major),
                    matrix.triangular.Dense(numeric_type, .upper, .unit, .col_major),
                    matrix.triangular.Dense(numeric_type, .lower, .non_unit, .col_major),
                    matrix.triangular.Dense(numeric_type, .lower, .unit, .col_major),
                    matrix.triangular.Dense(numeric_type, .upper, .non_unit, .row_major),
                    matrix.triangular.Dense(numeric_type, .upper, .unit, .row_major),
                    matrix.triangular.Dense(numeric_type, .lower, .non_unit, .row_major),
                    matrix.triangular.Dense(numeric_type, .lower, .unit, .row_major),
                    matrix.triangular.Sparse(numeric_type, .upper, .non_unit, .col_major),
                    matrix.triangular.Sparse(numeric_type, .upper, .unit, .col_major),
                    matrix.triangular.Sparse(numeric_type, .lower, .non_unit, .col_major),
                    matrix.triangular.Sparse(numeric_type, .lower, .unit, .col_major),
                    matrix.triangular.Sparse(numeric_type, .upper, .non_unit, .row_major),
                    matrix.triangular.Sparse(numeric_type, .upper, .unit, .row_major),
                    matrix.triangular.Sparse(numeric_type, .lower, .non_unit, .row_major),
                    matrix.triangular.Sparse(numeric_type, .lower, .unit, .row_major),
                    matrix.Diagonal(numeric_type),
                    matrix.Banded(numeric_type, .col_major),
                    matrix.Banded(numeric_type, .row_major),
                    matrix.Tridiagonal(numeric_type),
                    matrix.Permutation(numeric_type),
                } else .{
                    matrix.general.Dense(numeric_type, .col_major),
                    matrix.general.Dense(numeric_type, .row_major),
                    matrix.general.Sparse(numeric_type, .col_major),
                    matrix.general.Sparse(numeric_type, .row_major),
                    matrix.general.Block(numeric_type, .col_major, .col_major),
                    matrix.general.Block(numeric_type, .row_major, .col_major),
                    matrix.general.Block(numeric_type, .col_major, .row_major),
                    matrix.general.Block(numeric_type, .row_major, .row_major),
                    matrix.symmetric.Dense(numeric_type, .upper, .col_major),
                    matrix.symmetric.Dense(numeric_type, .lower, .col_major),
                    matrix.symmetric.Dense(numeric_type, .upper, .row_major),
                    matrix.symmetric.Dense(numeric_type, .lower, .row_major),
                    matrix.symmetric.Sparse(numeric_type, .upper, .col_major),
                    matrix.symmetric.Sparse(numeric_type, .lower, .col_major),
                    matrix.symmetric.Sparse(numeric_type, .upper, .row_major),
                    matrix.symmetric.Sparse(numeric_type, .lower, .row_major),
                    matrix.symmetric.Block(numeric_type, .upper, .col_major, .col_major),
                    matrix.symmetric.Block(numeric_type, .lower, .col_major, .col_major),
                    matrix.symmetric.Block(numeric_type, .upper, .row_major, .col_major),
                    matrix.symmetric.Block(numeric_type, .lower, .row_major, .col_major),
                    matrix.symmetric.Block(numeric_type, .upper, .col_major, .row_major),
                    matrix.symmetric.Block(numeric_type, .lower, .col_major, .row_major),
                    matrix.symmetric.Block(numeric_type, .upper, .row_major, .row_major),
                    matrix.symmetric.Block(numeric_type, .lower, .row_major, .row_major),
                    matrix.triangular.Dense(numeric_type, .upper, .non_unit, .col_major),
                    matrix.triangular.Dense(numeric_type, .upper, .unit, .col_major),
                    matrix.triangular.Dense(numeric_type, .lower, .non_unit, .col_major),
                    matrix.triangular.Dense(numeric_type, .lower, .unit, .col_major),
                    matrix.triangular.Dense(numeric_type, .upper, .non_unit, .row_major),
                    matrix.triangular.Dense(numeric_type, .upper, .unit, .row_major),
                    matrix.triangular.Dense(numeric_type, .lower, .non_unit, .row_major),
                    matrix.triangular.Dense(numeric_type, .lower, .unit, .row_major),
                    matrix.triangular.Sparse(numeric_type, .upper, .non_unit, .col_major),
                    matrix.triangular.Sparse(numeric_type, .upper, .unit, .col_major),
                    matrix.triangular.Sparse(numeric_type, .lower, .non_unit, .col_major),
                    matrix.triangular.Sparse(numeric_type, .lower, .unit, .col_major),
                    matrix.triangular.Sparse(numeric_type, .upper, .non_unit, .row_major),
                    matrix.triangular.Sparse(numeric_type, .upper, .unit, .row_major),
                    matrix.triangular.Sparse(numeric_type, .lower, .non_unit, .row_major),
                    matrix.triangular.Sparse(numeric_type, .lower, .unit, .row_major),
                    matrix.Diagonal(numeric_type),
                    matrix.Banded(numeric_type, .col_major),
                    matrix.Banded(numeric_type, .row_major),
                    matrix.Tridiagonal(numeric_type),
                    matrix.Permutation(numeric_type),
                };

                inline for (matrix_types) |matrix_type| {
                    if (T == matrix_type) return true;
                }
            }

            return false;
        },
        else => return false,
    }
}

pub fn isSquareMatrix(comptime T: type) bool {
    return isSymmetricDenseMatrix(T) or isHermitianDenseMatrix(T) or
        isTridiagonalMatrix(T) or isSymmetricSparseMatrix(T) or
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.general.Dense(numeric_type, .col_major) or
                    T == matrix.general.Dense(numeric_type, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.symmetric.Dense(numeric_type, .upper, .col_major) or
                    T == matrix.symmetric.Dense(numeric_type, .lower, .col_major) or
                    T == matrix.symmetric.Dense(numeric_type, .upper, .row_major) or
                    T == matrix.symmetric.Dense(numeric_type, .lower, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_complex_types) |numeric_type| {
                if (T == matrix.hermitian.Dense(numeric_type, .upper, .col_major) or
                    T == matrix.hermitian.Dense(numeric_type, .lower, .col_major) or
                    T == matrix.hermitian.Dense(numeric_type, .upper, .row_major) or
                    T == matrix.hermitian.Dense(numeric_type, .lower, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.triangular.Dense(numeric_type, .upper, .non_unit, .col_major) or
                    T == matrix.triangular.Dense(numeric_type, .upper, .unit, .col_major) or
                    T == matrix.triangular.Dense(numeric_type, .lower, .non_unit, .col_major) or
                    T == matrix.triangular.Dense(numeric_type, .lower, .unit, .col_major) or
                    T == matrix.triangular.Dense(numeric_type, .upper, .non_unit, .row_major) or
                    T == matrix.triangular.Dense(numeric_type, .upper, .unit, .row_major) or
                    T == matrix.triangular.Dense(numeric_type, .lower, .non_unit, .row_major) or
                    T == matrix.triangular.Dense(numeric_type, .lower, .unit, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.Diagonal(numeric_type)) return true;
            }

            return false;
        },
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.Banded`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.Banded`, `false` otherwise.
pub fn isBandedMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.Banded(numeric_type, .col_major) or
                    T == matrix.Banded(numeric_type, .row_major)) return true;
            }

            return false;
        },
        else => return false,
    }
}

/// Checks if the input type is an instance of a `matrix.Tridiagonal`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a `matrix.Tridiagonal`, `false` otherwise.
pub fn isTridiagonalMatrix(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.Tridiagonal(numeric_type)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.general.Sparse(numeric_type, .col_major) or
                    T == matrix.general.Sparse(numeric_type, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.symmetric.Sparse(numeric_type, .upper, .col_major) or
                    T == matrix.symmetric.Sparse(numeric_type, .lower, .col_major) or
                    T == matrix.symmetric.Sparse(numeric_type, .upper, .row_major) or
                    T == matrix.symmetric.Sparse(numeric_type, .lower, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_complex_types) |numeric_type| {
                if (T == matrix.hermitian.Sparse(numeric_type, .upper, .col_major) or
                    T == matrix.hermitian.Sparse(numeric_type, .lower, .col_major) or
                    T == matrix.hermitian.Sparse(numeric_type, .upper, .row_major) or
                    T == matrix.hermitian.Sparse(numeric_type, .lower, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.triangular.Sparse(numeric_type, .upper, .non_unit, .col_major) or
                    T == matrix.triangular.Sparse(numeric_type, .upper, .unit, .col_major) or
                    T == matrix.triangular.Sparse(numeric_type, .lower, .non_unit, .col_major) or
                    T == matrix.triangular.Sparse(numeric_type, .lower, .unit, .col_major) or
                    T == matrix.triangular.Sparse(numeric_type, .upper, .non_unit, .row_major) or
                    T == matrix.triangular.Sparse(numeric_type, .upper, .unit, .row_major) or
                    T == matrix.triangular.Sparse(numeric_type, .lower, .non_unit, .row_major) or
                    T == matrix.triangular.Sparse(numeric_type, .lower, .unit, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.general.Block(numeric_type, .col_major, .col_major) or
                    T == matrix.general.Block(numeric_type, .row_major, .col_major) or
                    T == matrix.general.Block(numeric_type, .col_major, .row_major) or
                    T == matrix.general.Block(numeric_type, .row_major, .row_major))
                    return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.symmetric.Block(numeric_type, .upper, .col_major, .col_major) or
                    T == matrix.symmetric.Block(numeric_type, .lower, .col_major, .col_major) or
                    T == matrix.symmetric.Block(numeric_type, .upper, .row_major, .col_major) or
                    T == matrix.symmetric.Block(numeric_type, .lower, .row_major, .col_major) or
                    T == matrix.symmetric.Block(numeric_type, .upper, .col_major, .row_major) or
                    T == matrix.symmetric.Block(numeric_type, .lower, .col_major, .row_major) or
                    T == matrix.symmetric.Block(numeric_type, .upper, .row_major, .row_major) or
                    T == matrix.symmetric.Block(numeric_type, .lower, .row_major, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_complex_types) |numeric_type| {
                if (T == matrix.hermitian.Block(numeric_type, .upper, .col_major, .col_major) or
                    T == matrix.hermitian.Block(numeric_type, .lower, .col_major, .col_major) or
                    T == matrix.hermitian.Block(numeric_type, .upper, .row_major, .col_major) or
                    T == matrix.hermitian.Block(numeric_type, .lower, .row_major, .col_major) or
                    T == matrix.hermitian.Block(numeric_type, .upper, .col_major, .row_major) or
                    T == matrix.hermitian.Block(numeric_type, .lower, .col_major, .row_major) or
                    T == matrix.hermitian.Block(numeric_type, .upper, .row_major, .row_major) or
                    T == matrix.hermitian.Block(numeric_type, .lower, .row_major, .row_major)) return true;
            }

            return false;
        },
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
    @setEvalBranchQuota(10000);

    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == matrix.Permutation(numeric_type)) return true;
            }

            return false;
        },
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
    return isGeneralDenseMatrix(T) or
        isSymmetricDenseMatrix(T) or
        isHermitianDenseMatrix(T) or
        isTriangularDenseMatrix(T);
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
    return isGeneralSparseMatrix(T) or
        isSymmetricSparseMatrix(T) or
        isHermitianSparseMatrix(T) or
        isTriangularSparseMatrix(T) or
        isGeneralBlockMatrix(T) or
        isSymmetricBlockMatrix(T) or
        isHermitianBlockMatrix(T);
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
    @setEvalBranchQuota(10000);
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                const array_types: [6]type = .{
                    array.Dense(numeric_type, .col_major),
                    array.Dense(numeric_type, .row_major),
                    array.Strided(numeric_type, .col_major),
                    array.Strided(numeric_type, .row_major),
                    array.Sparse(numeric_type, .col_major),
                    array.Sparse(numeric_type, .row_major),
                };

                inline for (array_types) |array_type| {
                    if (T == array_type) return true;
                }
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == array.Dense(numeric_type, .col_major) or
                    T == array.Dense(numeric_type, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == array.Strided(numeric_type, .col_major) or
                    T == array.Strided(numeric_type, .row_major)) return true;
            }

            return false;
        },
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
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            inline for (supported_numeric_types) |numeric_type| {
                if (T == array.Sparse(numeric_type, .col_major) or
                    T == array.Sparse(numeric_type, .row_major)) return true;
            }

            return false;
        },
        else => return false,
    }
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
        .expression => return true,
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
            },
            .numeric => unreachable,
        },
        .matrix => switch (comptime matrixType(X)) {
            .dense_general => switch (comptime domainType(Y)) {
                .numeric => return matrix.general.Dense(Coerce(Numeric(X), Y), orderOf(X)), // dense general matrix + numeric
                .vector => @compileError("Cannot coerce matrix and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense general matrix + vector
                .matrix => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense general matrix + matrix
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // general + array
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
                .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // symmetric + array
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
            },
            .numeric => unreachable,
        },
        .array => switch (comptime arrayType(X)) {
            .dense => switch (comptime domainType(Y)) {
                .numeric => return array.Dense(Coerce(Numeric(X), Y), orderOf(X)), // dense + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense + matrix
                .array => return array.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense + array
            },
            .strided => switch (comptime domainType(Y)) {
                .numeric => return array.Dense(Coerce(Numeric(X), Y), orderOf(X)), // strided + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // strided + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // strided + matrix
                .array => return array.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // strided + array
            },
            .sparse => switch (comptime domainType(Y)) {
                .numeric => return array.Sparse(Coerce(Numeric(X), Y), orderOf(X)), // sparse + numeric
                .vector => @compileError("Cannot coerce array and vector types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse + vector
                .matrix => @compileError("Cannot coerce array and matrix types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse + matrix
                .array => switch (comptime arrayType(Y)) {
                    .sparse => return array.Sparse(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse + sparse
                    else => return array.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // sparse + rest of arrays
                },
            },
            .numeric => unreachable,
        },
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
                .complex => {
                    if (Y == Complex(Integer)) {
                        return Complex(Rational);
                    } else {
                        return Y;
                    }
                },
                else => return Y,
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
            .complex => {
                if (Y == Complex(Integer)) {
                    return Complex(Rational);
                } else {
                    return Y;
                }
            },
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
                if (Y == Complex(Integer)) {
                    return Complex(Real);
                } else if (Y == Complex(Rational)) {
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
            .cfloat => {
                if (X == Complex(Integer)) {
                    return Complex(Rational);
                } else {
                    return X;
                }
            },
            .integer => return X,
            .rational => return X,
            .real => return X,
            .complex => {
                const x: X = .empty;
                const y: Y = .empty;
                return Complex(Coerce(@TypeOf(@field(x, "re")), @TypeOf(@field(y, "re"))));
            },
            .expression => return Y,
        },
        .expression => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return X,
            .integer => return X,
            .rational => return X,
            .real => return X,
            .complex => return X,
            .expression => return X,
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
            },
            .sparse => switch (comptime domainType(Y)) {
                .numeric => {}, // Same as Coerce
                .vector => Coerce(Numeric(X), Numeric(Y)), // sparse vector * vector
                .matrix => switch (comptime matrixType(Y)) {
                    .diagonal => return vector.Sparse(Coerce(Numeric(X), Numeric(Y))), // sparse vector * diagonal matrix
                    else => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // sparse vector * rest of matrices
                },
                .array => @compileError("Cannot coerce vector and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // sparse vector * array
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
                },
                .dense_symmetric => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense symmetric matrix * vector
                    .matrix => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense symmetric matrix * matrix
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense symmetric matrix * array
                },
                .dense_hermitian => switch (comptime domainType(Y)) {
                    .numeric => {}, // Same as Coerce
                    .vector => return vector.Dense(Coerce(Numeric(X), Numeric(Y))), // dense hermitian matrix * vector
                    .matrix => return matrix.general.Dense(Coerce(Numeric(X), Numeric(Y)), orderOf(X)), // dense hermitian matrix * matrix
                    .array => @compileError("Cannot coerce matrix and array types: " ++ @typeName(X) ++ " and " ++ @typeName(Y)), // dense hermitian matrix * array
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
            }
        },
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
                .complex => {
                    if (T2 == Complex(Integer)) {
                        T3 = Complex(Rational);
                    } else {
                        T3 = T2;
                    }
                },
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
            .complex => {
                if (T2 == Complex(Integer)) {
                    T3 = Complex(Rational);
                } else {
                    T3 = T2;
                }
            },
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
                if (T2 == Complex(Integer)) {
                    T3 = Complex(Real);
                } else if (T2 == Complex(Rational)) {
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
                } else if (T1 == Complex(Integer)) {
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
            .expression => T3 = T2,
            else => unreachable,
        },
        .expression => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = T1,
            .integer => T3 = T1,
            .rational => T3 = T1,
            .real => T3 = T1,
            .complex => T3 = T1,
            .expression => T3 = T1,
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
    if (isArray(X)) {
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
            .expression => return true,
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
            .expression => return true,
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
            .expression => return true,
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
            .expression => return true,
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
            .expression => return true,
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
            .expression => return true,
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
            .expression => return true,
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
            .expression => return true,
        },
        .expression => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting expression to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
        },
    }
}

/// Coerces the input type to a floating point type if it is not already a
/// higher range type.
pub fn EnsureFloat(comptime T: type) type {
    if (isArray(T)) {
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
        return vector.Vector(EnsureFloat(Numeric(T)));
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
        else => unreachable,
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

    const numeric = numericType(T);

    switch (comptime numeric) {
        .int => return T,
        .float => return T,
        .bool => return T,
        .cfloat => switch (T) {
            cf16 => return f16,
            cf32 => return f32,
            cf64 => return f64,
            cf80 => return f80,
            cf128 => return f128,
            comptime_complex => return comptime_complex,
            std.math.Complex(f16) => return f16,
            std.math.Complex(f32) => return f32,
            std.math.Complex(f64) => return f64,
            std.math.Complex(f80) => return f80,
            std.math.Complex(f128) => return f128,
            std.math.Complex(comptime_complex) => return comptime_complex,
            else => unreachable,
        },
        .integer => return Integer,
        .rational => return Rational,
        .real => return Real,
        .complex => switch (T) {
            Complex(Integer) => return Integer,
            Complex(Rational) => return Rational,
            Complex(Real) => return Real,
            else => unreachable,
        },
        //.expression => return Expression,
        else => unreachable,
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
    if (isArray(T) or isMatrix(T) or isVector(T)) {
        if (isPermutationMatrix(T)) {
            return T.tp();
        } else {
            const t: T = .empty;
            return @typeInfo(@TypeOf(@field(t, "data"))).pointer.child;
        }
    }

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
/// `complex` or `expression`): The value to cast.
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

    const onumeric = numericType(O);
    const inumeric = numericType(I);

    switch (comptime inumeric) {
        .bool => switch (comptime onumeric) {
            .bool => unreachable,
            .int => return if (value) 1 else 0,
            .float => return if (value) 1 else 0,
            .cfloat => return .{
                .re = if (value) 1 else 0,
                .im = 0,
            },
            else => unreachable,
        },
        .int => switch (comptime onumeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intCast(value),
            .float => return @floatFromInt(value),
            .cfloat => return .{
                .re = @floatFromInt(value),
                .im = 0,
            },
            else => unreachable,
        },
        .float => switch (comptime onumeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intFromFloat(value),
            .float => return @floatCast(value),
            .cfloat => return if (I == Scalar(O)) .{
                .re = value,
                .im = 0,
            } else .{
                .re = @floatCast(value),
                .im = 0,
            },
            else => unreachable,
        },
        .cfloat => switch (comptime onumeric) {
            .bool => return if (value.re != 0 or value.im != 0) true else false,
            .int => return @intFromFloat(value.re),
            .float => return if (Scalar(I) == O) value.re else @floatCast(value.re),
            .cfloat => return .{
                .re = @floatCast(value.re),
                .im = @floatCast(value.im),
            },
            else => unreachable,
        },
        else => unreachable,
    }
}

/// Casts a value of any numeric type to any other numeric type.
///
/// Parameters
/// ----------
/// comptime `T` (`type`): The type to cast to. Must be a supported numeric type.
///
/// `value` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`):The value to cast.
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

    comptime if (isArbitraryPrecision(O)) {
        if (I == O) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .copy = .{ .type = bool, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        }
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (I == O) {
        switch (comptime numericType(O)) {
            .bool, .int, .float, .cfloat => return value,
            .integer => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return integer.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
            .rational => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return rational.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
            .real => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return real.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
            .complex => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return complex.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
            .expression => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return expression.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
        }
        return value;
    }

    const onumeric = numericType(O);
    const inumeric = comptime numericType(I);

    switch (inumeric) {
        .bool => switch (comptime onumeric) {
            .bool => unreachable,
            .int => return if (value) 1 else 0,
            .float => return if (value) 1 else 0,
            .cfloat => return .{
                .re = if (value) 1 else 0,
                .im = 0,
            },
            .integer => @compileError("Not implemented yet: casting from bool to integer"),
            .rational => @compileError("Not implemented yet: casting from bool to rational"),
            .real => @compileError("Not implemented yet: casting from bool to real"),
            .complex => @compileError("Not implemented yet: casting from bool to complex"),
            .expression => @compileError("Not implemented yet: casting from bool to expression"),
        },
        .int => switch (comptime onumeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intCast(value),
            .float => return @floatFromInt(value),
            .cfloat => return .{
                .re = @floatFromInt(value),
                .im = 0,
            },
            .integer => @compileError("Not implemented yet: casting from int to integer"),
            .rational => @compileError("Not implemented yet: casting from int to rational"),
            .real => @compileError("Not implemented yet: casting from int to real"),
            .complex => @compileError("Not implemented yet: casting from int to complex"),
            .expression => @compileError("Not implemented yet: casting from int to expression"),
        },
        .float => switch (comptime onumeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intFromFloat(value),
            .float => return @floatCast(value),
            .cfloat => return if (I == Scalar(O)) .{
                .re = value,
                .im = 0,
            } else .{
                .re = @floatCast(value),
                .im = 0,
            },
            .integer => @compileError("Not implemented yet: casting from float to integer"),
            .rational => @compileError("Not implemented yet: casting from float to rational"),
            .real => @compileError("Not implemented yet: casting from float to real"),
            .complex => @compileError("Not implemented yet: casting from float to complex"),
            .expression => @compileError("Not implemented yet: casting from float to expression"),
        },
        .cfloat => switch (comptime onumeric) {
            .bool => return if (value.re != 0 or value.im != 0) true else false,
            .int => return @intFromFloat(value.re),
            .float => return if (Scalar(I) == O) value.re else @floatCast(value.re),
            .cfloat => return .{
                .re = @floatCast(value.re),
                .im = @floatCast(value.im),
            },
            .integer => @compileError("Not implemented yet: casting from cfloat to integer"),
            .rational => @compileError("Not implemented yet: casting from cfloat to rational"),
            .real => @compileError("Not implemented yet: casting from cfloat to real"),
            .complex => @compileError("Not implemented yet: casting from cfloat to complex"),
            .expression => @compileError("Not implemented yet: casting from cfloat to expression"),
        },
        .integer => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from integer to bool"),
            .int => @compileError("Not implemented yet: casting from integer to int"),
            .float => @compileError("Not implemented yet: casting from integer to float"),
            .cfloat => @compileError("Not implemented yet: casting from integer to cfloat"),
            .integer => unreachable,
            .rational => @compileError("Not implemented yet: casting from integer to rational"),
            .real => @compileError("Not implemented yet: casting from integer to real"),
            .complex => @compileError("Not implemented yet: casting from integer to complex"),
            .expression => @compileError("Not implemented yet: casting from integer to expression"),
        },
        .rational => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from rational to bool"),
            .int => @compileError("Not implemented yet: casting from rational to int"),
            .float => @compileError("Not implemented yet: casting from rational to float"),
            .cfloat => @compileError("Not implemented yet: casting from rational to cfloat"),
            .integer => @compileError("Not implemented yet: casting from rational to integer"),
            .rational => unreachable,
            .real => @compileError("Not implemented yet: casting from rational to real"),
            .complex => @compileError("Not implemented yet: casting from rational to complex"),
            .expression => @compileError("Not implemented yet: casting from rational to expression"),
        },
        .real => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from real to bool"),
            .int => @compileError("Not implemented yet: casting from real to int"),
            .float => @compileError("Not implemented yet: casting from real to float"),
            .cfloat => @compileError("Not implemented yet: casting from real to cfloat"),
            .integer => @compileError("Not implemented yet: casting from real to integer"),
            .rational => @compileError("Not implemented yet: casting from real to rational"),
            .real => unreachable,
            .complex => @compileError("Not implemented yet: casting from real to complex"),
            .expression => @compileError("Not implemented yet: casting from real to expression"),
        },
        .complex => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from complex to bool"),
            .int => @compileError("Not implemented yet: casting from complex to int"),
            .float => @compileError("Not implemented yet: casting from complex to float"),
            .cfloat => @compileError("Not implemented yet: casting from complex to cfloat"),
            .integer => @compileError("Not implemented yet: casting from complex to integer"),
            .rational => @compileError("Not implemented yet: casting from complex to rational"),
            .real => @compileError("Not implemented yet: casting from complex to real"),
            .complex => return value, // Will have to check the type held by the Complex
            .expression => @compileError("Not implemented yet: casting from complex to expression"),
        },
        .expression => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from expression to bool"),
            .int => @compileError("Not implemented yet: casting from expression to int"),
            .float => @compileError("Not implemented yet: casting from expression to float"),
            .cfloat => @compileError("Not implemented yet: casting from expression to cfloat"),
            .integer => @compileError("Not implemented yet: casting from expression to integer"),
            .rational => @compileError("Not implemented yet: casting from expression to rational"),
            .real => @compileError("Not implemented yet: casting from expression to real"),
            .complex => @compileError("Not implemented yet: casting from expression to complex"),
            .expression => unreachable,
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
/// required or not.
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

pub fn getFieldOrDefault(ctx: anytype, comptime field_name: []const u8, comptime FieldType: type, default_value: FieldType) FieldType {
    const T = @TypeOf(ctx);

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

    return default_value;
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

    if (finfo.@"struct".fields.len > info.@"struct".fields.len)
        @compileError("More fields to rename than fields in the struct");

    // Check that all old names exist and new names do not exist in the original struct
    inline for (finfo.@"struct".fields) |field| {
        if (!@hasField(S, field.name))
            @compileError("Field '" ++ field.name ++ "' does not exist in the struct");

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
