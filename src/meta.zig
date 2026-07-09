//! Namespace for type utilities.

const std = @import("std");

const matrix = @import("matrix.zig");

pub const standard_integer_types: [10]type = .{
    u8,   u16,
    u32,  u64,
    u128, i8,
    i16,  i32,
    i64,  i128,
};

pub const IterationOrder = enum {
    left_to_right,
    right_to_left,
};

/// `meta.NumericType` is an enum that represents the different numeric types
/// supported by the library. It is used to categorize types based on their
/// properties and capabilities, such as whether they are ints, floats, complex
/// numbers, etc.
///
/// ## Values
/// * `bool`: Represents the boolean type (`bool`).
/// * `int`: Represents integer types:
///   * `usize`, `u8`, `u16`, `u32`, `u64`, `u128`
///   * `isize`, `i8`, `i16`, `i32`, `i64`, `i128`
///   * `uX`, `iX` (where X is any bit size)
///   * `comptime_int`
/// * `float`: Represents floating point types:
///   * `f16`, `f32`, `f64`, `f80`, `f128`
///   * `comptime_float`
/// * `dyadic`: Represents dyadic rational types (`Dyadic(...)`).
/// * `complex`: Represents complex types:
///   * `cf16`, `cf32`, `cf64`, `cf80`, `cf128`
///   * `Complex(Dyadic(...))`
///   * `Complex(<custom>)`
///   * `comptime_complex`
/// * `custom`: Represents custom user-defined numeric types.
pub const NumericType = enum {
    bool,
    int,
    float,
    dyadic,
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
    static,
    dense,
    sparse,
    numeric, // Fallback for numeric types that are not vectors
};

pub const MatrixType = enum {
    general_static,
    general_dense,
    general_sparse,
    symmetric_static,
    symmetric_dense,
    symmetric_sparse,
    hermitian_static,
    hermitian_dense,
    hermitian_sparse,
    triangular_static,
    triangular_dense,
    triangular_sparse,
    diagonal_static,
    diagonal_sparse,
    permutation_static,
    permutation_sparse,
    builder_sparse,
    numeric, // Fallback for numeric types that are not matrices
};

pub const MatrixKind = enum {
    general,
    symmetric,
    hermitian,
    triangular,
    builder,
    diagonal,
    permutation,
    numeric,
};

pub const MatrixStorage = enum {
    static,
    dense,
    sparse,
    numeric,
};

pub const ArrayType = enum {
    dense,
    sparse,
    builder_sparse,
    numeric, // Fallback for numeric types that are not arrays
};

pub const Domain = enum {
    numeric,
    vector,
    matrix,
    array,
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

/// Checks the the input type `N` and returns the corresponding
/// `meta.NumericType`.
///
/// Checks that the input type is a supported numeric type and returns the
/// corresponding `meta.NumericType` enum value. If the type is not supported,
/// it will raise a compile error.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to check.
///
/// ## Returns
/// `meta.NumericType`: The corresponding `meta.NumericType` enum value.
pub inline fn numericType(comptime N: type) NumericType {
    switch (comptime @typeInfo(N)) {
        .bool => return .bool,
        .int, .comptime_int => return .int,
        .float, .comptime_float => return .float,
        .@"struct", .@"union" => {
            if (comptime !@hasDecl(N, "is_numeric") or !N.is_numeric)
                @compileError("zsl.meta.numericType: " ++ @typeName(N) ++ " is not a supported numeric type");

            if (comptime @hasDecl(N, "is_dyadic") and N.is_dyadic)
                return .dyadic;

            if (comptime (@hasDecl(N, "is_complex") and N.is_complex) or
                N == std.math.Complex(f16) or N == std.math.Complex(f32) or N == std.math.Complex(f64) or
                N == std.math.Complex(f80) or N == std.math.Complex(f128) or N == std.math.Complex(comptime_float))
                return .complex;

            return .custom;
        },
        else => @compileError("zsl.meta.numericType: " ++ @typeName(N) ++ " is not a supported numeric type"),
    }
}

/// Determines the vector type of the input type `V`.
///
/// ## Arguments
/// * `V` (`comptime type`): The type to check.
///
/// ## Returns
/// `meta.VectorType`: The corresponding `meta.VectorType` enum value.
pub fn vectorType(comptime V: type) VectorType {
    if (comptime isStaticVector(V))
        return .static;

    if (comptime isDenseVector(V))
        return .dense;

    if (comptime isSparseVector(V))
        return .sparse;

    return .numeric; // Fallback for numeric types that are not vectors
}

/// Determines the matrix type of the input type `M`.
///
/// ## Arguments
/// * `M` (`comptime type`): The type to check.
///
/// ## Returns
/// `meta.MatrixType`: The corresponding `meta.MatrixType` enum value.
pub fn matrixType(comptime M: type) MatrixType {
    if (comptime isGeneralStaticMatrix(M))
        return .general_static;

    if (comptime isSymmetricStaticMatrix(M))
        return .symmetric_static;

    if (comptime isHermitianStaticMatrix(M))
        return .hermitian_static;

    if (comptime isTriangularStaticMatrix(M))
        return .triangular_static;

    if (comptime isDiagonalStaticMatrix(M))
        return .diagonal_static;

    if (comptime isPermutationStaticMatrix(M))
        return .permutation_static;

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

    if (comptime isDiagonalSparseMatrix(M))
        return .diagonal_sparse;

    if (comptime isPermutationSparseMatrix(M))
        return .permutation_sparse;

    if (comptime isBuilderSparseMatrix(M))
        return .builder_sparse;

    return .numeric; // Fallback for numeric types that are not matrices
}

pub fn matrixKind(comptime M: type) MatrixKind {
    if (comptime isGeneralMatrix(M))
        return .general;

    if (comptime isSymmetricMatrix(M))
        return .symmetric;

    if (comptime isHermitianMatrix(M))
        return .hermitian;

    if (comptime isTriangularMatrix(M))
        return .triangular;

    if (comptime isDiagonalMatrix(M))
        return .diagonal;

    if (comptime isPermutationMatrix(M))
        return .permutation;

    return .numeric; // Fallback for numeric types that are not matrices
}

pub fn matrixStorage(comptime M: type) MatrixStorage {
    if (comptime isStaticMatrix(M))
        return .static;

    if (comptime isDenseMatrix(M))
        return .dense;

    if (comptime isSparseMatrix(M))
        return .sparse;

    return .numeric; // Fallback for numeric types that are not matrices
}

/// Determines the array type of the input type `A`.
///
/// ## Arguments
/// * `A` (`comptime type`): The type to check.
///
/// ## Returns
/// `meta.ArrayType`: The corresponding `meta.ArrayType` enum value.
pub fn arrayType(comptime A: type) ArrayType {
    if (comptime isDenseArray(A))
        return .dense;

    if (comptime isSparseArray(A))
        return .sparse;

    if (comptime isBuilderSparseArray(A))
        return .builder_sparse;

    return .numeric; // Fallback for numeric types that are not arrays
}

/// Determines the domain of the input type `T`.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to check.
///
/// ## Returns
/// `meta.Domain`: The corresponding `meta.Domain` enum value.
pub fn domain(comptime T: type) Domain {
    if (comptime isNumeric(T))
        return .numeric;

    if (comptime isVector(T))
        return .vector;

    if (comptime isMatrix(T))
        return .matrix;

    if (comptime isArray(T))
        return .array;

    @compileError("zsl.meta.domain: " ++ @typeName(T) ++ " does not belong to any supported domain");
}

const type_checks = @import("meta/type_checks.zig");
pub const isSupportedType = type_checks.isSupportedType;
pub const isPointer = type_checks.isPointer;
pub const isManyItemPointer = type_checks.isManyItemPointer;
pub const isConstPointer = type_checks.isConstPointer;
pub const isSlice = type_checks.isSlice;
pub const isSimdVector = type_checks.isSimdVector;
pub const isNumeric = type_checks.isNumeric;
pub const isCustomNumeric = type_checks.isCustomNumeric;
pub const isVector = type_checks.isVector;
pub const isStaticVector = type_checks.isStaticVector;
pub const isDenseVector = type_checks.isDenseVector;
pub const isSparseVector = type_checks.isSparseVector;
pub const isMatrix = type_checks.isMatrix;
pub const isSquareMatrix = type_checks.isSquareMatrix;
pub const isGeneralStaticMatrix = type_checks.isGeneralStaticMatrix;
pub const isSymmetricStaticMatrix = type_checks.isSymmetricStaticMatrix;
pub const isHermitianStaticMatrix = type_checks.isHermitianStaticMatrix;
pub const isTriangularStaticMatrix = type_checks.isTriangularStaticMatrix;
pub const isDiagonalStaticMatrix = type_checks.isDiagonalStaticMatrix;
pub const isPermutationStaticMatrix = type_checks.isPermutationStaticMatrix;
pub const isGeneralDenseMatrix = type_checks.isGeneralDenseMatrix;
pub const isSymmetricDenseMatrix = type_checks.isSymmetricDenseMatrix;
pub const isHermitianDenseMatrix = type_checks.isHermitianDenseMatrix;
pub const isTriangularDenseMatrix = type_checks.isTriangularDenseMatrix;
pub const isGeneralSparseMatrix = type_checks.isGeneralSparseMatrix;
pub const isSymmetricSparseMatrix = type_checks.isSymmetricSparseMatrix;
pub const isHermitianSparseMatrix = type_checks.isHermitianSparseMatrix;
pub const isTriangularSparseMatrix = type_checks.isTriangularSparseMatrix;
pub const isDiagonalSparseMatrix = type_checks.isDiagonalSparseMatrix;
pub const isPermutationSparseMatrix = type_checks.isPermutationSparseMatrix;
pub const isBuilderSparseMatrix = type_checks.isBuilderSparseMatrix;
pub const isGeneralMatrix = type_checks.isGeneralMatrix;
pub const isSymmetricMatrix = type_checks.isSymmetricMatrix;
pub const isHermitianMatrix = type_checks.isHermitianMatrix;
pub const isTriangularMatrix = type_checks.isTriangularMatrix;
pub const isDiagonalMatrix = type_checks.isDiagonalMatrix;
pub const isPermutationMatrix = type_checks.isPermutationMatrix;
pub const isBuilderMatrix = type_checks.isBuilderMatrix;
pub const isStaticMatrix = type_checks.isStaticMatrix;
pub const isDenseMatrix = type_checks.isDenseMatrix;
pub const isSparseMatrix = type_checks.isSparseMatrix;
pub const isArray = type_checks.isArray;
pub const isStaticArray = type_checks.isStaticArray;
pub const isDenseArray = type_checks.isDenseArray;
pub const isSparseArray = type_checks.isSparseArray;
pub const isBuilderSparseArray = type_checks.isBuilderSparseArray;
pub const isIntegral = type_checks.isIntegral;
pub const isNonIntegral = type_checks.isNonIntegral;
pub const isReal = type_checks.isReal;
pub const isComplex = type_checks.isComplex;
pub const isSigned = type_checks.isSigned;
pub const isUnsigned = type_checks.isUnsigned;

/// Determines the safe accumulator type for a given numeric type to prevent
/// overflow or catastrophic precision loss during loop summation.
///
/// ## Arguments
/// * `N` (`type`): The base type being accumulated.
///
/// ## Returns
/// `type`: The widened or stable type suitable for loop accumulation.
pub fn Accumulator(comptime N: type) type {
    comptime if (!isNumeric(N))
        @compileError("zsl.meta.Accumulator: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    switch (comptime numericType(N)) {
        .bool => return usize,
        .int => {
            const info = @typeInfo(N).int;
            return if (comptime info.bits < 64)
                (if (comptime info.signedness == .signed) i64 else u64)
            else
                N;
        },
        .float => {
            return if (comptime @typeInfo(N).float.bits <= 16) f32 else N;
        },
        .dyadic => return N.Accumulator,
        .complex => return N.Accumulator,
        .custom => {
            if (comptime !@hasDecl(N, "Accumulator"))
                @compileError("zsl.meta.Accumulator: custom numeric type " ++ @typeName(N) ++ " must have an `Accumulator` declaration");

            return N.Accumulator;
        },
    }
}

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

/// Returns the scalar type of a given numeric type.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to get the scalar type of. Must be a
///   numeric type.
///
/// ## Returns
/// `type`: The scalar type of the input type.
pub fn Scalar(comptime N: type) type {
    comptime if (!isNumeric(N))
        @compileError("zsl.meta.Scalar: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    switch (comptime numericType(N)) {
        .bool => return N,
        .int => return N,
        .float => return N,
        .dyadic => return N,
        .complex => switch (comptime N) {
            std.math.Complex(f16) => return f16,
            std.math.Complex(f32) => return f32,
            std.math.Complex(f64) => return f64,
            std.math.Complex(f80) => return f80,
            std.math.Complex(f128) => return f128,
            std.math.Complex(comptime_float) => return comptime_float,
            else => return N.Scalar,
        },
        .custom => {
            if (comptime !@hasDecl(N, "Scalar"))
                @compileError("zsl.meta.Scalar: custom numeric type " ++ @typeName(N) ++ " must have a `Scalar` declaration");

            return N.Scalar;
        },
    }
}

pub fn Real(comptime N: type) type {
    comptime if (!isNumeric(N))
        @compileError("zsl.meta.Real: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    switch (comptime numericType(N)) {
        .bool => return N,
        .int => return N,
        .float => return N,
        .dyadic => return N,
        .complex => switch (comptime N) {
            std.math.Complex(f16) => return f16,
            std.math.Complex(f32) => return f32,
            std.math.Complex(f64) => return f64,
            std.math.Complex(f80) => return f80,
            std.math.Complex(f128) => return f128,
            std.math.Complex(comptime_float) => return comptime_float,
            else => return N.Real,
        },
        .custom => {
            if (comptime !@hasDecl(N, "Real"))
                @compileError("zsl.meta.Real: custom numeric type " ++ @typeName(N) ++ " must have a `Real` declaration");

            return N.Real;
        },
    }
}

/// Returns the underlying numeric type of a given numeric type, vector, matrix
/// or array.
///
/// ## Arguments
/// * `T` (`comptime type`): The type to get the numeric type of. Must be a
///   numeric type, vector, matrix, or array.
///
/// ## Returns
/// `type`: The underlying numeric type of the input type.
pub fn Numeric(comptime T: type) type {
    if (comptime !isSupportedType(T))
        @compileError("zsl.meta.Numeric: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .numeric => return T,
        .vector => switch (comptime vectorType(T)) {
            else => return T.Numeric,
            .numeric => return T,
        },
        .matrix => switch (comptime matrixType(T)) {
            else => return T.Numeric,
            .numeric => return T,
        },
        .array => switch (comptime arrayType(T)) {
            else => return T.Numeric,
            .numeric => return T,
        },
    }
}

pub fn layoutOf(comptime T: type) ?matrix.Layout {
    if (comptime !isSupportedType(T))
        @compileError("zsl.meta.layoutOf: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .matrix => return T.storage_layout,
        else => return null,
    }
}

pub fn uploOf(comptime T: type) ?matrix.Uplo {
    if (comptime !isSupportedType(T))
        @compileError("zsl.meta.uploOf: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .matrix => switch (comptime matrixType(T)) {
            else => return T.storage_uplo,
            .numeric => unreachable,
        },
        else => @compileError("zsl.meta.uploOf: T must be a matrix type, got " ++ @typeName(T)),
    }
}

pub fn diagOf(comptime T: type) ?matrix.Diag {
    if (comptime !isSupportedType(T))
        @compileError("zsl.meta.diagOf: " ++ @typeName(T) ++ " is not a supported type");

    switch (comptime domain(T)) {
        .matrix => switch (comptime matrixType(T)) {
            else => return T.storage_diag,
            .numeric => unreachable,
        },
        else => @compileError("zsl.meta.diagOf: T must be a matrix type, got " ++ @typeName(T)),
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

/// Returns the return type of a function when called with the given input
/// types. If the return type does not depend on the input types, it is returned
/// directly.
///
/// ## Arguments
/// * `func` (`comptime anytype`): The function to get the return type of.
/// * `input_types` (`comptime []const type`): The types of the inputs to the
///   function, used to determine the return type for functions with inferred
///   return types. If a parameter in the function type is of type `type`, the
///   corresponding input should be passed directly as the parameter type. For
///   instance, for a function `fn (..., T: type, ...) ...`, `input_types`
///   should be `&.{..., T, ...}` instead of `&.{..., type, ...}`.
///
/// ## Returns
/// `type`: The return type of the function when called with the given input
/// types.
pub fn ReturnTypeFromInputs(
    comptime func: anytype,
    comptime input_types: []const type,
) type {
    const func_params = @typeInfo(@TypeOf(func)).@"fn".params;

    // Correct the input types:
    // - If a parameter in the function type is generic, keep the input type as is without any checks, since it can be any type
    // - If a parameter in the function type is of type `type`, keep `type`
    // - If a parameter in the function type is optional, keep non-optional type
    comptime var corrected_input_types: [input_types.len]type = undefined;
    inline for (func_params, 0..) |func_param, i| {
        if (func_param.is_generic) {
            corrected_input_types[i] = input_types[i];
            continue;
        }

        const info_param = @typeInfo(func_param.type.?);
        if (func_param.type.? == type) {
            corrected_input_types[i] = type;
        } else if (info_param == .optional) {
            if (info_param.optional.child == type)
                corrected_input_types[i] = type;

            if (info_param.optional.child != input_types[i])
                @compileError("zsl.meta.ReturnTypeFromInputs: input type " ++ @typeName(input_types[i]) ++ " does not match the non-optional type " ++ @typeName(info_param.optional.child) ++ " of the corresponding parameter in the function type");

            corrected_input_types[i] = input_types[i];
        } else {
            if (func_param.type.? != input_types[i])
                @compileError("zsl.meta.ReturnTypeFromInputs: input type " ++ @typeName(input_types[i]) ++ " does not match the type " ++ @typeName(func_param.type.?) ++ " of the corresponding parameter in the function type");

            corrected_input_types[i] = input_types[i];
        }
    }

    // Generate the inputs to pass to the function based on the corrected input types:
    // - If a parameter in the function type is of type `type`, pass the input type directly
    // - If a parameter in the function type is of type `std.mem.Allocator`, pass the `useless_allocator`
    // - Otherwise, pass an empty value of the corresponding input type
    comptime var inputs: std.meta.Tuple(&corrected_input_types) = undefined;
    inline for (func_params, input_types, 0..) |func_param, input_type, i| {
        inputs[i] = if (comptime func_param.type == type)
            input_type // The type is passed directly
        else if (comptime input_type == std.mem.Allocator)
            useless_allocator
        else
            undefined;
    }

    switch (comptime input_types.len) {
        0 => return @TypeOf(func()),
        1 => return @TypeOf(func(inputs[0])),
        2 => return @TypeOf(func(inputs[0], inputs[1])),
        3 => return @TypeOf(func(inputs[0], inputs[1], inputs[2])),
        4 => return @TypeOf(func(inputs[0], inputs[1], inputs[2], inputs[3])),
        5 => return @TypeOf(func(inputs[0], inputs[1], inputs[2], inputs[3], inputs[4])),
        else => @compileError("zsl.meta.ReturnTypeFromInputs: functions with more than 5 parameters are not supported"),
    }
}

/// Checks if the type `T` has a method with the given name and type. `anytype`
/// parameters are counted as matching any type.
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
///   types. If a parameter in `method_type` is of type `type`, the
///   corresponding type in `input_types` should be the input type for that
///   parameter.
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
        @compileError("zsl.meta.hasMethod: " ++ @typeName(T) ++ " is not a supported type");

    if (comptime !@hasDecl(T, method_name))
        return false;

    // Test that the method has the correct type
    const info_spec = @typeInfo(method_type);
    if (comptime info_spec != .@"fn")
        @compileError("zsl.meta.hasMethod: method_type must be a function type");

    const info_method = @typeInfo(@TypeOf(@field(T, method_name)));
    if (comptime info_method != .@"fn")
        return false;

    const spec_params = info_spec.@"fn".params;
    comptime var spec_return = info_spec.@"fn".return_type.?;
    const method_params = info_method.@"fn".params;
    comptime var method_return = info_method.@"fn".return_type orelse
        ReturnTypeFromInputs(@field(T, method_name), input_types);

    switch (comptime domain(T)) {
        .numeric => {
            if (method_params.len != spec_params.len)
                return false;

            // Check parameter types
            inline for (spec_params, method_params) |spec_param, method_param| {
                if (comptime spec_param.is_generic or method_param.is_generic)
                    continue;

                if (comptime spec_param.type == type)
                    continue; // The parameter type is `type`, so it matches any type

                if (comptime spec_param.type.? != method_param.type.?)
                    return false;
            }

            // Check return type
            const spec_return_info = @typeInfo(spec_return);
            const method_return_info = @typeInfo(method_return);
            if (comptime spec_return_info == .error_union) {
                if (method_return_info != .error_union)
                    return false;

                spec_return = spec_return_info.error_union.payload;
                method_return = method_return_info.error_union.payload;
            }

            if (comptime spec_return != method_return)
                return false;

            return true;
        },
        else => @compileError("zsl.meta.hasMethod: only implemented for numeric types so far"),
    }
}

/// Checks, in order, if any of the types in `types` has a method with the given
/// name and type. `anytype` parameters are counted as matching any type.
///
/// ## Arguments
/// * `types` (`comptime []const type`): The types to check.
/// * `method_name` (`comptime []const u8`): The name of the method to check.
/// * `method_type` (`comptime type`): The expected type of the method,
///   including parameters and return type. Must not include the allocator
///   parameter for allocated types, and the return type must not be an error
///   union.
/// * `input_types` (`comptime []const type`): The types of the inputs to the
///   method, used to determine the return type for methods with inferred return
///   types. If a parameter in `method_type` is of type `type`, the
///   corresponding type in `input_types` should be the input type for that
///   parameter.
///
/// ## Returns
/// `?type`: The first type in `types` that has the method with the correct
/// type, or `null` if no type has the method.
pub fn anyHasMethod(
    comptime types: []const type,
    comptime method_name: []const u8,
    comptime method_type: type,
    comptime input_types: []const type,
) ?type {
    inline for (types) |T| {
        if (hasMethod(T, method_name, method_type, input_types))
            return T;
    }

    return null;
}

pub fn structToArrayOfTypes(comptime Struct: type) [@typeInfo(Struct).@"struct".fields.len]type {
    const info = @typeInfo(Struct);

    if (info != .@"struct")
        @compileError("zsl.meta.structToArrayOfTypes: Struct must be a struct type");

    const fields = info.@"struct".fields;
    comptime var types: [fields.len]type = undefined;

    for (fields, 0..) |field, i| {
        types[i] = field.type;
    }

    return types;
}
