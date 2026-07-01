const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const matrix = @import("../../matrix.zig");

const linalg = @import("../../linalg.zig");

const matops = @import("../ops.zig");

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if ((!meta.isMatrix(X) and !meta.isNumeric(X)) or (!meta.isMatrix(Y) and !meta.isNumeric(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)))
        @compileError("zsl.matrix.Mul: at least one of X or Y must be a matrix type, the other must be a numeric or matrix type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return if (comptime meta.isMatrix(X) and meta.isMatrix(Y))
        linalg.Matmul(X, Y)
    else
        matops.Apply2(X, Y, numeric.mul);
}

/// Performs multiplication between two matrices, or a matrix and a numeric.
///
/// For two input matrices, the result inherits its memory layout from the
/// inputs, i.e., if the input layouts mismatch, the left operand (`x`) strictly
/// dictates the output layout, unless it provides no layout information. For
/// more control over layouts, use `matrix.mulInto`.
///
/// This function is intended for when the result matrix's dimensions is known
/// at compile time, i.e., at least one of the inputs is a static matrix. For
/// two static matrices, dimension checks are performed at compile time.
///
/// ## Signature
/// ```zig
/// matrix.mul(x: X, y: Y) !matrix.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left matrix or numeric operand.
/// * `y` (`anytype`): The right matrix or numeric operand.
///
/// ## Returns
/// `matrix.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
///
/// ## Errors
/// * `linalg.Error.DimensionMismatch`: If the two matrices do not have
///   compatible dimensions.
pub fn mul(x: anytype, y: anytype) matrix.Mul(@TypeOf(x), @TypeOf(y)) {
    return if (comptime meta.isMatrix(@TypeOf(x)) and meta.isMatrix(@TypeOf(y)))
        linalg.matmul(x, y)
    else
        matops.apply2(x, y, numeric.mul);
}

/// Performs multiplication between two matrices, or a matrix and a numeric,
/// without performing any dimension checks.
///
/// For two input matrices, the result inherits its memory layout from the
/// inputs, i.e., if the input layouts mismatch, the left operand (`x`) strictly
/// dictates the output layout, unless it provides no layout information. For
/// more control over layouts, use `matrix.mulInto`.
///
/// This function is intended for when the result matrix's dimensions is known
/// at compile time, i.e., at least one of the inputs is a static matrix.
///
/// ## Signature
/// ```zig
/// matrix.mulUnchecked(x: X, y: Y) matrix.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left matrix or numeric operand.
/// * `y` (`anytype`): The right matrix or numeric operand.
///
/// ## Returns
/// `matrix.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
pub fn mulUnchecked(x: anytype, y: anytype) matrix.Mul(@TypeOf(x), @TypeOf(y)) {
    return if (comptime meta.isMatrix(@TypeOf(x)) and meta.isMatrix(@TypeOf(y)))
        linalg.matmulUnchecked(x, y)
    else
        matops.apply2Unchecked(x, y, numeric.mul);
}

/// Performs multiplication between two matrices, or a matrix and a numeric,
/// dynamically allocating memory for the result.
///
/// For two input matrices, the result inherits its memory layout from the
/// inputs, i.e., if the input layouts mismatch, the left operand (`x`) strictly
/// dictates the output layout, unless it provides no layout information. For
/// more control over layouts, use `matrix.mulInto`.
///
/// This function is intended for when the result matrix's dimensions is known
/// at runtime, i.e., none of the inputs is a static matrix.
///
/// ## Signature
/// ```zig
/// matrix.mulAlloc(allocator: std.mem.Allocator, x: X, y: Y) !matrix.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left matrix or numeric operand.
/// * `y` (`anytype`): The right matrix or numeric operand.
///
/// ## Returns
/// `matrix.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `linalg.Error.DimensionMismatch`: If the two matrices do not have
///   compatible dimensions.
pub fn mulAlloc(allocator: std.mem.Allocator, x: anytype, y: anytype) !matrix.Mul(@TypeOf(x), @TypeOf(y)) {
    return if (comptime meta.isMatrix(@TypeOf(x)) and meta.isMatrix(@TypeOf(y)))
        linalg.matmulAlloc(allocator, x, y)
    else
        matops.apply2Alloc(allocator, x, y, numeric.mul);
}

/// Performs computation of the multiplication two matrices, or a matrix and a
/// numeric, `x` and `y`, into a matrix `o`.
///
/// Exact aliasing (in-place modification) between the output and an input, when
/// only one input is a matrix, is permitted. Any other form of memory overlap
/// might yield incorrect results.
///
/// For a static output matrix, two input static matrices, or an input matrix
/// and an input numeric, the function cannot return an error.
///
/// ## Signature
/// ```zig
/// matrix.mulInto(o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output matrix operand.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `matrix.Error.DimensionMismatch`: If the two matrices do not have the same
///   dimensions.
/// * `linalg.Error.DimensionMismatch`: If the three matrices do not have
///   compatible dimensions.
pub fn mulInto(o: anytype, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if ((!meta.isMatrix(X) and !meta.isNumeric(X)) or (!meta.isMatrix(Y) and !meta.isNumeric(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)))
        @compileError("zsl.matrix.mulInto: at least one of x or y must be a matrix, the other must be a numeric or a matrix, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime meta.isMatrix(@TypeOf(x)) and meta.isMatrix(@TypeOf(y)))
        linalg.matmulInto(o, x, y)
    else
        matops.apply2Into(o, x, y, numeric.mulInto);
}

/// Performs computation of the multiplication two matrices, or a matrix and a
/// numeric, `x` and `y`, into a matrix `o`, without performing any dimension
/// checks.
///
/// Exact aliasing (in-place modification) between the output and an input, when
/// only one input is a matrix, is permitted. Any other form of memory overlap
/// might yield incorrect results.
///
/// ## Signature
/// ```zig
/// matrix.mulIntoUnchecked(o: *O, x: X, y: Y) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output matrix operand.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
pub fn mulIntoUnchecked(o: anytype, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if ((!meta.isMatrix(X) and !meta.isNumeric(X)) or (!meta.isMatrix(Y) and !meta.isNumeric(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)))
        @compileError("zsl.matrix.mulIntoUnchecked: at least one of x or y must be a matrix, the other must be a numeric or a matrix, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime meta.isMatrix(@TypeOf(x)) and meta.isMatrix(@TypeOf(y)))
        linalg.matmulIntoUnchecked(o, x, y)
    else
        matops.apply2IntoUnchecked(o, x, y, numeric.mulInto);
}
