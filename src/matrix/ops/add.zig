const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const matrix = @import("../../matrix.zig");

const matops = @import("../ops.zig");

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isMatrix(X) or !meta.isMatrix(Y))
        @compileError("zsl.matrix.Add: X and Y must be matrix types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return matops.Apply2(X, Y, numeric.add);
}

/// Performs addition between two matrices.
///
/// The result inherits its memory layout from the inputs, i.e., if the input
/// layouts mismatch, the left operand (`x`) strictly dictates the output
/// layout, unless it provides no layout information. For more control over
/// layouts, use `matrix.addInto`.
///
/// This function is intended for when the result matrix's dimensions is known
/// at compile time, i.e., at least one of the inputs is a static matrix. For
/// two static matrices, dimension checks are performed at compile time, for any
/// other combination, dimension checks are performed at runtime throught
/// `std.debug.assert`.
///
/// ## Signature
/// ```zig
/// matrix.add(x: X, y: Y) matrix.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right matrix operand.
///
/// ## Returns
/// `matrix.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
pub fn add(x: anytype, y: anytype) matrix.Add(@TypeOf(x), @TypeOf(y)) {
    return matops.apply2(x, y, numeric.add);
}

/// Performs addition between two matrices, dynamically allocating memory for
/// the result.
///
/// The result inherits its memory layout from the inputs, i.e., if the input
/// layouts mismatch, the left operand (`x`) strictly dictates the output
/// layout, unless it provides no layout information. For more control over
/// layouts, use `matrix.addInto`.
///
/// This function is intended for when the result matrix's dimensions is known
/// at runtime, i.e., none of the inputs is a static matrix.
///
/// ## Signature
/// ```zig
/// matrix.addAlloc(allocator: std.mem.Allocator, x: X, y: Y) !matrix.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right matrix operand.
///
/// ## Returns
/// `matrix.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `matrix.Error.DimensionMismatch`: If the two matrices do not have the same
///   dimensions.
pub fn addAlloc(allocator: std.mem.Allocator, x: anytype, y: anytype) !matrix.Add(@TypeOf(x), @TypeOf(y)) {
    return matops.apply2Alloc(allocator, x, y, numeric.add);
}

/// Performs computation of the addition of two matrices `x` and `y` into a
/// matrix `o`.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted for static or dense outputs, or sparse outputs when both inputs
/// are not matrices. Any other form of memory overlap might yield incorrect
/// results.
///
/// For three static matrices, or a static output matrix, an input static matrix
/// and an input numeric, the function cannot return an error.
///
/// ## Signature
/// ```zig
/// matrix.addInto(o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output matrix operand.
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right matrix operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `matrix.Error.DimensionMismatch`: If the three matrices do not have the
///   same dimension.
pub fn addInto(o: anytype, x: anytype, y: anytype) !void {
    return matops.apply2Into(o, x, y, numeric.addInto);
}
