const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const matrix = @import("../../matrix.zig");

const matops = @import("../ops.zig");

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isMatrix(X) or !meta.isNumeric(Y))
        @compileError("zsl.matrix.Div: X must be a matrix type and Y must be a numeric type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return matops.Apply2(X, Y, numeric.div);
}

/// Performs division of a matrix by a numeric.
///
/// This function is intended for when the result matrix's dimensions is known
/// at compile time, i.e., `x` is a static matrix.
///
/// ## Signature
/// ```zig
/// matrix.div(x: X, y: Y) matrix.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `matrix.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
pub fn div(x: anytype, y: anytype) matrix.Div(@TypeOf(x), @TypeOf(y)) {
    return matops.apply2(x, y, numeric.div);
}

/// Performs division of matrix by a numeric, dynamically allocating memory for
/// the result.
///
/// This function is intended for when the result matrix's dimensions is known
/// at runtime, i.e., `x` is not a static matrix.
///
/// ## Signature
/// ```zig
/// matrix.divAlloc(allocator: std.mem.Allocator, x: X, y: Y) !matrix.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `matrix.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn divAlloc(allocator: std.mem.Allocator, x: anytype, y: anytype) !matrix.Div(@TypeOf(x), @TypeOf(y)) {
    return matops.apply2Alloc(allocator, x, y, numeric.div);
}

/// Performs computation of the division of a matrix `x` and a numeric `y` into
/// a matrix `o`.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted. Any other form of memory overlap might yield incorrect results.
///
/// For a static output matrix, an input static matrix and an input numeric, the
/// function cannot return an error.
///
/// ## Signature
/// ```zig
/// matrix.divInto(o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output matrix operand.
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `matrix.Error.DimensionMismatch`: If the two matrices do not have the
///   same dimensions.
pub fn divInto(o: anytype, x: anytype, y: anytype) !void {
    return matops.apply2Into(o, x, y, numeric.divInto);
}

/// Performs computation of the division of a matrix `x` and a numeric `y` into
/// a matrix `o`, without performing any dimension checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted. Any other form of memory overlap might yield incorrect results.
///
/// ## Signature
/// ```zig
/// matrix.divIntoUnchecked(o: *O, x: X, y: Y) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output matrix operand.
/// * `x` (`anytype`): The left matrix operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `void`
pub fn divIntoUnchecked(o: anytype, x: anytype, y: anytype) !void {
    return matops.apply2IntoUnchecked(o, x, y, numeric.divInto);
}
