const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if ((!meta.isVector(X) and !meta.isNumeric(X)) or (!meta.isVector(Y) and !meta.isNumeric(Y)) or
        (!meta.isVector(X) and !meta.isVector(Y)) or (meta.isVector(X) and meta.isVector(Y)))
        @compileError("zsl.vector.Mul: at least one of X or Y must be a vector type, the other must be a numeric type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.mul);
}

/// Performs multiplication between a vector and a numeric.
///
/// This function is intended for when the result vector's length is known at
/// compile time, i.e., the vector input is a static vector.
///
/// ## Signature
/// ```zig
/// vector.mul(x: X, y: Y) vector.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left vector or numeric operand.
/// * `y` (`anytype`): The right numeric or vector operand.
///
/// ## Returns
/// `vector.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
pub fn mul(x: anytype, y: anytype) vector.Mul(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2(x, y, numeric.mul);
}

/// Performs multiplication between a vector and a numeric, dynamically
/// allocating memory for the result.
///
/// This function is intended for when the result vector's length is known at
/// runtime, i.e., the vector input is not a static vector.
///
/// ## Signature
/// ```zig
/// vector.mulAlloc(allocator: std.mem.Allocator, x: X, y: Y) !vector.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left vector or numeric operand.
/// * `y` (`anytype`): The right numeric or vector operand.
///
/// ## Returns
/// `vector.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn mulAlloc(allocator: std.mem.Allocator, x: anytype, y: anytype) !vector.Mul(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2Alloc(allocator, x, y, numeric.mul);
}

/// Performs computation of the multiplication of a vectors and a numeric, `x`
/// and `y`, into a vector `o`.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted . Any other form of memory overlap might yield incorrect results.
///
/// For a static output vector, an input static vector and an input numeric, the
/// function cannot return an error.
///
/// ## Signature
/// ```zig
/// vector.mulInto(o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output vector operand.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `vector.Error.DimensionMismatch`: If the two vectors do not have the same
///   length.
pub fn mulInto(o: anytype, x: anytype, y: anytype) !void {
    return vecops.apply2Into(o, x, y, numeric.mulInto);
}

/// Performs computation of the multiplication of a vectors and a numeric, `x`
/// and `y`, into a vector `o`, without performing any dimension checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted . Any other form of memory overlap might yield incorrect results.
///
/// ## Signature
/// ```zig
/// vector.mulIntoUnchecked(o: *O, x: X, y: Y) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output vector operand.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
pub fn mulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    return vecops.apply2IntoUnchecked(o, x, y, numeric.mulInto);
}
