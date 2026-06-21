const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isVector(Y))
        @compileError("zsl.vector.Sub: X and Y must be vector types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.sub);
}

/// Performs subtraction between two vectors.
///
/// This function is intended for when the result vector's length is known at
/// compile time, i.e., at least one of the inputs is a static vector. For two
/// static vectors, dimension checks are performed at compile time, for any
/// other combination, dimension checks are performed at runtime throught
/// `std.debug.assert`.
///
/// ## Signature
/// ```zig
/// vector.sub(x: X, y: Y) vector.Sub(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right vector operand.
///
/// ## Returns
/// `vector.Sub(@TypeOf(x), @TypeOf(y))`: The result of the subtraction.
pub fn sub(x: anytype, y: anytype) vector.Sub(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2(x, y, numeric.sub);
}

/// Performs subtraction between two vectors, dynamically allocating memory for
/// the result.
///
/// This function is intended for when the result vector's length is known at
/// runtime, i.e., none of the inputs is a static vector.
///
/// ## Signature
/// ```zig
/// vector.subAlloc(allocator: std.mem.Allocator, x: X, y: Y) !vector.Sub(X, Y)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right vector operand.
///
/// ## Returns
/// `vector.Sub(@TypeOf(x), @TypeOf(y))`: The result of the subtraction.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `vector.Error.DimensionMismatch`: If the two vectors do not have the same
///   length.
pub fn subAlloc(allocator: std.mem.Allocator, x: anytype, y: anytype) !vector.Sub(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2Alloc(allocator, x, y, numeric.sub);
}

/// Performs computation of the subtraction of two vectors `x` and `y` into a
/// vector `o`.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted for static or dense outputs. Any other form of memory overlap
/// might yield incorrect results.
///
/// For three static vectors, or a static output vector, an input static vector
/// and an input numeric, the function cannot return an error.
///
/// ## Signature
/// ```zig
/// vector.subInto(o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output vector operand.
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right vector operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `vector.Error.DimensionMismatch`: If the three vectors do not have the
///   same length.
pub fn subInto(o: anytype, x: anytype, y: anytype) !void {
    return vecops.apply2Into(o, x, y, numeric.subInto);
}

/// Performs computation of the subtraction of two vectors `x` and `y` into a
/// vector `o`, without performing any dimension checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted for static or dense outputs. Any other form of memory overlap
/// might yield incorrect results.
///
/// ## Signature
/// ```zig
/// vector.subIntoUnchecked(o: *O, x: X, y: Y) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output vector operand.
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right vector operand.
///
/// ## Returns
/// `void`
pub fn subIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    return vecops.apply2IntoUnchecked(o, x, y, numeric.subInto);
}
