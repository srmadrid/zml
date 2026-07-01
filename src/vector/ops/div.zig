const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isNumeric(Y))
        @compileError("zsl.vector.Div: X must be a vector type and Y must be a numeric type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.div);
}

/// Performs division of a vector by a numeric.
///
/// This function is intended for when the result vector's length is known at
/// compile time, i.e., `x` is a static vector.
///
/// ## Signature
/// ```zig
/// vector.div(x: X, y: Y) vector.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `vector.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
pub fn div(x: anytype, y: anytype) vector.Div(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2Unchecked(x, y, numeric.div);
}

/// Performs division of vector by a numeric, dynamically allocating memory for
/// the result.
///
/// This function is intended for when the result vector's length is known at
/// runtime, i.e., `x` is not a static vector.
///
/// ## Signature
/// ```zig
/// vector.divAlloc(allocator: std.mem.Allocator, x: X, y: Y) !vector.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations.
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `vector.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn divAlloc(allocator: std.mem.Allocator, x: anytype, y: anytype) !vector.Div(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2Alloc(allocator, x, y, numeric.div);
}

/// Performs computation of the division of a vector `x` and a numeric `y` into
/// a vector `o`.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted . Any other form of memory overlap might yield incorrect results.
///
/// For a static output vector, an input static vector and an input numeric, the
/// function cannot return an error.
///
/// ## Signature
/// ```zig
/// vector.divInto(o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output vector operand.
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `vector.Error.DimensionMismatch`: If the two vectors do not have the
///   same length.
pub fn divInto(o: anytype, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isVector(X) or !meta.isNumeric(Y))
        @compileError("zsl.vector.divInto: X must be a vector type and Y must be a numeric type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return vecops.apply2Into(o, x, y, numeric.divInto);
}

/// Performs computation of the division of a vector `x` and a numeric `y` into
/// a vector `o`, without performing any dimension checks.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted . Any other form of memory overlap might yield incorrect results.
///
/// ## Signature
/// ```zig
/// vector.divIntoUnchecked(o: *O, x: X, y: Y) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output vector operand.
/// * `x` (`anytype`): The left vector operand.
/// * `y` (`anytype`): The right numeric operand.
///
/// ## Returns
/// `void`
pub fn divIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isVector(X) or !meta.isNumeric(Y))
        @compileError("zsl.vector.divIntoUnchecked: X must be a vector type and Y must be a numeric type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return vecops.apply2IntoUnchecked(o, x, y, numeric.divInto);
}
