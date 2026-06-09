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

/// Performs division of a static vector by a numeric.
///
/// This function is intended for stack-allocated vector types
/// (`vector.Static`).
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
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = vecops.Div(X, Y);

    if (comptime meta.isDenseVector(X) or meta.isSparseVector(X))
        @compileError("zsl.vector.div: the result cannot be a heap-allocated vector type, i.e., x must be a static vector, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tresult: " ++ @typeName(R) ++ "\nFor these inputs use zsl.vector.divAlloc instead.");

    return vecops.apply2(x, y, numeric.div);
}

/// Performs division of vector by a numeric, dynamically allocating memory for
/// the result.
///
/// This function is intended for heap-allocated vector types (`vector.Dense`
/// and `vector.Sparse`), and performs dimension checks at runtime.
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
/// permitted and often more efficient. Any other form of memory overlap might
/// yield incorrect results.
///
/// For three static vectors, or a static output vector, an input static vector
/// and an input numeric, the function cannot return an error.
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
    const O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or
        !meta.isVector(X) or !meta.isNumeric(Y))
        @compileError("zsl.vector.divInto: o must be a mutable one-itme pointer to a vector, x must be a vector, and y must be a numeric, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return vecops.apply2Into(o, x, y, numeric.divInto);
}
