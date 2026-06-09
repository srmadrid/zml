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

/// Performs multiplication between a static vector and a numeric.
///
/// This function is intended for stack-allocated vector types
/// (`vector.Static`).
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
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = vecops.Mul(X, Y);

    if (comptime meta.isDenseVector(X) or meta.isSparseVector(X) or
        meta.isDenseVector(Y) or meta.isSparseVector(Y))
        @compileError("zsl.vector.mul: the result cannot be a heap-allocated vector type, i.e., the vector input must be a static vector, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tresult: " ++ @typeName(R) ++ "\nFor these inputs use zsl.vector.mulAlloc instead.");

    return vecops.apply2(x, y, numeric.mul) catch unreachable;
}

/// Performs multiplication between a vector and a numeric, dynamically
/// allocating memory for the result.
///
/// This function is intended for heap-allocated vector types (`vector.Dense`
/// and `vector.Sparse`), and performs dimension checks at runtime.
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
/// permitted and often more efficient. Any other form of memory overlap might
/// yield incorrect results.
///
/// For three static vectors, or a static output vector, an input static vector
/// and an input numeric, the function cannot return an error.
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
    const O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or
        (!meta.isVector(X) and !meta.isNumeric(X)) or (!meta.isVector(Y) and !meta.isNumeric(Y)) or
        (!meta.isVector(X) and !meta.isVector(Y)))
        @compileError("zsl.vector.mulInto: o must be a mutable one-itme pointer to a vector, and at least one of x or y must be a vector, the other must be a vector or a numeric, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return vecops.apply2Into(o, x, y, numeric.mulInto);
}
