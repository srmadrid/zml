const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isVector(Y))
        @compileError("zsl.vector.Sub: X and T must be vector types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.sub);
}

/// Performs subtraction between two static vectors.
///
/// This function is intended for stack-allocated vector types
/// (`vector.Static`), and performs dimension checks at compile time.
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
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = vecops.Sub(X, Y);

    if (comptime meta.isDenseVector(X) or meta.isSparseVector(X) or
        meta.isDenseVector(Y) or meta.isSparseVector(Y))
        @compileError("zsl.vector.sub: the result cannot be a heap-allocated vector type, i.e., both inputs must be static vectors, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tresult: " ++ @typeName(R) ++ "\n");

    return vecops.apply2(x, y, numeric.sub);
}

/// Performs subtraction between two vectors, dynamically allocating memory for
/// the result.
///
/// This function is intended for heap-allocated vector types (`vector.Dense`
/// and `vector.Sparse`), and performs dimension checks at runtime.
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
/// permitted and often more efficient. Any other form of memory overlap might
/// yield incorrect results.
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
    const O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or
        !meta.isVector(X) or !meta.isVector(Y))
        @compileError("zsl.vector.subInto: o must be a mutable one-itme pointer to a vector, and x and y must be vectors, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return vecops.apply2Into(o, x, y, numeric.subInto);
}
