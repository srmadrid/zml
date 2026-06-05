const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isVector(Y))
        @compileError("zsl.vector.sub: x and y must be vectors, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.sub);
}

/// Performs subtraction between two vectors.
///
/// For two static vectors, the allocator is not used and can be set to
/// undefined, and the function cannot return an error.
///
/// ## Signature
/// ```zig
/// vector.sub(allocator: std.mem.Allocator, x: X, y: Y) !vector.Sub(X, Y)
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
pub fn sub(allocator: std.mem.Allocator, x: anytype, y: anytype) !vector.Sub(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2(allocator, x, y, numeric.sub);
}
