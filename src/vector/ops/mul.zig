const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if ((!meta.isVector(X) and !meta.isNumeric(X)) or (!meta.isVector(Y) and !meta.isNumeric(Y)) or
        (!meta.isVector(X) and !meta.isVector(Y)) or (meta.isVector(X) and meta.isVector(Y)))
        @compileError("zsl.vector.mul: at least one of x or y must be a vector, the other must be a numeric, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.mul);
}

/// Performs multiplication between a vector and a numeric.
///
/// ## Signature
/// ```zig
/// vector.mul(allocator: std.mem.Allocator, x: X, y: Y) !vector.Mul(X, Y)
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
pub fn mul(allocator: std.mem.Allocator, x: anytype, y: anytype) !vector.Mul(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2(allocator, x, y, numeric.mul);
}
