const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const vector = @import("../../vector.zig");

const vecops = @import("../ops.zig");

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isNumeric(Y))
        @compileError("zsl.vector.div: x must be a vector and y must be a numeric, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return vecops.Apply2(X, Y, numeric.div);
}

/// Performs division of vector by a numeric.
///
/// For two static vectors, or a static vector and a numeric, the allocator is
/// not used and can be set to undefined, and the function cannot return an
/// error.
///
/// ## Signature
/// ```zig
/// vector.div(allocator: std.mem.Allocator, x: X, y: Y) !vector.Div(X, Y)
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
pub fn div(allocator: std.mem.Allocator, x: anytype, y: anytype) !vector.Div(@TypeOf(x), @TypeOf(y)) {
    return vecops.apply2(allocator, x, y, numeric.div);
}
