const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs addition between two operands of any numeric type in integer
/// precision. The operation is performed by casting both operands to integer,
/// then adding them.
///
/// ## Signature
/// ```zig
/// integer.add(x: X, y: Y) !Integer
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `Integer`: The result of the addition.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn add(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.integer.add: x and y must be numerics, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    var result: Integer = try .init(allocator, 0);
    errdefer result.deinit(allocator);

    try integer.add_(allocator, &result, x, y);

    return result;
}
