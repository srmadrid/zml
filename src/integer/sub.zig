const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs subtraction between two operands of integer, cfloat, dyadic, float,
/// int or bool types, where at least one operand must be of integer type. The
/// operation is performed by casting both operands to integer, then subtracting
/// them.
///
/// ## Signature
/// ```zig
/// integer.sub(allocator: std.mem.Allocator, x: X, y: Y) !Integer
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `Integer`: The result of the subtraction.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn sub(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.integer) or !types.numericType(Y).le(.integer) or
        (types.numericType(X) != .integer and types.numericType(Y) != .integer))
        @compileError("zml.integer.sub: at least one of x or y must be an integer, the other must be a bool, an int, a float, a dyadic, a cfloat or an integer, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    var result: Integer = try .init(allocator, 0);
    errdefer result.deinit(allocator);

    try integer.sub_(allocator, &result, x, y);

    return result;
}
