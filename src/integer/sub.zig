const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs subtraction between two operands of any numeric type in integer
/// precision. The operation is performed by casting both operands to integer,
/// then subtracting them.
///
/// If either `x` or `y` is of custom numeric type, that type must implement the
/// required `copyToInteger` method. The expected signature and behavior of
/// `copyToInteger` are as follows:
/// * `fn copyToInteger(self: *const @This(), allocator: std.mem.Allocator) !Integer`:
///   Initializes and returns a new integer representing the value of the
///   instance.
///
/// ## Signature
/// ```zig
/// integer.sub(x: X, y: Y) !Integer
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

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.integer.sub: x and y must be numerics, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    var result: Integer = try .init(allocator, 0);
    errdefer result.deinit(allocator);

    try integer.sub_(allocator, &result, x, y);

    return result;
}
