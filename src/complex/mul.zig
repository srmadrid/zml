const std = @import("std");

const types = @import("../types.zig");
const complex = @import("../complex.zig");

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.complex) or !types.numericType(Y).le(.complex) or
        (types.numericType(X) != .complex and types.numericType(Y) != .complex))
        @compileError("zml.complex.mul: at least one of x or y must be a complex, the other must be a bool, an int, a float, a dyadic, a cfloat, an integer, a rational, a real or a complex, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

/// Performs multiplication between two operands of complex, real, rational,
/// integer, cfloat, dyadic, float, int or bool types, where at least one
/// operand must be of complex type. The result type is determined by coercing
/// the operand types, and the operation is performed by casting both operands
/// to the result type, then multiplying them.
///
/// ## Signature
/// ```zig
/// complex.mul(allocator: std.mem.Allocator, x: X, y: Y) !complex.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `complex.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn mul(allocator: std.mem.Allocator, x: anytype, y: anytype) !complex.Mul(@TypeOf(x), @TypeOf(y)) {
    const R: type = complex.Mul(@TypeOf(x), @TypeOf(y));

    var result: R = try .init(allocator, 0, 0);
    errdefer result.deinit(allocator);

    try complex.mul_(allocator, &result, x, y);

    return result;
}
