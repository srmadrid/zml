const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Performs subtraction between two operands of rational, integer, cfloat,
/// dyadic, float, int or bool types, where at least one operand must be of
/// rational type. The operation is performed by casting both operands to
/// rational, then subtracting them.
///
/// ## Signature
/// ```zig
/// rational.sub(allocator: std.mem.Allocator, x: X, y: Y) !Rational
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `Rational`: The result of the subtraction.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn sub(allocator: std.mem.Allocator, x: anytype, y: anytype) !Rational {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.rational) or !types.numericType(Y).le(.rational) or
        (types.numericType(X) != .rational and types.numericType(Y) != .rational))
        @compileError("zml.rational.sub: at least one of x or y must be a rational, the other must be a bool, an int, a float, a dyadic, a cfloat, an integer or a rational, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    var result: Rational = try .init(allocator, 0, 0);
    errdefer result.deinit(allocator);

    try rational.sub_(allocator, &result, x, y);

    return result;
}
