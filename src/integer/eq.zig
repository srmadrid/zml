const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Compares two operands of integer, dyadic, float, int or bool types, where at
/// least one operand must be of integer type, for equality. The operation is
/// performed by casting both operands to integer, then comparing them.
///
/// ## Signature
/// ```zig
/// integer.eq(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are equal, `false` otherwise.
pub fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.integer) or !types.numericType(Y).le(.integer) or
        types.numericType(X) == .cfloat or types.numericType(Y) == .cfloat or
        (types.numericType(X) != .integer and types.numericType(Y) != .integer))
        @compileError("zml.integer.eq: at least one of x or y must be an integer, the other must be a bool, an int, a float, a dyadic or an integer, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const cmp: types.Cmp = integer.cmp(x, y);
    return cmp == .eq;
}
