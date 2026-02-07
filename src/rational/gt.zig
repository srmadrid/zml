const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");

/// Compares two operands of rational, integer, dyadic, float, int or bool
/// types, where at least one operand must be of rational type, for greater-than
/// ordering. The operation is performed by casting both operands to rational,
/// then comparing them.
///
/// ## Signature
/// ```zig
/// rational.gt(x: X, y: Y) Cmp
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is greater than `y`, `false` otherwise.
pub fn gt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.rational) or !types.numericType(Y).le(.rational) or
        types.numericType(X) == .cfloat or types.numericType(Y) == .cfloat or
        (types.numericType(X) != .rational and types.numericType(Y) != .rational))
        @compileError("zml.rational.gt: at least one of x or y must be a rational, the other must be a bool, an int, a float, a dyadic, an integer or a rational, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const cmp: types.Cmp = rational.cmp(x, y);
    return cmp == .gt;
}
