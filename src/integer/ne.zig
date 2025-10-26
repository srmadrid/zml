const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Compares an `Integer` with another numeric type for non-equality.
///
/// Signature
/// ---------
/// ```zig
/// fn ne(x: X, y: Y) bool
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// Returns
/// -------
/// `bool`:
/// The comparison result.
pub fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .integer and types.numericType(Y) != .integer)
        @compileError("integer.ne requires at least x or y to be of integer type, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const cmp: types.Cmp = integer.cmp(x, y);
    return cmp != .eq;
}
