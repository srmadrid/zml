const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Compares an `Integer` with another lower or equal precision numeric type for
/// equality or less than.
///
/// Signature
/// ---------
/// ```zig
/// fn le(x: X, y: Y) bool
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
pub fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .integer and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .float) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .int) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .bool) and
        !(types.numericType(X) == .float and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .int and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .bool and types.numericType(Y) == .integer))
        @compileError("integer.le requires x or y to be an integer type, the other must be an integer, float, int or bool type, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const cmp: types.Cmp = integer.cmp(x, y);
    return cmp == .lt or cmp == .eq;
}
