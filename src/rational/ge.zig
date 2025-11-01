const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");

/// Compares a `Rational` with another lower or equal precision numeric type for
/// equality or greater than.
///
/// Signature
/// ---------
/// ```zig
/// fn ge(x: X, y: Y) bool
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
pub fn ge(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .rational and types.numericType(Y) == .rational) and
        !(types.numericType(X) == .rational and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .rational and types.numericType(Y) == .float) and
        !(types.numericType(X) == .rational and types.numericType(Y) == .int) and
        !(types.numericType(X) == .rational and types.numericType(Y) == .bool) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .rational) and
        !(types.numericType(X) == .float and types.numericType(Y) == .rational) and
        !(types.numericType(X) == .int and types.numericType(Y) == .rational) and
        !(types.numericType(X) == .bool and types.numericType(Y) == .rational))
        @compileError("rational.ge requires x or y to be a rational type, the other must be a rational, integer, float, int or bool type, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const cmp: types.Cmp = rational.cmp(x, y);
    return cmp == .gt or cmp == .eq;
}
