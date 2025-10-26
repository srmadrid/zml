const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Performs division between two operands of any numeric type in `Rational`
/// precision. For cfloat or complex types, only the real part is considered.
///
/// Signature
/// ---------
/// ```zig
/// fn div(allocator: std.mem.Allocator, x: X, y: Y) !Rational
/// ```
///
/// Parameters
/// ----------
/// `allocator` (`std.mem.Allocator`):
/// The allocator to use for memory allocations.
///
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// Returns
/// -------
/// `Rational`:
/// The result of the division.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
///
/// `rational.ZeroDivision`:
/// If `y` is zero.
pub fn div(allocator: std.mem.Allocator, x: anytype, y: anytype) !Rational {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("rational.div requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    var result: Rational = try .init(allocator, 0, 0);
    errdefer result.deinit(allocator);

    try rational.div_(allocator, &result, x, y);

    return result;
}
