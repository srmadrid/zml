const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Performs subtraction between two operands of any numeric type in `Rational`
/// precision. For cfloat or complex types, only the real part is considered.
///
/// Signature
/// ---------
/// ```zig
/// fn sub(allocator: std.mem.Allocator, x: X, y: Y) !Rational
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
/// The result of the subtraction.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
pub fn sub(allocator: std.mem.Allocator, x: anytype, y: anytype) !Rational {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("rational.sub requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    var result: Rational = try .init(allocator, 0, 0);
    errdefer result.deinit(allocator);

    try rational.sub_(allocator, &result, x, y);

    return result;
}
