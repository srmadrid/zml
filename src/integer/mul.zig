const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs multiplication between two operands of any numeric type in
/// `Integer` precision. Float, rational or real types are truncated towards
/// zero, and for cfloat or complex types, only the real part is considered.
///
/// Signature
/// ---------
/// ```zig
/// fn mul(allocator: std.mem.Allocator, x: X, y: Y) !Rational
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
/// `Integer`:
/// The result of the multiplication.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
pub fn mul(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("integer.mul requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    var result: Integer = try .init(allocator, 0);
    errdefer result.deinit(allocator);

    try integer.mul_(allocator, &result, x, y);

    return result;
}
