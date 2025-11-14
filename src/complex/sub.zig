const std = @import("std");

const types = @import("../types.zig");
const Scalar = types.Scalar;
const Coerce = types.Coerce;
const rational = @import("../rational.zig");
const Rational = rational.Rational;
const complex = @import("../complex.zig");
const Complex = complex.Complex;

/// Performs subtraction between two operands of any numeric type in `Complex`
/// precision.
///
/// Signature
/// ---------
/// ```zig
/// fn sub(allocator: std.mem.Allocator, x: X, y: Y) !Complex
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
/// `Complex`:
/// The result of the subtraction.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
pub fn sub(allocator: std.mem.Allocator, x: anytype, y: anytype) !Complex(Scalar(Coerce(Rational, Coerce(@TypeOf(x), @TypeOf(y))))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Complex(Scalar(Coerce(Rational, Coerce(X, Y))));

    comptime if (types.numericType(X) != .complex and types.numericType(X) != .real and types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .complex and types.numericType(X) != .real and types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.sub_ requires x and y to be an int, float, integer, rational, real or complex, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    var result: C = try .init(allocator, 0, 0);
    errdefer result.deinit(allocator);

    try complex.sub_(allocator, &result, x, y);

    return result;
}
