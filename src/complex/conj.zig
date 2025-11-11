const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");
const real = @import("../real.zig");
const complex = @import("../complex.zig");
const Complex = complex.Complex;

/// Returns the complex conjugate of a `Complex`. If an allocator is provided,
/// the result will be a new `Complex` allocated with the given allocator. If no
/// allocator is provided, the result will be a view of the `Complex` with the
/// sign of both the imaginary set to its negation.
///
/// Parameters
/// ----------
/// `allocator` (`?std.mem.Allocator`):
/// The allocator to use for memory allocations.
///
/// `x` (`Complex`):
/// The input `Complex`.
///
/// Returns
/// -------
/// `Complex`:
/// The complex conjugate of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only occur if an allocator is provided.
pub fn conj(allocator: ?std.mem.Allocator, x: anytype) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (types.numericType(X) != .complex)
        @compileError("complex.conj requires x to be a complex, got " ++ @typeName(X));

    var result: X = undefined;
    if (allocator) |a| {
        result = try x.copy(a);
    } else {
        result = x;
    }

    if (comptime types.Scalar(X) == rational.Rational) {
        result.im.num.positive = !result.im.num.positive;
    } else { // real
        result.im.rational.num.positive = !result.im.rational.num.positive;
    }

    result.flags.owns_data = false;
    return result;
}
