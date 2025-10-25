const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Returns the absolute value of a `Rational`. If an allocator is provided,
/// the result will be a new `Rational` allocated with the given allocator. If
/// no allocator is provided, the result will be a view of the `Rational` with
/// the sign set to positive.
///
/// Parameters
/// ----------
/// `allocator` (`?std.mem.Allocator`):
/// The allocator to use for memory allocations. If provided it must be the same
/// allocator used to initialize `o`.
///
/// `x` (Rational):
/// The input `Rational`.
///
/// Returns
/// -------
/// `Rational`:
/// The absolute value of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only occur if an allocator is provided.
pub fn abs(allocator: ?std.mem.Allocator, x: Rational) !Rational {
    if (allocator) |a| {
        var result: Rational = try x.copy(a);
        result.num.positive = true;
        return result;
    }

    var result: Rational = x;
    result.num.positive = true;
    result.flags.owns_data = false;
    return result;
}
