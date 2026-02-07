const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Returns the absolute value of a rational `x`.
///
/// If an allocator is provided, the result will be a new rational allocated
/// with the given allocator. If no allocator is provided, the result will be a
/// view of the rational with the sign set to positive.
///
/// ## Arguments
/// * `x` (`Rational`): The rational value to get the absolute value of.
///
/// ## Returns
/// `Rational`: The absolute value of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if an allocator is provided.
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
