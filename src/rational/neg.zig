const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Returns the negation of a rational `x`.
///
/// If an allocator is provided, the result will be a new rational allocated
/// with the given allocator. If no allocator is provided, the result will be a
/// view of the rational with the sign flipped.
///
/// ## Arguments
/// * `x` (`Rational`): The rational value to get the negation of.
///
/// ## Returns
/// `Rational`: The negation of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if an allocator is provided.
pub fn neg(allocator: ?std.mem.Allocator, x: Rational) !Rational {
    if (allocator) |a| {
        var result: Rational = try x.copy(a);
        result.num.positive = !result.num.positive;
        return result;
    }

    var result: Rational = x;
    result.num.positive = !result.num.positive;
    result.flags.owns_data = false;
    return result;
}
