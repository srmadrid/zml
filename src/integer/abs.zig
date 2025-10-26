const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Returns the absolute value of an `Integer`. If an allocator is provided,
/// the result will be a new `Integer` allocated with the given allocator. If
/// no allocator is provided, the result will be a view of the `Integer` with
/// the sign set to positive.
///
/// Parameters
/// ----------
/// `allocator` (`?std.mem.Allocator`):
/// The allocator to use for memory allocations.
///
/// `x` (`Integer`):
/// The input `Integer`.
///
/// Returns
/// -------
/// `Integer`:
/// The absolute value of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only occur if an allocator is provided.
pub fn abs(allocator: ?std.mem.Allocator, x: Integer) !Integer {
    if (allocator) |a| {
        var result: Integer = try x.copy(a);
        result.positive = true;
        return result;
    }

    var result: Integer = x;
    result.positive = true;
    result.flags.owns_data = false;
    return result;
}
