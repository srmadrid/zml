const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Returns the absolute value of an integer `x`.
///
/// If an allocator is provided, the result will be a new integer allocated with
/// the given allocator. If no allocator is provided, the result will be a view
/// of the integer with the sign set to positive.
///
/// ## Arguments
/// * `allocator` (`?std.mem.Allocator`): The allocator to use for the result.
/// * `x` (`Integer`): The integer value to get the absolute value of.
///
/// ## Returns
/// `Integer`: The absolute value of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if an allocator is provided.
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
