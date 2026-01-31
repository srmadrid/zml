const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Returns the negation of an integer `x`.
///
/// If an allocator is provided, the result will be a new integer allocated with
/// the given allocator. If no allocator is provided, the result will be a view
/// of the integer with the sign flipped.
///
/// ## Signature
/// ```zig
/// integer.neg(allocator: ?std.mem.Allocator, x: X) X
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The integer value to get the negation of.
///
/// ## Returns
/// `@TypeOf(x)`: The negation of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if an allocator is provided.
pub fn neg(allocator: ?std.mem.Allocator, x: Integer) !Integer {
    if (allocator) |a| {
        var result: Integer = try x.copy(a);
        result.positive = !result.positive;
        return result;
    }

    var result: Integer = x;
    result.positive = !result.positive;
    result.flags.owns_data = false;
    return result;
}
