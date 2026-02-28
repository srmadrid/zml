const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

const check_aliasing = @import("check_aliasing.zig").check_aliasing;

/// Performs in-place computation of the absolute value of an integer.
///
/// If the allocator is not provided, `o` and `x` must be the same integer, and
/// the operation will simply set the `positive` field of `o` to `true`.
///
/// ## Arguments
/// * `allocator` (`?std.mem.Allocator`): The allocator to use for memory
///   allocations. If provided, must be the same allocator used to initialize
///   `o`.
/// * `o` (`*Integer`): The output operand.
/// * `x` (`Integer`): The integer value to get the absolute value of.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if an allocator is provided and `o` and `x` are different
///   integers.
/// * `integer.Error.AllocatorRequired`: If `o` and `x` are different integers
///   and no allocator is provided.
pub fn abs_(allocator: ?std.mem.Allocator, o: *Integer, x: Integer) !void {
    if (check_aliasing(o, x)) {
        o.positive = true;
    } else {
        if (allocator) |a|
            try o.set(a, integer.abs(null, x) catch unreachable)
        else
            return integer.Error.AllocatorRequired;
    }
}
