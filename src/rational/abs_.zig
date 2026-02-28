const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

const check_aliasing = @import("check_aliasing.zig").check_aliasing;

/// Performs in-place computation of the absolute value of a rational.
///
/// If the allocator is not provided, `o` and `x` must be the same rational, and
/// the operation will simply set the `positive` field of `o` to `true`.
///
/// ## Arguments
/// * `allocator` (`?std.mem.Allocator`): The allocator to use for memory
///   allocations. If provided, must be the same allocator used to initialize
///   `o`.
/// * `o` (`*Rational`): The output operand.
/// * `x` (`Rational`): The rational value to get the absolute value of.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if an allocator is provided and `o` and `x` are different
///   rationals.
/// * `rational.Error.AllocatorRequired`: If `o` and `x` are different rationals
///   and no allocator is provided.
pub fn abs_(allocator: ?std.mem.Allocator, o: *Rational, x: Rational) !void {
    if (check_aliasing(o, x)) {
        o.num.positive = true;
    } else {
        if (allocator) |a|
            try o.set(a, rational.abs(null, x) catch unreachable)
        else
            return rational.Error.AllocatorRequired;
    }
}
