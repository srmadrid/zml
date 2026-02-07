const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const complex = @import("../complex.zig");
const Complex = complex.Complex;

/// Returns the negation of a complex `x`.
///
/// If an allocator is provided, the result will be a new complex allocated
/// with the given allocator. If no allocator is provided, the result will be a
/// view of the complex with the sign flipped.
///
/// ## Signature
/// ```zig
/// complex.neg(allocator: ?std.mem.Allocator, x: X) !X
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The complex value to get the negation of.
///
/// ## Returns
/// `@TypeOf(x)`: The negation of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if an allocator is provided.
pub fn neg(allocator: ?std.mem.Allocator, x: anytype) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!types.isNumeric(X) or types.numericType(X) != .complex)
        @compileError("zml.complex.neg: x must be a complex, got \n\tx: " ++ @typeName(X) ++ "\n");

    if (allocator) |a| {
        var result: X = try x.copy(a);
        ops.neg_(&result.re, result.re, .{}) catch unreachable;
        ops.neg_(&result.im, result.im, .{}) catch unreachable;
        return result;
    }

    var result: X = x;
    result.re = ops.neg(result.re, .{}) catch unreachable;
    result.im = ops.neg(result.im, .{}) catch unreachable;
    result.flags.owns_data = false;
    return result;
}
