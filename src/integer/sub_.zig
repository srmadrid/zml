const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs in-place subtraction between two operands of any numeric type in
/// `Integer` precision. Float, rational or real types are truncated towards
/// zero, and for cfloat or complex types, only the real part is considered.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// Signature
/// ---------
/// ```zig
/// fn sub_(allocator: std.mem.Allocator, o: *Integer, x: X, y: Y, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `allocator` (`std.mem.Allocator`):
/// The allocator to use for memory allocations. Must be the same allocator used
/// to initialize `o`.
///
/// `o` (`*Integer`):
/// A pointer to the output operand where the result will be stored.
///
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// Returns
/// -------
/// `void`
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
///
/// `integer.Error.NotWritable`:
/// If the output operand `o` is not writable.
pub fn sub_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("integer.sub_ requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return integer.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                return integer.add_(allocator, o, x, integer.neg(null, y) catch unreachable);
            },
            .float, .int => {
                var temp: Integer = try .initSet(allocator, y);
                defer temp.deinit(allocator);
                return integer.add_(allocator, o, x, integer.neg(null, temp) catch unreachable);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .float, .int => {
                var tx: Integer = try .initSet(allocator, x);
                defer tx.deinit(allocator);
                var ty: Integer = try .initSet(allocator, y);
                defer ty.deinit(allocator);
                return integer.add_(allocator, o, tx, integer.neg(null, ty) catch unreachable);
            },
            .integer => {
                var temp: Integer = try .initSet(allocator, x);
                defer temp.deinit(allocator);
                return integer.add_(allocator, o, temp, integer.neg(null, y) catch unreachable);
            },
            else => unreachable,
        },

        else => unreachable,
    }
}
