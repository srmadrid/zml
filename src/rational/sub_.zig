const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Performs in-place subtraction between two operands of any numeric type in
/// `Rational` precision. For complex types, only the real part is considered.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// Signature
/// ---------
/// ```zig
/// fn sub_(allocator: std.mem.Allocator, o: *Rational, x: X, y: Y, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `allocator` (`std.mem.Allocator`):
/// The allocator to use for memory allocations. Must be the same allocator used
/// to initialize `o`.
///
/// `o` (`*Rational`):
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
/// `rational.Error.NotWritable`:
/// If the output operand `o` is not writable, or if its numerator or
/// denominator are not writable when they need to be modified.
pub fn sub_(allocator: std.mem.Allocator, o: *Rational, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("rational.sub_ requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .rational => switch (comptime types.numericType(Y)) {
            .rational => {
                return rational.add_(allocator, o, x, rational.neg(null, y) catch unreachable);
            },
            .integer => {
                return rational.add_(allocator, o, x, (integer.neg(null, y) catch unreachable).asRational());
            },
            .float, .int => {
                var temp: Rational = try .initSet(allocator, y, 1);
                defer temp.deinit(allocator);
                return rational.add_(allocator, o, x, rational.neg(null, temp) catch unreachable);
            },
            else => unreachable,
        },
        .integer => switch (comptime types.numericType(Y)) {
            .rational => {
                return rational.add_(allocator, o, x.asRational(), rational.neg(null, y) catch unreachable);
            },
            .integer => {
                return rational.add_(allocator, o, x.asRational(), (integer.neg(null, y) catch unreachable).asRational());
            },
            .float, .int => {
                var temp: Rational = try .initSet(allocator, y, 1);
                defer temp.deinit(allocator);
                return rational.add_(allocator, o, x.asRational(), rational.neg(null, temp) catch unreachable);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .rational => {
                var temp: Rational = try .initSet(allocator, x, 1);
                defer temp.deinit(allocator);
                return rational.add_(allocator, o, temp, rational.neg(null, y) catch unreachable);
            },
            .integer => {
                var temp: Rational = try .initSet(allocator, x, 1);
                defer temp.deinit(allocator);
                return rational.add_(allocator, o, temp, (integer.neg(null, y) catch unreachable).asRational());
            },
            .float, .int => {
                var tx: Rational = try .initSet(allocator, x, 1);
                defer tx.deinit(allocator);
                var ty: Rational = try .initSet(allocator, y, 1);
                defer ty.deinit(allocator);
                return rational.add_(allocator, o, tx, rational.neg(null, ty) catch unreachable);
            },
            else => unreachable,
        },

        else => unreachable,
    }
}
