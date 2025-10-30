const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Performs in-place subtraction between two operands of any numeric type in
/// `Rational` precision. For cfloat or complex types, only the real part is
/// considered.
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
        .expression => @compileError("rational.sub_ not implemented for Expression yet"),
        .complex => switch (comptime types.numericType(Y)) {
            .expression => @compileError("rational.sub_ not implemented for Complex + Expression yet"),
            .complex => return sub_(allocator, o, x.re, y.re),
            .real => @compileError("rational.sub_ not implemented for Complex + Real yet"),
            .rational => return sub_(allocator, o, x.re, y),
            .integer => return sub_(allocator, o, x.re, y),
            .cfloat => return sub_(allocator, o, x.re, y.re),
            .float => return sub_(allocator, o, x.re, y),
            .int => return sub_(allocator, o, x.re, y),
            .bool => return sub_(allocator, o, x.re, y),
        },
        .real => @compileError("rational.sub_ not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .expression => @compileError("rational.sub_ not implemented for Rational + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("rational.sub_ not implemented for Rational + Real yet"),
            .rational => return rational.add_(allocator, o, x, rational.neg(null, y) catch unreachable),
            .integer => return rational.add_(allocator, o, x, rational.neg(null, y.asRational()) catch unreachable),
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return rational.add_(allocator, o, x, rational.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return rational.add_(allocator, o, x, rational.neg(null, ty[0]) catch unreachable);
            },
            .bool => return rational.add_(allocator, o, x, rational.neg(null, types.cast(Rational, y, .{}) catch unreachable) catch unreachable),
        },
        .integer => switch (comptime types.numericType(Y)) {
            .expression => @compileError("rational.sub_ not implemented for Integer + Expression yet"),
            .complex => sub_(allocator, o, x, y.re),
            .real => @compileError("rational.sub_ not implemented for Integer + Real yet"),
            .rational => return rational.add_(allocator, o, x.asRational(), rational.neg(null, y) catch unreachable),
            .integer => return rational.add_(allocator, o, x.asRational(), rational.neg(null, y.asRational()) catch unreachable),
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return rational.add_(allocator, o, x.asRational(), rational.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return rational.add_(allocator, o, x.asRational(), rational.neg(null, ty[0]) catch unreachable);
            },
            .bool => return rational.add_(allocator, o, x.asRational(), rational.neg(null, types.cast(Rational, y, .{}) catch unreachable) catch unreachable),
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .expression => @compileError("rational.sub_ not implemented for CFloat + Expression yet"),
            .complex => return sub_(allocator, o, x.re, y.re),
            .real => @compileError("rational.sub_ not implemented for CFloat + Real yet"),
            .rational => return sub_(allocator, o, x.re, y),
            .integer => return sub_(allocator, o, x.re, y),
            .cfloat => return sub_(allocator, o, x.re, y.re),
            .float => return sub_(allocator, o, x.re, y),
            .int => return sub_(allocator, o, x.re, y),
            .bool => return sub_(allocator, o, x.re, y),
        },
        .float => switch (comptime types.numericType(Y)) {
            .expression => @compileError("rational.sub_ not implemented for Float + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("rational.sub_ not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, y) catch unreachable);
            },
            .integer => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, y.asRational()) catch unreachable);
            },
            .cfloat => return rational.add_(allocator, o, x, y.re),
            .float => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, types.cast(Rational, y, .{}) catch unreachable) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .expression => @compileError("rational.sub_ not implemented for Int + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("rational.sub_ not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, y) catch unreachable);
            },
            .integer => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, y.asRational()) catch unreachable);
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                return rational.add_(allocator, o, tx[0], rational.neg(null, types.cast(Rational, y, .{}) catch unreachable) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .expression => @compileError("rational.sub_ not implemented for Bool + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("rational.sub_ not implemented for Bool + Real yet"),
            .rational => return rational.add_(allocator, o, types.cast(Rational, x, .{}) catch unreachable, rational.neg(null, y) catch unreachable),
            .integer => return rational.add_(allocator, o, types.cast(Rational, x, .{}) catch unreachable, rational.neg(null, y.asRational()) catch unreachable),
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return rational.add_(allocator, o, types.cast(Rational, x, .{}) catch unreachable, rational.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return rational.add_(allocator, o, types.cast(Rational, x, .{}) catch unreachable, rational.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                return rational.add_(
                    allocator,
                    o,
                    types.cast(Rational, x, .{}) catch unreachable,
                    rational.neg(null, types.cast(Rational, y, .{}) catch unreachable) catch unreachable,
                );
            },
        },
    }
}
