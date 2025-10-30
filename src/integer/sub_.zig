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
        .expression => @compileError("integer.sub_ not implemented for Expression yet"),
        .complex => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.sub_ not implemented for Complex + Expression yet"),
            .complex => return sub_(allocator, o, x.re, y.re),
            .real => @compileError("integer.sub_ not implemented for Complex + Real yet"),
            .rational => return sub_(allocator, o, x.re, y),
            .integer => return sub_(allocator, o, x.re, y),
            .cfloat => return sub_(allocator, o, x.re, y.re),
            .float => return sub_(allocator, o, x.re, y),
            .int => return sub_(allocator, o, x.re, y),
            .bool => return sub_(allocator, o, x.re, y),
        },
        .real => @compileError("integer.sub_ not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.sub_ not implemented for Rational + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("integer.sub_ not implemented for Rational + Real yet"),
            .rational => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return integer.add_(allocator, o, tx, integer.neg(null, ty) catch unreachable);
            },
            .integer => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                return integer.add_(allocator, o, tx, integer.neg(null, y) catch unreachable);
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, tx, integer.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, tx, integer.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                return integer.add_(allocator, o, tx, integer.neg(null, types.cast(Integer, y, .{}) catch unreachable) catch unreachable);
            },
        },
        .integer => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.sub_ not implemented for Integer + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("integer.sub_ not implemented for Integer + Real yet"),
            .rational => {
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return integer.add_(allocator, o, x, integer.neg(null, ty) catch unreachable);
            },
            .integer => {
                return integer.add_(allocator, o, x, integer.neg(null, y) catch unreachable);
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, x, integer.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, x, integer.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                return integer.add_(allocator, o, x, integer.neg(null, types.cast(Integer, y, .{}) catch unreachable) catch unreachable);
            },
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.sub_ not implemented for CFloat + Expression yet"),
            .complex => return sub_(allocator, o, x.re, y.re),
            .real => @compileError("integer.sub_ not implemented for CFloat + Real yet"),
            .rational => return sub_(allocator, o, x.re, y),
            .integer => return sub_(allocator, o, x.re, y),
            .cfloat => return sub_(allocator, o, x.re, y.re),
            .float => return sub_(allocator, o, x.re, y),
            .int => return sub_(allocator, o, x.re, y),
            .bool => return sub_(allocator, o, x.re, y),
        },
        .float => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.sub_ not implemented for Float + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("integer.sub_ not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return integer.add_(allocator, o, tx[0], integer.neg(null, ty) catch unreachable);
            },
            .integer => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return integer.add_(allocator, o, tx[0], integer.neg(null, y) catch unreachable);
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, tx[0], integer.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, tx[0], integer.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return integer.add_(allocator, o, tx[0], integer.neg(null, types.cast(Integer, y, .{}) catch unreachable) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.sub_ not implemented for Int + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("integer.sub_ not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return integer.add_(allocator, o, tx[0], integer.neg(null, ty) catch unreachable);
            },
            .integer => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return integer.add_(allocator, o, tx[0], integer.neg(null, y) catch unreachable);
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, tx[0], integer.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, tx[0], integer.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return integer.add_(allocator, o, tx[0], integer.neg(null, types.cast(Integer, y, .{}) catch unreachable) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.sub_ not implemented for Bool + Expression yet"),
            .complex => return sub_(allocator, o, x, y.re),
            .real => @compileError("integer.sub_ not implemented for Bool + Real yet"),
            .rational => {
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return integer.add_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, integer.neg(null, ty) catch unreachable);
            },
            .integer => {
                return integer.add_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, integer.neg(null, y) catch unreachable);
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, integer.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return integer.add_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, integer.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                return integer.add_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    integer.neg(null, types.cast(Integer, y, .{}) catch unreachable) catch unreachable,
                );
            },
        },
    }
}
