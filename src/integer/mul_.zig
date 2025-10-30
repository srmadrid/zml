const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs in-place multiplication between two operands of any numeric type in
/// `Integer` precision. Float, rational or real types are truncated towards
/// zero, and for cfloat or complex types, only the real part is considered.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// Signature
/// ---------
/// ```zig
/// fn mul_(allocator: std.mem.Allocator, o: *Integer, x: X, y: Y, ctx: anytype) !void
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
pub fn mul_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("integer.mul_ requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return integer.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .expression => @compileError("integer.mul_ not implemented for Expression yet"),
        .complex => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.mul_ not implemented for Complex + Expression yet"),
            .complex => return mul_(allocator, o, x.re, y.re),
            .real => @compileError("integer.mul_ not implemented for Complex + Real yet"),
            .rational => return mul_(allocator, o, x.re, y),
            .integer => return mul_(allocator, o, x.re, y),
            .cfloat => return mul_(allocator, o, x.re, y.re),
            .float => return mul_(allocator, o, x.re, y),
            .int => return mul_(allocator, o, x.re, y),
            .bool => return mul_(allocator, o, x.re, y),
        },
        .real => @compileError("integer.mul_ not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.mul_ not implemented for Rational + Expression yet"),
            .complex => return mul_(allocator, o, x, y.re),
            .real => @compileError("integer.mul_ not implemented for Rational + Real yet"),
            .rational => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return mul_(allocator, o, tx, ty);
            },
            .integer => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                return mul_(allocator, o, tx, y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .float => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, tx, ty[0]);
            },
            .int => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, tx, ty[0]);
            },
            .bool => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                return mul_(allocator, o, tx, types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .integer => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.mul_ not implemented for Integer + Expression yet"),
            .complex => return mul_(allocator, o, x, y.re),
            .real => @compileError("integer.mul_ not implemented for Integer + Real yet"),
            .rational => {
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return mul_(allocator, o, x, ty);
            },
            .integer => {
                if (x.size == 0 or y.size == 0) return o.set(allocator, 0);

                // Aliasing check.
                var tx: Integer = if (o.limbs == x.limbs)
                    try x.copy(allocator)
                else blk: {
                    var tmp: Integer = x;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                defer tx.deinit(allocator);
                var ty: Integer = if (o.limbs == y.limbs)
                    try y.copy(allocator)
                else blk: {
                    var tmp: Integer = y;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                defer ty.deinit(allocator);

                try o.reserve(allocator, tx.size + ty.size);

                var i: u32 = 0;
                while (i < tx.size + ty.size) : (i += 1) {
                    o.limbs[i] = 0;
                }

                i = 0;
                while (i < tx.size) : (i += 1) {
                    var carry: u128 = 0;
                    var j: u32 = 0;
                    while (j < ty.size) : (j += 1) {
                        const acc: u128 = types.scast(u128, o.limbs[i + j]) +
                            types.scast(u128, tx.limbs[i]) * types.scast(u128, ty.limbs[j]) + carry;
                        o.limbs[i + j] = @truncate(acc & 0xFFFFFFFF);
                        carry = acc >> 32;
                    }

                    // Propagate leftover carry.
                    var k: u32 = i + ty.size;
                    while (carry != 0) : (k += 1) {
                        const acc: u128 = types.scast(u128, o.limbs[k]) + carry;
                        o.limbs[k] = @truncate(acc & 0xFFFFFFFF);
                        carry = acc >> 32;
                    }
                }

                o.size = tx.size + ty.size;
                o.positive = tx.positive == ty.positive;
                o.truncate();
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, x, ty[0]);
            },
            .bool => {
                return mul_(allocator, o, x, types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.mul_ not implemented for CFloat + Expression yet"),
            .complex => return mul_(allocator, o, x.re, y.re),
            .real => @compileError("integer.mul_ not implemented for CFloat + Real yet"),
            .rational => return mul_(allocator, o, x.re, y),
            .integer => return mul_(allocator, o, x.re, y),
            .cfloat => return mul_(allocator, o, x.re, y.re),
            .float => return mul_(allocator, o, x.re, y),
            .int => return mul_(allocator, o, x.re, y),
            .bool => return mul_(allocator, o, x.re, y),
        },
        .float => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.mul_ not implemented for Float + Expression yet"),
            .complex => return mul_(allocator, o, x, y.re),
            .real => @compileError("integer.mul_ not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return mul_(allocator, o, tx[0], ty);
            },
            .integer => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return mul_(allocator, o, tx[0], y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .float => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return mul_(allocator, o, tx[0], types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.mul_ not implemented for Int + Expression yet"),
            .complex => return mul_(allocator, o, x, y.re),
            .real => @compileError("integer.mul_ not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return mul_(allocator, o, tx[0], ty);
            },
            .integer => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return mul_(allocator, o, tx[0], y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .float => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return mul_(allocator, o, tx[0], types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer.mul_ not implemented for Bool + Expression yet"),
            .complex => return mul_(allocator, o, x, y.re),
            .real => @compileError("integer.mul_ not implemented for Bool + Real yet"),
            .rational => {
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return mul_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, ty);
            },
            .integer => {
                return mul_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return mul_(allocator, o, types.cast(Integer, x, .{}) catch unreachable, ty[0]);
            },
            .bool => {
                return mul_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
    }
}
