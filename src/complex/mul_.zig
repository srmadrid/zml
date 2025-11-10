const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

const check_aliasing_alloc = @import("check_aliasing_alloc.zig").check_aliasing_alloc;

/// Performs in-place multiplication between two operands of any numeric type in
/// `Complex` precision.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// Signature
/// ---------
/// ```zig
/// fn mul_(allocator: std.mem.Allocator, o: *Complex, x: X, y: Y, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `allocator` (`std.mem.Allocator`):
/// The allocator to use for memory allocations. Must be the same allocator used
/// to initialize `o`.
///
/// `o` (`*Complex`):
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
pub fn mul_(allocator: std.mem.Allocator, o: anytype, x: anytype, y: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("complex.mul_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.numericType(O) != .complex)
        @compileError("complex.mul_ requires o to be a complex, got " ++ @typeName(O));

    comptime if (types.numericType(X) != .complex and types.numericType(X) != .real and types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .complex and types.numericType(X) != .real and types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.mul_ requires x and y to be an int, float, integer, rational, real or complex, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .expression => @compileError("complex.mul_ not implemented for Expression yet"),
        .complex => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.mul_ not implemented for Complex + Expression yet"),
            .complex => {
                if (x.num.size == 0) {
                    try o.re.set(allocator, y.re);
                    try o.im.set(allocator, y.im);

                    return;
                }

                if (y.num.size == 0) {
                    try o.re.set(allocator, x.re);
                    try o.im.set(allocator, x.im);

                    return;
                }

                // Aliasing checks
                var tx: X = try check_aliasing_alloc(allocator, o, x);
                defer tx.deinit(allocator);
                var ty: Y = try check_aliasing_alloc(allocator, o, y);
                defer ty.deinit(allocator);

                // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
                var ac = try ops.mul(x.re, y.re, .{ .allocator = allocator });
                defer ac.deinit(allocator);
                var bd = try ops.mul(x.im, y.im, .{ .allocator = allocator });
                defer bd.deinit(allocator);

                try ops.sub_(&o.re, ac, bd, .{ .allocator = allocator });

                try ops.mul_(&ac, x.re, y.im, .{ .allocator = allocator });
                try ops.mul_(&bd, x.im, y.re, .{ .allocator = allocator });

                try ops.mul_(&o.im, ac, bd, .{ .allocator = allocator });

                return;
            },
            .real => @compileError("complex.mul_ not implemented for Complex + Real yet"),
            .rational => return mul_(allocator, o, x, y.asComplex()),
            .integer => return mul_(allocator, o, x, y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return mul_(allocator, o, x, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return mul_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return mul_(allocator, o, x, ty[0]);
            },
            .bool => return mul_(allocator, o, x, types.cast(X, y, .{}) catch unreachable),
        },
        .real => @compileError("complex.mul_ not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.mul_ not implemented for Rational + Expression yet"),
            .complex => return mul_(allocator, o, x.asComplex(), y),
            .real => @compileError("complex.mul_ not implemented for Rational + Real yet"),
            .rational => return mul_(allocator, o, x.asComplex(), y.asComplex()),
            .integer => return mul_(allocator, o, x.asComplex(), y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return mul_(allocator, o, x.asComplex(), ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return mul_(allocator, o, x.asComplex(), ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return mul_(allocator, o, x.asComplex(), ty[0]);
            },
            .bool => return mul_(allocator, o, x.asComplex(), types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable),
        },
        .integer => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.mul_ not implemented for Integer + Expression yet"),
            .complex => return mul_(allocator, o, x.asComplex(), y),
            .real => @compileError("complex.mul_ not implemented for Integer + Real yet"),
            .rational => return mul_(allocator, o, x.asComplex(), y.asComplex()),
            .integer => return mul_(allocator, o, x.asComplex(), y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return mul_(allocator, o, x.asComplex(), ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return mul_(allocator, o, x.asComplex(), ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return mul_(allocator, o, x.asComplex(), ty[0]);
            },
            .bool => {
                return mul_(allocator, o, x.asComplex(), types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable);
            },
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.mul_ not implemented for CFloat + Expression yet"),
            .complex => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return mul_(allocator, o, tx[0], y);
            },
            .real => @compileError("complex.mul_ not implemented for CFloat + Real yet"),
            .rational => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return mul_(allocator, o, tx[0], y.asComplex());
            },
            .integer => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return mul_(allocator, o, tx[0], y.asComplex());
            },
            .cfloat => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .float => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return mul_(allocator, o, tx[0], types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable);
            },
        },
        .float => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.mul_ not implemented for Float + Expression yet"),
            .complex => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return mul_(allocator, o, tx[0], y);
            },
            .real => @compileError("complex.mul_ not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return mul_(allocator, o, tx[0], y.asComplex());
            },
            .integer => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return mul_(allocator, o, tx[0], y.asComplex());
            },
            .cfloat => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .float => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return mul_(allocator, o, tx[0], types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.mul_ not implemented for Int + Expression yet"),
            .complex => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return mul_(allocator, o, tx[0], y);
            },
            .real => @compileError("complex.mul_ not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return mul_(allocator, o, tx[0], y.asComplex());
            },
            .integer => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return mul_(allocator, o, tx[0], y.asComplex());
            },
            .cfloat => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .float => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return mul_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return mul_(allocator, o, tx[0], types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.mul_ not implemented for Bool + Expression yet"),
            .complex => return mul_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, y),
            .real => @compileError("complex.mul_ not implemented for Bool + Real yet"),
            .rational => return mul_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, y.asComplex()),
            .integer => return mul_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return mul_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return mul_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return mul_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, ty[0]);
            },
            .bool => {
                return mul_(
                    allocator,
                    o,
                    types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable,
                    types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable,
                );
            },
        },
    }
}
