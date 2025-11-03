const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

const check_aliasing = @import("check_aliasing.zig").check_aliasing;

/// Performs in-place addition between two operands of any numeric type in
/// `Complex` precision.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// Signature
/// ---------
/// ```zig
/// fn add_(allocator: std.mem.Allocator, o: *Complex, x: X, y: Y, ctx: anytype) !void
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
pub fn add_(allocator: std.mem.Allocator, o: anytype, x: anytype, y: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("complex.add_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.numericType(O) != .complex)
        @compileError("complex.add_ requires o to be a complex, got " ++ @typeName(O));

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("complex.add_ requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return complex.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .expression => @compileError("complex.add_ not implemented for Expression yet"),
        .complex => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.add_ not implemented for Complex + Expression yet"),
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
                var tx: X = try check_aliasing(allocator, o, x);
                defer tx.deinit(allocator);
                var ty: Y = try check_aliasing(allocator, o, y);
                defer ty.deinit(allocator);

                // (a + bi) + (c + di) = (a + c) + (b + d)i
                try ops.add_(&o.re, x.re, y.re, .{ .allocator = allocator });
                try ops.add_(&o.im, x.im, y.im, .{ .allocator = allocator });

                return;
            },
            .real => @compileError("complex.add_ not implemented for Complex + Real yet"),
            .rational => return add_(allocator, o, x, y.asComplex()),
            .integer => return add_(allocator, o, x, y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return add_(allocator, o, x, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return add_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return add_(allocator, o, x, ty[0]);
            },
            .bool => return add_(allocator, o, x, types.cast(X, y, .{}) catch unreachable),
        },
        .real => @compileError("complex.add_ not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.add_ not implemented for Rational + Expression yet"),
            .complex => return add_(allocator, o, x.asComplex(), y),
            .real => @compileError("complex.add_ not implemented for Rational + Real yet"),
            .rational => return add_(allocator, o, x.asComplex(), y.asComplex()),
            .integer => return add_(allocator, o, x.asComplex(), y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return add_(allocator, o, x.asComplex(), ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return add_(allocator, o, x.asComplex(), ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return add_(allocator, o, x.asComplex(), ty[0]);
            },
            .bool => return add_(allocator, o, x.asComplex(), types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable),
        },
        .integer => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.add_ not implemented for Integer + Expression yet"),
            .complex => return add_(allocator, o, x.asComplex(), y),
            .real => @compileError("complex.add_ not implemented for Integer + Real yet"),
            .rational => return add_(allocator, o, x.asComplex(), y.asComplex()),
            .integer => return add_(allocator, o, x.asComplex(), y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return add_(allocator, o, x.asComplex(), ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return add_(allocator, o, x.asComplex(), ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return add_(allocator, o, x.asComplex(), ty[0]);
            },
            .bool => {
                return add_(allocator, o, x.asComplex(), types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable);
            },
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.add_ not implemented for CFloat + Expression yet"),
            .complex => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return add_(allocator, o, tx[0], y);
            },
            .real => @compileError("complex.add_ not implemented for CFloat + Real yet"),
            .rational => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return add_(allocator, o, tx[0], y.asComplex());
            },
            .integer => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return add_(allocator, o, tx[0], y.asComplex());
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
                return add_(allocator, o, tx[0], ty[0]);
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
                return add_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return add_(allocator, o, tx[0], types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable);
            },
        },
        .float => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.add_ not implemented for Float + Expression yet"),
            .complex => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return add_(allocator, o, tx[0], y);
            },
            .real => @compileError("complex.add_ not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return add_(allocator, o, tx[0], y.asComplex());
            },
            .integer => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return add_(allocator, o, tx[0], y.asComplex());
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
                return add_(allocator, o, tx[0], ty[0]);
            },
            .float => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return add_(allocator, o, tx[0], types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.add_ not implemented for Int + Expression yet"),
            .complex => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return add_(allocator, o, tx[0], y);
            },
            .real => @compileError("complex.add_ not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return add_(allocator, o, tx[0], y.asComplex());
            },
            .integer => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return add_(allocator, o, tx[0], y.asComplex());
            },
            .cfloat => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .float => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return add_(allocator, o, tx[0], types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.add_ not implemented for Bool + Expression yet"),
            .complex => return add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, y),
            .real => @compileError("complex.add_ not implemented for Bool + Real yet"),
            .rational => return add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, y.asComplex()),
            .integer => return add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, ty[0]);
            },
            .bool => {
                return add_(
                    allocator,
                    o,
                    types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable,
                    types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable,
                );
            },
        },
    }
}
