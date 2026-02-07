const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

const check_aliasing_alloc = @import("check_aliasing_alloc.zig").check_aliasing_alloc;

/// Performs in-place division between two operands of complex, real, rational,
/// integer, cfloat, dyadic, float, int or bool types, where at least one
/// operand must be of complex type. The operation is performed by casting both
/// operands to complex, then diving them in-place.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// ## Signature
/// ```zig
/// complex.div_(allocator: std.mem.Allocator, o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations. Must be the same allocator used to initialize `o`.
/// * `o` (`anytype`): A pointer to the output operand where the result will be
///   stored.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `complex.Error.ZeroDivision`: If `y` is zero.
/// * `complex.Error.NotWritable`: If the output operand `o` is not writable.
/// * `complex.Error.DataNotOwned`: If the output operand `o` does not own its
///   data and resizing is needed.
/// * `real.Error.NotWritable`: If the output operand `o` holds reals and its
///   real or imaginary part is not writable when it needs to be modified.
/// * `real.Error.DataNotOwned`: If the output operand `o` holds reals and its
///   real or imaginary part does not own its data and resizing is needed.
/// * `rational.Error.NotWritable`: If the output operand `o` holds rationals
///   and its real or imaginary part is not writable when it needs to be
///   modified.
/// * `rational.Error.DataNotOwned`: If the output operand `o` holds rationals
///   and its real or imaginary part does not own its data and resizing is
///   needed.
pub fn div_(allocator: std.mem.Allocator, o: anytype, x: anytype, y: anytype) !void {
    const O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O) or
        !types.isNumeric(types.Child(O)) or types.numericType(types.Child(O)) != .complex or
        !types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.complex) or !types.numericType(Y).le(.complex) or
        (types.numericType(X) != .complex and types.numericType(Y) != .complex))
        @compileError("zml.complex.div_: o must be a mutable pointer to a complex, and at least one of x or y must be a complex, the other must be a bool, an int, a float, a dyadic, a cfloat, an integer, a rational, a real or a complex, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (!o.flags.writable)
        return complex.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .complex => switch (comptime types.numericType(Y)) {
            .complex => {
                if (ops.eq(x, 0, .{}) catch unreachable) {
                    try o.re.set(allocator, x.re);
                    try o.im.set(allocator, x.im);

                    return;
                }

                if (ops.eq(y, 0, .{}) catch unreachable)
                    return complex.Error.ZeroDivision;

                // Aliasing checks
                var tx: X = try check_aliasing_alloc(allocator, o, x);
                defer tx.deinit(allocator);
                var ty: Y = try check_aliasing_alloc(allocator, o, y);
                defer ty.deinit(allocator);

                // (a + bi) / (c + di) = ((ac + bd) / (c² + d²)) + ((bc - ad) / (c² + d²))i
                var ac = try ops.div(x.re, y.re, .{ .allocator = allocator });
                defer ac.deinit(allocator);
                var bd = try ops.div(x.im, y.im, .{ .allocator = allocator });
                defer bd.deinit(allocator);

                var c2 = try ops.mul(y.re, y.re, .{ .allocator = allocator });
                defer c2.deinit(allocator);
                var d2 = try ops.mul(y.im, y.im, .{ .allocator = allocator });
                defer d2.deinit(allocator);

                try ops.div_(&c2, c2, d2, .{ .allocator = allocator });

                try ops.div_(&ac, ac, bd, .{ .allocator = allocator });
                try ops.div_(&o.re, ac, c2, .{ .allocator = allocator });

                try ops.mul_(&ac, x.im, y.re, .{ .allocator = allocator });
                try ops.mul_(&bd, x.re, y.im, .{ .allocator = allocator });

                try ops.sub_(&ac, ac, bd, .{ .allocator = allocator });
                try ops.div_(&o.im, ac, c2, .{ .allocator = allocator });

                return;
            },
            .real => return div_(allocator, o, x, y.asComplex()),
            .rational => return div_(allocator, o, x, y.asComplex()),
            .integer => return div_(allocator, o, x, y.asComplex()),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];

                return div_(allocator, o, x, ty[0]);
            },
            .dyadic => {
                var ty = try @import("../dyadic/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];

                return div_(allocator, o, x, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];

                return div_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];

                return div_(allocator, o, x, ty[0]);
            },
            .bool => return div_(
                allocator,
                o,
                x,
                types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable,
            ),
            else => unreachable,
        },
        .real => return div_(allocator, o, x.asComplex(), y),
        .rational => return div_(allocator, o, x.asComplex(), y),
        .integer => return div_(allocator, o, x.asComplex(), y),
        .cfloat => {
            var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
            tx[0].re.num.limbs = &tx[1][0];
            tx[0].re.den.limbs = &tx[1][1];
            tx[0].im.num.limbs = &tx[1][2];
            tx[0].im.den.limbs = &tx[1][3];

            return div_(allocator, o, tx[0], y);
        },
        .dyadic => {
            var tx = try @import("../dyadic/asComplex.zig").asComplex(x);
            tx[0].re.num.limbs = &tx[1][0];
            tx[0].re.den.limbs = &tx[1][1];

            return div_(allocator, o, tx[0], y);
        },
        .float => {
            var tx = try @import("../float/asComplex.zig").asComplex(x);
            tx[0].re.num.limbs = &tx[1][0];
            tx[0].re.den.limbs = &tx[1][1];

            return div_(allocator, o, tx[0], y);
        },
        .int => {
            var tx = @import("../int/asComplex.zig").asComplex(x);
            tx[0].re.num.limbs = &tx[1];

            return div_(allocator, o, tx[0], y);
        },
        .bool => return div_(
            allocator,
            o,
            types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable,
            y,
        ),
        else => unreachable,
    }
}
