const std = @import("std");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

const check_aliasing_alloc = @import("check_aliasing_alloc.zig").check_aliasing_alloc;

/// Performs in-place multiplication between two operands of rational, integer,
/// cfloat, dyadic, float, int or bool types, where at least one operand must be
/// of rational type. The operation is performed by casting both operands to
/// integer, then multiplying them in-place.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// ## Signature
/// ```zig
/// rational.mul_(allocator: std.mem.Allocator, o: *Rational, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations. Must be the same allocator used to initialize `o`.
/// * `o` (`*Rational`): A pointer to the output operand where the result will
///   be stored.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `rational.Error.NotWritable`: If the output operand `o` is not writable.
/// * `rational.Error.DataNotOwned`: If the output operand `o` does not own its
///   data and resizing is needed.
/// * `integer.Error.NotWritable`: If the numerator or denominator of the output
///   operand `o` is not writable when it needs to be modified.
/// * `integer.Error.DataNotOwned`: If the numerator or denominator of the
///   output operand `o` does not own its data and resizing is needed.
pub fn mul_(allocator: std.mem.Allocator, o: *Rational, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.rational) or !types.numericType(Y).le(.rational) or
        (types.numericType(X) != .rational and types.numericType(Y) != .rational))
        @compileError("zml.rational.mul_: at least one of x or y must be a rational, the other must be a bool, an int, a float, a dyadic, a cfloat, an integer or a rational, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .rational => switch (comptime types.numericType(Y)) {
            .rational => {
                if (x.num.size == 0 or y.num.size == 0) {
                    if (integer.ne(o.num, y.num)) {
                        if (!o.num.flags.writable)
                            return rational.Error.NotWritable;

                        try o.num.set(allocator, 0);
                    }

                    if (integer.ne(o.den, y.den)) {
                        if (!o.den.flags.writable)
                            return rational.Error.NotWritable;

                        try o.den.set(allocator, 1);
                    }

                    return;
                }

                // Aliasing checks
                var tx: Rational = try check_aliasing_alloc(allocator, o, x);
                defer tx.deinit(allocator);
                var ty: Rational = try check_aliasing_alloc(allocator, o, y);
                defer ty.deinit(allocator);

                // a/b * c/d = (a*c)/(b*d)
                var temp: integer.Integer = try .init(allocator, 0);
                defer temp.deinit(allocator);

                if (integer.eq(tx.num, constants.one(integer.Integer, .{}) catch unreachable)) {
                    if (integer.ne(o.num, ty.num)) {
                        if (!o.num.flags.writable)
                            return rational.Error.NotWritable;

                        try o.num.set(allocator, ty.num);
                    }
                } else if (integer.eq(ty.num, constants.one(integer.Integer, .{}) catch unreachable)) {
                    if (integer.ne(o.num, tx.num)) {
                        if (!o.num.flags.writable)
                            return rational.Error.NotWritable;

                        try o.num.set(allocator, tx.num);
                    }
                } else {
                    try integer.mul_(allocator, &temp, tx.num, ty.num);

                    if (integer.ne(o.num, temp)) {
                        if (!o.num.flags.writable)
                            return rational.Error.NotWritable;

                        try o.num.set(allocator, temp);
                    }
                }

                if (integer.eq(tx.den, constants.one(integer.Integer, .{}) catch unreachable)) {
                    if (integer.ne(o.den, ty.den)) {
                        if (!o.den.flags.writable)
                            return rational.Error.NotWritable;

                        try o.den.set(allocator, ty.den);
                    }
                } else if (integer.eq(ty.den, constants.one(integer.Integer, .{}) catch unreachable)) {
                    if (integer.ne(o.den, tx.den)) {
                        if (!o.den.flags.writable)
                            return rational.Error.NotWritable;

                        try o.den.set(allocator, tx.den);
                    }
                } else {
                    try integer.mul_(allocator, &temp, tx.den, ty.den);

                    if (integer.ne(o.den, temp)) {
                        if (!o.den.flags.writable)
                            return rational.Error.NotWritable;

                        try o.den.set(allocator, temp);
                    }
                }

                o.num.positive = (tx.num.positive == ty.num.positive);

                if (o.num.flags.writable and o.den.flags.writable)
                    try o.reduce(allocator);

                return;
            },
            .integer => return mul_(allocator, o, x, y.asRational()),
            .cfloat => return mul_(allocator, o, x, y.re),
            .dyadic => {
                var ty = @import("../dyadic/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];

                return mul_(allocator, o, x, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];

                return mul_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];

                return mul_(allocator, o, x, ty[0]);
            },
            .bool => return mul_(
                allocator,
                o,
                x,
                types.cast(Rational, y, .{}) catch unreachable,
            ),
            else => unreachable,
        },
        .integer => return mul_(allocator, o, x.asRational(), y),
        .cfloat => return mul_(allocator, o, x.re, y),
        .dyadic => {
            var tx = @import("../dyadic/asRational.zig").asRational(x);
            tx[0].num.limbs = &tx[1][0];
            tx[0].den.limbs = &tx[1][1];

            return mul_(allocator, o, tx[0], y);
        },
        .float => {
            var tx = try @import("../float/asRational.zig").asRational(x);
            tx[0].num.limbs = &tx[1][0];
            tx[0].den.limbs = &tx[1][1];

            return mul_(allocator, o, tx[0], y);
        },
        .int => {
            var tx = @import("../int/asRational.zig").asRational(x);
            tx[0].num.limbs = &tx[1];

            return mul_(allocator, o, tx[0], y);
        },
        .bool => return mul_(
            allocator,
            o,
            types.cast(Rational, x, .{}) catch unreachable,
            y,
        ),
        else => unreachable,
    }
}
