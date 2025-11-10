const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Performs in-place addition between two operands of any numeric type in
/// `Rational` precision. For cfloat or complex types, only the real part is
/// considered.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// Signature
/// ---------
/// ```zig
/// fn add_(allocator: std.mem.Allocator, o: *Rational, x: X, y: Y, ctx: anytype) !void
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
pub fn add_(allocator: std.mem.Allocator, o: *Rational, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("rational.add_ requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .complex => switch (comptime types.numericType(Y)) {
            .complex => return add_(allocator, o, x.re, y.re),
            .real => @compileError("rational.add_ not implemented for Complex + Real yet"),
            .rational => return add_(allocator, o, x.re, y),
            .integer => return add_(allocator, o, x.re, y),
            .cfloat => return add_(allocator, o, x.re, y.re),
            .float => return add_(allocator, o, x.re, y),
            .int => return add_(allocator, o, x.re, y),
            .bool => return add_(allocator, o, x.re, y),
        },
        .real => @compileError("rational.add_ not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .complex => return add_(allocator, o, x, y.re),
            .real => @compileError("rational.add_ not implemented for Rational + Real yet"),
            .rational => {
                if (x.num.size == 0) {
                    if (integer.ne(o.num, y.num)) {
                        if (!o.num.flags.writable)
                            return rational.Error.NotWritable;

                        try o.num.set(allocator, y.num);
                    }

                    if (integer.ne(o.den, y.den)) {
                        if (!o.den.flags.writable)
                            return rational.Error.NotWritable;

                        try o.den.set(allocator, y.den);
                    }

                    return;
                }

                if (y.num.size == 0) {
                    if (integer.ne(o.num, x.num)) {
                        if (!o.num.flags.writable)
                            return rational.Error.NotWritable;

                        try o.num.set(allocator, x.num);
                    }

                    if (integer.ne(o.den, x.den)) {
                        if (!o.den.flags.writable)
                            return rational.Error.NotWritable;

                        try o.den.set(allocator, x.den);
                    }

                    return;
                }

                const cmpox: types.Cmp = integer.cmp(o.den, x.den);
                const cmpxy: types.Cmp = integer.cmp(x.den, y.den);

                if (cmpox != .eq or cmpxy != .eq) // If the result denominator won't change, it can be non-writable
                    if (!o.den.flags.writable)
                        return rational.Error.NotWritable;

                // Aliasing checks
                var tx: Rational = undefined;
                tx.num = if (o.num.limbs == x.num.limbs or o.den.limbs == x.num.limbs)
                    try x.num.copy(allocator)
                else blk: {
                    var tmp: integer.Integer = x.num;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                tx.den = if (o.num.limbs == x.den.limbs or o.den.limbs == x.den.limbs)
                    try x.den.copy(allocator)
                else blk: {
                    var tmp: integer.Integer = x.den;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                tx.flags = .{ .owns_data = true, .writable = false };
                defer tx.deinit(allocator);
                var ty: Rational = undefined;
                ty.num = if (o.num.limbs == y.num.limbs or o.den.limbs == y.num.limbs)
                    try y.num.copy(allocator)
                else blk: {
                    var tmp: integer.Integer = y.num;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                ty.den = if (o.num.limbs == y.den.limbs or o.den.limbs == y.den.limbs)
                    try y.den.copy(allocator)
                else blk: {
                    var tmp: integer.Integer = y.den;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                ty.flags = .{ .owns_data = true, .writable = false };
                defer ty.deinit(allocator);

                if (cmpxy == .eq) {
                    try integer.add_(allocator, &o.num, tx.num, ty.num);

                    if (cmpox != .eq) { // o.den is guaranteed writable
                        try o.den.set(allocator, tx.den);
                        try o.reduce(allocator);
                    }

                    return;
                }

                // a/b + c/d = (a*d + b*c) / (b*d)
                var ad: integer.Integer = try integer.mul(allocator, tx.num, ty.den);
                defer ad.deinit(allocator);

                try integer.mul_(allocator, &o.den, tx.den, ty.num);
                try integer.add_(allocator, &o.den, ad, o.den);

                if (integer.ne(o.num, o.den)) {
                    if (!o.num.flags.writable)
                        return rational.Error.NotWritable;

                    try o.num.set(allocator, o.den);
                }

                try integer.mul_(allocator, &o.num, tx.den, ty.num);
                try integer.add_(allocator, &o.num, o.num, ad);

                try integer.mul_(allocator, &o.den, tx.den, ty.den);

                if (o.num.flags.writable) // o.den is guaranteed writable
                    try o.reduce(allocator);

                return;
            },
            .integer => return add_(allocator, o, x, y.asRational()),
            .cfloat => return add_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return add_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return add_(allocator, o, x, ty[0]);
            },
            .bool => return add_(allocator, o, x, types.cast(Rational, y, .{}) catch unreachable),
        },
        .integer => switch (comptime types.numericType(Y)) {
            .complex => return add_(allocator, o, x, y.re),
            .real => @compileError("rational.add_ not implemented for Integer + Real yet"),
            .rational => return add_(allocator, o, x.asRational(), y),
            .integer => return add_(allocator, o, x.asRational(), y.asRational()),
            .cfloat => return add_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return add_(allocator, o, x.asRational(), ty[0]);
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return add_(allocator, o, x.asRational(), ty[0]);
            },
            .bool => {
                return add_(allocator, o, x.asRational(), types.cast(Rational, y, .{}) catch unreachable);
            },
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .complex => return add_(allocator, o, x.re, y.re),
            .real => @compileError("rational.add_ not implemented for CFloat + Real yet"),
            .rational => return add_(allocator, o, x.re, y),
            .integer => return add_(allocator, o, x.re, y),
            .cfloat => return add_(allocator, o, x.re, y.re),
            .float => return add_(allocator, o, x.re, y),
            .int => return add_(allocator, o, x.re, y),
            .bool => return add_(allocator, o, x.re, y),
        },
        .float => switch (comptime types.numericType(Y)) {
            .complex => return add_(allocator, o, x, y.re),
            .real => @compileError("rational.add_ not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                return add_(allocator, o, tx[0], y);
            },
            .integer => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                return add_(allocator, o, tx[0], y.asRational());
            },
            .cfloat => return add_(allocator, o, x, y.re),
            .float => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../float/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1][0];
                tx[0].den.limbs = &tx[1][1];
                return add_(allocator, o, tx[0], types.cast(Rational, y, .{}) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .complex => return add_(allocator, o, x, y.re),
            .real => @compileError("rational.add_ not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                return add_(allocator, o, tx[0], y);
            },
            .integer => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                return add_(allocator, o, tx[0], y.asRational());
            },
            .cfloat => return add_(allocator, o, x, y.re),
            .float => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return add_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = @import("../int/asRational.zig").asRational(x);
                tx[0].num.limbs = &tx[1];
                return add_(allocator, o, tx[0], types.cast(Rational, y, .{}) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .complex => return add_(allocator, o, x, y.re),
            .real => @compileError("rational.add_ not implemented for Bool + Real yet"),
            .rational => return add_(allocator, o, types.cast(Rational, x, .{}) catch unreachable, y),
            .integer => return add_(allocator, o, types.cast(Rational, x, .{}) catch unreachable, y.asRational()),
            .cfloat => return add_(allocator, o, x, y.re),
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];
                return add_(allocator, o, types.cast(Rational, x, .{}) catch unreachable, ty[0]);
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];
                return add_(allocator, o, types.cast(Rational, x, .{}) catch unreachable, ty[0]);
            },
            .bool => {
                return add_(
                    allocator,
                    o,
                    types.cast(Rational, x, .{}) catch unreachable,
                    types.cast(Rational, y, .{}) catch unreachable,
                );
            },
        },
    }
}
