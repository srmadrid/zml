const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Compute the greatest common divisor of two numbers in `Integer` precision.
/// Float, rational or real types are truncated towards zero, and for cfloat or
/// complex types, only the real part is considered.
///
/// Signature
/// ---------
/// ```zig
/// fn gcd(allocator: std.mem.Allocator, x: X, y: Y) !Integer
/// ```
///
/// Parameters
/// ----------
/// `allocator` (`std.mem.Allocator`):
/// The allocator to use for memory allocations.
///
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// Returns
/// -------
/// `Integer`:
/// The greatest common divisor of `x` and `y`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
pub fn gcd(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("integer.gcd requires x and y to be numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(X)) {
        .expression => @compileError("integer. gcd not implemented for Expression yet"),
        .complex => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer. gcd not implemented for Complex + Expression yet"),
            .complex => return gcd(allocator, x.re, y.re),
            .real => @compileError("integer. gcd not implemented for Complex + Real yet"),
            .rational => return gcd(allocator, x.re, y),
            .integer => return gcd(allocator, x.re, y),
            .cfloat => return gcd(allocator, x.re, y.re),
            .float => return gcd(allocator, x.re, y),
            .int => return gcd(allocator, x.re, y),
            .bool => return gcd(allocator, x.re, y),
        },
        .real => @compileError("integer. gcd not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer. gcd not implemented for Rational + Expression yet"),
            .complex => return gcd(allocator, x, y.re),
            .real => @compileError("integer. gcd not implemented for Rational + Real yet"),
            .rational => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return gcd(allocator, tx, ty);
            },
            .integer => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                return gcd(allocator, tx, y);
            },
            .cfloat => return gcd(allocator, x, y.re),
            .float => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, tx, ty[0]);
            },
            .int => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, tx, ty[0]);
            },
            .bool => {
                var tx: Integer = try integer.div(allocator, x.num, x.den);
                defer tx.deinit(allocator);
                return gcd(allocator, tx, types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .integer => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer. gcd not implemented for Integer + Expression yet"),
            .complex => return gcd(allocator, x, y.re),
            .real => @compileError("integer. gcd not implemented for Integer + Real yet"),
            .rational => {
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return gcd(allocator, x, ty);
            },
            .integer => {
                if (x.size == 0) return integer.abs(allocator, y);
                if (y.size == 0) return integer.abs(allocator, x);

                var a = try integer.abs(allocator, x);
                errdefer a.deinit(allocator);

                var b = try integer.abs(allocator, y);
                defer b.deinit(allocator);

                // Remove common factors of 2 from both a and b, store in `shift`.
                const shift: u32 = int.min(trailingZeroBits(a), trailingZeroBits(b));

                var s: u32 = shift;
                while (s >= 32) : (s -= 32) {
                    removeFirstLimb(&a);
                    removeFirstLimb(&b);
                }

                if (s > 0) {
                    @import("div_.zig").shiftRightInPlace(a.limbs[0..a.size], @intCast(s));
                    a.truncate();
                    @import("div_.zig").shiftRightInPlace(b.limbs[0..b.size], @intCast(s));
                    b.truncate();
                }

                // Make sure a is odd.
                while (a.size != 0 and (a.limbs[0] & 1) == 0) {
                    @import("div_.zig").shiftRightInPlace(a.limbs[0..a.size], 1);
                    a.truncate();
                }

                while (b.size != 0) {
                    while (b.size != 0 and (b.limbs[0] & 1) == 0) {
                        @import("div_.zig").shiftRightInPlace(b.limbs[0..b.size], 1);
                        b.truncate();
                    }

                    if (b.size == 0) break;

                    a.truncate();
                    b.truncate();

                    if (integer.gt(a, b)) {
                        const tmp = a;
                        a = b;
                        b = tmp;
                    }

                    try integer.sub_(allocator, &b, b, a);
                }

                // Restore common factors of 2.
                const limb_shift = shift / 32;
                const bit_shift: u5 = @intCast(shift % 32);

                if (limb_shift > 0) try addZeroLimbs(allocator, &a, limb_shift);
                if (bit_shift > 0) {
                    if (@import("div_.zig").shiftLeftInPlace(a.limbs[0..a.size], &a.limbs[a.size], bit_shift))
                        a.size += 1;

                    a.truncate();
                }

                return a;
            },
            .cfloat => return gcd(allocator, x, y.re),
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, x, ty[0]);
            },
            .bool => {
                return gcd(allocator, x, types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer. gcd not implemented for CFloat + Expression yet"),
            .complex => return gcd(allocator, x.re, y.re),
            .real => @compileError("integer. gcd not implemented for CFloat + Real yet"),
            .rational => return gcd(allocator, x.re, y),
            .integer => return gcd(allocator, x.re, y),
            .cfloat => return gcd(allocator, x.re, y.re),
            .float => return gcd(allocator, x.re, y),
            .int => return gcd(allocator, x.re, y),
            .bool => return gcd(allocator, x.re, y),
        },
        .float => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer. gcd not implemented for Float + Expression yet"),
            .complex => return gcd(allocator, x, y.re),
            .real => @compileError("integer. gcd not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return gcd(allocator, tx[0], ty);
            },
            .integer => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return gcd(allocator, tx[0], y);
            },
            .cfloat => return gcd(allocator, x, y.re),
            .float => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return gcd(allocator, tx[0], types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer. gcd not implemented for Int + Expression yet"),
            .complex => return gcd(allocator, x, y.re),
            .real => @compileError("integer. gcd not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return gcd(allocator, tx[0], ty);
            },
            .integer => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return gcd(allocator, tx[0], y);
            },
            .cfloat => return gcd(allocator, x, y.re),
            .float => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, tx[0], ty[0]);
            },
            .int => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, tx[0], ty[0]);
            },
            .bool => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return gcd(allocator, tx[0], types.cast(Integer, y, .{}) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .expression => @compileError("integer. gcd not implemented for Bool + Expression yet"),
            .complex => return gcd(allocator, x, y.re),
            .real => @compileError("integer. gcd not implemented for Bool + Real yet"),
            .rational => {
                var ty: Integer = try integer.div(allocator, y.num, y.den);
                defer ty.deinit(allocator);
                return gcd(allocator, types.cast(Integer, x, .{}) catch unreachable, ty);
            },
            .integer => {
                return gcd(allocator, types.cast(Integer, x, .{}) catch unreachable, y);
            },
            .cfloat => return gcd(allocator, x, y.re),
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, types.cast(Integer, x, .{}) catch unreachable, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];
                return gcd(allocator, types.cast(Integer, x, .{}) catch unreachable, ty[0]);
            },
            .bool => {
                return gcd(
                    allocator,

                    types.cast(Integer, x, .{}) catch unreachable,
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
    }
}

fn trailingZeroBits(self: Integer) u32 {
    var count: u32 = 0;
    var i: u32 = 0;
    while (i < self.size) : (i += 1) {
        const limb = self.limbs[i];
        if (limb == 0) {
            count += 32;
        } else {
            count += @ctz(limb);
            break;
        }
    }
    return count;
}

fn removeFirstLimb(self: *Integer) void {
    if (self.size == 0) return;
    var i: u32 = 0;
    while (i + 1 < self.size) : (i += 1) {
        self.limbs[i] = self.limbs[i + 1];
    }
    self.size -= 1;
    if (self.size > 0) self.limbs[self.size] = 0;
    self.truncate();
}

fn addZeroLimbs(allocator: std.mem.Allocator, self: *Integer, n: u32) !void {
    if (n == 0 or self.size == 0) return;

    try self.reserve(allocator, self.size + n);

    var i: u32 = self.size;
    while (i > 0) : (i -= 1) {
        self.limbs[i + n - 1] = self.limbs[i - 1];
    }

    var j: u32 = 0;
    while (j < n) : (j += 1) {
        self.limbs[j] = 0;
    }

    self.size += n;
}
