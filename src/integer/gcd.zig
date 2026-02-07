const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Computes the greatest common divisor between two operands of integer,
/// cfloat, dyadic, float, int or bool types, where at least one operand must be
/// of integer type. The operation is performed by casting both operands to
/// integer, then applying the binary GCD algorithm.
///
/// ## Signature
/// ```zig
/// integer.gcd(allocator: std.mem.Allocator, x: X, y: Y) !Integer
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `Integer`: The greatest common divisor of `x` and `y`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
pub fn gcd(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.integer) or !types.numericType(Y).le(.integer) or
        (types.numericType(X) != .integer and types.numericType(Y) != .integer))
        @compileError("zml.integer.gcd: at least one of x or y must be an integer, the other must be a bool, an int, a float, a dyadic, a cfloat or an integer, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
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
            .dyadic => {
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return gcd(allocator, x, ty[0]);
            },
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
            .bool => return gcd(
                allocator,

                x,
                types.cast(Integer, y, .{}) catch unreachable,
            ),
            else => unreachable,
        },
        .cfloat => return gcd(allocator, x.re, y),
        .dyadic => {
            var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return gcd(allocator, tx[0], y);
        },
        .float => {
            var tx = try @import("../float/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return gcd(allocator, tx[0], y);
        },
        .int => {
            var tx = @import("../int/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return gcd(allocator, tx[0], y);
        },
        .bool => return gcd(
            allocator,

            types.cast(Integer, x, .{}) catch unreachable,
            y,
        ),
        else => unreachable,
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
