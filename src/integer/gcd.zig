const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn gcd(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .integer and types.numericType(Y) == .int) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .float) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .int and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .float and types.numericType(Y) == .integer))
        @compileError("integer.gcd requires at least one of x or y to be an integer, the other must be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

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
            .float, .int => {
                var temp: Integer = try .initSet(allocator, y);
                defer temp.deinit(allocator);
                return gcd(allocator, x, temp);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .integer => {
                var temp: Integer = try .initSet(allocator, x);
                defer temp.deinit(allocator);
                return gcd(allocator, temp, y);
            },
            else => unreachable,
        },
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
