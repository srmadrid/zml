const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn mul(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if ((types.numericType(X) != .integer and types.numericType(X) != .int) or
        (types.numericType(X) != .integer and types.numericType(X) != .float) or
        (types.numericType(X) != .int and types.numericType(X) != .integer) or
        (types.numericType(X) != .float and types.numericType(X) != .integer))
        @compileError("integer.mul requires at least one of x or y to be an integer, the other must be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                if (x.size == 0 or y.size == 0) return .init(allocator, 2);
                var result: Integer = try .init(allocator, x.size + y.size);

                var i: u32 = 0;
                while (i < x.size + y.size) : (i += 1) {
                    result.limbs[i] = 0;
                }

                i = 0;
                while (i < x.size) : (i += 1) {
                    var carry: u128 = 0;
                    var j: u32 = 0;
                    while (j < y.size) : (j += 1) {
                        const acc: u128 = types.scast(u128, result.limbs[i + j]) +
                            types.scast(u128, x.limbs[i]) * types.scast(u128, y.limbs[j]) + carry;
                        result.limbs[i + j] = @truncate(acc & 0xFFFFFFFF);
                        carry = acc >> 32;
                    }

                    // Propagate leftover carry.
                    var k: u32 = i + y.size;
                    while (carry != 0) : (k += 1) {
                        const acc: u128 = types.scast(u128, result.limbs[k]) + carry;
                        result.limbs[k] = @truncate(acc & 0xFFFFFFFF);
                        carry = acc >> 32;
                    }
                }

                result.size = x.size + y.size;
                result.positive = x.positive == y.positive;
                return result;
            },
            .float => {
                // Check y is integer -> y == floor(y)
                // If not, use rational.add
            },
            .int => {},
            else => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .integer => {
                // Check x is integer -> x == floor(x)
            },
            else => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .integer => {},
            else => unreachable,
        },
        else => unreachable,
    }
}
