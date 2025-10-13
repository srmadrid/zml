const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn mul_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("integer.mul_ requires x and y to be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                if (x.size == 0 or y.size == 0) return o.set(allocator, 0);

                try o.reserve(allocator, x.size + y.size);

                var i: u32 = 0;
                while (i < x.size + y.size) : (i += 1) {
                    o.limbs[i] = 0;
                }

                i = 0;
                while (i < x.size) : (i += 1) {
                    var carry: u128 = 0;
                    var j: u32 = 0;
                    while (j < y.size) : (j += 1) {
                        const acc: u128 = types.scast(u128, o.limbs[i + j]) +
                            types.scast(u128, x.limbs[i]) * types.scast(u128, y.limbs[j]) + carry;
                        o.limbs[i + j] = @truncate(acc & 0xFFFFFFFF);
                        carry = acc >> 32;
                    }

                    // Propagate leftover carry.
                    var k: u32 = i + y.size;
                    while (carry != 0) : (k += 1) {
                        const acc: u128 = types.scast(u128, o.limbs[k]) + carry;
                        o.limbs[k] = @truncate(acc & 0xFFFFFFFF);
                        carry = acc >> 32;
                    }
                }

                o.size = x.size + y.size;
                o.trimSize();
                o.positive = x.positive == y.positive;
            },
            .float, .int => {
                var temp: Integer = try .initSet(allocator, y);
                defer temp.deinit(allocator);
                return mul_(allocator, o, x, temp);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .float, .int => {
                var tx: Integer = try .initSet(allocator, x);
                defer tx.deinit(allocator);
                var ty: Integer = try .initSet(allocator, y);
                defer ty.deinit(allocator);
                return mul_(allocator, o, tx, ty);
            },
            .integer => {
                var temp: Integer = try .initSet(allocator, x);
                defer temp.deinit(allocator);
                return mul_(allocator, o, temp, y);
            },
            else => unreachable,
        },
        else => unreachable,
    }
}
