const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn add(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if ((types.numericType(X) != .integer and types.numericType(X) != .int) or
        (types.numericType(X) != .integer and types.numericType(X) != .float) or
        (types.numericType(X) != .int and types.numericType(X) != .integer) or
        (types.numericType(X) != .float and types.numericType(X) != .integer))
        @compileError("integer.add requires at least one of x or y to be an integer, the other must be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                if (x.size == 0) return y.copy(allocator);
                if (y.size == 0) return x.copy(allocator);

                if (x.positive == y.positive) {
                    var result: Integer = try add_abs(allocator, x, y);
                    result.positive = x.positive;
                    result.trimSize();
                    return result;
                }

                const cmp_abs: types.Cmp = integer.cmp(integer.abs(null, x) catch unreachable, integer.abs(null, y) catch unreachable);
                if (cmp_abs == .eq) return .init(allocator, 2); // zero
                if (cmp_abs == .gt) {
                    var result: Integer = try sub_abs_pos_geq(allocator, x, y);
                    result.positive = x.positive;
                    result.trimSize();
                    return result;
                } else {
                    var result: Integer = try sub_abs_pos_geq(allocator, y, x);
                    result.positive = y.positive;
                    result.trimSize();
                    return result;
                }
            },
            .float, .int => {
                var temp: Integer = try .initSet(allocator, y);
                defer temp.deinit(allocator);
                return add(allocator, x, temp);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .integer => {
                var temp: Integer = try .initSet(allocator, x);
                defer temp.deinit(allocator);
                return add(allocator, temp, y);
            },
            else => unreachable,
        },
        else => unreachable,
    }
}

fn add_abs(allocator: std.mem.Allocator, x: Integer, y: Integer) !Integer {
    var result: Integer = try .init(allocator, (if (x.size > y.size) x.size else y.size) + 1);
    var carry: u64 = 0;
    var i: u32 = 0;
    while (i < x.size or i < y.size or carry != 0) : (i += 1) {
        const xi: u64 = if (i < x.size) x.limbs[i] else 0;
        const yi: u64 = if (i < y.size) y.limbs[i] else 0;

        const s: u64 = xi + yi + carry;
        result.limbs[i] = @truncate(s & 0xFFFFFFFF);
        carry = (s >> 32);
    }

    result.size = i;

    return result;
}

fn sub_abs_pos_geq(allocator: std.mem.Allocator, x: Integer, y: Integer) !Integer {
    // Assumes a >= b >= 0 in absolute value.
    var result: Integer = try .init(allocator, x.size);
    var borrow: u64 = 0;
    var i: u32 = 0;
    while (i < x.size) : (i += 1) {
        const xi: u64 = x.limbs[i];
        const yi: u64 = if (i < y.size) y.limbs[i] else 0;

        var s: u64 = undefined;
        if (xi < yi + borrow) {
            s = xi + (@as(u64, 1) << 32) - yi - borrow;
            borrow = 1;
        } else {
            s = xi - yi - borrow;
            borrow = 0;
        }

        result.limbs[result.size] = @truncate(s & 0xFFFFFFFF);
        result.size += 1;
    }

    return result;
}
