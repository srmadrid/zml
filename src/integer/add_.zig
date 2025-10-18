const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn add_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("integer.add_ requires x and y to be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return integer.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                if (x.size == 0) {
                    try o.set(allocator, y);
                    return;
                }

                if (y.size == 0) {
                    try o.set(allocator, x);
                    return;
                }

                // Aliasing checks
                var tx: Integer = if (std.meta.eql(o.*, x))
                    try x.copy(allocator)
                else blk: {
                    var tmp: Integer = x;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                defer tx.deinit(allocator);
                var ty: Integer = if (std.meta.eql(o.*, y))
                    try y.copy(allocator)
                else blk: {
                    var tmp: Integer = y;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                defer ty.deinit(allocator);

                if (tx.positive == ty.positive) {
                    try add_abs(allocator, o, tx, ty);
                    o.positive = tx.positive;
                    return;
                }

                // Reset o's size
                o.size = 0;

                const cmp_abs: types.Cmp =
                    integer.cmp(
                        integer.abs(null, tx) catch unreachable,
                        integer.abs(null, ty) catch unreachable,
                    );

                if (cmp_abs == .eq) {
                    return o.set(allocator, 0);
                } else if (cmp_abs == .gt) {
                    try sub_abs(allocator, o, tx, ty);
                    o.positive = tx.positive;
                } else {
                    try sub_abs(allocator, o, ty, tx);
                    o.positive = ty.positive;
                }

                return;
            },
            .float, .int => {
                var temp: Integer = try .initSet(allocator, y);
                defer temp.deinit(allocator);
                return add_(allocator, o, x, temp);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .float, .int => {
                var tx: Integer = try .initSet(allocator, x);
                defer tx.deinit(allocator);
                var ty: Integer = try .initSet(allocator, y);
                defer ty.deinit(allocator);
                return add_(allocator, o, tx, ty);
            },
            .integer => {
                var temp: Integer = try .initSet(allocator, x);
                defer temp.deinit(allocator);
                return add_(allocator, o, temp, y);
            },
            else => unreachable,
        },
        else => unreachable,
    }
}

fn add_abs(allocator: std.mem.Allocator, o: *Integer, x: Integer, y: Integer) !void {
    // No aliasing allowed.
    try o.reserve(allocator, int.max(x.size, y.size) + 1);
    var carry: u64 = 0;
    var i: u32 = 0;
    while (i < x.size or i < y.size or carry != 0) : (i += 1) {
        const xi: u64 = if (i < x.size) x.limbs[i] else 0;
        const yi: u64 = if (i < y.size) y.limbs[i] else 0;

        const s: u64 = xi + yi + carry;
        o.limbs[i] = @truncate(s & 0xFFFFFFFF);
        carry = (s >> 32);
    }

    o.size = i;
    o.trimSize();

    return;
}

fn sub_abs(allocator: std.mem.Allocator, o: *Integer, x: Integer, y: Integer) !void {
    // Assumes a >= b >= 0 in absolute value. No aliasing allowed.
    try o.reserve(allocator, x.size);
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

        o.limbs[i] = @truncate(s & 0xFFFFFFFF);
    }

    o.size = i;
    o.trimSize();

    return;
}
