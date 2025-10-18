const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

pub fn add_(allocator: std.mem.Allocator, o: *Rational, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.add_ requires x and y to be an int, float, integer or rational, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .rational => switch (comptime types.numericType(Y)) {
            .rational => {
                if (x.num.size == 0) {
                    try o.num.set(allocator, y.num);
                    try o.den.set(allocator, y.den);
                    return;
                }

                if (y.num.size == 0) {
                    try o.num.set(allocator, x.num);
                    try o.den.set(allocator, x.den);
                    return;
                }

                // Aliasing checks
                var tx: Rational = if (std.meta.eql(o.*, x))
                    try x.copy(allocator)
                else blk: {
                    var tmp: Rational = x;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                defer tx.deinit(allocator);
                var ty: Rational = if (std.meta.eql(o.*, y))
                    try y.copy(allocator)
                else blk: {
                    var tmp: Rational = y;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                defer ty.deinit(allocator);

                if (integer.eq(tx.den, ty.den)) {
                    try integer.add_(allocator, &o.num, tx.num, ty.num);

                    if (integer.ne(o.den, tx.den))
                        try o.den.set(allocator, tx.den);

                    o.num.positive = tx.num.positive;
                    return;
                }

                // a/b + c/d = (a*d + b*c) / (b*d)
                var ad: integer.Integer = try integer.mul(allocator, tx.num, ty.den);
                defer ad.deinit(allocator);

                try integer.mul_(allocator, &o.num, tx.den, ty.num);
                try integer.add_(allocator, &o.num, o.num, ad);

                try integer.mul_(allocator, &o.den, tx.den, ty.den);

                try o.reduce(allocator);

                return;
            },
            .integer => {
                return add_(allocator, o, x, y.asRational());
            },
            .float, .int => {
                var temp: Rational = try .initSet(allocator, y, 1);
                defer temp.deinit(allocator);
                return add_(allocator, o, x, temp);
            },
            else => unreachable,
        },
        .integer => switch (comptime types.numericType(Y)) {
            .rational => {
                return add_(allocator, o, x.asRational(), y);
            },
            .integer => {
                return add_(allocator, o, x.asRational(), y.asRational());
            },
            .float, .int => {
                var temp: Rational = try .initSet(allocator, y, 1);
                defer temp.deinit(allocator);
                return add_(allocator, o, x.asRational(), temp);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .rational => {
                var temp: Rational = try .initSet(allocator, x, 1);
                defer temp.deinit(allocator);
                return add_(allocator, o, temp, y);
            },
            .integer => {
                var temp: Rational = try .initSet(allocator, x, 1);
                defer temp.deinit(allocator);
                return add_(allocator, o, temp, y.asRational());
            },
            .float, .int => {
                var tx: Rational = try .initSet(allocator, x, 1);
                defer tx.deinit(allocator);
                var ty: Rational = try .initSet(allocator, y, 1);
                defer ty.deinit(allocator);
                return add_(allocator, o, tx, ty);
            },
            else => unreachable,
        },

        else => unreachable,
    }
}
