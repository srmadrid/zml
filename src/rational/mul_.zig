const std = @import("std");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

pub fn mul_(allocator: std.mem.Allocator, o: *Rational, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.mul_ requires x and y to be an int, float, integer or rational, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

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

                if (o.num.flags.writable and o.den.flags.writable)
                    try o.reduce(allocator);

                return;
            },
            .integer => {
                return mul_(allocator, o, x, y.asRational());
            },
            .float, .int => {
                var temp: Rational = try .initSet(allocator, y, 1);
                defer temp.deinit(allocator);
                return mul_(allocator, o, x, temp);
            },
            else => unreachable,
        },
        .integer => switch (comptime types.numericType(Y)) {
            .rational => {
                return mul_(allocator, o, x.asRational(), y);
            },
            .integer => {
                return mul_(allocator, o, x.asRational(), y.asRational());
            },
            .float, .int => {
                var temp: Rational = try .initSet(allocator, y, 1);
                defer temp.deinit(allocator);
                return mul_(allocator, o, x.asRational(), temp);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .rational => {
                var temp: Rational = try .initSet(allocator, x, 1);
                defer temp.deinit(allocator);
                return mul_(allocator, o, temp, y);
            },
            .integer => {
                var temp: Rational = try .initSet(allocator, x, 1);
                defer temp.deinit(allocator);
                return mul_(allocator, o, temp, y.asRational());
            },
            .float, .int => {
                var tx: Rational = try .initSet(allocator, x, 1);
                defer tx.deinit(allocator);
                var ty: Rational = try .initSet(allocator, y, 1);
                defer ty.deinit(allocator);
                return mul_(allocator, o, tx, ty);
            },
            else => unreachable,
        },

        else => unreachable,
    }
}
