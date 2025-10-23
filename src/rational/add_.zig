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
