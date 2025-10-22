const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

pub fn add_(allocator: std.mem.Allocator, o: anytype, x: anytype, y: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("complex.add_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.numericType(O) != .complex)
        @compileError("complex.add_ requires o to be a complex, got " ++ @typeName(O));

    comptime if (types.numericType(X) != .complex and types.numericType(X) != .real and types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .complex and types.numericType(X) != .real and types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.add_ requires x and y to be an int, float, integer, rational, real or complex, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .complex => switch (comptime types.numericType(Y)) {
            .complex => {
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

                // (a + bi) + (c + di) = (a + c) + (b + d)i
                try ops.add_(allocator, &o.re, x.re, y.re);
                try ops.add_(allocator, &o.im, x.im, y.im);

                try o.reduce(allocator);

                return;
            },
            .integer => {
                return add_(allocator, o, x, y.asComplex());
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
