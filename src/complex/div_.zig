const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

const check_aliasing = @import("check_aliasing.zig").check_aliasing;

pub fn div_(allocator: std.mem.Allocator, o: anytype, x: anytype, y: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("complex.div_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.numericType(O) != .complex)
        @compileError("complex.div_ requires o to be a complex, got " ++ @typeName(O));

    comptime if (types.numericType(X) != .complex and types.numericType(X) != .real and types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .complex and types.numericType(X) != .real and types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.div_ requires x and y to be an int, float, integer, rational, real or complex, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .complex => switch (comptime types.numericType(Y)) {
            .complex => {
                if (x.num.size == 0) {
                    try o.re.set(allocator, y.re);
                    try o.im.set(allocator, y.im);

                    return;
                }

                if (y.num.size == 0) {
                    try o.re.set(allocator, x.re);
                    try o.im.set(allocator, x.im);

                    return;
                }

                // Aliasing checks
                var tx: X = try check_aliasing(allocator, o, x);
                defer tx.deinit(allocator);
                var ty: Y = try check_aliasing(allocator, o, y);
                defer ty.deinit(allocator);

                // (a + bi) / (c + di) = ((ac + bd) / (c² + d²)) + ((bc - ad) / (c² + d²))i
                var ac = try ops.div(x.re, y.re, .{ .allocator = allocator });
                defer ac.deinit(allocator);
                var bd = try ops.div(x.im, y.im, .{ .allocator = allocator });
                defer bd.deinit(allocator);

                var c2 = try ops.mul(y.re, y.re, .{ .allocator = allocator });
                defer c2.deinit(allocator);
                var d2 = try ops.mul(y.im, y.im, .{ .allocator = allocator });
                defer d2.deinit(allocator);

                try ops.add_(&c2, c2, d2, .{ .allocator = allocator });

                try ops.add_(&ac, ac, bd, .{ .allocator = allocator });
                try ops.div_(&o.re, ac, c2, .{ .allocator = allocator });

                try ops.mul_(&ac, x.im, y.re, .{ .allocator = allocator });
                try ops.mul_(&bd, x.re, y.im, .{ .allocator = allocator });

                try ops.sub_(&ac, ac, bd, .{ .allocator = allocator });
                try ops.div_(&o.im, ac, c2, .{ .allocator = allocator });

                return;
            },
            .rational, .integer => {
                return div_(allocator, o, x, y.asComplex());
            },
            .float, .int => {
                var temp: complex.Complex(rational.Rational) = try .initSet(allocator, y, 0);
                defer temp.deinit(allocator);
                return div_(allocator, o, x, temp);
            },
            else => unreachable,
        },
        .rational, .integer => switch (comptime types.numericType(Y)) {
            .complex => {
                return div_(allocator, o, x.asComplex(), y);
            },
            .rational, .integer => {
                return div_(allocator, o, x.asComplex(), y.asComplex());
            },
            .float, .int => {
                var temp: complex.Complex(rational.Rational) = try .initSet(allocator, y, 0);
                defer temp.deinit(allocator);
                return div_(allocator, o, x.asComplex(), temp);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .complex => {
                var temp: complex.Complex(rational.Rational) = try .initSet(allocator, x, 0);
                defer temp.deinit(allocator);
                return div_(allocator, o, temp, y);
            },
            .rational, .integer => {
                var temp: complex.Complex(rational.Rational) = try .initSet(allocator, x, 0);
                defer temp.deinit(allocator);
                return div_(allocator, o, temp, y.asComplex());
            },
            .float, .int => {
                var tx: complex.Complex(rational.Rational) = try .initSet(allocator, x, 0);
                defer tx.deinit(allocator);
                var ty: complex.Complex(rational.Rational) = try .initSet(allocator, y, 0);
                defer ty.deinit(allocator);
                return div_(allocator, o, tx, ty);
            },
            else => unreachable,
        },

        else => unreachable,
    }
}
