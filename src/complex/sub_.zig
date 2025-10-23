const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

const check_aliasing = @import("check_aliasing.zig").check_aliasing;

pub fn sub_(allocator: std.mem.Allocator, o: anytype, x: anytype, y: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("complex.sub_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.numericType(O) != .complex)
        @compileError("complex.sub_ requires o to be a complex, got " ++ @typeName(O));

    comptime if (types.numericType(X) != .complex and types.numericType(X) != .real and types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .complex and types.numericType(X) != .real and types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.sub_ requires x and y to be an int, float, integer, rational, real or complex, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .complex => switch (comptime types.numericType(Y)) {
            .complex => {
                return complex.add_(allocator, o, x, complex.neg(null, y) catch unreachable);
            },
            .rational, .integer => {
                return complex.add_(allocator, o, x, complex.neg(null, y.asComplex()) catch unreachable);
            },
            .float, .int => {
                var temp: complex.Complex(rational.Rational) = try .initSet(allocator, y, 0);
                defer temp.deinit(allocator);
                return complex.add_(allocator, o, x, complex.neg(null, temp) catch unreachable);
            },
            else => unreachable,
        },
        .rational, .integer => switch (comptime types.numericType(Y)) {
            .complex => {
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y) catch unreachable);
            },
            .rational, .integer => {
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y.asComplex()) catch unreachable);
            },
            .float, .int => {
                var temp: complex.Complex(rational.Rational) = try .initSet(allocator, y, 0);
                defer temp.deinit(allocator);
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, temp) catch unreachable);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .complex => {
                var temp: complex.Complex(rational.Rational) = try .initSet(allocator, x, 0);
                defer temp.deinit(allocator);
                return complex.add_(allocator, o, temp, complex.neg(null, y) catch unreachable);
            },
            .rational, .integer => {
                var temp: complex.Complex(rational.Rational) = try .initSet(allocator, x, 0);
                defer temp.deinit(allocator);
                return complex.add_(allocator, o, temp, complex.neg(null, y.asComplex()) catch unreachable);
            },
            .float, .int => {
                var tx: complex.Complex(rational.Rational) = try .initSet(allocator, x, 0);
                defer tx.deinit(allocator);
                var ty: complex.Complex(rational.Rational) = try .initSet(allocator, y, 0);
                defer ty.deinit(allocator);
                return complex.add_(allocator, o, tx, complex.neg(null, ty) catch unreachable);
            },
            else => unreachable,
        },

        else => unreachable,
    }
}
