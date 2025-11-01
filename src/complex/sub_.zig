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
        .expression => @compileError("complex.sub_ not implemented for Expression yet"),
        .complex => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.sub_ not implemented for Complex + Expression yet"),
            .complex => return complex.add_(allocator, o, x, complex.neg(null, y) catch unreachable),
            .real => @compileError("complex.sub_ not implemented for Complex + Real yet"),
            .rational => return complex.add_(allocator, o, x, complex.neg(null, y.asComplex()) catch unreachable),
            .integer => return complex.add_(allocator, o, x, complex.neg(null, y.asComplex()) catch unreachable),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return complex.add_(allocator, o, x, complex.neg(null, ty[0]) catch unreachable);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return complex.add_(allocator, o, x, complex.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return complex.add_(allocator, o, x, ty[0]);
            },
            .bool => return complex.add_(allocator, o, x, complex.neg(null, types.cast(X, y, .{}) catch unreachable) catch unreachable),
        },
        .real => @compileError("complex.sub_ not implemented for Real yet"),
        .rational => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.sub_ not implemented for Rational + Expression yet"),
            .complex => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y) catch unreachable),
            .real => @compileError("complex.sub_ not implemented for Rational + Real yet"),
            .rational => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y.asComplex()) catch unreachable),
            .integer => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y.asComplex()) catch unreachable),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, ty[0]) catch unreachable);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, ty[0]) catch unreachable);
            },
            .bool => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable) catch unreachable),
        },
        .integer => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.sub_ not implemented for Integer + Expression yet"),
            .complex => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y) catch unreachable),
            .real => @compileError("complex.sub_ not implemented for Integer + Real yet"),
            .rational => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y.asComplex()) catch unreachable),
            .integer => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y.asComplex()) catch unreachable),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, ty[0]) catch unreachable);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return complex.add_(allocator, o, x.asComplex(), complex.neg(null, ty[0]) catch unreachable);
            },
            .bool => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable) catch unreachable),
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.sub_ not implemented for CFloat + Expression yet"),
            .complex => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y) catch unreachable);
            },
            .real => @compileError("complex.sub_ not implemented for CFloat + Real yet"),
            .rational => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y.asComplex()) catch unreachable);
            },
            .integer => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y.asComplex()) catch unreachable);
            },
            .cfloat => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .float => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                tx[0].im.num.limbs = &tx[1][2];
                tx[0].im.den.limbs = &tx[1][3];
                return complex.add_(allocator, o, tx[0], complex.neg(null, types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable) catch unreachable);
            },
        },
        .float => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.sub_ not implemented for Float + Expression yet"),
            .complex => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y) catch unreachable);
            },
            .real => @compileError("complex.sub_ not implemented for Float + Real yet"),
            .rational => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y.asComplex()) catch unreachable);
            },
            .integer => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y.asComplex()) catch unreachable);
            },
            .cfloat => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .float => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                var tx = try @import("../float/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1][0];
                tx[0].re.den.limbs = &tx[1][1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable) catch unreachable);
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.sub_ not implemented for Int + Expression yet"),
            .complex => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y) catch unreachable);
            },
            .real => @compileError("complex.sub_ not implemented for Int + Real yet"),
            .rational => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y.asComplex()) catch unreachable);
            },
            .integer => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, y.asComplex()) catch unreachable);
            },
            .cfloat => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .float => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                var tx = @import("../int/asComplex.zig").asComplex(x);
                tx[0].re.num.limbs = &tx[1];
                return complex.add_(allocator, o, tx[0], complex.neg(null, types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable) catch unreachable);
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .expression => @compileError("complex.sub_ not implemented for Bool + Expression yet"),
            .complex => return complex.add_(allocator, o, complex.neg(null, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, y) catch unreachable),
            .real => @compileError("complex.sub_ not implemented for Bool + Real yet"),
            .rational => return complex.add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, complex.neg(null, y.asComplex()) catch unreachable),
            .integer => return complex.add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, complex.neg(null, y.asComplex()) catch unreachable),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];
                return complex.add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, complex.neg(null, ty[0]) catch unreachable);
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                return complex.add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, complex.neg(null, ty[0]) catch unreachable);
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];
                return complex.add_(allocator, o, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, complex.neg(null, ty[0]) catch unreachable);
            },
            .bool => {
                return complex.add_(
                    allocator,
                    o,
                    types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable,
                    complex.neg(null, types.cast(complex.Complex(rational.Rational), y, .{}) catch unreachable) catch unreachable,
                );
            },
        },
    }
}
