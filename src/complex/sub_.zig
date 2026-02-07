const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

/// Performs in-place subtraction between two operands of complex, real,
/// rational, integer, cfloat, dyadic, float, int or bool types, where at least
/// one operand must be of complex type. The operation is performed by casting
/// both operands to complex, then subtracting them in-place.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// ## Signature
/// ```zig
/// complex.sub_(allocator: std.mem.Allocator, o: *O, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations. Must be the same allocator used to initialize `o`.
/// * `o` (`anytype`): A pointer to the output operand where the result will be
///   stored.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `complex.Error.NotWritable`: If the output operand `o` is not writable.
/// * `complex.Error.DataNotOwned`: If the output operand `o` does not own its
///   data and resizing is needed.
/// * `real.Error.NotWritable`: If the output operand `o` holds reals and its
///   real or imaginary part is not writable when it needs to be modified.
/// * `real.Error.DataNotOwned`: If the output operand `o` holds reals and its
///   real or imaginary part does not own its data and resizing is needed.
/// * `rational.Error.NotWritable`: If the output operand `o` holds rationals
///   and its real or imaginary part is not writable when it needs to be
///   modified.
/// * `rational.Error.DataNotOwned`: If the output operand `o` holds rationals
///   and its real or imaginary part does not own its data and resizing is
///   needed.
pub fn sub_(allocator: std.mem.Allocator, o: anytype, x: anytype, y: anytype) !void {
    const O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O) or
        !types.isNumeric(types.Child(O)) or types.numericType(types.Child(O)) != .complex or
        !types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.complex) or !types.numericType(Y).le(.complex) or
        (types.numericType(X) != .complex and types.numericType(Y) != .complex))
        @compileError("zml.complex.sub_: o must be a mutable pointer to a complex, and at least one of x or y must be a complex, the other must be a bool, an int, a float, a dyadic, a cfloat, an integer, a rational, a real or a complex, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (!o.flags.writable)
        return complex.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .complex => switch (comptime types.numericType(Y)) {
            .complex => return complex.add_(allocator, o, x, complex.neg(null, y) catch unreachable),
            .real => return complex.add_(allocator, o, x, complex.neg(null, y.asComplex()) catch unreachable),
            .rational => return complex.add_(allocator, o, x, complex.neg(null, y.asComplex()) catch unreachable),
            .integer => return complex.add_(allocator, o, x, complex.neg(null, y.asComplex()) catch unreachable),
            .cfloat => {
                var ty = try @import("../cfloat/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];
                ty[0].im.num.limbs = &ty[1][2];
                ty[0].im.den.limbs = &ty[1][3];

                return complex.add_(
                    allocator,
                    o,
                    x,
                    complex.neg(null, ty[0]) catch unreachable,
                );
            },
            .dyadic => {
                var ty = try @import("../dyadic/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];

                return complex.add_(
                    allocator,
                    o,
                    x,
                    complex.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var ty = try @import("../float/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1][0];
                ty[0].re.den.limbs = &ty[1][1];

                return complex.add_(
                    allocator,
                    o,
                    x,
                    complex.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var ty = @import("../int/asComplex.zig").asComplex(y);
                ty[0].re.num.limbs = &ty[1];

                return complex.add_(allocator, o, x, ty[0]);
            },
            .bool => return complex.add_(
                allocator,
                o,
                x,
                complex.neg(null, types.cast(X, y, .{}) catch unreachable) catch unreachable,
            ),
            else => unreachable,
        },
        .real => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y) catch unreachable),
        .rational => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y) catch unreachable),
        .integer => return complex.add_(allocator, o, x.asComplex(), complex.neg(null, y) catch unreachable),
        .cfloat => {
            var tx = try @import("../cfloat/asComplex.zig").asComplex(x);
            tx[0].re.num.limbs = &tx[1][0];
            tx[0].re.den.limbs = &tx[1][1];
            tx[0].im.num.limbs = &tx[1][2];
            tx[0].im.den.limbs = &tx[1][3];

            return complex.add_(
                allocator,
                o,
                tx[0],
                complex.neg(null, y) catch unreachable,
            );
        },
        .dyadic => {
            var tx = try @import("../dyadic/asComplex.zig").asComplex(x);
            tx[0].re.num.limbs = &tx[1][0];
            tx[0].re.den.limbs = &tx[1][1];

            return complex.add_(
                allocator,
                o,
                tx[0],
                complex.neg(null, y) catch unreachable,
            );
        },
        .float => {
            var tx = try @import("../float/asComplex.zig").asComplex(x);
            tx[0].re.num.limbs = &tx[1][0];
            tx[0].re.den.limbs = &tx[1][1];

            return complex.add_(
                allocator,
                o,
                tx[0],
                complex.neg(null, y) catch unreachable,
            );
        },
        .int => {
            var tx = @import("../int/asComplex.zig").asComplex(x);
            tx[0].re.num.limbs = &tx[1];

            return complex.add_(
                allocator,
                o,
                tx[0],
                complex.neg(null, y) catch unreachable,
            );
        },
        .bool => return complex.add_(
            allocator,
            o,
            complex.neg(null, types.cast(complex.Complex(rational.Rational), x, .{}) catch unreachable, y) catch unreachable,
        ),
        else => unreachable,
    }
}
