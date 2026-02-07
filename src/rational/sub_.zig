const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

/// Performs in-place subtraction between two operands of rational, integer,
/// cfloat, dyadic, float, int or bool types, where at least one operand must be
/// of rational type. The operation is performed by casting both operands to
/// integer, then subtracting them in-place.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// ## Signature
/// ```zig
/// rational.sub_(allocator: std.mem.Allocator, o: *Rational, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations. Must be the same allocator used to initialize `o`.
/// * `o` (`*Rational`): A pointer to the output operand where the result will
///   be stored.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `rational.Error.NotWritable`: If the output operand `o` is not writable.
/// * `rational.Error.DataNotOwned`: If the output operand `o` does not own its
///   data and resizing is needed.
/// * `integer.Error.NotWritable`: If the numerator or denominator of the output
///   operand `o` is not writable when it needs to be modified.
/// * `integer.Error.DataNotOwned`: If the numerator or denominator of the
///   output operand `o` does not own its data and resizing is needed.
pub fn sub_(allocator: std.mem.Allocator, o: *Rational, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.rational) or !types.numericType(Y).le(.rational) or
        (types.numericType(X) != .rational and types.numericType(Y) != .rational))
        @compileError("zml.rational.sub_: at least one of x or y must be a rational, the other must be a bool, an int, a float, a dyadic, a cfloat, an integer or a rational, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .rational => switch (comptime types.numericType(Y)) {
            .rational => return rational.add_(
                allocator,
                o,
                x,
                rational.neg(null, y) catch unreachable,
            ),
            .integer => return rational.add_(
                allocator,
                o,
                x,
                rational.neg(null, y.asRational()) catch unreachable,
            ),
            .cfloat => return sub_(allocator, o, x, y.re),
            .dyadic => {
                var ty = @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];

                return rational.add_(
                    allocator,
                    o,
                    x,
                    rational.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var ty = try @import("../float/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1][0];
                ty[0].den.limbs = &ty[1][1];

                return rational.add_(
                    allocator,
                    o,
                    x,
                    rational.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var ty = @import("../int/asRational.zig").asRational(y);
                ty[0].num.limbs = &ty[1];

                return rational.add_(
                    allocator,
                    o,
                    x,
                    rational.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => return rational.add_(
                allocator,
                o,
                x,
                rational.neg(null, types.cast(Rational, y, .{}) catch unreachable) catch unreachable,
            ),
            else => unreachable,
        },
        .integer => return rational.add_(
            allocator,
            o,
            x.asRational(),
            rational.neg(null, y) catch unreachable,
        ),
        .cfloat => return sub_(allocator, o, x.re, y),
        .dyadic => {
            var tx = @import("../float/asRational.zig").asRational(x);
            tx[0].num.limbs = &tx[1][0];
            tx[0].den.limbs = &tx[1][1];

            return rational.add_(
                allocator,
                o,
                tx[0],
                rational.neg(null, y) catch unreachable,
            );
        },
        .float => {
            var tx = try @import("../float/asRational.zig").asRational(x);
            tx[0].num.limbs = &tx[1][0];
            tx[0].den.limbs = &tx[1][1];

            return rational.add_(
                allocator,
                o,
                tx[0],
                rational.neg(null, y) catch unreachable,
            );
        },
        .int => {
            var tx = @import("../int/asRational.zig").asRational(x);
            tx[0].num.limbs = &tx[1];

            return rational.add_(
                allocator,
                o,
                tx[0],
                rational.neg(null, y) catch unreachable,
            );
        },
        .bool => return rational.add_(
            allocator,
            o,
            types.cast(Rational, x, .{}) catch unreachable,
            rational.neg(null, y) catch unreachable,
        ),
        else => unreachable,
    }
}
