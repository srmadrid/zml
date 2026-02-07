const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs in-place subtraction between two operands of integer, cfloat,
/// dyadic, float, int or bool types, where at least one operand must be of
/// integer type. The operation is performed by casting both operands to
/// integer, then subtracting them in-place.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// ## Signature
/// ```zig
/// integer.sub_(allocator: std.mem.Allocator, o: *Integer, x: X, y: Y) !void
/// ```
///
/// ## Arguments
/// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
///   allocations. Must be the same allocator used to initialize `o`.
/// * `o` (`*Integer`): A pointer to the output operand where the result will be
///   stored.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
/// * `integer.Error.NotWritable`: If the output operand `o` is not writable.
/// * `integer.Error.DataNotOwned`: If the output operand `o` does not own its
///   data and resizing is needed.
pub fn sub_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.integer) or !types.numericType(Y).le(.integer) or
        (types.numericType(X) != .integer and types.numericType(Y) != .integer))
        @compileError("zml.integer.sub_: at least one of x or y must be an integer, the other must be a bool, an int, a float, a dyadic, a cfloat or an integer, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (!o.flags.writable)
        return integer.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => return integer.add_(
                allocator,
                o,
                x,
                integer.neg(null, y) catch unreachable,
            ),
            .cfloat => return sub_(allocator, o, x, y.re),
            .dyadic => {
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    x,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    x,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    x,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => return integer.add_(
                allocator,
                o,
                x,
                integer.neg(
                    null,
                    types.cast(Integer, y, .{}) catch unreachable,
                ) catch unreachable,
            ),
            else => unreachable,
        },
        .cfloat => return sub_(allocator, o, x.re, y),
        .dyadic => {
            var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return integer.add_(
                allocator,
                o,
                tx[0],
                integer.neg(null, y) catch unreachable,
            );
        },
        .float => {
            var tx = try @import("../float/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return integer.add_(
                allocator,
                o,
                tx[0],
                integer.neg(null, y) catch unreachable,
            );
        },
        .int => {
            var tx = @import("../int/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return integer.add_(
                allocator,
                o,
                tx[0],
                integer.neg(null, y) catch unreachable,
            );
        },
        .bool => return integer.add_(
            allocator,
            o,
            types.cast(Integer, x, .{}) catch unreachable,
            integer.neg(null, y) catch unreachable,
        ),
        else => unreachable,
    }
}
