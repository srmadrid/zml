const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs in-place multiplication between two operands of any numeric type in
/// integer precision. The operation is performed by casting both operands to
/// integer, then multiplying them in-place.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// If either `x` or `y` is of custom numeric type, that type must implement the
/// required `copyToInteger` method. The expected signature and behavior of
/// `copyToInteger` are as follows:
/// * `fn copyToInteger(self: *const @This(), allocator: std.mem.Allocator) !Integer`:
///   Initializes and returns a new integer representing the value of the
///   instance.
///
/// ## Signature
/// ```zig
/// integer.mul_(allocator: std.mem.Allocator, o: *Integer, x: X, y: Y) !void
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
pub fn mul_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.integer.mul_: x and y must be numerics, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (!o.flags.writable)
        return integer.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .custom => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ops.deinit(&ty, .{ .allocator = allocator });

                return mul_(allocator, o, tx, ty);
            },
            .complex => return mul_(allocator, o, x.re, y.re),
            .real => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ops.deinit(&ty, .{ .allocator = allocator });

                return mul_(allocator, o, tx, ty);
            },
            .rational => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ops.deinit(&ty, .{ .allocator = allocator });

                return mul_(allocator, o, tx, ty);
            },
            .integer => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });

                return mul_(allocator, o, tx, y);
            },
            .cfloat => return mul_(allocator, o, x.re, y.re),
            .dyadic => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .float => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .int => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .bool => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });

                return mul_(
                    allocator,
                    o,
                    tx,
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
        .complex => switch (comptime types.numericType(Y)) {
            .custom => return mul_(allocator, o, x.re, y),
            .complex => return mul_(allocator, o, x.re, y.re),
            .real => return mul_(allocator, o, x.re, y),
            .rational => return mul_(allocator, o, x.re, y),
            .integer => return mul_(allocator, o, x.re, y),
            .cfloat => return mul_(allocator, o, x.re, y.re),
            .dyadic => return mul_(allocator, o, x.re, y),
            .float => return mul_(allocator, o, x.re, y),
            .int => return mul_(allocator, o, x.re, y),
            .bool => return mul_(allocator, o, x.re, y),
        },
        .real => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx, ty);
            },
            .complex => return mul_(allocator, o, x, y.re),
            .real => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx, ty);
            },
            .rational => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx, ty);
            },
            .integer => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);

                return mul_(allocator, o, tx, y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .dyadic => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .float => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .int => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .bool => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);

                return mul_(
                    allocator,
                    o,
                    tx,
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
        .rational => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx, ty);
            },
            .complex => return mul_(allocator, o, x, y.re),
            .real => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx, ty);
            },
            .rational => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx, ty);
            },
            .integer => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);

                return mul_(allocator, o, tx, y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .dyadic => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .float => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .int => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx, ty[0]);
            },
            .bool => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);

                return mul_(
                    allocator,
                    o,
                    tx,
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
        .integer => switch (comptime types.numericType(Y)) {
            .custom => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, x, ty);
            },
            .complex => return mul_(allocator, o, x, y.re),
            .real => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, x, ty);
            },
            .rational => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, x, ty);
            },
            .integer => {
                if (x.size == 0 or y.size == 0)
                    return o.set(allocator, 0);

                // Aliasing check.
                var tx: Integer = if (o.limbs == x.limbs)
                    try x.copy(allocator)
                else blk: {
                    var tmp: Integer = x;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                defer tx.deinit(allocator);

                var ty: Integer = if (o.limbs == y.limbs)
                    try y.copy(allocator)
                else blk: {
                    var tmp: Integer = y;
                    tmp.flags.owns_data = false;
                    break :blk tmp;
                };
                defer ty.deinit(allocator);

                try o.reserve(allocator, tx.size + ty.size);

                var i: u32 = 0;
                while (i < tx.size + ty.size) : (i += 1) {
                    o.limbs[i] = 0;
                }

                i = 0;
                while (i < tx.size) : (i += 1) {
                    var carry: u128 = 0;
                    var j: u32 = 0;
                    while (j < ty.size) : (j += 1) {
                        const acc: u128 = types.scast(u128, o.limbs[i + j]) +
                            types.scast(u128, tx.limbs[i]) * types.scast(u128, ty.limbs[j]) + carry;
                        o.limbs[i + j] = @truncate(acc & 0xFFFFFFFF);
                        carry = acc >> 32;
                    }

                    // Propagate leftover carry.
                    var k: u32 = i + ty.size;
                    while (carry != 0) : (k += 1) {
                        const acc: u128 = types.scast(u128, o.limbs[k]) + carry;
                        o.limbs[k] = @truncate(acc & 0xFFFFFFFF);
                        carry = acc >> 32;
                    }
                }

                o.size = tx.size + ty.size;
                o.positive = tx.positive == ty.positive;
                o.truncate();
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .dyadic => {
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, x, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, x, ty[0]);
            },
            .bool => return mul_(
                allocator,
                o,
                x,
                types.cast(Integer, y, .{}) catch unreachable,
            ),
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .custom => return mul_(allocator, o, x.re, y),
            .complex => return mul_(allocator, o, x.re, y.re),
            .real => return mul_(allocator, o, x.re, y),
            .rational => return mul_(allocator, o, x.re, y),
            .integer => return mul_(allocator, o, x.re, y),
            .cfloat => return mul_(allocator, o, x.re, y.re),
            .dyadic => return mul_(allocator, o, x.re, y),
            .float => return mul_(allocator, o, x.re, y),
            .int => return mul_(allocator, o, x.re, y),
            .bool => return mul_(allocator, o, x.re, y),
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .complex => return mul_(allocator, o, x, y.re),
            .real => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .rational => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .integer => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return mul_(allocator, o, tx[0], y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .dyadic => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .float => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return mul_(
                    allocator,
                    o,
                    tx[0],
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
        .float => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .complex => return mul_(allocator, o, x, y.re),
            .real => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .rational => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .integer => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return mul_(allocator, o, tx[0], y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .dyadic => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .float => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return mul_(
                    allocator,
                    o,
                    tx[0],
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .complex => return mul_(allocator, o, x, y.re),
            .real => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .rational => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(allocator, o, tx[0], ty);
            },
            .integer => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return mul_(allocator, o, tx[0], y);
            },
            .cfloat => return mul_(allocator, o, x, y.re),
            .dyadic => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .float => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .int => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(allocator, o, tx[0], ty[0]);
            },
            .bool => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return mul_(
                    allocator,
                    o,
                    tx[0],
                    types.cast(Integer, y, .{}) catch unreachable,
                );
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .custom => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    ty,
                );
            },
            .complex => return mul_(allocator, o, x, y.re),
            .real => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    ty,
                );
            },
            .rational => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return mul_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    ty,
                );
            },
            .integer => return mul_(
                allocator,
                o,
                types.cast(Integer, x, .{}) catch unreachable,
                y,
            ),
            .cfloat => return mul_(allocator, o, x, y.re),
            .dyadic => {
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    ty[0],
                );
            },
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    ty[0],
                );
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return mul_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    ty[0],
                );
            },
            .bool => return mul_(
                allocator,
                o,
                types.cast(Integer, x, .{}) catch unreachable,
                types.cast(Integer, y, .{}) catch unreachable,
            ),
        },
    }
}
