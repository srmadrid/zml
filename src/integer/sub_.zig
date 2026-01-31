const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Performs in-place subtraction between two operands of any numeric type in
/// integer precision. The operation is performed by casting both operands to
/// integer, then subtracting them in-place.
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
pub fn sub_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.integer.sub_: x and y must be numerics, got\n\tx: " ++
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

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .complex => return sub_(allocator, o, x.re, y.re),
            .real => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ops.deinit(&ty, .{ .allocator = allocator });

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .rational => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ops.deinit(&ty, .{ .allocator = allocator });

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .integer => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, y) catch unreachable,
                );
            },
            .cfloat => return sub_(allocator, o, x.re, y.re),
            .dyadic => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer ops.deinit(&tx, .{ .allocator = allocator });

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(
                        null,
                        types.cast(Integer, y) catch unreachable,
                    ) catch unreachable,
                );
            },
        },
        .complex => switch (comptime types.numericType(Y)) {
            .custom => return sub_(allocator, o, x.re, y),
            .complex => return sub_(allocator, o, x.re, y.re),
            .real => return sub_(allocator, o, x.re, y),
            .rational => return sub_(allocator, o, x.re, y),
            .integer => return sub_(allocator, o, x.re, y),
            .cfloat => return sub_(allocator, o, x.re, y.re),
            .dyadic => return sub_(allocator, o, x.re, y),
            .float => return sub_(allocator, o, x.re, y),
            .int => return sub_(allocator, o, x.re, y),
            .bool => return sub_(allocator, o, x.re, y),
        },
        .real => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .complex => return sub_(allocator, o, x, y.re),
            .real => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .rational => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .integer => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, y) catch unreachable,
                );
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .dyadic => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(
                        null,
                        types.cast(Integer, y, .{}) catch unreachable,
                    ) catch unreachable,
                );
            },
        },
        .rational => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .complex => return sub_(allocator, o, x, y.re),
            .real => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .rational => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .integer => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, y) catch unreachable,
                );
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .dyadic => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => {
                var tx: Integer = try types.cast(Integer, x, .{ .allocator = allocator });
                defer tx.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx,
                    integer.neg(
                        null,
                        types.cast(Integer, y, .{}) catch unreachable,
                    ) catch unreachable,
                );
            },
        },
        .integer => switch (comptime types.numericType(Y)) {
            .custom => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    x,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .complex => return sub_(allocator, o, x, y.re),
            .real => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    x,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .rational => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    x,
                    integer.neg(null, ty) catch unreachable,
                );
            },
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
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .custom => return sub_(allocator, o, x.re, y),
            .complex => return sub_(allocator, o, x.re, y.re),
            .real => return sub_(allocator, o, x.re, y),
            .rational => return sub_(allocator, o, x.re, y),
            .integer => return sub_(allocator, o, x.re, y),
            .cfloat => return sub_(allocator, o, x.re, y.re),
            .dyadic => return sub_(allocator, o, x.re, y),
            .float => return sub_(allocator, o, x.re, y),
            .int => return sub_(allocator, o, x.re, y),
            .bool => return sub_(allocator, o, x.re, y),
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .complex => return sub_(allocator, o, x, y.re),
            .real => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .rational => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .integer => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, y) catch unreachable,
                );
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .dyadic => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => {
                var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(
                        null,
                        types.cast(Integer, y, .{}) catch unreachable,
                    ) catch unreachable,
                );
            },
        },
        .float => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .complex => return sub_(allocator, o, x, y.re),
            .real => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .rational => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .integer => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, y) catch unreachable,
                );
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .dyadic => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(
                        null,
                        types.cast(Integer, y, .{}) catch unreachable,
                    ) catch unreachable,
                );
            },
        },
        .int => switch (comptime types.numericType(Y)) {
            .custom => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .complex => return sub_(allocator, o, x, y.re),
            .real => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .rational => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .integer => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, y) catch unreachable,
                );
            },
            .cfloat => return sub_(allocator, o, x, y.re),
            .dyadic => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];

                return integer.add_(
                    allocator,
                    o,
                    tx[0],
                    integer.neg(
                        null,
                        types.cast(Integer, y, .{}) catch unreachable,
                    ) catch unreachable,
                );
            },
        },
        .bool => switch (comptime types.numericType(Y)) {
            .custom => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .complex => return sub_(allocator, o, x, y.re),
            .real => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .rational => {
                var ty: Integer = try types.cast(Integer, y, .{ .allocator = allocator });
                defer ty.deinit(allocator);

                return integer.add_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    integer.neg(null, ty) catch unreachable,
                );
            },
            .integer => return integer.add_(
                allocator,
                o,
                types.cast(Integer, x, .{}) catch unreachable,
                integer.neg(null, y) catch unreachable,
            ),
            .cfloat => return sub_(allocator, o, x, y.re),
            .dyadic => {
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return integer.add_(
                    allocator,
                    o,
                    types.cast(Integer, x, .{}) catch unreachable,
                    integer.neg(null, ty[0]) catch unreachable,
                );
            },
            .bool => return integer.add_(
                allocator,
                o,
                types.cast(Integer, x, .{}) catch unreachable,
                integer.neg(
                    null,
                    types.cast(Integer, y, .{}) catch unreachable,
                ) catch unreachable,
            ),
        },
    }
}
