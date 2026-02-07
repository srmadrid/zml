const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

const ops = @import("../ops.zig");
const constants = @import("../constants.zig");

const check_aliasing_alloc = @import("check_aliasing_alloc.zig").check_aliasing_alloc;

/// Performs in-place addition between two operands of integer, cfloat, dyadic,
/// float, int or bool types, where at least one operand must be of integer
/// type. The operation is performed by casting both operands to integer, then
/// adding them in-place.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
///
/// ## Signature
/// ```zig
/// integer.add_(allocator: std.mem.Allocator, o: *Integer, x: X, y: Y) !void
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
pub fn add_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.integer) or !types.numericType(Y).le(.integer) or
        (types.numericType(X) != .integer and types.numericType(Y) != .integer))
        @compileError("zml.integer.add_: at least one of x or y must be an integer, the other must be a bool, an int, a float, a dyadic, a cfloat or an integer, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (!o.flags.writable)
        return integer.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                if (x.size == 0)
                    return o.set(allocator, y);

                if (y.size == 0)
                    return o.set(allocator, x);

                // Aliasing checks
                var tx: Integer = try check_aliasing_alloc(allocator, o, x);
                defer tx.deinit(allocator);
                var ty: Integer = try check_aliasing_alloc(allocator, o, y);
                defer ty.deinit(allocator);

                if (tx.positive == ty.positive) {
                    try _addAbs_(allocator, o, tx, ty);
                    o.positive = tx.positive;
                    return;
                }

                // Reset o's size
                o.size = 0;

                const cmp_abs: types.Cmp =
                    integer.cmp(
                        integer.abs(null, tx) catch unreachable,
                        integer.abs(null, ty) catch unreachable,
                    );

                if (cmp_abs == .eq) {
                    return o.set(allocator, 0);
                } else if (cmp_abs == .gt) {
                    try _subAbs_(allocator, o, tx, ty);
                    o.positive = tx.positive;
                } else {
                    try _subAbs_(allocator, o, ty, tx);
                    o.positive = ty.positive;
                }

                return;
            },
            .cfloat => return add_(allocator, o, x, y.re),
            .dyadic => {
                var ty = try @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return add_(allocator, o, x, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return add_(allocator, o, x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return add_(allocator, o, x, ty[0]);
            },
            .bool => return add_(
                allocator,
                o,
                x,
                types.cast(Integer, y, .{}) catch unreachable,
            ),
            else => unreachable,
        },
        .cfloat => return add_(allocator, o, x.re, y),
        .dyadic => {
            var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return add_(allocator, o, tx[0], y);
        },
        .float => {
            var tx = try @import("../float/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return add_(allocator, o, tx[0], y);
        },
        .int => {
            var tx = @import("../int/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return add_(allocator, o, tx[0], y);
        },
        .bool => return add_(
            allocator,
            o,
            types.cast(Integer, x, .{}) catch unreachable,
            types.cast(Integer, y, .{}) catch unreachable,
        ),
        else => unreachable,
    }
}

fn _addAbs_(allocator: std.mem.Allocator, o: *Integer, x: Integer, y: Integer) !void {
    // No aliasing allowed.
    try o.reserve(allocator, int.max(x.size, y.size) + 1);
    var carry: u64 = 0;
    var i: u32 = 0;
    while (i < x.size or i < y.size or carry != 0) : (i += 1) {
        const xi: u64 = if (i < x.size) x.limbs[i] else 0;
        const yi: u64 = if (i < y.size) y.limbs[i] else 0;

        const s: u64 = xi + yi + carry;
        o.limbs[i] = @truncate(s & 0xFFFFFFFF);
        carry = (s >> 32);
    }

    o.size = i;
    o.truncate();

    return;
}

fn _subAbs_(allocator: std.mem.Allocator, o: *Integer, x: Integer, y: Integer) !void {
    // Assumes x >= y >= 0 in absolute value. No aliasing allowed.
    try o.reserve(allocator, x.size);
    var borrow: u64 = 0;
    var i: u32 = 0;
    while (i < x.size) : (i += 1) {
        const xi: u64 = x.limbs[i];
        const yi: u64 = if (i < y.size) y.limbs[i] else 0;

        var s: u64 = undefined;
        if (xi < yi + borrow) {
            s = xi + (@as(u64, 1) << 32) - yi - borrow;
            borrow = 1;
        } else {
            s = xi - yi - borrow;
            borrow = 0;
        }

        o.limbs[i] = @truncate(s & 0xFFFFFFFF);
    }

    o.size = i;
    o.truncate();

    return;
}
