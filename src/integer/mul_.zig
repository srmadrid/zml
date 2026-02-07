const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

const check_aliasing_alloc = @import("check_aliasing_alloc.zig").check_aliasing_alloc;

/// Performs in-place multiplication between two operands of integer, cfloat,
/// dyadic, float, int or bool types, where at least one operand must be of
/// integer type. The operation is performed by casting both operands to
/// integer, then multiplying them in-place.
///
/// Aliasing between the output operand `o` and the input operands `x` or `y` is
/// allowed.
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
/// * `integer.Error.DataNotOwned`: If the output operand `o` does not own its
///   data and resizing is needed.
pub fn mul_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.integer) or !types.numericType(Y).le(.integer) or
        (types.numericType(X) != .integer and types.numericType(Y) != .integer))
        @compileError("zml.integer.mul_: at least one of x or y must be an integer, the other must be a bool, an int, a float, a dyadic, a cfloat or an integer, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (!o.flags.writable)
        return integer.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                if (x.size == 0 or y.size == 0)
                    return o.set(allocator, 0);

                // Aliasing check.
                var tx: Integer = try check_aliasing_alloc(allocator, o, x);
                defer tx.deinit(allocator);
                var ty: Integer = try check_aliasing_alloc(allocator, o, y);
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
            else => unreachable,
        },
        .cfloat => mul_(allocator, o, x.re, y),
        .dyadic => {
            var tx = try @import("../dyadic/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return mul_(allocator, o, tx[0], y);
        },
        .float => {
            var tx = try @import("../float/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return mul_(allocator, o, tx[0], y);
        },
        .int => {
            var tx = @import("../int/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return mul_(allocator, o, tx[0], y);
        },
        .bool => return mul_(
            allocator,
            o,
            types.cast(Integer, x, .{}) catch unreachable,
            y,
        ),
        else => unreachable,
    }
}
