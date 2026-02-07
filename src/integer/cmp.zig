const std = @import("std");

const types = @import("../types.zig");
const Cmp = types.Cmp;
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Compares two operands of integer, dyadic, float, int or bool types, where at
/// least one operand must be of integer type, for ordering. The operation is
/// performed by casting both operands to integer, then comparing them.
///
/// ## Signature
/// ```zig
/// integer.cmp(x: X, y: Y) Cmp
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `Cmp`: The result of the comparison.
pub fn cmp(x: anytype, y: anytype) Cmp {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.integer) or !types.numericType(Y).le(.integer) or
        types.numericType(X) == .cfloat or types.numericType(Y) == .cfloat or
        (types.numericType(X) != .integer and types.numericType(Y) != .integer))
        @compileError("zml.integer.cmp: at least one of x or y must be an integer, the other must be a bool, an int, a float, a dyadic or an integer, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                if (x.positive != y.positive)
                    return if (!x.positive and y.positive) .lt else .gt;

                if (x.size == 0 and y.size == 0) return .eq;

                var cmp_abs: Cmp = .eq;
                if (x.size != y.size) {
                    cmp_abs = if (x.size < y.size) .lt else .gt;
                } else {
                    var i: i32 = types.scast(i32, x.size - 1);
                    while (i >= 0) : (i -= 1) {
                        if (x.limbs[types.scast(u32, i)] != y.limbs[types.scast(u32, i)]) {
                            cmp_abs = if (x.limbs[types.scast(u32, i)] < y.limbs[types.scast(u32, i)]) .lt else .gt;
                            break;
                        }
                    }
                }

                return if (x.positive) cmp_abs else cmp_abs.invert();
            },
            .dyadic => {
                var ty = @import("../dyadic/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return cmp(x, ty[0]);
            },
            .float => {
                var ty = try @import("../float/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return cmp(x, ty[0]);
            },
            .int => {
                var ty = @import("../int/asInteger.zig").asInteger(y);
                ty[0].limbs = &ty[1];

                return cmp(x, ty[0]);
            },
            .bool => return cmp(x, types.cast(Integer, y, .{}) catch unreachable),
            else => unreachable,
        },
        .dyadic => {
            var tx = @import("../dyadic/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return cmp(tx[0], y);
        },
        .float => {
            var tx = try @import("../float/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return cmp(tx[0], y);
        },
        .int => {
            var tx = @import("../int/asInteger.zig").asInteger(x);
            tx[0].limbs = &tx[1];

            return cmp(tx[0], y);
        },
        .bool => return cmp(types.cast(Integer, x, .{}) catch unreachable, y),
        else => unreachable,
    }
}
