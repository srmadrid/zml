const std = @import("std");

const types = @import("../types.zig");
const Cmp = types.Cmp;
const integer = @import("../integer.zig");
const Integer = integer.Integer;

/// Compares an `Integer` with another lower or equal precision numeric type for
/// ordering.
///
/// Signature
/// ---------
/// ```zig
/// fn cmp(x: X, y: Y) Cmp
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// Returns
/// -------
/// `Cmp`:
/// The comparison result.
pub fn cmp(x: anytype, y: anytype) Cmp {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .integer and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .float) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .int) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .bool) and
        !(types.numericType(X) == .float and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .int and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .bool and types.numericType(Y) == .integer))
        @compileError("integer.cmp requires x or y to be an integer type, the other must be an integer, float, int or bool type, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

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
            .bool => {
                return cmp(x, types.cast(Integer, y, .{}) catch unreachable);
            },
            else => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .integer => {
                var tx = try @import("../float/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return cmp(tx[0], y);
            },
            else => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .integer => {
                var tx = @import("../int/asInteger.zig").asInteger(x);
                tx[0].limbs = &tx[1];
                return cmp(tx[0], y);
            },
            else => unreachable,
        },
        .bool => switch (comptime types.numericType(Y)) {
            .integer => {
                return cmp(types.cast(Integer, x, .{}) catch unreachable, y);
            },
            else => unreachable,
        },
        else => unreachable,
    }
}
