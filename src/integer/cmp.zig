const std = @import("std");

const types = @import("../types.zig");
const Cmp = types.Cmp;
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn cmp(x: anytype, y: anytype) Cmp {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if ((types.numericType(X) != .integer and types.numericType(X) != .int) or
        (types.numericType(X) != .integer and types.numericType(X) != .float) or
        (types.numericType(X) != .int and types.numericType(X) != .integer) or
        (types.numericType(X) != .float and types.numericType(X) != .integer))
        @compileError("integer.gt requires at least one of x or y to be an integer, the other must be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

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
            .float => {},
            .int => {},
            else => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .integer => {},
            else => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .integer => {},
            else => unreachable,
        },
        else => unreachable,
    }
}
