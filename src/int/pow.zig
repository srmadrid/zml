const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const Coerce = types.Coerce;

pub inline fn pow(x: anytype, y: anytype) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!(types.numericType(@TypeOf(x)) != .bool and types.numericType(@TypeOf(y)) != .int) and
        !(types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(y)) != .bool) and
        !(types.numericType(@TypeOf(x)) == .int and types.numericType(@TypeOf(y)) == .int))
        @compileError("int.pow requires at least one of x or y to be an int, the other must be a bool or an int, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (comptime C == comptime_int) {
        comptime var result: C = 1;
        comptime var base: C = types.scast(C, x);
        comptime var exponent: C = types.scast(C, y);

        if (exponent < 0)
            return error.NegativeExponent;

        inline while (exponent != 0) : (exponent >>= 1) {
            if ((exponent & 1) != 0)
                result *= base;

            base *= base;
        }

        return result;
    } else {
        var result: C = 1;
        var base: C = types.scast(C, x);
        var exponent: C = types.scast(C, y);

        if (exponent < 0)
            return error.NegativeExponent;

        while (exponent != 0) : (exponent >>= 1) {
            if ((exponent & 1) != 0)
                result *= base;

            base *= base;
        }

        return result;
    }
}
