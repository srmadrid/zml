const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");

pub fn Pow(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.pow: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

/// Performs exponentiation between two operands of int or bool types, where at
/// least one operand must be of int type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then using exponentiation by squaring for
/// efficient computation.
///
/// The operation to be performed is $x^y$.
///
/// ## Arguments
/// * `x` (`anytype`): The base value.
/// * `y` (`anytype`): The exponent value.
///
/// ## Returns
/// `Pow(@TypeOf(x), @TypeOf(y))`: The result of raising `x` to the power of
/// `y`.
///
/// ## Errors
/// * `error.NegativeExponent`: If `y` is negative.
pub inline fn pow(x: anytype, y: anytype) !Pow(@TypeOf(x), @TypeOf(y)) {
    const R: type = Pow(@TypeOf(x), @TypeOf(y));

    if (comptime R == comptime_int) {
        comptime var result: R = 1;
        comptime var base: R = types.scast(R, x);
        comptime var exponent: R = types.scast(R, y);

        if (exponent < 0)
            return error.NegativeExponent;

        inline while (exponent != 0) : (exponent >>= 1) {
            if ((exponent & 1) != 0)
                result *= base;

            base *= base;
        }

        return result;
    } else {
        var result: R = 1;
        var base: R = types.scast(R, x);
        var exponent: R = types.scast(R, y);

        if (exponent < 0)
            return error.NegativeExponent;

        while (exponent != 0) : (exponent >>= 1) {
            if ((exponent & 1) != 0)
                result = int.mul(result, base);

            base = int.mul(base, base);
        }

        return result;
    }
}
