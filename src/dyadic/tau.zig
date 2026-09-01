const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

// Michin's formula, 2π = 32 arctan(1/5) − 8 arctan(1/239)
pub fn tau(comptime Dyadic: type) Dyadic {
    comptime if (!meta.isNumeric(Dyadic) or meta.numericType(Dyadic) != .dyadic)
        @compileError("zsl.dyadic.tau: Dyadic must be a dyadic type, got \n\nDyadic = " ++ @typeName(Dyadic) ++ "\n");

    const mantissa_bits = @typeInfo(Dyadic.Mantissa).int.bits;

    comptime var term1: Dyadic.WideMantissa = (1 << (2 * mantissa_bits - 1)) / 5;
    comptime var sum1: Dyadic.WideMantissa = term1;
    comptime var divisor: Dyadic.WideMantissa = 3;
    comptime var subtract: bool = true;
    inline while (true) {
        term1 /= 25; // 5² = 25
        if (term1 == 0)
            break;

        if (subtract)
            sum1 -= term1 / divisor
        else
            sum1 += term1 / divisor;

        divisor += 2;
        subtract = !subtract;
    }

    comptime var term2: Dyadic.WideMantissa = (1 << (2 * mantissa_bits - 3)) / 239;
    comptime var sum2: Dyadic.WideMantissa = term2;
    divisor = 3;
    subtract = true;
    inline while (true) {
        term2 /= 57121; // 239² = 57121
        if (term2 == 0)
            break;

        if (subtract)
            sum2 -= term2 / divisor
        else
            sum2 += term2 / divisor;

        divisor += 2;
        subtract = !subtract;
    }

    return .{
        .mantissa = @truncate((sum1 - sum2 + (1 << (mantissa_bits - 3 - 1))) >> (mantissa_bits - 3)),
        .exponent = -numeric.cast(Dyadic.Exponent, mantissa_bits - 3),
        .positive = true,
    };
}
