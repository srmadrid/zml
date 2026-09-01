const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

// Taylor expansion of eˣ at x = 1, e = 1 + 1/1 + 1/2 + 1/6 + 1/24 …
pub fn e(comptime Dyadic: type) Dyadic {
    comptime if (!meta.isNumeric(Dyadic) or meta.numericType(Dyadic) != .dyadic)
        @compileError("zsl.dyadic.e: Dyadic must be a dyadic type, got \n\nDyadic = " ++ @typeName(Dyadic) ++ "\n");

    const mantissa_bits = @typeInfo(Dyadic.Mantissa).int.bits;

    comptime var term: Dyadic.WideMantissa = 1 << (2 * mantissa_bits - 2);
    comptime var sum: Dyadic.WideMantissa = term;
    comptime var divisor: Dyadic.WideMantissa = 1;
    inline while (true) {
        term /= divisor;
        if (term == 0)
            break;

        sum += term;

        divisor += 1;
    }

    return .{
        .mantissa = @truncate((sum + (1 << (mantissa_bits - 1))) >> mantissa_bits),
        .exponent = -numeric.cast(Dyadic.Exponent, mantissa_bits - 2),
        .positive = true,
    };
}
