pub fn extractComps(f: f64) struct { positive: bool, mantissa: u53, exponent: i11 } {
    const bits: u64 = @bitCast(f);

    const positive = (bits >> 63) == 0;

    const biased_exp: u11 = @truncate((bits >> 52) & 0x7FF);
    const exponent: i11 = if (biased_exp == 0)
        -1022
    else
        @intCast(@as(i16, biased_exp) - 1023);

    const explicit_bit: u53 = if (biased_exp == 0) 0 else (1 << 52);
    const mantissa: u53 = @as(u53, @truncate(bits)) | explicit_bit;

    return .{
        .positive = positive,
        .mantissa = mantissa,
        .exponent = exponent,
    };
}

test {
    // Override test flags
    const test_all = true;

    // Individual test flags
    const test_casts = true;
    const test_arithmetic = false;

    if (test_all or test_casts) {
        _ = @import("dyadic/initBool.zig");
        _ = @import("dyadic/initInt.zig");
        _ = @import("dyadic/initFloat.zig");

        _ = @import("dyadic/toFloat.zig");
    }

    if (test_all or test_arithmetic) {
        _ = @import("dyadic/add.zig");
    }
}
