const std = @import("std");

const zsl = @import("zsl");

test "initInt: 0" {
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(i8, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(u8, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(i16, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(u16, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(i32, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(u32, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(i64, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(u64, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(i128, 0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(u128, 0)));
}

test "initInt: 1 -> 0x8000 * 2^-15" {
    const r = zsl.Dyadic(16, 16).initValue(@as(i32, 1));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initInt: -1 -> -0x8000 * 2^-15" {
    const r = zsl.Dyadic(16, 16).initValue(@as(i32, -1));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(false, r.positive);
}

test "initInt: 3 -> 0xC000 * 2^-14" {
    const r = zsl.Dyadic(16, 16).initValue(@as(i32, 3));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initInt: 5 -> 0xA000 * 2^-13" {
    const r = zsl.Dyadic(16, 16).initValue(@as(i32, 5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xA000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -13), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initInt: -100 -> -0xC800 * 2^-9" {
    const r = zsl.Dyadic(16, 16).initValue(@as(i32, -100));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC800), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -9), r.exponent);
    try std.testing.expectEqual(false, r.positive);
}

test "initInt: 2^n -> 0x8000 * 2^(-15 + n)" {
    inline for (.{
        .{ @as(u64, 1) << 0, @as(zsl.Dyadic(16, 16).Exponent, -15) },
        .{ @as(u64, 1) << 1, @as(zsl.Dyadic(16, 16).Exponent, -14) },
        .{ @as(u64, 1) << 7, @as(zsl.Dyadic(16, 16).Exponent, -8) },
        .{ @as(u64, 1) << 15, @as(zsl.Dyadic(16, 16).Exponent, 0) },
        .{ @as(u64, 1) << 16, @as(zsl.Dyadic(16, 16).Exponent, 1) },
        .{ @as(u64, 1) << 31, @as(zsl.Dyadic(16, 16).Exponent, 16) },
        .{ @as(u64, 1) << 50, @as(zsl.Dyadic(16, 16).Exponent, 35) },
    }) |tc| {
        const r = zsl.Dyadic(16, 16).initValue(tc[0]);
        try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r.mantissa);
        try std.testing.expectEqual(tc[1], r.exponent);
        try std.testing.expectEqual(true, r.positive);
    }
}

test "initInt: 65535 -> 0xFFFF * 2^0" {
    const r = zsl.Dyadic(16, 16).initValue(@as(u32, 0xFFFF));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xFFFF), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, 0), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initInt: round down, 513 = 0x201 -> 0x80 * 2^2" {
    // 513 = 0x201, shift = 2 -> shifted = 0x80, discarded = 0b01, half = 0b10.
    // discarded < half -> no round.
    // r = 512.
    const r = zsl.Dyadic(8, 16).initValue(@as(u32, 513));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 2), r.exponent);
}

test "initInt: round up, 515 = 0x203 -> 0x81 * 2^2" {
    // 515 = 0x203, shift = 2 -> shifted = 0x80, discarded = 0b11, half = 0b10.
    // discarded > half -> round up to 0x81.
    // r = 516.
    const r = zsl.Dyadic(8, 16).initValue(@as(u32, 515));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x81), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 2), r.exponent);
}

test "initInt: round down (tie to even), 257 = 0x101 -> 0x80 * 2^1" {
    // 257 = 0x101, shift = 1 -> shifted = 0x80 (even), discarded = 0b1, half = 0b1.
    // discarded == half and shifted is even -> no round.
    // r = 256.
    const r = zsl.Dyadic(8, 16).initValue(@as(u32, 257));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 1), r.exponent);
}

test "initInt: round up (tie to even), 259 = 0x103 -> 0x82 * 2^1" {
    // 259 = 0x103, shift = 1 -> shifted = 0x81 (odd), discarded = 0b1, half = 0b1.
    // discarded == half and shifted is odd -> round up to 0x82.
    // r = 260.
    const r = zsl.Dyadic(8, 16).initValue(@as(u32, 259));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 1), r.exponent);
}

test "initInt: round (carry), 511 = 0x1FF -> 0x80 * 2^2" {
    // 511 = 0x1FF, shift = 1 -> shifted = 0xFF, discarded = 0b1, half = 0b1.
    // discarded == half -> round up to 0x100 (even).
    // overflow, shift right one -> 0x80, shift = 2.
    // r = 512.
    const r = zsl.Dyadic(8, 16).initValue(@as(u32, 511));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 2), r.exponent);
}

test "initInt: round (carry), 65535 = 0xFFFF -> 0x80 * 2^9" {
    // 65535 = 0xFFFF, shift = 8 -> shifted = 0xFF, discarded = 0xFF, half = 0x80).
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 9.
    // r = 65536.
    const r = zsl.Dyadic(8, 16).initValue(@as(u32, 0xFFFF));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 9), r.exponent);
}

test "initInt: maxInt(u8) = 0xFF -> 0xFF00000000000000 * 2^-56" {
    const r = zsl.Dyadic(64, 16).initValue(@as(u8, zsl.int.maxVal(u8)));
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0xFF00000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, -56), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initInt: maxInt(u64) = 0xFF -> 0xFFFFFFFFFFFFFFFF * 2^0" {
    const r = zsl.Dyadic(64, 16).initValue(@as(u64, zsl.int.maxVal(u64)));
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0xFFFFFFFFFFFFFFFF), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, 0), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initInt: minInt(i8) = -0x80 -> -0x80000000 * 2^-24" {
    const r = zsl.Dyadic(32, 16).initValue(@as(i8, zsl.int.minVal(i8)));
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0x80000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, -24), r.exponent);
    try std.testing.expectEqual(false, r.positive);
}

test "initInt: minInt(i32) = -0x80000000 -> -0x80000000 * 2^0" {
    const r = zsl.Dyadic(32, 16).initValue(@as(i32, zsl.int.minVal(i32)));
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0x80000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, 0), r.exponent);
    try std.testing.expectEqual(false, r.positive);
}

test "initInt: minInt(i64) = -0x8000000000000000 -> -0x80000000 * 2^32" {
    const r = zsl.Dyadic(32, 16).initValue(@as(i64, zsl.int.minVal(i64)));
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0x80000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, 32), r.exponent);
    try std.testing.expectEqual(false, r.positive);
}

test "initInt: minInt(i128) = -0x80000000000000000000000000000000 -> -0x80000000 * 2^96" {
    const r = zsl.Dyadic(32, 16).initValue(@as(i128, zsl.int.minVal(i128)));
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0x80000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, 96), r.exponent);
    try std.testing.expectEqual(false, r.positive);
}

test "initInt: positive exponent overflow -> +inf" {
    const r = zsl.Dyadic(8, 4).initValue(@as(u64, 1) << 20);
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf, r);
}

test "initInt: negative exponent overflow -> -inf" {
    const r = zsl.Dyadic(8, 4).initValue(-(@as(i64, 1) << 20));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf.neg(), r);
}

test "initInt: exponent underflow -> 0.0" {
    const r = zsl.Dyadic(32, 4).initValue(@as(i32, 1));
    try std.testing.expectEqual(zsl.Dyadic(32, 4).zero, r);
}

test "initInt: exponent underflow -> -0.0" {
    const r = zsl.Dyadic(32, 4).initValue(@as(i32, -1));
    try std.testing.expectEqual(zsl.Dyadic(32, 4).zero.neg(), r);
}

test "initInt: mantissa msb is always set" {
    inline for (.{
        @as(i64, 1),      @as(i64, 2),       @as(i64, 3),                   @as(i64, 7),
        @as(i64, 100),    @as(i64, 12345),   @as(i64, 0xDEADBEEF),          @as(i64, -1),
        @as(i64, -12345), @as(i64, -0xCAFE), @as(i64, zsl.int.maxVal(i32)), @as(i64, zsl.int.minVal(i32)),
    }) |v| {
        const r = zsl.Dyadic(16, 16).initValue(v);
        if (r.mantissa != 0)
            try std.testing.expect((r.mantissa & (@as(zsl.Dyadic(16, 16).Mantissa, 1) << 15)) != 0);
    }
}

test "initInt: unsigned values always produce positive results" {
    inline for (.{
        @as(u8, 0),  @as(u8, 1),                    @as(u8, zsl.int.maxVal(u8)),
        @as(u32, 1), @as(u32, zsl.int.maxVal(u32)), @as(u64, zsl.int.maxVal(u64)),
    }) |v| {
        try std.testing.expectEqual(true, zsl.Dyadic(16, 16).initValue(v).positive);
    }
}

test "initInt: positive * mantissa * 2^exponent exactly reconstructs the unsigned input (if no rounding)" {
    inline for (.{
        @as(u64, 1),                       @as(u64, 2),          @as(u64, 100),
        @as(u64, 12345),                   @as(u64, 0xDEADBEEF), @as(u64, 1) << 30,
        @as(u64, (@as(u64, 1) << 32) - 1), @as(u64, 1) << 50,    @as(u64, zsl.int.maxVal(u64)),
    }) |v| {
        const r = zsl.Dyadic(64, 32).initValue(v);
        const e: i64 = r.exponent;
        const m: u128 = r.mantissa;
        const reconstructed: u128 = if (e >= 0)
            m << @intCast(e)
        else
            m >> @intCast(-e);
        try std.testing.expectEqual(zsl.numeric.cast(u128, v), reconstructed);
        try std.testing.expectEqual(true, r.positive);
    }
}

test "initInt: positive * mantissa * 2^exponent exactly reconstructs the signed input (if no rounding)" {
    inline for (.{
        @as(i64, 1),      @as(i64, -1),                  @as(i64, 12345),
        @as(i64, -12345), @as(i64, zsl.int.maxVal(i64)), @as(i64, zsl.int.minVal(i64) + 1),
    }) |v| {
        const r = zsl.Dyadic(64, 32).initValue(v);
        const abs_v: u64 = @abs(v);
        const e: i64 = r.exponent;
        const m: u128 = r.mantissa;
        const reconstructed: u128 = if (e >= 0)
            m << @intCast(e)
        else
            m >> @intCast(-e);
        try std.testing.expectEqual(@as(u128, abs_v), reconstructed);
        try std.testing.expectEqual(v >= 0, r.positive);
    }
}

test "initInt: largest finite e (e = maxVal - 1)" {
    // 1 << 13 -> msb_pos = 13, shift = 6, m = 0x80, e = 6 (= maxVal - 1 for i4).
    const r = zsl.Dyadic(8, 4).initValue(@as(u64, 1) << 13);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 4).Mantissa, 0x80), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 4).Exponent, 6), r.exponent);
}

test "initInt: e exactly maxVal overflows to +inf" {
    // 1 << 14 -> msb_pos = 14, shift = 7, m = 0x80, e = 7 (= maxVal) -> inf.
    const r = zsl.Dyadic(8, 4).initValue(@as(u64, 1) << 14);
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf, r);
}

test "initInt: carry rounding pushes e to maxVal -> +inf" {
    // 0x3FFF: msb_pos = 13, shift = 6 -> shifted = 0xFF, discarded = 0x3F, half = 0x20.
    // discarded > half -> round up to 0x100.
    // overflow, shift = 7, m = 0x80, e = 7 (= maxVal) -> inf.
    const r = zsl.Dyadic(8, 4).initValue(@as(u64, 0x3FFF));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf, r);
}

test "initInt: carry rounding pushes e to maxVal -> -inf" {
    const r = zsl.Dyadic(8, 4).initValue(-@as(i64, 0x3FFF));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf.neg(), r);
}

test "initInt: 1-bit mantissa, 1 -> 1 * 2^0" {
    // msb_pos = 0 = msb_pos_result, equal branch, shift = 0. m = 1, e = 0.
    const r = zsl.Dyadic(1, 4).initValue(@as(u32, 1));
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Mantissa, 1), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Exponent, 0), r.exponent);
}

test "initInt: 1-bit mantissa, 2 -> 1 * 2^1" {
    // 2 = 0b10, msb_pos = 1, shift = 1, shifted = 1, discarded = 0. No round. m = 1, e = 1.
    const r = zsl.Dyadic(1, 4).initValue(@as(u32, 2));
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Mantissa, 1), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Exponent, 1), r.exponent);
}

test "initInt: 1-bit mantissa, 3 rounds up to 4 (1 * 2^2)" {
    // 3 = 0b11, msb_pos = 1, shift = 1, shifted = 1 (odd), discarded = 1, half = 1.
    // discarded == half and shifted is odd -> round up to 0b10.
    // overflow, shift = 2, m = 1. e = 2.
    const r = zsl.Dyadic(1, 4).initValue(@as(u32, 3));
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Mantissa, 1), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Exponent, 2), r.exponent);
}
