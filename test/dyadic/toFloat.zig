const std = @import("std");

const zsl = @import("zsl");

test "toFloat: nan" {
    try std.testing.expect(std.math.isNan(zsl.Dyadic(16, 16).nan.toFloat(f16)));
    try std.testing.expect(std.math.isNan(zsl.Dyadic(16, 16).nan.toFloat(f32)));
    try std.testing.expect(std.math.isNan(zsl.Dyadic(16, 16).nan.toFloat(f64)));
    try std.testing.expect(std.math.isNan(zsl.Dyadic(16, 16).nan.toFloat(f128)));
}

test "toFloat: +inf" {
    try std.testing.expectEqual(std.math.inf(f16), zsl.Dyadic(16, 16).inf.toFloat(f16));
    try std.testing.expectEqual(std.math.inf(f32), zsl.Dyadic(16, 16).inf.toFloat(f32));
    try std.testing.expectEqual(std.math.inf(f64), zsl.Dyadic(16, 16).inf.toFloat(f64));
    try std.testing.expectEqual(std.math.inf(f128), zsl.Dyadic(16, 16).inf.toFloat(f128));
}

test "toFloat: -inf" {
    try std.testing.expectEqual(-std.math.inf(f16), zsl.Dyadic(16, 16).inf.neg().toFloat(f16));
    try std.testing.expectEqual(-std.math.inf(f32), zsl.Dyadic(16, 16).inf.neg().toFloat(f32));
    try std.testing.expectEqual(-std.math.inf(f64), zsl.Dyadic(16, 16).inf.neg().toFloat(f64));
    try std.testing.expectEqual(-std.math.inf(f128), zsl.Dyadic(16, 16).inf.neg().toFloat(f128));
}

test "toFloat: 0.0" {
    try std.testing.expectEqual(@as(f16, 0.0), zsl.Dyadic(16, 16).zero.toFloat(f16));
    try std.testing.expectEqual(@as(f32, 0.0), zsl.Dyadic(16, 16).zero.toFloat(f32));
    try std.testing.expectEqual(@as(f64, 0.0), zsl.Dyadic(16, 16).zero.toFloat(f64));
    try std.testing.expectEqual(@as(f128, 0.0), zsl.Dyadic(16, 16).zero.toFloat(f128));
}

test "toFloat: -0.0" {
    try std.testing.expectEqual(@as(f16, -0.0), zsl.Dyadic(16, 16).zero.neg().toFloat(f16));
    try std.testing.expectEqual(@as(f32, -0.0), zsl.Dyadic(16, 16).zero.neg().toFloat(f32));
    try std.testing.expectEqual(@as(f64, -0.0), zsl.Dyadic(16, 16).zero.neg().toFloat(f64));
    try std.testing.expectEqual(@as(f128, -0.0), zsl.Dyadic(16, 16).zero.neg().toFloat(f128));
}

test "toFloat: 1.0" {
    try std.testing.expectEqual(@as(f16, 1.0), zsl.Dyadic(16, 16).one.toFloat(f16));
    try std.testing.expectEqual(@as(f32, 1.0), zsl.Dyadic(16, 16).one.toFloat(f32));
    try std.testing.expectEqual(@as(f64, 1.0), zsl.Dyadic(16, 16).one.toFloat(f64));
    try std.testing.expectEqual(@as(f128, 1.0), zsl.Dyadic(16, 16).one.toFloat(f128));
}

test "toFloat: 0x8000 * 2^(-15 + n) -> 2^n" {
    @setEvalBranchQuota(10_000);

    inline for (.{
        .{ @as(f16, 1.0), @as(zsl.Dyadic(16, 16).Exponent, -15) },
        .{ @as(f16, 2.0), @as(zsl.Dyadic(16, 16).Exponent, -14) },
        .{ @as(f16, 4.0), @as(zsl.Dyadic(16, 16).Exponent, -13) },
        .{ @as(f16, 0.5), @as(zsl.Dyadic(16, 16).Exponent, -16) },
        .{ @as(f16, 0.25), @as(zsl.Dyadic(16, 16).Exponent, -17) },
        .{ @as(f16, 1024.0), @as(zsl.Dyadic(16, 16).Exponent, -5) },
        .{ @as(f32, 1.0), @as(zsl.Dyadic(16, 16).Exponent, -15) },
        .{ @as(f32, 2.0), @as(zsl.Dyadic(16, 16).Exponent, -14) },
        .{ @as(f32, 4.0), @as(zsl.Dyadic(16, 16).Exponent, -13) },
        .{ @as(f32, 0.5), @as(zsl.Dyadic(16, 16).Exponent, -16) },
        .{ @as(f32, 0.25), @as(zsl.Dyadic(16, 16).Exponent, -17) },
        .{ @as(f32, 1024.0), @as(zsl.Dyadic(16, 16).Exponent, -5) },
        .{ @as(f64, 1.0), @as(zsl.Dyadic(16, 16).Exponent, -15) },
        .{ @as(f64, 2.0), @as(zsl.Dyadic(16, 16).Exponent, -14) },
        .{ @as(f64, 4.0), @as(zsl.Dyadic(16, 16).Exponent, -13) },
        .{ @as(f64, 0.5), @as(zsl.Dyadic(16, 16).Exponent, -16) },
        .{ @as(f64, 0.25), @as(zsl.Dyadic(16, 16).Exponent, -17) },
        .{ @as(f64, 1024.0), @as(zsl.Dyadic(16, 16).Exponent, -5) },
        .{ @as(f128, 1.0), @as(zsl.Dyadic(16, 16).Exponent, -15) },
        .{ @as(f128, 2.0), @as(zsl.Dyadic(16, 16).Exponent, -14) },
        .{ @as(f128, 4.0), @as(zsl.Dyadic(16, 16).Exponent, -13) },
        .{ @as(f128, 0.5), @as(zsl.Dyadic(16, 16).Exponent, -16) },
        .{ @as(f128, 0.25), @as(zsl.Dyadic(16, 16).Exponent, -17) },
        .{ @as(f128, 1024.0), @as(zsl.Dyadic(16, 16).Exponent, -5) },
    }) |tc| {
        const FloatType = @TypeOf(tc[0]);
        const d = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = tc[1], .positive = true };
        try std.testing.expectEqual(tc[0], d.toFloat(FloatType));
    }
}

test "toFloat: 0xC000 × 2^-15 -> 1.5" {
    const d = zsl.Dyadic(16, 16){ .mantissa = 0xC000, .exponent = -15, .positive = true };
    try std.testing.expectEqual(@as(f16, 1.5), d.toFloat(f16));
    try std.testing.expectEqual(@as(f32, 1.5), d.toFloat(f32));
    try std.testing.expectEqual(@as(f64, 1.5), d.toFloat(f64));
    try std.testing.expectEqual(@as(f128, 1.5), d.toFloat(f128));
}

test "toFloat: 0xC000 × 2^-14 -> 3.0" {
    const d = zsl.Dyadic(16, 16){ .mantissa = 0xC000, .exponent = -14, .positive = true };
    try std.testing.expectEqual(@as(f16, 3.0), d.toFloat(f16));
    try std.testing.expectEqual(@as(f32, 3.0), d.toFloat(f32));
    try std.testing.expectEqual(@as(f64, 3.0), d.toFloat(f64));
    try std.testing.expectEqual(@as(f128, 3.0), d.toFloat(f128));
}

test "toFloat: -0xA000 × 2^-14 -> -2.5" {
    const d = zsl.Dyadic(16, 16){ .mantissa = 0xA000, .exponent = -14, .positive = false };
    try std.testing.expectEqual(@as(f16, -2.5), d.toFloat(f16));
    try std.testing.expectEqual(@as(f32, -2.5), d.toFloat(f32));
    try std.testing.expectEqual(@as(f64, -2.5), d.toFloat(f64));
    try std.testing.expectEqual(@as(f128, -2.5), d.toFloat(f128));
}

test "toFloat: smallest positive normal" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = -29, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x0400))), d16.toFloat(f16)); // 2^-14

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000000, .exponent = -157, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x00800000))), d32.toFloat(f32)); // 2^-126

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000000, .exponent = -1085, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x0010000000000000))), d64.toFloat(f64)); // 2^-1022

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000000000, .exponent = -16509, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x00010000000000000000000000000000))), d128.toFloat(f128)); // 2^-16382
}

test "toFloat: largest positive subnormal" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFC0, .exponent = -30, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x03FF))), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFE00, .exponent = -158, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x007FFFFF))), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFF000, .exponent = -1086, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x000FFFFFFFFFFFFF))), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000, .exponent = -16510, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x0000FFFFFFFFFFFFFFFFFFFFFFFFFFFF))), d128.toFloat(f128));
}

test "toFloat: smallest positive subnormal" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = -39, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x0001))), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000000, .exponent = -180, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x00000001))), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000000, .exponent = -1137, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x0000000000000001))), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000000000, .exponent = -16621, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x00000000000000000000000000000001))), d128.toFloat(f128));
}

test "toFloat: round down" {
    // m = 0x800A, shift = 5, shifted = 0x0400, discarded = 0x0A, half = 0x10.
    // discarded < half -> no round.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x800A, .exponent = -15, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x000))), d16.toFloat(f16));

    // m = 0x80000040, shift = 8, shifted = 0x00800000, discarded = 0x40, half = 0x80.
    // discarded < half -> no round.
    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000040, .exponent = -31, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x000000))), d32.toFloat(f32));

    // m = 0x8000000000000200, shift = 11, shifted = 0x0010000000000000, discarded = 0x200, half = 0x400.
    // discarded < half -> no round.
    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000200, .exponent = -63, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x0000000000000))), d64.toFloat(f64));

    // m = 0x80000000000000000000000000002000, shift = 15, shifted = 0x00010000000000000000000000000000, discarded = 0x2000, half = 0x4000.
    // discarded < half -> no round.
    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000002000, .exponent = -127, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 0))), d128.toFloat(f128));
}

test "toFloat: round up" {
    // m = 0x8012, shift = 5, shifted = 0x0400, discarded = 0x12, half = 0x10.
    // discarded > half -> round up.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8012, .exponent = -15, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x001))), d16.toFloat(f16));

    // m = 0x800000C0, shift = 8, shifted = 0x00800000, discarded = 0xC0, half = 0x80.
    // discarded > half -> round up.
    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x800000C0, .exponent = -31, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x000001))), d32.toFloat(f32));

    // m = 0x8000000000000600, shift = 11, shifted = 0x0010000000000000, discarded = 0x600, half = 0x400.
    // discarded > half -> round up.
    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000600, .exponent = -63, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x0000000000001))), d64.toFloat(f64));

    // m = 0x80000000000000000000000000006000, shift = 15, shifted = 0x00010000000000000000000000000000, discarded = 0x6000, half = 0x4000.
    // discarded > half -> round up.
    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000006000, .exponent = -127, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 1))), d128.toFloat(f128));
}

test "toFloat: round down (tie to even)" {
    // m = 0x8010, shift = 5, shifted = 0x0400 (even), discarded = 0x10, half = 0x10.
    // discarded == half and shifted is even -> no round.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8010, .exponent = -15, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x000))), d16.toFloat(f16));

    // m = 0x80000080, shift = 8, shifted = 0x00800000 (even), discarded = 0x80, half = 0x80.
    // discarded == half and shifted is even -> no round.
    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000080, .exponent = -31, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x000000))), d32.toFloat(f32));

    // m = 0x8000000000000400, shift = 11, shifted = 0x0010000000000000 (even), discarded = 0x400, half = 0x400.
    // discarded == half and shifted is even -> no round.
    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000400, .exponent = -63, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x0000000000000))), d64.toFloat(f64));

    // m = 0x80000000000000000000000000004000, shift = 15, shifted = 0x00010000000000000000000000000000 (even), discarded = 0x4000, half = 0x4000.
    // discarded == half and shifted is even -> no round.
    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000004000, .exponent = -127, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 0))), d128.toFloat(f128));
}

test "toFloat: round up (tie to even)" {
    // m = 0x8030, shift = 5, shifted = 0x0401 (odd), discarded = 0x10, half = 0x10.
    // discarded == half and shifted is odd -> round up.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8030, .exponent = -15, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x002))), d16.toFloat(f16));

    // m = 0x80000180, shift = 8, shifted = 0x00800001 (odd), discarded = 0x80, half = 0x80.
    // discarded == half and shifted is odd -> round up.
    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000180, .exponent = -31, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x000002))), d32.toFloat(f32));

    // m = 0x8000000000000C00, shift = 11, shifted = 0x0010000000000001 (odd), discarded = 0x400, half = 0x400.
    // discarded == half and shifted is odd -> round up.
    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000C00, .exponent = -63, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x0000000000002))), d64.toFloat(f64));

    // m = 0x8000000000000000000000000000C000, shift = 15, shifted = 0x00010000000000000000000000000001 (odd), discarded = 0x4000, half = 0x4000.
    // discarded == half and shifted is odd -> round up.
    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x8000000000000000000000000000C000, .exponent = -127, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 2))), d128.toFloat(f128));
}

test "toFloat: round (carry)" {
    // m = 0xFFF2, shift = 5, shifted = 0x07FF, discarded = 0x12, half = 0x10.
    // discarded > half -> round up to 0x0800.
    // overflow, shift right one -> 0x0400, exponent carry.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFF2, .exponent = -15, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, (16 << 10) | 0x000))), d16.toFloat(f16));

    // m = 0xFFFFFFC0, shift = 8, shifted = 0x00FFFFFF, discarded = 0xC0, half = 0x80.
    // discarded > half -> round up to 0x01000000.
    // overflow, shift right one -> 0x00800000, exponent carry.
    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFFC0, .exponent = -31, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, (128 << 23) | 0x000000))), d32.toFloat(f32));

    // m = 0xFFFFFFFFFFFFFE00, shift = 11, shifted = 0x001FFFFFFFFFFFFF, discarded = 0x600, half = 0x400.
    // discarded > half -> round up to 0x0020000000000000.
    // overflow, shift right one -> 0x0010000000000000, exponent carry.
    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFFE00, .exponent = -63, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, (@as(u64, 1024) << 52) | 0x0000000000000))), d64.toFloat(f64));

    // m = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFE000, shift = 15, shifted = 0x0001FFFFFFFFFFFFFFFFFFFFFFFFFFFF, discarded = 0x6000, half = 0x4000.
    // discarded > half -> round up to 0x00020000000000000000000000000000.
    // overflow, shift right one -> 0x00010000000000000000000000000000, exponent carry.
    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFE000, .exponent = -127, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, (@as(u128, 16384) << 112) | 0))), d128.toFloat(f128));
}

test "toFloat: (2^n - 1) rounds up to 2^n in narrow float" {
    // m = 0xFFFF, shift = 5, shifted = 0x07FF, discarded = 0x1F, half = 0x10.
    // discarded > half -> round up to 0x0800, exponent carry.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFFF, .exponent = -5, .positive = true };
    try std.testing.expectEqual(@as(f16, 2048.0), d16.toFloat(f16));

    // m = 0xFFFFFFFF, shift = 8, shifted = 0x00FFFFFF, discarded = 0xFF, half = 0x80.
    // discarded > half -> round up to 0x01000000, exponent carry.
    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFFFF, .exponent = -8, .positive = true };
    try std.testing.expectEqual(@as(f32, 16777216.0), d32.toFloat(f32));

    // m = 0xFFFFFFFFFFFFFFFF, shift = 11, shifted = 0x001FFFFFFFFFFFFF, discarded = 0x7FF, half = 0x400.
    // discarded > half -> round up to 0x0020000000000000, exponent carry.
    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFFFFF, .exponent = -11, .positive = true };
    try std.testing.expectEqual(@as(f64, 9007199254740992.0), d64.toFloat(f64));

    // m = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF, shift = 15, shifted = 0x0001FFFFFFFFFFFFFFFFFFFFFFFFFFFF, discarded = 0x7FFF, half = 0x4000.
    // discarded > half -> round up to 0x00020000000000000000000000000000, exponent carry.
    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF, .exponent = -15, .positive = true };
    try std.testing.expectEqual(@as(f128, 10384593717069655257060992658440192.0), d128.toFloat(f128));
}

test "toFloat: positive exponent overflow -> +inf" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = 16, .positive = true };
    try std.testing.expectEqual(std.math.inf(f16), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000000, .exponent = 128, .positive = true };
    try std.testing.expectEqual(std.math.inf(f32), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000000, .exponent = 1024, .positive = true };
    try std.testing.expectEqual(std.math.inf(f64), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000000000, .exponent = 16384, .positive = true };
    try std.testing.expectEqual(std.math.inf(f128), d128.toFloat(f128));
}

test "toFloat: negative exponent overflow -> -inf" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = 16, .positive = false };
    try std.testing.expectEqual(-std.math.inf(f16), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000000, .exponent = 128, .positive = false };
    try std.testing.expectEqual(-std.math.inf(f32), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000000, .exponent = 1024, .positive = false };
    try std.testing.expectEqual(-std.math.inf(f64), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000000000, .exponent = 16384, .positive = false };
    try std.testing.expectEqual(-std.math.inf(f128), d128.toFloat(f128));
}

test "toFloat: positive exponent underflow -> 0.0" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = -40, .positive = true };
    try std.testing.expectEqual(@as(f16, 0.0), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000000, .exponent = -181, .positive = true };
    try std.testing.expectEqual(@as(f32, 0.0), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000000, .exponent = -1138, .positive = true };
    try std.testing.expectEqual(@as(f64, 0.0), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000000000, .exponent = -16622, .positive = true };
    try std.testing.expectEqual(@as(f128, 0.0), d128.toFloat(f128));
}

test "toFloat: negative exponent underflow -> -0.0" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = -40, .positive = false };
    try std.testing.expectEqual(@as(f16, -0.0), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000000, .exponent = -181, .positive = false };
    try std.testing.expectEqual(@as(f32, -0.0), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000000, .exponent = -1138, .positive = false };
    try std.testing.expectEqual(@as(f64, -0.0), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000000000, .exponent = -16622, .positive = false };
    try std.testing.expectEqual(@as(f128, -0.0), d128.toFloat(f128));
}

test "toFloat: exact roundtrip (float -> dyadic -> float)" {
    @setEvalBranchQuota(10_000);

    inline for (.{
        @as(f16, 1.0),                        @as(f16, 2.0),                        @as(f16, 3.0),
        @as(f16, 0.5),                        @as(f16, 0.25),                       @as(f16, 1.5),
        @as(f16, 2.5),                        @as(f16, 3.14),                       @as(f16, -2.71),
        @as(f16, 100.0),                      @as(f16, -1000.0),                    @as(f16, 0.001),
        @as(f16, @bitCast(@as(u16, 0x0001))), @as(f16, @bitCast(@as(u16, 0x03FF))), @as(f16, @bitCast(@as(u16, 0x0400))),
    }) |v| {
        const d = zsl.Dyadic(11, 16).initValue(v);
        try std.testing.expectEqual(v, d.toFloat(f16));
    }

    inline for (.{
        @as(f32, 1.0),                            @as(f32, 2.0),                            @as(f32, 3.0),
        @as(f32, 0.5),                            @as(f32, 0.25),                           @as(f32, 1.5),
        @as(f32, 2.5),                            @as(f32, 3.14159),                        @as(f32, -2.71828),
        @as(f32, 100.0),                          @as(f32, -1000.0),                        @as(f32, 0.001),
        @as(f32, @bitCast(@as(u32, 0x00000001))), @as(f32, @bitCast(@as(u32, 0x007FFFFF))), @as(f32, @bitCast(@as(u32, 0x00800000))),
    }) |v| {
        const d = zsl.Dyadic(24, 16).initValue(v);
        try std.testing.expectEqual(v, d.toFloat(f32));
    }

    inline for (.{
        @as(f64, 1.0),                                    @as(f64, 2.0),                                    @as(f64, 3.0),
        @as(f64, 0.5),                                    @as(f64, 0.1),                                    @as(f64, 0.2),
        @as(f64, 1.5),                                    @as(f64, 3.141592653589793),                      @as(f64, -2.718281828459045),
        @as(f64, 1e10),                                   @as(f64, 1e-10),                                  @as(f64, -1e100),
        @as(f64, @bitCast(@as(u64, 0x0000000000000001))), @as(f64, @bitCast(@as(u64, 0x000FFFFFFFFFFFFF))), @as(f64, @bitCast(@as(u64, 0x0010000000000000))),
    }) |v| {
        const d = zsl.Dyadic(53, 16).initValue(v);
        try std.testing.expectEqual(v, d.toFloat(f64));
    }

    const values: []const f128 = &.{
        @as(f128, 1.0),                                                     @as(f128, 2.0),                                                     @as(f128, 3.0),
        @as(f128, 0.5),                                                     @as(f128, 0.1),                                                     @as(f128, 0.2),
        @as(f128, 1.5),                                                     @as(f128, 3.141592653589793),                                       @as(f128, -2.718281828459045),
        @as(f128, 1e10),                                                    @as(f128, 1e-10),                                                   @as(f128, -1e100),
        @as(f128, @bitCast(@as(u128, 0x00000000000000000000000000000001))), @as(f128, @bitCast(@as(u128, 0x0000FFFFFFFFFFFFFFFFFFFFFFFFFFFF))), @as(f128, @bitCast(@as(u128, 0x00010000000000000000000000000000))),
    };
    for (values) |v| {
        const d = zsl.Dyadic(113, 16).initValue(v);
        try std.testing.expectEqual(v, d.toFloat(f128));
    }
}

test "toFloat: very narrow mantissa maps correctly" {
    const d = zsl.Dyadic(4, 16){ .mantissa = 0x8, .exponent = -3, .positive = true };
    try std.testing.expectEqual(@as(f16, 1.0), d.toFloat(f16));
    try std.testing.expectEqual(@as(f32, 1.0), d.toFloat(f32));
    try std.testing.expectEqual(@as(f64, 1.0), d.toFloat(f64));
    try std.testing.expectEqual(@as(f128, 1.0), d.toFloat(f128));
}

test "toFloat: wider target mantissa pads with zeros" {
    const d = zsl.Dyadic(8, 16){ .mantissa = 0xC0, .exponent = -7, .positive = true };
    try std.testing.expectEqual(@as(f16, 1.5), d.toFloat(f16));
    try std.testing.expectEqual(@as(f32, 1.5), d.toFloat(f32));
    try std.testing.expectEqual(@as(f64, 1.5), d.toFloat(f64));
    try std.testing.expectEqual(@as(f128, 1.5), d.toFloat(f128));
}

test "toFloat: smallest negative subnormal" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = -39, .positive = false };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x8001))), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000000, .exponent = -180, .positive = false };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x80000001))), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000000, .exponent = -1137, .positive = false };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x8000000000000001))), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000000000, .exponent = -16621, .positive = false };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x80000000000000000000000000000001))), d128.toFloat(f128));
}

test "toFloat: largest negative subnormal" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFC0, .exponent = -30, .positive = false };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x83FF))), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFE00, .exponent = -158, .positive = false };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x807FFFFF))), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFF000, .exponent = -1086, .positive = false };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x800FFFFFFFFFFFFF))), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000, .exponent = -16510, .positive = false };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x8000FFFFFFFFFFFFFFFFFFFFFFFFFFFF))), d128.toFloat(f128));
}

test "toFloat: smallest negative normal" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0x8000, .exponent = -29, .positive = false };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x8400))), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0x80000000, .exponent = -157, .positive = false };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x80800000))), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0x8000000000000000, .exponent = -1085, .positive = false };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x8010000000000000))), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0x80000000000000000000000000000000, .exponent = -16509, .positive = false };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x80010000000000000000000000000000))), d128.toFloat(f128));
}

test "toFloat: largest finite normal" {
    // self.exponent + 15 + 15 = 30, m has top 11 bits set.
    // raw_e = 15, biased_exp = 30. mantissa_shift = -5.
    // m >> 5 = 0x7FF, discarded = 0, no round. frac = 0x3FF, exp_field = 30 -> 0x7BFF.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFE0, .exponent = 0, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x7BFF))), d16.toFloat(f16));

    // self.exponent + 31 + 127 = 254 -> self.exponent = 96. m has top 24 bits set.
    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFF00, .exponent = 96, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x7F7FFFFF))), d32.toFloat(f32));

    // self.exponent + 63 + 1023 = 2046 -> self.exponent = 960. m has top 53 bits set.
    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFF800, .exponent = 960, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x7FEFFFFFFFFFFFFF))), d64.toFloat(f64));

    // self.exponent + 127 + 16383 = 32766 -> self.exponent = 16256. m has top 113 bits set.
    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF8000, .exponent = 16256, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF))), d128.toFloat(f128));
}

test "toFloat: carry rounding pushes biased_exp to overflow -> +inf" {
    // self.exponent = 0, m = 0xFFFF.
    // raw_e = 15, biased_exp = 30. mantissa_shift = -5.
    // shifted = 0x7FF, discarded = 0x1F, half = 0x10 -> round up to 0x800.
    // m == normal_threshold (1 << 11), carry: m = 0x400, biased_exp = 31 (= max_biased) -> +inf.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFFF, .exponent = 0, .positive = true };
    try std.testing.expectEqual(std.math.inf(f16), d16.toFloat(f16));

    // self.exponent = 96, m = 0xFFFFFFFF.
    // mantissa_shift = -8, shifted = 0xFFFFFF, discarded = 0xFF, half = 0x80 -> round up to 0x1000000.
    // m == normal_threshold (1 << 24), carry: biased_exp = 255 (= max_biased) -> +inf.
    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFFFF, .exponent = 96, .positive = true };
    try std.testing.expectEqual(std.math.inf(f32), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFFFFF, .exponent = 960, .positive = true };
    try std.testing.expectEqual(std.math.inf(f64), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF, .exponent = 16256, .positive = true };
    try std.testing.expectEqual(std.math.inf(f128), d128.toFloat(f128));
}

test "toFloat: carry rounding pushes biased_exp to overflow -> -inf" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFFF, .exponent = 0, .positive = false };
    try std.testing.expectEqual(-std.math.inf(f16), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFFFF, .exponent = 96, .positive = false };
    try std.testing.expectEqual(-std.math.inf(f32), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFFFFF, .exponent = 960, .positive = false };
    try std.testing.expectEqual(-std.math.inf(f64), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF, .exponent = 16256, .positive = false };
    try std.testing.expectEqual(-std.math.inf(f128), d128.toFloat(f128));
}

test "toFloat: subnormal rounds up to smallest positive normal" {
    // self.exponent = -30, m = 0xFFFF.
    // biased_exp = 0 pre-adjust, extra_shift = 1, mantissa_shift = -6.
    // shifted = 0x3FF, discarded = 0x3F, half = 0x20 -> round up to 0x400 = subnormal_threshold.
    // biased_exp = 1, m = 0 -> smallest normal.
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFFF, .exponent = -30, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x0400))), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFFFF, .exponent = -158, .positive = true };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x00800000))), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFFFFF, .exponent = -1086, .positive = true };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x0010000000000000))), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF, .exponent = -16510, .positive = true };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x00010000000000000000000000000000))), d128.toFloat(f128));
}

test "toFloat: subnormal rounds up to smallest negative normal" {
    const d16 = zsl.Dyadic(16, 16){ .mantissa = 0xFFFF, .exponent = -30, .positive = false };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x8400))), d16.toFloat(f16));

    const d32 = zsl.Dyadic(32, 16){ .mantissa = 0xFFFFFFFF, .exponent = -158, .positive = false };
    try std.testing.expectEqual(@as(f32, @bitCast(@as(u32, 0x80800000))), d32.toFloat(f32));

    const d64 = zsl.Dyadic(64, 16){ .mantissa = 0xFFFFFFFFFFFFFFFF, .exponent = -1086, .positive = false };
    try std.testing.expectEqual(@as(f64, @bitCast(@as(u64, 0x8010000000000000))), d64.toFloat(f64));

    const d128 = zsl.Dyadic(128, 16){ .mantissa = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF, .exponent = -16510, .positive = false };
    try std.testing.expectEqual(@as(f128, @bitCast(@as(u128, 0x80010000000000000000000000000000))), d128.toFloat(f128));
}

test "toFloat: right_shift_full == mantissa_bits, sticky bit one -> smallest subnormal" {
    // self.exponent = -40, mantissa = 0x8001.
    // biased_exp = -10, extra_shift = 11, mantissa_shift = -16 = -mantissa_bits.
    // right_shift_full == mantissa_bits. (mantissa & lower_mask) = 1 != 0 -> m = 1.
    // Result: biased_exp = 0, frac = 1 -> smallest positive subnormal.
    const d = zsl.Dyadic(16, 16){ .mantissa = 0x8001, .exponent = -40, .positive = true };
    try std.testing.expectEqual(@as(f16, @bitCast(@as(u16, 0x0001))), d.toFloat(f16));
}

test "toFloat: right_shift_full > mantissa_bits -> 0.0 regardless of sticky bits" {
    // self.exponent = -50, mantissa = 0xFFFF.
    // biased_exp = -20, extra_shift = 21, mantissa_shift = -26 -> right_shift_full = 26 > mantissa_bits = 16.
    // m = 0 directly (no sticky check). Result: 0.0.
    const d_pos = zsl.Dyadic(16, 16){ .mantissa = 0xFFFF, .exponent = -50, .positive = true };
    try std.testing.expectEqual(@as(f16, 0.0), d_pos.toFloat(f16));

    const d_neg = zsl.Dyadic(16, 16){ .mantissa = 0xFFFF, .exponent = -50, .positive = false };
    try std.testing.expectEqual(@as(f16, -0.0), d_neg.toFloat(f16));
}
