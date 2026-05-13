const std = @import("std");

const zsl = @import("zsl");

test "initFloat: nan" {
    try std.testing.expectEqual(zsl.Dyadic(16, 16).nan, zsl.Dyadic(16, 16).initValue(std.math.nan(f16)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).nan, zsl.Dyadic(16, 16).initValue(std.math.nan(f32)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).nan, zsl.Dyadic(16, 16).initValue(std.math.nan(f64)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).nan, zsl.Dyadic(16, 16).initValue(std.math.nan(f128)));
}

test "initFloat: +inf" {
    try std.testing.expectEqual(zsl.Dyadic(16, 16).inf, zsl.Dyadic(16, 16).initValue(std.math.inf(f16)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).inf, zsl.Dyadic(16, 16).initValue(std.math.inf(f32)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).inf, zsl.Dyadic(16, 16).initValue(std.math.inf(f64)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).inf, zsl.Dyadic(16, 16).initValue(std.math.inf(f128)));
}

test "initFloat: -inf" {
    try std.testing.expectEqual(zsl.Dyadic(16, 16).inf.neg(), zsl.Dyadic(16, 16).initValue(-std.math.inf(f16)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).inf.neg(), zsl.Dyadic(16, 16).initValue(-std.math.inf(f32)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).inf.neg(), zsl.Dyadic(16, 16).initValue(-std.math.inf(f64)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).inf.neg(), zsl.Dyadic(16, 16).initValue(-std.math.inf(f128)));
}

test "initFloat: 0.0" {
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(f16, 0.0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(f32, 0.0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(f64, 0.0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero, zsl.Dyadic(16, 16).initValue(@as(f128, 0.0)));
}

test "initFloat: -0.0" {
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero.neg(), zsl.Dyadic(16, 16).initValue(@as(f16, -0.0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero.neg(), zsl.Dyadic(16, 16).initValue(@as(f32, -0.0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero.neg(), zsl.Dyadic(16, 16).initValue(@as(f64, -0.0)));
    try std.testing.expectEqual(zsl.Dyadic(16, 16).zero.neg(), zsl.Dyadic(16, 16).initValue(@as(f128, -0.0)));
}

test "initFloat: 1.0" {
    var r = zsl.Dyadic(16, 16).initValue(@as(f16, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f32, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f64, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f128, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initFloat: 2^n -> 0x8000 * 2^(-15 + n)" {
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
        const r = zsl.Dyadic(16, 16).initValue(tc[0]);
        try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r.mantissa);
        try std.testing.expectEqual(tc[1], r.exponent);
        try std.testing.expectEqual(true, r.positive);
    }
}

test "initFloat: 1.5 -> 0xC000 × 2^-15" {
    var r = zsl.Dyadic(16, 16).initValue(@as(f16, 1.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f32, 1.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f64, 1.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f128, 1.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -15), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initFloat: 3.0 -> 0xC000 × 2^-14" {
    var r = zsl.Dyadic(16, 16).initValue(@as(f16, 3.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f32, 3.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f64, 3.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f128, 3.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xC000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initFloat: -2.5 = -0xA000 × 2^-14" {
    var r = zsl.Dyadic(16, 16).initValue(@as(f16, -2.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xA000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(false, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f32, -2.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xA000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(false, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f64, -2.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xA000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(false, r.positive);

    r = zsl.Dyadic(16, 16).initValue(@as(f128, -2.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xA000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -14), r.exponent);
    try std.testing.expectEqual(false, r.positive);
}

test "initFloat: smallest positive normal" {
    const r16 = zsl.Dyadic(16, 16).initValue(@as(f16, @bitCast(@as(u16, 0x0400)))); // 2^-14
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -29), r16.exponent);
    try std.testing.expectEqual(true, r16.positive);

    const r32 = zsl.Dyadic(32, 16).initValue(@as(f32, @bitCast(@as(u32, 0x00800000)))); // 2^-126
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Mantissa, 0x80000000), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Exponent, -157), r32.exponent);
    try std.testing.expectEqual(true, r32.positive);

    const r64 = zsl.Dyadic(64, 16).initValue(@as(f64, @bitCast(@as(u64, 0x0010000000000000)))); // 2^-1022
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0x8000000000000000), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, -1085), r64.exponent);
    try std.testing.expectEqual(true, r64.positive);

    const r128 = zsl.Dyadic(128, 16).initValue(@as(f128, @bitCast(@as(u128, 0x00010000000000000000000000000000)))); // 2^-16382
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Mantissa, 0x80000000000000000000000000000000), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Exponent, -16509), r128.exponent);
    try std.testing.expectEqual(true, r128.positive);
}

test "initFloat: largest positive subnormal" {
    const r16 = zsl.Dyadic(16, 16).initValue(@as(f16, @bitCast(@as(u16, 0x03FF)))); // (2^10 - 1) * 2^-24
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0xFFC0), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -30), r16.exponent);
    try std.testing.expectEqual(true, r16.positive);

    const r32 = zsl.Dyadic(32, 16).initValue(@as(f32, @bitCast(@as(u32, 0x007FFFFF)))); // (2^23 - 1) * 2^-149
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Mantissa, 0xFFFFFE00), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Exponent, -158), r32.exponent);
    try std.testing.expectEqual(true, r32.positive);

    const r64 = zsl.Dyadic(64, 16).initValue(@as(f64, @bitCast(@as(u64, 0x000FFFFFFFFFFFFF)))); // (2^52 - 1) * 2^-1074
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0xFFFFFFFFFFFFF000), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, -1086), r64.exponent);
    try std.testing.expectEqual(true, r64.positive);

    const r128 = zsl.Dyadic(128, 16).initValue(@as(f128, @bitCast(@as(u128, 0x0000FFFFFFFFFFFFFFFFFFFFFFFFFFFF)))); // (2^112 - 1) * 2^-16494
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Mantissa, 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Exponent, -16510), r128.exponent);
    try std.testing.expectEqual(true, r128.positive);
}

test "initFloat: smallest positive subnormal" {
    const r16 = zsl.Dyadic(16, 16).initValue(@as(f16, @bitCast(@as(u16, 0x0001)))); // 2^-24
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -39), r16.exponent);
    try std.testing.expectEqual(true, r16.positive);

    const r32 = zsl.Dyadic(32, 16).initValue(@as(f32, @bitCast(@as(u32, 0x00000001)))); // 2^-149
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Mantissa, 0x80000000), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Exponent, -180), r32.exponent);
    try std.testing.expectEqual(true, r32.positive);

    const r64 = zsl.Dyadic(64, 16).initValue(@as(f64, @bitCast(@as(u64, 0x0000000000000001)))); // 2^-1074
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0x8000000000000000), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, -1137), r64.exponent);
    try std.testing.expectEqual(true, r64.positive);

    const r128 = zsl.Dyadic(128, 16).initValue(@as(f128, @bitCast(@as(u128, 0x00000000000000000000000000000001)))); // 2^-16494
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Mantissa, 0x80000000000000000000000000000000), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Exponent, -16621), r128.exponent);
    try std.testing.expectEqual(true, r128.positive);
}

test "initFloat: round down" {
    // frac = 0x00A -> raw_m = 0x40A, shift = 3, shifted = 0x81, discarded = 0x02, half = 0x04.
    // discarded < half -> no round.
    const r16 = zsl.Dyadic(8, 16).initValue(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x00A))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x81), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r16.exponent);
    try std.testing.expectEqual(true, r16.positive);

    // frac = 0x14000 -> raw_m = 0x814000, shift = 16, shifted = 0x81, discarded = 0x4000, half = 0x8000.
    // discarded < half -> no round.
    const r32 = zsl.Dyadic(8, 16).initValue(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x14000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x81), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r32.exponent);
    try std.testing.expectEqual(true, r32.positive);

    // frac = 0x2080000000000 -> raw_m = 0x10208000000000, shift = 45, shifted = 0x81, discarded = 0x80000000000 (1 << 43), half = 0x100000000000 (1 << 44).
    // discarded < half -> no round.
    const r64 = zsl.Dyadic(8, 16).initValue(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x00208000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x81), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r64.exponent);
    try std.testing.expectEqual(true, r64.positive);

    // frac = 0x0280...0000 -> raw_m = 0x1028...0000, shift = 105, shifted = 0x81, discarded = (1 << 103), half = (1 << 104).
    // discarded < half -> no round.
    const r128 = zsl.Dyadic(8, 16).initValue(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 0x0280000000000000000000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x81), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r128.exponent);
    try std.testing.expectEqual(true, r128.positive);
}

test "initFloat: round up" {
    // frac = 0x00E -> raw_m = 0x40E, shift = 3, shifted = 0x81, discarded = 0x6, half = 0x4.
    // discarded > half -> round up to 0x82.
    const r16 = zsl.Dyadic(8, 16).initValue(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x00E))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r16.exponent);
    try std.testing.expectEqual(true, r16.positive);

    // frac = 0x1A000 -> raw_m = 0x81A000, shift = 16, shifted = 0x81, discarded = 0xA000, half = 0x8000.
    // discarded > half -> round up to 0x82.
    const r32 = zsl.Dyadic(8, 16).initValue(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x1A000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r32.exponent);
    try std.testing.expectEqual(true, r32.positive);

    // frac = 0x0380... -> raw_m = 0x1038..., shift = 45, shifted = 0x81, discarded = 0x180000000000 (1.5 * 2^44), half = 0x100000000000 (2^44).
    // discarded > half -> round up to 0x82.
    const r64 = zsl.Dyadic(8, 16).initValue(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x0380000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r64.exponent);
    try std.testing.expectEqual(true, r64.positive);

    // frac = 0x0380... -> raw_m = 0x1038..., shift = 105, shifted = 0x81, discarded = 0x01800... (1.5 * 2^104), half = 0x01000... (2^104).
    // discarded > half -> round up to 0x82.
    const r128 = zsl.Dyadic(8, 16).initValue(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 0x0380000000000000000000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r128.exponent);
    try std.testing.expectEqual(true, r128.positive);
}

test "initFloat: round down (tie to even)" {
    // frac = 0x004 -> raw_m = 0x404, shift = 3, shifted = 0x80 (even), discarded = 0x4, half = 0x4.
    // discarded == half and shifted is even -> no round.
    const r16 = zsl.Dyadic(8, 16).initValue(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x004))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r16.exponent);

    // frac = 0x008000 -> raw_m = 0x808000, shift = 16, shifted = 0x80 (even), discarded = 0x8000, half = 0x8000.
    // discarded == half and shifted is even -> no round.
    const r32 = zsl.Dyadic(8, 16).initValue(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x008000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r32.exponent);

    // frac = 0x0000100000000000 -> raw_m = 0x0010100000000000, shift = 45, shifted = 0x80 (even), discarded = 0x100000000000 (2^44), half = 0x100000000000 (2^44).
    // discarded == half and shifted is even -> no round.
    const r64 = zsl.Dyadic(8, 16).initValue(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x0000100000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r64.exponent);

    // frac = 0x00000100... -> raw_m = 0x00010100..., shift = 105, shifted = 0x80 (even), discarded = 0x00000100... (2^104), half = 0x00000100... (2^104).
    // discarded == half and shifted is even -> no round.
    const r128 = zsl.Dyadic(8, 16).initValue(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 0x00000100000000000000000000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r128.exponent);
}

test "initFloat: round up (tie to even)" {
    // frac = 0x00C -> raw_m = 0x40C, shift = 3, shifted = 0x81 (odd), discarded = 0x4, half = 0x4.
    // discarded == half and shifted is odd -> round up to 0x82.
    const r16 = zsl.Dyadic(8, 16).initValue(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x00C))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r16.exponent);

    // frac = 0x018000 -> raw_m = 0x818000, shift = 16, shifted = 0x81 (odd), discarded = 0x8000, half = 0x8000.
    // discarded == half and shifted is odd -> round up to 0x82.
    const r32 = zsl.Dyadic(8, 16).initValue(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x018000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r32.exponent);

    // frac = 0x0000300000000000 -> raw_m = 0x0010003000000000, shift = 45, shifted = 0x81 (odd), discarded = 0x0000100000000000 (2^44), half = 0x0000100000000000 (2^44).
    // discarded == half and shifted is odd -> round up to 0x82.
    const r64 = zsl.Dyadic(8, 16).initValue(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x0000300000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r64.exponent);

    // frac = 0x00000300... -> raw_m = 0x00010300..., shift = 105, shifted = 0x81 (odd), discarded = 0x00000100... (2^104), half = 0x00000100... (2^104).
    // discarded == half and shifted is odd -> round up to 0x82.
    const r128 = zsl.Dyadic(8, 16).initValue(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 0x00000300000000000000000000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x82), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -7), r128.exponent);
}

test "initFloat: round (carry)" {
    // frac = 0x03FE -> raw_m = 0x07FE (bits 10..1 set), shift = 3, shifted = 0xFF, discarded = 0x6, half = 0x4.
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 4.
    const r16 = zsl.Dyadic(8, 16).initValue(@as(f16, @bitCast(@as(u16, (15 << 10) | 0x03FE))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -6), r16.exponent);

    // frac = 0x7FC000 -> raw_m = 0xFFC000 (bits 23..14 set), shift = 16, shifted = 0xFF, discarded = 0xC000, half = 0x8000.
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 17.
    const r32 = zsl.Dyadic(8, 16).initValue(@as(f32, @bitCast(@as(u32, (127 << 23) | 0x7FC000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -6), r32.exponent);

    // frac = 0x000FFF8000000000 -> raw_m = 0x001FFF8000000000 (bits 52..43 set), shift = 45, shifted = 0xFF, discarded = 0x0000180000000000, half = 0x0000100000000000.
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 46.
    const r64 = zsl.Dyadic(8, 16).initValue(@as(f64, @bitCast(@as(u64, (@as(u64, 1023) << 52) | 0x000FFF8000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -6), r64.exponent);

    // frac = 0x0000FF80... -> raw_m = 0x0001FF80... (bits 112..103 set) shift = 105, shifted = 0xFF, discarded = 0x00000180..., half = 0x00000100....
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 106.
    const r128 = zsl.Dyadic(8, 16).initValue(@as(f128, @bitCast(@as(u128, (@as(u128, 16383) << 112) | 0x0000FF80000000000000000000000000))));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, -6), r128.exponent);
}

test "initFloat: (2^n - 1) rounds up to 2^n in narrow dyadic" {
    // frac = 0x03FF -> raw_m = 0x07FF (bits 10..0 set), shift = 3, shifted = 0xFF, discarded = 0x7, half = 0x4.
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 4.
    const r16 = zsl.Dyadic(8, 16).initValue(@as(f16, @bitCast(@as(u16, (25 << 10) | 0x03FF)))); // 2^11 - 1
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 4), r16.exponent);

    // frac = 0x7FFFFF -> raw_m = 0xFFFFFF (bits 23..0 set), shift = 16, shifted = 0xFF, discarded = 0xFFFF, half = 0x8000.
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 17.
    const r32 = zsl.Dyadic(8, 16).initValue(@as(f32, @bitCast(@as(u32, (150 << 23) | 0x007FFFFF)))); // 2^24 - 1
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 17), r32.exponent);

    // frac = 0x000FFFFFFFFFFFFF -> raw_m = 0x001FFFFFFFFFFFFF (bits 52..0 set), shift = 45, shifted = 0xFF, discarded = 0x00001FFFFFFFFFFF, half = 0x0000100000000000.
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 46.
    const r64 = zsl.Dyadic(8, 16).initValue(@as(f64, @bitCast(@as(u64, (@as(u64, 1075) << 52) | 0x000FFFFFFFFFFFFF)))); // 2^53 - 1
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 46), r64.exponent);

    // frac = 0x0000FFFF... -> raw_m = 0x0001FFFF... (bits 112..0 set), shift = 105, shifted = 0xFF, discarded = 0x000001FFFF..., half = 0x00000100....
    // discarded > half -> round up to 0x100.
    // overflow, shift right one -> 0x80, shift = 106.
    const r128 = zsl.Dyadic(8, 16).initValue(@as(f128, @bitCast(@as(u128, (@as(u128, 16495) << 112) | 0x0000FFFFFFFFFFFFFFFFFFFFFFFFFFFF)))); // 2^113 - 1
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Mantissa, 0x80), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 16).Exponent, 106), r128.exponent);
}

test "initFloat: positive exponent overflow -> +inf" {
    const r16 = zsl.Dyadic(8, 4).initValue(@as(f16, @as(f16, 1) * @as(f16, 1 << 14)));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf, r16);

    const r32 = zsl.Dyadic(16, 4).initValue(@as(f32, @as(f32, 1) * @as(f32, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).inf, r32);

    const r64 = zsl.Dyadic(16, 4).initValue(@as(f64, @as(f64, 1) * @as(f64, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).inf, r64);

    const r128 = zsl.Dyadic(16, 4).initValue(@as(f128, @as(f128, 1) * @as(f128, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).inf, r128);
}

test "initFloat: negative exponent overflow -> -inf" {
    const r16 = zsl.Dyadic(8, 4).initValue(@as(f16, @as(f16, -1) * @as(f16, 1 << 14)));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf.neg(), r16);

    const r32 = zsl.Dyadic(16, 4).initValue(@as(f32, @as(f32, -1) * @as(f32, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).inf.neg(), r32);

    const r64 = zsl.Dyadic(16, 4).initValue(@as(f64, @as(f64, -1) * @as(f64, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).inf.neg(), r64);

    const r128 = zsl.Dyadic(16, 4).initValue(@as(f128, @as(f128, -1) * @as(f128, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).inf.neg(), r128);
}

test "initFloat: positive exponent underflow -> 0.0" {
    const r16 = zsl.Dyadic(8, 4).initValue(@as(f16, @as(f16, 1) / @as(f16, 1 << 14)));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).zero, r16);

    const r32 = zsl.Dyadic(16, 4).initValue(@as(f32, @as(f32, 1) / @as(f32, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero, r32);

    const r64 = zsl.Dyadic(16, 4).initValue(@as(f64, @as(f64, 1) / @as(f64, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero, r64);

    const r128 = zsl.Dyadic(16, 4).initValue(@as(f128, @as(f128, 1) / @as(f128, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero, r128);
}

test "initFloat: negative exponent underflow -> -0.0" {
    const r16 = zsl.Dyadic(8, 4).initValue(@as(f16, @as(f16, -1) / @as(f16, 1 << 14)));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).zero.neg(), r16);

    const r32 = zsl.Dyadic(16, 4).initValue(@as(f32, @as(f32, -1) / @as(f32, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero.neg(), r32);

    const r64 = zsl.Dyadic(16, 4).initValue(@as(f64, @as(f64, -1) / @as(f64, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero.neg(), r64);

    const r128 = zsl.Dyadic(16, 4).initValue(@as(f128, @as(f128, -1) / @as(f128, 1 << 22)));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero.neg(), r128);
}

test "initFloat: tiny subnormal underflows narrow exponent type" {
    const r16 = zsl.Dyadic(16, 4).initValue(@as(f16, @bitCast(@as(u16, 0x0001))));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero, r16);

    const r32 = zsl.Dyadic(16, 4).initValue(@as(f32, @bitCast(@as(u32, 0x00000001))));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero, r32);

    const r64 = zsl.Dyadic(16, 4).initValue(@as(f64, @bitCast(@as(u64, 0x0000000000000001))));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero, r64);

    const r128 = zsl.Dyadic(16, 4).initValue(@as(f128, @bitCast(@as(u128, 0x00000000000000000000000000000001))));
    try std.testing.expectEqual(zsl.Dyadic(16, 4).zero, r128);
}

test "initFloat: mantissa msb is always set" {
    inline for (.{
        @as(f16, 1.0),                                                      @as(f16, -2.5),                       @as(f16, 100.0),
        @as(f16, -0.001),                                                   @as(f16, @bitCast(@as(u16, 0x0001))), @as(f32, 12345.0),
        @as(f32, -0.12345),                                                 @as(f32, 1e20),                       @as(f32, -1e-20),
        @as(f32, @bitCast(@as(u32, 0x00000001))),                           @as(f64, 2.718281828),                @as(f64, -3735928559.0),
        @as(f64, 1e100),                                                    @as(f64, -1e-100),                    @as(f64, @bitCast(@as(u64, 0x0000000000000001))),
        @as(f128, 1.6180339887),                                            @as(f128, 1e4000),                    @as(f128, -1e-4000),
        @as(f128, @bitCast(@as(u128, 0x00000000000000000000000000000001))),
    }) |v| {
        const r = zsl.Dyadic(16, 16).initValue(v);
        if (r.mantissa != 0)
            try std.testing.expect((r.mantissa & (@as(zsl.Dyadic(16, 16).Mantissa, 1) << 15)) != 0);
    }
}

test "initFloat: positive * mantissa * 2^exponent exactly reconstructs the input" {
    @setEvalBranchQuota(10_000);

    inline for (.{
        @as(f16, 1.0),                        @as(f16, 2.0),                        @as(f16, 3.0),
        @as(f16, 0.5),                        @as(f16, 0.25),                       @as(f16, 1.5),
        @as(f16, 2.5),                        @as(f16, 3.14),                       @as(f16, -2.71),
        @as(f16, 100.0),                      @as(f16, -1000.0),                    @as(f16, 0.001),
        @as(f16, @bitCast(@as(u16, 0x0001))), @as(f16, @bitCast(@as(u16, 0x03FF))), @as(f16, @bitCast(@as(u16, 0x0400))),
    }) |v| {
        const r = zsl.Dyadic(11, 16).initValue(v);
        const abs_v: f64 = zsl.numeric.cast(f64, zsl.float.abs(v));
        const m_f: f64 = zsl.numeric.cast(f64, r.mantissa);
        const reconstructed: f64 = zsl.float.ldexp(m_f, r.exponent);
        try std.testing.expectEqual(abs_v, reconstructed);
        try std.testing.expectEqual(v >= 0, r.positive);
    }

    inline for (.{
        @as(f32, 1.0),                            @as(f32, 2.0),                            @as(f32, 3.0),
        @as(f32, 0.5),                            @as(f32, 0.25),                           @as(f32, 1.5),
        @as(f32, 2.5),                            @as(f32, 3.14159),                        @as(f32, -2.71828),
        @as(f32, 100.0),                          @as(f32, -1000.0),                        @as(f32, 0.001),
        @as(f32, @bitCast(@as(u32, 0x00000001))), @as(f32, @bitCast(@as(u32, 0x007FFFFF))), @as(f32, @bitCast(@as(u32, 0x00800000))),
    }) |v| {
        const r = zsl.Dyadic(24, 16).initValue(v);
        const abs_v: f64 = zsl.numeric.cast(f64, zsl.float.abs(v));
        const m_f: f64 = zsl.numeric.cast(f64, r.mantissa);
        const reconstructed: f64 = zsl.float.ldexp(m_f, r.exponent);
        try std.testing.expectEqual(abs_v, reconstructed);
        try std.testing.expectEqual(v >= 0, r.positive);
    }

    inline for (.{
        @as(f64, 1.0),                                    @as(f64, 2.0),                                    @as(f64, 3.0),
        @as(f64, 0.5),                                    @as(f64, 0.1),                                    @as(f64, 0.2),
        @as(f64, 1.5),                                    @as(f64, 3.141592653589793),                      @as(f64, -2.718281828459045),
        @as(f64, 1e10),                                   @as(f64, 1e-10),                                  @as(f64, -1e100),
        @as(f64, @bitCast(@as(u64, 0x0000000000000001))), @as(f64, @bitCast(@as(u64, 0x000FFFFFFFFFFFFF))), @as(f64, @bitCast(@as(u64, 0x0010000000000000))),
    }) |v| {
        const r = zsl.Dyadic(53, 16).initValue(v);
        const abs_v: f64 = zsl.float.abs(v);
        const m_f: f64 = zsl.numeric.cast(f64, r.mantissa);
        const reconstructed: f64 = zsl.float.ldexp(m_f, r.exponent);
        try std.testing.expectEqual(abs_v, reconstructed);
        try std.testing.expectEqual(v >= 0, r.positive);
    }

    const values: []const f128 = &.{
        @as(f128, 1.0),                                                     @as(f128, 2.0),
        @as(f128, 3.0),                                                     @as(f128, 0.5),
        @as(f128, 0.1),                                                     @as(f128, 0.2),
        @as(f128, 1.5),                                                     @as(f128, 3.141592653589793),
        @as(f128, -2.718281828459045),                                      @as(f128, 1e10),
        @as(f128, 1e-10),                                                   @as(f128, -1e100),
        @as(f128, @bitCast(@as(u128, 0x00000000000000000000000000000001))), @as(f128, @bitCast(@as(u128, 0x0000FFFFFFFFFFFFFFFFFFFFFFFFFFFF))),
        @as(f128, @bitCast(@as(u128, 0x00010000000000000000000000000000))),
    };
    for (values) |v| {
        const r = zsl.Dyadic(113, 16).initValue(v);
        const abs_v: f128 = zsl.float.abs(v);
        const m_f: f128 = zsl.numeric.cast(f128, r.mantissa);
        const reconstructed: f128 = zsl.float.ldexp(m_f, r.exponent);
        try std.testing.expectEqual(abs_v, reconstructed);
        try std.testing.expectEqual(v >= 0, r.positive);
    }
}

test "initFloat: very narrow mantissa rounds correctly" {
    const r16 = zsl.Dyadic(4, 16).initValue(@as(f16, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Mantissa, 0x8), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Exponent, -3), r16.exponent);

    const r32 = zsl.Dyadic(4, 16).initValue(@as(f32, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Mantissa, 0x8), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Exponent, -3), r32.exponent);

    const r64 = zsl.Dyadic(4, 16).initValue(@as(f64, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Mantissa, 0x8), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Exponent, -3), r64.exponent);

    const r128 = zsl.Dyadic(4, 16).initValue(@as(f128, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Mantissa, 0x8), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Exponent, -3), r128.exponent);
}

test "initFloat: wider target mantissa pads with zeros" {
    const r16 = zsl.Dyadic(32, 16).initValue(@as(f16, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Mantissa, 0x80000000), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Exponent, -31), r16.exponent);

    const r32 = zsl.Dyadic(64, 16).initValue(@as(f32, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0x8000000000000000), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, -63), r32.exponent);

    const r64 = zsl.Dyadic(128, 16).initValue(@as(f64, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Mantissa, 0x80000000000000000000000000000000), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Exponent, -127), r64.exponent);

    const r128 = zsl.Dyadic(128, 16).initValue(@as(f128, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Mantissa, 0x80000000000000000000000000000000), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Exponent, -127), r128.exponent);
}

test "initFloat: largest finite e" {
    // 8192 = 2^13 (f32). raw_m = 1 << 23, raw_e = 140 - 127 - 23 = -10.
    // msb_pos = 23, msb_pos_result = 7, shift = 16. m = 0x80, e = -10 + 16 = 6 (= maxVal - 1 for i4).
    const r = zsl.Dyadic(8, 4).initValue(@as(f32, 8192.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 4).Mantissa, 0x80), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(8, 4).Exponent, 6), r.exponent);
}

test "initFloat: e exactly maxVal overflows to +inf" {
    // 16384 = 2^14 (f32). raw_m = 1 << 23, raw_e = -9.
    // msb_pos = 23, shift = 16, m = 0x80, e = 7 (= maxVal) -> inf.
    const r = zsl.Dyadic(8, 4).initValue(@as(f32, 16384.0));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf, r);
}

test "initFloat: carry rounding pushes e to maxVal -> +inf" {
    // 16383.0 (f32) = (2^14 - 1). biased_exp = 140, frac = 0x7FFE00, raw_m = 0xFFFC00, raw_e = -10.
    // msb_pos = 23, shift = 16, shifted = 0xFF, discarded = 0xFC00, half = 0x8000.
    // discarded > half -> round up to 0x100.
    // overflow, shift = 17, m = 0x80, e = -10 + 17 = 7 (= maxVal) -> inf.
    const r = zsl.Dyadic(8, 4).initValue(@as(f32, 16383.0));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf, r);
}

test "initFloat: carry rounding pushes e to maxVal -> -inf" {
    const r = zsl.Dyadic(8, 4).initValue(@as(f32, -16383.0));
    try std.testing.expectEqual(zsl.Dyadic(8, 4).inf.neg(), r);
}

test "initFloat: negative subnormals" {
    const r16 = zsl.Dyadic(16, 16).initValue(-@as(f16, @bitCast(@as(u16, 0x0001))));
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Mantissa, 0x8000), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(16, 16).Exponent, -39), r16.exponent);
    try std.testing.expectEqual(false, r16.positive);

    const r32 = zsl.Dyadic(32, 16).initValue(-@as(f32, @bitCast(@as(u32, 0x00000001))));
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Mantissa, 0x80000000), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(32, 16).Exponent, -180), r32.exponent);
    try std.testing.expectEqual(false, r32.positive);

    const r64 = zsl.Dyadic(64, 16).initValue(-@as(f64, @bitCast(@as(u64, 0x0000000000000001))));
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Mantissa, 0x8000000000000000), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(64, 16).Exponent, -1137), r64.exponent);
    try std.testing.expectEqual(false, r64.positive);

    const r128 = zsl.Dyadic(128, 16).initValue(-@as(f128, @bitCast(@as(u128, 0x00000000000000000000000000000001))));
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Mantissa, 0x80000000000000000000000000000000), r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(128, 16).Exponent, -16621), r128.exponent);
    try std.testing.expectEqual(false, r128.positive);
}

test "initFloat: subnormal with rounding" {
    // frac = 0x03FF, biased_exp = 0 -> raw_m = 0x03FF, raw_e = -24.
    // msb_pos = 9, msb_pos_result = 3, shift = 6.
    // shifted = 0xF, discarded = 0x3F, half = 0x20.
    // discarded > half -> round up to 0x10.
    // overflow, shift = 7, m = 0x8. e = -24 + 7 = -17.
    const r = zsl.Dyadic(4, 16).initValue(@as(f16, @bitCast(@as(u16, 0x03FF))));
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Mantissa, 0x8), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(4, 16).Exponent, -17), r.exponent);
    try std.testing.expectEqual(true, r.positive);
}

test "initFloat: exact-equal msb branch" {
    // biased_exp = 15, frac = 0. raw_m = 0x400 (bit 10 set). raw_e = -10.
    // msb_pos = 10 = msb_pos_result. Equal branch. m = 0x400, e = -10.
    const r16 = zsl.Dyadic(11, 16).initValue(@as(f16, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(11, 16).Mantissa, 0x400), r16.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(11, 16).Exponent, -10), r16.exponent);
    try std.testing.expectEqual(true, r16.positive);

    // raw_m = 0x00800000, raw_e = -23.
    // msb_pos = 23 = msb_pos_result. Equal branch. m = 0x00800000, e = -23.
    const r32 = zsl.Dyadic(24, 16).initValue(@as(f32, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(24, 16).Mantissa, 0x00800000), r32.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(24, 16).Exponent, -23), r32.exponent);
    try std.testing.expectEqual(true, r32.positive);

    // raw_m = 0x0010000000000000, raw_e = -52.
    // msb_pos = 52 = msb_pos_result. Equal branch.
    const r64 = zsl.Dyadic(53, 16).initValue(@as(f64, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 16).Mantissa, 0x0010000000000000), r64.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 16).Exponent, -52), r64.exponent);
    try std.testing.expectEqual(true, r64.positive);

    // raw_m = 1 << 112, raw_e = -112.
    // msb_pos = 112 = msb_pos_result. Equal branch.
    const r128 = zsl.Dyadic(113, 16).initValue(@as(f128, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(113, 16).Mantissa, 1) << 112, r128.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(113, 16).Exponent, -112), r128.exponent);
    try std.testing.expectEqual(true, r128.positive);
}

test "initFloat: 1-bit mantissa, 1.0 -> 1 * 2^0" {
    // raw_m = 0x00800000, raw_e = -23.
    // msb_pos = 23, msb_pos_result = 0, shift = 23.
    // shifted = 1, discarded = 0. No round. m = 1, e = 0.
    const r = zsl.Dyadic(1, 4).initValue(@as(f32, 1.0));
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Mantissa, 1), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Exponent, 0), r.exponent);
}

test "initFloat: 1-bit mantissa, 1.5 rounds up to 2.0" {
    // raw_m = 0xC00000, raw_e = -23.
    // msb_pos = 23, shift = 23, shifted = 1 (odd), discarded = 0x400000, half = 0x400000.
    // discarded == half and shifted is odd -> round up to 0b10.
    // overflow, shift = 24, m = 1. e = 1.
    const r = zsl.Dyadic(1, 4).initValue(@as(f32, 1.5));
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Mantissa, 1), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(1, 4).Exponent, 1), r.exponent);
}
