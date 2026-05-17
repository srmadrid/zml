const std = @import("std");

const zsl = @import("zsl");

// Verifies that y is the correctly-rounded square root of x in Dyadic D.
//
// For finite positive x = m * 2^e and finite positive normalized
// y = m' * 2^e':
//   y is correctly rounded  iff  (2m' - 1)² ≤ m * 2^(e - 2e' + 2) ≤ (2m' + 1)²
//
// Derivation: y is correctly rounded  iff  sqrt(x) ∈ [y - 0.5 ulp, y + 0.5 ulp]
// iff  x ∈ [(y - 0.5 ulp)², (y + 0.5 ulp)²]. Substituting y = m' * 2^e' and
// ulp = 2^e' and cancelling the common 2^(2e' - 2) factor gives the
// inequality above. Ties (sqrt landing exactly on y ± 0.5 ulp) cannot occur
// for the integer n = m * 2^shift the algorithm operates on, since they
// would require shifted_m = (2m' ± 1)² to be odd while
// shifted_m = m * 2^(≥1) is always even.
fn isCorrectlyRoundedSqrt(comptime D: type, x: D, y: D) bool {
    if (x.isNan())
        return y.isNan();

    if (!x.positive) {
        if (x.isZero())
            return y.isZero() and !y.positive;

        return y.isNan();
    }

    if (x.isInf())
        return y.isInf() and y.positive;

    if (x.isZero())
        return y.isZero() and y.positive;

    if (y.isNan() or y.isInf() or y.isZero() or !y.positive)
        return false;

    const N: u16 = @typeInfo(D.Mantissa).int.bits;
    const Wide = @Int(.unsigned, 2 * N + 4);

    const shift_signed: i64 =
        @as(i64, x.exponent) - 2 * @as(i64, y.exponent) + 2;

    // For a valid result of the algorithm shift ∈ {N-1, N, N+1, N+2}; anything
    // else is a sign the algorithm produced an unreachable representation.
    if (shift_signed < 0 or shift_signed > N + 2)
        return false;
    const shift: u16 = @intCast(shift_signed);

    const shifted_m: Wide = @as(Wide, x.mantissa) << @intCast(shift);

    const m_prime: Wide = @as(Wide, y.mantissa);
    const lower: Wide = (2 * m_prime - 1) * (2 * m_prime - 1);
    const upper: Wide = (2 * m_prime + 1) * (2 * m_prime + 1);

    return lower <= shifted_m and shifted_m <= upper;
}

fn expectCorrectlyRoundedSqrt(comptime D: type, x: D, y: D) !void {
    if (!isCorrectlyRoundedSqrt(D, x, y)) {
        std.debug.print(
            "\nsqrt result not correctly rounded:\n  x:        {}\n  computed: {}\n",
            .{ x, y },
        );
        return error.IncorrectSqrt;
    }
}

test "sqrt: nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).nan,
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).nan),
    );
}

test "sqrt: +inf -> +inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).inf,
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).inf),
    );
}

test "sqrt: -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).nan,
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "sqrt: +0.0 -> +0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).zero,
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).zero),
    );
}

test "sqrt: -0.0 -> -0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).zero.neg(),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "sqrt: -1.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).nan,
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).one.neg()),
    );
}

test "sqrt: -4.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).nan,
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, -4.0))),
    );
}

test "sqrt: large negative -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).nan,
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, -1e10))),
    );
}

test "sqrt: smallest negative finite -> nan" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent) + 1, .positive = false };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).nan,
        zsl.dyadic.sqrt(x),
    );
}

test "sqrt: 1.0 -> 1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).one,
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).one),
    );
}

test "sqrt: 4.0 -> 2.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, 4.0))),
    );
}

test "sqrt: 9.0 -> 3.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 3.0)),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, 9.0))),
    );
}

test "sqrt: 16.0 -> 4.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 4.0)),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, 16.0))),
    );
}

test "sqrt: 25.0 -> 5.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 5.0)),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, 25.0))),
    );
}

test "sqrt: 49.0 -> 7.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 7.0)),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, 49.0))),
    );
}

test "sqrt: 256.0 -> 16.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 16.0)),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, 256.0))),
    );
}

test "sqrt: 0.25 -> 0.5" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.5)),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, 0.25))),
    );
}

test "sqrt: 0.0625 -> 0.25" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.25)),
        zsl.dyadic.sqrt(zsl.Dyadic(53, 11).initValue(@as(f64, 0.0625))),
    );
}

test "sqrt: 2^104 -> 2^52" {
    // x = 2^52 * 2^52. y = 2^52 * 2^0.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 52, .positive = true };
    const expected = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    try std.testing.expectEqual(expected, zsl.dyadic.sqrt(x));
}

test "sqrt: 2^-104 -> 2^-52" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -156, .positive = true };
    const expected = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -104, .positive = true };
    try std.testing.expectEqual(expected, zsl.dyadic.sqrt(x));
}

test "sqrt: e even, shift = N - 1 (parity mismatch, perfect square)" {
    // m = 2^52, e = 0 (even). N = 53 (odd). Parity differ -> shift = 52.
    // n = m << 52 = 2^104. q = isqrt(2^104) = 2^52. rem = 0.
    // m' = 2^52, e' = (0 - 52) / 2 = -26. value = 2^26.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -26), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
}

test "sqrt: e odd, shift = N (parity match, perfect square)" {
    // m = 9 * 2^49 = 0x12000000000000, e = -49 (odd). N = 53 (odd). Same parity -> shift = 53.
    // n = m << 53 = 9 * 2^102. q = isqrt(9 * 2^102) = 3 * 2^51 = 0x18000000000000. rem = 0.
    // m' = 0x18000000000000, e' = (-49 - 53) / 2 = -51. value = 3.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x12000000000000, .exponent = -49, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x18000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -51), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
}

test "sqrt: round down boundary (rem == q)" {
    // m = 2^52 + 1, e = 0 (even). shift = 52.
    // n = (2^52 + 1) << 52 = 2^104 + 2^52. q = 2^52 (next square is 2^104 + 2^53 + 1).
    // rem = n - q² = 2^52 = q exactly. Algorithm rounds down.
    // m' = 2^52 = 0x10000000000000, e' = -26.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = 0, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -26), r.exponent);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
}

test "sqrt: round up boundary (rem > q)" {
    // m = 2^52 + 2, e = 0. shift = 52.
    // n = (2^52 + 2) << 52 = 2^104 + 2^53. q = 2^52, rem = 2^53 = 2q > q. Round UP.
    // m' = 2^52 + 1 = 0x10000000000001, e' = -26.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000002, .exponent = 0, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000001), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -26), r.exponent);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
}

test "sqrt: max mantissa, e odd (largest q achievable for shift = N)" {
    // m = 2^53 - 1 = 0x1FFFFFFFFFFFFF, e = -1 (odd). shift = 53.
    // n = (2^53 - 1) * 2^53 = 2^106 - 2^53.
    // q = 2^53 - 1 (since (2^53)² = 2^106 > n).
    // rem = n - (2^53 - 1)² = 2^53 - 1 = q. Round DOWN.
    // m' = 2^53 - 1, e' = (-1 - 53) / 2 = -27.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = -1, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x1FFFFFFFFFFFFF), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -27), r.exponent);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
}

test "sqrt: max mantissa, e even (mid-range q after shift = N - 1)" {
    // m = 2^53 - 1, e = 0. shift = 52.
    // n = (2^53 - 1) << 52. q lands in the middle of the normalized range.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    // Expected q ≈ floor(sqrt(2) * 2^52) -> 0x16A09E667F3BCC after rounding.
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x16A09E667F3BCC), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -26), r.exponent);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
}

test "sqrt: min mantissa, e odd (shift = N, rounding required)" {
    // m = 2^52, e = -1 (odd). shift = 53.
    // n = 2^52 << 53 = 2^105. q = floor(sqrt(2) * 2^52) → 0x16A09E667F3BCD.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -1, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x16A09E667F3BCD), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -27), r.exponent);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
}

test "sqrt: arbitrary mid-range mantissa" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x12345678ABCDEF, .exponent = 3, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
}

test "sqrt: max representable exponent" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = max_e - 1, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
    // sqrt halves the exponent, so the result never overflows.
    try std.testing.expect(!r.isInf());
}

test "sqrt: min representable exponent" {
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
    try std.testing.expect(!r.isZero());
}

test "sqrt: Dyadic(6, 4) path, sqrt(2)" {
    // 2 = 32 * 2^-4. m=32, e=-4 (even). N=6 (even). shift = 6.
    // n = 32 << 6 = 2048. q = 45 (45² = 2025), rem = 23. rem ≤ q -> round down.
    // m' = 45, e' = (-4 - 6) / 2 = -5. value = 45/32 = 1.40625.
    const x = zsl.Dyadic(6, 4){ .mantissa = 32, .exponent = -4, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try std.testing.expectEqual(@as(zsl.Dyadic(6, 4).Mantissa, 45), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(6, 4).Exponent, -5), r.exponent);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(6, 4), x, r);
}

test "sqrt: Dyadic(6, 4) path, sqrt(5) (rounds up)" {
    // 5 = 40 * 2^-3. m=40, e=-3 (odd). N=6 (even). Parity differ -> shift = 5.
    // n = 40 << 5 = 1280. q = 35 (35² = 1225), rem = 55. rem > q -> round UP.
    // m' = 36, e' = (-3 - 5) / 2 = -4. value = 36/16 = 2.25.
    const x = zsl.Dyadic(6, 4){ .mantissa = 40, .exponent = -3, .positive = true };
    const r = zsl.dyadic.sqrt(x);
    try std.testing.expectEqual(@as(zsl.Dyadic(6, 4).Mantissa, 36), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(6, 4).Exponent, -4), r.exponent);

    try expectCorrectlyRoundedSqrt(zsl.Dyadic(6, 4), x, r);
}

test "sqrt: monotonicity" {
    inline for (.{
        .{ @as(f64, 0.001), @as(f64, 0.01) },
        .{ @as(f64, 0.5), @as(f64, 1.0) },
        .{ @as(f64, 1.0), @as(f64, 2.0) },
        .{ @as(f64, 2.0), @as(f64, 4.0) },
        .{ @as(f64, 100.0), @as(f64, 200.0) },
        .{ @as(f64, 1e10), @as(f64, 1e20) },
    }) |pair| {
        const x1 = zsl.Dyadic(53, 11).initValue(pair[0]);
        const x2 = zsl.Dyadic(53, 11).initValue(pair[1]);
        const y1 = zsl.dyadic.sqrt(x1);
        const y2 = zsl.dyadic.sqrt(x2);

        const order = y1.cmp(y2);
        try std.testing.expect(order == .lt or order == .eq);
    }
}

test "sqrt: positive output for positive input" {
    inline for (.{
        @as(f64, 0.001),
        @as(f64, 0.5),
        @as(f64, 1.0),
        @as(f64, 2.0),
        @as(f64, 100.0),
        @as(f64, 1e10),
        @as(f64, 1e100),
    }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        const y = zsl.dyadic.sqrt(x);
        try std.testing.expectEqual(true, y.positive);
        try std.testing.expect(!y.isZero());
        try std.testing.expect(!y.isNan());
        try std.testing.expect(!y.isInf());
    }
}

test "sqrt: 0ulp for many specific values" {
    inline for (.{
        @as(f64, 2.0),   @as(f64, 3.0),       @as(f64, 5.0),
        @as(f64, 7.0),   @as(f64, 10.0),      @as(f64, 0.1),
        @as(f64, 0.5),   @as(f64, 0.7),       @as(f64, 0.9),
        @as(f64, 1.5),   @as(f64, 2.5),       @as(f64, 1e-100),
        @as(f64, 1e100), @as(f64, 1234.5678), @as(f64, 0.000123),
    }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        const r = zsl.dyadic.sqrt(x);
        try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, r);
    }
}

test "sqrt: randomized testing" {
    const n = 1_000_000;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    for (0..n) |_| {
        const f = rand.floatNorm(f64) * if (rand.boolean())
            zsl.numeric.cast(f64, rand.int(u32))
        else
            @as(f64, 1.0) / zsl.numeric.max(1.0, zsl.numeric.cast(f64, rand.int(u32)));
        const x = zsl.Dyadic(53, 11).initValue(f);
        const y = zsl.dyadic.sqrt(x);

        try expectCorrectlyRoundedSqrt(zsl.Dyadic(53, 11), x, y);
    }
}

test "sqrt: exhaustive testing" {
    // const D = zsl.Dyadic(15, 6); // ~4M values
    const D = zsl.Dyadic(11, 5); // ~61K values
    // const D = zsl.Dyadic(6, 4); // ~900 values
    const allocator = std.testing.allocator;

    var values: std.ArrayList(D) = .empty;
    defer values.deinit(allocator);

    try values.appendSlice(allocator, &.{ D.nan, D.inf, D.inf.neg(), D.zero, D.zero.neg() });

    const m_min: D.Mantissa = @as(D.Mantissa, 1) << (@typeInfo(D.Mantissa).int.bits - 1);
    const m_max: D.Mantissa = std.math.maxInt(D.Mantissa);
    const e_min: D.Exponent = std.math.minInt(D.Exponent) + 1;
    const e_max: D.Exponent = std.math.maxInt(D.Exponent) - 1;

    var m: D.Mantissa = m_min;
    while (true) : (m +%= 1) {
        var e: D.Exponent = e_min;
        while (true) : (e +%= 1) {
            try values.append(allocator, .{ .mantissa = m, .exponent = e, .positive = true });
            try values.append(allocator, .{ .mantissa = m, .exponent = e, .positive = false });

            if (e == e_max)
                break;
        }
        if (m == m_max)
            break;
    }

    for (values.items) |x| {
        const y = zsl.dyadic.sqrt(x);

        if (!isCorrectlyRoundedSqrt(D, x, y)) {
            std.debug.print("x: {}\ny: {}\n", .{ x, y });
            return error.Fail;
        }
    }
}
