const std = @import("std");

const zsl = @import("zsl");

test "add: nan + finite -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) + @as(f64, 1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) + @as(f64, -1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one.neg()),
    );
}

test "add: finite + nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + std.math.nan(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) + std.math.nan(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "add: nan + nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) + std.math.nan(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).nan),
    );
}

test "add: nan + inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) + std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).inf),
    );
}

test "add: nan + -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) + -std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "add: inf + nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) + std.math.nan(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).nan),
    );
}

test "add: -inf + nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) + std.math.nan(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "add: nan + 0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) + @as(f64, 0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).zero),
    );
}

test "add: nan + -0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) + @as(f64, -0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "add: 0.0 + nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) + std.math.nan(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).nan),
    );
}

test "add: -0.0 + nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) + std.math.nan(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "add: inf + inf -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) + std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).inf),
    );
}

test "add: -inf + -inf -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) + -std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "add: inf + -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) + -std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "add: -inf + inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) + std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).inf),
    );
}

test "add: inf + finite -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) + @as(f64, 1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) + @as(f64, -1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one.neg()),
    );
}

test "add: -inf + finite -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) + @as(f64, 1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) + @as(f64, -1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one.neg()),
    );
}

test "add: finite + inf -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) + std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).inf),
    );
}

test "add: finite + -inf -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + -std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf.neg()),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) + -std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "add: inf + 0.0 -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) + @as(f64, 0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).zero),
    );
}

test "add: inf + -0.0 -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) + @as(f64, -0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "add: -inf + 0.0 -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) + @as(f64, 0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).zero),
    );
}

test "add: -inf + -0.0 -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) + @as(f64, -0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "add: 0.0 + inf -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) + std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).inf),
    );
}

test "add: -0.0 + inf -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) + std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).inf),
    );
}

test "add: 0.0 + -inf -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) + -std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "add: -0.0 + -inf -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) + -std.math.inf(f64)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "add: 0.0 + x -> x" {
    inline for (.{
        @as(f64, 1.0),
        @as(f64, -42.0),
        @as(f64, 51966.0),
    }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) + v),
            zsl.dyadic.add(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).initValue(v)),
        );
    }
}

test "add: -0.0 + x -> x" {
    inline for (.{
        @as(f64, 1.0),
        @as(f64, -42.0),
        @as(f64, 51966.0),
    }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) + v),
            zsl.dyadic.add(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).initValue(v)),
        );
    }
}

test "add: x + 0.0 -> x" {
    inline for (.{
        @as(f64, 1.0),
        @as(f64, -42.0),
        @as(f64, 51966.0),
    }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v + @as(f64, 0.0)),
            zsl.dyadic.add(zsl.Dyadic(53, 11).initValue(v), zsl.Dyadic(53, 11).zero),
        );
    }
}

test "add: x + -0.0 -> x" {
    inline for (.{
        @as(f64, 1.0),
        @as(f64, -42.0),
        @as(f64, 51966.0),
    }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v + @as(f64, -0.0)),
            zsl.dyadic.add(zsl.Dyadic(53, 11).initValue(v), zsl.Dyadic(53, 11).zero.neg()),
        );
    }
}

test "add: 0.0 + 0.0 -> 0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) + @as(f64, 0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero),
    );
}

test "add: 0.0 + -0.0 -> 0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) + @as(f64, -0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "add: -0.0 + 0.0 -> 0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) + @as(f64, 0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero),
    );
}

test "add: -0.0 + -0.0 -> -0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) + @as(f64, -0.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "add: 1.0 + 1.0 -> 2.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + @as(f64, 1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one),
    );
}

test "add: 1.0 + 2.0 -> 3.0 (|x| < |y|)" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + @as(f64, 2.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).initValue(@as(f64, 2.0))),
    );
}

test "add: 2.0 + 1.0 -> 3.0 (|x| > |y|)" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 2.0) + @as(f64, 1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)), zsl.Dyadic(53, 11).one),
    );
}

test "add: 1.0 + 0.5 -> 1.5" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + @as(f64, 0.5)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).initValue(@as(f64, 0.5))),
    );
}

test "add: -1.0 + -1.0 -> -2.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) + @as(f64, -1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).one.neg()),
    );
}

test "add: -3.0 + -5.0 -> -8.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -3.0) + @as(f64, -5.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).initValue(@as(f64, -3.0)), zsl.Dyadic(53, 11).initValue(@as(f64, -5.0))),
    );
}

test "add: same exponent forces carry-from-addition" {
    // x_m = 0x10000000000000, y_m = 0x10000000000000, exp_diff = 0.
    // sum = 0x20000000000000, carry. After carry-handling: m = 0x10000000000000, e += 1.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 1), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: max mantissa same exponent" {
    // x_m = 0x1FFFFFFFFFFFFF, y_m = 0x1FFFFFFFFFFFFF, exp_diff = 0.
    // sum = 0x3FFFFFFFFFFFFE, carry. After carry-handling: m = 0x1FFFFFFFFFFFFF, e += 1.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x1FFFFFFFFFFFFF), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 1), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: small exponent diff, no rounding" {
    // x_m = 0x10000000000000, y_m = 0x10000000000000, exp_diff = 5. y_shifted = 0x10000000000000 >> 5 = 0x00800000000000.
    // sum = 0x10800000000000. m = 0x10800000000000, e = 0. remainder = 0.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -5, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10800000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 0), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: round down" {
    // x_m = 0x10000000000000, y_m = 0x10000000000000, exp_diff = 54. y_shifted = 0x10000000000000 >> 54 = 0. sticky bits = 0.
    // sum = 0x10000000000000. m = 0x10000000000000. remainder < halfway. No round.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -54, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 0), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: round up" {
    // x_m = 0x10000000000000, y_m = 0x18000000000000, exp_diff = 53. y_shifted = 0x18000000000000 >> 53 = 0. sticky = 0.
    // sum = 0x10000000000000. m = 0x10000000000000. remainder > halfway. Round up to 0x10000000000001.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x18000000000000, .exponent = -53, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000001), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 0), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: tie to even, mantissa even -> no round" {
    // x_m = 0x10000000000000, y_m = 0x10000000000000, exp_diff = 53. y_shifted = 0x10000000000000 >> 53 = 0. sticky = 0.
    // sum = 0x10000000000000. m = 0x10000000000000 (even). remainder = 0x10000000000000 = halfway. Tie, no round.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -53, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 0), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: tie to even, mantissa odd -> round up" {
    // x_m = 0x10000000000001, y_m = 0x10000000000000, exp_diff = 53. y_shifted = 0x10000000000000 >> 53 = 0. sticky = 0.
    // sum = 0x10000000000001. m = 0x10000000000001 (odd). remainder = 0x10000000000000 = halfway. Tie + odd -> round up to 0x10000000000002.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -53, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000002), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 0), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: rounding overflows mantissa, carries exponent" {
    // x_m = 0x1FFFFFFFFFFFFF, y_m = 0x10000000000000, exp_diff = 53. y_shifted = 0x10000000000000 >> 53 = 0. sticky = 0.
    // sum = 0x1FFFFFFFFFFFFF. m = 0x1FFFFFFFFFFFFF (odd). remainder = halfway. Round up.
    // mantissa wraps to 0x20000000000000 -> reset to 0x10000000000000 and bump exponent. m = 0x10000000000000, e = 1.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -53, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 1), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: exp_diff >= 2 * mantissa_bits (y entirely sticky)" {
    // x_m = 0x10000000000001, y_m = 0x10000000000000, exp_diff = 106 >= 2 * 53. y_shifted = 0, sticky = 1.
    // sum = x_wide. remainder = 0 (no tie). Result == x.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -106, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 0), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: extreme exp_diff" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = -1000, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(x, r);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: exponent overflow returns inf-shape (positive)" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf, r);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: exponent overflow returns inf-shape (negative)" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = false };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = false };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf.neg(), r);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: exponent overflow via rounding carry" {
    // x.m = 0x1FFFFFFFFFFFFF, y.m = 0x10000000000000 with exp_diff = 53 at max_e - 1: rounds up, mantissa carries, exp bumps to max_e.
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = max_e - 1 - 53, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf, r);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: 1.0 + -1.0 -> 0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + @as(f64, -1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one.neg()),
    );
}

test "add: 2.0 + -1.0 -> 1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 2.0) + @as(f64, -1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)), zsl.Dyadic(53, 11).one.neg()),
    );
}

test "add: 1.0 + -2.0 -> -1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + @as(f64, -2.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).initValue(@as(f64, -2.0))),
    );
}

test "add: -2.0 + 1.0 -> -1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -2.0) + @as(f64, 1.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).initValue(@as(f64, -2.0)), zsl.Dyadic(53, 11).one),
    );
}

test "add: -1.0 + 2.0 -> 1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) + @as(f64, 2.0)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).initValue(@as(f64, 2.0))),
    );
}

test "add: 1.0 + -0.5 -> 0.5" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) + @as(f64, -0.5)),
        zsl.dyadic.add(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).initValue(@as(f64, -0.5))),
    );
}

test "add: equal magnitudes opposite signs -> 0.0" {
    inline for (.{
        @as(f64, 1.0),
        @as(f64, 42.0),
        @as(f64, 3735928559.0),
    }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v + -v),
            zsl.dyadic.add(zsl.Dyadic(53, 11).initValue(v), zsl.Dyadic(53, 11).initValue(v).neg()),
        );
    }
}

test "add: cancellation requires single-bit normalization" {
    // x = 1.0, y = -(2^53 - 2) * 2^-53 = -0x1FFFFFFFFFFFFE * 2^-53. |x| > |y|.
    // exp_diff = 1, y_shifted = 0x0FFFFFFFFFFFFF..., diff = 0x10000000000000..., lz = 52.
    // Result: 0x10000000000000 * 2^-67 = 2^-15.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -15, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFE, .exponent = -16, .positive = false };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -67), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: cancellation with sub-rounding (round up)" {
    // x.m = 0x10000000000000, y = -0x10000000000001 at e - 53. exp_diff = 53, y_shifted = 0x10000000000001, sticky = 0.
    // diff = ..., lz = 1. m = 0x1FFFFFFFFFFFFE, remainder = 0x1FFFFFFFFFFFFE > halfway -> round up to 0x1FFFFFFFFFFFFF.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = -53, .positive = false };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x1FFFFFFFFFFFFF), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -1), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: cancellation with sticky from y (sticky decrements diff)" {
    // exp_diff = 54 > 53: y.mantissa lsb falls below cutoff -> sticky = 1.
    // y.m = 0x10000000000001 -> y_shifted = 0x08000000000000, sticky = 1.
    // diff -= 1 due to sticky. lz = 1. m = 0x1FFFFFFFFFFFFF, remainder = 0x0FFFFFFFFFFFFE < halfway. No round.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = -54, .positive = false };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x1FFFFFFFFFFFFF), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -1), r.exponent);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: result sign tracks the larger magnitude" {
    // |y| > |x|, signs differ -> result sign = y.positive.
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 1.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, -10.0));
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(false, r.positive);

    const x2 = zsl.Dyadic(53, 11).initValue(@as(f64, -1.0));
    const y2 = zsl.Dyadic(53, 11).initValue(@as(f64, 10.0));
    const r2 = zsl.dyadic.add(x2, y2);
    try std.testing.expectEqual(true, r2.positive);
}

test "add: cancellation underflow returns 0.0" {
    // Near min_e, |x| - |y| = 1 ULP -> after normalization exponent <= min_e -> .zero.
    // Caller sets result.positive = x.positive = true.
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = min_e + 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = false };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero, r);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: cancellation underflow with negative dominant -> -0.0" {
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = min_e + 1, .positive = false };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero.neg(), r);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: commutativity (x + y == y + x)" {
    inline for (.{
        .{ @as(f64, 1.0), @as(f64, 2.0) },      .{ @as(f64, -5.0), @as(f64, 7.0) },
        .{ @as(f64, 100.0), @as(f64, -100.0) }, .{ @as(f64, 12345.0), @as(f64, -6789.0) },
        .{ @as(f64, 0.0), @as(f64, 42.0) },     .{ @as(f64, -1000.0), @as(f64, -2000.0) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const y = zsl.Dyadic(53, 11).initValue(tc[1]);
        try std.testing.expectEqual(zsl.dyadic.add(x, y), zsl.dyadic.add(y, x));

        try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(tc[0] + tc[1]), zsl.dyadic.add(x, y));
        try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(tc[1] + tc[0]), zsl.dyadic.add(y, x));
    }
}

test "add: identity x + 0.0 == 0.0 + x == x" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0), @as(f64, -48879.0) }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        try std.testing.expectEqual(x, zsl.dyadic.add(x, zsl.Dyadic(53, 11).zero));
        try std.testing.expectEqual(x, zsl.dyadic.add(zsl.Dyadic(53, 11).zero, x));
    }
}

test "add: deep cancellation requires massive left shift" {
    // x_m = 0x10000000000000, e = 0. Value = 2^52
    // y_m = 0x1FFFFFFFFFFFFE, e = -1. Value = -(2^53 - 2) * 2^-1 = -(2^52 - 1)
    // exp_diff = 1. y_shifted = 0x0FFFFFFFFFFFFF.
    // diff = 1. Normalizing requires shifting left by 52 bits to restore the hidden/highest bit.
    // Result: m = 0x10000000000000, e = x.e - 52.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFE, .exponent = -1, .positive = false };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000000), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, -52), r.exponent);
    try std.testing.expectEqual(true, r.positive);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: exact cancellation always yields 0.0" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = false };

    const r1 = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero, r1);

    const r2 = zsl.dyadic.add(y, x);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero, r2);
}

test "add: absolute boundary of underflow" {
    // If x and y cancel out to exactly 1 ULP at min_e, normalizing requires
    // left shifting which immediately drops the exponent below min_e, triggering 0.0.
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = min_e, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e, .positive = false };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero, r);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: double rounding illusion (halfway + distant sticky bit)" {
    // exp_diff = 53. y_m = 0x10000000000001. y_shifted = 0.
    // remainder is 0x10000000000000 (exactly halfway) BUT sticky = 1 because of the LSB of y_m.
    // Because sticky = 1, discarded > halfway. It MUST round UP.
    // x_m is even, so a naïve tie-to-even check that forgets the sticky bit would fail and not round.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = -53, .positive = true };
    const r = zsl.dyadic.add(x, y);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Mantissa, 0x10000000000001), r.mantissa);
    try std.testing.expectEqual(@as(zsl.Dyadic(53, 11).Exponent, 0), r.exponent);

    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) + y.toFloat(f64)), r);
}

test "add: randomized testing" {
    const n = 1_000_000;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();

    for (0..n) |_| {
        const x = rand.floatNorm(f64) * if (rand.boolean())
            zsl.numeric.cast(f64, rand.int(u32))
        else
            @as(f64, 1.0) / zsl.numeric.max(1.0, zsl.numeric.cast(f64, rand.int(u32)));
        const y = rand.floatNorm(f64) * if (rand.boolean())
            zsl.numeric.cast(f64, rand.int(u32))
        else
            @as(f64, 1.0) / zsl.numeric.max(1.0, zsl.numeric.cast(f64, rand.int(u32)));

        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(x + y),
            zsl.dyadic.add(zsl.Dyadic(53, 11).initValue(x), zsl.Dyadic(53, 11).initValue(y)),
        );
    }
}

test "add: exhaustive testing" {
    // const D = zsl.Dyadic(11, 5); // ~3.77B pairs
    const D = zsl.Dyadic(6, 4); // ~812K pairs
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
        for (values.items) |y| {
            const expected = D.initValue(x.toFloat(f64) + y.toFloat(f64));
            const actual = zsl.dyadic.add(x, y);

            if (expected.isNan())
                try std.testing.expect(actual.isNan())
            else
                std.testing.expectEqual(expected, actual) catch {
                    std.debug.print("x:        {}\ny:        {}\nexpected: {}\nactual:   {}\n", .{ x, y, expected, actual });

                    return error.Fail;
                };
        }
    }
}
