const std = @import("std");

const zsl = @import("zsl");

test "div: nan / finite -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) / @as(f64, 1.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) / @as(f64, -1.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one.neg()),
    );
}

test "div: finite / nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) / std.math.nan(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) / std.math.nan(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "div: nan / nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) / std.math.nan(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).nan),
    );
}

test "div: nan / inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) / std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).inf),
    );
}

test "div: nan / -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) / -std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "div: inf / nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) / std.math.nan(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).nan),
    );
}

test "div: -inf / nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) / std.math.nan(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "div: nan / 0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) / @as(f64, 0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).zero),
    );
}

test "div: nan / -0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) / @as(f64, -0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "div: 0.0 / nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) / std.math.nan(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).nan),
    );
}

test "div: -0.0 / nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) / std.math.nan(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "div: inf / inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) / std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).inf),
    );
}

test "div: -inf / -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) / -std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "div: inf / -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) / -std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "div: -inf / inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) / std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).inf),
    );
}

test "div: inf / finite -> ±inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) / @as(f64, 1.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) / @as(f64, -1.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one.neg()),
    );
}

test "div: -inf / finite -> ±inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) / @as(f64, 1.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) / @as(f64, -1.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one.neg()),
    );
}

test "div: finite / inf -> ±0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) / std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) / std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).inf),
    );
}

test "div: finite / -inf -> ±0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) / -std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf.neg()),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) / -std.math.inf(f64)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "div: inf / 0.0 -> +inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) / @as(f64, 0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).zero),
    );
}

test "div: inf / -0.0 -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) / @as(f64, -0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "div: -inf / 0.0 -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) / @as(f64, 0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).zero),
    );
}

test "div: -inf / -0.0 -> +inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) / @as(f64, -0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "div: 0.0 / 0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) / @as(f64, 0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero),
    );
}

test "div: 0.0 / -0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) / @as(f64, -0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "div: -0.0 / 0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) / @as(f64, 0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero),
    );
}

test "div: -0.0 / -0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) / @as(f64, -0.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "div: finite / 0.0 -> ±inf" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v / @as(f64, 0.0)),
            zsl.dyadic.div(zsl.Dyadic(53, 11).initValue(v), zsl.Dyadic(53, 11).zero),
        );
    }
}

test "div: finite / -0.0 -> ∓inf" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v / @as(f64, -0.0)),
            zsl.dyadic.div(zsl.Dyadic(53, 11).initValue(v), zsl.Dyadic(53, 11).zero.neg()),
        );
    }
}

test "div: 0.0 / x -> ±0.0" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) / v),
            zsl.dyadic.div(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).initValue(v)),
        );
    }
}

test "div: -0.0 / x -> ±0.0" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) / v),
            zsl.dyadic.div(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).initValue(v)),
        );
    }
}

test "div: 1.0 / 1.0 -> 1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) / @as(f64, 1.0)),
        zsl.dyadic.div(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one),
    );
}

test "div: 6.0 / 2.0 -> 3.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 6.0) / @as(f64, 2.0)),
        zsl.dyadic.div(
            zsl.Dyadic(53, 11).initValue(@as(f64, 6.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
        ),
    );
}

test "div: 6.0 / 3.0 -> 2.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 6.0) / @as(f64, 3.0)),
        zsl.dyadic.div(
            zsl.Dyadic(53, 11).initValue(@as(f64, 6.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 3.0)),
        ),
    );
}

test "div: 1.0 / 2.0 -> 0.5" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) / @as(f64, 2.0)),
        zsl.dyadic.div(
            zsl.Dyadic(53, 11).one,
            zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
        ),
    );
}

test "div: 1.0 / 3.0 (inexact, rounding required)" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) / @as(f64, 3.0)),
        zsl.dyadic.div(
            zsl.Dyadic(53, 11).one,
            zsl.Dyadic(53, 11).initValue(@as(f64, 3.0)),
        ),
    );
}

test "div: 2.0 / 4.0 -> 0.5" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 2.0) / @as(f64, 4.0)),
        zsl.dyadic.div(
            zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 4.0)),
        ),
    );
}

test "div: 7.0 / 7.0 -> 1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 7.0) / @as(f64, 7.0)),
        zsl.dyadic.div(
            zsl.Dyadic(53, 11).initValue(@as(f64, 7.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 7.0)),
        ),
    );
}

test "div: 1.0 / x (reciprocals)" {
    inline for (.{ @as(f64, 2.0), @as(f64, 3.0), @as(f64, 7.0), @as(f64, -5.0), @as(f64, 0.5) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) / v),
            zsl.dyadic.div(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).initValue(v)),
        );
    }
}

test "div: (+)/(+) -> (+)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 15.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, 3.0));
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(true, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(@as(f64, 15.0) / @as(f64, 3.0)), r);
}

test "div: (-)/(-) -> (+)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, -15.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, -3.0));
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(true, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(@as(f64, -15.0) / @as(f64, -3.0)), r);
}

test "div: (+)/(-) -> (-)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 15.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, -3.0));
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(false, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(@as(f64, 15.0) / @as(f64, -3.0)), r);
}

test "div: (-)/(+) -> (-)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, -15.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, 3.0));
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(false, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(@as(f64, -15.0) / @as(f64, 3.0)), r);
}

test "div: x.mantissa == y.mantissa, same exponent -> 1.0" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 0, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)),
        zsl.dyadic.div(x, y),
    );
}

test "div: division by power of two (exact)" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 3, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)),
        zsl.dyadic.div(x, y),
    );
}

test "div: small mantissa / large mantissa (quotient renormalization)" {
    // m_x < m_y -> quotient mantissa < 2^(N-1), needs left shift and exponent decrement.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)),
        zsl.dyadic.div(x, y),
    );
}

test "div: large mantissa / small mantissa (quotient already normalized)" {
    // m_x > m_y -> quotient mantissa MSB lands at the expected position.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)),
        zsl.dyadic.div(x, y),
    );
}

test "div: mid-range mantissas (rounding likely)" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1234567890ABCD, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 0, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)),
        zsl.dyadic.div(x, y),
    );
}

test "div: mantissas where rounding carries into exponent" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 5, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = 0, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)),
        zsl.dyadic.div(x, y),
    );
}

test "div: exponent subtraction only (powers of two)" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 100, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -200, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)),
        zsl.dyadic.div(x, y),
    );
}

test "div: exponent overflow -> +inf" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf, r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)), r);
}

test "div: exponent overflow with mixed sign -> -inf" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = false };
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf.neg(), r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)), r);
}

test "div: exponent overflow with both negative -> +inf" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = false };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = false };
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf, r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)), r);
}

test "div: exponent underflow -> +0.0" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = max_e - 1, .positive = true };
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero, r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)), r);
}

test "div: exponent underflow with mixed sign -> -0.0" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = max_e - 1, .positive = false };
    const r = zsl.dyadic.div(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero.neg(), r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) / y.toFloat(f64)), r);
}

test "div: identity x / 1.0 == x" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0), @as(f64, -48879.0), @as(f64, 0.001) }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        try std.testing.expectEqual(x, zsl.dyadic.div(x, zsl.Dyadic(53, 11).one));
    }
}

test "div: x / x -> 1.0" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0), @as(f64, 3.14) }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v / v),
            zsl.dyadic.div(x, x),
        );
    }
}

test "div: x / -1.0 == -x" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0), @as(f64, 3.14) }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        try std.testing.expectEqual(x.neg(), zsl.dyadic.div(x, zsl.Dyadic(53, 11).one.neg()));
    }
}

test "div: -x / y == -(x / y) and x / -y == -(x / y)" {
    @setEvalBranchQuota(10_000);

    inline for (.{
        .{ @as(f64, 6.0), @as(f64, 2.0) },   .{ @as(f64, 1.0), @as(f64, 3.0) },
        .{ @as(f64, 100.0), @as(f64, 7.0) }, .{ @as(f64, 12345.0), @as(f64, 89.0) },
        .{ @as(f64, 3.14), @as(f64, 2.71) }, .{ @as(f64, 1e10), @as(f64, 1e-10) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const y = zsl.Dyadic(53, 11).initValue(tc[1]);
        const xy = zsl.dyadic.div(x, y);
        try std.testing.expectEqual(xy.neg(), zsl.dyadic.div(x.neg(), y));
        try std.testing.expectEqual(xy.neg(), zsl.dyadic.div(x, y.neg()));
        try std.testing.expectEqual(xy, zsl.dyadic.div(x.neg(), y.neg()));
    }
}

test "div: randomized testing" {
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
            zsl.Dyadic(53, 11).initValue(x / y),
            zsl.dyadic.div(zsl.Dyadic(53, 11).initValue(x), zsl.Dyadic(53, 11).initValue(y)),
        );
    }
}

test "div: exhaustive testing" {
    // const D = zsl.Dyadic(11, 5);
    const D = zsl.Dyadic(6, 4);
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
            const expected = D.initValue(x.toFloat(f64) / y.toFloat(f64));
            const actual = zsl.dyadic.div(x, y);

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
