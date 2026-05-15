const std = @import("std");

const zsl = @import("zsl");

test "mul: nan * finite -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) * @as(f64, 1.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) * @as(f64, -1.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one.neg()),
    );
}

test "mul: finite * nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) * std.math.nan(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) * std.math.nan(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "mul: nan * nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) * std.math.nan(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).nan),
    );
}

test "mul: nan * inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) * std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).inf),
    );
}

test "mul: nan * -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) * -std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "mul: inf * nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) * std.math.nan(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).nan),
    );
}

test "mul: -inf * nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) * std.math.nan(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "mul: nan * 0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) * @as(f64, 0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).zero),
    );
}

test "mul: nan * -0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.nan(f64) * @as(f64, -0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "mul: 0.0 * nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) * std.math.nan(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).nan),
    );
}

test "mul: -0.0 * nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) * std.math.nan(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "mul: inf * inf -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) * std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).inf),
    );
}

test "mul: -inf * -inf -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) * -std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "mul: inf * -inf -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) * -std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "mul: -inf * inf -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) * std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).inf),
    );
}

test "mul: inf * finite -> ±inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) * @as(f64, 1.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) * @as(f64, -1.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one.neg()),
    );
}

test "mul: -inf * finite -> ±inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) * @as(f64, 1.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) * @as(f64, -1.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one.neg()),
    );
}

test "mul: finite * inf -> ±inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) * std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) * std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).inf),
    );
}

test "mul: finite * -inf -> ±inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) * -std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf.neg()),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -1.0) * -std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "mul: inf * 0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) * @as(f64, 0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).zero),
    );
}

test "mul: inf * -0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(std.math.inf(f64) * @as(f64, -0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "mul: -inf * 0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) * @as(f64, 0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).zero),
    );
}

test "mul: -inf * -0.0 -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(-std.math.inf(f64) * @as(f64, -0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "mul: 0.0 * inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) * std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).inf),
    );
}

test "mul: -0.0 * inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) * std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).inf),
    );
}

test "mul: 0.0 * -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) * -std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "mul: -0.0 * -inf -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) * -std.math.inf(f64)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "mul: 0.0 * x -> ±0.0" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) * v),
            zsl.dyadic.mul(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).initValue(v)),
        );
    }
}

test "mul: -0.0 * x -> ±0.0" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) * v),
            zsl.dyadic.mul(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).initValue(v)),
        );
    }
}

test "mul: x * 0.0 -> ±0.0" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v * @as(f64, 0.0)),
            zsl.dyadic.mul(zsl.Dyadic(53, 11).initValue(v), zsl.Dyadic(53, 11).zero),
        );
    }
}

test "mul: x * -0.0 -> ±0.0" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |v| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v * @as(f64, -0.0)),
            zsl.dyadic.mul(zsl.Dyadic(53, 11).initValue(v), zsl.Dyadic(53, 11).zero.neg()),
        );
    }
}

test "mul: 0.0 * 0.0 -> 0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) * @as(f64, 0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero),
    );
}

test "mul: 0.0 * -0.0 -> -0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) * @as(f64, -0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "mul: -0.0 * 0.0 -> -0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) * @as(f64, 0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero),
    );
}

test "mul: -0.0 * -0.0 -> 0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, -0.0) * @as(f64, -0.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "mul: 1.0 * 1.0 -> 1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 1.0) * @as(f64, 1.0)),
        zsl.dyadic.mul(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one),
    );
}

test "mul: 2.0 * 3.0 -> 6.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 2.0) * @as(f64, 3.0)),
        zsl.dyadic.mul(
            zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 3.0)),
        ),
    );
}

test "mul: 3.0 * 2.0 -> 6.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 3.0) * @as(f64, 2.0)),
        zsl.dyadic.mul(
            zsl.Dyadic(53, 11).initValue(@as(f64, 3.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
        ),
    );
}

test "mul: 0.5 * 0.5 -> 0.25" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 0.5) * @as(f64, 0.5)),
        zsl.dyadic.mul(
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.5)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.5)),
        ),
    );
}

test "mul: 2.0 * 0.5 -> 1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 2.0) * @as(f64, 0.5)),
        zsl.dyadic.mul(
            zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.5)),
        ),
    );
}

test "mul: 7.0 * 7.0 -> 49.0 (odd squaring)" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(@as(f64, 7.0) * @as(f64, 7.0)),
        zsl.dyadic.mul(
            zsl.Dyadic(53, 11).initValue(@as(f64, 7.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 7.0)),
        ),
    );
}

test "mul: (+)(+) -> (+)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 3.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, 5.0));
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(true, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(@as(f64, 3.0) * @as(f64, 5.0)), r);
}

test "mul: (-)(-) -> (+)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, -3.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, -5.0));
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(true, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(@as(f64, -3.0) * @as(f64, -5.0)), r);
}

test "mul: (+)(-) -> (-)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 3.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, -5.0));
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(false, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(@as(f64, 3.0) * @as(f64, -5.0)), r);
}

test "mul: (-)(+) -> (-)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, -3.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, 5.0));
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(false, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(@as(f64, -3.0) * @as(f64, 5.0)), r);
}

test "mul: smallest * smallest mantissa (product msb at 2N-2)" {
    // m_x = m_y = 2^(N-1). Product = 2^(2N-2), msb at bit 2N-2 -> shift by N-1.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)),
        zsl.dyadic.mul(x, y),
    );
}

test "mul: largest * largest mantissa (product msb at 2N-1)" {
    // m_x = m_y = 2^N - 1. Product ≈ 2^(2N), msb at bit 2N-1 -> shift by N.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)),
        zsl.dyadic.mul(x, y),
    );
}

test "mul: mid-range mantissas (rounding likely)" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1234567890ABCD, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 0, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)),
        zsl.dyadic.mul(x, y),
    );
}

test "mul: mantissas where rounding carries into exponent" {
    // Product just below the next power of 2, rounding flips it over.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 5, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 5, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)),
        zsl.dyadic.mul(x, y),
    );
}

test "mul: exponent addition only (powers of 2)" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 100, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -200, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)),
        zsl.dyadic.mul(x, y),
    );
}

test "mul: exponent overflow -> +inf" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 100, .positive = true };
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf, r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)), r);
}

test "mul: exponent overflow with mixed sign -> -inf" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 100, .positive = false };
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf.neg(), r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)), r);
}

test "mul: exponent overflow with both negative -> +inf" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = false };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 100, .positive = false };
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf, r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)), r);
}

test "mul: exponent underflow -> +0.0" {
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -100, .positive = true };
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero, r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)), r);
}

test "mul: exponent underflow with mixed sign -> -0.0" {
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -100, .positive = false };
    const r = zsl.dyadic.mul(x, y);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).zero.neg(), r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(x.toFloat(f64) * y.toFloat(f64)), r);
}

test "mul: commutativity (x * y == y * x)" {
    @setEvalBranchQuota(10_000);

    inline for (.{
        .{ @as(f64, 1.0), @as(f64, 2.0) },    .{ @as(f64, -5.0), @as(f64, 7.0) },
        .{ @as(f64, 100.0), @as(f64, -0.5) }, .{ @as(f64, 12345.0), @as(f64, -6789.0) },
        .{ @as(f64, 0.0), @as(f64, 42.0) },   .{ @as(f64, -1000.0), @as(f64, -2000.0) },
        .{ @as(f64, 3.14), @as(f64, 2.71) },  .{ @as(f64, 1e10), @as(f64, 1e-10) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const y = zsl.Dyadic(53, 11).initValue(tc[1]);
        try std.testing.expectEqual(zsl.dyadic.mul(x, y), zsl.dyadic.mul(y, x));

        try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(tc[0] * tc[1]), zsl.dyadic.mul(x, y));
        try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(tc[1] * tc[0]), zsl.dyadic.mul(y, x));
    }
}

test "mul: identity x * 1.0 == 1.0 * x == x" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0), @as(f64, -48879.0), @as(f64, 0.001) }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        try std.testing.expectEqual(x, zsl.dyadic.mul(x, zsl.Dyadic(53, 11).one));
        try std.testing.expectEqual(x, zsl.dyadic.mul(zsl.Dyadic(53, 11).one, x));
    }
}

test "mul: annihilator x * 0.0 -> ±0.0" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0), @as(f64, -48879.0) }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(v * @as(f64, 0.0)),
            zsl.dyadic.mul(x, zsl.Dyadic(53, 11).zero),
        );
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.0) * v),
            zsl.dyadic.mul(zsl.Dyadic(53, 11).zero, x),
        );
    }
}

test "mul: x * -1.0 == -x" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0), @as(f64, 3.14) }) |v| {
        const x = zsl.Dyadic(53, 11).initValue(v);
        try std.testing.expectEqual(x.neg(), zsl.dyadic.mul(x, zsl.Dyadic(53, 11).one.neg()));
        try std.testing.expectEqual(x.neg(), zsl.dyadic.mul(zsl.Dyadic(53, 11).one.neg(), x));
    }
}

test "mul: randomized testing" {
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
            zsl.Dyadic(53, 11).initValue(x * y),
            zsl.dyadic.mul(zsl.Dyadic(53, 11).initValue(x), zsl.Dyadic(53, 11).initValue(y)),
        );
    }
}

test "mul: exhaustive testing" {
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
            const expected = D.initValue(x.toFloat(f64) * y.toFloat(f64));
            const actual = zsl.dyadic.mul(x, y);

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
