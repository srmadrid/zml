const std = @import("std");

const zsl = @import("zsl");

test "fma: nan * y + z -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.nan(f64), @as(f64, 1.0), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.nan(f64), @as(f64, -1.0), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).one),
    );
}

test "fma: x * nan + z -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), std.math.nan(f64), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -1.0), std.math.nan(f64), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one),
    );
}

test "fma: x * y + nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), @as(f64, 1.0), std.math.nan(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), @as(f64, -1.0), std.math.nan(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).nan),
    );
}

test "fma: nan * nan + nan -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.nan(f64), std.math.nan(f64), std.math.nan(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).nan),
    );
}

test "fma: nan with inf in any position -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.nan(f64), std.math.inf(f64), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.nan(f64), @as(f64, 1.0), std.math.inf(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.inf(f64), std.math.nan(f64), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), std.math.nan(f64), std.math.inf(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).inf),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.inf(f64), @as(f64, 1.0), std.math.nan(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), std.math.inf(f64), std.math.nan(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).nan),
    );
}

test "fma: nan with zero in any position -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.nan(f64), @as(f64, 0.0), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.nan(f64), @as(f64, 1.0), @as(f64, 0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).zero),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.0), std.math.nan(f64), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).one),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), std.math.nan(f64), @as(f64, 0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan, zsl.Dyadic(53, 11).zero),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.0), @as(f64, 1.0), std.math.nan(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).nan),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), @as(f64, 0.0), std.math.nan(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).nan),
    );
}

test "fma: inf * 0.0 + z -> nan (regardless of z)" {
    inline for (.{
        @as(f64, 1.0),     @as(f64, -1.0),
        @as(f64, 0.0),     @as(f64, -0.0),
        std.math.inf(f64), -std.math.inf(f64),
        @as(f64, 42.0),    @as(f64, -3735928559.0),
    }) |z_f| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.inf(f64), @as(f64, 0.0), z_f)),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).initValue(z_f)),
        );
    }
}

test "fma: -inf * 0.0 + z -> nan (regardless of z)" {
    inline for (.{
        @as(f64, 1.0),     @as(f64, -1.0),
        @as(f64, 0.0),     @as(f64, -0.0),
        std.math.inf(f64), -std.math.inf(f64),
    }) |z_f| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(-std.math.inf(f64), @as(f64, 0.0), z_f)),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).initValue(z_f)),
        );
    }
}

test "fma: 0.0 * inf + z -> nan (regardless of z)" {
    inline for (.{
        @as(f64, 1.0),     @as(f64, -1.0),
        std.math.inf(f64), -std.math.inf(f64),
        @as(f64, 42.0),
    }) |z_f| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.0), std.math.inf(f64), z_f)),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).initValue(z_f)),
        );
    }
}

test "fma: 0.0 * -inf + z -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.0), -std.math.inf(f64), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one),
    );
}

test "fma: -0.0 * inf + z -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -0.0), std.math.inf(f64), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one),
    );
}

test "fma: inf * -0.0 + z -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.inf(f64), @as(f64, -0.0), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).one),
    );
}

test "fma: (+inf) + (-inf) -> nan (product inf, z opposite inf)" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.inf(f64), @as(f64, 1.0), -std.math.inf(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf.neg()),
    );

    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -1.0), std.math.inf(f64), std.math.inf(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one.neg(), zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).inf),
    );
}

test "fma: (-inf) + (+inf) -> nan" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(-std.math.inf(f64), @as(f64, 1.0), std.math.inf(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf),
    );
}

test "fma: inf + inf (same sign) -> inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.inf(f64), @as(f64, 1.0), std.math.inf(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf),
    );
}

test "fma: -inf + -inf -> -inf" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(-std.math.inf(f64), @as(f64, 1.0), -std.math.inf(f64))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).inf.neg()),
    );
}

test "fma: inf * finite + finite -> ±inf" {
    inline for (.{
        .{ @as(f64, 2.0), @as(f64, 3.0) },
        .{ @as(f64, -2.0), @as(f64, 3.0) },
        .{ @as(f64, 2.0), @as(f64, -3.0) },
        .{ @as(f64, -2.0), @as(f64, -3.0) },
    }) |tc| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(std.math.inf(f64), tc[0], tc[1])),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).initValue(tc[0]), zsl.Dyadic(53, 11).initValue(tc[1])),
        );

        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(-std.math.inf(f64), tc[0], tc[1])),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).inf.neg(), zsl.Dyadic(53, 11).initValue(tc[0]), zsl.Dyadic(53, 11).initValue(tc[1])),
        );
    }
}

test "fma: finite * inf + finite -> ±inf" {
    inline for (.{
        .{ @as(f64, 2.0), @as(f64, 3.0) },
        .{ @as(f64, -2.0), @as(f64, -3.0) },
    }) |tc| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(tc[0], std.math.inf(f64), tc[1])),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).initValue(tc[0]), zsl.Dyadic(53, 11).inf, zsl.Dyadic(53, 11).initValue(tc[1])),
        );
    }
}

test "fma: finite * finite + ±inf -> ±inf" {
    inline for (.{
        .{ @as(f64, 2.0), @as(f64, 3.0) },
        .{ @as(f64, -2.0), @as(f64, -3.0) },
    }) |tc| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(tc[0], tc[1], std.math.inf(f64))),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).initValue(tc[0]), zsl.Dyadic(53, 11).initValue(tc[1]), zsl.Dyadic(53, 11).inf),
        );

        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(tc[0], tc[1], -std.math.inf(f64))),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).initValue(tc[0]), zsl.Dyadic(53, 11).initValue(tc[1]), zsl.Dyadic(53, 11).inf.neg()),
        );
    }
}

test "fma: 0.0 * y + z -> z" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 0.5), @as(f64, 51966.0) }) |z_f| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.0), @as(f64, 5.0), z_f)),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).initValue(@as(f64, 5.0)), zsl.Dyadic(53, 11).initValue(z_f)),
        );
    }
}

test "fma: x * 0.0 + z -> z" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 0.5), @as(f64, 51966.0) }) |z_f| {
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 5.0), @as(f64, 0.0), z_f)),
            zsl.dyadic.fma(zsl.Dyadic(53, 11).initValue(@as(f64, 5.0)), zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).initValue(z_f)),
        );
    }
}

test "fma: x * y + 0.0 -> mul(x, y)" {
    inline for (.{
        .{ @as(f64, 2.0), @as(f64, 3.0) },
        .{ @as(f64, -5.0), @as(f64, 7.0) },
        .{ @as(f64, 0.5), @as(f64, 0.5) },
        .{ @as(f64, 3.14), @as(f64, 2.71) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const y = zsl.Dyadic(53, 11).initValue(tc[1]);
        try std.testing.expectEqual(zsl.dyadic.mul(x, y), zsl.dyadic.fma(x, y, zsl.Dyadic(53, 11).zero));
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(tc[0], tc[1], @as(f64, 0.0))),
            zsl.dyadic.fma(x, y, zsl.Dyadic(53, 11).zero),
        );
    }
}

test "fma: 0.0 * 0.0 + 0.0 -> 0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.0), @as(f64, 0.0), @as(f64, 0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero),
    );
}

test "fma: 0.0 * 0.0 + -0.0 -> 0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.0), @as(f64, 0.0), @as(f64, -0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "fma: -0.0 * -0.0 + 0.0 -> 0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -0.0), @as(f64, -0.0), @as(f64, 0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero),
    );
}

test "fma: -0.0 * -0.0 + -0.0 -> -0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -0.0), @as(f64, -0.0), @as(f64, -0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "fma: 0.0 * -0.0 + -0.0 -> -0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.0), @as(f64, -0.0), @as(f64, -0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero.neg()),
    );
}

test "fma: -0.0 * 0.0 + 0.0 -> 0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -0.0), @as(f64, 0.0), @as(f64, 0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).zero.neg(), zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).zero),
    );
}

test "fma: 1.0 * 1.0 + 0.0 -> 1.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), @as(f64, 1.0), @as(f64, 0.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).zero),
    );
}

test "fma: 1.0 * 1.0 + 1.0 -> 2.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 1.0), @as(f64, 1.0), @as(f64, 1.0))),
        zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).one),
    );
}

test "fma: 2.0 * 3.0 + 1.0 -> 7.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 2.0), @as(f64, 3.0), @as(f64, 1.0))),
        zsl.dyadic.fma(
            zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 3.0)),
            zsl.Dyadic(53, 11).one,
        ),
    );
}

test "fma: 0.5 * 0.5 + 1.0 -> 1.25" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 0.5), @as(f64, 0.5), @as(f64, 1.0))),
        zsl.dyadic.fma(
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.5)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 0.5)),
            zsl.Dyadic(53, 11).one,
        ),
    );
}

test "fma: 2.0 * 3.0 + -6.0 -> 0.0 (exact cancellation)" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 2.0), @as(f64, 3.0), @as(f64, -6.0))),
        zsl.dyadic.fma(
            zsl.Dyadic(53, 11).initValue(@as(f64, 2.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 3.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, -6.0)),
        ),
    );
}

test "fma: -2.0 * 3.0 + 6.0 -> 0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -2.0), @as(f64, 3.0), @as(f64, 6.0))),
        zsl.dyadic.fma(
            zsl.Dyadic(53, 11).initValue(@as(f64, -2.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 3.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 6.0)),
        ),
    );
}

test "fma: 7.0 * 7.0 + -49.0 -> 0.0" {
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 7.0), @as(f64, 7.0), @as(f64, -49.0))),
        zsl.dyadic.fma(
            zsl.Dyadic(53, 11).initValue(@as(f64, 7.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, 7.0)),
            zsl.Dyadic(53, 11).initValue(@as(f64, -49.0)),
        ),
    );
}

test "fma: (+)(+) + (+) -> (+)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 2.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, 3.0));
    const z = zsl.Dyadic(53, 11).initValue(@as(f64, 1.0));
    const r = zsl.dyadic.fma(x, y, z);
    try std.testing.expectEqual(true, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 2.0), @as(f64, 3.0), @as(f64, 1.0))), r);
}

test "fma: (-)(-) + (+) -> (+)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, -2.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, -3.0));
    const z = zsl.Dyadic(53, 11).initValue(@as(f64, 1.0));
    const r = zsl.dyadic.fma(x, y, z);
    try std.testing.expectEqual(true, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -2.0), @as(f64, -3.0), @as(f64, 1.0))), r);
}

test "fma: (+)(-) + (+) -> sign of larger magnitude" {
    // x * y = -6.0, z = 1.0. -6.0 + 1.0 = -5.0 (negative).
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 2.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, -3.0));
    const z = zsl.Dyadic(53, 11).initValue(@as(f64, 1.0));
    const r = zsl.dyadic.fma(x, y, z);
    try std.testing.expectEqual(false, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 2.0), @as(f64, -3.0), @as(f64, 1.0))), r);
}

test "fma: (+)(+) + (-) -> sign of larger magnitude" {
    // x * y = 6.0, z = -1.0. 6.0 + (-1.0) = 5.0 (positive).
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 2.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, 3.0));
    const z = zsl.Dyadic(53, 11).initValue(@as(f64, -1.0));
    const r = zsl.dyadic.fma(x, y, z);
    try std.testing.expectEqual(true, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, 2.0), @as(f64, 3.0), @as(f64, -1.0))), r);
}

test "fma: (-)(+) + (-) -> (-)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, -2.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, 3.0));
    const z = zsl.Dyadic(53, 11).initValue(@as(f64, -1.0));
    const r = zsl.dyadic.fma(x, y, z);
    try std.testing.expectEqual(false, r.positive);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(zsl.float.fma(@as(f64, -2.0), @as(f64, 3.0), @as(f64, -1.0))), r);
}

test "fma: product msb at 2 * mantissa_bits - 1 (large mantissas)" {
    // m_x * m_y near 2^(2 * mantissa_bits) -> product msb lands at 2 * mantissa_bits - 1, no left shift needed.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const z = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 50, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: product msb at 2 * mantissa_bits - 2 (small mantissas, requires left shift)" {
    // m_x = m_y = 2^(mantissa_bits - 1) -> product = 2^(2 * mantissa_bits - 2), msb at 2 * mantissa_bits - 2 -> shift left by 1.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const z = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 50, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: |z| >> |x * y| (z dominates, x * y contributes sticky)" {
    const x = zsl.Dyadic(53, 11).initValue(@as(f64, 1.0));
    const y = zsl.Dyadic(53, 11).initValue(@as(f64, 1.0));
    const z = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 50, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: |x * y| >> |z| (product dominates, z contributes sticky)" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 50, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 50, .positive = true };
    const z = zsl.Dyadic(53, 11).one;
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: z fully shifted out (exp_diff >= 2 * mantissa_bits, z is sticky only)" {
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 100, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const z = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -300, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: addition with carry (sum crosses 2 * mantissa_bits bit boundary)" {
    // Choose x, y, z such that |x * y| ≈ |z| and same sign -> sum carries into bit 2N.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const z = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 53, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: cancellation z ≈ -x * y (small residual)" {
    // fma's characteristic case: x * y + (-x * y_rounded) recovers the rounding error.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1234567890ABCD, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1ABCDEF0123456, .exponent = 0, .positive = true };
    const z = zsl.dyadic.mul(x, y).neg();
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: cancellation requires deep renormalization (large left shift)" {
    // Subtraction near-equal magnitudes -> diff has many leading zeros, exp decreases sharply.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = 0, .positive = true };
    // x * y = 1.0. z = -1.0 + tiny -> almost cancels.
    const z = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFE, .exponent = -53, .positive = false };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: exponent overflow -> +inf" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 100, .positive = true };
    const z = zsl.Dyadic(53, 11).one;
    const r = zsl.dyadic.fma(x, y, z);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf, r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))), r);
}

test "fma: exponent overflow with mixed sign -> -inf" {
    const max_e = zsl.int.maxVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = max_e - 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 100, .positive = false };
    const z = zsl.Dyadic(53, 11).one;
    const r = zsl.dyadic.fma(x, y, z);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).inf.neg(), r);
    try std.testing.expectEqual(zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))), r);
}

test "fma: exponent underflow -> +0.0" {
    const min_e = zsl.int.minVal(zsl.Dyadic(53, 11).Exponent);
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = -100, .positive = true };
    const z = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000000, .exponent = min_e + 1, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: rounding tie - sticky from product affects tie-break" {
    // Construct x, y where the product has a halfway-rounding case, with z contributing a sticky bit.
    const x = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const y = zsl.Dyadic(53, 11){ .mantissa = 0x1FFFFFFFFFFFFF, .exponent = 0, .positive = true };
    const z = zsl.Dyadic(53, 11){ .mantissa = 0x10000000000001, .exponent = -53, .positive = true };
    try std.testing.expectEqual(
        zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64))),
        zsl.dyadic.fma(x, y, z),
    );
}

test "fma: fma(x, 1, z) == add(x, z)" {
    inline for (.{
        .{ @as(f64, 1.0), @as(f64, 2.0) },
        .{ @as(f64, -5.0), @as(f64, 7.0) },
        .{ @as(f64, 100.0), @as(f64, -0.5) },
        .{ @as(f64, 12345.0), @as(f64, -6789.0) },
        .{ @as(f64, 0.0), @as(f64, 42.0) },
        .{ @as(f64, 3.14), @as(f64, 2.71) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const z = zsl.Dyadic(53, 11).initValue(tc[1]);
        try std.testing.expectEqual(zsl.dyadic.add(x, z), zsl.dyadic.fma(x, zsl.Dyadic(53, 11).one, z));
    }
}

test "fma: fma(1, y, z) == add(y, z)" {
    inline for (.{
        .{ @as(f64, 1.0), @as(f64, 2.0) },
        .{ @as(f64, -5.0), @as(f64, 7.0) },
        .{ @as(f64, 0.5), @as(f64, -0.25) },
    }) |tc| {
        const y = zsl.Dyadic(53, 11).initValue(tc[0]);
        const z = zsl.Dyadic(53, 11).initValue(tc[1]);
        try std.testing.expectEqual(zsl.dyadic.add(y, z), zsl.dyadic.fma(zsl.Dyadic(53, 11).one, y, z));
    }
}

test "fma: fma(x, y, 0) == mul(x, y)" {
    inline for (.{
        .{ @as(f64, 2.0), @as(f64, 3.0) },
        .{ @as(f64, -5.0), @as(f64, 7.0) },
        .{ @as(f64, 0.5), @as(f64, 0.5) },
        .{ @as(f64, 3.14), @as(f64, 2.71) },
        .{ @as(f64, 1e10), @as(f64, 1e-10) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const y = zsl.Dyadic(53, 11).initValue(tc[1]);
        try std.testing.expectEqual(zsl.dyadic.mul(x, y), zsl.dyadic.fma(x, y, zsl.Dyadic(53, 11).zero));
    }
}

test "fma: fma(0, y, z) == z (for finite y, z)" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |z_f| {
        const z = zsl.Dyadic(53, 11).initValue(z_f);
        try std.testing.expectEqual(z, zsl.dyadic.fma(zsl.Dyadic(53, 11).zero, zsl.Dyadic(53, 11).one, z));
    }
}

test "fma: fma(x, 0, z) == z (for finite x, z)" {
    inline for (.{ @as(f64, 1.0), @as(f64, -42.0), @as(f64, 51966.0) }) |z_f| {
        const z = zsl.Dyadic(53, 11).initValue(z_f);
        try std.testing.expectEqual(z, zsl.dyadic.fma(zsl.Dyadic(53, 11).one, zsl.Dyadic(53, 11).zero, z));
    }
}

test "fma: commutativity in x, y (fma(x, y, z) == fma(y, x, z))" {
    @setEvalBranchQuota(10_000);

    inline for (.{
        .{ @as(f64, 2.0), @as(f64, 3.0), @as(f64, 1.0) },
        .{ @as(f64, -5.0), @as(f64, 7.0), @as(f64, -2.0) },
        .{ @as(f64, 100.0), @as(f64, -0.5), @as(f64, 1.5) },
        .{ @as(f64, 3.14), @as(f64, 2.71), @as(f64, 0.0) },
        .{ @as(f64, 1e10), @as(f64, 1e-10), @as(f64, 42.0) },
        .{ @as(f64, 12345.0), @as(f64, -6789.0), @as(f64, 0.001) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const y = zsl.Dyadic(53, 11).initValue(tc[1]);
        const z = zsl.Dyadic(53, 11).initValue(tc[2]);
        try std.testing.expectEqual(zsl.dyadic.fma(x, y, z), zsl.dyadic.fma(y, x, z));
    }
}

test "fma: sign symmetry fma(-x, y, -z) == fma(x, -y, -z) == -fma(x, y, z)" {
    @setEvalBranchQuota(10_000);

    inline for (.{
        .{ @as(f64, 2.0), @as(f64, 3.0), @as(f64, 1.0) },
        .{ @as(f64, -5.0), @as(f64, 7.0), @as(f64, -2.0) },
        .{ @as(f64, 0.5), @as(f64, 0.5), @as(f64, 1.0) },
        .{ @as(f64, 3.14), @as(f64, 2.71), @as(f64, -0.5) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const y = zsl.Dyadic(53, 11).initValue(tc[1]);
        const z = zsl.Dyadic(53, 11).initValue(tc[2]);
        const r = zsl.dyadic.fma(x, y, z);
        try std.testing.expectEqual(r.neg(), zsl.dyadic.fma(x.neg(), y, z.neg()));
        try std.testing.expectEqual(r.neg(), zsl.dyadic.fma(x, y.neg(), z.neg()));
        try std.testing.expectEqual(r, zsl.dyadic.fma(x.neg(), y.neg(), z));
    }
}

test "fma: error term fma(x, y, -mul(x, y)) recovers rounding residual exactly" {
    // The classic application: with true fma, this expression gives the
    // exact rounding error of mul(x, y). With sequential mul + add, it always returns 0.
    inline for (.{
        .{ @as(f64, 1.0) + @as(f64, 1.0) / @as(f64, 7.0), @as(f64, 1.0) + @as(f64, 1.0) / @as(f64, 13.0) },
        .{ @as(f64, 1.1), @as(f64, 3.7) },
        .{ @as(f64, 1.234567890), @as(f64, 9.876543210) },
    }) |tc| {
        const x = zsl.Dyadic(53, 11).initValue(tc[0]);
        const y = zsl.Dyadic(53, 11).initValue(tc[1]);
        const xy = zsl.dyadic.mul(x, y);
        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), -xy.toFloat(f64))),
            zsl.dyadic.fma(x, y, xy.neg()),
        );
    }
}

test "fma: randomized testing" {
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
        const z = rand.floatNorm(f64) * if (rand.boolean())
            zsl.numeric.cast(f64, rand.int(u32))
        else
            @as(f64, 1.0) / zsl.numeric.max(1.0, zsl.numeric.cast(f64, rand.int(u32)));

        try std.testing.expectEqual(
            zsl.Dyadic(53, 11).initValue(zsl.float.fma(x, y, z)),
            zsl.dyadic.fma(
                zsl.Dyadic(53, 11).initValue(x),
                zsl.Dyadic(53, 11).initValue(y),
                zsl.Dyadic(53, 11).initValue(z),
            ),
        );
    }
}

test "fma: exhaustive testing" {
    // const D = zsl.Dyadic(11, 5); // ~232T triplets
    // const D = zsl.Dyadic(7, 4); // ~5.8B triples
    // const D = zsl.Dyadic(6, 4); // ~730M triples
    // const D = zsl.Dyadic(5, 3); // ~7.6M triples
    const D = zsl.Dyadic(4, 3); // ~1M triples
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
            for (values.items) |z| {
                const expected = D.initValue(zsl.float.fma(x.toFloat(f64), y.toFloat(f64), z.toFloat(f64)));
                const actual = zsl.dyadic.fma(x, y, z);

                if (expected.isNan())
                    try std.testing.expect(actual.isNan())
                else
                    std.testing.expectEqual(expected, actual) catch {
                        std.debug.print("x:        {}\ny:        {}\nz:        {}\nexpected: {}\nactual:   {}\n", .{ x, y, z, expected, actual });

                        return error.Fail;
                    };
            }
        }
    }
}
