const std = @import("std");

const zsl = @import("zsl");

test "sub: exhaustive testing" {
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
            const expected = D.initValue(x.toFloat(f64) - y.toFloat(f64));
            const actual = zsl.dyadic.sub(x, y);

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
