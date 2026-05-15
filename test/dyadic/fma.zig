const std = @import("std");

const zsl = @import("zsl");

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
            zsl.dyadic.fma(zsl.Dyadic(53, 11).initValue(x), zsl.Dyadic(53, 11).initValue(y), zsl.Dyadic(53, 11).initValue(z)),
        );
    }
}
