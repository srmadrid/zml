const std = @import("std");

const zsl = @import("zsl");
const mul = zsl.dyadic.mul;
const tzsl = @import("../zsl.zig");

test mul {
    const n = 10_000_000;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.Io.failing).toNanoseconds()));
    const rand = prng.random();

    for (0..n) |i| {
        const x_f = rand.floatNorm(f64) * rand.floatNorm(f64) * if (rand.boolean())
            zsl.numeric.cast(f64, rand.intRangeAtMost(u32, 1, 1_000_000_000))
        else
            1 / zsl.numeric.cast(f64, rand.intRangeAtMost(u32, 1, 1_000_000_000));
        const x_d = zsl.numeric.cast(zsl.Dyadic(53, 11), x_f);

        const y_f = rand.floatNorm(f64) * rand.floatNorm(f64) * if (rand.boolean())
            zsl.numeric.cast(f64, rand.intRangeAtMost(u32, 1, 1_000_000_000))
        else
            1 / zsl.numeric.cast(f64, rand.intRangeAtMost(u32, 1, 1_000_000_000));
        const y_d = zsl.numeric.cast(zsl.Dyadic(53, 11), y_f);

        const r_f = zsl.float.mul(x_f, y_f);
        const r_d = zsl.dyadic.mul(x_d, y_d);

        const f_comps = tzsl.dyadic.extractComps(r_f);

        const r_f_positive: bool = f_comps.positive;
        const r_d_positive: bool = r_d.positive;

        const r_f_exponent: i11 = f_comps.exponent - 52;
        const r_d_exponent: i11 = r_d.exponent;

        const r_f_mantissa: u53 = f_comps.mantissa;
        const r_d_mantissa: u53 = r_d.mantissa;

        if (r_f_positive != r_d_positive or
            r_f_exponent != r_d_exponent or
            r_f_mantissa != r_d_mantissa)
        {
            const x_comps = tzsl.dyadic.extractComps(x_f);
            const y_comps = tzsl.dyadic.extractComps(y_f);

            std.debug.print(
                "Test number: {}\nx: {d},  positive ({})  mantissa ({b:0>53})  exponent ({b:0>11})\ny: {d},  positive ({})  mantissa ({b:0>53})  exponent ({b:0>11})\n--mul--\nf: {d},  positive ({})  mantissa ({b:0>53})  exponent ({b:0>11})\nd: {d},  positive ({})  mantissa ({b:0>53})  exponent ({b:0>11})\n",
                .{
                    i,
                    x_f,
                    x_comps.positive,
                    x_comps.mantissa,
                    @as(u11, @bitCast(x_comps.exponent)),
                    y_f,
                    y_comps.positive,
                    y_comps.mantissa,
                    @as(u11, @bitCast(y_comps.exponent)),
                    r_f,
                    r_f_positive,
                    r_f_mantissa,
                    @as(u11, @bitCast(r_f_exponent)),
                    zsl.numeric.cast(f64, r_d),
                    r_d_positive,
                    r_d_mantissa,
                    @as(u11, @bitCast(r_d_exponent)),
                },
            );

            return error.MultiplicationMismatch;
        }
    }
}
