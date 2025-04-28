const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const exp2_data = @import("exp2_data.zig");
const exp_data = @import("exp_data.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub fn exp(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return exp(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, exp32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_expf.c
                    return exp32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_exp.c
                    return exp64(x);
                },
                f80 => return cast(f80, exp128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_expl.c
                    return exp128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

inline fn top12(x: f32) u32 {
    return @as(u32, @bitCast(x)) >> 20;
}

fn exp32(x: f32) f32 {
    const xd: f64 = cast(f64, x, .{});
    const abstop: u32 = (@as(u32, @bitCast(x)) >> 20) & 0x7ff;
    if (abstop >= @as(u32, @bitCast(@as(f32, (88)))) >> 20) {
        @branchHint(.unlikely);
        // |x| >= 88 or x is nan.
        if (@as(u32, @bitCast(x)) == @as(u32, @bitCast(-std.math.inf(f32)))) {
            return 0;
        }
        if (abstop >= @as(u32, @bitCast(std.math.inf(f32))) >> 20) {
            return x + x;
        }
        if (x > 0x1.62e42ep6) { // x > log(0x1p128) ~= 88.72
            return std.math.inf(f32);
        }
        if (x < -0x1.9fe368p6) { // x < log(0x1p-150) ~= -103.97
            return 0;
        }
        if (x < -0x1.9d1d9ep6) { // x < log(0x1p-149) ~= -103.28
            return 0;
        }
    }

    // x*N/Ln2 = k + r with r in [-1/2, 1/2] and int k.
    var z: f64 = exp2_data.InvLn2N_32 * xd;

    // Round and convert z to int, the result is in [-150*32, 128*32] and
    // ideally ties-to-even rule is used, otherwise the magnitude of r
    // can be bigger which gives larger approximation error.
    var kd: f64 = z + exp2_data.SHIFT_32;
    const ki: u64 = @bitCast(kd);
    kd -= exp2_data.SHIFT_32;
    const r: f64 = z - kd;

    // exp(x) = 2^(k/N) * 2^(r/N) ~= s * (C0*r^3 + C1*r^2 + C2*r + 1)
    var t: u64 = exp2_data.T_32[ki % 32];
    t +%= ki << 47;
    const s: f64 = @bitCast(t);
    z = exp2_data.poly_scaled_32[0] * r + exp2_data.poly_scaled_32[1];
    const r2: f64 = r * r;
    var y: f64 = exp2_data.poly_scaled_32[2] * r + 1;
    y = z * r2 + y;
    y = y * s;
    return cast(f32, y, .{});
}

// Handle cases that may overflow or underflow when computing the result that
// is scale*(1+TMP) without intermediate rounding.  The bit representation of
// scale is in SBITS, however it has a computed exponent that may have
// overflown into the sign bit so that needs to be adjusted before using it as
// a double.  (int32_t)KI is the k used in the argument reduction and exponent
// adjustment of scale, positive k here means the result may overflow and
// negative k means the result may underflow.  */
inline fn specialcase(tmp: f64, sbits: u64, ki: u64) f64 {
    var sb: u64 = sbits;
    if ((ki & 0x80000000) == 0) {
        // k > 0, the exponent of scale might have overflowed by <= 460.
        sb -= 1009 << 52;
        const scale: f64 = @bitCast(sb);
        const y: f64 = 0x1p1009 * (scale + scale * tmp);
        return y;
    }
    // k < 0, need special care in the subnormal range.
    sb +%= 1022 << 52;
    const scale: f64 = @bitCast(sb);
    var y: f64 = scale + scale * tmp;
    if (y < 1) {
        // Round y to the right precision before scaling it into the subnormal
        // range to avoid double rounding that can cause 0.5+E/2 ulp error where
        // E is the worst-case ulp error outside the subnormal range.  So this
        // is only useful if the goal is better than 1 ulp worst-case error.
        var lo: f64 = scale - y + scale * tmp;
        const hi: f64 = 1 + y;
        lo = 1 - hi + y + lo;
        y = hi + lo - 1;
        // Avoid -0.0 with downward rounding.
        if (y == 0)
            y = 0;
        // The underflow exception needs to be signaled explicitly.
        std.mem.doNotOptimizeAway(0x1p-1022 * 0x1p-1022);
    }
    y = 0x1p-1022 * y;
    return y;
}

fn exp64(x: f64) f64 {
    var abstop: u32 = @intCast((@as(u64, @bitCast(x)) >> 52) & 0x7ff);
    if (abstop -% (@as(u64, @bitCast(@as(f64, 0x1p-54))) >> 52) >= ((@as(u64, @bitCast(@as(f64, 512))) >> 52) - (@as(u64, @bitCast(@as(f64, 0x1p-54))) >> 52))) {
        @branchHint(.unlikely);
        if (abstop -% (@as(u64, @bitCast(@as(f64, 0x1p-54))) >> 52) >= 0x80000000) {
            // Avoid spurious underflow for tiny x.
            // Note: 0 is common input.
            return 1 + x;
        }

        if (abstop >= ((@as(u64, @bitCast(@as(f64, 1024)))) >> 52)) {
            if (@as(u64, @bitCast(x)) == @as(u64, @bitCast(-std.math.inf(f64))))
                return 0;

            if (abstop >= (@as(u64, @bitCast(@as(f64, std.math.inf(f64)))) >> 52))
                return 1 + x;

            if ((@as(u64, @bitCast(x)) >> 63) != 0) {
                return 0;
            } else {
                return std.math.inf(f64);
            }
        }
        // Large x is special cased below.
        abstop = 0;
    }

    // exp(x) = 2^(k/128) * exp(r), with exp(r) in [2^(-1/256),2^(1/256)].
    // x = ln2/128*k + r, with int k and r in [-ln2/256, ln2/256].  */
    const z: f64 = exp_data.InvLn2N_64 * x;
    // z - kd is in [-1, 1] in non-nearest rounding modes.
    var kd: f64 = z + exp_data.Shift_64;
    const ki: u64 = @bitCast(kd);
    kd -= exp_data.Shift_64;
    const r: f64 = x + kd * exp_data.NegLn2hiN_64 + kd * exp_data.NegLn2loN_64;
    // 2^(k/128) ~= scale * (1 + tail).
    const idx: u64 = 2 * (ki % 128);
    const top: u64 = ki << 45;
    const tail: f64 = @bitCast(exp_data.T_64[idx]);
    // This is only a valid scale when -1023*128 < k < 1024*128.
    const sbits: u64 = exp_data.T_64[idx + 1] +% top;
    // exp(x) = 2^(k/128) * exp(r) ~= scale + scale * (tail + exp(r) - 1).
    // Evaluation is optimized assuming superscalar pipelined execution.
    const r2: f64 = r * r;
    // Without fma the worst case error is 0.25/N ulp larger.
    // Worst case error is less than 0.5+1.11/N+(abs poly error * 2^53) ulp.
    const tmp: f64 = tail + r + r2 * (exp_data.poly[0] + r * exp_data.poly[1]) + r2 * r2 * (exp_data.poly[2] + r * exp_data.poly[3]);
    if (abstop == 0) {
        @branchHint(.unlikely);
        return specialcase(tmp, sbits, ki);
    }
    const scale: f64 = @bitCast(sbits);
    // Note: tmp == 0 or |tmp| > 2^-65 and scale > 2^-739, so there
    // is no spurious underflow here even without fma.
    return scale + scale * tmp;
}

fn exp128(x: f128) f128 {
    // Smallest integer x for which e^x overflows.
    const himark: f128 = 11356.523406294143949491931077970765;
    // Largest integer x for which e^x underflows.
    const lomark: f128 = -11433.4627433362978788372438434526231;
    // 3x2^96
    const THREEp96: f128 = 59421121885698253195157962752;
    // 3x2^103
    const THREEp103: f128 = 30423614405477505635920876929024;
    // 3x2^111
    const THREEp111: f128 = 7788445287802241442795744493830144;
    // 1/ln(2)
    const M_1_LN2: f128 = 1.44269504088896340735992468100189204;
    // first 93 bits of ln(2)
    const M_LN2_0: f128 = 0.693147180559945309417232121457981864;
    // ln2_0 - ln(2)
    const M_LN2_1: f128 = -1.94704509238074995158795957333327386e-31;
    // very small number
    const TINY: f128 = 1.0e-4900;
    // 2^16383
    const TWO16383: f128 = 5.94865747678615882542879663314003565e+4931;
    // 256
    const TWO8: f128 = 256;
    // 32768
    const TWO15: f128 = 32768;
    // Chebyshev polynom coefficients for (exp(x)-1)/x
    const P1: f128 = 0.5;
    const P2: f128 = 1.66666666666666666666666666666666683e-01;
    const P3: f128 = 4.16666666666666666666654902320001674e-02;
    const P4: f128 = 8.33333333333333333333314659767198461e-03;
    const P5: f128 = 1.38888888889899438565058018857254025e-03;
    const P6: f128 = 1.98412698413981650382436541785404286e-04;

    // Check for usual case.
    if (x < himark and x > lomark) {

        // Calculate n.
        var n: f128 = x * M_1_LN2 + THREEp111;
        n -= THREEp111;
        var xx: f128 = x - n * M_LN2_0;
        var xl: f128 = n * M_LN2_1;

        // Calculate t/256.
        var t: f128 = xx + THREEp103;
        t -= THREEp103;

        // Compute tval1 = t.
        const tval1: i32 = cast(i32, t * TWO8, .{});

        xx -= exp_data.__expl_table[@intCast(2 * 89 + 2 * tval1)];
        xl -= exp_data.__expl_table[@intCast(2 * 89 + 2 * tval1 + 1)];

        // Calculate t/32768.
        t = xx + THREEp96;
        t -= THREEp96;

        // Compute tval2 = t.
        const tval2: i32 = cast(i32, t * TWO15, .{});

        xx -= exp_data.__expl_table[@intCast((2 * (2 * 89) + 2 + 2 * 65) + 2 * tval2)];
        xl -= exp_data.__expl_table[@intCast((2 * (2 * 89) + 2 + 2 * 65) + 2 * tval2 + 1)];

        xx = xx + xl;

        // Compute ex2 = 2^n_0 e^(argtable[tval1]) e^(argtable[tval2]).
        var ex2_u: ldbl128.ieee_f128_shape = @bitCast(exp_data.__expl_table[@intCast(((2 * (2 * 89) + 2 + 2 * 65) + 2 + 2 * 65 + 89) + tval1)] * exp_data.__expl_table[@intCast((((2 * (2 * 89) + 2 + 2 * 65) + 2 + 2 * 65 + 89) + 1 + 89 + 65) + tval2)]);
        const n_i: i32 = cast(i32, n, .{});
        // 'unsafe' is 1 iff n_1 != 0.
        const unsafe: i32 = cast(i32, math.abs(n_i) >= 15000, .{});
        {
            @setRuntimeSafety(false);
            ex2_u.exponent += @intCast(n_i >> @intCast(unsafe));
        }

        // Compute scale = 2^n_1.
        var scale_u: ldbl128.ieee_f128_shape = @bitCast(@as(f128, 1));
        {
            @setRuntimeSafety(false);
            scale_u.exponent += @intCast(n_i - (n_i >> @intCast(unsafe)));
        }

        // Approximate e^x2 - 1, using a seventh-degree polynomial,
        // with maximum error in [-2^-16-2^-53,2^-16+2^-53]
        // less than 4.8e-39.
        const x22: f128 = xx + xx * xx * (P1 + xx * (P2 + xx * (P3 + xx * (P4 + xx * (P5 + xx * P6)))));
        std.mem.doNotOptimizeAway(x22);

        // Return result.
        var result: f128 = x22 * @as(f128, @bitCast(ex2_u)) + @as(f128, @bitCast(ex2_u));

        // Now we can test whether the result is ultimate or if we are unsure.
        // In the later case we should probably call a mpn based routine to give
        // the ultimate result.
        // Empirically, this routine is already ultimate in about 99.9986% of
        // cases, the test below for the round to nearest case will be false
        // in ~ 99.9963% of cases.
        // Without proc2 routine maximum error which has been seen is
        // 0.5000262 ulp.

        // union ieee854_long_double ex3_u;
        // ex3_u.d = (result - ex2_u.d) - x22 * ex2_u.d;
        // ex2_u.d = result;
        // ex3_u.exponent += LDBL_MANT_DIG + 15 + IEEE854_LONG_DOUBLE_BIAS
        //          - ex2_u.exponent;
        // n_i = abs (ex3_u.d);
        // n_i = (n_i + 1) / 2;
        // fesetenv (&oldenv);
        // #ifdef FE_TONEAREST
        // if (fegetround () == FE_TONEAREST)
        //   n_i -= 0x4000;
        // #endif
        // if (!n_i) {
        //   return __ieee754_expl_proc2 (origx);
        // }
        if (unsafe == 0) {
            return result;
        } else {
            result *= @bitCast(scale_u);

            if (result < std.math.floatMin(f128)) {
                const vresult: f128 = result * result;
                std.mem.doNotOptimizeAway(vresult);
            }

            return result;
        }
    }
    // Exceptional cases:
    else if (x < himark) {
        if (std.math.isInf(x)) {
            // e^-inf == 0, with no error.
            return 0;
        } else {
            // Underflow
            return TINY * TINY;
        }
    } else {
        // Return x, if x is a NaN or Inf; or overflow, otherwise.
        return TWO16383 * x;
    }
}

test exp {
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x2.b7e15p+0, exp(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x7.63993p+0, exp(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(0x1.415e5cp+4, exp(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(0x2.1df3b8p+0, exp(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x1.19103ep+72, exp(@as(f32, 0x3.2p+4)));
    try std.testing.expectEqual(0xf.ff684p+124, exp(@as(f32, 0x5.8b9028p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5cp+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x3.e8p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c6p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x4.d2p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.e870a4p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.e870a8p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.ebe224p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.ebe228p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c4edp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c469d8p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c469dcp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46d94p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46d98p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46724p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46728p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c469dcp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c469ep+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46c04p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46adcp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46aep+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c471bp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c471b4p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c4699p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46994p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c49fap+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c4ac1p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c4d89p+8)));
    try std.testing.expectEqual(0x1.004008p+0, exp(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0xf.fc008p-4, exp(@as(f32, -0x4p-12)));
    try std.testing.expectEqual(0x1.00001p+0, exp(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0xf.ffffp-4, exp(@as(f32, -0x1p-20)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x4p-32)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x1p-40)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0xf.fff84p+124, exp(@as(f32, 0x5.8b90b8p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x5.8b90cp+4)));
    try std.testing.expectEqual(0x3.ffff3p-128, exp(@as(f32, -0x5.75628p+4)));
    try std.testing.expectEqual(0x4.00013p-128, exp(@as(f32, -0x5.756278p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c86p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c86p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46578p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c4657cp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c46578p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c4657cp+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c86p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c86p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.9fa8dcp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.9fa8ep+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.9fa8dcp+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.9fa8ep+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c86p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c86p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5bd48p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5bd4cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5bd48p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5bd4cp+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c86p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c86p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp(@as(f32, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x3.b629a4p-4, exp(@as(f32, -0x1.760cdp+0)));
    try std.testing.expectEqual(0x3.b6299cp-4, exp(@as(f32, -0x1.760cd2p+0)));
    try std.testing.expectEqual(0x3.a823dp+0, exp(@as(f32, 0x1.4bed28p+0)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x1.f1cf36p+8)));
    try std.testing.expectEqual(0x1.72a52cp+4, exp(@as(f32, 0x3.248524p+0)));
    try std.testing.expectEqual(0x6.f5dcep+0, exp(@as(f32, 0x1.f0b362p+0)));
    try std.testing.expectEqual(0xb.8c7b8p+16, exp(@as(f32, 0xd.89747p+0)));
    try std.testing.expectEqual(0xb.8c7adp+16, exp(@as(f32, 0xd.89746p+0)));
    try std.testing.expectEqual(0xa.c2d26p-4, exp(@as(f32, -0x6.58b64p-4)));
    // try std.testing.expectEqual(0x1.0001fep+0, exp(@as(f32, 0x1.fefe02p-16)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.1895ep+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f32, -0x2.1895e4p+8)));

    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x2.b7e151628aed2p+0, exp(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x7.63992e35376b8p+0, exp(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(0x1.415e5bf6fb106p+4, exp(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(0x2.1df3b68cfb9fp+0, exp(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x1.19103e4080b45p+72, exp(@as(f64, 0x3.2p+4)));
    try std.testing.expectEqual(0xf.ff6844410e1f8p+124, exp(@as(f64, 0x5.8b9028p+4)));
    try std.testing.expectEqual(0xf.7c2d08f39f968p+1020, exp(@as(f64, 0x2.c5cp+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x3.e8p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c6p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x4.d2p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4p-1076, exp(@as(f64, -0x2.e870a4p+8)));
    try std.testing.expectEqual(0x4p-1076, exp(@as(f64, -0x2.e870a8p+8)));
    try std.testing.expectEqual(0x4p-1076, exp(@as(f64, -0x2.e870a7e5e88cp+8)));
    try std.testing.expectEqual(0x4p-1076, exp(@as(f64, -0x2.e870a7e5e88c2p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.ebe224p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.ebe228p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.ebe227861639p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c4edp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x3.eeb48df211be4p-1024, exp(@as(f64, -0x2.c469d8p+8)));
    try std.testing.expectEqual(0x3.eea4d33f4f708p-1024, exp(@as(f64, -0x2.c469dcp+8)));
    try std.testing.expectEqual(0x3.eeb09f3f7b25cp-1024, exp(@as(f64, -0x2.c469d9p+8)));
    try std.testing.expectEqual(0x3.e0206d364fe24p-1024, exp(@as(f64, -0x2.c46d94p+8)));
    try std.testing.expectEqual(0x3.e010ecd39be3p-1024, exp(@as(f64, -0x2.c46d98p+8)));
    try std.testing.expectEqual(0x3.e018acfd35b14p-1024, exp(@as(f64, -0x2.c46d96p+8)));
    try std.testing.expectEqual(0x3.f96438ed17acp-1024, exp(@as(f64, -0x2.c46724p+8)));
    try std.testing.expectEqual(0x3.f954537bfeefp-1024, exp(@as(f64, -0x2.c46728p+8)));
    try std.testing.expectEqual(0x3.f9584cd24f15cp-1024, exp(@as(f64, -0x2.c46727p+8)));
    try std.testing.expectEqual(0x3.eea4d33f4f708p-1024, exp(@as(f64, -0x2.c469dcp+8)));
    try std.testing.expectEqual(0x3.ee9518cb776fcp-1024, exp(@as(f64, -0x2.c469ep+8)));
    try std.testing.expectEqual(0x3.ee9cf5fd86364p-1024, exp(@as(f64, -0x2.c469dep+8)));
    try std.testing.expectEqual(0x3.e6335d7047b5cp-1024, exp(@as(f64, -0x2.c46c04p+8)));
    try std.testing.expectEqual(0x3.eab82516dd868p-1024, exp(@as(f64, -0x2.c46adcp+8)));
    try std.testing.expectEqual(0x3.eaa87a559ec28p-1024, exp(@as(f64, -0x2.c46aep+8)));
    try std.testing.expectEqual(0x3.eab04fae68c4p-1024, exp(@as(f64, -0x2.c46adep+8)));
    try std.testing.expectEqual(0x3.d053f45176d6cp-1024, exp(@as(f64, -0x2.c471bp+8)));
    try std.testing.expectEqual(0x3.d044b3202807cp-1024, exp(@as(f64, -0x2.c471b4p+8)));
    try std.testing.expectEqual(0x3.d0488366c34bp-1024, exp(@as(f64, -0x2.c471b3p+8)));
    try std.testing.expectEqual(0x3.efcfd88e9dc9p-1024, exp(@as(f64, -0x2.c4699p+8)));
    try std.testing.expectEqual(0x3.efc0196eb9e34p-1024, exp(@as(f64, -0x2.c46994p+8)));
    try std.testing.expectEqual(0x3.efc40930cb32cp-1024, exp(@as(f64, -0x2.c46993p+8)));
    try std.testing.expectEqual(0x3.2ff3a3a879b1p-1024, exp(@as(f64, -0x2.c49fap+8)));
    try std.testing.expectEqual(0x3.0941d1ff351b8p-1024, exp(@as(f64, -0x2.c4ac1p+8)));
    try std.testing.expectEqual(0x2.8d3d2f65d464cp-1024, exp(@as(f64, -0x2.c4d89p+8)));
    try std.testing.expectEqual(0x1.00400800aab55p+0, exp(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0xf.fc007ff556p-4, exp(@as(f64, -0x4p-12)));
    try std.testing.expectEqual(0x1.00001000008p+0, exp(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0xf.ffff000008p-4, exp(@as(f64, -0x1p-20)));
    try std.testing.expectEqual(0x1.00000004p+0, exp(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffffcp-4, exp(@as(f64, -0x4p-32)));
    try std.testing.expectEqual(0x1.0000000001p+0, exp(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0xf.fffffffffp-4, exp(@as(f64, -0x1p-40)));
    try std.testing.expectEqual(0x1.0000000000004p+0, exp(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffffcp-4, exp(@as(f64, -0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0xf.fff8417363ff8p+124, exp(@as(f64, 0x5.8b90b8p+4)));
    try std.testing.expectEqual(0x1.00000417184b8p+128, exp(@as(f64, 0x5.8b90cp+4)));
    try std.testing.expectEqual(0x3.ffff2fe5259dp-128, exp(@as(f64, -0x5.75628p+4)));
    try std.testing.expectEqual(0x4.00012fe53d8f8p-128, exp(@as(f64, -0x5.756278p+4)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f8p+1020, exp(@as(f64, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0xf.ffffffffff95p+1020, exp(@as(f64, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f8p+1020, exp(@as(f64, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0x4.000ebd79918d4p-1024, exp(@as(f64, -0x2.c46578p+8)));
    try std.testing.expectEqual(0x3.fffebd5e9bf24p-1024, exp(@as(f64, -0x2.c4657cp+8)));
    try std.testing.expectEqual(0x3.ffffffffff9fp-1024, exp(@as(f64, -0x2.c4657baf579a6p+8)));
    try std.testing.expectEqual(0x4.000ebd79918d4p-1024, exp(@as(f64, -0x2.c46578p+8)));
    try std.testing.expectEqual(0x3.fffebd5e9bf24p-1024, exp(@as(f64, -0x2.c4657cp+8)));
    try std.testing.expectEqual(0x4.00000000001fp-1024, exp(@as(f64, -0x2.c4657baf579a4p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f8p+1020, exp(@as(f64, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0xf.ffffffffff95p+1020, exp(@as(f64, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f8p+1020, exp(@as(f64, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0xf.ffffffffff95p+1020, exp(@as(f64, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(0x8.0005c84b6996p-972, exp(@as(f64, -0x2.9fa8dcp+8)));
    try std.testing.expectEqual(0x7.ffe5c8744841p-972, exp(@as(f64, -0x2.9fa8ep+8)));
    try std.testing.expectEqual(0x8.00000000009c8p-972, exp(@as(f64, -0x2.9fa8dcb9092a4p+8)));
    try std.testing.expectEqual(0x7.ffffffffff9c4p-972, exp(@as(f64, -0x2.9fa8dcb9092a6p+8)));
    try std.testing.expectEqual(0x8.0005c84b6996p-972, exp(@as(f64, -0x2.9fa8dcp+8)));
    try std.testing.expectEqual(0x7.ffe5c8744841p-972, exp(@as(f64, -0x2.9fa8ep+8)));
    try std.testing.expectEqual(0x8.00000000009c8p-972, exp(@as(f64, -0x2.9fa8dcb9092a4p+8)));
    try std.testing.expectEqual(0x7.ffffffffff9c4p-972, exp(@as(f64, -0x2.9fa8dcb9092a6p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c86p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c86p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5bd48p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5bd4cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5bd48bdc7c0cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5bd48bdc7c0ep+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5bd48p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5bd4cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5bd48bdc7c0cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5bd48bdc7c0ep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c86p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c86p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp(@as(f64, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f64, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x3.b629a25f7a4c8p-4, exp(@as(f64, -0x1.760cdp+0)));
    try std.testing.expectEqual(0x3.b6299af3270f2p-4, exp(@as(f64, -0x1.760cd2p+0)));
    try std.testing.expectEqual(0x3.b6299da019d6cp-4, exp(@as(f64, -0x1.760cd14774bd9p+0)));
    try std.testing.expectEqual(0x3.a823cf4b14606p+0, exp(@as(f64, 0x1.4bed28p+0)));
    try std.testing.expectEqual(0x3.8366d35e29fb4p-720, exp(@as(f64, -0x1.f1cf36p+8)));
    try std.testing.expectEqual(0x1.72a52c383a488p+4, exp(@as(f64, 0x3.248524p+0)));
    try std.testing.expectEqual(0x6.f5dcdfffff3ccp+0, exp(@as(f64, 0x1.f0b362p+0)));
    try std.testing.expectEqual(0xb.8c7b86075631p+16, exp(@as(f64, 0xd.89747p+0)));
    try std.testing.expectEqual(0xb.8c7acd3fa397p+16, exp(@as(f64, 0xd.89746p+0)));
    try std.testing.expectEqual(0xb.8c7b4638d385p+16, exp(@as(f64, 0xd.89746a799ac5p+0)));
    try std.testing.expectEqual(0xb.8c7b4638d37fp+16, exp(@as(f64, 0xd.89746a799ac48p+0)));
    try std.testing.expectEqual(0xa.c2d2580088708p-4, exp(@as(f64, -0x6.58b64p-4)));
    try std.testing.expectEqual(0x1.0001fefffffdep+0, exp(@as(f64, 0x1.fefe02p-16)));
    try std.testing.expectEqual(0x3.a84ddaee8cc56p-776, exp(@as(f64, -0x2.1895ep+8)));
    try std.testing.expectEqual(0x3.a83f39d46353p-776, exp(@as(f64, -0x2.1895e4p+8)));
    // try std.testing.expectEqual(0x3.a84196ad208cep-776, exp(@as(f64, -0x2.1895e35a9dc6cp+8)));

    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x2.b7e151628aed2a6cp+0, exp(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x7.63992e35376b731p+0, exp(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(0x1.415e5bf6fb105f2ep+4, exp(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(0x2.1df3b68cfb9ef7a8p+0, exp(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x1.19103e4080b45664p+72, exp(@as(f80, 0x3.2p+4)));
    try std.testing.expectEqual(0xf.ff6844410e1f547p+124, exp(@as(f80, 0x5.8b9028p+4)));
    try std.testing.expectEqual(0xf.7c2d08f39f969a2p+1020, exp(@as(f80, 0x2.c5cp+8)));
    try std.testing.expectEqual(0x6.79c8de6bb5ceb6p+1440, exp(@as(f80, 0x3.e8p+8)));
    try std.testing.expectEqual(0x1.3e21a464507f94ap+1024, exp(@as(f80, 0x2.c6p+8)));
    try std.testing.expectEqual(0xd.202c22e749b3087p-1784, exp(@as(f80, -0x4.d2p+8)));
    try std.testing.expectEqual(0x2.0004118603e6de38p+1024, exp(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0x1.fffc1185bdda0562p+1024, exp(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x2.000000000013ae58p+1024, exp(@as(f80, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0x1.ffffffffffd3ae5ap+1024, exp(@as(f80, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.fffffffffffffe5ap+1024, exp(@as(f80, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4.001236e259a704p-1076, exp(@as(f80, -0x2.e870a4p+8)));
    try std.testing.expectEqual(0x4.000236b97e84a93p-1076, exp(@as(f80, -0x2.e870a8p+8)));
    try std.testing.expectEqual(0x4.00029f178d98fa1p-1076, exp(@as(f80, -0x2.e870a7e5e88cp+8)));
    try std.testing.expectEqual(0x4.00029f178d18f9b8p-1076, exp(@as(f80, -0x2.e870a7e5e88c2p+8)));
    try std.testing.expectEqual(0x4.00029f178d1cc9cp-1076, exp(@as(f80, -0x2.e870a7e5e88c1f0cp+8)));
    try std.testing.expectEqual(0x4.00029f178d1cb9cp-1076, exp(@as(f80, -0x2.e870a7e5e88c1f1p+8)));
    try std.testing.expectEqual(0x2.0b9f4f64aed595b8p-1080, exp(@as(f80, -0x2.ebe224p+8)));
    try std.testing.expectEqual(0x2.0b9720f7ce27845p-1080, exp(@as(f80, -0x2.ebe228p+8)));
    try std.testing.expectEqual(0x2.0b981a509bab7298p-1080, exp(@as(f80, -0x2.ebe227861639p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0xd.be48e2532594eccp-16368, exp(@as(f80, -0x2.c4edp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c28p-16384, exp(@as(f80, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce9588p-16384, exp(@as(f80, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843accp+12)));
    try std.testing.expectEqual(0x3.eeb48df211be3694p-1024, exp(@as(f80, -0x2.c469d8p+8)));
    try std.testing.expectEqual(0x3.eea4d33f4f706d24p-1024, exp(@as(f80, -0x2.c469dcp+8)));
    try std.testing.expectEqual(0x3.eeb09f3f7b25c404p-1024, exp(@as(f80, -0x2.c469d9p+8)));
    try std.testing.expectEqual(0x3.e0206d364fe236d8p-1024, exp(@as(f80, -0x2.c46d94p+8)));
    try std.testing.expectEqual(0x3.e010ecd39be30a7cp-1024, exp(@as(f80, -0x2.c46d98p+8)));
    try std.testing.expectEqual(0x3.e018acfd35b146bp-1024, exp(@as(f80, -0x2.c46d96p+8)));
    try std.testing.expectEqual(0x3.f96438ed17abe504p-1024, exp(@as(f80, -0x2.c46724p+8)));
    try std.testing.expectEqual(0x3.f954537bfeeee9b8p-1024, exp(@as(f80, -0x2.c46728p+8)));
    try std.testing.expectEqual(0x3.f9584cd24f15bbf4p-1024, exp(@as(f80, -0x2.c46727p+8)));
    try std.testing.expectEqual(0x3.eea4d33f4f706d24p-1024, exp(@as(f80, -0x2.c469dcp+8)));
    try std.testing.expectEqual(0x3.ee9518cb776fd7fcp-1024, exp(@as(f80, -0x2.c469ep+8)));
    try std.testing.expectEqual(0x3.ee9cf5fd86363694p-1024, exp(@as(f80, -0x2.c469dep+8)));
    try std.testing.expectEqual(0x3.e6335d7047b5b9b8p-1024, exp(@as(f80, -0x2.c46c04p+8)));
    try std.testing.expectEqual(0x3.eab82516dd8695e8p-1024, exp(@as(f80, -0x2.c46adcp+8)));
    try std.testing.expectEqual(0x3.eaa87a559ec28104p-1024, exp(@as(f80, -0x2.c46aep+8)));
    try std.testing.expectEqual(0x3.eab04fae68c3ec14p-1024, exp(@as(f80, -0x2.c46adep+8)));
    try std.testing.expectEqual(0x3.d053f45176d6b17p-1024, exp(@as(f80, -0x2.c471bp+8)));
    try std.testing.expectEqual(0x3.d044b3202807caap-1024, exp(@as(f80, -0x2.c471b4p+8)));
    try std.testing.expectEqual(0x3.d0488366c34aeefp-1024, exp(@as(f80, -0x2.c471b3p+8)));
    try std.testing.expectEqual(0x3.efcfd88e9dc8fe0cp-1024, exp(@as(f80, -0x2.c4699p+8)));
    try std.testing.expectEqual(0x3.efc0196eb9e34d88p-1024, exp(@as(f80, -0x2.c46994p+8)));
    try std.testing.expectEqual(0x3.efc40930cb32bc18p-1024, exp(@as(f80, -0x2.c46993p+8)));
    try std.testing.expectEqual(0x3.2ff3a3a879b0f404p-1024, exp(@as(f80, -0x2.c49fap+8)));
    try std.testing.expectEqual(0x3.0941d1ff351b759cp-1024, exp(@as(f80, -0x2.c4ac1p+8)));
    try std.testing.expectEqual(0x2.8d3d2f65d464b28p-1024, exp(@as(f80, -0x2.c4d89p+8)));
    try std.testing.expectEqual(0x1.00400800aab555dep+0, exp(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0xf.fc007ff555fff77p-4, exp(@as(f80, -0x4p-12)));
    try std.testing.expectEqual(0x1.0000100000800002p+0, exp(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0xf.ffff000007ffffdp-4, exp(@as(f80, -0x1p-20)));
    try std.testing.expectEqual(0x1.0000000400000008p+0, exp(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffffc00000008p-4, exp(@as(f80, -0x4p-32)));
    try std.testing.expectEqual(0x1.0000000001p+0, exp(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0xf.fffffffffp-4, exp(@as(f80, -0x1p-40)));
    try std.testing.expectEqual(0x1.0000000000004p+0, exp(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffffcp-4, exp(@as(f80, -0x4p-52)));
    try std.testing.expectEqual(0x1.000000000000001p+0, exp(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0xf.ffffffffffffffp-4, exp(@as(f80, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x1p-10000)));
    try std.testing.expectEqual(0xf.fff8417363ff68fp+124, exp(@as(f80, 0x5.8b90b8p+4)));
    try std.testing.expectEqual(0x1.00000417184b8786p+128, exp(@as(f80, 0x5.8b90cp+4)));
    try std.testing.expectEqual(0x3.ffff2fe5259d01c8p-128, exp(@as(f80, -0x5.75628p+4)));
    try std.testing.expectEqual(0x4.00012fe53d8f8fe8p-128, exp(@as(f80, -0x5.756278p+4)));
    try std.testing.expectEqual(0x1.000020b8c430abdep+1024, exp(@as(f80, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f64p+1020, exp(@as(f80, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0xf.ffffffffff950d8p+1020, exp(@as(f80, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(0x1.000020b8c430abdep+1024, exp(@as(f80, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f64p+1020, exp(@as(f80, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0x1.00000000001950d8p+1024, exp(@as(f80, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0x4.000ebd79918d4c2p-1024, exp(@as(f80, -0x2.c46578p+8)));
    try std.testing.expectEqual(0x3.fffebd5e9bf2469cp-1024, exp(@as(f80, -0x2.c4657cp+8)));
    try std.testing.expectEqual(0x3.ffffffffff9eef4p-1024, exp(@as(f80, -0x2.c4657baf579a6p+8)));
    try std.testing.expectEqual(0x4.000ebd79918d4c2p-1024, exp(@as(f80, -0x2.c46578p+8)));
    try std.testing.expectEqual(0x3.fffebd5e9bf2469cp-1024, exp(@as(f80, -0x2.c4657cp+8)));
    try std.testing.expectEqual(0x4.00000000001eef4p-1024, exp(@as(f80, -0x2.c4657baf579a4p+8)));
    try std.testing.expectEqual(0x1.000020b8c430abdep+1024, exp(@as(f80, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f64p+1020, exp(@as(f80, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0x1.00000000001950d8p+1024, exp(@as(f80, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0xf.ffffffffff950d8p+1020, exp(@as(f80, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffcd8p+1020, exp(@as(f80, 0x2.c5c85fdf473de6acp+8)));
    try std.testing.expectEqual(0xf.ffffffffffff8d8p+1020, exp(@as(f80, 0x2.c5c85fdf473de6a8p+8)));
    try std.testing.expectEqual(0x1.000020b8c430abdep+1024, exp(@as(f80, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f64p+1020, exp(@as(f80, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0x1.00000000001950d8p+1024, exp(@as(f80, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0xf.ffffffffff950d8p+1020, exp(@as(f80, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffcd8p+1020, exp(@as(f80, 0x2.c5c85fdf473de6acp+8)));
    try std.testing.expectEqual(0xf.ffffffffffff8d8p+1020, exp(@as(f80, 0x2.c5c85fdf473de6a8p+8)));
    try std.testing.expectEqual(0x8.0005c84b6995d8ep-972, exp(@as(f80, -0x2.9fa8dcp+8)));
    try std.testing.expectEqual(0x7.ffe5c87448411fa8p-972, exp(@as(f80, -0x2.9fa8ep+8)));
    try std.testing.expectEqual(0x8.00000000009c5ap-972, exp(@as(f80, -0x2.9fa8dcb9092a4p+8)));
    try std.testing.expectEqual(0x7.ffffffffff9c59f8p-972, exp(@as(f80, -0x2.9fa8dcb9092a6p+8)));
    try std.testing.expectEqual(0x8.0000000000001ap-972, exp(@as(f80, -0x2.9fa8dcb9092a5388p+8)));
    try std.testing.expectEqual(0x7.fffffffffffff9f8p-972, exp(@as(f80, -0x2.9fa8dcb9092a538cp+8)));
    try std.testing.expectEqual(0x8.0005c84b6995d8ep-972, exp(@as(f80, -0x2.9fa8dcp+8)));
    try std.testing.expectEqual(0x7.ffe5c87448411fa8p-972, exp(@as(f80, -0x2.9fa8ep+8)));
    try std.testing.expectEqual(0x8.00000000009c5ap-972, exp(@as(f80, -0x2.9fa8dcb9092a4p+8)));
    try std.testing.expectEqual(0x7.ffffffffff9c59f8p-972, exp(@as(f80, -0x2.9fa8dcb9092a6p+8)));
    try std.testing.expectEqual(0x8.0000000000001ap-972, exp(@as(f80, -0x2.9fa8dcb9092a5388p+8)));
    try std.testing.expectEqual(0x7.fffffffffffff9f8p-972, exp(@as(f80, -0x2.9fa8dcb9092a538cp+8)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c86p+12)));
    try std.testing.expectEqual(0xf.fc2130abb1e4163p+16380, exp(@as(f80, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(0xf.fffffffff950d87p+16380, exp(@as(f80, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(0xf.fffffffffffcd87p+16380, exp(@as(f80, 0x2.c5c85fdf473de6acp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c86p+12)));
    try std.testing.expectEqual(0xf.fc2130abb1e4163p+16380, exp(@as(f80, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(0xf.fffffffff950d87p+16380, exp(@as(f80, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c85fdf473de6bp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c28p-16384, exp(@as(f80, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce9588p-16384, exp(@as(f80, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843accp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c28p-16384, exp(@as(f80, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce9588p-16384, exp(@as(f80, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x4.000000000000fc88p-16384, exp(@as(f80, -0x2.c5b2319c4843acbcp+12)));
    try std.testing.expectEqual(0x2.0017b984cbf1bfep-16384, exp(@as(f80, -0x2.c5bd48p+12)));
    try std.testing.expectEqual(0x1.ff97c395d33b18b8p-16384, exp(@as(f80, -0x2.c5bd4cp+12)));
    try std.testing.expectEqual(0x2.000000000136f198p-16384, exp(@as(f80, -0x2.c5bd48bdc7c0cp+12)));
    try std.testing.expectEqual(0x1.fffffffffd36f198p-16384, exp(@as(f80, -0x2.c5bd48bdc7c0ep+12)));
    try std.testing.expectEqual(0x1.fffffffffffff198p-16384, exp(@as(f80, -0x2.c5bd48bdc7c0c9b8p+12)));
    try std.testing.expectEqual(0x2.0017b984cbf1bfep-16384, exp(@as(f80, -0x2.c5bd48p+12)));
    try std.testing.expectEqual(0x1.ff97c395d33b18b8p-16384, exp(@as(f80, -0x2.c5bd4cp+12)));
    try std.testing.expectEqual(0x2.000000000136f198p-16384, exp(@as(f80, -0x2.c5bd48bdc7c0cp+12)));
    try std.testing.expectEqual(0x1.fffffffffd36f198p-16384, exp(@as(f80, -0x2.c5bd48bdc7c0ep+12)));
    try std.testing.expectEqual(0x2.0000000000007198p-16384, exp(@as(f80, -0x2.c5bd48bdc7c0c9b4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c86p+12)));
    try std.testing.expectEqual(0xf.fc2130abb1e4163p+16380, exp(@as(f80, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(0xf.fffffffff950d87p+16380, exp(@as(f80, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c85fdf473de6bp+12)));
    try std.testing.expectEqual(0xf.fffffffffffcd87p+16380, exp(@as(f80, 0x2.c5c85fdf473de6acp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c86p+12)));
    try std.testing.expectEqual(0xf.fc2130abb1e4163p+16380, exp(@as(f80, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(0xf.fffffffff950d87p+16380, exp(@as(f80, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp(@as(f80, 0x2.c5c85fdf473de6bp+12)));
    try std.testing.expectEqual(0xf.fffffffffffcd87p+16380, exp(@as(f80, 0x2.c5c85fdf473de6acp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c28p-16384, exp(@as(f80, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce9588p-16384, exp(@as(f80, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x4.000000000000fc88p-16384, exp(@as(f80, -0x2.c5b2319c4843acbcp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843accp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c28p-16384, exp(@as(f80, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce9588p-16384, exp(@as(f80, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x4.000000000000fc88p-16384, exp(@as(f80, -0x2.c5b2319c4843acbcp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffc88p-16384, exp(@as(f80, -0x2.c5b2319c4843accp+12)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x3.b629a25f7a4c82ap-4, exp(@as(f80, -0x1.760cdp+0)));
    try std.testing.expectEqual(0x3.b6299af3270f2ffcp-4, exp(@as(f80, -0x1.760cd2p+0)));
    try std.testing.expectEqual(0x3.b6299da019d6b33p-4, exp(@as(f80, -0x1.760cd14774bd9p+0)));
    try std.testing.expectEqual(0x3.a823cf4b14605f3cp+0, exp(@as(f80, 0x1.4bed28p+0)));
    try std.testing.expectEqual(0x3.8366d35e29fb32e8p-720, exp(@as(f80, -0x1.f1cf36p+8)));
    try std.testing.expectEqual(0x1.72a52c383a487ffep+4, exp(@as(f80, 0x3.248524p+0)));
    try std.testing.expectEqual(0x6.f5dcdfffff3ca048p+0, exp(@as(f80, 0x1.f0b362p+0)));
    try std.testing.expectEqual(0xb.8c7b8607563104p+16, exp(@as(f80, 0xd.89747p+0)));
    try std.testing.expectEqual(0xb.8c7acd3fa396cc4p+16, exp(@as(f80, 0xd.89746p+0)));
    try std.testing.expectEqual(0xb.8c7b4638d384df1p+16, exp(@as(f80, 0xd.89746a799ac5p+0)));
    try std.testing.expectEqual(0xb.8c7b4638d37f18dp+16, exp(@as(f80, 0xd.89746a799ac48p+0)));
    try std.testing.expectEqual(0xb.8c7b4638d384189p+16, exp(@as(f80, 0xd.89746a799ac4eedp+0)));
    try std.testing.expectEqual(0xa.c2d2580088709f3p-4, exp(@as(f80, -0x6.58b64p-4)));
    try std.testing.expectEqual(0x1.0001fefffffdd954p+0, exp(@as(f80, 0x1.fefe02p-16)));
    try std.testing.expectEqual(0x3.a84ddaee8cc5565cp-776, exp(@as(f80, -0x2.1895ep+8)));
    try std.testing.expectEqual(0x3.a83f39d46352f7a4p-776, exp(@as(f80, -0x2.1895e4p+8)));
    try std.testing.expectEqual(0x3.a84196ad208cefb4p-776, exp(@as(f80, -0x2.1895e35a9dc6cp+8)));

    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x2.b7e151628aed2a6abf7158809cf4p+0, exp(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x7.63992e35376b730ce8ee881ada2cp+0, exp(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(0x1.415e5bf6fb105f2d4bdfc53744c4p+4, exp(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(0x2.1df3b68cfb9ef7a986addc7dcee2p+0, exp(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x1.19103e4080b45664d7740cf8c5d9p+72, exp(@as(f128, 0x3.2p+4)));
    try std.testing.expectEqual(0xf.ff6844410e1f547369129d530d48p+124, exp(@as(f128, 0x5.8b9028p+4)));
    try std.testing.expectEqual(0xf.7c2d08f39f969a25d99164d121d8p+1020, exp(@as(f128, 0x2.c5cp+8)));
    try std.testing.expectEqual(0x6.79c8de6bb5ceb60158acfea8d148p+1440, exp(@as(f128, 0x3.e8p+8)));
    try std.testing.expectEqual(0x1.3e21a464507f94a0ae03700b899dp+1024, exp(@as(f128, 0x2.c6p+8)));
    try std.testing.expectEqual(0xd.202c22e749b30873a3228b398b5p-1784, exp(@as(f128, -0x4.d2p+8)));
    try std.testing.expectEqual(0x2.0004118603e6de38929bc069c7p+1024, exp(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0x1.fffc1185bdda0561d3753d38842fp+1024, exp(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x2.000000000013ae594e9bd9113664p+1024, exp(@as(f128, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0x1.ffffffffffd3ae594e9bda9b6b3bp+1024, exp(@as(f128, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.fffffffffffffe594e9bd8b06065p+1024, exp(@as(f128, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4.001236e259a703fe35e4572c21b8p-1076, exp(@as(f128, -0x2.e870a4p+8)));
    try std.testing.expectEqual(0x4.000236b97e84a932aa555f5e8ce4p-1076, exp(@as(f128, -0x2.e870a8p+8)));
    try std.testing.expectEqual(0x4.00029f178d98fa0e72d456f74c28p-1076, exp(@as(f128, -0x2.e870a7e5e88cp+8)));
    try std.testing.expectEqual(0x4.00029f178d18f9ba8fe2abd80f98p-1076, exp(@as(f128, -0x2.e870a7e5e88c2p+8)));
    try std.testing.expectEqual(0x4.00029f178d1cc9bd0f851e55aecp-1076, exp(@as(f128, -0x2.e870a7e5e88c1f0cp+8)));
    try std.testing.expectEqual(0x4.00029f178d1cb9bd0508c0213bb8p-1076, exp(@as(f128, -0x2.e870a7e5e88c1f1p+8)));
    try std.testing.expectEqual(0x4.00029f178d1cbba1a34fc1f1ad4cp-1076, exp(@as(f128, -0x2.e870a7e5e88c1f0f86d8bda5cef2p+8)));
    try std.testing.expectEqual(0x4.00029f178d1cbba1a34fc1f1a54cp-1076, exp(@as(f128, -0x2.e870a7e5e88c1f0f86d8bda5cef4p+8)));
    try std.testing.expectEqual(0x4.00029f178d1cbba1a34fc1f5755p-1076, exp(@as(f128, -0x2.e870a7e5e88c1f0f86d8bda5cep+8)));
    try std.testing.expectEqual(0x4.00029f178d1cbba1a34fc1f1754cp-1076, exp(@as(f128, -0x2.e870a7e5e88c1f0f86d8bda5cfp+8)));
    try std.testing.expectEqual(0x2.0b9f4f64aed595b7b1e41fe97b6ep-1080, exp(@as(f128, -0x2.ebe224p+8)));
    try std.testing.expectEqual(0x2.0b9720f7ce27844ea9674284d868p-1080, exp(@as(f128, -0x2.ebe228p+8)));
    try std.testing.expectEqual(0x2.0b981a509bab72997118df69cd8p-1080, exp(@as(f128, -0x2.ebe227861639p+8)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, exp(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0xd.be48e2532594ecc1a3b8f7ce2038p-16368, exp(@as(f128, -0x2.c4edp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c27d0a36c18105cp-16384, exp(@as(f128, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce958ba803f3e779bp-16384, exp(@as(f128, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc85647bac501164p-16384, exp(@as(f128, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc85647d4c57069cp-16384, exp(@as(f128, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffc85647a6732d718p-16384, exp(@as(f128, -0x2.c5b2319c4843accp+12)));
    try std.testing.expectEqual(0x3.eeb48df211be3692ba54c47bcc8p-1024, exp(@as(f128, -0x2.c469d8p+8)));
    try std.testing.expectEqual(0x3.eea4d33f4f706d23cb49d8c61f8cp-1024, exp(@as(f128, -0x2.c469dcp+8)));
    try std.testing.expectEqual(0x3.eeb09f3f7b25c404f4404aac7ae8p-1024, exp(@as(f128, -0x2.c469d9p+8)));
    try std.testing.expectEqual(0x3.e0206d364fe236d9c5d77c55c0cep-1024, exp(@as(f128, -0x2.c46d94p+8)));
    try std.testing.expectEqual(0x3.e010ecd39be30a7d887797e1f9dcp-1024, exp(@as(f128, -0x2.c46d98p+8)));
    try std.testing.expectEqual(0x3.e018acfd35b146aea65f341569fp-1024, exp(@as(f128, -0x2.c46d96p+8)));
    try std.testing.expectEqual(0x3.f96438ed17abe50216a8a6561002p-1024, exp(@as(f128, -0x2.c46724p+8)));
    try std.testing.expectEqual(0x3.f954537bfeeee9b89a29977f429ep-1024, exp(@as(f128, -0x2.c46728p+8)));
    try std.testing.expectEqual(0x3.f9584cd24f15bbf3bba7324c725ap-1024, exp(@as(f128, -0x2.c46727p+8)));
    try std.testing.expectEqual(0x3.eea4d33f4f706d23cb49d8c61f8cp-1024, exp(@as(f128, -0x2.c469dcp+8)));
    try std.testing.expectEqual(0x3.ee9518cb776fd7fdb6575a15d89cp-1024, exp(@as(f128, -0x2.c469ep+8)));
    try std.testing.expectEqual(0x3.ee9cf5fd863636931550dd9e73d4p-1024, exp(@as(f128, -0x2.c469dep+8)));
    try std.testing.expectEqual(0x3.e6335d7047b5b9b96c422bdab03ap-1024, exp(@as(f128, -0x2.c46c04p+8)));
    try std.testing.expectEqual(0x3.eab82516dd8695e83929b0026a4ap-1024, exp(@as(f128, -0x2.c46adcp+8)));
    try std.testing.expectEqual(0x3.eaa87a559ec281025eff0982481ap-1024, exp(@as(f128, -0x2.c46aep+8)));
    try std.testing.expectEqual(0x3.eab04fae68c3ec15de16fa21fe0cp-1024, exp(@as(f128, -0x2.c46adep+8)));
    try std.testing.expectEqual(0x3.d053f45176d6b170998260fa99e4p-1024, exp(@as(f128, -0x2.c471bp+8)));
    try std.testing.expectEqual(0x3.d044b3202807caa00a6f8afedc28p-1024, exp(@as(f128, -0x2.c471b4p+8)));
    try std.testing.expectEqual(0x3.d0488366c34aeeee2f9bbbae7086p-1024, exp(@as(f128, -0x2.c471b3p+8)));
    try std.testing.expectEqual(0x3.efcfd88e9dc8fe0b263c320689ep-1024, exp(@as(f128, -0x2.c4699p+8)));
    try std.testing.expectEqual(0x3.efc0196eb9e34d87beedb9c2a6c2p-1024, exp(@as(f128, -0x2.c46994p+8)));
    try std.testing.expectEqual(0x3.efc40930cb32bc17ecf6f28f4976p-1024, exp(@as(f128, -0x2.c46993p+8)));
    try std.testing.expectEqual(0x3.2ff3a3a879b0f4051f06a63de3f2p-1024, exp(@as(f128, -0x2.c49fap+8)));
    try std.testing.expectEqual(0x3.0941d1ff351b759b1317765f7692p-1024, exp(@as(f128, -0x2.c4ac1p+8)));
    try std.testing.expectEqual(0x2.8d3d2f65d464b281f702d1fadab4p-1024, exp(@as(f128, -0x2.c4d89p+8)));
    try std.testing.expectEqual(0x1.00400800aab555dde38e6ce86e92p+0, exp(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0xf.fc007ff555fff777d279e7b87adp-4, exp(@as(f128, -0x4p-12)));
    try std.testing.expectEqual(0x1.0000100000800002aaaab5555577p+0, exp(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0xf.ffff000007ffffd55555fffffdep-4, exp(@as(f128, -0x1p-20)));
    try std.testing.expectEqual(0x1.00000004000000080000000aaaabp+0, exp(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffffc00000007fffffff555558p-4, exp(@as(f128, -0x4p-32)));
    try std.testing.expectEqual(0x1.000000000100000000008p+0, exp(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0xf.fffffffff00000000008p-4, exp(@as(f128, -0x1p-40)));
    try std.testing.expectEqual(0x1.00000000000040000000000008p+0, exp(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffffc0000000000008p-4, exp(@as(f128, -0x4p-52)));
    try std.testing.expectEqual(0x1.000000000000001p+0, exp(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0xf.ffffffffffffffp-4, exp(@as(f128, -0x1p-60)));
    try std.testing.expectEqual(0x1.0000000000000000000000001p+0, exp(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffp-4, exp(@as(f128, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x1p-10000)));
    try std.testing.expectEqual(0xf.fff8417363ff68e9feae7fbf0fbp+124, exp(@as(f128, 0x5.8b90b8p+4)));
    try std.testing.expectEqual(0x1.00000417184b8786d5655cf85948p+128, exp(@as(f128, 0x5.8b90cp+4)));
    try std.testing.expectEqual(0x3.ffff2fe5259d01c9284638c92dd4p-128, exp(@as(f128, -0x5.75628p+4)));
    try std.testing.expectEqual(0x4.00012fe53d8f8fe9a18888523c4cp-128, exp(@as(f128, -0x5.756278p+4)));
    try std.testing.expectEqual(0x1.000020b8c430abde3653c03623p+1024, exp(@as(f128, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f6404b50edab3e208p+1020, exp(@as(f128, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0xf.ffffffffff950d87131a0068afep+1020, exp(@as(f128, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(0x1.000020b8c430abde3653c03623p+1024, exp(@as(f128, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f6404b50edab3e208p+1020, exp(@as(f128, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0x1.00000000001950d87131a130a60cp+1024, exp(@as(f128, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0x4.000ebd79918d4c20f5b15cc47ba8p-1024, exp(@as(f128, -0x2.c46578p+8)));
    try std.testing.expectEqual(0x3.fffebd5e9bf2469b17acee55c738p-1024, exp(@as(f128, -0x2.c4657cp+8)));
    try std.testing.expectEqual(0x3.ffffffffff9eef3fec1be818c628p-1024, exp(@as(f128, -0x2.c4657baf579a6p+8)));
    try std.testing.expectEqual(0x4.000ebd79918d4c20f5b15cc47ba8p-1024, exp(@as(f128, -0x2.c46578p+8)));
    try std.testing.expectEqual(0x3.fffebd5e9bf2469b17acee55c738p-1024, exp(@as(f128, -0x2.c4657cp+8)));
    try std.testing.expectEqual(0x4.00000000001eef3fec1be3f6ae24p-1024, exp(@as(f128, -0x2.c4657baf579a4p+8)));
    try std.testing.expectEqual(0x1.000020b8c430abde3653c03623p+1024, exp(@as(f128, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f6404b50edab3e208p+1020, exp(@as(f128, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0x1.00000000001950d87131a130a60cp+1024, exp(@as(f128, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0xf.ffffffffff950d87131a0068afep+1020, exp(@as(f128, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffcd871319ff0342ap+1020, exp(@as(f128, 0x2.c5c85fdf473de6acp+8)));
    try std.testing.expectEqual(0xf.ffffffffffff8d871319ff0343fp+1020, exp(@as(f128, 0x2.c5c85fdf473de6a8p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffbffffffffff342d8p+1020, exp(@as(f128, 0x2.c5c85fdf473de6ab278ece600fp+8)));
    try std.testing.expectEqual(0x1.000020b8c430abde3653c03623p+1024, exp(@as(f128, 0x2.c5c86p+8)));
    try std.testing.expectEqual(0xf.ffc20c04143f6404b50edab3e208p+1020, exp(@as(f128, 0x2.c5c85cp+8)));
    try std.testing.expectEqual(0x1.00000000001950d87131a130a60cp+1024, exp(@as(f128, 0x2.c5c85fdf473ep+8)));
    try std.testing.expectEqual(0xf.ffffffffff950d87131a0068afep+1020, exp(@as(f128, 0x2.c5c85fdf473dep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffcd871319ff0342ap+1020, exp(@as(f128, 0x2.c5c85fdf473de6acp+8)));
    try std.testing.expectEqual(0xf.ffffffffffff8d871319ff0343fp+1020, exp(@as(f128, 0x2.c5c85fdf473de6a8p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffc0000000000342d8p+1020, exp(@as(f128, 0x2.c5c85fdf473de6ab278ece601p+8)));
    try std.testing.expectEqual(0x8.0005c84b6995d8e0ebcf71f99d9p-972, exp(@as(f128, -0x2.9fa8dcp+8)));
    try std.testing.expectEqual(0x7.ffe5c87448411fa72724e261227p-972, exp(@as(f128, -0x2.9fa8ep+8)));
    try std.testing.expectEqual(0x8.00000000009c59f97716592f71b8p-972, exp(@as(f128, -0x2.9fa8dcb9092a4p+8)));
    try std.testing.expectEqual(0x7.ffffffffff9c59f9771655a43288p-972, exp(@as(f128, -0x2.9fa8dcb9092a6p+8)));
    try std.testing.expectEqual(0x8.00000000000019f9771653379568p-972, exp(@as(f128, -0x2.9fa8dcb9092a5388p+8)));
    try std.testing.expectEqual(0x7.fffffffffffff9f977165337954p-972, exp(@as(f128, -0x2.9fa8dcb9092a538cp+8)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffff953cp-972, exp(@as(f128, -0x2.9fa8dcb9092a538b3f2ee2ca67p+8)));
    try std.testing.expectEqual(0x8.0005c84b6995d8e0ebcf71f99d9p-972, exp(@as(f128, -0x2.9fa8dcp+8)));
    try std.testing.expectEqual(0x7.ffe5c87448411fa72724e261227p-972, exp(@as(f128, -0x2.9fa8ep+8)));
    try std.testing.expectEqual(0x8.00000000009c59f97716592f71b8p-972, exp(@as(f128, -0x2.9fa8dcb9092a4p+8)));
    try std.testing.expectEqual(0x7.ffffffffff9c59f9771655a43288p-972, exp(@as(f128, -0x2.9fa8dcb9092a6p+8)));
    try std.testing.expectEqual(0x8.00000000000019f9771653379568p-972, exp(@as(f128, -0x2.9fa8dcb9092a5388p+8)));
    try std.testing.expectEqual(0x7.fffffffffffff9f977165337954p-972, exp(@as(f128, -0x2.9fa8dcb9092a538cp+8)));
    try std.testing.expectEqual(0x8.000000000000000000000007954p-972, exp(@as(f128, -0x2.9fa8dcb9092a538b3f2ee2ca66p+8)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c86p+12)));
    try std.testing.expectEqual(0xf.fc2130abb1e4162a1031a29a3918p+16380, exp(@as(f128, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(0xf.fffffffff950d87131a155a1b3a8p+16380, exp(@as(f128, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(0xf.fffffffffffcd871319ff03474ep+16380, exp(@as(f128, 0x2.c5c85fdf473de6acp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c86p+12)));
    try std.testing.expectEqual(0xf.fc2130abb1e4162a1031a29a3918p+16380, exp(@as(f128, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(0xf.fffffffff950d87131a155a1b3a8p+16380, exp(@as(f128, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473de6bp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c27d0a36c18105cp-16384, exp(@as(f128, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce958ba803f3e779bp-16384, exp(@as(f128, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc85647bac501164p-16384, exp(@as(f128, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc85647d4c57069cp-16384, exp(@as(f128, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffc85647a6732d718p-16384, exp(@as(f128, -0x2.c5b2319c4843accp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c27d0a36c18105cp-16384, exp(@as(f128, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce958ba803f3e779bp-16384, exp(@as(f128, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc85647bac501164p-16384, exp(@as(f128, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc85647d4c57069cp-16384, exp(@as(f128, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x4.000000000000fc85647a6732f63cp-16384, exp(@as(f128, -0x2.c5b2319c4843acbcp+12)));
    // try std.testing.expectEqual(0x2.0017b984cbf1bfdc86bdf5ca5c0cp-16384, exp(@as(f128, -0x2.c5bd48p+12)));
    try std.testing.expectEqual(0x1.ff97c395d33b18b4e8d5b82a339p-16384, exp(@as(f128, -0x2.c5bd4cp+12)));
    try std.testing.expectEqual(0x2.000000000136f19a4604f934f4dcp-16384, exp(@as(f128, -0x2.c5bd48bdc7c0cp+12)));
    try std.testing.expectEqual(0x1.fffffffffd36f19a46068b51c05p-16384, exp(@as(f128, -0x2.c5bd48bdc7c0ep+12)));
    // try std.testing.expectEqual(0x1.fffffffffffff19a46049ac973a4p-16384, exp(@as(f128, -0x2.c5bd48bdc7c0c9b8p+12)));
    // try std.testing.expectEqual(0x2.0017b984cbf1bfdc86bdf5ca5c0cp-16384, exp(@as(f128, -0x2.c5bd48p+12)));
    try std.testing.expectEqual(0x1.ff97c395d33b18b4e8d5b82a339p-16384, exp(@as(f128, -0x2.c5bd4cp+12)));
    try std.testing.expectEqual(0x2.000000000136f19a4604f934f4dcp-16384, exp(@as(f128, -0x2.c5bd48bdc7c0cp+12)));
    try std.testing.expectEqual(0x1.fffffffffd36f19a46068b51c05p-16384, exp(@as(f128, -0x2.c5bd48bdc7c0ep+12)));
    try std.testing.expectEqual(0x2.000000000000719a46049ac9800cp-16384, exp(@as(f128, -0x2.c5bd48bdc7c0c9b4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c86p+12)));
    try std.testing.expectEqual(0xf.fc2130abb1e4162a1031a29a3918p+16380, exp(@as(f128, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(0xf.fffffffff950d87131a155a1b3a8p+16380, exp(@as(f128, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473de6bp+12)));
    try std.testing.expectEqual(0xf.fffffffffffcd871319ff03474ep+16380, exp(@as(f128, 0x2.c5c85fdf473de6acp+12)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffe254p+16380, exp(@as(f128, 0x2.c5c85fdf473de6af278ece600fcap+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473de6af278ece601p+12)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffff34254p+16380, exp(@as(f128, 0x2.c5c85fdf473de6af278ece600fp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c86p+12)));
    try std.testing.expectEqual(0xf.fc2130abb1e4162a1031a29a3918p+16380, exp(@as(f128, 0x2.c5c85cp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473ep+12)));
    try std.testing.expectEqual(0xf.fffffffff950d87131a155a1b3a8p+16380, exp(@as(f128, 0x2.c5c85fdf473dep+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473de6bp+12)));
    try std.testing.expectEqual(0xf.fffffffffffcd871319ff03474ep+16380, exp(@as(f128, 0x2.c5c85fdf473de6acp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473de6af278ece600fccp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp(@as(f128, 0x2.c5c85fdf473de6af278ece601p+12)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffff34254p+16380, exp(@as(f128, 0x2.c5c85fdf473de6af278ece600fp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c27d0a36c18105cp-16384, exp(@as(f128, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce958ba803f3e779bp-16384, exp(@as(f128, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc85647bac501164p-16384, exp(@as(f128, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc85647d4c57069cp-16384, exp(@as(f128, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x4.000000000000fc85647a6732f63cp-16384, exp(@as(f128, -0x2.c5b2319c4843acbcp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffc85647a6732d718p-16384, exp(@as(f128, -0x2.c5b2319c4843accp+12)));
    try std.testing.expectEqual(0x3.ffffffffffffffffffffffffd718p-16384, exp(@as(f128, -0x2.c5b2319c4843acbff21591e99cccp+12)));
    try std.testing.expectEqual(0x4.000000000000000000000032d718p-16384, exp(@as(f128, -0x2.c5b2319c4843acbff21591e99cp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffffffffffff2d718p-16384, exp(@as(f128, -0x2.c5b2319c4843acbff21591e99dp+12)));
    try std.testing.expectEqual(0x4.00671741091b8c27d0a36c18105cp-16384, exp(@as(f128, -0x2.c5b23p+12)));
    try std.testing.expectEqual(0x3.ff671d7bc6ce958ba803f3e779bp-16384, exp(@as(f128, -0x2.c5b234p+12)));
    try std.testing.expectEqual(0x4.00000000032ffc85647bac501164p-16384, exp(@as(f128, -0x2.c5b2319c4843ap+12)));
    try std.testing.expectEqual(0x3.fffffffffb2ffc85647d4c57069cp-16384, exp(@as(f128, -0x2.c5b2319c4843cp+12)));
    try std.testing.expectEqual(0x4.000000000000fc85647a6732f63cp-16384, exp(@as(f128, -0x2.c5b2319c4843acbcp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffc85647a6732d718p-16384, exp(@as(f128, -0x2.c5b2319c4843accp+12)));
    try std.testing.expectEqual(0x4.0000000000000000000000005718p-16384, exp(@as(f128, -0x2.c5b2319c4843acbff21591e99ccap+12)));
    try std.testing.expectEqual(0x4.000000000000000000000032d718p-16384, exp(@as(f128, -0x2.c5b2319c4843acbff21591e99cp+12)));
    try std.testing.expectEqual(0x3.fffffffffffffffffffffff2d718p-16384, exp(@as(f128, -0x2.c5b2319c4843acbff21591e99dp+12)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, exp(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.b629a25f7a4c829f1c5a4743e518p-4, exp(@as(f128, -0x1.760cdp+0)));
    try std.testing.expectEqual(0x3.b6299af3270f2ffdc3211b78da96p-4, exp(@as(f128, -0x1.760cd2p+0)));
    try std.testing.expectEqual(0x3.b6299da019d6b3318f1b14156892p-4, exp(@as(f128, -0x1.760cd14774bd9p+0)));
    try std.testing.expectEqual(0x3.a823cf4b14605f3bc47d07a323p+0, exp(@as(f128, 0x1.4bed28p+0)));
    try std.testing.expectEqual(0x3.8366d35e29fb32e7aebc60e8c72cp-720, exp(@as(f128, -0x1.f1cf36p+8)));
    try std.testing.expectEqual(0x1.72a52c383a487ffd4852b11d665p+4, exp(@as(f128, 0x3.248524p+0)));
    try std.testing.expectEqual(0x6.f5dcdfffff3ca04a93d557465a7p+0, exp(@as(f128, 0x1.f0b362p+0)));
    try std.testing.expectEqual(0xb.8c7b8607563103fb90d5d7b02058p+16, exp(@as(f128, 0xd.89747p+0)));
    try std.testing.expectEqual(0xb.8c7acd3fa396cc3cb84da240fed8p+16, exp(@as(f128, 0xd.89746p+0)));
    try std.testing.expectEqual(0xb.8c7b4638d384df0dd0342350a048p+16, exp(@as(f128, 0xd.89746a799ac5p+0)));
    try std.testing.expectEqual(0xb.8c7b4638d37f18d02d17b98fa25p+16, exp(@as(f128, 0xd.89746a799ac48p+0)));
    try std.testing.expectEqual(0xb.8c7b4638d384188f496d329d2b38p+16, exp(@as(f128, 0xd.89746a799ac4eedp+0)));
    try std.testing.expectEqual(0xa.c2d2580088709f3262612e0cb808p-4, exp(@as(f128, -0x6.58b64p-4)));
    try std.testing.expectEqual(0x1.0001fefffffdd953027f9648617p+0, exp(@as(f128, 0x1.fefe02p-16)));
    try std.testing.expectEqual(0x3.a84ddaee8cc5565aff6debc74164p-776, exp(@as(f128, -0x2.1895ep+8)));
    try std.testing.expectEqual(0x3.a83f39d46352f7a29ab0ace8657ap-776, exp(@as(f128, -0x2.1895e4p+8)));
    try std.testing.expectEqual(0x3.a84196ad208cefb3b0c82f6b5b38p-776, exp(@as(f128, -0x2.1895e35a9dc6cp+8)));
}
