const std = @import("std");
const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const exp2_data = @import("exp2_data.zig");
const exp_data = @import("exp_data.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn exp(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.exp: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return scast(f16, exp32(scast(f32, x))),
        f32 => {
            // glibc/sysdeps/ieee754/flt-32/e_expf.c
            return exp32(scast(f32, x));
        },
        f64 => {
            // glibc/sysdeps/ieee754/dbl-64/e_exp.c
            return exp64(scast(f64, x));
        },
        f80 => return scast(f80, exp128(scast(f128, x))),
        f128 => {
            // glibc/sysdeps/ieee754/ldbl-128/e_expl.c
            return exp128(scast(f128, x));
        },
        else => unreachable,
    }
}

inline fn top12(x: f32) u32 {
    return @as(u32, @bitCast(x)) >> 20;
}

fn exp32(x: f32) f32 {
    const xd: f64 = scast(f64, x);
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
    return scast(f32, y);
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
        const tval1: i32 = scast(i32, t * TWO8);

        xx -= exp_data.__expl_table[@intCast(2 * 89 + 2 * tval1)];
        xl -= exp_data.__expl_table[@intCast(2 * 89 + 2 * tval1 + 1)];

        // Calculate t/32768.
        t = xx + THREEp96;
        t -= THREEp96;

        // Compute tval2 = t.
        const tval2: i32 = scast(i32, t * TWO15);

        xx -= exp_data.__expl_table[@intCast((2 * (2 * 89) + 2 + 2 * 65) + 2 * tval2)];
        xl -= exp_data.__expl_table[@intCast((2 * (2 * 89) + 2 + 2 * 65) + 2 * tval2 + 1)];

        xx = xx + xl;

        // Compute ex2 = 2^n_0 e^(argtable[tval1]) e^(argtable[tval2]).
        var ex2_u: ldbl128.ieee_f128_shape = @bitCast(exp_data.__expl_table[@intCast(((2 * (2 * 89) + 2 + 2 * 65) + 2 + 2 * 65 + 89) + tval1)] * exp_data.__expl_table[@intCast((((2 * (2 * 89) + 2 + 2 * 65) + 2 + 2 * 65 + 89) + 1 + 89 + 65) + tval2)]);
        const n_i: i32 = scast(i32, n);
        // 'unsafe' is 1 iff n_1 != 0.
        const unsafe: i32 = scast(i32, int.abs(n_i) >= 15000);
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
