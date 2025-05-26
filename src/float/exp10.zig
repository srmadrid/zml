const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const exp2_data = @import("exp2_data.zig");
const exp_data = @import("exp_data.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn exp10(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.exp10: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return scast(f16, exp10_32(scast(f32, x))),
        f32 => {
            // glibc/sysdeps/ieee754/flt-32/e_exp10f.c
            return exp10_32(scast(f32, x));
        },
        f64 => {
            // glibc/sysdeps/ieee754/dbl-64/e_exp10.c
            return exp10_64(scast(f64, x));
        },
        f80 => return scast(f80, exp10_128(scast(f128, x))),
        f128 => {
            // glibc/sysdeps/ieee754/ldbl-128/e_exp10l.c
            return exp10_128(scast(f128, x));
        },
        else => unreachable,
    }
}

inline fn top13(x: f32) u32 {
    return @as(u32, @bitCast(x)) >> 19;
}

fn exp10_32(x: f32) f32 {
    const xd: f64 = scast(f64, x);
    const abstop: u32 = top13(x) & 0xfff; // Ignore sign.
    if (abstop >= top13(38)) {
        @branchHint(.unlikely);
        // |x| >= 38 or x is nan.
        if (@as(u32, @bitCast(x)) == @as(u32, @bitCast(-std.math.inf(f32))))
            return 0;

        if (abstop >= top13(std.math.inf(f32)))
            return x + x;

        // 0x26.8826ap0 is the largest value such that 10^x < 2^128.
        if (x > 0x26.8826ap0)
            return 0x1p97 * 0x1p97;

        // -0x2d.278d4p0 is the smallest value such that 10^x > 2^-150.
        if (x < -0x2d.278d4p0)
            return 0x1p-95 * 0x1p-95;

        if (x < -0x2c.da7cfp0)
            return 0x1.4p-75 * 0x1.4p-75;

        // the smallest value such that 10^x >= 2^-126 (normal range)
        // is x = -0x25.ee060p0
        // we go through here for 2014929 values out of 2060451840
        // (not counting NaN and infinities, i.e., about 0.1%
    }

    // x*N*Ln10/Ln2 = k + r with r in [-1/2, 1/2] and int k.
    var z: f64 = (0x3.5269e12f346e2p0 * 32) * xd;
    // |xd| < 38 thus |z| < 1216
    var kd: f64 = z + exp2_data.SHIFT_32;
    const ki: u64 = @bitCast(kd);
    kd -= exp2_data.SHIFT_32;
    const r: f64 = z - kd;

    // 10^x = 10^(k/N) * 10^(r/N) ~= s * (C0*r^3 + C1*r^2 + C2*r + 1)
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

fn special_case(sbits: u64, tmp: f64, ki: u64) f64 {
    if ((ki & 0x80000000) == 0) {
        // The exponent of scale might have overflowed by 1.
        const scale: f64 = @bitCast(sbits - (1 << 52));
        const y: f64 = 2 * (scale + scale * tmp);
        return y;
    }

    // n < 0, need special care in the subnormal range.
    const scale: f64 = @bitCast(sbits + (1022 << 52));
    var y: f64 = scale + scale * tmp;

    if (y < 1) {
        // Round y to the right precision before scaling it into the subnormal
        // range to avoid double rounding that can cause 0.5+E/2 ulp error where
        // E is the worst-case ulp error outside the subnormal range.  So this
        // is only useful if the goal is better than 1 ulp worst-case error.
        var lo: f64 = scale - y + scale * tmp;
        const hi: f64 = 1 + y;
        lo = 1 - hi + y + lo;
        y = (hi + lo) - 1;
        // Avoid -0.0 with downward rounding.
        if (y == 0)
            y = 0;

        // The underflow exception needs to be signaled explicitly.
        std.mem.doNotOptimizeAway(0x1p-1022 * 0x1p-1022);
    }

    y = 0x1p-1022 * y;

    return y;
}

// Double-precision 10^x approximation. Largest observed error is ~0.513 ULP.
fn exp10_64(x: f64) f64 {
    const ix: u64 = @bitCast(x);
    var abstop: u32 = @truncate((ix >> 52) & 0x7ff);
    if (abstop -% 0x3c6 >= 0x41) {
        @branchHint(.unlikely);
        if (abstop -% 0x3c6 >= 0x80000000) {
            // Avoid spurious underflow for tiny x.
            // Note: 0 is common input.
            return x + 1;
        }

        if (abstop == 0x7ff)
            return if (ix == @as(u64, @bitCast(-std.math.inf(f64)))) 0 else x + 1;

        if (x >= 0x1.34413509f79ffp8)
            return 0x1p769 * 0x1p769;

        if (x < -0x1.5ep+8)
            return 0x1p-767 * 0x1p-767;

        // Large x is special-cased below.
        abstop = 0;
    }

    // Reduce x: z = x * N / log10(2), k = round(z).
    const z: f64 = exp_data.invlog10_2N * x;
    var kd: f64 = z + exp_data.Shift_64;
    const ki: u64 = @bitCast(kd);
    kd -= exp_data.Shift_64;

    // r = x - k * log10(2), r in [-0.5, 0.5].
    var r: f64 = x;
    r = exp_data.neglog10_2hiN * kd + r;
    r = exp_data.neglog10_2loN * kd + r;

    // exp10(x) = 2^(k/N) * 2^(r/N).
    // Approximate the two components separately.

    // s = 2^(k/N), using lookup table.
    const e: u64 = ki << 45;
    const i: u64 = (ki & 127) * 2;
    const u: u64 = exp_data.T_64[i + 1];
    const sbits: u64 = u +% e;

    const tail: f64 = @bitCast(exp_data.T_64[i]);

    // 2^(r/N) ~= 1 + r * Poly(r).
    const r2: f64 = r * r;
    const p: f64 = exp_data.exp10_poly[0] + r * exp_data.exp10_poly[1];
    var y: f64 = exp_data.exp10_poly[2] + r * exp_data.exp10_poly[3];
    y = y + r2 * exp_data.exp10_poly[4];
    y = p + r2 * y;
    y = tail + y * r;

    if (abstop == 0) {
        @branchHint(.unlikely);
        return special_case(sbits, y, ki);
    }

    // Assemble components:
    // y  = 2^(r/N) * 2^(k/N)
    //   ~= (y + 1) * s.  */
    const s: f64 = @bitCast(sbits);
    return s * y + s;
}

fn exp10_128(arg: f128) f128 {
    const log10_high: f128 = 0x2.4d763776aaa2bp0;
    const log10_low: f128 = 0x5.ba95b58ae0b4c28a38a3fb3e7698p-60;

    if (!std.math.isFinite(arg))
        return float.exp(arg);

    if (arg < -4931 - 18 - 10) {
        return 3.36210314311209350626e-4932 * 3.36210314311209350626e-4932;
    } else if (arg > 4932 + 1) {
        return 1.18973149535723176502e+4932 * 1.18973149535723176502e+4932;
    } else if (float.abs(arg) < 0x1p-116) {
        return 1;
    }

    var u: ldbl128.ieee_f128_shape64 = @bitCast(arg);
    u.lsw &= 0xfe00000000000000;
    const arg_high: f128 = @bitCast(u);
    const arg_low: f128 = arg - arg_high;
    const exp_high: f128 = arg_high * log10_high;
    const exp_low: f128 = arg_high * log10_low + arg_low * 2.302585092994045684017991454684364208;
    return float.exp(exp_high) * float.exp(exp_low);
}
