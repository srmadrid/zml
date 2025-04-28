const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const exp2_data = @import("exp2_data.zig");
const exp_data = @import("exp_data.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub fn exp10(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return exp10(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, exp10_32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_exp10f.c
                    return exp10_32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_exp10.c
                    return exp10_64(x);
                },
                f80 => return cast(f80, exp10_128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_exp10l.c
                    return exp10_128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

inline fn top13(x: f32) u32 {
    return @as(u32, @bitCast(x)) >> 19;
}

fn exp10_32(x: f32) f32 {
    const xd: f64 = cast(f64, x, .{});
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
    return cast(f32, y, .{});
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
        return math.exp(arg);

    if (arg < -4931 - 18 - 10) {
        return 3.36210314311209350626e-4932 * 3.36210314311209350626e-4932;
    } else if (arg > 4932 + 1) {
        return 1.18973149535723176502e+4932 * 1.18973149535723176502e+4932;
    } else if (math.abs(arg) < 0x1p-116) {
        return 1;
    }

    var u: ldbl128.ieee_f128_shape64 = @bitCast(arg);
    u.lsw &= 0xfe00000000000000;
    const arg_high: f128 = @bitCast(u);
    const arg_low: f128 = arg - arg_high;
    const exp_high: f128 = arg_high * log10_high;
    const exp_low: f128 = arg_high * log10_low + arg_low * 2.302585092994045684017991454684364208;
    return math.exp(exp_high) * math.exp(exp_low);
}

test exp10 {
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x3.e8p+8, exp10(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(0x1.99999ap-4, exp10(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0xc.097cep+116, exp10(@as(f32, 0x2.4p+4)));
    try std.testing.expectEqual(0x1.54484ap-120, exp10(@as(f32, -0x2.4p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.31p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.31p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.344p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343794p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.86ap+16)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0xf.424p+16)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0xf.424p+16)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x5.9f98p+0, exp10(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.348e46p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.348e44p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33aa02p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33aa04p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33ad16p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33ad18p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33afcap+8)));
    try std.testing.expectEqual(0x1.009388p+0, exp10(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0xf.f6ccdp-4, exp10(@as(f32, -0x4p-12)));
    try std.testing.expectEqual(0x1.000024p+0, exp10(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffdbp-4, exp10(@as(f32, -0x1p-20)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x4p-32)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x1p-40)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0xf.fffb3p+124, exp10(@as(f32, 0x2.688268p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x2.68826cp+4)));
    try std.testing.expectEqual(0x3.fffdfp-128, exp10(@as(f32, -0x2.5ee064p+4)));
    try std.testing.expectEqual(0x4.00004p-128, exp10(@as(f32, -0x2.5ee06p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344136p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344134p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344136p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344134p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33a714p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33a716p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33a714p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.33a716p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344136p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344134p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344136p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344134p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.23b2b4p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.23b2b6p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.23b2b4p+8)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.23b2b6p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344136p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344136p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344134p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343794p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343794p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343c64p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343c66p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343c64p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343c66p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344136p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344136p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0x1.344134p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343794p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f32, -0x1.343794p+12)));
    try std.testing.expectEqual(0x2.86b32cp-36, exp10(@as(f32, -0xa.6f431p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x7.764c68p+0, exp10(@as(f32, 0xd.f73d6p-4)));
    try std.testing.expectEqual(0x3.edf194p+4, exp10(@as(f32, 0x1.cc6776p+0)));
    try std.testing.expectEqual(0x7.6f01bp+16, exp10(@as(f32, 0x5.b00bdp+0)));
    try std.testing.expectEqual(0x7.6f012p+16, exp10(@as(f32, 0x5.b00bc8p+0)));
    try std.testing.expectEqual(std.math.inf(f32), exp10(@as(f32, 0xe.8b349p+4)));
    try std.testing.expectEqual(0x7.8e7e4p+8, exp10(@as(f32, 0x3.495c78p+0)));
    try std.testing.expectEqual(0x1.fad592p+52, exp10(@as(f32, 0xf.f33f6p+0)));

    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x3.e8p+8, exp10(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(0x1.999999999999ap-4, exp10(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0xc.097ce7bc90718p+116, exp10(@as(f64, 0x2.4p+4)));
    try std.testing.expectEqual(0x1.54484932d2e72p-120, exp10(@as(f64, -0x2.4p+4)));
    try std.testing.expectEqual(0x2.474a2dd05b374p+1012, exp10(@as(f64, 0x1.31p+8)));
    try std.testing.expectEqual(0x7.05b171494d5d4p-1016, exp10(@as(f64, -0x1.31p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.344p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343794p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.86ap+16)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0xf.424p+16)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0xf.424p+16)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x5.9f9802c8d1898p+0, exp10(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.348e46p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.348e44p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.348e45573a1dep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.348e45573a1ddp+8)));
    try std.testing.expectEqual(0x3.e5625e7a4219cp-1024, exp10(@as(f64, -0x1.33aa02p+8)));
    try std.testing.expectEqual(0x3.e5506d83c44fp-1024, exp10(@as(f64, -0x1.33aa04p+8)));
    try std.testing.expectEqual(0x3.e55965f4af484p-1024, exp10(@as(f64, -0x1.33aa03p+8)));
    try std.testing.expectEqual(0x3.ca263994bd44p-1024, exp10(@as(f64, -0x1.33ad16p+8)));
    try std.testing.expectEqual(0x3.ca14c60907b7p-1024, exp10(@as(f64, -0x1.33ad18p+8)));
    try std.testing.expectEqual(0x3.ca1d7fc4d6c3cp-1024, exp10(@as(f64, -0x1.33ad17p+8)));
    try std.testing.expectEqual(0x3.b2d8a908d0634p-1024, exp10(@as(f64, -0x1.33afcap+8)));
    try std.testing.expectEqual(0x1.009388004be7ep+0, exp10(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0xf.f6cccd4498ccp-4, exp10(@as(f64, -0x4p-12)));
    try std.testing.expectEqual(0x1.000024d7661e1p+0, exp10(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffdb289f2f38p-4, exp10(@as(f64, -0x1p-20)));
    try std.testing.expectEqual(0x1.0000000935d8ep+0, exp10(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffff6ca272p-4, exp10(@as(f64, -0x4p-32)));
    try std.testing.expectEqual(0x1.00000000024d7p+0, exp10(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0xf.ffffffffdb288p-4, exp10(@as(f64, -0x1p-40)));
    try std.testing.expectEqual(0x1.0000000000009p+0, exp10(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffff7p-4, exp10(@as(f64, -0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0xf.fffb372d9da8p+124, exp10(@as(f64, 0x2.688268p+4)));
    try std.testing.expectEqual(0x1.000046d066117p+128, exp10(@as(f64, 0x2.68826cp+4)));
    try std.testing.expectEqual(0x3.fffdf07e3d25p-128, exp10(@as(f64, -0x2.5ee064p+4)));
    try std.testing.expectEqual(0x4.00003df3ee9c4p-128, exp10(@as(f64, -0x2.5ee06p+4)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc598p+1020, exp10(@as(f64, 0x1.344134p+8)));
    try std.testing.expectEqual(0xf.fffffffffdd08p+1020, exp10(@as(f64, 0x1.34413509f79fep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc598p+1020, exp10(@as(f64, 0x1.344134p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0x4.0004027aecdd8p-1024, exp10(@as(f64, -0x1.33a714p+8)));
    try std.testing.expectEqual(0x3.fff196e1243d8p-1024, exp10(@as(f64, -0x1.33a716p+8)));
    try std.testing.expectEqual(0x3.fffffffffff8cp-1024, exp10(@as(f64, -0x1.33a7146f72a42p+8)));
    try std.testing.expectEqual(0x4.0004027aecdd8p-1024, exp10(@as(f64, -0x1.33a714p+8)));
    try std.testing.expectEqual(0x3.fff196e1243d8p-1024, exp10(@as(f64, -0x1.33a716p+8)));
    try std.testing.expectEqual(0x4.00000000008c4p-1024, exp10(@as(f64, -0x1.33a7146f72a41p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc598p+1020, exp10(@as(f64, 0x1.344134p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0xf.fffffffffdd08p+1020, exp10(@as(f64, 0x1.34413509f79fep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc598p+1020, exp10(@as(f64, 0x1.344134p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0xf.fffffffffdd08p+1020, exp10(@as(f64, 0x1.34413509f79fep+8)));
    try std.testing.expectEqual(0x8.00081bb1a65e8p-972, exp10(@as(f64, -0x1.23b2b4p+8)));
    try std.testing.expectEqual(0x7.ffe3447dac6e4p-972, exp10(@as(f64, -0x1.23b2b6p+8)));
    try std.testing.expectEqual(0x8.000000000095p-972, exp10(@as(f64, -0x1.23b2b470ae931p+8)));
    try std.testing.expectEqual(0x7.ffffffffff6e8p-972, exp10(@as(f64, -0x1.23b2b470ae932p+8)));
    try std.testing.expectEqual(0x8.00081bb1a65e8p-972, exp10(@as(f64, -0x1.23b2b4p+8)));
    try std.testing.expectEqual(0x7.ffe3447dac6e4p-972, exp10(@as(f64, -0x1.23b2b6p+8)));
    try std.testing.expectEqual(0x8.000000000095p-972, exp10(@as(f64, -0x1.23b2b470ae931p+8)));
    try std.testing.expectEqual(0x7.ffffffffff6e8p-972, exp10(@as(f64, -0x1.23b2b470ae932p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344136p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344136p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343794p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343794p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343c64p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343c66p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343c640523781p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343c640523782p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343c64p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343c66p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343c640523781p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343c640523782p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344136p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344136p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp10(@as(f64, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343794p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343792p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343794p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f64, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x2.86b32a000000ep-36, exp10(@as(f64, -0xa.6f431p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x7.764c69914e79cp+0, exp10(@as(f64, 0xd.f73d6p-4)));
    try std.testing.expectEqual(0x3.edf195be93514p+4, exp10(@as(f64, 0x1.cc6776p+0)));
    try std.testing.expectEqual(0x7.6f01ac1f6639cp+16, exp10(@as(f64, 0x5.b00bdp+0)));
    try std.testing.expectEqual(0x7.6f012330be264p+16, exp10(@as(f64, 0x5.b00bc8p+0)));
    try std.testing.expectEqual(0x7.6f0181f100c4cp+16, exp10(@as(f64, 0x5.b00bcd891ffe8p+0)));
    try std.testing.expectEqual(0x7.6f0181f100c08p+16, exp10(@as(f64, 0x5.b00bcd891ffe4p+0)));
    try std.testing.expectEqual(0x2.04e945593f42p+772, exp10(@as(f64, 0xe.8b349p+4)));
    try std.testing.expectEqual(0x7.8e7e436efa1d4p+8, exp10(@as(f64, 0x3.495c78p+0)));
    try std.testing.expectEqual(0x1.fad59245e4f68p+52, exp10(@as(f64, 0xf.f33f6p+0)));

    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x3.e8p+8, exp10(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(0x1.999999999999999ap-4, exp10(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0xc.097ce7bc90715b3p+116, exp10(@as(f80, 0x2.4p+4)));
    try std.testing.expectEqual(0x1.54484932d2e725a6p-120, exp10(@as(f80, -0x2.4p+4)));
    try std.testing.expectEqual(0x2.474a2dd05b3749f8p+1012, exp10(@as(f80, 0x1.31p+8)));
    try std.testing.expectEqual(0x7.05b171494d5d41ep-1016, exp10(@as(f80, -0x1.31p+8)));
    try std.testing.expectEqual(0xd.72cb2a95c7ef6cdp+16380, exp10(@as(f80, 0x1.344p+12)));
    try std.testing.expectEqual(0x1.30923e47949abf8p-16384, exp10(@as(f80, -0x1.344p+12)));
    try std.testing.expectEqual(0x4.009395d78ebc9b68p-16384, exp10(@as(f80, -0x1.343792p+12)));
    try std.testing.expectEqual(0x3.ff6cdaadaae05f2p-16384, exp10(@as(f80, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f5p-16384, exp10(@as(f80, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x3.fffffffff80d767p-16384, exp10(@as(f80, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffcbfp-16384, exp10(@as(f80, -0x1.343793004f503232p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.86ap+16)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f80, -0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0xf.424p+16)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f80, -0xf.424p+16)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x5.9f9802c8d1896578p+0, exp10(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x2.0003093cc02bf7cp+1024, exp10(@as(f80, 0x1.348e46p+8)));
    try std.testing.expectEqual(0x1.fff9d36b1c2656fp+1024, exp10(@as(f80, 0x1.348e44p+8)));
    try std.testing.expectEqual(0x2.000000000028a374p+1024, exp10(@as(f80, 0x1.348e45573a1dep+8)));
    try std.testing.expectEqual(0x1.ffffffffffdef4acp+1024, exp10(@as(f80, 0x1.348e45573a1ddp+8)));
    try std.testing.expectEqual(0x1.fffffffffffffbc4p+1024, exp10(@as(f80, 0x1.348e45573a1dd72cp+8)));
    try std.testing.expectEqual(0x3.e5625e7a4219b1f4p-1024, exp10(@as(f80, -0x1.33aa02p+8)));
    try std.testing.expectEqual(0x3.e5506d83c44ee174p-1024, exp10(@as(f80, -0x1.33aa04p+8)));
    try std.testing.expectEqual(0x3.e55965f4af4844bp-1024, exp10(@as(f80, -0x1.33aa03p+8)));
    try std.testing.expectEqual(0x3.ca263994bd441e7cp-1024, exp10(@as(f80, -0x1.33ad16p+8)));
    try std.testing.expectEqual(0x3.ca14c60907b717bp-1024, exp10(@as(f80, -0x1.33ad18p+8)));
    try std.testing.expectEqual(0x3.ca1d7fc4d6c3bc58p-1024, exp10(@as(f80, -0x1.33ad17p+8)));
    try std.testing.expectEqual(0x3.b2d8a908d0634328p-1024, exp10(@as(f80, -0x1.33afcap+8)));
    try std.testing.expectEqual(0x1.009388004be7e55ap+0, exp10(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0xf.f6cccd4498cbd18p-4, exp10(@as(f80, -0x4p-12)));
    try std.testing.expectEqual(0x1.000024d7661e0f64p+0, exp10(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffdb289f2f39cep-4, exp10(@as(f80, -0x1p-20)));
    try std.testing.expectEqual(0x1.0000000935d8de06p+0, exp10(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffff6ca27225p-4, exp10(@as(f80, -0x4p-32)));
    try std.testing.expectEqual(0x1.00000000024d7638p+0, exp10(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0xf.ffffffffdb289c9p-4, exp10(@as(f80, -0x1p-40)));
    try std.testing.expectEqual(0x1.000000000000935ep+0, exp10(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffff6ca2p-4, exp10(@as(f80, -0x4p-52)));
    try std.testing.expectEqual(0x1.0000000000000024p+0, exp10(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0xf.fffffffffffffdbp-4, exp10(@as(f80, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x1p-10000)));
    try std.testing.expectEqual(0xf.fffb372d9da7f63p+124, exp10(@as(f80, 0x2.688268p+4)));
    try std.testing.expectEqual(0x1.000046d066116874p+128, exp10(@as(f80, 0x2.68826cp+4)));
    try std.testing.expectEqual(0x3.fffdf07e3d250a68p-128, exp10(@as(f80, -0x2.5ee064p+4)));
    try std.testing.expectEqual(0x4.00003df3ee9c5b1p-128, exp10(@as(f80, -0x2.5ee06p+4)));
    try std.testing.expectEqual(0x1.0002368555061c26p+1024, exp10(@as(f80, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc5990cp+1020, exp10(@as(f80, 0x1.344134p+8)));
    try std.testing.expectEqual(0xf.fffffffffdd04f8p+1020, exp10(@as(f80, 0x1.34413509f79fep+8)));
    try std.testing.expectEqual(0x1.0002368555061c26p+1024, exp10(@as(f80, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc5990cp+1020, exp10(@as(f80, 0x1.344134p+8)));
    try std.testing.expectEqual(0x1.000000000001dc5cp+1024, exp10(@as(f80, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0x4.0004027aecdd7a7p-1024, exp10(@as(f80, -0x1.33a714p+8)));
    try std.testing.expectEqual(0x3.fff196e1243d704cp-1024, exp10(@as(f80, -0x1.33a716p+8)));
    try std.testing.expectEqual(0x3.fffffffffff8dbfcp-1024, exp10(@as(f80, -0x1.33a7146f72a42p+8)));
    try std.testing.expectEqual(0x4.0004027aecdd7a7p-1024, exp10(@as(f80, -0x1.33a714p+8)));
    try std.testing.expectEqual(0x3.fff196e1243d704cp-1024, exp10(@as(f80, -0x1.33a716p+8)));
    try std.testing.expectEqual(0x4.00000000008c3988p-1024, exp10(@as(f80, -0x1.33a7146f72a41p+8)));
    try std.testing.expectEqual(0x1.0002368555061c26p+1024, exp10(@as(f80, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc5990cp+1020, exp10(@as(f80, 0x1.344134p+8)));
    try std.testing.expectEqual(0x1.000000000001dc5cp+1024, exp10(@as(f80, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0xf.fffffffffdd04f8p+1020, exp10(@as(f80, 0x1.34413509f79fep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffd6bp+1020, exp10(@as(f80, 0x1.34413509f79fef3p+8)));
    try std.testing.expectEqual(0xf.ffffffffffff8dp+1020, exp10(@as(f80, 0x1.34413509f79fef2ep+8)));
    try std.testing.expectEqual(0x1.0002368555061c26p+1024, exp10(@as(f80, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc5990cp+1020, exp10(@as(f80, 0x1.344134p+8)));
    try std.testing.expectEqual(0x1.000000000001dc5cp+1024, exp10(@as(f80, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0xf.fffffffffdd04f8p+1020, exp10(@as(f80, 0x1.34413509f79fep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffd6bp+1020, exp10(@as(f80, 0x1.34413509f79fef3p+8)));
    try std.testing.expectEqual(0xf.ffffffffffff8dp+1020, exp10(@as(f80, 0x1.34413509f79fef2ep+8)));
    try std.testing.expectEqual(0x8.00081bb1a65e9f4p-972, exp10(@as(f80, -0x1.23b2b4p+8)));
    try std.testing.expectEqual(0x7.ffe3447dac6e5168p-972, exp10(@as(f80, -0x1.23b2b6p+8)));
    try std.testing.expectEqual(0x8.0000000000951bfp-972, exp10(@as(f80, -0x1.23b2b470ae931p+8)));
    try std.testing.expectEqual(0x7.ffffffffff6e60d8p-972, exp10(@as(f80, -0x1.23b2b470ae932p+8)));
    try std.testing.expectEqual(0x8.000000000000045p-972, exp10(@as(f80, -0x1.23b2b470ae931818p+8)));
    try std.testing.expectEqual(0x7.ffffffffffffdf7p-972, exp10(@as(f80, -0x1.23b2b470ae93181ap+8)));
    try std.testing.expectEqual(0x8.00081bb1a65e9f4p-972, exp10(@as(f80, -0x1.23b2b4p+8)));
    try std.testing.expectEqual(0x7.ffe3447dac6e5168p-972, exp10(@as(f80, -0x1.23b2b6p+8)));
    try std.testing.expectEqual(0x8.0000000000951bfp-972, exp10(@as(f80, -0x1.23b2b470ae931p+8)));
    try std.testing.expectEqual(0x7.ffffffffff6e60d8p-972, exp10(@as(f80, -0x1.23b2b470ae932p+8)));
    try std.testing.expectEqual(0x8.000000000000045p-972, exp10(@as(f80, -0x1.23b2b470ae931818p+8)));
    try std.testing.expectEqual(0x7.ffffffffffffdf7p-972, exp10(@as(f80, -0x1.23b2b470ae93181ap+8)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.344136p+12)));
    try std.testing.expectEqual(0xf.fd9bc4394211e89p+16380, exp10(@as(f80, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(0xf.ffffffffdd04f79p+16380, exp10(@as(f80, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(0xf.fffffffffffd6bp+16380, exp10(@as(f80, 0x1.34413509f79fef3p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.344136p+12)));
    try std.testing.expectEqual(0xf.fd9bc4394211e89p+16380, exp10(@as(f80, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(0xf.ffffffffdd04f79p+16380, exp10(@as(f80, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.34413509f79fef32p+12)));
    try std.testing.expectEqual(0x4.009395d78ebc9b68p-16384, exp10(@as(f80, -0x1.343792p+12)));
    try std.testing.expectEqual(0x3.ff6cdaadaae05f2p-16384, exp10(@as(f80, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f5p-16384, exp10(@as(f80, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x3.fffffffff80d767p-16384, exp10(@as(f80, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffcbfp-16384, exp10(@as(f80, -0x1.343793004f503232p+12)));
    try std.testing.expectEqual(0x4.009395d78ebc9b68p-16384, exp10(@as(f80, -0x1.343792p+12)));
    try std.testing.expectEqual(0x3.ff6cdaadaae05f2p-16384, exp10(@as(f80, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f5p-16384, exp10(@as(f80, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x3.fffffffff80d767p-16384, exp10(@as(f80, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x4.000000000000f2a8p-16384, exp10(@as(f80, -0x1.343793004f50323p+12)));
    try std.testing.expectEqual(0x2.00017a9fe296ea9p-16384, exp10(@as(f80, -0x1.343c64p+12)));
    try std.testing.expectEqual(0x1.ff6e31d8368f1dep-16384, exp10(@as(f80, -0x1.343c66p+12)));
    try std.testing.expectEqual(0x2.0000000000330e2p-16384, exp10(@as(f80, -0x1.343c640523781p+12)));
    try std.testing.expectEqual(0x1.fffffffffb9821b8p-16384, exp10(@as(f80, -0x1.343c640523782p+12)));
    try std.testing.expectEqual(0x1.ffffffffffffd2ap-16384, exp10(@as(f80, -0x1.343c6405237810b2p+12)));
    try std.testing.expectEqual(0x2.00017a9fe296ea9p-16384, exp10(@as(f80, -0x1.343c64p+12)));
    try std.testing.expectEqual(0x1.ff6e31d8368f1dep-16384, exp10(@as(f80, -0x1.343c66p+12)));
    try std.testing.expectEqual(0x2.0000000000330e2p-16384, exp10(@as(f80, -0x1.343c640523781p+12)));
    try std.testing.expectEqual(0x1.fffffffffb9821b8p-16384, exp10(@as(f80, -0x1.343c640523782p+12)));
    try std.testing.expectEqual(0x2.00000000000065f8p-16384, exp10(@as(f80, -0x1.343c6405237810bp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.344136p+12)));
    try std.testing.expectEqual(0xf.fd9bc4394211e89p+16380, exp10(@as(f80, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(0xf.ffffffffdd04f79p+16380, exp10(@as(f80, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.34413509f79fef32p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd6bp+16380, exp10(@as(f80, 0x1.34413509f79fef3p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.344136p+12)));
    try std.testing.expectEqual(0xf.fd9bc4394211e89p+16380, exp10(@as(f80, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(0xf.ffffffffdd04f79p+16380, exp10(@as(f80, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp10(@as(f80, 0x1.34413509f79fef32p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd6bp+16380, exp10(@as(f80, 0x1.34413509f79fef3p+12)));
    try std.testing.expectEqual(0x4.009395d78ebc9b68p-16384, exp10(@as(f80, -0x1.343792p+12)));
    try std.testing.expectEqual(0x3.ff6cdaadaae05f2p-16384, exp10(@as(f80, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f5p-16384, exp10(@as(f80, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x3.fffffffff80d767p-16384, exp10(@as(f80, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x4.000000000000f2a8p-16384, exp10(@as(f80, -0x1.343793004f50323p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffcbfp-16384, exp10(@as(f80, -0x1.343793004f503232p+12)));
    try std.testing.expectEqual(0x4.009395d78ebc9b68p-16384, exp10(@as(f80, -0x1.343792p+12)));
    try std.testing.expectEqual(0x3.ff6cdaadaae05f2p-16384, exp10(@as(f80, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f5p-16384, exp10(@as(f80, -0x1.343793004f503p+12)));
    try std.testing.expectEqual(0x3.fffffffff80d767p-16384, exp10(@as(f80, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x4.000000000000f2a8p-16384, exp10(@as(f80, -0x1.343793004f50323p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffcbfp-16384, exp10(@as(f80, -0x1.343793004f503232p+12)));
    try std.testing.expectEqual(0x2.86b32a000000da34p-36, exp10(@as(f80, -0xa.6f431p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x7.764c69914e79a9d8p+0, exp10(@as(f80, 0xd.f73d6p-4)));
    try std.testing.expectEqual(0x3.edf195be935146fp+4, exp10(@as(f80, 0x1.cc6776p+0)));
    try std.testing.expectEqual(0x7.6f01ac1f6639ae28p+16, exp10(@as(f80, 0x5.b00bdp+0)));
    try std.testing.expectEqual(0x7.6f012330be263708p+16, exp10(@as(f80, 0x5.b00bc8p+0)));
    try std.testing.expectEqual(0x7.6f0181f100c4cea8p+16, exp10(@as(f80, 0x5.b00bcd891ffe8p+0)));
    try std.testing.expectEqual(0x7.6f0181f100c0873p+16, exp10(@as(f80, 0x5.b00bcd891ffe4p+0)));
    try std.testing.expectEqual(0x7.6f0181f100c20fdp+16, exp10(@as(f80, 0x5.b00bcd891ffe56fp+0)));
    try std.testing.expectEqual(0x2.04e945593f41f0c8p+772, exp10(@as(f80, 0xe.8b349p+4)));
    try std.testing.expectEqual(0x7.8e7e436efa1d5b18p+8, exp10(@as(f80, 0x3.495c78p+0)));
    try std.testing.expectEqual(0x1.fad59245e4f68156p+52, exp10(@as(f80, 0xf.f33f6p+0)));

    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x3.e8p+8, exp10(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(0x1.999999999999999999999999999ap-4, exp10(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0xc.097ce7bc90715b34b9f1p+116, exp10(@as(f128, 0x2.4p+4)));
    // try std.testing.expectEqual(0x1.54484932d2e725a5bbca17a3aba1p-120, exp10(@as(f128, -0x2.4p+4)));
    try std.testing.expectEqual(0x2.474a2dd05b3749f93370cc755feap+1012, exp10(@as(f128, 0x1.31p+8)));
    try std.testing.expectEqual(0x7.05b171494d5d41e198d66d5ff4a8p-1016, exp10(@as(f128, -0x1.31p+8)));
    // try std.testing.expectEqual(0xd.72cb2a95c7ef6cce81bf1e825ba8p+16380, exp10(@as(f128, 0x1.344p+12)));
    try std.testing.expectEqual(0x1.30923e47949abf816b7d38ebc01p-16384, exp10(@as(f128, -0x1.344p+12)));
    // try std.testing.expectEqual(0x4.009395d78ebc9b64a0aa93fc93dp-16384, exp10(@as(f128, -0x1.343792p+12)));
    // try std.testing.expectEqual(0x3.ff6cdaadaae05f1d9410e8bb22f8p-16384, exp10(@as(f128, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f4e77c5e6c4ab5p-16384, exp10(@as(f128, -0x1.343793004f503p+12)));
    // try std.testing.expectEqual(0x3.fffffffff80d76709d230e22dc24p-16384, exp10(@as(f128, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffcbed07250d70bb4p-16384, exp10(@as(f128, -0x1.343793004f503232p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.86ap+16)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f128, -0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0xf.424p+16)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f128, -0xf.424p+16)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, exp10(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x5.9f9802c8d189657416ee3fd818c4p+0, exp10(@as(f128, 0xcp-4)));
    // try std.testing.expectEqual(0x2.0003093cc02bf7be0dd170fd426p+1024, exp10(@as(f128, 0x1.348e46p+8)));
    try std.testing.expectEqual(0x1.fff9d36b1c2656ef7dd26d07ce3fp+1024, exp10(@as(f128, 0x1.348e44p+8)));
    try std.testing.expectEqual(0x2.000000000028a3736b9d8e05898ep+1024, exp10(@as(f128, 0x1.348e45573a1dep+8)));
    try std.testing.expectEqual(0x1.ffffffffffdef4ac7cc8392399ffp+1024, exp10(@as(f128, 0x1.348e45573a1ddp+8)));
    // try std.testing.expectEqual(0x1.fffffffffffffbc4285657a030a5p+1024, exp10(@as(f128, 0x1.348e45573a1dd72cp+8)));
    // try std.testing.expectEqual(0x3.e5625e7a4219b1f23b7f41e1933ep-1024, exp10(@as(f128, -0x1.33aa02p+8)));
    try std.testing.expectEqual(0x3.e5506d83c44ee174a1cd22369ecep-1024, exp10(@as(f128, -0x1.33aa04p+8)));
    try std.testing.expectEqual(0x3.e55965f4af4844b0187da80e25ep-1024, exp10(@as(f128, -0x1.33aa03p+8)));
    try std.testing.expectEqual(0x3.ca263994bd441e7c46ea7c3f2964p-1024, exp10(@as(f128, -0x1.33ad16p+8)));
    try std.testing.expectEqual(0x3.ca14c60907b717ae36dc1f6cac46p-1024, exp10(@as(f128, -0x1.33ad18p+8)));
    try std.testing.expectEqual(0x3.ca1d7fc4d6c3bc586b2b65fe7e66p-1024, exp10(@as(f128, -0x1.33ad17p+8)));
    // try std.testing.expectEqual(0x3.b2d8a908d063432616cd82f6a4eep-1024, exp10(@as(f128, -0x1.33afcap+8)));
    try std.testing.expectEqual(0x1.009388004be7e5592e3f8d6c8273p+0, exp10(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0xf.f6cccd4498cbd185346978c830b8p-4, exp10(@as(f128, -0x4p-12)));
    // try std.testing.expectEqual(0x1.000024d7661e0f63a0af573a6217p+0, exp10(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffdb289f2f39ce2e8d9a96332d8p-4, exp10(@as(f128, -0x1p-20)));
    try std.testing.expectEqual(0x1.0000000935d8de0514d4506ab26bp+0, exp10(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffff6ca27224fbfbecc88f737p-4, exp10(@as(f128, -0x4p-32)));
    try std.testing.expectEqual(0x1.00000000024d763776ad4954f491p+0, exp10(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0xf.ffffffffdb289c8895803f43d3ep-4, exp10(@as(f128, -0x1p-40)));
    try std.testing.expectEqual(0x1.000000000000935d8dddaaa8d681p+0, exp10(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffff6ca2722255577e538p-4, exp10(@as(f128, -0x4p-52)));
    try std.testing.expectEqual(0x1.0000000000000024d763776aaa2bp+0, exp10(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0xf.fffffffffffffdb289c889555d5p-4, exp10(@as(f128, -0x1p-60)));
    try std.testing.expectEqual(0x1.00000000000000000000000024d7p+0, exp10(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffdb288p-4, exp10(@as(f128, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x1p-10000)));
    try std.testing.expectEqual(0xf.fffb372d9da7f632ed2ce00c161p+124, exp10(@as(f128, 0x2.688268p+4)));
    try std.testing.expectEqual(0x1.000046d0661168747a32f41ab7b4p+128, exp10(@as(f128, 0x2.68826cp+4)));
    try std.testing.expectEqual(0x3.fffdf07e3d250a69f129cbd543a6p-128, exp10(@as(f128, -0x2.5ee064p+4)));
    // try std.testing.expectEqual(0x4.00003df3ee9c5b137b1279114d6p-128, exp10(@as(f128, -0x2.5ee06p+4)));
    // try std.testing.expectEqual(0x1.0002368555061c26d904a4f4deb5p+1024, exp10(@as(f128, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc5990b8fb9ade558848p+1020, exp10(@as(f128, 0x1.344134p+8)));
    try std.testing.expectEqual(0xf.fffffffffdd04f7930e506f06848p+1020, exp10(@as(f128, 0x1.34413509f79fep+8)));
    // try std.testing.expectEqual(0x1.0002368555061c26d904a4f4deb5p+1024, exp10(@as(f128, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc5990b8fb9ade558848p+1020, exp10(@as(f128, 0x1.344134p+8)));
    try std.testing.expectEqual(0x1.000000000001dc5b0a78f837f53dp+1024, exp10(@as(f128, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0x4.0004027aecdd7a6d329b9571e458p-1024, exp10(@as(f128, -0x1.33a714p+8)));
    // try std.testing.expectEqual(0x3.fff196e1243d704d262ae46d3424p-1024, exp10(@as(f128, -0x1.33a716p+8)));
    try std.testing.expectEqual(0x3.fffffffffff8dbfb531fe67239fap-1024, exp10(@as(f128, -0x1.33a7146f72a42p+8)));
    try std.testing.expectEqual(0x4.0004027aecdd7a6d329b9571e458p-1024, exp10(@as(f128, -0x1.33a714p+8)));
    // try std.testing.expectEqual(0x3.fff196e1243d704d262ae46d3424p-1024, exp10(@as(f128, -0x1.33a716p+8)));
    try std.testing.expectEqual(0x4.00000000008c398930ca98b1d098p-1024, exp10(@as(f128, -0x1.33a7146f72a41p+8)));
    // try std.testing.expectEqual(0x1.0002368555061c26d904a4f4deb5p+1024, exp10(@as(f128, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc5990b8fb9ade558848p+1020, exp10(@as(f128, 0x1.344134p+8)));
    try std.testing.expectEqual(0x1.000000000001dc5b0a78f837f53dp+1024, exp10(@as(f128, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0xf.fffffffffdd04f7930e506f06848p+1020, exp10(@as(f128, 0x1.34413509f79fep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffd6afd688d920ac5p+1020, exp10(@as(f128, 0x1.34413509f79fef3p+8)));
    // try std.testing.expectEqual(0xf.ffffffffffff8d010f9a03cc57a8p+1020, exp10(@as(f128, 0x1.34413509f79fef2ep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffbffffffffffa35a8p+1020, exp10(@as(f128, 0x1.34413509f79fef2f625b0205a88p+8)));
    // try std.testing.expectEqual(0x1.0002368555061c26d904a4f4deb5p+1024, exp10(@as(f128, 0x1.344136p+8)));
    try std.testing.expectEqual(0xf.ffd9b994fc5990b8fb9ade558848p+1020, exp10(@as(f128, 0x1.344134p+8)));
    try std.testing.expectEqual(0x1.000000000001dc5b0a78f837f53dp+1024, exp10(@as(f128, 0x1.34413509f79ffp+8)));
    try std.testing.expectEqual(0xf.fffffffffdd04f7930e506f06848p+1020, exp10(@as(f128, 0x1.34413509f79fep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffd6afd688d920ac5p+1020, exp10(@as(f128, 0x1.34413509f79fef3p+8)));
    // try std.testing.expectEqual(0xf.ffffffffffff8d010f9a03cc57a8p+1020, exp10(@as(f128, 0x1.34413509f79fef2ep+8)));
    // try std.testing.expectEqual(0xf.ffffffffffffc0000000000ca16p+1020, exp10(@as(f128, 0x1.34413509f79fef2f625b0205a9p+8)));
    try std.testing.expectEqual(0x8.00081bb1a65e9f3993115f203f68p-972, exp10(@as(f128, -0x1.23b2b4p+8)));
    try std.testing.expectEqual(0x7.ffe3447dac6e5164756078e4d228p-972, exp10(@as(f128, -0x1.23b2b6p+8)));
    try std.testing.expectEqual(0x8.0000000000951bf165b2bc2da89p-972, exp10(@as(f128, -0x1.23b2b470ae931p+8)));
    // try std.testing.expectEqual(0x7.ffffffffff6e60d5aa5d6a953d2cp-972, exp10(@as(f128, -0x1.23b2b470ae932p+8)));
    try std.testing.expectEqual(0x8.000000000000044ade6f0e19f4cp-972, exp10(@as(f128, -0x1.23b2b470ae931818p+8)));
    // try std.testing.expectEqual(0x7.ffffffffffffdf737af7a36fcap-972, exp10(@as(f128, -0x1.23b2b470ae93181ap+8)));
    // try std.testing.expectEqual(0x7.fffffffffffffffffffffffb3284p-972, exp10(@as(f128, -0x1.23b2b470ae9318183ba772361cp+8)));
    try std.testing.expectEqual(0x8.00081bb1a65e9f3993115f203f68p-972, exp10(@as(f128, -0x1.23b2b4p+8)));
    try std.testing.expectEqual(0x7.ffe3447dac6e5164756078e4d228p-972, exp10(@as(f128, -0x1.23b2b6p+8)));
    try std.testing.expectEqual(0x8.0000000000951bf165b2bc2da89p-972, exp10(@as(f128, -0x1.23b2b470ae931p+8)));
    // try std.testing.expectEqual(0x7.ffffffffff6e60d5aa5d6a953d2cp-972, exp10(@as(f128, -0x1.23b2b470ae932p+8)));
    try std.testing.expectEqual(0x8.000000000000044ade6f0e19f4cp-972, exp10(@as(f128, -0x1.23b2b470ae931818p+8)));
    // try std.testing.expectEqual(0x7.ffffffffffffdf737af7a36fcap-972, exp10(@as(f128, -0x1.23b2b470ae93181ap+8)));
    // try std.testing.expectEqual(0x8.000000000000000000000004686p-972, exp10(@as(f128, -0x1.23b2b470ae9318183ba772361b8p+8)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.344136p+12)));
    // try std.testing.expectEqual(0xf.fd9bc4394211e88a367d0df4a918p+16380, exp10(@as(f128, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(0xf.ffffffffdd04f7930e744857852p+16380, exp10(@as(f128, 0x1.34413509f79fep+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd6afd688d920af6ep+16380, exp10(@as(f128, 0x1.34413509f79fef3p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.344136p+12)));
    // try std.testing.expectEqual(0xf.fd9bc4394211e88a367d0df4a918p+16380, exp10(@as(f128, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(0xf.ffffffffdd04f7930e744857852p+16380, exp10(@as(f128, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79fef32p+12)));
    // try std.testing.expectEqual(0x4.009395d78ebc9b64a0aa93fc93dp-16384, exp10(@as(f128, -0x1.343792p+12)));
    // try std.testing.expectEqual(0x3.ff6cdaadaae05f1d9410e8bb22f8p-16384, exp10(@as(f128, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f4e77c5e6c4ab5p-16384, exp10(@as(f128, -0x1.343793004f503p+12)));
    // try std.testing.expectEqual(0x3.fffffffff80d76709d230e22dc24p-16384, exp10(@as(f128, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffcbed07250d70bb4p-16384, exp10(@as(f128, -0x1.343793004f503232p+12)));
    // try std.testing.expectEqual(0x4.009395d78ebc9b64a0aa93fc93dp-16384, exp10(@as(f128, -0x1.343792p+12)));
    // try std.testing.expectEqual(0x3.ff6cdaadaae05f1d9410e8bb22f8p-16384, exp10(@as(f128, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f4e77c5e6c4ab5p-16384, exp10(@as(f128, -0x1.343793004f503p+12)));
    // try std.testing.expectEqual(0x3.fffffffff80d76709d230e22dc24p-16384, exp10(@as(f128, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x4.000000000000f2a822e062c22edcp-16384, exp10(@as(f128, -0x1.343793004f50323p+12)));
    try std.testing.expectEqual(0x2.00017a9fe296ea91c2392281858p-16384, exp10(@as(f128, -0x1.343c64p+12)));
    try std.testing.expectEqual(0x1.ff6e31d8368f1ddd2bedf091d16p-16384, exp10(@as(f128, -0x1.343c66p+12)));
    try std.testing.expectEqual(0x2.0000000000330e22f661ec16a8fcp-16384, exp10(@as(f128, -0x1.343c640523781p+12)));
    try std.testing.expectEqual(0x1.fffffffffb9821b409117e700974p-16384, exp10(@as(f128, -0x1.343c640523782p+12)));
    // try std.testing.expectEqual(0x1.ffffffffffffd29ca45194e72e5cp-16384, exp10(@as(f128, -0x1.343c6405237810b2p+12)));
    try std.testing.expectEqual(0x2.00017a9fe296ea91c2392281858p-16384, exp10(@as(f128, -0x1.343c64p+12)));
    try std.testing.expectEqual(0x1.ff6e31d8368f1ddd2bedf091d16p-16384, exp10(@as(f128, -0x1.343c66p+12)));
    try std.testing.expectEqual(0x2.0000000000330e22f661ec16a8fcp-16384, exp10(@as(f128, -0x1.343c640523781p+12)));
    try std.testing.expectEqual(0x1.fffffffffb9821b409117e700974p-16384, exp10(@as(f128, -0x1.343c640523782p+12)));
    // try std.testing.expectEqual(0x2.00000000000065fa322f3f8fe298p-16384, exp10(@as(f128, -0x1.343c6405237810bp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.344136p+12)));
    // try std.testing.expectEqual(0xf.fd9bc4394211e88a367d0df4a918p+16380, exp10(@as(f128, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(0xf.ffffffffdd04f7930e744857852p+16380, exp10(@as(f128, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79fef32p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd6afd688d920af6ep+16380, exp10(@as(f128, 0x1.34413509f79fef3p+12)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffaf9p+16380, exp10(@as(f128, 0x1.34413509f79fef311f12b35816f9p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79fef311f12b35817p+12)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffee912bp+16380, exp10(@as(f128, 0x1.34413509f79fef311f12b358168p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.344136p+12)));
    // try std.testing.expectEqual(0xf.fd9bc4394211e88a367d0df4a918p+16380, exp10(@as(f128, 0x1.344134p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79ffp+12)));
    try std.testing.expectEqual(0xf.ffffffffdd04f7930e744857852p+16380, exp10(@as(f128, 0x1.34413509f79fep+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79fef32p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd6afd688d920af6ep+16380, exp10(@as(f128, 0x1.34413509f79fef3p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79fef311f12b35816fap+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp10(@as(f128, 0x1.34413509f79fef311f12b35817p+12)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffee912bp+16380, exp10(@as(f128, 0x1.34413509f79fef311f12b358168p+12)));
    // try std.testing.expectEqual(0x4.009395d78ebc9b64a0aa93fc93dp-16384, exp10(@as(f128, -0x1.343792p+12)));
    // try std.testing.expectEqual(0x3.ff6cdaadaae05f1d9410e8bb22f8p-16384, exp10(@as(f128, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f4e77c5e6c4ab5p-16384, exp10(@as(f128, -0x1.343793004f503p+12)));
    // try std.testing.expectEqual(0x3.fffffffff80d76709d230e22dc24p-16384, exp10(@as(f128, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x4.000000000000f2a822e062c22edcp-16384, exp10(@as(f128, -0x1.343793004f50323p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffcbed07250d70bb4p-16384, exp10(@as(f128, -0x1.343793004f503232p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffffffffffffffaa5p-16384, exp10(@as(f128, -0x1.343793004f503231a589bac27c39p+12)));
    try std.testing.expectEqual(0x4.0000000000000000000000207a24p-16384, exp10(@as(f128, -0x1.343793004f503231a589bac27cp+12)));
    // try std.testing.expectEqual(0x3.ffffffffffffffffffffffd6cb6p-16384, exp10(@as(f128, -0x1.343793004f503231a589bac27c8p+12)));
    // try std.testing.expectEqual(0x4.009395d78ebc9b64a0aa93fc93dp-16384, exp10(@as(f128, -0x1.343792p+12)));
    // try std.testing.expectEqual(0x3.ff6cdaadaae05f1d9410e8bb22f8p-16384, exp10(@as(f128, -0x1.343794p+12)));
    try std.testing.expectEqual(0x4.0000000001434f4e77c5e6c4ab5p-16384, exp10(@as(f128, -0x1.343793004f503p+12)));
    // try std.testing.expectEqual(0x3.fffffffff80d76709d230e22dc24p-16384, exp10(@as(f128, -0x1.343793004f504p+12)));
    try std.testing.expectEqual(0x4.000000000000f2a822e062c22edcp-16384, exp10(@as(f128, -0x1.343793004f50323p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffcbed07250d70bb4p-16384, exp10(@as(f128, -0x1.343793004f503232p+12)));
    // try std.testing.expectEqual(0x4.0000000000000000000000003dbp-16384, exp10(@as(f128, -0x1.343793004f503231a589bac27c38p+12)));
    try std.testing.expectEqual(0x4.0000000000000000000000207a24p-16384, exp10(@as(f128, -0x1.343793004f503231a589bac27cp+12)));
    // try std.testing.expectEqual(0x3.ffffffffffffffffffffffd6cb6p-16384, exp10(@as(f128, -0x1.343793004f503231a589bac27c8p+12)));
    try std.testing.expectEqual(0x2.86b32a000000da34970abbb69d3ap-36, exp10(@as(f128, -0xa.6f431p+0)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, exp10(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x7.764c69914e79a9d63d384048a754p+0, exp10(@as(f128, 0xd.f73d6p-4)));
    try std.testing.expectEqual(0x3.edf195be935146f144ca1eb390e2p+4, exp10(@as(f128, 0x1.cc6776p+0)));
    try std.testing.expectEqual(0x7.6f01ac1f6639ae2a1bbe346dd88p+16, exp10(@as(f128, 0x5.b00bdp+0)));
    try std.testing.expectEqual(0x7.6f012330be26370bdca477f8f5fp+16, exp10(@as(f128, 0x5.b00bc8p+0)));
    // try std.testing.expectEqual(0x7.6f0181f100c4cea8583757f40d14p+16, exp10(@as(f128, 0x5.b00bcd891ffe8p+0)));
    // try std.testing.expectEqual(0x7.6f0181f100c08733087a227b3b38p+16, exp10(@as(f128, 0x5.b00bcd891ffe4p+0)));
    // try std.testing.expectEqual(0x7.6f0181f100c20fcf53ce3264ffecp+16, exp10(@as(f128, 0x5.b00bcd891ffe56fp+0)));
    try std.testing.expectEqual(0x2.04e945593f41f0c960f8e9467afcp+772, exp10(@as(f128, 0xe.8b349p+4)));
    // try std.testing.expectEqual(0x7.8e7e436efa1d5b19bda5590a91c4p+8, exp10(@as(f128, 0x3.495c78p+0)));
    try std.testing.expectEqual(0x1.fad59245e4f681552bf17b541a5fp+52, exp10(@as(f128, 0xf.f33f6p+0)));
}
