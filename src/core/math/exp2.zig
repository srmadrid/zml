const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const exp2_data = @import("exp2_data.zig");
const exp_data = @import("exp_data.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub fn exp2(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return exp2(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, exp2_32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_exp2f.c
                    return exp2_32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_exp2.c
                    return exp2_64(x);
                },
                f80 => return cast(f80, exp2_128(cast(f128, x, .{})), .{}),
                f128 => {
                    // openlibm/ld128/s_exp2l.c
                    return exp2_128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

inline fn top12_32(x: f32) u32 {
    return @as(u32, @bitCast(x)) >> 20;
}

fn exp2_32(x: f32) f32 {
    const xd: f64 = cast(f64, x, .{});
    const abstop: u32 = top12_32(x) & 0x7ff;
    if (abstop >= top12_32(128)) {
        @branchHint(.unlikely);
        // |x| >= 128 or x is nan
        if (@as(u32, @bitCast(x)) == @as(u32, @bitCast(-std.math.inf(f32))))
            return 0;

        if (abstop >= top12_32(std.math.inf(f32)))
            return x + x;

        if (x > 0)
            return 0x1p97 * 0x1p97;

        if (x <= -150)
            return 0x1p-95 * 0x1p-95;

        if (x < -149)
            return 0x1.4p-75 * 0x1.4p-75;
    }

    // x = k/N + r with r in [-1/(2N), 1/(2N)] and int k.
    var kd: f64 = xd + exp2_data.shift_scaled_32;
    const ki: u64 = @bitCast(kd);
    kd -= exp2_data.shift_scaled_32; // k/N for int k.
    const r: f64 = xd - kd;

    // exp2(x) = 2^(k/N) * 2^r ~= s * (C0*r^3 + C1*r^2 + C2*r + 1)
    var t: u64 = exp2_data.T_32[ki % 32];
    t +%= ki << 47;
    const s: f64 = @bitCast(t);
    const z: f64 = exp2_data.poly_32[0] * r + exp2_data.poly_32[1];
    const r2: f64 = r * r;
    var y: f64 = exp2_data.poly_32[2] * r + 1;
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
// negative k means the result may underflow.
inline fn specialcase(tmp: f64, sbits: u64, ki: u64) f64 {
    var scale: f64 = undefined;
    var y: f64 = undefined;
    if ((ki & 0x80000000) == 0) {
        // k > 0, the exponent of scale might have overflowed by 1.
        scale = @bitCast(sbits - (1 << 52));
        y = 2 * (scale + scale * tmp);
        return y;
    }

    // k < 0, need special care in the subnormal range.
    scale = @bitCast(sbits +% (1022 << 52));
    y = scale + scale * tmp;
    if (y < 1) {
        // Round y to the right precision before scaling it into the subnormal
        // range to avoid double rounding that can cause 0.5+E/2 ulp error where
        // E is the worst-case ulp error outside the subnormal range.  So this
        // is only useful if the goal is better than 1 ulp worst-case error.
        var lo: f64 = scale - y + scale * tmp;
        const hi: f64 = 1.0 + y;
        lo = 1.0 - hi + y + lo;
        y = (hi + lo) - 1.0;
        // Avoid -0.0 with downward rounding.
        if (y == 0)
            y = 0;

        // The underflow exception needs to be signaled explicitly.
        std.mem.doNotOptimizeAway(0x1p-1022 * 0x1p-1022);
    }

    y = 0x1p-1022 * y;
    return y;
}

// Top 12 bits of a double (sign and exponent bits).
inline fn top12_64(x: f64) u32 {
    return @truncate(@as(u64, @bitCast(x)) >> 52);
}

fn exp2_64(x: f64) f64 {
    var abstop: u32 = @truncate(top12_64(x) & 0x7ff);
    if (abstop -% top12_64(@as(f64, 0x1p-54)) >= top12_64(@as(f64, 512)) -% top12_64(@as(f64, 0x1p-54))) {
        @branchHint(.unlikely);
        if (abstop -% top12_64(@as(f64, 0x1p-54)) >= 0x80000000) {
            // Avoid spurious underflow for tiny x.
            // Note: 0 is common input.
            return 1 + x;
        }

        if (abstop >= top12_64(@as(f64, 1024))) {
            if (@as(u64, @bitCast(x)) == @as(u64, @bitCast(-std.math.inf(f64))))
                return 0;

            if (abstop >= top12_64(@as(f64, @bitCast(std.math.inf(f64)))))
                return 1 + x;

            if (@as(u64, @bitCast(x)) >> 63 == 0) {
                return 0x1p769 * 0x1p769;
            } else if (@as(u64, @bitCast(x)) >= @as(u64, @bitCast(@as(f64, -1075)))) {
                return 0x1p-767 * 0x1p-767;
            }
        }

        if (2 *% @as(u64, @bitCast(x)) > 2 *% @as(u64, @bitCast(@as(f64, 928)))) {
            // Large x is special cased below.
            abstop = 0;
        }
    }

    // exp2(x) = 2^(k/N) * 2^r, with 2^r in [2^(-1/2N),2^(1/2N)].
    // x = k/N + r, with int k and r in [-1/2N, 1/2N].
    var kd: f64 = x + exp_data.exp2_shift;
    const ki: u64 = @bitCast(kd); // k.
    kd -= exp_data.exp2_shift; // k/N for int k.
    const r: f64 = x - kd;
    // 2^(k/N) ~= scale * (1 + tail).
    const idx: u64 = 2 * (ki % 128);
    const top: u64 = ki << 45;
    const tail: f64 = @bitCast(exp_data.T_64[idx]);
    // This is only a valid scale when -1023*N < k < 1024*N.
    const sbits: u64 = exp_data.T_64[idx + 1] +% top;
    // exp2(x) = 2^(k/N) * 2^r ~= scale + scale * (tail + 2^r - 1).
    // Evaluation is optimized assuming superscalar pipelined execution.
    const r2: f64 = r * r;
    // Without fma the worst case error is 0.5/N ulp larger.
    // Worst case error is less than 0.5+0.86/N+(abs poly error * 2^53) ulp.
    const tmp: f64 = tail + r * exp_data.exp2_poly[0] + r2 * (exp_data.exp2_poly[1] + r * exp_data.exp2_poly[2]) + r2 * r2 * (exp_data.exp2_poly[3] + r * exp_data.exp2_poly[4]);

    if (abstop == 0) {
        @branchHint(.unlikely);
        return specialcase(tmp, sbits, ki);
    }

    const scale: f64 = @bitCast(sbits);
    // Note: tmp == 0 or |tmp| > 2^-65 and scale > 2^-928, so there
    // is no spurious underflow here even without fma.
    return scale + scale * tmp;
}

fn exp2_128(x: f128) f128 {
    var u: ldbl128.ieee_f128_shape2 = @bitCast(x);

    // Filter out exceptional cases.
    const hx: u32 = @intCast(u.expsign);
    const ix: u32 = hx & ((16384 - 1) + 16384);
    if (ix >= (16384 - 1) + 14) { // |x| >= 16384
        if (ix == (16384 - 1) + 16384) {
            if (u.manh != 0 or u.manl != 0 or (hx & 0x8000) == 0) {
                return x + x; // x is NaN or +Inf
            } else {
                return 0; // x is -Inf
            }
        }

        if (x >= 16384)
            return exp2_data.huge_128 * exp2_data.huge_128; // overflow

        if (x <= -16495)
            return exp2_data.twom10000_128 * exp2_data.twom10000_128; // underflow

    } else if (ix <= (16384 - 1) - 115) { // |x| < 0x1p-115
        return 1 + x;
    }

    // Reduce x, computing z, I0, and k. The low bits of x + redux
    // contain the 16-bit integer part of the exponent (k) followed by
    // TBLBITS fractional bits (I0). We use bit tricks to extract these
    // as integers, then set z to the remainder.
    //
    // Example: Suppose x is 0xabc.123456p0 and TBLBITS is 8.
    // Then the low-order word of x + redux is 0x000abc12,
    // We split this into k = 0xabc and I0 = 0x12 (adjusted to
    // index into the table), then we compute z = 0x0.003456p0.
    //
    // XXX If the exponent is negative, the computation of k depends on
    //     '>>' doing sign extension.
    u = @bitCast(x + exp2_data.redux_128);
    var I0: u32 = @truncate((u.manl & 0xffffffff) + 128 / 2);
    var k: i32 = undefined;
    {
        @setRuntimeSafety(false);
        k = @as(i32, @intCast(I0)) >> 7;
    }
    I0 = I0 & (128 - 1);
    u = @bitCast(@as(f128, @bitCast(u)) - exp2_data.redux_128);
    var z: f128 = x - @as(f128, @bitCast(u));
    var v: ldbl128.ieee_f128_shape2 = undefined;
    v.manh = 0;
    v.manl = 0;
    var twopk: f128 = undefined;
    var twopkp10000: f128 = undefined;
    if (k >= -16381) {
        @setRuntimeSafety(false);
        v.expsign = @intCast(16384 - 1 + k);
        twopk = @bitCast(v);
    } else {
        @setRuntimeSafety(false);
        v.expsign = @intCast(16384 - 1 + k + 10000);
        twopkp10000 = @bitCast(v);
    }

    // Compute r = exp2(y) = exp2t[I0] * p(z - eps[i]).
    const t: f128 = exp2_data.tbl_128[I0]; // exp2t[I0]
    z -= exp2_data.eps_128[I0]; // eps[I0]
    const r: f128 = t + t * z * (exp2_data.P1_128 + z * (exp2_data.P2_128 + z * (exp2_data.P3_128 + z * (exp2_data.P4_128 + z * (exp2_data.P5_128 + z * (exp2_data.P6_128 + z * (exp2_data.P7_128 + z * (exp2_data.P8_128 + z * (exp2_data.P9_128 + z * exp2_data.P10_128)))))))));

    // Scale by 2**k.
    if (k >= -16381) {
        if (k == 16384)
            return r * 2.0 * 0x1p16383;

        return r * twopk;
    } else {
        return r * twopkp10000 * exp2_data.twom10000_128;
    }
}

test exp2 {
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x4p+8, exp2(@as(f32, 0xap+0)));
    try std.testing.expectEqual(0x8p-4, exp2(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0xf.424p+16)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0xf.424p+16)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.ae89fap+0, exp2(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x1.6a09e6p+100, exp2(@as(f32, 0x6.48p+4)));
    try std.testing.expectEqual(0xb.504f3p-120, exp2(@as(f32, -0x7.48p+4)));
    try std.testing.expectEqual(0x1.6a09e6p-124, exp2(@as(f32, -0x7.b8p+4)));
    try std.testing.expectEqual(0xb.504f3p-128, exp2(@as(f32, -0x7.c8p+4)));
    try std.testing.expectEqual(0x5.a82798p-128, exp2(@as(f32, -0x7.d8p+4)));
    try std.testing.expectEqual(0x8p+124, exp2(@as(f32, 0x7.fp+4)));
    try std.testing.expectEqual(0x8p-152, exp2(@as(f32, -0x9.5p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x3.e84p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fb8p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fc8p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fd8p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x3.ffp+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x4.32p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x3.fffp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x4.01p+12)));
    try std.testing.expectEqual(0x3.ab0318p-128, exp2(@as(f32, -0x7.e2p+4)));
    try std.testing.expectEqual(0x3.5d13fp-128, exp2(@as(f32, -0x7.e4p+4)));
    try std.testing.expectEqual(0x3.159ca8p-128, exp2(@as(f32, -0x7.e6p+4)));
    try std.testing.expectEqual(0x2.d413dp-128, exp2(@as(f32, -0x7.e8p+4)));
    try std.testing.expectEqual(0x2.97fb58p-128, exp2(@as(f32, -0x7.eap+4)));
    try std.testing.expectEqual(0x2.60dfcp-128, exp2(@as(f32, -0x7.ecp+4)));
    try std.testing.expectEqual(0x2.2e5708p-128, exp2(@as(f32, -0x7.eep+4)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fe2p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fe4p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fe6p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fe8p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.feap+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fecp+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.feep+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fe4e8p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fe513p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffe2p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffe4p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffe6p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffe8p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffeap+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffecp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffeep+12)));
    try std.testing.expectEqual(0x1.002c6p+0, exp2(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0xf.fd3a7p-4, exp2(@as(f32, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000cp+0, exp2(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0xf.ffff5p-4, exp2(@as(f32, -0x1p-20)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x4p-32)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x1p-40)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0xf.fffa7p+124, exp2(@as(f32, 0x7.fffff8p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x8.00001p+4)));
    try std.testing.expectEqual(0x3.fffeap-128, exp2(@as(f32, -0x7.e00008p+4)));
    try std.testing.expectEqual(0x4.00016p-128, exp2(@as(f32, -0x7.dffff8p+4)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x3.fffffcp+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4.000008p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fep+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fe0004p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fdfffcp+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fep+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x3.fffffcp+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.c9p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.c90004p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.c8fffcp+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.c9p+8)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x3.fffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4.000008p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffep+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffe004p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffdffcp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffep+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fffp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fff004p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffeffcp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.fffp+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x3.fffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4.000008p+12)));
    try std.testing.expectEqual(std.math.inf(f32), exp2(@as(f32, 0x4p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffep+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffe004p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffdffcp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x3.ffep+12)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x1.a44722p+0, exp2(@as(f32, 0xb.71754p-4)));
    try std.testing.expectEqual(0x3.959e68p+12, exp2(@as(f32, 0xd.d77dp+0)));
    try std.testing.expectEqual(0x1.afdd74p+0, exp2(@as(f32, 0xc.122c4p-4)));
    try std.testing.expectEqual(0x6.546d6p-4, exp2(@as(f32, -0x1.567cc8p+0)));
    try std.testing.expectEqual(0x4.cfe008p-4, exp2(@as(f32, -0x1.bbbd76p+0)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f32, -0x1.3045fep+8)));
    try std.testing.expectEqual(0x5.c6bfd8p+8, exp2(@as(f32, 0xa.87b8bp+0)));
    try std.testing.expectEqual(0x8.a8745p-4, exp2(@as(f32, -0xe.2ce69p-4)));
    try std.testing.expectEqual(0xf.ff79bp-4, exp2(@as(f32, -0xc.1bf12p-16)));
    try std.testing.expectEqual(0xd.23272p-4, exp2(@as(f32, -0x4.8ce878p-4)));
    try std.testing.expectEqual(0x1.f6b64ap+0, exp2(@as(f32, 0xf.93d19p-4)));
    try std.testing.expectEqual(0x1.f6b64ap+0, exp2(@as(f32, 0xf.93d18p-4)));

    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x4p+8, exp2(@as(f64, 0xap+0)));
    try std.testing.expectEqual(0x8p-4, exp2(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0xf.424p+16)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0xf.424p+16)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.ae89f995ad3adp+0, exp2(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcdp+100, exp2(@as(f64, 0x6.48p+4)));
    try std.testing.expectEqual(0xb.504f333f9de68p-120, exp2(@as(f64, -0x7.48p+4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcdp-124, exp2(@as(f64, -0x7.b8p+4)));
    try std.testing.expectEqual(0xb.504f333f9de68p-128, exp2(@as(f64, -0x7.c8p+4)));
    try std.testing.expectEqual(0x5.a827999fcef34p-128, exp2(@as(f64, -0x7.d8p+4)));
    try std.testing.expectEqual(0x8p+124, exp2(@as(f64, 0x7.fp+4)));
    try std.testing.expectEqual(0x8p-152, exp2(@as(f64, -0x9.5p+4)));
    try std.testing.expectEqual(0x1.306fe0a31b715p+1000, exp2(@as(f64, 0x3.e84p+8)));
    try std.testing.expectEqual(0x1.6a09e667f3bcdp-1020, exp2(@as(f64, -0x3.fb8p+8)));
    try std.testing.expectEqual(0xb.504f333f9de68p-1024, exp2(@as(f64, -0x3.fc8p+8)));
    try std.testing.expectEqual(0x5.a827999fcef34p-1024, exp2(@as(f64, -0x3.fd8p+8)));
    try std.testing.expectEqual(0x8p+1020, exp2(@as(f64, 0x3.ffp+8)));
    try std.testing.expectEqual(0x4p-1076, exp2(@as(f64, -0x4.32p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x3.fffp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x4.01p+12)));
    try std.testing.expectEqual(0x3.ab031b9f7490ep-128, exp2(@as(f64, -0x7.e2p+4)));
    try std.testing.expectEqual(0x3.5d13f32b5a75ap-128, exp2(@as(f64, -0x7.e4p+4)));
    try std.testing.expectEqual(0x3.159ca845541b6p-128, exp2(@as(f64, -0x7.e6p+4)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-128, exp2(@as(f64, -0x7.e8p+4)));
    try std.testing.expectEqual(0x2.97fb5aa6c544ep-128, exp2(@as(f64, -0x7.eap+4)));
    try std.testing.expectEqual(0x2.60dfc14636e2ap-128, exp2(@as(f64, -0x7.ecp+4)));
    try std.testing.expectEqual(0x2.2e57078faa2f6p-128, exp2(@as(f64, -0x7.eep+4)));
    try std.testing.expectEqual(0x3.ab031b9f7491p-1024, exp2(@as(f64, -0x3.fe2p+8)));
    try std.testing.expectEqual(0x3.5d13f32b5a75cp-1024, exp2(@as(f64, -0x3.fe4p+8)));
    try std.testing.expectEqual(0x3.159ca845541b8p-1024, exp2(@as(f64, -0x3.fe6p+8)));
    try std.testing.expectEqual(0x2.d413cccfe7798p-1024, exp2(@as(f64, -0x3.fe8p+8)));
    try std.testing.expectEqual(0x2.97fb5aa6c545p-1024, exp2(@as(f64, -0x3.feap+8)));
    try std.testing.expectEqual(0x2.60dfc14636e2cp-1024, exp2(@as(f64, -0x3.fecp+8)));
    try std.testing.expectEqual(0x2.2e57078faa2f4p-1024, exp2(@as(f64, -0x3.feep+8)));
    try std.testing.expectEqual(0x3.3bed4179f82bcp-1024, exp2(@as(f64, -0x3.fe4e8p+8)));
    try std.testing.expectEqual(0x3.35ec906f22fbcp-1024, exp2(@as(f64, -0x3.fe513p+8)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffe2p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffe4p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffe6p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffe8p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffeap+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffecp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffeep+12)));
    try std.testing.expectEqual(0x1.002c605e2e8cfp+0, exp2(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0xf.fd3a751c0f7ep-4, exp2(@as(f64, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000b1721bdp+0, exp2(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0xf.ffff4e8debep-4, exp2(@as(f64, -0x1p-20)));
    try std.testing.expectEqual(0x1.00000002c5c86p+0, exp2(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffffd3a37ap-4, exp2(@as(f64, -0x4p-32)));
    try std.testing.expectEqual(0x1.0000000000b17p+0, exp2(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0xf.fffffffff4e9p-4, exp2(@as(f64, -0x1p-40)));
    try std.testing.expectEqual(0x1.0000000000003p+0, exp2(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffffdp-4, exp2(@as(f64, -0x4p-52)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0xf.fffa7470363f8p+124, exp2(@as(f64, 0x7.fffff8p+4)));
    try std.testing.expectEqual(0x1.0000b17255776p+128, exp2(@as(f64, 0x8.00001p+4)));
    try std.testing.expectEqual(0x3.fffe9d1c0d8fep-128, exp2(@as(f64, -0x7.e00008p+4)));
    try std.testing.expectEqual(0x4.000162e46d6f4p-128, exp2(@as(f64, -0x7.dffff8p+4)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814e8p+1020, exp2(@as(f64, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9dp+1020, exp2(@as(f64, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4.000008p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4p+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4.0000000000004p+8)));
    try std.testing.expectEqual(0x4p-1024, exp2(@as(f64, -0x3.fep+8)));
    try std.testing.expectEqual(0x3.fff4e8ede053cp-1024, exp2(@as(f64, -0x3.fe0004p+8)));
    try std.testing.expectEqual(0x3.ffffffffffa74p-1024, exp2(@as(f64, -0x3.fe00000000002p+8)));
    try std.testing.expectEqual(0x4.000b1730df6a4p-1024, exp2(@as(f64, -0x3.fdfffcp+8)));
    try std.testing.expectEqual(0x4p-1024, exp2(@as(f64, -0x3.fep+8)));
    try std.testing.expectEqual(0x4.000000000058cp-1024, exp2(@as(f64, -0x3.fdffffffffffep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814e8p+1020, exp2(@as(f64, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9dp+1020, exp2(@as(f64, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814e8p+1020, exp2(@as(f64, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9dp+1020, exp2(@as(f64, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(0x8p-972, exp2(@as(f64, -0x3.c9p+8)));
    try std.testing.expectEqual(0x7.ffe9d1dbc0a74p-972, exp2(@as(f64, -0x3.c90004p+8)));
    try std.testing.expectEqual(0x7.ffffffffff4e8p-972, exp2(@as(f64, -0x3.c900000000002p+8)));
    try std.testing.expectEqual(0x8.00162e61bed48p-972, exp2(@as(f64, -0x3.c8fffcp+8)));
    try std.testing.expectEqual(0x8p-972, exp2(@as(f64, -0x3.c9p+8)));
    try std.testing.expectEqual(0x8.0000000000b18p-972, exp2(@as(f64, -0x3.c8ffffffffffep+8)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x3.fffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x3.ffffffffffffep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4.000008p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4.0000000000004p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffep+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffe004p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffe0000000002p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffdffcp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffep+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffdfffffffffep+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.fffp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.fff004p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.fff0000000002p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffeffcp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.fffp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffefffffffffep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x3.fffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x3.ffffffffffffep+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4.000008p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f64), exp2(@as(f64, 0x4.0000000000004p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffep+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffe004p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffe0000000002p+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffdffcp+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffep+12)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f64, -0x3.ffdfffffffffep+12)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x1.a44722ff862d7p+0, exp2(@as(f64, 0xb.71754p-4)));
    try std.testing.expectEqual(0x3.959e67fd7ff86p+12, exp2(@as(f64, 0xd.d77dp+0)));
    // try std.testing.expectEqual(0x1.afdd736c287abp+0, exp2(@as(f64, 0xc.122c4p-4)));
    try std.testing.expectEqual(0x6.546d5ccd21bacp-4, exp2(@as(f64, -0x1.567cc8p+0)));
    try std.testing.expectEqual(0x4.cfe0085ef004cp-4, exp2(@as(f64, -0x1.bbbd76p+0)));
    try std.testing.expectEqual(0xd.3ce16388003dp-308, exp2(@as(f64, -0x1.3045fep+8)));
    try std.testing.expectEqual(0x5.c6bfd7fd625f8p+8, exp2(@as(f64, 0xa.87b8bp+0)));
    try std.testing.expectEqual(0x8.a8744fff686fp-4, exp2(@as(f64, -0xe.2ce69p-4)));
    try std.testing.expectEqual(0xf.ff79b6bee6bd8p-4, exp2(@as(f64, -0xc.1bf12p-16)));
    try std.testing.expectEqual(0xd.23271e170998p-4, exp2(@as(f64, -0x4.8ce878p-4)));
    try std.testing.expectEqual(0x1.f6b64a6870e6bp+0, exp2(@as(f64, 0xf.93d19p-4)));
    try std.testing.expectEqual(0x1.f6b6490bfcd17p+0, exp2(@as(f64, 0xf.93d18p-4)));
    try std.testing.expectEqual(0x1.f6b64a10a015fp+0, exp2(@as(f64, 0xf.93d18bf7be8d8p-4)));
    try std.testing.expectEqual(0x1.f6b64a10a015ep+0, exp2(@as(f64, 0xf.93d18bf7be8dp-4)));

    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x4p+8, exp2(@as(f80, 0xap+0)));
    try std.testing.expectEqual(0x8p-4, exp2(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0xf.424p+16)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f80, -0xf.424p+16)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.ae89f995ad3ad5e8p+0, exp2(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908p+100, exp2(@as(f80, 0x6.48p+4)));
    try std.testing.expectEqual(0xb.504f333f9de6484p-120, exp2(@as(f80, -0x7.48p+4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908p-124, exp2(@as(f80, -0x7.b8p+4)));
    try std.testing.expectEqual(0xb.504f333f9de6484p-128, exp2(@as(f80, -0x7.c8p+4)));
    try std.testing.expectEqual(0x5.a827999fcef3242p-128, exp2(@as(f80, -0x7.d8p+4)));
    try std.testing.expectEqual(0x8p+124, exp2(@as(f80, 0x7.fp+4)));
    try std.testing.expectEqual(0x8p-152, exp2(@as(f80, -0x9.5p+4)));
    try std.testing.expectEqual(0x1.306fe0a31b7152dep+1000, exp2(@as(f80, 0x3.e84p+8)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908p-1020, exp2(@as(f80, -0x3.fb8p+8)));
    try std.testing.expectEqual(0xb.504f333f9de6484p-1024, exp2(@as(f80, -0x3.fc8p+8)));
    try std.testing.expectEqual(0x5.a827999fcef3242p-1024, exp2(@as(f80, -0x3.fd8p+8)));
    try std.testing.expectEqual(0x8p+1020, exp2(@as(f80, 0x3.ffp+8)));
    try std.testing.expectEqual(0x4p-1076, exp2(@as(f80, -0x4.32p+8)));
    try std.testing.expectEqual(0x8p+16380, exp2(@as(f80, 0x3.fffp+12)));
    try std.testing.expectEqual(0x1p-16400, exp2(@as(f80, -0x4.01p+12)));
    try std.testing.expectEqual(0x3.ab031b9f7490e4bcp-128, exp2(@as(f80, -0x7.e2p+4)));
    try std.testing.expectEqual(0x3.5d13f32b5a75abdp-128, exp2(@as(f80, -0x7.e4p+4)));
    try std.testing.expectEqual(0x3.159ca845541b6b74p-128, exp2(@as(f80, -0x7.e6p+4)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-128, exp2(@as(f80, -0x7.e8p+4)));
    try std.testing.expectEqual(0x2.97fb5aa6c544e3a8p-128, exp2(@as(f80, -0x7.eap+4)));
    try std.testing.expectEqual(0x2.60dfc14636e2a5bcp-128, exp2(@as(f80, -0x7.ecp+4)));
    try std.testing.expectEqual(0x2.2e57078faa2f5b9cp-128, exp2(@as(f80, -0x7.eep+4)));
    try std.testing.expectEqual(0x3.ab031b9f7490e4bcp-1024, exp2(@as(f80, -0x3.fe2p+8)));
    try std.testing.expectEqual(0x3.5d13f32b5a75abdp-1024, exp2(@as(f80, -0x3.fe4p+8)));
    try std.testing.expectEqual(0x3.159ca845541b6b74p-1024, exp2(@as(f80, -0x3.fe6p+8)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-1024, exp2(@as(f80, -0x3.fe8p+8)));
    try std.testing.expectEqual(0x2.97fb5aa6c544e3a8p-1024, exp2(@as(f80, -0x3.feap+8)));
    try std.testing.expectEqual(0x2.60dfc14636e2a5bcp-1024, exp2(@as(f80, -0x3.fecp+8)));
    try std.testing.expectEqual(0x2.2e57078faa2f5b9cp-1024, exp2(@as(f80, -0x3.feep+8)));
    try std.testing.expectEqual(0x3.3bed4179f82bc004p-1024, exp2(@as(f80, -0x3.fe4e8p+8)));
    try std.testing.expectEqual(0x3.35ec906f22fbcp-1024, exp2(@as(f80, -0x3.fe513p+8)));
    try std.testing.expectEqual(0x3.ab031b9f7490e4b8p-16384, exp2(@as(f80, -0x3.ffe2p+12)));
    try std.testing.expectEqual(0x3.5d13f32b5a75abdp-16384, exp2(@as(f80, -0x3.ffe4p+12)));
    try std.testing.expectEqual(0x3.159ca845541b6b78p-16384, exp2(@as(f80, -0x3.ffe6p+12)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-16384, exp2(@as(f80, -0x3.ffe8p+12)));
    try std.testing.expectEqual(0x2.97fb5aa6c544e3a8p-16384, exp2(@as(f80, -0x3.ffeap+12)));
    try std.testing.expectEqual(0x2.60dfc14636e2a5cp-16384, exp2(@as(f80, -0x3.ffecp+12)));
    try std.testing.expectEqual(0x2.2e57078faa2f5b98p-16384, exp2(@as(f80, -0x3.ffeep+12)));
    try std.testing.expectEqual(0x1.002c605e2e8cec5p+0, exp2(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0xf.fd3a751c0f7e10cp-4, exp2(@as(f80, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000b1721bcfc9ap+0, exp2(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0xf.ffff4e8debe025ep-4, exp2(@as(f80, -0x1p-20)));
    try std.testing.expectEqual(0x1.00000002c5c85fe4p+0, exp2(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffffd3a37a025p-4, exp2(@as(f80, -0x4p-32)));
    try std.testing.expectEqual(0x1.0000000000b17218p+0, exp2(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0xf.fffffffff4e8de8p-4, exp2(@as(f80, -0x1p-40)));
    try std.testing.expectEqual(0x1.0000000000002c5cp+0, exp2(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffffd3a3p-4, exp2(@as(f80, -0x4p-52)));
    try std.testing.expectEqual(0x1.000000000000000cp+0, exp2(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0xf.ffffffffffffff5p-4, exp2(@as(f80, -0x1p-60)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x1p-10000)));
    try std.testing.expectEqual(0xf.fffa7470363f451p+124, exp2(@as(f80, 0x7.fffff8p+4)));
    try std.testing.expectEqual(0x1.0000b17255775c04p+128, exp2(@as(f80, 0x8.00001p+4)));
    try std.testing.expectEqual(0x3.fffe9d1c0d8fd144p-128, exp2(@as(f80, -0x7.e00008p+4)));
    try std.testing.expectEqual(0x4.000162e46d6f26b8p-128, exp2(@as(f80, -0x7.dffff8p+4)));
    try std.testing.expectEqual(0x1p+1024, exp2(@as(f80, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814eb54p+1020, exp2(@as(f80, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d1bdp+1020, exp2(@as(f80, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(0x1.00058ba01fb9f96ep+1024, exp2(@as(f80, 0x4.000008p+8)));
    try std.testing.expectEqual(0x1p+1024, exp2(@as(f80, 0x4p+8)));
    try std.testing.expectEqual(0x1.00000000002c5c86p+1024, exp2(@as(f80, 0x4.0000000000004p+8)));
    try std.testing.expectEqual(0x4p-1024, exp2(@as(f80, -0x3.fep+8)));
    try std.testing.expectEqual(0x3.fff4e8ede053ad5p-1024, exp2(@as(f80, -0x3.fe0004p+8)));
    try std.testing.expectEqual(0x3.ffffffffffa746f4p-1024, exp2(@as(f80, -0x3.fe00000000002p+8)));
    try std.testing.expectEqual(0x4.000b1730df6a5248p-1024, exp2(@as(f80, -0x3.fdfffcp+8)));
    try std.testing.expectEqual(0x4p-1024, exp2(@as(f80, -0x3.fep+8)));
    try std.testing.expectEqual(0x4.000000000058b908p-1024, exp2(@as(f80, -0x3.fdffffffffffep+8)));
    try std.testing.expectEqual(0x1p+1024, exp2(@as(f80, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814eb54p+1020, exp2(@as(f80, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d1bdp+1020, exp2(@as(f80, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffd3ap+1020, exp2(@as(f80, 0x3.fffffffffffffffcp+8)));
    try std.testing.expectEqual(0xf.ffffffffffffa74p+1020, exp2(@as(f80, 0x3.fffffffffffffff8p+8)));
    try std.testing.expectEqual(0x1p+1024, exp2(@as(f80, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814eb54p+1020, exp2(@as(f80, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d1bdp+1020, exp2(@as(f80, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffd3ap+1020, exp2(@as(f80, 0x3.fffffffffffffffcp+8)));
    try std.testing.expectEqual(0xf.ffffffffffffa74p+1020, exp2(@as(f80, 0x3.fffffffffffffff8p+8)));
    try std.testing.expectEqual(0x8p-972, exp2(@as(f80, -0x3.c9p+8)));
    try std.testing.expectEqual(0x7.ffe9d1dbc0a75aap-972, exp2(@as(f80, -0x3.c90004p+8)));
    try std.testing.expectEqual(0x7.ffffffffff4e8de8p-972, exp2(@as(f80, -0x3.c900000000002p+8)));
    try std.testing.expectEqual(0x7.ffffffffffffe9dp-972, exp2(@as(f80, -0x3.c900000000000004p+8)));
    try std.testing.expectEqual(0x8.00162e61bed4a49p-972, exp2(@as(f80, -0x3.c8fffcp+8)));
    try std.testing.expectEqual(0x8p-972, exp2(@as(f80, -0x3.c9p+8)));
    try std.testing.expectEqual(0x8.0000000000b1721p-972, exp2(@as(f80, -0x3.c8ffffffffffep+8)));
    try std.testing.expectEqual(0x8.000000000000163p-972, exp2(@as(f80, -0x3.c8fffffffffffffcp+8)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4p+12)));
    try std.testing.expectEqual(0xf.fd3a751c0f7e10cp+16380, exp2(@as(f80, 0x3.fffffcp+12)));
    try std.testing.expectEqual(0xf.ffffffffe9d1bdp+16380, exp2(@as(f80, 0x3.ffffffffffffep+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3a3p+16380, exp2(@as(f80, 0x3.fffffffffffffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4.000008p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4.0000000000004p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4.0000000000000008p+12)));
    try std.testing.expectEqual(0x4p-16384, exp2(@as(f80, -0x3.ffep+12)));
    try std.testing.expectEqual(0x3.ff4e9d4703df843p-16384, exp2(@as(f80, -0x3.ffe004p+12)));
    try std.testing.expectEqual(0x3.fffffffffa746f4p-16384, exp2(@as(f80, -0x3.ffe0000000002p+12)));
    try std.testing.expectEqual(0x3.ffffffffffff4e9p-16384, exp2(@as(f80, -0x3.ffe0000000000004p+12)));
    try std.testing.expectEqual(0x4.00b18178ba33b14p-16384, exp2(@as(f80, -0x3.ffdffcp+12)));
    try std.testing.expectEqual(0x4p-16384, exp2(@as(f80, -0x3.ffep+12)));
    try std.testing.expectEqual(0x4.00000000058b90cp-16384, exp2(@as(f80, -0x3.ffdfffffffffep+12)));
    try std.testing.expectEqual(0x4.000000000000b17p-16384, exp2(@as(f80, -0x3.ffdffffffffffffcp+12)));
    try std.testing.expectEqual(0x2p-16384, exp2(@as(f80, -0x3.fffp+12)));
    try std.testing.expectEqual(0x1.ffa74ea381efc218p-16384, exp2(@as(f80, -0x3.fff004p+12)));
    try std.testing.expectEqual(0x1.fffffffffd3a37ap-16384, exp2(@as(f80, -0x3.fff0000000002p+12)));
    try std.testing.expectEqual(0x1.ffffffffffffa748p-16384, exp2(@as(f80, -0x3.fff0000000000004p+12)));
    try std.testing.expectEqual(0x2.0058c0bc5d19d8ap-16384, exp2(@as(f80, -0x3.ffeffcp+12)));
    try std.testing.expectEqual(0x2p-16384, exp2(@as(f80, -0x3.fffp+12)));
    try std.testing.expectEqual(0x2.0000000002c5c86p-16384, exp2(@as(f80, -0x3.ffefffffffffep+12)));
    try std.testing.expectEqual(0x2.00000000000058b8p-16384, exp2(@as(f80, -0x3.ffeffffffffffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4p+12)));
    try std.testing.expectEqual(0xf.fd3a751c0f7e10cp+16380, exp2(@as(f80, 0x3.fffffcp+12)));
    try std.testing.expectEqual(0xf.ffffffffe9d1bdp+16380, exp2(@as(f80, 0x3.ffffffffffffep+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3a3p+16380, exp2(@as(f80, 0x3.fffffffffffffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4.000008p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4.0000000000004p+12)));
    try std.testing.expectEqual(std.math.inf(f80), exp2(@as(f80, 0x4.0000000000000008p+12)));
    try std.testing.expectEqual(0x4p-16384, exp2(@as(f80, -0x3.ffep+12)));
    try std.testing.expectEqual(0x3.ff4e9d4703df843p-16384, exp2(@as(f80, -0x3.ffe004p+12)));
    try std.testing.expectEqual(0x3.fffffffffa746f4p-16384, exp2(@as(f80, -0x3.ffe0000000002p+12)));
    try std.testing.expectEqual(0x3.ffffffffffff4e9p-16384, exp2(@as(f80, -0x3.ffe0000000000004p+12)));
    try std.testing.expectEqual(0x4.00b18178ba33b14p-16384, exp2(@as(f80, -0x3.ffdffcp+12)));
    try std.testing.expectEqual(0x4p-16384, exp2(@as(f80, -0x3.ffep+12)));
    try std.testing.expectEqual(0x4.00000000058b90cp-16384, exp2(@as(f80, -0x3.ffdfffffffffep+12)));
    try std.testing.expectEqual(0x4.000000000000b17p-16384, exp2(@as(f80, -0x3.ffdffffffffffffcp+12)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x1.a44722ff862d7436p+0, exp2(@as(f80, 0xb.71754p-4)));
    try std.testing.expectEqual(0x3.959e67fd7ff858c4p+12, exp2(@as(f80, 0xd.d77dp+0)));
    try std.testing.expectEqual(0x1.afdd736c287aa8p+0, exp2(@as(f80, 0xc.122c4p-4)));
    try std.testing.expectEqual(0x6.546d5ccd21bad058p-4, exp2(@as(f80, -0x1.567cc8p+0)));
    try std.testing.expectEqual(0x4.cfe0085ef004d24p-4, exp2(@as(f80, -0x1.bbbd76p+0)));
    try std.testing.expectEqual(0xd.3ce16388003d33ap-308, exp2(@as(f80, -0x1.3045fep+8)));
    try std.testing.expectEqual(0x5.c6bfd7fd625f812p+8, exp2(@as(f80, 0xa.87b8bp+0)));
    try std.testing.expectEqual(0x8.a8744fff686ede8p-4, exp2(@as(f80, -0xe.2ce69p-4)));
    try std.testing.expectEqual(0xf.ff79b6bee6bd8p-4, exp2(@as(f80, -0xc.1bf12p-16)));
    try std.testing.expectEqual(0xd.23271e170998p-4, exp2(@as(f80, -0x4.8ce878p-4)));
    try std.testing.expectEqual(0x1.f6b64a6870e6ae12p+0, exp2(@as(f80, 0xf.93d19p-4)));
    try std.testing.expectEqual(0x1.f6b6490bfcd17676p+0, exp2(@as(f80, 0xf.93d18p-4)));
    try std.testing.expectEqual(0x1.f6b64a10a015e998p+0, exp2(@as(f80, 0xf.93d18bf7be8d8p-4)));
    try std.testing.expectEqual(0x1.f6b64a10a015deb4p+0, exp2(@as(f80, 0xf.93d18bf7be8dp-4)));
    try std.testing.expectEqual(0x1.f6b64a10a015e208p+0, exp2(@as(f80, 0xf.93d18bf7be8d272p-4)));

    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x4p+8, exp2(@as(f128, 0xap+0)));
    try std.testing.expectEqual(0x8p-4, exp2(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0xf.424p+16)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f128, -0xf.424p+16)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, exp2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.ae89f995ad3ad5e8734d1773205ap+0, exp2(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908b2fb1366ea95p+100, exp2(@as(f128, 0x6.48p+4)));
    try std.testing.expectEqual(0xb.504f333f9de6484597d89b3754a8p-120, exp2(@as(f128, -0x7.48p+4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908b2fb1366ea95p-124, exp2(@as(f128, -0x7.b8p+4)));
    try std.testing.expectEqual(0xb.504f333f9de6484597d89b3754a8p-128, exp2(@as(f128, -0x7.c8p+4)));
    try std.testing.expectEqual(0x5.a827999fcef32422cbec4d9baa54p-128, exp2(@as(f128, -0x7.d8p+4)));
    try std.testing.expectEqual(0x8p+124, exp2(@as(f128, 0x7.fp+4)));
    try std.testing.expectEqual(0x8p-152, exp2(@as(f128, -0x9.5p+4)));
    try std.testing.expectEqual(0x1.306fe0a31b7152de8d5a46305c86p+1000, exp2(@as(f128, 0x3.e84p+8)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908b2fb1366ea95p-1020, exp2(@as(f128, -0x3.fb8p+8)));
    try std.testing.expectEqual(0xb.504f333f9de6484597d89b3754a8p-1024, exp2(@as(f128, -0x3.fc8p+8)));
    try std.testing.expectEqual(0x5.a827999fcef32422cbec4d9baa54p-1024, exp2(@as(f128, -0x3.fd8p+8)));
    try std.testing.expectEqual(0x8p+1020, exp2(@as(f128, 0x3.ffp+8)));
    try std.testing.expectEqual(0x4p-1076, exp2(@as(f128, -0x4.32p+8)));
    try std.testing.expectEqual(0x8p+16380, exp2(@as(f128, 0x3.fffp+12)));
    try std.testing.expectEqual(0x1p-16400, exp2(@as(f128, -0x4.01p+12)));
    try std.testing.expectEqual(0x3.ab031b9f7490e4bb40b5d6cdc1bap-128, exp2(@as(f128, -0x7.e2p+4)));
    try std.testing.expectEqual(0x3.5d13f32b5a75abd0e69a2ee640b4p-128, exp2(@as(f128, -0x7.e4p+4)));
    try std.testing.expectEqual(0x3.159ca845541b6b74f8ab43259376p-128, exp2(@as(f128, -0x7.e6p+4)));
    try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-128, exp2(@as(f128, -0x7.e8p+4)));
    try std.testing.expectEqual(0x2.97fb5aa6c544e3a872f5fd885c42p-128, exp2(@as(f128, -0x7.eap+4)));
    try std.testing.expectEqual(0x2.60dfc14636e2a5bd1ab48c60b90cp-128, exp2(@as(f128, -0x7.ecp+4)));
    try std.testing.expectEqual(0x2.2e57078faa2f5b9bef918a1d6294p-128, exp2(@as(f128, -0x7.eep+4)));
    try std.testing.expectEqual(0x3.ab031b9f7490e4bb40b5d6cdc1bap-1024, exp2(@as(f128, -0x3.fe2p+8)));
    try std.testing.expectEqual(0x3.5d13f32b5a75abd0e69a2ee640b4p-1024, exp2(@as(f128, -0x3.fe4p+8)));
    try std.testing.expectEqual(0x3.159ca845541b6b74f8ab43259376p-1024, exp2(@as(f128, -0x3.fe6p+8)));
    try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-1024, exp2(@as(f128, -0x3.fe8p+8)));
    try std.testing.expectEqual(0x2.97fb5aa6c544e3a872f5fd885c42p-1024, exp2(@as(f128, -0x3.feap+8)));
    try std.testing.expectEqual(0x2.60dfc14636e2a5bd1ab48c60b90cp-1024, exp2(@as(f128, -0x3.fecp+8)));
    try std.testing.expectEqual(0x2.2e57078faa2f5b9bef918a1d6294p-1024, exp2(@as(f128, -0x3.feep+8)));
    try std.testing.expectEqual(0x3.3bed4179f82bc002979648b91cfap-1024, exp2(@as(f128, -0x3.fe4e8p+8)));
    try std.testing.expectEqual(0x3.35ec906f22fbbffeffc0d272938p-1024, exp2(@as(f128, -0x3.fe513p+8)));
    try std.testing.expectEqual(0x3.ab031b9f7490e4bb40b5d6cdc1b8p-16384, exp2(@as(f128, -0x3.ffe2p+12)));
    try std.testing.expectEqual(0x3.5d13f32b5a75abd0e69a2ee640b4p-16384, exp2(@as(f128, -0x3.ffe4p+12)));
    try std.testing.expectEqual(0x3.159ca845541b6b74f8ab43259378p-16384, exp2(@as(f128, -0x3.ffe6p+12)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52cp-16384, exp2(@as(f128, -0x3.ffe8p+12)));
    try std.testing.expectEqual(0x2.97fb5aa6c544e3a872f5fd885c4p-16384, exp2(@as(f128, -0x3.ffeap+12)));
    try std.testing.expectEqual(0x2.60dfc14636e2a5bd1ab48c60b90cp-16384, exp2(@as(f128, -0x3.ffecp+12)));
    try std.testing.expectEqual(0x2.2e57078faa2f5b9bef918a1d6294p-16384, exp2(@as(f128, -0x3.ffeep+12)));
    try std.testing.expectEqual(0x1.002c605e2e8cec506d21bfc89a24p+0, exp2(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0xf.fd3a751c0f7e10bd3b9f8ae012f8p-4, exp2(@as(f128, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000b1721bcfc99d9f890ea0691p+0, exp2(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0xf.ffff4e8debe025e24128a3d4607p-4, exp2(@as(f128, -0x1p-20)));
    try std.testing.expectEqual(0x1.00000002c5c85fe31f35a6a30da2p+0, exp2(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffffd3a37a02490b9d93da3c2p-4, exp2(@as(f128, -0x4p-32)));
    try std.testing.expectEqual(0x1.0000000000b17217f7d20cf927c9p+0, exp2(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0xf.fffffffff4e8de8082e6e05d035p-4, exp2(@as(f128, -0x1p-40)));
    try std.testing.expectEqual(0x1.0000000000002c5c85fdf473e243p+0, exp2(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffffd3a37a020b8c256dp-4, exp2(@as(f128, -0x4p-52)));
    try std.testing.expectEqual(0x1.000000000000000b17217f7d1cf8p+0, exp2(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0xf.ffffffffffffff4e8de8082e3088p-4, exp2(@as(f128, -0x1p-60)));
    try std.testing.expectEqual(0x1.0000000000000000000000000b17p+0, exp2(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffff4e9p-4, exp2(@as(f128, -0x1p-100)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x1p-600)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x1p-10000)));
    try std.testing.expectEqual(0xf.fffa7470363f4515426d76c762b8p+124, exp2(@as(f128, 0x7.fffff8p+4)));
    try std.testing.expectEqual(0x1.0000b17255775c040618bf4a4adfp+128, exp2(@as(f128, 0x8.00001p+4)));
    try std.testing.expectEqual(0x3.fffe9d1c0d8fd145509b5db1d8aep-128, exp2(@as(f128, -0x7.e00008p+4)));
    try std.testing.expectEqual(0x4.000162e46d6f26b8bbb607a6df54p-128, exp2(@as(f128, -0x7.dffff8p+4)));
    try std.testing.expectEqual(0x1p+1024, exp2(@as(f128, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814eb53cd7629d70feap+1020, exp2(@as(f128, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d1bd0105c706c8768p+1020, exp2(@as(f128, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(0x1.00058ba01fb9f96d6cacd4b18091p+1024, exp2(@as(f128, 0x4.000008p+8)));
    try std.testing.expectEqual(0x1p+1024, exp2(@as(f128, 0x4p+8)));
    try std.testing.expectEqual(0x1.00000000002c5c85fdf477b662b2p+1024, exp2(@as(f128, 0x4.0000000000004p+8)));
    try std.testing.expectEqual(0x4p-1024, exp2(@as(f128, -0x3.fep+8)));
    try std.testing.expectEqual(0x3.fff4e8ede053ad4f35d8a75c3fa8p-1024, exp2(@as(f128, -0x3.fe0004p+8)));
    try std.testing.expectEqual(0x3.ffffffffffa746f404171c1b21dap-1024, exp2(@as(f128, -0x3.fe00000000002p+8)));
    try std.testing.expectEqual(0x4.000b1730df6a5247426170d231a4p-1024, exp2(@as(f128, -0x3.fdfffcp+8)));
    try std.testing.expectEqual(0x4p-1024, exp2(@as(f128, -0x3.fep+8)));
    try std.testing.expectEqual(0x4.000000000058b90bfbe8eb94cda4p-1024, exp2(@as(f128, -0x3.fdffffffffffep+8)));
    try std.testing.expectEqual(0x1p+1024, exp2(@as(f128, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814eb53cd7629d70feap+1020, exp2(@as(f128, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d1bd0105c706c8768p+1020, exp2(@as(f128, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffd3a37a020b8c21dp+1020, exp2(@as(f128, 0x3.fffffffffffffffcp+8)));
    try std.testing.expectEqual(0xf.ffffffffffffa746f4041718442p+1020, exp2(@as(f128, 0x3.fffffffffffffff8p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffbffffffffffd357p+1020, exp2(@as(f128, 0x3.fffffffffffffffa3aae26b51fp+8)));
    try std.testing.expectEqual(0x1p+1024, exp2(@as(f128, 0x4p+8)));
    try std.testing.expectEqual(0xf.ffd3a3b7814eb53cd7629d70feap+1020, exp2(@as(f128, 0x3.fffffcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d1bd0105c706c8768p+1020, exp2(@as(f128, 0x3.ffffffffffffep+8)));
    try std.testing.expectEqual(0xf.ffffffffffffd3a37a020b8c21dp+1020, exp2(@as(f128, 0x3.fffffffffffffffcp+8)));
    try std.testing.expectEqual(0xf.ffffffffffffa746f4041718442p+1020, exp2(@as(f128, 0x3.fffffffffffffff8p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffc000000000084c9p+1020, exp2(@as(f128, 0x3.fffffffffffffffa3aae26b52p+8)));
    try std.testing.expectEqual(0x8p-972, exp2(@as(f128, -0x3.c9p+8)));
    try std.testing.expectEqual(0x7.ffe9d1dbc0a75a9e6bb14eb87f5p-972, exp2(@as(f128, -0x3.c90004p+8)));
    try std.testing.expectEqual(0x7.ffffffffff4e8de8082e383643b4p-972, exp2(@as(f128, -0x3.c900000000002p+8)));
    try std.testing.expectEqual(0x7.ffffffffffffe9d1bd0105c610e8p-972, exp2(@as(f128, -0x3.c900000000000004p+8)));
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffa747p-972, exp2(@as(f128, -0x3.c9000000000000000000000001p+8)));
    try std.testing.expectEqual(0x8.00162e61bed4a48e84c2e1a46348p-972, exp2(@as(f128, -0x3.c8fffcp+8)));
    try std.testing.expectEqual(0x8p-972, exp2(@as(f128, -0x3.c9p+8)));
    try std.testing.expectEqual(0x8.0000000000b17217f7d1d7299b48p-972, exp2(@as(f128, -0x3.c8ffffffffffep+8)));
    try std.testing.expectEqual(0x8.000000000000162e42fefa39ef58p-972, exp2(@as(f128, -0x3.c8fffffffffffffcp+8)));
    try std.testing.expectEqual(0x8.0000000000000000000000058b9p-972, exp2(@as(f128, -0x3.c8ffffffffffffffffffffffffp+8)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4p+12)));
    try std.testing.expectEqual(0xf.fd3a751c0f7e10bd3b9f8ae012f8p+16380, exp2(@as(f128, 0x3.fffffcp+12)));
    try std.testing.expectEqual(0xf.ffffffffe9d1bd0105d570a98688p+16380, exp2(@as(f128, 0x3.ffffffffffffep+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3a37a020b8c256dp+16380, exp2(@as(f128, 0x3.fffffffffffffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4.000008p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4.0000000000004p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4.0000000000000008p+12)));
    try std.testing.expectEqual(0x4p-16384, exp2(@as(f128, -0x3.ffep+12)));
    try std.testing.expectEqual(0x3.ff4e9d4703df842f4ee7e2b804cp-16384, exp2(@as(f128, -0x3.ffe004p+12)));
    try std.testing.expectEqual(0x3.fffffffffa746f4041755c2a61ap-16384, exp2(@as(f128, -0x3.ffe0000000002p+12)));
    try std.testing.expectEqual(0x3.ffffffffffff4e8de8082e3095b4p-16384, exp2(@as(f128, -0x3.ffe0000000000004p+12)));
    try std.testing.expectEqual(0x4.00b18178ba33b141b486ff22689p-16384, exp2(@as(f128, -0x3.ffdffcp+12)));
    try std.testing.expectEqual(0x4p-16384, exp2(@as(f128, -0x3.ffep+12)));
    try std.testing.expectEqual(0x4.00000000058b90bfbe9253c51e4p-16384, exp2(@as(f128, -0x3.ffdfffffffffep+12)));
    try std.testing.expectEqual(0x4.000000000000b17217f7d1cf890cp-16384, exp2(@as(f128, -0x3.ffdffffffffffffcp+12)));
    try std.testing.expectEqual(0x2p-16384, exp2(@as(f128, -0x3.fffp+12)));
    try std.testing.expectEqual(0x1.ffa74ea381efc217a773f15c026p-16384, exp2(@as(f128, -0x3.fff004p+12)));
    try std.testing.expectEqual(0x1.fffffffffd3a37a020baae1530dp-16384, exp2(@as(f128, -0x3.fff0000000002p+12)));
    // try std.testing.expectEqual(0x1.ffffffffffffa746f40417184adcp-16384, exp2(@as(f128, -0x3.fff0000000000004p+12)));
    try std.testing.expectEqual(0x2.0058c0bc5d19d8a0da437f913448p-16384, exp2(@as(f128, -0x3.ffeffcp+12)));
    try std.testing.expectEqual(0x2p-16384, exp2(@as(f128, -0x3.fffp+12)));
    try std.testing.expectEqual(0x2.0000000002c5c85fdf4929e28f2p-16384, exp2(@as(f128, -0x3.ffefffffffffep+12)));
    // try std.testing.expectEqual(0x2.00000000000058b90bfbe8e7c484p-16384, exp2(@as(f128, -0x3.ffeffffffffffffcp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4p+12)));
    try std.testing.expectEqual(0xf.fd3a751c0f7e10bd3b9f8ae012f8p+16380, exp2(@as(f128, 0x3.fffffcp+12)));
    try std.testing.expectEqual(0xf.ffffffffe9d1bd0105d570a98688p+16380, exp2(@as(f128, 0x3.ffffffffffffep+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3a37a020b8c256dp+16380, exp2(@as(f128, 0x3.fffffffffffffffcp+12)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffe9d18p+16380, exp2(@as(f128, 0x3.fffffffffffffffffffffffffffep+12)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffff4e8de8p+16380, exp2(@as(f128, 0x3.ffffffffffffffffffffffffffp+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4.000008p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4.0000000000004p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4.0000000000000008p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4.0000000000000000000000000004p+12)));
    try std.testing.expectEqual(std.math.inf(f128), exp2(@as(f128, 0x4.00000000000000000000000002p+12)));
    try std.testing.expectEqual(0x4p-16384, exp2(@as(f128, -0x3.ffep+12)));
    try std.testing.expectEqual(0x3.ff4e9d4703df842f4ee7e2b804cp-16384, exp2(@as(f128, -0x3.ffe004p+12)));
    try std.testing.expectEqual(0x3.fffffffffa746f4041755c2a61ap-16384, exp2(@as(f128, -0x3.ffe0000000002p+12)));
    try std.testing.expectEqual(0x3.ffffffffffff4e8de8082e3095b4p-16384, exp2(@as(f128, -0x3.ffe0000000000004p+12)));
    try std.testing.expectEqual(0x3.ffffffffffffffffffffffffa748p-16384, exp2(@as(f128, -0x3.ffe0000000000000000000000002p+12)));
    //  try std.testing.expectEqual(0x3.ffffffffffffffffffffffd3a37cp-16384, exp2(@as(f128, -0x3.ffe00000000000000000000001p+12)));
    try std.testing.expectEqual(0x4.00b18178ba33b141b486ff22689p-16384, exp2(@as(f128, -0x3.ffdffcp+12)));
    try std.testing.expectEqual(0x4p-16384, exp2(@as(f128, -0x3.ffep+12)));
    try std.testing.expectEqual(0x4.00000000058b90bfbe9253c51e4p-16384, exp2(@as(f128, -0x3.ffdfffffffffep+12)));
    try std.testing.expectEqual(0x4.000000000000b17217f7d1cf890cp-16384, exp2(@as(f128, -0x3.ffdffffffffffffcp+12)));
    try std.testing.expectEqual(0x4.00000000000000000000000058b8p-16384, exp2(@as(f128, -0x3.ffdffffffffffffffffffffffffep+12)));
    try std.testing.expectEqual(0x4.00000000000000000000002c5c84p-16384, exp2(@as(f128, -0x3.ffdfffffffffffffffffffffffp+12)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, exp2(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1.a44722ff862d74360c22ab25d2cdp+0, exp2(@as(f128, 0xb.71754p-4)));
    try std.testing.expectEqual(0x3.959e67fd7ff858c3dda97946a1a2p+12, exp2(@as(f128, 0xd.d77dp+0)));
    try std.testing.expectEqual(0x1.afdd736c287aa8000406087bccf5p+0, exp2(@as(f128, 0xc.122c4p-4)));
    try std.testing.expectEqual(0x6.546d5ccd21bad0545e3ae48d3b3p-4, exp2(@as(f128, -0x1.567cc8p+0)));
    try std.testing.expectEqual(0x4.cfe0085ef004d24004a566c1b274p-4, exp2(@as(f128, -0x1.bbbd76p+0)));
    try std.testing.expectEqual(0xd.3ce16388003d339d8e42c2ed7088p-308, exp2(@as(f128, -0x1.3045fep+8)));
    try std.testing.expectEqual(0x5.c6bfd7fd625f811d85ee0f45e71p+8, exp2(@as(f128, 0xa.87b8bp+0)));
    try std.testing.expectEqual(0x8.a8744fff686ede7e5204943f8a98p-4, exp2(@as(f128, -0xe.2ce69p-4)));
    try std.testing.expectEqual(0xf.ff79b6bee6bd7ffc6db60f67e948p-4, exp2(@as(f128, -0xc.1bf12p-16)));
    try std.testing.expectEqual(0xd.23271e170997ffff8d5111790ddp-4, exp2(@as(f128, -0x4.8ce878p-4)));
    try std.testing.expectEqual(0x1.f6b64a6870e6ae124dad946fb894p+0, exp2(@as(f128, 0xf.93d19p-4)));
    try std.testing.expectEqual(0x1.f6b6490bfcd17676f008c989d53ap+0, exp2(@as(f128, 0xf.93d18p-4)));
    try std.testing.expectEqual(0x1.f6b64a10a015e99701a69e715b1fp+0, exp2(@as(f128, 0xf.93d18bf7be8d8p-4)));
    try std.testing.expectEqual(0x1.f6b64a10a015deb360fb026c4e21p+0, exp2(@as(f128, 0xf.93d18bf7be8dp-4)));
    try std.testing.expectEqual(0x1.f6b64a10a015e20774d776dcd953p+0, exp2(@as(f128, 0xf.93d18bf7be8d272p-4)));
}
