const std = @import("std");
const types = @import("../types.zig");
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
