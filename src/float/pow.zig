const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const log2_data = @import("log2_data.zig");
const exp2_data = @import("exp2_data.zig");
const pow_data = @import("pow_data.zig");
const exp_data = @import("exp_data.zig");
const ldbl128 = @import("ldbl128.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn pow(x: anytype, y: anytype) EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.pow: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    comptime if (types.numericType(@TypeOf(y)) != .int and types.numericType(@TypeOf(y)) != .float)
        @compileError("float.pow: y must be an int or float, got " ++ @typeName(@TypeOf(y)));

    switch (EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y)))) {
        f16 => return scast(f16, pow32(scast(f32, x), scast(f32, y))),
        f32 => {
            // glibc/sysdeps/ieee754/flt-32/e_powf.c
            return pow32(scast(f32, x), scast(f32, y));
        },
        f64 => {
            // glibc/sysdeps/ieee754/dbl-64/e_pow.c
            return pow64(scast(f64, x), scast(f64, y));
        },
        f80 => return scast(f80, pow128(scast(f128, x), scast(f128, y))),
        f128 => {
            // glibc/sysdeps/ieee754/ldbl-128/e_powl.c
            return pow128(scast(f128, x), scast(f128, y));
        },
        else => unreachable,
    }
}

// Subnormal input is normalized so ix has negative biased exponent.
// Output is multiplied by N (POWF_SCALE) if TOINT_INTRINICS is set.
inline fn log2_inline32(ix: u32) f64 {
    // x = 2^k z; where z is in range [0x3f330000,2*0x3f330000] and exact.
    // The range is split into N subintervals.
    // The ith subinterval contains z and c is near its center.
    const tmp: u32 = ix -% 0x3f330000;
    const i: i32 = scast(i32, @mod(tmp >> 19, 16));
    const top: u32 = tmp & 0xff800000;
    const iz: u32 = ix -% top;
    var k: i32 = undefined;
    {
        @setRuntimeSafety(false);
        k = @as(i32, @intCast(top)) >> 23; // arithmetic shift
    }
    const invc: f64 = log2_data.tab_32[@intCast(i)].invc;
    const logc: f64 = log2_data.tab_32[@intCast(i)].logc;
    const z: f64 = scast(f64, @as(f32, @bitCast(iz)));

    // log2(x) = log1p(z/c-1)/ln2 + log2(c) + k
    const r: f64 = z * invc - 1;
    const y0: f64 = logc + scast(f64, k);

    // Pipelined polynomial evaluation to approximate log1p(r)/ln2.
    const r2: f64 = r * r;
    var y: f64 = log2_data.poly_pow_32[0] * r + log2_data.poly_pow_32[1];
    const p: f64 = log2_data.poly_pow_32[2] * r + log2_data.poly_pow_32[3];
    const r4: f64 = r2 * r2;
    var q: f64 = log2_data.poly_pow_32[4] * r + y0;
    q = p * r2 + q;
    y = y * r4 + q;
    return y;
}

// The output of log2 and thus the input of exp2 is either scaled by N
// (in case of fast toint intrinsics) or not.  The unscaled xd must be
// in [-1021,1023], sign_bias sets the sign of the result.
inline fn exp2_inline32(xd: f64, sign_bias: u32) f64 {
    // x = k/N + r with r in [-1/(64), 1/(64)]
    var kd: f64 = xd + exp2_data.shift_scaled_32;
    const ki: u64 = @bitCast(kd);
    kd -= exp2_data.shift_scaled_32; // k/32
    const r: f64 = xd - kd;

    // exp2(x) = 2^(k/32) * 2^r ~= s * (C0*r^3 + C1*r^2 + C2*r + 1)
    var t: u64 = exp2_data.T_32[ki % 32];
    const ski: u64 = ki + scast(u64, sign_bias);
    t +%= ski << 47;
    const s: f64 = @bitCast(t);
    const z: f64 = exp2_data.poly_32[0] * r + exp2_data.poly_32[1];
    const r2: f64 = r * r;
    var y: f64 = exp2_data.poly_32[2] * r + 1;
    y = z * r2 + y;
    y = y * s;
    return y;
}

// Returns 0 if not int, 1 if odd int, 2 if even int.  The argument is
// the bit representation of a non-zero finite floating-point value.
inline fn checkint32(iy: u32) i32 {
    const e: i32 = scast(i32, iy >> 23 & 0xff);
    if (e < 0x7f)
        return 0;

    if (e > 0x7f + 23)
        return 2;

    if ((iy & ((@as(u32, 1) << @as(u5, @intCast(0x7f + 23 - e))) - 1)) != 0)
        return 0;

    if ((iy & (@as(u32, 1) << @as(u5, @intCast(0x7f + 23 - e)))) != 0)
        return 1;

    return 2;
}

inline fn zeroinfnan32(ix: u32) bool {
    return 2 *% ix -% 1 >= 2 * 0x7f800000 -% 1;
}

fn pow32(x: f32, y: f32) f32 {
    var sign_bias: u32 = 0;
    var ix: u32 = @bitCast(x);
    const iy: u32 = @bitCast(y);
    if (ix -% 0x00800000 >= 0x7f800000 -% 0x00800000 or zeroinfnan32(iy)) {
        @branchHint(.unlikely);
        // Either (x < 0x1p-126 or inf or nan) or (y is 0 or inf or nan).
        if (zeroinfnan32(iy)) {
            @branchHint(.unlikely);
            if (2 *% iy == 0)
                return if (std.math.isSignalNan(x)) x + y else 1;

            if (ix == 0x3f800000)
                return if (std.math.isSignalNan(y)) x + y else 1;

            if (2 *% ix > 2 *% 0x7f800000 or 2 *% iy > 2 * 0x7f800000)
                return x + y;

            if (2 *% ix == 2 *% 0x3f800000)
                return 1;

            if ((2 *% ix < 2 *% 0x3f800000) == (iy & 0x80000000 == 0))
                return 0; // |x|<1 && y==inf or |x|>1 && y==-inf.

            return y * y;
        }
        if (zeroinfnan32(ix)) {
            @branchHint(.unlikely);
            var x2: f32 = x * x;
            if (ix & 0x80000000 != 0 and checkint32(iy) == 1) {
                x2 = -x2;
                sign_bias = 1;
            }

            if (2 *% ix == 0 and iy & 0x80000000 != 0)
                return @as(f32, (if (sign_bias != 0) -1 else 1)) / @as(f32, 0);

            return if (iy & 0x80000000 != 0) 1 / x2 else x2;
        }
        // x and y are non-zero finite.
        if ((ix & 0x80000000) != 0) {
            // Finite x < 0.
            const yint: i32 = checkint32(iy);
            if (yint == 0)
                return (x - x) / (x - x);

            if (yint == 1)
                sign_bias = 65536;

            ix &= 0x7fffffff;
        }

        if (ix < 0x00800000) {
            // Normalize subnormal x so exponent becomes negative.
            ix = @bitCast(x * 0x1p23);
            ix &= 0x7fffffff;
            ix -%= 23 << 23;
        }
    }

    const logx: f64 = log2_inline32(ix);
    const ylogx: f64 = scast(f64, y) * logx; // Note: cannot overflow, y is single prec.
    if ((@as(u64, @bitCast(ylogx)) >> 47 & 0xffff) >= @as(u64, @bitCast(@as(f64, 126))) >> 47) {
        @branchHint(.unlikely);
        // |y*log(x)| >= 126.
        if (ylogx > 0x1.fffffffd1d571p+6) {
            // |x^y| > 0x1.ffffffp127.
            return @as(f32, (if (sign_bias != 0) -0x1p97 else 0x1p97)) * 0x1p97;
        }

        if (ylogx > 0x1.fffffffa3aae2p+6) {
            // |x^y| > 0x1.fffffep127, check if we round away from 0.
            if ((sign_bias == 0 and 1 + 0x1p-25 != 1) or (sign_bias != 0 and -1 - 0x1p-25 != -1))
                return @as(f32, (if (sign_bias != 0) -0x1p97 else 0x1p97)) * 0x1p97;
        }

        if (ylogx <= -150)
            return @as(f32, (if (sign_bias != 0) -0x1p-95 else 0x1p-95)) * 0x1p-95;

        if (ylogx < -149)
            return @as(f32, (if (sign_bias != 0) -0x1.4p-75 else 0x1.4p-75)) * 0x1.4p-75;
    }
    return scast(f32, exp2_inline32(ylogx, sign_bias));
}

// Top 12 bits of a double (sign and exponent bits).
inline fn top12(x: f64) u32 {
    return @truncate(@as(u64, @bitCast(x)) >> 52);
}

// Compute y+TAIL = log(x) where the rounded result is y and TAIL has about
// additional 15 bits precision.  IX is the bit representation of x, but
// normalized in the subnormal range using the sign bit for the exponent.
inline fn log_inline64(ix: u64, tail: *f64) f64 {
    // x = 2^k z; where z is in range [0x3fe6955500000000,2*0x3fe6955500000000) and exact.
    // The range is split into N subintervals.
    // The ith subinterval contains z and c is near its center.
    const tmp: u64 = ix -% 0x3fe6955500000000;
    const i: i32 = scast(i32, (tmp >> 45) % 128);
    var k: i32 = undefined;
    {
        @setRuntimeSafety(false);
        k = scast(i32, @as(i64, @intCast(tmp)) >> 52); // arithmetic shift
    }
    const iz: u64 = ix -% (tmp & 0xfff << 52);
    const z: f64 = @bitCast(iz);
    const kd: f64 = scast(f64, k);

    // log(x) = k*Ln2 + log(c) + log1p(z/c-1).
    const invc: f64 = pow_data.tab_64[@intCast(i)].invc;
    const logc: f64 = pow_data.tab_64[@intCast(i)].logc;
    const logctail: f64 = pow_data.tab_64[@intCast(i)].logctail;

    // Note: 1/c is j/N or j/N/2 where j is an integer in [N,2N) and
    // |z/c - 1| < 1/N, so r = z/c - 1 is exactly representible.
    var r: f64 = undefined;
    var rhi: f64 = undefined;
    var rlo: f64 = undefined;
    if (true) {
        r = @mulAdd(f64, z, invc, -1);
    } else {
        // Split z such that rhi, rlo and rhi*rhi are exact and |rlo| <= |r|.
        const zhi: f64 = @bitCast((iz + (1 << 31)) & (-1 << 32));
        const zlo: f64 = z - zhi;
        rhi = zhi * invc - 1.0;
        rlo = zlo * invc;
        r = rhi + rlo;
    }

    // k*Ln2 + log(c) + r.
    const t1: f64 = kd * pow_data.ln2hi_64 + logc;
    const t2: f64 = t1 + r;
    const lo1: f64 = kd * pow_data.ln2lo_64 + logctail;
    const lo2: f64 = t1 - t2 + r;

    // Evaluation is optimized assuming superscalar pipelined execution.
    const ar: f64 = pow_data.poly_64[0] * r; // pow_data.poly_64[0] = -0.5.
    const ar2: f64 = r * ar;
    const ar3: f64 = r * ar2;
    // k*Ln2 + log(c) + r + pow_data.poly_64[0]*r*r.
    var hi: f64 = undefined;
    var lo3: f64 = undefined;
    var lo4: f64 = undefined;
    if (true) {
        hi = t2 + ar2;
        lo3 = @mulAdd(f64, ar, r, -ar2);
        lo4 = t2 - hi + ar2;
    } else {
        const arhi: f64 = pow_data.poly_64[0] * rhi;
        const arhi2: f64 = rhi * arhi;
        hi = t2 + arhi2;
        lo3 = rlo * (ar + arhi);
        lo4 = t2 - hi + arhi2;
    }

    // p = log1p(r) - r - pow_data.poly_64[0]*r*r.
    const p: f64 = (ar3 * (pow_data.poly_64[1] + r * pow_data.poly_64[2] + ar2 * (pow_data.poly_64[3] + r * pow_data.poly_64[4] + ar2 * (pow_data.poly_64[5] + r * pow_data.poly_64[6]))));
    const lo: f64 = lo1 + lo2 + lo3 + lo4 + p;
    const y: f64 = hi + lo;
    tail.* = hi - y + lo;
    return y;
}

// Handle cases that may overflow or underflow when computing the result that
// is scale*(1+TMP) without intermediate rounding.  The bit representation of
// scale is in SBITS, however it has a computed exponent that may have
// overflown into the sign bit so that needs to be adjusted before using it as
// a double.  (int32_t)KI is the k used in the argument reduction and exponent
// adjustment of scale, positive k here means the result may overflow and
// negative k means the result may underflow.  */
inline fn specialcase(tmp: f64, sbits: u64, ki: u64) f64 {
    var ssbits: u64 = sbits;
    if ((ki & 0x80000000) == 0) {
        // k > 0, the exponent of scale might have overflowed by <= 460.
        ssbits -%= 1009 << 52;
        const scale: f64 = @bitCast(ssbits);
        const y: f64 = 0x1p1009 * (scale + scale * tmp);
        return y;
    }

    // k < 0, need special care in the subnormal range.
    ssbits += 1022 << 52;
    // Note: ssbits is signed scale.
    const scale: f64 = @bitCast(ssbits);
    var y: f64 = scale + scale * tmp;
    if (float.abs(y) < 1) {
        // Round y to the right precision before scaling it into the subnormal
        // range to avoid double rounding that can cause 0.5+E/2 ulp error where
        // E is the worst-case ulp error outside the subnormal range.  So this
        // is only useful if the goal is better than 1 ulp worst-case error.
        var one: f64 = 1;
        if (y < 0)
            one = -1;

        var lo: f64 = scale - y + scale * tmp;
        const hi: f64 = one + y;
        lo = one - hi + y + lo;
        y = hi + lo - one;

        // Fix the sign of 0.
        if (y == 0.0)
            y = @bitCast(ssbits & 0x8000000000000000);

        // The underflow exception needs to be signaled explicitly.
        std.mem.doNotOptimizeAway(0x1p-1022 * 0x1p-1022);
    }

    y = 0x1p-1022 * y;
    return y;
}

// Computes sign*exp(x+xtail) where |xtail| < 2^-8/N and |xtail| <= |x|.
// The sign_bias argument is SIGN_BIAS or 0 and sets the sign to -1 or 1.
inline fn exp_inline64(x: f64, xtail: f64, sign_bias: u32) f64 {
    var abstop: u32 = top12(x) & 0x7ff;
    if (abstop -% top12(0x1p-54) >= top12(512) -% top12(0x1p-54)) {
        @branchHint(.unlikely);
        if (abstop -% top12(0x1p-54) >= 0x80000000) {
            // Avoid spurious underflow for tiny x.
            // Note: 0 is common input.
            const one: f64 = 1 + x;
            return if (sign_bias != 0) -one else one;
        }
        if (abstop >= top12(1024.0)) {
            // Note: inf and nan are already handled.
            if ((@as(u64, @bitCast(x)) >> 63) != 0) {
                return @as(f64, (if (sign_bias != 0) -0x1p-767 else 0x1p-767)) * 0x1p-767;
            } else {
                return @as(f64, (if (sign_bias != 0) -0x1p769 else 0x1p769)) * 0x1p769;
            }
        }
        // Large x is special cased below.
        abstop = 0;
    }

    // exp(x) = 2^(k/N) * exp(r), with exp(r) in [2^(-1/2N),2^(1/2N)].
    // x = ln2/N*k + r, with int k and r in [-ln2/2N, ln2/2N].
    const z: f64 = exp_data.InvLn2N_64 * x;
    // z - kd is in [-1, 1] in non-nearest rounding modes.
    var kd: f64 = z + exp_data.Shift_64;
    const ki: u64 = @bitCast(kd);
    kd -= exp_data.Shift_64;
    var r: f64 = x + kd * exp_data.NegLn2hiN_64 + kd * exp_data.NegLn2loN_64;
    // The code assumes 2^-200 < |xtail| < 2^-8/N.
    r += xtail;
    // 2^(k/N) ~= scale * (1 + tail).
    const idx: u64 = 2 * (ki % 128);
    const top: u64 = (ki + scast(u64, sign_bias)) << 45;
    const tail: f64 = @bitCast(exp_data.T_64[idx]);
    // This is only a valid scale when -1023*N < k < 1024*N.
    const sbits: u64 = exp_data.T_64[idx + 1] +% top;
    // exp(x) = 2^(k/N) * exp(r) ~= scale + scale * (tail + exp(r) - 1).
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
    // Note: tmp == 0 or |tmp| > 2^-200 and scale > 2^-739, so there
    // is no spurious underflow here even without fma.
    return scale + scale * tmp;
}

// Returns 0 if not int, 1 if odd int, 2 if even int.  The argument is
// the bit representation of a non-zero finite floating-point value.
inline fn checkint64(iy: u64) i32 {
    const e: i32 = scast(i32, iy >> 52 & 0x7ff);
    if (e < 0x3ff)
        return 0;

    if (e > 0x3ff + 52)
        return 2;

    if ((iy & ((@as(u64, 1) << @as(u6, @intCast(0x3ff + 52 - e))) - 1)) != 0)
        return 0;

    if ((iy & (@as(u64, 1) << @as(u6, @intCast(0x3ff + 52 - e)))) != 0)
        return 1;

    return 2;
}

// Returns 1 if input is the bit representation of 0, infinity or nan.
inline fn zeroinfnan64(i: u64) bool {
    return 2 *% i -% 1 >= 2 *% @as(u64, @bitCast(std.math.inf(f64))) -% 1;
}

fn pow64(x: f64, y: f64) f64 {
    var sign_bias: u32 = 0;
    var ix: u64 = @bitCast(x);
    const iy: u64 = @bitCast(y);
    var topx: u32 = top12(x);
    const topy: u32 = top12(y);
    if (topx -% 0x001 >= 0x7ff - 0x001 or (topy & 0x7ff) -% 0x3be >= 0x43e -% 0x3be) {
        @branchHint(.unlikely);
        // Note: if |y| > 1075 * ln2 * 2^53 ~= 0x1.749p62 then pow(x,y) = inf/0
        // and if |y| < 2^-54 / 1075 ~= 0x1.e7b6p-65 then pow(x,y) = +-1.
        // Special cases: (x < 0x1p-126 or inf or nan) or
        // (|y| < 0x1p-65 or |y| >= 0x1p63 or nan).
        if (zeroinfnan64(iy)) {
            @branchHint(.unlikely);
            if (2 *% iy == 0)
                return if (std.math.isSignalNan(x)) x + y else 1;

            if (ix == @as(u64, @bitCast(@as(f64, 1))))
                return if (std.math.isSignalNan(y)) x + y else 1;

            if (2 *% ix > 2 *% @as(u64, @bitCast(std.math.inf(f64))) or 2 *% iy > 2 *% @as(u64, @bitCast(std.math.inf(f64))))
                return x + y;
            if (2 *% ix == 2 *% @as(u64, @bitCast(@as(f64, 1))))
                return 1;
            if ((2 *% ix < 2 *% @as(u64, @bitCast(@as(f64, 1)))) == ((iy >> 63) == 0))
                return 0; // |x|<1 && y==inf or |x|>1 && y==-inf.

            return y * y;
        }

        if (zeroinfnan64(ix)) {
            @branchHint(.unlikely);
            var x2: f64 = x * x;
            if ((ix >> 63) != 0 and checkint64(iy) == 1) {
                x2 = -x2;
                sign_bias = 1;
            }
            if (2 *% ix == 0 and (iy >> 63) != 0)
                return @as(f64, (if (sign_bias != 0) -1 else 1)) / @as(f64, 0);

            return if ((iy >> 63) != 0) 1 / x2 else x2;
        }

        // Here x and y are non-zero finite.
        if ((ix >> 63) != 0) {
            // Finite x < 0.
            const yint: i32 = checkint64(iy);
            if (yint == 0)
                return (x - x) / (x - x);

            if (yint == 1)
                sign_bias = (0x800 << 7);

            ix &= 0x7fffffffffffffff;
            topx &= 0x7ff;
        }

        if ((topy & 0x7ff) -% 0x3be >= 0x43e - 0x3be) {
            // Note: sign_bias == 0 here because y is not odd.
            if (ix == @as(u64, @bitCast(@as(f64, 1))))
                return 1;

            if ((topy & 0x7ff) < 0x3be) {
                // |y| < 2^-65, x^y ~= 1 + y*log(x).
                return if (ix > @as(u64, @bitCast(@as(f64, 1)))) 1 + y else 1 - y;
            }

            return if ((ix > @as(u64, @bitCast(@as(f64, 1)))) == (topy < 0x800)) 0x1p769 * 0x1p769 else 0x1p-767 * 0x1p-767;
        }
        if (topx == 0) {
            // Normalize subnormal x so exponent becomes negative.
            ix = @bitCast(x * 0x1p52);
            ix &= 0x7fffffffffffffff;
            ix -%= 52 << 52;
        }
    }

    var lo: f64 = undefined;
    const hi: f64 = log_inline64(ix, &lo);
    var ehi: f64 = undefined;
    var elo: f64 = undefined;
    if (true) {
        ehi = y * hi;
        elo = y * lo + @mulAdd(f64, y, hi, -ehi);
    } else {
        const yhi: f64 = @bitCast(iy & -1 << 27);
        const ylo: f64 = y - yhi;
        const lhi: f64 = @bitCast(@as(u64, @bitCast(hi)) & -1 << 27);
        const llo: f64 = hi - lhi + lo;
        ehi = yhi * lhi;
        elo = ylo * lhi + y * llo; // |elo| < |ehi| * 2^-25.
    }

    return exp_inline64(ehi, elo, sign_bias);
}

fn pow128(x: f128, y: f128) f128 {
    const bp: [2]f128 = .{
        1,
        1.5,
    };
    // log_2(1.5)
    const dp_h: [2]f128 = .{
        0,
        5.8496250072115607565592654282227158546448e-1,
    };
    // Low part of log_2(1.5)
    const dp_l: [2]f128 = .{
        0,
        1.0579781240112554492329533686862998106046e-16,
    };
    const two113: f128 = 1.0384593717069655257060992658440192e34;
    const huge: f128 = 1.0e3000;
    const tiny: f128 = 1.0e-3000;
    // 3/2 log x = 3 z + z^3 + z^3 (z^2 R(z^2))
    // z = (x-1)/(x+1)
    // 1 <= x <= 1.25
    // Peak relative error 2.3e-37
    const LN: [5]f128 = .{
        -3.0779177200290054398792536829702930623200e1,
        6.5135778082209159921251824580292116201640e1,
        -4.6312921812152436921591152809994014413540e1,
        1.2510208195629420304615674658258363295208e1,
        -9.9266909031921425609179910128531667336670e-1,
    };
    const LD: [5]f128 = .{
        -5.129862866715009066465422805058933131960e1,
        1.452015077564081884387441590064272782044e2,
        -1.524043275549860505277434040464085593165e2,
        7.236063513651544224319663428634139768808e1,
        -1.494198912340228235853027849917095580053e1,
    };
    // exp(x) = 1 + x - x / (1 - 2 / (x - x^2 R(x^2)))
    // 0 <= x <= 0.5
    // Peak relative error 5.7e-38
    const PN: [5]f128 = .{
        5.081801691915377692446852383385968225675e8,
        9.360895299872484512023336636427675327355e6,
        4.213701282274196030811629773097579432957e4,
        5.201006511142748908655720086041570288182e1,
        9.088368420359444263703202925095675982530e-3,
    };
    const PD: [4]f128 = .{
        3.049081015149226615468111430031590411682e9,
        1.069833887183886839966085436512368982758e8,
        8.259257717868875207333991924545445705394e5,
        1.872583833284143212651746812884298360922e3,
    };
    // ln 2
    const lg2: f128 = 6.9314718055994530941723212145817656807550e-1;
    const lg2_h: f128 = 6.9314718055994528622676398299518041312695e-1;
    const lg2_l: f128 = 2.3190468138462996154948554638754786504121e-17;
    const ovt: f128 = 8.0085662595372944372e-0017;
    // 2/(3*log(2))
    const cp: f128 = 9.6179669392597560490661645400126142495110e-1;
    const cp_h: f128 = 9.6179669392597555432899980587535537779331e-1;
    const cp_l: f128 = 5.0577616648125906047157785230014751039424e-17;

    const p: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const hx: i32 = @bitCast(p.w0);
    var ix: u32 = @bitCast(hx & 0x7fffffff);

    const q: ldbl128.ieee_f128_shape32 = @bitCast(y);
    const hy: i32 = @bitCast(q.w0);
    const iy: u32 = @bitCast(hy & 0x7fffffff);

    // y==0: x**0 = 1
    if ((iy | q.w1 | q.w2 | q.w3) == 0 and !std.math.isSignalNan(x))
        return 1;

    // 1.0**y = 1; -1.0**+-Inf = 1
    if (x == 1 and !std.math.isSignalNan(y))
        return 1;

    if (x == -1 and iy == 0x7fff0000 and (q.w1 | q.w2 | q.w3) == 0)
        return 1;

    // +-NaN return x+y
    if ((ix > 0x7fff0000) or ((ix == 0x7fff0000) and ((p.w1 | p.w2 | p.w3) != 0)) or (iy > 0x7fff0000) or ((iy == 0x7fff0000) and ((q.w1 | q.w2 | q.w3) != 0)))
        return x + y;

    // determine if y is an odd int when x < 0
    // yisint = 0       ... y is not an integer
    // yisint = 1       ... y is an odd int
    // yisint = 2       ... y is an even int
    var yisint: i32 = 0;
    if (hx < 0) {
        if (iy >= 0x40700000) { // 2^113
            yisint = 2; // even integer y
        } else if (iy >= 0x3fff0000) { // 1.0
            if (float.floor(y) == y) {
                const z: f128 = 0.5 * y;
                if (float.floor(z) == z) {
                    yisint = 2;
                } else {
                    yisint = 1;
                }
            }
        }
    }

    // special value of y
    if ((q.w1 | q.w2 | q.w3) == 0) {
        if (iy == 0x7fff0000) { // y is +-inf
            if (((ix - 0x3fff0000) | p.w1 | p.w2 | p.w3) == 0) {
                return y - y; // +-1**inf is NaN
            } else if (ix >= 0x3fff0000) { // (|x|>1)**+-inf = inf,0
                return if (hy >= 0) y else 0;
            } else { // (|x|<1)**-,+inf = inf,0
                return if (hy < 0) -y else 0;
            }
        }

        if (iy == 0x3fff0000) { // y is  +-1
            if (hy < 0) {
                return 1 / x;
            } else {
                return x;
            }
        }

        if (hy == 0x40000000)
            return x * x; // y is  2

        if (hy == 0x3ffe0000) { // y is  0.5
            if (hx >= 0) // x >= +0
                return float.sqrt(x);
        }
    }

    var ax: f128 = float.abs(x);
    // special value of x
    if ((p.w1 | p.w2 | p.w3) == 0) {
        if (ix == 0x7fff0000 or ix == 0 or ix == 0x3fff0000) {
            var z: f128 = ax; // x is +-0,+-inf,+-1
            if (hy < 0)
                z = 1 / z; // z = (1/|x|)

            if (hx < 0) {
                if (((ix -% 0x3fff0000) | @as(u32, @bitCast(yisint))) == 0) {
                    z = (z - z) / (z - z); // (-1)**non-int is NaN
                } else if (yisint == 1) {
                    z = -z; // (x<0)**odd = -(|x|**odd)
                }
            }
            return z;
        }
    }

    // (x<0)**(non-int) is NaN
    {
        @setRuntimeSafety(false);
        if ((((@as(u32, @intCast(hx)) >> 31) -% 1) | @as(u32, @bitCast(yisint))) == 0)
            return (x - x) / (x - x);
    }

    // sgn (sign of result -ve**odd) = -1 else = 1
    var sgn: f128 = 1;
    {
        @setRuntimeSafety(false);
        if ((((@as(u32, @intCast(hx)) >> 31) -% 1) | @as(u32, @bitCast(yisint -% 1))) == 0)
            sgn = -1; // (-ve)**(odd int)
    }

    // |y| is huge.
    // 2^-16495 = 1/2 of smallest representable value.
    // If (1 - 1/131072)^y underflows, y > 1.4986e9
    if (iy > 0x401d654b) {
        // if (1 - 2^-113)^y underflows, y > 1.1873e38
        if (iy > 0x407d654b) {
            if (ix <= 0x3ffeffff)
                return if (hy < 0) huge * huge else tiny * tiny;

            if (ix >= 0x3fff0000)
                return if (hy > 0) huge * huge else tiny * tiny;
        }
        // over/underflow if x is not close to 1
        if (ix < 0x3ffeffff)
            return if (hy < 0) sgn * huge * huge else sgn * tiny * tiny;

        if (ix > 0x3fff0000)
            return if (hy > 0) sgn * huge * huge else sgn * tiny * tiny;
    }

    const ay: f128 = if (y > 0) y else -y;
    var yy: f128 = y;
    if (ay < 0x1p-128)
        yy = if (y < 0) -0x1p-128 else 0x1p-128;

    var n: i32 = 0;
    // take care subnormal number
    if (ix < 0x00010000) {
        ax *= two113;
        n -= 113;
        const o: ldbl128.ieee_f128_shape32 = @bitCast(ax);
        ix = o.w0;
    }

    {
        @setRuntimeSafety(false);
        n +%= @intCast(((ix) >> 16) -% 0x3fff);
    }
    var j: i32 = @bitCast(ix & 0x0000ffff);
    // determine interval
    ix = @bitCast(j | 0x3fff0000); // normalize ix
    var k: i32 = undefined;
    if (j <= 0x3988) {
        k = 0; // |x|<sqrt(3/2)
    } else if (j < 0xbb67) {
        k = 1; // |x|<sqrt(3)
    } else {
        k = 0;
        n += 1;
        ix -= 0x00010000;
    }

    var o: ldbl128.ieee_f128_shape32 = @bitCast(ax);
    o.w0 = ix;
    ax = @bitCast(o);

    // compute s = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5)
    var u: f128 = ax - bp[@intCast(k)]; // bp[0]=1.0, bp[1]=1.5
    var v: f128 = 1 / (ax + bp[@intCast(k)]);
    const s: f128 = u * v;
    var s_h: f128 = s;

    o = @bitCast(s_h);
    o.w3 = 0;
    o.w2 &= 0xf8000000;
    s_h = @bitCast(o);
    // t_h=ax+bp[k] High
    var t_h: f128 = ax + bp[@intCast(k)];
    o = @bitCast(t_h);
    o.w3 = 0;
    o.w2 &= 0xf8000000;
    t_h = @bitCast(o);
    var t_l: f128 = ax - (t_h - bp[@intCast(k)]);
    const s_l: f128 = v * ((u - s_h * t_h) - s_h * t_l);
    // compute log(ax)
    var s2: f128 = s * s;
    u = LN[0] + s2 * (LN[1] + s2 * (LN[2] + s2 * (LN[3] + s2 * LN[4])));
    v = LD[0] + s2 * (LD[1] + s2 * (LD[2] + s2 * (LD[3] + s2 * (LD[4] + s2))));
    var r: f128 = s2 * s2 * u / v;
    r += s_l * (s_h + s);
    s2 = s_h * s_h;
    t_h = 3.0 + s2 + r;
    o = @bitCast(t_h);
    o.w3 = 0;
    o.w2 &= 0xf8000000;
    t_h = @bitCast(o);
    t_l = r - ((t_h - 3.0) - s2);
    // u+v = s*(1+...)
    u = s_h * t_h;
    v = s_l * t_h + t_l * s;
    // 2/(3log2)*(s+...)
    var p_h: f128 = u + v;
    o = @bitCast(p_h);
    o.w3 = 0;
    o.w2 &= 0xf8000000;
    p_h = @bitCast(o);
    var p_l: f128 = v - (p_h - u);
    const z_h: f128 = cp_h * p_h; // cp_h+cp_l = 2/(3*log2)
    const z_l: f128 = cp_l * p_h + p_l * cp + dp_l[@intCast(k)];
    // log2(ax) = (s+..)*2/(3*log2) = n + dp_h + z_h + z_l
    var t: f128 = scast(f128, n);
    var t1: f128 = (((z_h + z_l) + dp_h[@intCast(k)]) + t);
    o = @bitCast(t1);
    o.w3 = 0;
    o.w2 &= 0xf8000000;
    t1 = @bitCast(o);
    const t2: f128 = z_l - (((t1 - t) - dp_h[@intCast(k)]) - z_h);

    // split up y into y1+y2 and compute (y1+y2)*(t1+t2)
    var y1: f128 = yy;
    o = @bitCast(y1);
    o.w3 = 0;
    o.w2 &= 0xf8000000;
    y1 = @bitCast(o);
    p_l = (yy - y1) * t1 + yy * t2;
    p_h = y1 * t1;
    var z: f128 = p_l + p_h;
    o = @bitCast(z);
    j = @bitCast(o.w0);
    if (j >= 0x400d0000) { // z >= 16384
        // if z > 16384
        if ((@as(u32, @bitCast(j - 0x400d0000)) | o.w1 | o.w2 | o.w3) != 0) {
            return sgn * huge * huge; // overflow
        } else {
            if (p_l + ovt > z - p_h)
                return sgn * huge * huge; // overflow
        }
    } else if ((j & 0x7fffffff) >= 0x400d01b9) { // z <= -16495
        // z < -16495
        if (((@as(u32, @bitCast(j)) - 0xc00d01bc) | o.w1 | o.w2 | o.w3) != 0) {
            return sgn * tiny * tiny; // underflow
        } else {
            if (p_l <= z - p_h)
                return sgn * tiny * tiny; // underflow
        }
    }

    // compute 2**(p_h+p_l)
    const i: i32 = j & 0x7fffffff;
    k = (i >> 16) - 0x3fff;
    n = 0;
    if (i > 0x3ffe0000) { // if |z| > 0.5, set n = [z+0.5]
        n = scast(i32, float.floor(z + 0.5));
        t = scast(f128, n);
        p_h -= t;
    }
    t = p_l + p_h;
    o = @bitCast(t);
    o.w3 = 0;
    o.w2 &= 0xf8000000;
    t = @bitCast(o);
    u = t * lg2_h;
    v = (p_l - (t - p_h)) * lg2 + t * lg2_l;
    z = u + v;
    const w: f128 = v - (z - u);
    // exp(z)
    t = z * z;
    u = PN[0] + t * (PN[1] + t * (PN[2] + t * (PN[3] + t * PN[4])));
    v = PD[0] + t * (PD[1] + t * (PD[2] + t * (PD[3] + t)));
    t1 = z - t * u / v;
    r = (z * t1) / (t1 - 2) - (w + z * w);
    z = 1 - (r - z);
    o = @bitCast(z);
    j = @bitCast(o.w0);
    j += (n << 16);
    if ((j >> 16) <= 0) {
        z = float.scalbn(z, n); // subnormal output
        const force_underflow: f128 = z * z;
        std.mem.doNotOptimizeAway(force_underflow);
    } else {
        o.w0 = @bitCast(j);
        z = @bitCast(o);
    }

    return sgn * z;
}
