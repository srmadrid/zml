const cast = @import("../types.zig").cast;
const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const log_data = @import("log_data.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;

pub fn log(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return log(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, log32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_logf.c
                    return log32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_log.c
                    return log64(x);
                },
                f80 => return cast(f80, log128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_logl.c
                    return log128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn log32(x: f32) f32 {
    var ix: u32 = @bitCast(x);

    if (ix == 0x3f800000) {
        @branchHint(.unlikely);
        return 0;
    }

    if (ix -% 0x00800000 >= 0x7f800000 - 0x00800000) {
        @branchHint(.unlikely);
        // x < 0x1p-126 or inf or nan.
        if (ix * 2 == 0)
            return std.math.inf(f32);

        if (ix == 0x7f800000) // log(inf) == inf.
            return x;

        if ((ix & 0x80000000) != 0 or ix * 2 >= 0xff000000)
            return std.math.nan(f32);

        // x is subnormal, normalize it.
        ix = @bitCast(cast(f32, x, .{}) * 0x1p23);
        ix -%= 23 << 23;
    }

    // x = 2^k z; where z is in range [0x3f330000,2*0x3f330000] and exact.
    // The range is split into 16 subintervals.
    // The ith subinterval contains z and c is near its center.
    const tmp: u32 = ix -% 0x3f330000;
    const i: i32 = @bitCast((tmp >> 19) % 16);
    const k: i32 = @as(i32, @bitCast(tmp)) >> 23; // arithmetic shift
    const iz: u32 = ix -% (tmp & (0x1ff << 23));
    const invc: f64 = log_data.T_32[@intCast(i)][0];
    const logc: f64 = log_data.T_32[@intCast(i)][1];
    const z: f64 = cast(f64, @as(f32, @bitCast(iz)), .{});

    // log(x) = log1p(z/c-1) + log(c) + k*Ln2
    const r: f64 = z * invc - 1;
    const y0: f64 = logc + cast(f64, k, .{}) * log_data.Ln2_32;

    // Pipelined polynomial evaluation to approximate log1p(r).
    const r2: f64 = r * r;
    var y: f64 = log_data.A_32[1] * r + log_data.A_32[2];
    y = log_data.A_32[0] * r2 + y;
    y = y * r2 + (y0 + r);
    return cast(f32, y, .{});
}

fn log64(x: f64) f64 {
    var ix: u64 = @bitCast(x);
    const top: u32 = @intCast(@as(u64, @bitCast(x)) >> 48);

    const LO: u64 = @bitCast(1.0 - @as(f64, 0x1p-4));
    const HI: u64 = @bitCast(1.0 + @as(f64, 0x1.09p-4));
    if (ix -% LO < HI - LO) {
        @branchHint(.unlikely);
        // Handle close to 1.0 inputs separately.
        if (ix == 1)
            return 0;

        const r: f64 = x - 1;
        const r2: f64 = r * r;
        const r3: f64 = r * r2;
        var y: f64 = r3 * (log_data.B_64[1] + r * log_data.B_64[2] + r2 * log_data.B_64[3] + r3 * (log_data.B_64[4] + r * log_data.B_64[5] + r2 * log_data.B_64[6] + r3 * (log_data.B_64[7] + r * log_data.B_64[8] + r2 * log_data.B_64[9] + r3 * log_data.B_64[10])));
        // Worst-case error is around 0.507 ULP.
        var w: f64 = r * 0x1p27;
        const rhi: f64 = r + w - w;
        const rlo: f64 = r - rhi;
        w = rhi * rhi * log_data.B_64[0]; // log_data.B_64[0] == -0.5.
        const hi: f64 = r + w;
        var lo: f64 = r - hi + w;
        lo += log_data.B_64[0] * rlo * (rhi + r);
        y += lo;
        y += hi;
        return y;
    }

    if (top -% 0x0010 >= 0x7ff0 - 0x0010) {
        @branchHint(.unlikely);
        // x < 0x1p-1022 or inf or nan.
        if (ix * 2 == 0)
            return std.math.inf(f64);

        if (ix == @as(u64, @bitCast(std.math.inf(f64)))) // log(inf) == inf.
            return x;

        if ((top & 0x8000) != 0 or (top & 0x7ff0) == 0x7ff0)
            return std.math.nan(f64);

        // x is subnormal, normalize it.
        ix = @bitCast(x * 0x1p52);
        ix -%= 52 << 52;
    }

    // x = 2^k z; where z is in range [0x3fe6000000000000,2*0x3fe6000000000000) and exact.
    // The range is split into 128 subintervals.
    // The ith subinterval contains z and c is near its center.
    const tmp: u64 = ix -% 0x3fe6000000000000;
    const i: u64 = (tmp >> 45) % 128;
    const k: i32 = @intCast(@as(i64, @bitCast(tmp)) >> 52); // arithmetic shift
    const iz: u64 = ix -% (tmp & (0xfff << 52));
    const invc: f64 = log_data.T_64[i][0];
    const logc: f64 = log_data.T_64[i][1];
    const z: f64 = @bitCast(iz);

    // log(x) = log1p(z/c-1) + log(c) + k*Ln2.
    // r ~= z/c - 1, |r| < 1/(2*128).
    var r: f64 = undefined;
    if (true) { // if fma
        // rounding error: 0x1p-55/128.
        r = @mulAdd(f64, z, invc, -1);
    } else {
        // rounding error: 0x1p-55/128 + 0x1p-66.
        r = (z - log_data.T2_64[i][0] - log_data.T2_64[i][1]) * invc;
    }
    const kd: f64 = cast(f64, k, .{});

    // hi + lo = r + log(c) + k*Ln2.
    const w: f64 = kd * log_data.Ln2hi_64 + logc;
    const hi: f64 = w + r;
    const lo: f64 = w - hi + r + kd * log_data.Ln2lo_64;

    // log(x) = lo + (log1p(r) - r) + hi.
    const r2: f64 = r * r; // rounding error: 0x1p-54/128^2.
    // Worst case error if |y| > 0x1p-4: 0.519 ULP (0.520 ULP without fma).
    // 0.5 + 2.06/128 + abs-poly-error*2^56 ULP (+ 0.001 ULP without fma).
    const y: f64 = lo + r2 * log_data.A_64[0] + r * r2 * (log_data.A_64[1] + r * log_data.A_64[2] + r2 * (log_data.A_64[3] + r * log_data.A_64[4])) + hi;
    return y;
}

fn log128(x: f128) f128 {
    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    var m: u32 = u.w0;

    // Check for IEEE special cases.
    var k: i32 = cast(i32, m & 0x7fffffff, .{});
    // log(0) = -infinity.
    if ((@as(u32, @bitCast(k)) | u.w1 | u.w2 | u.w3) == 0)
        return -0.5 / @as(f128, 0);

    // log ( x < 0 ) = NaN
    if ((m & 0x80000000) != 0)
        return (x - x) / 0;

    // log (infinity or NaN)
    if (k >= 0x7fff0000)
        return x + x;

    // Extract exponent and reduce domain to 0.703125 <= u < 1.40625
    var e: i32 = undefined;
    u = @bitCast(math.frexp(x, &e));
    m = u.w0 & 0xffff;
    m |= 0x10000;
    // Find lookup table index k from high order bits of the significand.
    var t: ldbl128.ieee_f128_shape32 = undefined;
    if (m < 0x16800) {
        k = cast(i32, (m - 0xff00) >> 9, .{});
        // t is the argument 0.5 + (k+26)/128
        // of the nearest item to u in the lookup table.
        t.w0 = @bitCast(0x3fff0000 + (k << 9));
        t.w1 = 0;
        t.w2 = 0;
        t.w3 = 0;
        u.w0 += 0x10000;
        e -= 1;
        k += 64;
    } else {
        k = cast(i32, (m - 0xfe00) >> 10, .{});
        t.w0 = @bitCast(0x3ffe0000 + (k << 10));
        t.w1 = 0;
        t.w2 = 0;
        t.w3 = 0;
    }

    // On this interval the table is not used due to cancellation error.
    var z: f128 = undefined;
    if ((x <= 1.0078125) and (x >= 0.9921875)) {
        if (x == 1)
            return 0;

        z = x - 1;
        k = 64;
        t = @bitCast(@as(f128, 1));
        e = 0;
    } else {
        // log(u) = log( t u/t ) = log(t) + log(u/t)
        // log(t) is tabulated in the lookup table.
        // Express log(u/t) = log(1+z),  where z = u/t - 1 = (u-t)/t.
        // cf. Cody & Waite.
        z = (@as(f128, @bitCast(u)) - @as(f128, @bitCast(t))) / @as(f128, @bitCast(t));
    }

    // Series expansion of log(1+z).
    const w: f128 = z * z;
    var y: f128 = ((((((((((((log_data.l15_128 * z + log_data.l14_128) * z + log_data.l13_128) * z + log_data.l12_128) * z + log_data.l11_128) * z + log_data.l10_128) * z + log_data.l9_128) * z + log_data.l8_128) * z + log_data.l7_128) * z + log_data.l6_128) * z + log_data.l5_128) * z + log_data.l4_128) * z + log_data.l3_128) * z * w;
    y -= 0.5 * w;
    y += cast(f128, e, .{}) * log_data.ln2b_128; // Base 2 exponent offset times ln(2).
    y += z;
    y += log_data.logtbl_128[@intCast(k - 26)]; // log(t) - (t-1)
    y += (@as(f128, @bitCast(t)) - 1);
    y += cast(f128, e, .{}) * log_data.ln2a_128;
    return y;
}

test log {
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1p+0, log(@as(f32, 0x2.b7e154p+0)));
    try std.testing.expectEqual(0xf.fffffp-4, log(@as(f32, 0x2.b7e15p+0)));
    try std.testing.expectEqual(-0x1p+0, log(@as(f32, 0x5.e2d59p-4)));
    try std.testing.expectEqual(-0x1p+0, log(@as(f32, 0x5.e2d588p-4)));
    try std.testing.expectEqual(0xb.17218p-4, log(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(0x2.4d7638p+0, log(@as(f32, 0xap+0)));
    try std.testing.expectEqual(-0x4.9a5888p-4, log(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x1.fffffep-24, log(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1.fffffep-24, log(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.fffffep-24, log(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.fffffep-24, log(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.fffffep-24, log(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x1p-24, log(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x1p-24, log(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x1p-24, log(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x1p-24, log(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x1p-24, log(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x5.75628p+4, log(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x6.74768p+4, log(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x5.8b90cp+4, log(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x5.eb58f8p-4, log(@as(f32, 0xb.0d5dfp-4)));
    try std.testing.expectEqual(0x5.a47afp-4, log(@as(f32, 0x1.6c3f6p+0)));
    try std.testing.expectEqual(-0x6.772d38p-4, log(@as(f32, 0xa.ae688p-4)));
    try std.testing.expectEqual(0x1.e811a8p+4, log(@as(f32, 0x1.017f8ap+44)));
    try std.testing.expectEqual(0x1.8ff28cp+4, log(@as(f32, 0x1.0b5c1ep+36)));
    try std.testing.expectEqual(0x5.8b90cp+4, log(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x7.3eb07p-4, log(@as(f32, 0x1.929d9cp+0)));
    try std.testing.expectEqual(0x6.1ba94p-4, log(@as(f32, 0x1.770072p+0)));
    try std.testing.expectEqual(-0x1.6fe0dp-4, log(@as(f32, 0xe.a0289p-4)));
    try std.testing.expectEqual(-0x1.6fe0e2p-4, log(@as(f32, 0xe.a0288p-4)));

    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.000000f647926p+0, log(@as(f64, 0x2.b7e154p+0)));
    try std.testing.expectEqual(0xf.fffff7d922f5p-4, log(@as(f64, 0x2.b7e15p+0)));
    try std.testing.expectEqual(0x1p+0, log(@as(f64, 0x2.b7e151628aed4p+0)));
    try std.testing.expectEqual(0x1p+0, log(@as(f64, 0x2.b7e151628aed2p+0)));
    try std.testing.expectEqual(-0xf.fffff952d5f5p-4, log(@as(f64, 0x5.e2d59p-4)));
    try std.testing.expectEqual(-0x1.000000f11e086p+0, log(@as(f64, 0x5.e2d588p-4)));
    try std.testing.expectEqual(-0x1p+0, log(@as(f64, 0x5.e2d58d8b3bcep-4)));
    try std.testing.expectEqual(-0x1.0000000000001p+0, log(@as(f64, 0x5.e2d58d8b3bcdcp-4)));
    try std.testing.expectEqual(0xb.17217f7d1cf78p-4, log(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(0x2.4d763776aaa2cp+0, log(@as(f64, 0xap+0)));
    try std.testing.expectEqual(-0x4.9a58844d36e48p-4, log(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x1.fffffe000002bp-24, log(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1.fffffe000002bp-24, log(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-56, log(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffe000002bp-24, log(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-56, log(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffe000002bp-24, log(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-56, log(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffe000002bp-24, log(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-56, log(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(-0x1.0000008000005p-24, log(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005p-24, log(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8p-56, log(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005p-24, log(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8p-56, log(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005p-24, log(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8p-56, log(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005p-24, log(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8p-56, log(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x5.75627cbf9441cp+4, log(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x2.c4657baf579a4p+8, log(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x2.9fa8dcb9092a6p+8, log(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x6.74767f33d1dcp+4, log(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x2.e870a88dae386p+8, log(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bcp+4, log(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c5c85fdf473dep+8, log(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x5.eb58f885a32d8p-4, log(@as(f64, 0xb.0d5dfp-4)));
    try std.testing.expectEqual(0x5.a47aee2b5c35p-4, log(@as(f64, 0x1.6c3f6p+0)));
    try std.testing.expectEqual(-0x6.772d36f0dd28cp-4, log(@as(f64, 0xa.ae688p-4)));
    try std.testing.expectEqual(0x1.e811a8a66aa57p+4, log(@as(f64, 0x1.017f8ap+44)));
    try std.testing.expectEqual(0x1.8ff28cfed79a1p+4, log(@as(f64, 0x1.0b5c1ep+36)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bcp+4, log(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.83fc3c37fc59p+8, log(@as(f64, 0x2.1b17c2887e938p+928)));
    try std.testing.expectEqual(0x7.3eb06c60714c4p-4, log(@as(f64, 0x1.929d9cp+0)));
    try std.testing.expectEqual(0x6.1ba943bb20434p-4, log(@as(f64, 0x1.770072p+0)));
    try std.testing.expectEqual(-0x1.6fe0d0a6311e3p-4, log(@as(f64, 0xe.a0289p-4)));
    try std.testing.expectEqual(-0x1.6fe0e22718ad7p-4, log(@as(f64, 0xe.a0288p-4)));
    // try std.testing.expectEqual(-0x1.6fe0d4c400979p-4, log(@as(f64, 0xe.a0288c3cb5ecp-4)));

    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.000000f647925f34p+0, log(@as(f80, 0x2.b7e154p+0)));
    try std.testing.expectEqual(0xf.fffff7d922f51a3p-4, log(@as(f80, 0x2.b7e15p+0)));
    try std.testing.expectEqual(0x1.00000000000007fp+0, log(@as(f80, 0x2.b7e151628aed4p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffc2bp-4, log(@as(f80, 0x2.b7e151628aed2p+0)));
    try std.testing.expectEqual(0x1p+0, log(@as(f80, 0x2.b7e151628aed2a6cp+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffp-4, log(@as(f80, 0x2.b7e151628aed2a68p+0)));
    try std.testing.expectEqual(-0xf.fffff952d5f52b9p-4, log(@as(f80, 0x5.e2d59p-4)));
    try std.testing.expectEqual(-0x1.000000f11e085f42p+0, log(@as(f80, 0x5.e2d588p-4)));
    try std.testing.expectEqual(-0xf.ffffffffffffd91p-4, log(@as(f80, 0x5.e2d58d8b3bcep-4)));
    try std.testing.expectEqual(-0x1.000000000000087p+0, log(@as(f80, 0x5.e2d58d8b3bcdcp-4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffp-4, log(@as(f80, 0x5.e2d58d8b3bcdf1bp-4)));
    try std.testing.expectEqual(-0x1p+0, log(@as(f80, 0x5.e2d58d8b3bcdf1a8p-4)));
    try std.testing.expectEqual(0xb.17217f7d1cf79acp-4, log(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(0x2.4d763776aaa2b05cp+0, log(@as(f80, 0xap+0)));
    try std.testing.expectEqual(-0x4.9a58844d36e49e1p-4, log(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x1.fffffe000002aaaap-24, log(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1.fffffe000002aaaap-24, log(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-56, log(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffe000002aaaap-24, log(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-56, log(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffffffffffffep-64, log(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x1.fffffe000002aaaap-24, log(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-56, log(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffffffffffffep-64, log(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x1.fffffe000002aaaap-24, log(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-56, log(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffffffffffffep-64, log(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(-0x1.0000008000005556p-24, log(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005556p-24, log(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.0000000000002p-56, log(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005556p-24, log(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.0000000000002p-56, log(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1p-64, log(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005556p-24, log(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.0000000000002p-56, log(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1p-64, log(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005556p-24, log(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.0000000000002p-56, log(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1p-64, log(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x5.75627cbf9441de28p+4, log(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x2.c4657baf579a47bcp+8, log(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x2.c5b2319c4843accp+12, log(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x2.c5bd48bdc7c0c9b8p+12, log(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x2.9fa8dcb9092a538cp+8, log(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x6.74767f33d1dc1d1p+4, log(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x2.e870a88dae386c74p+8, log(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x2.c86ce2daa80dcdbp+12, log(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bc56p+4, log(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6a8p+8, log(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6bp+12, log(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x5.eb58f885a32d894p-4, log(@as(f80, 0xb.0d5dfp-4)));
    try std.testing.expectEqual(0x5.a47aee2b5c34f8p-4, log(@as(f80, 0x1.6c3f6p+0)));
    try std.testing.expectEqual(-0x6.772d36f0dd28c27p-4, log(@as(f80, 0xa.ae688p-4)));
    try std.testing.expectEqual(0x1.e811a8a66aa56988p+4, log(@as(f80, 0x1.017f8ap+44)));
    try std.testing.expectEqual(0x1.8ff28cfed79a1002p+4, log(@as(f80, 0x1.0b5c1ep+36)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bc56p+4, log(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.83fc3c37fc59p+8, log(@as(f80, 0x2.1b17c2887e938p+928)));
    try std.testing.expectEqual(0x7.3eb06c60714c5ff8p-4, log(@as(f80, 0x1.929d9cp+0)));
    try std.testing.expectEqual(0x6.1ba943bb20434dc8p-4, log(@as(f80, 0x1.770072p+0)));
    try std.testing.expectEqual(-0x1.6fe0d0a6311e31f2p-4, log(@as(f80, 0xe.a0289p-4)));
    try std.testing.expectEqual(-0x1.6fe0e22718ad7752p-4, log(@as(f80, 0xe.a0288p-4)));
    try std.testing.expectEqual(-0x1.6fe0d4c40097884ep-4, log(@as(f80, 0xe.a0288c3cb5ecp-4)));

    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.000000f647925f34d03716a8b6ccp+0, log(@as(f128, 0x2.b7e154p+0)));
    try std.testing.expectEqual(0xf.fffff7d922f51a2d208d1c4e821p-4, log(@as(f128, 0x2.b7e15p+0)));
    try std.testing.expectEqual(0x1.00000000000007f0a06e4ddb0222p+0, log(@as(f128, 0x2.b7e151628aed4p+0)));
    // try std.testing.expectEqual(0xf.ffffffffffffc2af553376366578p-4, log(@as(f128, 0x2.b7e151628aed2p+0)));
    try std.testing.expectEqual(0x1.000000000000000075ed29d49ac4p+0, log(@as(f128, 0x2.b7e151628aed2a6cp+0)));
    // try std.testing.expectEqual(0xf.ffffffffffffffefd37c671cbd08p-4, log(@as(f128, 0x2.b7e151628aed2a68p+0)));
    try std.testing.expectEqual(0x1p+0, log(@as(f128, 0x2.b7e151628aed2a6abf7158809cf6p+0)));
    // try std.testing.expectEqual(0xf.fffffffffffffffffffffffffff8p-4, log(@as(f128, 0x2.b7e151628aed2a6abf7158809cf4p+0)));
    try std.testing.expectEqual(0x1.0000000000000000000000000004p+0, log(@as(f128, 0x2.b7e151628aed2a6abf7158809dp+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffa6p-4, log(@as(f128, 0x2.b7e151628aed2a6abf7158809cp+0)));
    try std.testing.expectEqual(-0xf.fffff952d5f52b972627765c7b8p-4, log(@as(f128, 0x5.e2d59p-4)));
    try std.testing.expectEqual(-0x1.000000f11e085f422347d5acdb97p+0, log(@as(f128, 0x5.e2d588p-4)));
    // try std.testing.expectEqual(-0xf.ffffffffffffd90c7882a506a588p-4, log(@as(f128, 0x5.e2d58d8b3bcep-4)));
    // try std.testing.expectEqual(-0x1.00000000000008704ccdb47c1f23p+0, log(@as(f128, 0x5.e2d58d8b3bcdcp-4)));
    // try std.testing.expectEqual(-0xf.fffffffffffffff4415f776b07c8p-4, log(@as(f128, 0x5.e2d58d8b3bcdf1bp-4)));
    try std.testing.expectEqual(-0x1.0000000000000000a006a027f5f3p+0, log(@as(f128, 0x5.e2d58d8b3bcdf1a8p-4)));
    try std.testing.expectEqual(-0x1p+0, log(@as(f128, 0x5.e2d58d8b3bcdf1abadec7829055p-4)));
    // try std.testing.expectEqual(-0x1.0000000000000000000000000001p+0, log(@as(f128, 0x5.e2d58d8b3bcdf1abadec7829054cp-4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffffe2p-4, log(@as(f128, 0x5.e2d58d8b3bcdf1abadec782906p-4)));
    try std.testing.expectEqual(-0x1.0000000000000000000000000039p+0, log(@as(f128, 0x5.e2d58d8b3bcdf1abadec782904p-4)));
    try std.testing.expectEqual(0xb.17217f7d1cf79abc9e3b39803f3p-4, log(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(0x2.4d763776aaa2b05ba95b58ae0b4cp+0, log(@as(f128, 0xap+0)));
    try std.testing.expectEqual(-0x4.9a58844d36e49e0efadd9db02aa8p-4, log(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x1.fffffe000002aaaaa6aaaab11111p-24, log(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1.fffffe000002aaaaa6aaaab11111p-24, log(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8000000000000558p-56, log(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffe000002aaaaa6aaaab11111p-24, log(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8000000000000558p-56, log(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffffffffffffep-64, log(@as(f128, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x1.fffffe000002aaaaa6aaaab11111p-24, log(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8000000000000558p-56, log(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffffffffffffep-64, log(@as(f128, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffffep-108, log(@as(f128, 0x1.000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x1.fffffe000002aaaaa6aaaab11111p-24, log(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xf.ffffffffffff8000000000000558p-56, log(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.fffffffffffffffep-64, log(@as(f128, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffff8p-116, log(@as(f128, 0x1.0000000000000000000000000001p+0)));
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffffep-108, log(@as(f128, 0x1.000000000000000000000000008p+0)));
    try std.testing.expectEqual(-0x1.0000008000005555559555558889p-24, log(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005555559555558889p-24, log(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.00000000000020000000000000a8p-56, log(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005555559555558889p-24, log(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.00000000000020000000000000a8p-56, log(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.00000000000000008p-64, log(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005555559555558889p-24, log(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.00000000000020000000000000a8p-56, log(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.00000000000000008p-64, log(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x4.000000000000000000000000008p-108, log(@as(f128, 0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(0x0p+0, log(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.0000008000005555559555558889p-24, log(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.00000000000020000000000000a8p-56, log(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.00000000000000008p-64, log(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x8p-116, log(@as(f128, 0xf.fffffffffffffffffffffffffff8p-4)));
    try std.testing.expectEqual(-0x4.000000000000000000000000008p-108, log(@as(f128, 0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(-0x5.75627cbf9441de28d5e1264d1f18p+4, log(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x2.c4657baf579a47bbcffb06f8dfc4p+8, log(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x2.c5b2319c4843acbff21591e99cccp+12, log(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x2.c5bd48bdc7c0c9b78cd23024d64cp+12, log(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x2.9fa8dcb9092a538b3f2ee2ca66f2p+8, log(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x6.74767f33d1dc1d0fc8187877a4c8p+4, log(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x2.e870a88dae386c72b4fd4773c092p+8, log(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x2.c86ce2daa80dcdaf0680827cc35ap+12, log(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x2.c877f9fc278aeaa6a13d20b7fcdcp+12, log(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x2.ca8c50440f005913a49acbd2c4e8p+12, log(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bc55e4f18476ac644p+4, log(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6a7278ece600fccp+8, log(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6af277ece600fccp+12, log(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6af278ece600fccp+12, log(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6ab278ece600fccp+8, log(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x5.eb58f885a32d893cb272dbe106p-4, log(@as(f128, 0xb.0d5dfp-4)));
    try std.testing.expectEqual(0x5.a47aee2b5c34f7fed8c38bb73dc4p-4, log(@as(f128, 0x1.6c3f6p+0)));
    try std.testing.expectEqual(-0x6.772d36f0dd28c26cc42127335304p-4, log(@as(f128, 0xa.ae688p-4)));
    try std.testing.expectEqual(0x1.e811a8a66aa569880c5e8ea2ec2p+4, log(@as(f128, 0x1.017f8ap+44)));
    try std.testing.expectEqual(0x1.8ff28cfed79a1001419ce243f3acp+4, log(@as(f128, 0x1.0b5c1ep+36)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bc55e4f18476ac644p+4, log(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.83fc3c37fc58ffff8f99749ff1e6p+8, log(@as(f128, 0x2.1b17c2887e938p+928)));
    // try std.testing.expectEqual(0x7.3eb06c60714c5ffbcdb915af2664p-4, log(@as(f128, 0x1.929d9cp+0)));
    try std.testing.expectEqual(0x6.1ba943bb20434dc4abd932bca664p-4, log(@as(f128, 0x1.770072p+0)));
    try std.testing.expectEqual(-0x1.6fe0d0a6311e31f19855212ae415p-4, log(@as(f128, 0xe.a0289p-4)));
    try std.testing.expectEqual(-0x1.6fe0e22718ad77516665df92f5dbp-4, log(@as(f128, 0xe.a0288p-4)));
    try std.testing.expectEqual(-0x1.6fe0d4c40097884d86068c297d0cp-4, log(@as(f128, 0xe.a0288c3cb5ecp-4)));
}
