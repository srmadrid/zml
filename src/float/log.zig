const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const log_data = @import("log_data.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn log(x: anytype) EnsureFloat(@TypeOf(x)) {
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
            return 1 / x;

        if (ix == 0x7f800000) // log(inf) == inf.
            return x;

        if ((ix & 0x80000000) != 0 or ix * 2 >= 0xff000000)
            return (x - x) / (x - x);

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
    u = @bitCast(float.frexp(x, &e));
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
