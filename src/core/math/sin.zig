const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const sincos_data = @import("sincos_data.zig");
const atnat = @import("atnat.zig");
const usncs = @import("usncs.zig");
const cos = @import("cos.zig");
const branred = @import("branred.zig");
const ldbl128 = @import("ldbl128.zig");
const rem_pio2 = @import("rem_pio2.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn sin(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return sin(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, sin32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_sinf.c
                    return sin32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_sin.c
                    return sin64(x);
                },
                f80 => return cast(f80, sin128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_sinl.c
                    return sin128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

// Fast sinf implementation.  Worst-case ULP is 0.5607, maximum relative
// error is 0.5303 * 2^-23.  A single-step range reduction is used for
// small values.  Large inputs have their range reduced using fast integer
// arithmetic.
fn sin32(y: f32) f32 {
    var x: f64 = cast(f64, y, .{});
    var p: *const sincos_data.Sincos32 = &sincos_data.__sincos32_table[0];

    if (sincos_data.abstop12(y) < sincos_data.abstop12(sincos_data.pio4)) {
        const s: f64 = x * x;

        if (sincos_data.abstop12(y) < sincos_data.abstop12(0x1p-12)) {
            @branchHint(.unlikely);
            // Force underflow for tiny y.
            if (sincos_data.abstop12(y) < sincos_data.abstop12(0x1p-126)) {
                @branchHint(.unlikely);
                std.mem.doNotOptimizeAway(cast(f32, s, .{}));
            }

            return y;
        }

        return sincos_data.sin32_poly(x, s, p, 0);
    } else if (sincos_data.abstop12(y) < sincos_data.abstop12(120.0)) {
        @branchHint(.likely);
        var n: i32 = undefined;
        x = sincos_data.reduce_fast(x, p, &n);

        // Setup the signs for sin and cos.
        const s: f64 = p.sign[@intCast(n & 3)];

        if ((n & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        return sincos_data.sin32_poly(x * s, x * x, p, n);
    } else if (sincos_data.abstop12(y) < sincos_data.abstop12(std.math.inf(f32))) {
        const xi: u32 = @bitCast(y);
        const sign: i32 = @intCast(xi >> 31);

        var n: i32 = undefined;
        x = sincos_data.reduce_large(xi, &n);

        // Setup signs for sin and cos - include original sign.
        const s: f64 = p.sign[@intCast((n + sign) & 3)];

        if (((n + sign) & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        return sincos_data.sin32_poly(x * s, x * x, p, n);
    } else {
        return (y - y) / (y - y);
    }
}

// Given a number partitioned into X and DX, this function computes the sine of
// the number by combining the sin and cos of X (as computed by a variation of
// the Taylor series) with the values looked up from the sin/cos table to get
// the result.
pub inline fn do_sin64(x: f64, dx: f64) f64 {
    // Max ULP is 0.501 if |x| < 0.126, otherwise ULP is 0.518.
    if (math.abs(x) < 0.126) {
        const t: f64 = ((((((usncs.s5 * x * x + usncs.s4) * x * x + usncs.s3) * x * x + usncs.s2) * x * x) + usncs.s1) * x - 0.5 * dx) * x * x + dx;
        return x + t;
    }

    var dxnew: f64 = dx;
    if (x <= 0)
        dxnew = -dxnew;

    const u: [2]i32 = @bitCast(usncs.big + math.abs(x));
    var xnew: f64 = x;
    xnew = math.abs(xnew) - (@as(f64, @bitCast(u)) - usncs.big);

    const xx: f64 = xnew * xnew;
    const s: f64 = xnew + (dxnew + xnew * xx * (sincos_data.sn3 + xx * sincos_data.sn5));
    const c: f64 = xnew * dxnew + xx * (sincos_data.cs2 + xx * (sincos_data.cs4 + xx * sincos_data.cs6));
    const k: i32 = u[atnat.LOW_HALF] << 2;
    const sn: f64 = @bitCast(sincos_data.__sincostab[@intCast(k)]);
    const ssn: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 1)]);
    const cs: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 2)]);
    const ccs: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 3)]);
    const cor: f64 = (ssn + s * ccs - sn * c) + cs * s;
    return math.copysign(sn + cor, x);
}

//******************************************************************
// An ultimate sin routine. Given an IEEE double machine number x
// it computes the rounded value of sin(x).
//******************************************************************
fn sin64(x: f64) f64 {
    const u: [2]i32 = @bitCast(x);
    const m: i32 = u[atnat.HIGH_HALF];
    const k: i32 = 0x7fffffff & m; // no sign
    if (k < 0x3e500000) { // if x->0 =>sin(x)=x
        if (math.abs(x) < std.math.floatMin(f64)) {
            const vx: f64 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }
        return x;
    }
    //--------------------------- 2^-26<|x|< 0.855469----------------------
    else if (k < 0x3feb6000) {
        // Max ULP is 0.548.
        return do_sin64(x, 0);
    }
    //----------------------- 0.855469  <|x|<2.426265  ----------------------
    else if (k < 0x400368fd) {
        const t: f64 = usncs.hp0 - math.abs(x);
        // Max ULP is 0.51.
        return math.copysign(cos.do_cos64(t, usncs.hp1), x);
    }
    //-------------------------- 2.426265<|x|< 105414350 ----------------------
    else if (k < 0x419921fb) {
        var a: f64 = undefined;
        var da: f64 = undefined;
        const n: i32 = sincos_data.reduce_sincos64(x, &a, &da);
        return sincos_data.do_sincos64(a, da, n);
    }
    //--------------------105414350 <|x| <2^1024------------------------------
    else if (k < 0x7ff00000) {
        var a: f64 = undefined;
        var da: f64 = undefined;
        const n: i32 = branred.branred(x, &a, &da);
        return sincos_data.do_sincos64(a, da, n);
    }
    //--------------------- |x| > 2^1024 ----------------------------------
    else {
        return x / x;
    }
}

pub fn kernel_sin128(x: f128, y: f128, iy: i32) f128 {
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);
    var tix: u32 = @truncate(@as(u64, @bitCast(ix)) >> 32);
    tix &= ~@as(u32, 0x80000000); // tix = |x|'s high 32 bits
    if (tix < 0x3ffc3000) { // |x| < 0.1484375
        // Argument is small enough to approximate it by a Chebyshev
        // polynomial of degree 17.
        if (tix < 0x3fc60000) { // |x| < 2^-57
            if (math.abs(x) < std.math.floatMin(f128)) {
                const vx: f128 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            if (cast(i32, x, .{}) == 0)
                return x; // generate inexact
        }
        const z: f128 = x * x;
        return x + (x * (z * (sincos_data.SIN1 + z * (sincos_data.SIN2 + z * (sincos_data.SIN3 + z * (sincos_data.SIN4 +
            z * (sincos_data.SIN5 + z * (sincos_data.SIN6 + z * (sincos_data.SIN7 + z * sincos_data.SIN8)))))))));
    } else {
        // So that we don't have to use too large polynomial,  we find
        // l and h such that x = l + h,  where fabsl(l) <= 1.0/256 with 83
        // possible values for h.  We look up cosl(h) and sinl(h) in
        // pre-computed tables,  compute cosl(l) and sinl(l) using a
        // Chebyshev polynomial of degree 10(11) and compute
        // sinl(h+l) = sinl(h)cosl(l) + cosl(h)sinl(l).
        var index: u32 = 0x3ffe - (tix >> 16);
        const hix: u32 = (tix + (@as(u32, 0x200) << @as(u5, @intCast(index)))) & (@as(u32, 0xfffffc00) << @as(u5, @intCast(index)));
        const xx: f128 = math.abs(x);
        switch (index) {
            0 => index = ((45 << 10) + hix - 0x3ffe0000) >> 8,
            1 => index = ((13 << 11) + hix - 0x3ffd0000) >> 9,
            else => index = (hix - 0x3ffc3000) >> 10,
        }

        var h: f128 = undefined;
        ldbl128.setWords(&h, cast(u64, hix, .{}) << 32, @as(u64, 0));
        var l: f128 = undefined;
        if (iy != 0) {
            l = (if (ix < 0) -y else y) - (h - xx);
        } else {
            l = xx - h;
        }

        var z: f128 = l * l;
        const sin_l: f128 = l * (1 + z * (sincos_data.SSIN1 + z * (sincos_data.SSIN2 + z * (sincos_data.SSIN3 + z * (sincos_data.SSIN4 + z * sincos_data.SSIN5)))));
        const cos_l_m1: f128 = z * (sincos_data.SCOS1 + z * (sincos_data.SCOS2 + z * (sincos_data.SCOS3 + z * (sincos_data.SCOS4 + z * sincos_data.SCOS5))));
        z = usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_LO] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] * cos_l_m1) + (usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] * sin_l));
        return if (ix < 0) -z else z;
    }
}

fn sin128(x: f128) f128 {
    // High word of x.
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);

    // |x| ~< pi/4
    ix &= 0x7fffffffffffffff;
    if (ix <= 0x3ffe921fb54442d1) {
        return kernel_sin128(x, 0, 0);
    } else if (ix >= 0x7fff000000000000) { // sin(Inf or NaN) is NaN
        return x - x;
    } else { // argument reduction needed
        var y: [2]f128 = undefined;
        const n: i32 = rem_pio2.rem_pio2_128(x, &y);
        switch (n & 3) {
            0 => return kernel_sin128(y[0], y[1], 1),
            1 => return cos.kernel_cos128(y[0], y[1]),
            2 => return -kernel_sin128(y[0], y[1], 1),
            else => return -cos.kernel_cos128(y[0], y[1]),
        }
    }
}

test sin {
    try std.testing.expectEqual(0x0p+0, sin(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sin(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x8p-4, sin(@as(f32, 0x8.60a92p-4)));
    try std.testing.expectEqual(0x7.fffff8p-4, sin(@as(f32, 0x8.60a91p-4)));
    try std.testing.expectEqual(-0x7.fffff8p-4, sin(@as(f32, -0x8.60a91p-4)));
    try std.testing.expectEqual(-0x8p-4, sin(@as(f32, -0x8.60a92p-4)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f32, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f32, 0x1.921fb4p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f32, -0x1.921fb6p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f32, -0x1.921fb4p+0)));
    try std.testing.expectEqual(-0x1.777a5cp-24, sin(@as(f32, 0x3.243f6cp+0)));
    try std.testing.expectEqual(0x2.8885a4p-24, sin(@as(f32, 0x3.243f68p+0)));
    try std.testing.expectEqual(0x1.777a5cp-24, sin(@as(f32, -0x3.243f6cp+0)));
    try std.testing.expectEqual(-0x2.8885a4p-24, sin(@as(f32, -0x3.243f68p+0)));
    try std.testing.expectEqual(0xa.e7fe1p-4, sin(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(-0xc.143e1p-8, sin(@as(f32, 0x2p+64)));
    try std.testing.expectEqual(0xc.143e1p-8, sin(@as(f32, -0x2p+64)));
    try std.testing.expectEqual(-0x1.1e7cfap-24, sin(@as(f32, 0xb.fa09ap+100)));
    try std.testing.expectEqual(0xb.7fb6p-4, sin(@as(f32, 0xc.d4967p-4)));
    try std.testing.expectEqual(0xb.7fb5fp-4, sin(@as(f32, 0xc.d4966p-4)));
    try std.testing.expectEqual(0x3.fe478p-4, sin(@as(f32, 0x4.093388p-4)));
    try std.testing.expectEqual(0x3.fe4778p-4, sin(@as(f32, 0x4.09338p-4)));
    try std.testing.expectEqual(-0x4.cd7e88p-4, sin(@as(f32, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(-0xb.becc4p-4, sin(@as(f32, 0x2.1e19ep+72)));
    try std.testing.expectEqual(-0x8.599b3p-4, sin(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x8.599b3p-4, sin(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x6.0b8d18p-4, sin(@as(f32, 0x1p+120)));
    try std.testing.expectEqual(0x9.f9631p-4, sin(@as(f32, 0x8p+124)));
    try std.testing.expectEqual(0xc.6fa5cp-8, sin(@as(f32, 0xf.ffffcp+124)));
    try std.testing.expectEqual(-0x8.599b3p-4, sin(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x7.f13d78p-4, sin(@as(f32, 0x4p+48)));
    try std.testing.expectEqual(-0xf.c777cp-4, sin(@as(f32, 0x1p+28)));
    try std.testing.expectEqual(0xc.dbc1ap-4, sin(@as(f32, 0xe.ef3bp-4)));
    try std.testing.expectEqual(0xc.dbc19p-4, sin(@as(f32, 0xe.ef3afp-4)));
    try std.testing.expectEqual(0xb.93256p-4, sin(@as(f32, 0x2.553538p+0)));
    try std.testing.expectEqual(0xb.93258p-4, sin(@as(f32, 0x2.553534p+0)));
    try std.testing.expectEqual(-0x9.10bb4p-4, sin(@as(f32, 0x3.be736p+0)));
    try std.testing.expectEqual(-0x9.10bb1p-4, sin(@as(f32, 0x3.be735cp+0)));
    try std.testing.expectEqual(-0xb.4352p-4, sin(@as(f32, 0x3.ec2a04p+0)));
    try std.testing.expectEqual(-0xb.4351dp-4, sin(@as(f32, 0x3.ec2ap+0)));
    try std.testing.expectEqual(-0xc.d263ap-4, sin(@as(f32, 0x4.1237e8p+0)));
    try std.testing.expectEqual(-0xc.d2635p-4, sin(@as(f32, 0x4.1237ep+0)));
    try std.testing.expectEqual(-0xf.f4f47p-4, sin(@as(f32, 0x4.c92d1p+0)));
    try std.testing.expectEqual(-0xf.f4f47p-4, sin(@as(f32, 0x4.c92d08p+0)));
    try std.testing.expectEqual(-0x4.b6f61p-4, sin(@as(f32, 0x5.fbec78p+0)));
    try std.testing.expectEqual(-0x4.b6f688p-4, sin(@as(f32, 0x5.fbec7p+0)));
    try std.testing.expectEqual(0xd.76aa4p-4, sin(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0xe.8c7b7p-4, sin(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(0x2.42070cp-4, sin(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(-0xc.1bdcfp-4, sin(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(-0xf.57c1p-4, sin(@as(f32, 0x5p+0)));
    try std.testing.expectEqual(-0x4.787c6p-4, sin(@as(f32, 0x6p+0)));
    try std.testing.expectEqual(0xa.83046p-4, sin(@as(f32, 0x7p+0)));
    try std.testing.expectEqual(0xf.d4695p-4, sin(@as(f32, 0x8p+0)));
    try std.testing.expectEqual(0x6.98099p-4, sin(@as(f32, 0x9p+0)));
    try std.testing.expectEqual(-0x8.b44f8p-4, sin(@as(f32, 0xap+0)));
    try std.testing.expectEqual(-0x5.595d8p-4, sin(@as(f32, 0x1.200148p+32)));
    try std.testing.expectEqual(0x4.220ffp-4, sin(@as(f32, 0x1.200146p+32)));
    try std.testing.expectEqual(0x8.599b3p-4, sin(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.773a3p-4, sin(@as(f32, 0xe.9f1e6p+112)));
    try std.testing.expectEqual(0x7.76d6p-4, sin(@as(f32, 0xe.9f1e5p+112)));
    try std.testing.expectEqual(-0x1.ffb67ap-4, sin(@as(f32, 0x4.7857dp+68)));
    try std.testing.expectEqual(-0x1.fecbp-4, sin(@as(f32, 0x6.287cdp+0)));
    try std.testing.expectEqual(-0x1.fecb7ep-4, sin(@as(f32, 0x6.287cc8p+0)));
    try std.testing.expectEqual(-0xd.8f692p-4, sin(@as(f32, -0x1.02e34cp+0)));
    try std.testing.expectEqual(-0x8.3beep-4, sin(@as(f32, 0xf.f0274p+4)));
    try std.testing.expectEqual(0x1.ffc6dap-4, sin(@as(f32, 0x3.042d88p+0)));
    // try std.testing.expectEqual(0x1.d12edp-12, sin(@as(f32, 0x1.d12ed2p-12)));
    try std.testing.expectEqual(-0x1.bf207cp-4, sin(@as(f32, -0x6.e23688p+16)));
    try std.testing.expectEqual(-0x2.3e1f7cp-4, sin(@as(f32, -0x6.e2369p+16)));
    try std.testing.expectEqual(-0x1.455068p-4, sin(@as(f32, 0x5.6a5008p+64)));
    try std.testing.expectEqual(-0x1.ee01dcp-4, sin(@as(f32, 0x5.6a5p+64)));
    try std.testing.expectEqual(-0x8.599b3p-4, sin(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x8.599b3p-4, sin(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x4p-128, sin(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-128, sin(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, sin(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, sin(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0xf.fa2aep-4, sin(@as(f32, 0x1.8475e6p+0)));
    try std.testing.expectEqual(0xf.fa2aep-4, sin(@as(f32, 0x1.8475e4p+0)));

    try std.testing.expectEqual(0x0p+0, sin(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sin(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x8.0000036321168p-4, sin(@as(f64, 0x8.60a92p-4)));
    try std.testing.expectEqual(0x7.fffff587e3a04p-4, sin(@as(f64, 0x8.60a91p-4)));
    try std.testing.expectEqual(0x8p-4, sin(@as(f64, 0x8.60a91c16b9b3p-4)));
    try std.testing.expectEqual(0x7.ffffffffffffcp-4, sin(@as(f64, 0x8.60a91c16b9b28p-4)));
    try std.testing.expectEqual(-0x7.fffff587e3a04p-4, sin(@as(f64, -0x8.60a91p-4)));
    try std.testing.expectEqual(-0x8.0000036321168p-4, sin(@as(f64, -0x8.60a92p-4)));
    try std.testing.expectEqual(-0x7.ffffffffffffcp-4, sin(@as(f64, -0x8.60a91c16b9b28p-4)));
    try std.testing.expectEqual(-0x8p-4, sin(@as(f64, -0x8.60a91c16b9b3p-4)));
    try std.testing.expectEqual(0xf.fffffffffffb8p-4, sin(@as(f64, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0xf.fffffffffff3p-4, sin(@as(f64, 0x1.921fb4p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f64, 0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f64, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(-0xf.fffffffffffb8p-4, sin(@as(f64, -0x1.921fb6p+0)));
    try std.testing.expectEqual(-0xf.fffffffffff3p-4, sin(@as(f64, -0x1.921fb4p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f64, -0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f64, -0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(-0x1.777a5cf72cec6p-24, sin(@as(f64, 0x3.243f6cp+0)));
    try std.testing.expectEqual(0x2.8885a308d3106p-24, sin(@as(f64, 0x3.243f68p+0)));
    try std.testing.expectEqual(-0x1.72cece675d1fdp-52, sin(@as(f64, 0x3.243f6a8885a32p+0)));
    try std.testing.expectEqual(0x8.d313198a2e038p-56, sin(@as(f64, 0x3.243f6a8885a3p+0)));
    try std.testing.expectEqual(0x1.777a5cf72cec6p-24, sin(@as(f64, -0x3.243f6cp+0)));
    try std.testing.expectEqual(-0x2.8885a308d3106p-24, sin(@as(f64, -0x3.243f68p+0)));
    try std.testing.expectEqual(0x1.72cece675d1fdp-52, sin(@as(f64, -0x3.243f6a8885a32p+0)));
    try std.testing.expectEqual(-0x8.d313198a2e038p-56, sin(@as(f64, -0x3.243f6a8885a3p+0)));
    try std.testing.expectEqual(0xa.e7fe0b5fc7868p-4, sin(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(-0xc.143e153b0702p-8, sin(@as(f64, 0x2p+64)));
    try std.testing.expectEqual(0xc.143e153b0702p-8, sin(@as(f64, -0x2p+64)));
    try std.testing.expectEqual(-0x1.1e7cf9ec10917p-24, sin(@as(f64, 0xb.fa09ap+100)));
    try std.testing.expectEqual(0xb.7fb6002758778p-4, sin(@as(f64, 0xc.d4967p-4)));
    try std.testing.expectEqual(0xb.7fb5f50739fa8p-4, sin(@as(f64, 0xc.d4966p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776979p-4, sin(@as(f64, 0xc.d4966d92d171p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776979p-4, sin(@as(f64, 0xc.d4966d92d1708p-4)));
    try std.testing.expectEqual(0x3.fe4780403e808p-4, sin(@as(f64, 0x4.093388p-4)));
    // try std.testing.expectEqual(0x3.fe4778810e026p-4, sin(@as(f64, 0x4.09338p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7958p-4, sin(@as(f64, 0x4.093385688a2d4p-4)));
    // try std.testing.expectEqual(0x3.fe477dbdc7954p-4, sin(@as(f64, 0x4.093385688a2dp-4)));
    try std.testing.expectEqual(-0x4.cd7e86c4077cp-4, sin(@as(f64, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(-0xb.becc47ab1b8c8p-4, sin(@as(f64, 0x2.1e19ep+72)));
    try std.testing.expectEqual(-0xd.a29d5bb5f9cb8p-4, sin(@as(f64, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(-0x8.599b32844aba8p-4, sin(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x9.0292465edbbp-4, sin(@as(f64, 0x8p+1020)));
    try std.testing.expectEqual(-0x8.599b32844aba8p-4, sin(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.452fc98b34e97p-8, sin(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x6.0b8d19579bf2cp-4, sin(@as(f64, 0x1p+120)));
    try std.testing.expectEqual(0x9.f963166f215e8p-4, sin(@as(f64, 0x8p+124)));
    try std.testing.expectEqual(0xc.6fa5c5665985p-8, sin(@as(f64, 0xf.ffffcp+124)));
    try std.testing.expectEqual(-0x8.599b32844aba8p-4, sin(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x7.f13d78eabb76cp-4, sin(@as(f64, 0x4p+48)));
    try std.testing.expectEqual(-0xf.c777c6b36a75p-4, sin(@as(f64, 0x1p+28)));
    try std.testing.expectEqual(0xc.dbc19bb4a58a8p-4, sin(@as(f64, 0xe.ef3bp-4)));
    try std.testing.expectEqual(0xc.dbc1922f1d9fp-4, sin(@as(f64, 0xe.ef3afp-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3bp-4, sin(@as(f64, 0xe.ef3af1b5d8008p-4)));
    // try std.testing.expectEqual(0xc.dbc19333ad3a8p-4, sin(@as(f64, 0xe.ef3af1b5d8p-4)));
    try std.testing.expectEqual(0xb.93255854754ap-4, sin(@as(f64, 0x2.553538p+0)));
    try std.testing.expectEqual(0xb.932584840807p-4, sin(@as(f64, 0x2.553534p+0)));
    // try std.testing.expectEqual(0xb.93255eeda1028p-4, sin(@as(f64, 0x2.5535376715bap+0)));
    try std.testing.expectEqual(0xb.93255eeda1038p-4, sin(@as(f64, 0x2.5535376715b9ep+0)));
    try std.testing.expectEqual(-0x9.10bb448d3cbp-4, sin(@as(f64, 0x3.be736p+0)));
    try std.testing.expectEqual(-0x9.10bb0fd0c39d8p-4, sin(@as(f64, 0x3.be735cp+0)));
    try std.testing.expectEqual(-0x9.10bb11242ecap-4, sin(@as(f64, 0x3.be735c19beap+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec8p-4, sin(@as(f64, 0x3.be735c19be9fep+0)));
    try std.testing.expectEqual(-0xb.4351fdda3d818p-4, sin(@as(f64, 0x3.ec2a04p+0)));
    try std.testing.expectEqual(-0xb.4351d06546e7p-4, sin(@as(f64, 0x3.ec2ap+0)));
    try std.testing.expectEqual(-0xb.4351eaad09848p-4, sin(@as(f64, 0x3.ec2a0250032a2p+0)));
    try std.testing.expectEqual(-0xb.4351eaad0983p-4, sin(@as(f64, 0x3.ec2a0250032ap+0)));
    try std.testing.expectEqual(-0xc.d2639f1afc7f8p-4, sin(@as(f64, 0x4.1237e8p+0)));
    try std.testing.expectEqual(-0xc.d2635289f075p-4, sin(@as(f64, 0x4.1237ep+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf59p-4, sin(@as(f64, 0x4.1237e153f7084p+0)));
    // try std.testing.expectEqual(-0xc.d2635f3faf568p-4, sin(@as(f64, 0x4.1237e153f708p+0)));
    try std.testing.expectEqual(-0xf.f4f46a017cb88p-4, sin(@as(f64, 0x4.c92d1p+0)));
    try std.testing.expectEqual(-0xf.f4f4736648dcp-4, sin(@as(f64, 0x4.c92d08p+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f28p-4, sin(@as(f64, 0x4.c92d0ffa4bf04p+0)));
    // try std.testing.expectEqual(-0xf.f4f46a082f288p-4, sin(@as(f64, 0x4.c92d0ffa4bfp+0)));
    try std.testing.expectEqual(-0x4.b6f60ca8d415p-4, sin(@as(f64, 0x5.fbec78p+0)));
    try std.testing.expectEqual(-0x4.b6f686f9ea13p-4, sin(@as(f64, 0x5.fbec7p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a0cp-4, sin(@as(f64, 0x5.fbec7477d4a84p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a48p-4, sin(@as(f64, 0x5.fbec7477d4a8p+0)));
    try std.testing.expectEqual(0xd.76aa47848677p-4, sin(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xe.8c7b7568da23p-4, sin(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(0x2.42070db6daab6p-4, sin(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(-0xc.1bdceeee0f57p-4, sin(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(-0xf.57c0faf04c998p-4, sin(@as(f64, 0x5p+0)));
    try std.testing.expectEqual(-0x4.787c62ac28bp-4, sin(@as(f64, 0x6p+0)));
    try std.testing.expectEqual(0xa.830461368504p-4, sin(@as(f64, 0x7p+0)));
    try std.testing.expectEqual(0xf.d469501467bd8p-4, sin(@as(f64, 0x8p+0)));
    try std.testing.expectEqual(0x6.98098d830be44p-4, sin(@as(f64, 0x9p+0)));
    try std.testing.expectEqual(-0x8.b44f7af9a7a9p-4, sin(@as(f64, 0xap+0)));
    try std.testing.expectEqual(-0x5.595d7e536fe34p-4, sin(@as(f64, 0x1.200148p+32)));
    try std.testing.expectEqual(0x4.220ff25f5cf04p-4, sin(@as(f64, 0x1.200146p+32)));
    try std.testing.expectEqual(-0x6.444fda50019f8p-4, sin(@as(f64, 0x1.2001469775ce6p+32)));
    try std.testing.expectEqual(0x8.599b32844aba8p-4, sin(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.e00885042dd78p-4, sin(@as(f64, -0x3.3de320f6be87ep+1020)));
    try std.testing.expectEqual(0xc.773a2eac3001p-4, sin(@as(f64, 0xe.9f1e6p+112)));
    try std.testing.expectEqual(0x7.76d600e03152p-4, sin(@as(f64, 0xe.9f1e5p+112)));
    try std.testing.expectEqual(0xf.dfffd7bde0fb8p-4, sin(@as(f64, 0xe.9f1e5bc3bb88p+112)));
    try std.testing.expectEqual(-0x1.ffb679ba994b7p-4, sin(@as(f64, 0x4.7857dp+68)));
    try std.testing.expectEqual(-0x1.fecaff6878a11p-4, sin(@as(f64, 0x6.287cdp+0)));
    try std.testing.expectEqual(-0x1.fecb7e68ad6eap-4, sin(@as(f64, 0x6.287cc8p+0)));
    try std.testing.expectEqual(-0x1.fecb772e1b83p-4, sin(@as(f64, 0x6.287cc8749213p+0)));
    try std.testing.expectEqual(-0x1.fecb772e1b87p-4, sin(@as(f64, 0x6.287cc8749212cp+0)));
    try std.testing.expectEqual(-0xd.8f691a7a95428p-4, sin(@as(f64, -0x1.02e34cp+0)));
    try std.testing.expectEqual(-0x8.3bee07bc90768p-4, sin(@as(f64, 0xf.f0274p+4)));
    try std.testing.expectEqual(0x1.ffc6da9f1ffeep-4, sin(@as(f64, 0x3.042d88p+0)));
    try std.testing.expectEqual(0x1.d12ed0fffffep-12, sin(@as(f64, 0x1.d12ed2p-12)));
    try std.testing.expectEqual(-0x1.bf207c900d878p-4, sin(@as(f64, -0x6.e23688p+16)));
    try std.testing.expectEqual(-0x2.3e1f7a26f5944p-4, sin(@as(f64, -0x6.e2369p+16)));
    try std.testing.expectEqual(-0x1.feb6a3619e804p-4, sin(@as(f64, -0x6.e2368c006c018p+16)));
    try std.testing.expectEqual(-0x1.feb6a36596829p-4, sin(@as(f64, -0x6.e2368c006c01cp+16)));
    try std.testing.expectEqual(-0x1.4550689b93bbep-4, sin(@as(f64, 0x5.6a5008p+64)));
    try std.testing.expectEqual(-0x1.ee01db6bc8ef3p-4, sin(@as(f64, 0x5.6a5p+64)));
    try std.testing.expectEqual(0x6.5ea3351c9d9dcp-4, sin(@as(f64, 0x5.6a5005df4363cp+64)));
    try std.testing.expectEqual(0x2.f0e32ed649b32p-4, sin(@as(f64, 0x5.6a5005df43638p+64)));
    try std.testing.expectEqual(-0x8.599b32844aba8p-4, sin(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.452fc98b34e97p-8, sin(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x8.599b32844aba8p-4, sin(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.452fc98b34e97p-8, sin(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4p-128, sin(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, sin(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, sin(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, sin(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, sin(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, sin(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, sin(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, sin(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-152, sin(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sin(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0xf.fa2add3e58948p-4, sin(@as(f64, 0x1.8475e6p+0)));
    try std.testing.expectEqual(0xf.fa2adb8953aep-4, sin(@as(f64, 0x1.8475e4p+0)));
    try std.testing.expectEqual(0xf.fa2adcf9ea84p-4, sin(@as(f64, 0x1.8475e5afd4481p+0)));

    try std.testing.expectEqual(0x0p+0, sin(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sin(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x8.000003632116885p-4, sin(@as(f80, 0x8.60a92p-4)));
    try std.testing.expectEqual(0x7.fffff587e3a050dp-4, sin(@as(f80, 0x8.60a91p-4)));
    try std.testing.expectEqual(0x8.000000000000358p-4, sin(@as(f80, 0x8.60a91c16b9b3p-4)));
    try std.testing.expectEqual(0x7.ffffffffffffc6a8p-4, sin(@as(f80, 0x8.60a91c16b9b28p-4)));
    try std.testing.expectEqual(0x8.000000000000001p-4, sin(@as(f80, 0x8.60a91c16b9b2c24p-4)));
    try std.testing.expectEqual(0x8p-4, sin(@as(f80, 0x8.60a91c16b9b2c23p-4)));
    try std.testing.expectEqual(-0x7.fffff587e3a050dp-4, sin(@as(f80, -0x8.60a91p-4)));
    try std.testing.expectEqual(-0x8.000003632116885p-4, sin(@as(f80, -0x8.60a92p-4)));
    try std.testing.expectEqual(-0x7.ffffffffffffc6a8p-4, sin(@as(f80, -0x8.60a91c16b9b28p-4)));
    try std.testing.expectEqual(-0x8.000000000000358p-4, sin(@as(f80, -0x8.60a91c16b9b3p-4)));
    try std.testing.expectEqual(-0x8p-4, sin(@as(f80, -0x8.60a91c16b9b2c23p-4)));
    try std.testing.expectEqual(-0x8.000000000000001p-4, sin(@as(f80, -0x8.60a91c16b9b2c24p-4)));
    try std.testing.expectEqual(0xf.fffffffffffbb29p-4, sin(@as(f80, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0xf.fffffffffff32a3p-4, sin(@as(f80, 0x1.921fb4p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f80, 0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f80, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f80, 0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f80, 0x1.921fb54442d18468p+0)));
    try std.testing.expectEqual(-0xf.fffffffffffbb29p-4, sin(@as(f80, -0x1.921fb6p+0)));
    try std.testing.expectEqual(-0xf.fffffffffff32a3p-4, sin(@as(f80, -0x1.921fb4p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f80, -0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f80, -0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f80, -0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f80, -0x1.921fb54442d18468p+0)));
    try std.testing.expectEqual(-0x1.777a5cf72cec5fd6p-24, sin(@as(f80, 0x3.243f6cp+0)));
    try std.testing.expectEqual(0x2.8885a308d31063e4p-24, sin(@as(f80, 0x3.243f68p+0)));
    try std.testing.expectEqual(-0x1.72cece675d1fc8f8p-52, sin(@as(f80, 0x3.243f6a8885a32p+0)));
    try std.testing.expectEqual(0x8.d313198a2e03707p-56, sin(@as(f80, 0x3.243f6a8885a3p+0)));
    try std.testing.expectEqual(-0xe.ce675d1fc8f8cbbp-68, sin(@as(f80, 0x3.243f6a8885a308d4p+0)));
    try std.testing.expectEqual(0x3.13198a2e03707344p-64, sin(@as(f80, 0x3.243f6a8885a308dp+0)));
    try std.testing.expectEqual(0x1.777a5cf72cec5fd6p-24, sin(@as(f80, -0x3.243f6cp+0)));
    try std.testing.expectEqual(-0x2.8885a308d31063e4p-24, sin(@as(f80, -0x3.243f68p+0)));
    try std.testing.expectEqual(0x1.72cece675d1fc8f8p-52, sin(@as(f80, -0x3.243f6a8885a32p+0)));
    try std.testing.expectEqual(-0x8.d313198a2e03707p-56, sin(@as(f80, -0x3.243f6a8885a3p+0)));
    try std.testing.expectEqual(0xe.ce675d1fc8f8cbbp-68, sin(@as(f80, -0x3.243f6a8885a308d4p+0)));
    try std.testing.expectEqual(-0x3.13198a2e03707344p-64, sin(@as(f80, -0x3.243f6a8885a308dp+0)));
    try std.testing.expectEqual(0xa.e7fe0b5fc786b2ep-4, sin(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(-0xc.143e153b0701e8p-8, sin(@as(f80, 0x2p+64)));
    try std.testing.expectEqual(0xc.143e153b0701e8p-8, sin(@as(f80, -0x2p+64)));
    try std.testing.expectEqual(-0x1.1e7cf9ec10916c24p-24, sin(@as(f80, 0xb.fa09ap+100)));
    try std.testing.expectEqual(0xb.7fb600275877a6p-4, sin(@as(f80, 0xc.d4967p-4)));
    try std.testing.expectEqual(0xb.7fb5f50739fa5f9p-4, sin(@as(f80, 0xc.d4966p-4)));
    try std.testing.expectEqual(0xb.7fb5fe7769793e6p-4, sin(@as(f80, 0xc.d4966d92d171p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e56p-4, sin(@as(f80, 0xc.d4966d92d1708p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e74p-4, sin(@as(f80, 0xc.d4966d92d17082ap-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e73p-4, sin(@as(f80, 0xc.d4966d92d170829p-4)));
    try std.testing.expectEqual(0x3.fe4780403e8078ccp-4, sin(@as(f80, 0x4.093388p-4)));
    try std.testing.expectEqual(0x3.fe4778810e026fep-4, sin(@as(f80, 0x4.09338p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7958dccp-4, sin(@as(f80, 0x4.093385688a2d4p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7954fd4p-4, sin(@as(f80, 0x4.093385688a2dp-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7956438p-4, sin(@as(f80, 0x4.093385688a2d151p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc795643p-4, sin(@as(f80, 0x4.093385688a2d1508p-4)));
    try std.testing.expectEqual(-0x4.cd7e86c4077bf0ep-4, sin(@as(f80, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(-0xb.becc47ab1b8c708p-4, sin(@as(f80, 0x2.1e19ep+72)));
    try std.testing.expectEqual(-0xd.a29d5bb5f9cb87dp-4, sin(@as(f80, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(-0x8.599b32844aba907p-4, sin(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x9.0292465edbaff2dp-4, sin(@as(f80, 0x8p+1020)));
    try std.testing.expectEqual(-0x8.599b32844aba907p-4, sin(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.452fc98b34e96b62p-8, sin(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x6.3ad4b2136cc6882p-4, sin(@as(f80, 0x8p+16380)));
    try std.testing.expectEqual(0x6.0b8d19579bf2db6p-4, sin(@as(f80, 0x1p+120)));
    try std.testing.expectEqual(0x9.f963166f215eb89p-4, sin(@as(f80, 0x8p+124)));
    try std.testing.expectEqual(0xc.6fa5c5665984d89p-8, sin(@as(f80, 0xf.ffffcp+124)));
    try std.testing.expectEqual(-0x8.599b32844aba907p-4, sin(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x7.f13d78eabb76b8a8p-4, sin(@as(f80, 0x4p+48)));
    try std.testing.expectEqual(-0xf.c777c6b36a750a6p-4, sin(@as(f80, 0x1p+28)));
    try std.testing.expectEqual(0xc.dbc19bb4a58a819p-4, sin(@as(f80, 0xe.ef3bp-4)));
    try std.testing.expectEqual(0xc.dbc1922f1d9f2c7p-4, sin(@as(f80, 0xe.ef3afp-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3b0c3p-4, sin(@as(f80, 0xe.ef3af1b5d8008p-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3acp-4, sin(@as(f80, 0xe.ef3af1b5d8p-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3ac01p-4, sin(@as(f80, 0xe.ef3af1b5d800001p-4)));
    try std.testing.expectEqual(0xb.93255854754a36dp-4, sin(@as(f80, 0x2.553538p+0)));
    try std.testing.expectEqual(0xb.932584840806c61p-4, sin(@as(f80, 0x2.553534p+0)));
    try std.testing.expectEqual(0xb.93255eeda1024p-4, sin(@as(f80, 0x2.5535376715bap+0)));
    try std.testing.expectEqual(0xb.93255eeda103a18p-4, sin(@as(f80, 0x2.5535376715b9ep+0)));
    try std.testing.expectEqual(0xb.93255eeda102403p-4, sin(@as(f80, 0x2.5535376715b9fffcp+0)));
    try std.testing.expectEqual(-0x9.10bb448d3cb0167p-4, sin(@as(f80, 0x3.be736p+0)));
    try std.testing.expectEqual(-0x9.10bb0fd0c39d5ap-4, sin(@as(f80, 0x3.be735cp+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec9cp-4, sin(@as(f80, 0x3.be735c19beap+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec81a2p-4, sin(@as(f80, 0x3.be735c19be9fep+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec9bfdp-4, sin(@as(f80, 0x3.be735c19be9ffffcp+0)));
    try std.testing.expectEqual(-0xb.4351fdda3d81514p-4, sin(@as(f80, 0x3.ec2a04p+0)));
    try std.testing.expectEqual(-0xb.4351d06546e7181p-4, sin(@as(f80, 0x3.ec2ap+0)));
    try std.testing.expectEqual(-0xb.4351eaad0984abap-4, sin(@as(f80, 0x3.ec2a0250032a2p+0)));
    try std.testing.expectEqual(-0xb.4351eaad09834p-4, sin(@as(f80, 0x3.ec2a0250032ap+0)));
    try std.testing.expectEqual(-0xb.4351eaad0983403p-4, sin(@as(f80, 0x3.ec2a0250032a0004p+0)));
    try std.testing.expectEqual(-0xc.d2639f1afc7f46ap-4, sin(@as(f80, 0x4.1237e8p+0)));
    try std.testing.expectEqual(-0xc.d2635289f074d2bp-4, sin(@as(f80, 0x4.1237ep+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf59249p-4, sin(@as(f80, 0x4.1237e153f7084p+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf56cp-4, sin(@as(f80, 0x4.1237e153f708p+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf56c05p-4, sin(@as(f80, 0x4.1237e153f7080008p+0)));
    try std.testing.expectEqual(-0xf.f4f46a017cb884p-4, sin(@as(f80, 0x4.c92d1p+0)));
    try std.testing.expectEqual(-0xf.f4f4736648dc2ap-4, sin(@as(f80, 0x4.c92d08p+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f27f4ep-4, sin(@as(f80, 0x4.c92d0ffa4bf04p+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f284p-4, sin(@as(f80, 0x4.c92d0ffa4bfp+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f283ffp-4, sin(@as(f80, 0x4.c92d0ffa4bf00008p+0)));
    try std.testing.expectEqual(-0x4.b6f60ca8d4150bcp-4, sin(@as(f80, 0x5.fbec78p+0)));
    try std.testing.expectEqual(-0x4.b6f686f9ea12e8fp-4, sin(@as(f80, 0x5.fbec7p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a0cd78p-4, sin(@as(f80, 0x5.fbec7477d4a84p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a4ap-4, sin(@as(f80, 0x5.fbec7477d4a8p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a49f88p-4, sin(@as(f80, 0x5.fbec7477d4a80008p+0)));
    try std.testing.expectEqual(0xd.76aa47848677021p-4, sin(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xe.8c7b7568da22efdp-4, sin(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(0x2.42070db6daab69e4p-4, sin(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(-0xc.1bdceeee0f57386p-4, sin(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(-0xf.57c0faf04c99913p-4, sin(@as(f80, 0x5p+0)));
    try std.testing.expectEqual(-0x4.787c62ac28b00e98p-4, sin(@as(f80, 0x6p+0)));
    try std.testing.expectEqual(0xa.83046136850421ep-4, sin(@as(f80, 0x7p+0)));
    try std.testing.expectEqual(0xf.d469501467bd75p-4, sin(@as(f80, 0x8p+0)));
    try std.testing.expectEqual(0x6.98098d830be42488p-4, sin(@as(f80, 0x9p+0)));
    try std.testing.expectEqual(-0x8.b44f7af9a7a92cep-4, sin(@as(f80, 0xap+0)));
    try std.testing.expectEqual(-0x5.595d7e536fe35ed8p-4, sin(@as(f80, 0x1.200148p+32)));
    try std.testing.expectEqual(0x4.220ff25f5cf02a48p-4, sin(@as(f80, 0x1.200146p+32)));
    try std.testing.expectEqual(-0x6.444fda50019f9f58p-4, sin(@as(f80, 0x1.2001469775ce6p+32)));
    try std.testing.expectEqual(0x8.599b32844aba907p-4, sin(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.e00885042dd770dp-4, sin(@as(f80, -0x3.3de320f6be87ep+1020)));
    try std.testing.expectEqual(0xc.773a2eac3000ddfp-4, sin(@as(f80, 0xe.9f1e6p+112)));
    try std.testing.expectEqual(0x7.76d600e031521b8p-4, sin(@as(f80, 0xe.9f1e5p+112)));
    try std.testing.expectEqual(0xf.dfffd7bde0fb4ecp-4, sin(@as(f80, 0xe.9f1e5bc3bb88p+112)));
    try std.testing.expectEqual(-0x1.ffb679ba994b7618p-4, sin(@as(f80, 0x4.7857dp+68)));
    try std.testing.expectEqual(-0x1.fecaff6878a10ce6p-4, sin(@as(f80, 0x6.287cdp+0)));
    try std.testing.expectEqual(-0x1.fecb7e68ad6e9c4p-4, sin(@as(f80, 0x6.287cc8p+0)));
    try std.testing.expectEqual(-0x1.fecb772e1b8300e6p-4, sin(@as(f80, 0x6.287cc8749213p+0)));
    try std.testing.expectEqual(-0x1.fecb772e1b86f8e8p-4, sin(@as(f80, 0x6.287cc8749212cp+0)));
    try std.testing.expectEqual(-0x1.fecb772e1b848bcap-4, sin(@as(f80, 0x6.287cc8749212e72p+0)));
    try std.testing.expectEqual(-0xd.8f691a7a95426p-4, sin(@as(f80, -0x1.02e34cp+0)));
    try std.testing.expectEqual(-0x8.3bee07bc9076425p-4, sin(@as(f80, 0xf.f0274p+4)));
    try std.testing.expectEqual(0x1.ffc6da9f1ffed896p-4, sin(@as(f80, 0x3.042d88p+0)));
    try std.testing.expectEqual(0x1.d12ed0fffffdfe1p-12, sin(@as(f80, 0x1.d12ed2p-12)));
    try std.testing.expectEqual(-0x1.bf207c900d877cb8p-4, sin(@as(f80, -0x6.e23688p+16)));
    try std.testing.expectEqual(-0x2.3e1f7a26f594337p-4, sin(@as(f80, -0x6.e2369p+16)));
    try std.testing.expectEqual(-0x1.feb6a3619e803c4ap-4, sin(@as(f80, -0x6.e2368c006c018p+16)));
    try std.testing.expectEqual(-0x1.feb6a36596828878p-4, sin(@as(f80, -0x6.e2368c006c01cp+16)));
    try std.testing.expectEqual(-0x1.feb6a361c0bb501cp-4, sin(@as(f80, -0x6.e2368c006c018228p+16)));
    try std.testing.expectEqual(-0x1.4550689b93bbe14p-4, sin(@as(f80, 0x5.6a5008p+64)));
    try std.testing.expectEqual(-0x1.ee01db6bc8ef288cp-4, sin(@as(f80, 0x5.6a5p+64)));
    try std.testing.expectEqual(0x6.5ea3351c9d9da32p-4, sin(@as(f80, 0x5.6a5005df4363cp+64)));
    try std.testing.expectEqual(0x2.f0e32ed649b32644p-4, sin(@as(f80, 0x5.6a5005df43638p+64)));
    try std.testing.expectEqual(-0xa.8640e82e7924ecp-4, sin(@as(f80, 0x5.6a5005df43638338p+64)));
    try std.testing.expectEqual(0xd.7457bd2255e689fp-4, sin(@as(f80, 0x5.6a5005df4363833p+64)));
    try std.testing.expectEqual(-0x8.599b32844aba907p-4, sin(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.452fc98b34e96b62p-8, sin(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.dfd9d4b6d0e5f7cp-4, sin(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x8.599b32844aba907p-4, sin(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.452fc98b34e96b62p-8, sin(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xf.dfd9d4b6d0e5f7cp-4, sin(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4p-128, sin(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, sin(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, sin(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, sin(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, sin(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, sin(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, sin(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, sin(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, sin(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, sin(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, sin(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, sin(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, sin(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-152, sin(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sin(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, sin(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0xf.fa2add3e58948d1p-4, sin(@as(f80, 0x1.8475e6p+0)));
    try std.testing.expectEqual(0xf.fa2adb8953ae262p-4, sin(@as(f80, 0x1.8475e4p+0)));
    try std.testing.expectEqual(0xf.fa2adcf9ea83dbep-4, sin(@as(f80, 0x1.8475e5afd4481p+0)));

    try std.testing.expectEqual(0x0p+0, sin(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sin(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x8.0000036321168852c4130c64b4cp-4, sin(@as(f128, 0x8.60a92p-4)));
    try std.testing.expectEqual(0x7.fffff587e3a050cf967fba7bc728p-4, sin(@as(f128, 0x8.60a91p-4)));
    try std.testing.expectEqual(0x8.00000000000035858118fd5157ep-4, sin(@as(f128, 0x8.60a91c16b9b3p-4)));
    try std.testing.expectEqual(0x7.ffffffffffffc6ab95779c1eae0cp-4, sin(@as(f128, 0x8.60a91c16b9b28p-4)));
    try std.testing.expectEqual(0x8.000000000000000b5feca2ed673p-4, sin(@as(f128, 0x8.60a91c16b9b2c24p-4)));
    try std.testing.expectEqual(0x7.fffffffffffffffd84af2ec140dcp-4, sin(@as(f128, 0x8.60a91c16b9b2c23p-4)));
    try std.testing.expectEqual(0x8p-4, sin(@as(f128, 0x8.60a91c16b9b2c232dd99707ab3d8p-4)));
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffffffcp-4, sin(@as(f128, 0x8.60a91c16b9b2c232dd99707ab3dp-4)));
    try std.testing.expectEqual(0x8.000000000000000000000000002p-4, sin(@as(f128, 0x8.60a91c16b9b2c232dd99707ab4p-4)));
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffffcacp-4, sin(@as(f128, 0x8.60a91c16b9b2c232dd99707abp-4)));
    try std.testing.expectEqual(-0x7.fffff587e3a050cf967fba7bc728p-4, sin(@as(f128, -0x8.60a91p-4)));
    try std.testing.expectEqual(-0x8.0000036321168852c4130c64b4cp-4, sin(@as(f128, -0x8.60a92p-4)));
    try std.testing.expectEqual(-0x7.ffffffffffffc6ab95779c1eae0cp-4, sin(@as(f128, -0x8.60a91c16b9b28p-4)));
    try std.testing.expectEqual(-0x8.00000000000035858118fd5157ep-4, sin(@as(f128, -0x8.60a91c16b9b3p-4)));
    try std.testing.expectEqual(-0x7.fffffffffffffffd84af2ec140dcp-4, sin(@as(f128, -0x8.60a91c16b9b2c23p-4)));
    try std.testing.expectEqual(-0x8.000000000000000b5feca2ed673p-4, sin(@as(f128, -0x8.60a91c16b9b2c24p-4)));
    try std.testing.expectEqual(-0x7.fffffffffffffffffffffffffffcp-4, sin(@as(f128, -0x8.60a91c16b9b2c232dd99707ab3dp-4)));
    try std.testing.expectEqual(-0x8p-4, sin(@as(f128, -0x8.60a91c16b9b2c232dd99707ab3d8p-4)));
    try std.testing.expectEqual(-0x7.fffffffffffffffffffffffffcacp-4, sin(@as(f128, -0x8.60a91c16b9b2c232dd99707abp-4)));
    try std.testing.expectEqual(-0x8.000000000000000000000000002p-4, sin(@as(f128, -0x8.60a91c16b9b2c232dd99707ab4p-4)));
    try std.testing.expectEqual(0xf.fffffffffffbb290924e3a114988p-4, sin(@as(f128, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0xf.fffffffffff32a3661c108e136d8p-4, sin(@as(f128, 0x1.921fb4p+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffbdp-4, sin(@as(f128, 0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffff68p-4, sin(@as(f128, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f128, 0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f128, 0x1.921fb54442d18468p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f128, 0x1.921fb54442d18469898cc51701b9p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f128, 0x1.921fb54442d18469898cc51701b8p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f128, 0x1.921fb54442d18469898cc51702p+0)));
    try std.testing.expectEqual(0x1p+0, sin(@as(f128, 0x1.921fb54442d18469898cc517018p+0)));
    try std.testing.expectEqual(-0xf.fffffffffffbb290924e3a114988p-4, sin(@as(f128, -0x1.921fb6p+0)));
    try std.testing.expectEqual(-0xf.fffffffffff32a3661c108e136d8p-4, sin(@as(f128, -0x1.921fb4p+0)));
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffffbdp-4, sin(@as(f128, -0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(-0xf.ffffffffffffffffffffffffff68p-4, sin(@as(f128, -0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f128, -0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f128, -0x1.921fb54442d18468p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f128, -0x1.921fb54442d18469898cc51701b9p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f128, -0x1.921fb54442d18469898cc51701b8p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f128, -0x1.921fb54442d18469898cc51702p+0)));
    try std.testing.expectEqual(-0x1p+0, sin(@as(f128, -0x1.921fb54442d18469898cc517018p+0)));
    try std.testing.expectEqual(-0x1.777a5cf72cec5fd61896cb4f40d2p-24, sin(@as(f128, 0x3.243f6cp+0)));
    try std.testing.expectEqual(0x2.8885a308d31063e2b6c62b7f4d6cp-24, sin(@as(f128, 0x3.243f68p+0)));
    // try std.testing.expectEqual(-0x1.72cece675d1fc8f8cbb5bf6c7d5cp-52, sin(@as(f128, 0x3.243f6a8885a32p+0)));
    try std.testing.expectEqual(0x8.d313198a2e03707344a4093821b8p-56, sin(@as(f128, 0x3.243f6a8885a3p+0)));
    try std.testing.expectEqual(-0xe.ce675d1fc8f8cbb5bf6c7ddd661p-68, sin(@as(f128, 0x3.243f6a8885a308d4p+0)));
    try std.testing.expectEqual(0x3.13198a2e03707344a409382229ap-64, sin(@as(f128, 0x3.243f6a8885a308dp+0)));
    try std.testing.expectEqual(-0x1.8cbb5bf6c7ddd660ce2ff7d10567p-112, sin(@as(f128, 0x3.243f6a8885a308d313198a2e0372p+0)));
    try std.testing.expectEqual(0x7.344a4093822299f31d0082efa99p-116, sin(@as(f128, 0x3.243f6a8885a308d313198a2e037p+0)));
    try std.testing.expectEqual(-0x8.f8cbb5bf6c7ddd660ce2ff7d1058p-108, sin(@as(f128, 0x3.243f6a8885a308d313198a2e04p+0)));
    try std.testing.expectEqual(0x7.07344a4093822299f31d0082efa8p-108, sin(@as(f128, 0x3.243f6a8885a308d313198a2e03p+0)));
    try std.testing.expectEqual(0x1.777a5cf72cec5fd61896cb4f40d2p-24, sin(@as(f128, -0x3.243f6cp+0)));
    try std.testing.expectEqual(-0x2.8885a308d31063e2b6c62b7f4d6cp-24, sin(@as(f128, -0x3.243f68p+0)));
    // try std.testing.expectEqual(0x1.72cece675d1fc8f8cbb5bf6c7d5cp-52, sin(@as(f128, -0x3.243f6a8885a32p+0)));
    try std.testing.expectEqual(-0x8.d313198a2e03707344a4093821b8p-56, sin(@as(f128, -0x3.243f6a8885a3p+0)));
    try std.testing.expectEqual(0xe.ce675d1fc8f8cbb5bf6c7ddd661p-68, sin(@as(f128, -0x3.243f6a8885a308d4p+0)));
    try std.testing.expectEqual(-0x3.13198a2e03707344a409382229ap-64, sin(@as(f128, -0x3.243f6a8885a308dp+0)));
    try std.testing.expectEqual(0x1.8cbb5bf6c7ddd660ce2ff7d10567p-112, sin(@as(f128, -0x3.243f6a8885a308d313198a2e0372p+0)));
    try std.testing.expectEqual(-0x7.344a4093822299f31d0082efa99p-116, sin(@as(f128, -0x3.243f6a8885a308d313198a2e037p+0)));
    try std.testing.expectEqual(0x8.f8cbb5bf6c7ddd660ce2ff7d1058p-108, sin(@as(f128, -0x3.243f6a8885a308d313198a2e04p+0)));
    try std.testing.expectEqual(-0x7.07344a4093822299f31d0082efa8p-108, sin(@as(f128, -0x3.243f6a8885a308d313198a2e03p+0)));
    try std.testing.expectEqual(0xa.e7fe0b5fc786b2d966e1d6af1408p-4, sin(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(-0xc.143e153b0701e800f9b8a47b75b8p-8, sin(@as(f128, 0x2p+64)));
    try std.testing.expectEqual(0xc.143e153b0701e800f9b8a47b75b8p-8, sin(@as(f128, -0x2p+64)));
    try std.testing.expectEqual(-0x1.1e7cf9ec10916c247834e44dabf1p-24, sin(@as(f128, 0xb.fa09ap+100)));
    try std.testing.expectEqual(0xb.7fb600275877a60011766c8a3178p-4, sin(@as(f128, 0xc.d4967p-4)));
    try std.testing.expectEqual(0xb.7fb5f50739fa5f8acc8f4f3f1b3p-4, sin(@as(f128, 0xc.d4966p-4)));
    try std.testing.expectEqual(0xb.7fb5fe7769793e65c978bd3cef98p-4, sin(@as(f128, 0xc.d4966d92d171p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e564d5ae94f8cb08p-4, sin(@as(f128, 0xc.d4966d92d1708p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e7381aae7a4c30dp-4, sin(@as(f128, 0xc.d4966d92d17082ap-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e72cfa9001072848p-4, sin(@as(f128, 0xc.d4966d92d170829p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e732912810356318p-4, sin(@as(f128, 0xc.d4966d92d17082980965c1a663c8p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e732912810356318p-4, sin(@as(f128, 0xc.d4966d92d17082980965c1a663cp-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e73291281035634p-4, sin(@as(f128, 0xc.d4966d92d17082980965c1a664p-4)));
    try std.testing.expectEqual(0xb.7fb5fe776978e732912810356078p-4, sin(@as(f128, 0xc.d4966d92d17082980965c1a66p-4)));
    try std.testing.expectEqual(0x3.fe4780403e8078ca89790118cb8cp-4, sin(@as(f128, 0x4.093388p-4)));
    try std.testing.expectEqual(0x3.fe4778810e026fe1e37f141da492p-4, sin(@as(f128, 0x4.09338p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7958dcbe52ad14f88f2p-4, sin(@as(f128, 0x4.093385688a2d4p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7954fd2613bf1f96c24p-4, sin(@as(f128, 0x4.093385688a2dp-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7956437bee74ef98326p-4, sin(@as(f128, 0x4.093385688a2d151p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc795642fffb6d11d9862p-4, sin(@as(f128, 0x4.093385688a2d1508p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7956433e0084536148p-4, sin(@as(f128, 0x4.093385688a2d150c00bf42a08e84p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7956433e0084536147cp-4, sin(@as(f128, 0x4.093385688a2d150c00bf42a08e8p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7956433e008453615fp-4, sin(@as(f128, 0x4.093385688a2d150c00bf42a09p-4)));
    try std.testing.expectEqual(0x3.fe477dbdc7956433e008453614p-4, sin(@as(f128, 0x4.093385688a2d150c00bf42a08ep-4)));
    try std.testing.expectEqual(-0x4.cd7e86c4077bf0debc87d70d196p-4, sin(@as(f128, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(-0xb.becc47ab1b8c70793712c4ff2bcp-4, sin(@as(f128, 0x2.1e19ep+72)));
    try std.testing.expectEqual(-0xd.a29d5bb5f9cb87d14de41dc991fp-4, sin(@as(f128, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(-0x8.599b32844aba906cee446be04998p-4, sin(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x9.0292465edbaff2d2e64a2845e558p-4, sin(@as(f128, 0x8p+1020)));
    try std.testing.expectEqual(-0x8.599b32844aba906cee446be04998p-4, sin(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.452fc98b34e96b61139b09a7c84ap-8, sin(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x6.3ad4b2136cc6881f0ca607c7946p-4, sin(@as(f128, 0x8p+16380)));
    try std.testing.expectEqual(-0xe.f1a3e1dc468a921dddb4e37fbe6p-4, sin(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x6.0b8d19579bf2db5e5f1aa933f37cp-4, sin(@as(f128, 0x1p+120)));
    try std.testing.expectEqual(0x9.f963166f215eb89381836cafaa3p-4, sin(@as(f128, 0x8p+124)));
    try std.testing.expectEqual(0xc.6fa5c5665984d8892761be1537b8p-8, sin(@as(f128, 0xf.ffffcp+124)));
    try std.testing.expectEqual(-0x8.599b32844aba906cee446be04998p-4, sin(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x7.f13d78eabb76b8a986d98d6703e8p-4, sin(@as(f128, 0x4p+48)));
    try std.testing.expectEqual(-0xf.c777c6b36a750a5fdeb8807a156p-4, sin(@as(f128, 0x1p+28)));
    try std.testing.expectEqual(0xc.dbc19bb4a58a818888614bb1337p-4, sin(@as(f128, 0xe.ef3bp-4)));
    try std.testing.expectEqual(0xc.dbc1922f1d9f2c71bb4e0682653p-4, sin(@as(f128, 0xe.ef3afp-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3b0c2c3f7f22cb9868p-4, sin(@as(f128, 0xe.ef3af1b5d8008p-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3abffffffc0a0c496p-4, sin(@as(f128, 0xe.ef3af1b5d8p-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3ac0098583fa6f6148p-4, sin(@as(f128, 0xe.ef3af1b5d800001p-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3abffffffc0a0c4a18p-4, sin(@as(f128, 0xe.ef3af1b5d800000000000000014p-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3abffffffc0a0c4a18p-4, sin(@as(f128, 0xe.ef3af1b5d8000000000000000138p-4)));
    try std.testing.expectEqual(0xc.dbc19333ad3abffffffc0a0c4bcp-4, sin(@as(f128, 0xe.ef3af1b5d80000000000000004p-4)));
    try std.testing.expectEqual(0xb.93255854754a36d261581d7c0da8p-4, sin(@as(f128, 0x2.553538p+0)));
    try std.testing.expectEqual(0xb.932584840806c6090aefd5f2505p-4, sin(@as(f128, 0x2.553534p+0)));
    try std.testing.expectEqual(0xb.93255eeda10240000004f6fd44f8p-4, sin(@as(f128, 0x2.5535376715bap+0)));
    try std.testing.expectEqual(0xb.93255eeda103a17c97f0fb51252p-4, sin(@as(f128, 0x2.5535376715b9ep+0)));
    try std.testing.expectEqual(0xb.93255eeda102402c2f97f47dcf78p-4, sin(@as(f128, 0x2.5535376715b9fffcp+0)));
    try std.testing.expectEqual(0xb.93255eeda10240000004f6fd4acp-4, sin(@as(f128, 0x2.5535376715b9ffffffffffffff7ap+0)));
    try std.testing.expectEqual(0xb.93255eeda10240000004f6fd4ad8p-4, sin(@as(f128, 0x2.5535376715b9ffffffffffffff78p+0)));
    try std.testing.expectEqual(0xb.93255eeda10240000004f6fd5008p-4, sin(@as(f128, 0x2.5535376715b9ffffffffffffffp+0)));
    try std.testing.expectEqual(-0x9.10bb448d3cb0166e220f3af793c8p-4, sin(@as(f128, 0x3.be736p+0)));
    try std.testing.expectEqual(-0x9.10bb0fd0c39d59f8b7898e86413p-4, sin(@as(f128, 0x3.be735cp+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec9c000000287a188ap-4, sin(@as(f128, 0x3.be735c19beap+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec81a1c3545a17905f8p-4, sin(@as(f128, 0x3.be735c19be9fep+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec9bfcb43893004c39p-4, sin(@as(f128, 0x3.be735c19be9ffffcp+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec9c000000287a1878p-4, sin(@as(f128, 0x3.be735c19be9fffffffffffffffeap+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec9c000000287a18768p-4, sin(@as(f128, 0x3.be735c19be9fffffffffffffffe8p+0)));
    try std.testing.expectEqual(-0x9.10bb11242ec9c000000287a17b7p-4, sin(@as(f128, 0x3.be735c19be9fffffffffffffffp+0)));
    try std.testing.expectEqual(-0xb.4351fdda3d81513dedde4fd2cd28p-4, sin(@as(f128, 0x3.ec2a04p+0)));
    try std.testing.expectEqual(-0xb.4351d06546e7181306453a5b2ec8p-4, sin(@as(f128, 0x3.ec2ap+0)));
    try std.testing.expectEqual(-0xb.4351eaad0984aba7b4606b57ad68p-4, sin(@as(f128, 0x3.ec2a0250032a2p+0)));
    try std.testing.expectEqual(-0xb.4351eaad09833fffffff47a70dd8p-4, sin(@as(f128, 0x3.ec2a0250032ap+0)));
    try std.testing.expectEqual(-0xb.4351eaad0983402d74f5d3cb83fp-4, sin(@as(f128, 0x3.ec2a0250032a0004p+0)));
    try std.testing.expectEqual(-0xb.4351eaad09833fffffff47a712e8p-4, sin(@as(f128, 0x3.ec2a0250032a0000000000000072p+0)));
    try std.testing.expectEqual(-0xb.4351eaad09833fffffff47a712dp-4, sin(@as(f128, 0x3.ec2a0250032a000000000000007p+0)));
    try std.testing.expectEqual(-0xb.4351eaad09833fffffff47a71938p-4, sin(@as(f128, 0x3.ec2a0250032a00000000000001p+0)));
    try std.testing.expectEqual(-0xc.d2639f1afc7f4698649bd4cf58c8p-4, sin(@as(f128, 0x4.1237e8p+0)));
    try std.testing.expectEqual(-0xc.d2635289f074d2b2df33331ddd4p-4, sin(@as(f128, 0x4.1237ep+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf59248868df5425cd98p-4, sin(@as(f128, 0x4.1237e153f7084p+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf56bffffffb0fcac28p-4, sin(@as(f128, 0x4.1237e153f708p+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf56c04c91082c534dfp-4, sin(@as(f128, 0x4.1237e153f7080008p+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf56bffffffb0fcac2a8p-4, sin(@as(f128, 0x4.1237e153f7080000000000000004p+0)));
    try std.testing.expectEqual(-0xc.d2635f3faf56bffffffb0fcad5ap-4, sin(@as(f128, 0x4.1237e153f70800000000000002p+0)));
    try std.testing.expectEqual(-0xf.f4f46a017cb883f95b64a197e97p-4, sin(@as(f128, 0x4.c92d1p+0)));
    try std.testing.expectEqual(-0xf.f4f4736648dc2a042045142724ap-4, sin(@as(f128, 0x4.c92d08p+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f27f4d98f09e3c6cf38p-4, sin(@as(f128, 0x4.c92d0ffa4bf04p+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f2840000005480a4ap-4, sin(@as(f128, 0x4.c92d0ffa4bfp+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f283ff69b37289dc1ap-4, sin(@as(f128, 0x4.c92d0ffa4bf00008p+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f2840000005480a4958p-4, sin(@as(f128, 0x4.c92d0ffa4bf0000000000000008cp+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f2840000005480a496p-4, sin(@as(f128, 0x4.c92d0ffa4bf00000000000000088p+0)));
    try std.testing.expectEqual(-0xf.f4f46a082f2840000005480a47ap-4, sin(@as(f128, 0x4.c92d0ffa4bf000000000000002p+0)));
    try std.testing.expectEqual(-0x4.b6f60ca8d4150bc1732b1580fc8p-4, sin(@as(f128, 0x5.fbec78p+0)));
    try std.testing.expectEqual(-0x4.b6f686f9ea12e8ec702d198cbfa8p-4, sin(@as(f128, 0x5.fbec7p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a0cd774f8372f49334p-4, sin(@as(f128, 0x5.fbec7477d4a84p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a49fffffff9a1c7cf8p-4, sin(@as(f128, 0x5.fbec7477d4a8p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a49f85aee98a9798p-4, sin(@as(f128, 0x5.fbec7477d4a80008p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a49fffffff9a1c73a8p-4, sin(@as(f128, 0x5.fbec7477d4a8000000000000009cp+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a49fffffff9a1c73e4p-4, sin(@as(f128, 0x5.fbec7477d4a80000000000000098p+0)));
    try std.testing.expectEqual(-0x4.b6f642a935a49fffffff9a1c5e64p-4, sin(@as(f128, 0x5.fbec7477d4a800000000000002p+0)));
    try std.testing.expectEqual(0xd.76aa47848677020c6e9e909c50fp-4, sin(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xe.8c7b7568da22efd5c240c4004e5p-4, sin(@as(f128, 0x2p+0)));
    // try std.testing.expectEqual(0x2.42070db6daab69e3902e8468315p-4, sin(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(-0xc.1bdceeee0f5738674c02ad072288p-4, sin(@as(f128, 0x4p+0)));
    try std.testing.expectEqual(-0xf.57c0faf04c999135789f2ab1de3p-4, sin(@as(f128, 0x5p+0)));
    try std.testing.expectEqual(-0x4.787c62ac28b00e989105113d2884p-4, sin(@as(f128, 0x6p+0)));
    try std.testing.expectEqual(0xa.83046136850421dd444208fd7788p-4, sin(@as(f128, 0x7p+0)));
    try std.testing.expectEqual(0xf.d469501467bd74fb15853828467p-4, sin(@as(f128, 0x8p+0)));
    try std.testing.expectEqual(0x6.98098d830be42487274a8291c67cp-4, sin(@as(f128, 0x9p+0)));
    try std.testing.expectEqual(-0x8.b44f7af9a7a92ce7fb22be024e2p-4, sin(@as(f128, 0xap+0)));
    try std.testing.expectEqual(-0x5.595d7e536fe35edbe2ad0df9d94p-4, sin(@as(f128, 0x1.200148p+32)));
    try std.testing.expectEqual(0x4.220ff25f5cf02a464dbb3a679ccp-4, sin(@as(f128, 0x1.200146p+32)));
    try std.testing.expectEqual(-0x6.444fda50019f9f5ba3779ca706p-4, sin(@as(f128, 0x1.2001469775ce6p+32)));
    try std.testing.expectEqual(0x8.599b32844aba906cee446be04998p-4, sin(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.e00885042dd770c93962abdb61f8p-4, sin(@as(f128, -0x3.3de320f6be87ep+1020)));
    try std.testing.expectEqual(0xc.773a2eac3000ddec0c69e7ddef68p-4, sin(@as(f128, 0xe.9f1e6p+112)));
    try std.testing.expectEqual(0x7.76d600e031521b7cc3cd579a135p-4, sin(@as(f128, 0xe.9f1e5p+112)));
    try std.testing.expectEqual(0xf.dfffd7bde0fb4ec139784e3b799p-4, sin(@as(f128, 0xe.9f1e5bc3bb88p+112)));
    try std.testing.expectEqual(-0x1.ffb679ba994b76173f9040637ff9p-4, sin(@as(f128, 0x4.7857dp+68)));
    // try std.testing.expectEqual(-0x1.fecaff6878a10ce5d42fde40e7p-4, sin(@as(f128, 0x6.287cdp+0)));
    try std.testing.expectEqual(-0x1.fecb7e68ad6e9c3f77c1915bc919p-4, sin(@as(f128, 0x6.287cc8p+0)));
    // try std.testing.expectEqual(-0x1.fecb772e1b8300e5ab16d9008ea9p-4, sin(@as(f128, 0x6.287cc8749213p+0)));
    try std.testing.expectEqual(-0x1.fecb772e1b86f8e74fbeae63ee4cp-4, sin(@as(f128, 0x6.287cc8749212cp+0)));
    try std.testing.expectEqual(-0x1.fecb772e1b848bca4e961470b22p-4, sin(@as(f128, 0x6.287cc8749212e72p+0)));
    try std.testing.expectEqual(-0xd.8f691a7a95425ffcb89dc2b97cep-4, sin(@as(f128, -0x1.02e34cp+0)));
    try std.testing.expectEqual(-0x8.3bee07bc9076424bef274717106p-4, sin(@as(f128, 0xf.f0274p+4)));
    // try std.testing.expectEqual(0x1.ffc6da9f1ffed895f9fa424ba91p-4, sin(@as(f128, 0x3.042d88p+0)));
    try std.testing.expectEqual(0x1.d12ed0fffffdfe0f0008c8b28233p-12, sin(@as(f128, 0x1.d12ed2p-12)));
    try std.testing.expectEqual(-0x1.bf207c900d877cb73f555829e3f1p-4, sin(@as(f128, -0x6.e23688p+16)));
    // try std.testing.expectEqual(-0x2.3e1f7a26f594336f9e583b26bbbep-4, sin(@as(f128, -0x6.e2369p+16)));
    // try std.testing.expectEqual(-0x1.feb6a3619e803c49fb3b778718d1p-4, sin(@as(f128, -0x6.e2368c006c018p+16)));
    // try std.testing.expectEqual(-0x1.feb6a365968288771150b6f6c51fp-4, sin(@as(f128, -0x6.e2368c006c01cp+16)));
    // try std.testing.expectEqual(-0x1.feb6a361c0bb501b009ef2c1f11ap-4, sin(@as(f128, -0x6.e2368c006c018228p+16)));
    try std.testing.expectEqual(-0x1.4550689b93bbe1406aee26103891p-4, sin(@as(f128, 0x5.6a5008p+64)));
    try std.testing.expectEqual(-0x1.ee01db6bc8ef288c92dcf3ee915cp-4, sin(@as(f128, 0x5.6a5p+64)));
    try std.testing.expectEqual(0x6.5ea3351c9d9da321a84877b1bf9cp-4, sin(@as(f128, 0x5.6a5005df4363cp+64)));
    try std.testing.expectEqual(0x2.f0e32ed649b326445c86bd0d5a5ep-4, sin(@as(f128, 0x5.6a5005df43638p+64)));
    try std.testing.expectEqual(-0xa.8640e82e7924ec0007c75179739p-4, sin(@as(f128, 0x5.6a5005df43638338p+64)));
    try std.testing.expectEqual(0xd.7457bd2255e689f0662a7ba85488p-4, sin(@as(f128, 0x5.6a5005df4363833p+64)));
    // try std.testing.expectEqual(-0xf.fdc3052396dd47c564b32734cc2p-8, sin(@as(f128, 0x5.6a5005df4363833413fa44f74ae8p+64)));
    // try std.testing.expectEqual(-0xf.fdc305227f69439ae6b7b4254f78p-8, sin(@as(f128, 0x5.6a5005df4363833413fa44f74cp+64)));
    // try std.testing.expectEqual(-0xf.fdc305247e694b390edb67a7fa7p-8, sin(@as(f128, 0x5.6a5005df4363833413fa44f74ap+64)));
    try std.testing.expectEqual(-0x8.599b32844aba906cee446be04998p-4, sin(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.452fc98b34e96b61139b09a7c84ap-8, sin(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.dfd9d4b6d0e5f7b9650cab0f5438p-4, sin(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0xf.3b0b11ed85b7fe43d110251580b8p-4, sin(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0xe.f1a3e1dc468a921dddb4e37fbe6p-4, sin(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x8.599b32844aba906cee446be04998p-4, sin(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.452fc98b34e96b61139b09a7c84ap-8, sin(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xf.dfd9d4b6d0e5f7b9650cab0f5438p-4, sin(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0xf.3b0b11ed85b7fe43d110251580b8p-4, sin(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0xe.f1a3e1dc468a921dddb4e37fbe6p-4, sin(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4p-128, sin(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, sin(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, sin(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, sin(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, sin(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, sin(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, sin(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, sin(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, sin(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, sin(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, sin(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, sin(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, sin(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, sin(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, sin(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-152, sin(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sin(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, sin(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, sin(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, sin(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0xf.fa2add3e58948d10238cc27b562p-4, sin(@as(f128, 0x1.8475e6p+0)));
    try std.testing.expectEqual(0xf.fa2adb8953ae26229c919ec8f6cp-4, sin(@as(f128, 0x1.8475e4p+0)));
    try std.testing.expectEqual(0xf.fa2adcf9ea83dbdd053ee455ea7p-4, sin(@as(f128, 0x1.8475e5afd4481p+0)));
}
