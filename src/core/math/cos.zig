const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const sincos_data = @import("sincos_data.zig");
const atnat = @import("atnat.zig");
const usncs = @import("usncs.zig");
const sin = @import("sin.zig");
const branred = @import("branred.zig");
const ldbl128 = @import("ldbl128.zig");
const rem_pio2 = @import("rem_pio2.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;
const ta = @import("builtin").target.cpu.arch;

pub inline fn cos(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return cos(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, cos32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_sinf.c
                    return cos32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_sin.c
                    return cos64(x);
                },
                f80 => return cast(f80, cos128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_sinl.c
                    return cos128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

// Fast cosf implementation.  Worst-case ULP is 0.5607, maximum relative
// error is 0.5303 * 2^-23.  A single-step range reduction is used for
// small values.  Large inputs have their range reduced using fast integer
// arithmetic.
fn cos32(y: f32) f32 {
    var x: f64 = cast(f64, y, .{});
    var p: *const sincos_data.Sincos32 = &sincos_data.__sincos32_table[0];

    if (sincos_data.abstop12(y) < sincos_data.abstop12(sincos_data.pio4)) {
        const x2: f64 = x * x;

        if (sincos_data.abstop12(y) < sincos_data.abstop12(0x1p-12)) {
            @branchHint(.unlikely);
            return 1.0;
        }

        return sincos_data.sin32_poly(x, x2, p, 1);
    } else if (sincos_data.abstop12(y) < sincos_data.abstop12(120.0)) {
        @branchHint(.likely);
        var n: i32 = undefined;
        x = sincos_data.reduce_fast(x, p, &n);

        // Setup the signs for sin and cos.
        const s: f64 = p.sign[@intCast(n & 3)];

        if ((n & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        return sincos_data.sin32_poly(x * s, x * x, p, n ^ 1);
    } else if (sincos_data.abstop12(y) < sincos_data.abstop12(std.math.inf(f32))) {
        const xi: u32 = @bitCast(y);
        const sign: i32 = @intCast(xi >> 31);

        var n: i32 = undefined;
        x = sincos_data.reduce_large(xi, &n);

        // Setup signs for sin and cos - include original sign.
        const s: f64 = p.sign[@intCast((n + sign) & 3)];

        if (((n + sign) & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        return sincos_data.sin32_poly(x * s, x * x, p, n ^ 1);
    } else return (y - y) / (y - y);
}

// Given a number partitioned into X and DX, this function computes the cosine
// of the number by combining the sin and cos of X (as computed by a variation
// of the Taylor series) with the values looked up from the sin/cos table to
// get the result.
pub inline fn do_cos64(x: f64, dx: f64) f64 {
    var dxnew: f64 = dx;
    if (x < 0)
        dxnew = -dxnew;

    const u: [2]i32 = @bitCast(usncs.big + math.abs(x));
    const xnew: f64 = math.abs(x) - (@as(f64, @bitCast(u)) - usncs.big) + dxnew;

    const xx: f64 = xnew * xnew;
    const s: f64 = xnew + xnew * xx * (sincos_data.sn3 + xx * sincos_data.sn5);
    const c: f64 = xx * (sincos_data.cs2 + xx * (sincos_data.cs4 + xx * sincos_data.cs6));
    const k: i32 = u[atnat.LOW_HALF] << 2;
    const sn: f64 = @bitCast(sincos_data.__sincostab[@intCast(k)]);
    const ssn: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 1)]);
    const cs: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 2)]);
    const ccs: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 3)]);
    const cor: f64 = (ccs - s * ssn - cs * c) - sn * s;
    return cs + cor;
}

//******************************************************************/
// An ultimate cos routine. Given an IEEE double machine number x
// it computes the rounded value of cos(x).
//******************************************************************/
fn cos64(x: f64) f64 {
    const u: [2]i32 = @bitCast(x);
    const m: i32 = u[atnat.HIGH_HALF];
    const k: i32 = 0x7fffffff & m;

    // |x|<2^-27 => cos(x)=1
    if (k < 0x3e400000) {
        return 1.0;
    } else if (k < 0x3feb6000) { // 2^-27 < |x| < 0.855469
        // Max ULP is 0.51.
        return do_cos64(x, 0);
    } else if (k < 0x400368fd) { // 0.855469  <|x|<2.426265
        const y: f64 = usncs.hp0 - math.abs(x);
        const a: f64 = y + usncs.hp1;
        const da: f64 = (y - a) + usncs.hp1;
        // Max ULP is 0.501 if xx < 0.01588 or 0.518 otherwise.
        // Range reduction uses 106 bits here which is sufficient.
        return sin.do_sin64(a, da);
    } else if (k < 0x419921fb) { // 2.426265<|x|< 105414350
        var a: f64 = undefined;
        var da: f64 = undefined;
        const n: i32 = sincos_data.reduce_sincos64(x, &a, &da);
        return sincos_data.do_sincos64(a, da, n + 1);
    } else if (k < 0x7ff00000) { // 105414350 <|x| <2^1024
        var a: f64 = undefined;
        var da: f64 = undefined;
        const n: i32 = branred.branred(x, &a, &da);
        return sincos_data.do_sincos64(a, da, n + 1);
    } else {
        return x / x; // |x| > 2^1024
    }
}

pub fn kernel_cos128(x: f128, y: f128) f128 {
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);
    var tix: u32 = @truncate(@as(u64, @bitCast(ix)) >> 32);
    tix &= ~@as(u32, 0x80000000); // tix = |x|'s high 32 bits
    if (tix < 0x3ffc3000) // |x| < 0.1484375
    {
        // Argument is small enough to approximate it by a Chebyshev
        // polynomial of degree 16.
        if (tix < 0x3fc60000) { // |x| < 2^-57
            if (cast(i32, x, .{}) == 1)
                return 1; // generate inexact
        }

        const z: f128 = x * x;
        return 1 + (z * (sincos_data.COS1 + z * (sincos_data.COS2 + z * (sincos_data.COS3 + z * (sincos_data.COS4 +
            z * (sincos_data.COS5 + z * (sincos_data.COS6 + z * (sincos_data.COS7 + z * sincos_data.COS8))))))));
    } else {
        // So that we don't have to use too large polynomial,  we find
        // l and h such that x = l + h,  where fabsl(l) <= 1.0/256 with 83
        // possible values for h.  We look up cosl(h) and sinl(h) in
        // pre-computed tables,  compute cosl(l) and sinl(l) using a
        // Chebyshev polynomial of degree 10(11) and compute
        // cosl(h+l) = cosl(h)cosl(l) - sinl(h)sinl(l).
        var index: u32 = 0x3ffe - (tix >> 16);
        const hix: u32 = (tix + (@as(u32, 0x200) << @as(u5, @intCast(index)))) & (@as(u32, 0xfffffc00) << @as(u5, @intCast(index)));
        var xx: f128 = x;
        var yy: f128 = y;
        if (std.math.signbit(x)) {
            xx = -x;
            yy = -y;
        }
        switch (index) {
            0 => index = ((45 << 10) + hix - 0x3ffe0000) >> 8,
            1 => index = ((13 << 11) + hix - 0x3ffd0000) >> 9,
            else => index = (hix - 0x3ffc3000) >> 10,
        }

        var h: f128 = undefined;
        ldbl128.setWords(&h, cast(u64, hix, .{}) << 32, @as(u64, 0));
        const l: f128 = yy - (h - xx);
        const z: f128 = l * l;
        const sin_l: f128 = l * (1 + z * (sincos_data.SSIN1 + z * (sincos_data.SSIN2 + z * (sincos_data.SSIN3 + z * (sincos_data.SSIN4 + z * sincos_data.SSIN5)))));
        const cos_l_m1: f128 = z * (sincos_data.SCOS1 + z * (sincos_data.SCOS2 + z * (sincos_data.SCOS3 + z * (sincos_data.SCOS4 + z * sincos_data.SCOS5))));
        return usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_LO] - (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] * sin_l - usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] * cos_l_m1));
    }
}

fn cos128(x: f128) f128 {
    // High word of x.
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);

    // |x| ~< pi/4
    ix &= 0x7fffffffffffffff;
    if (ix <= 0x3ffe921fb54442d1) {
        return kernel_cos128(x, 0);
    } else if (ix >= 0x7fff000000000000) { // cos(Inf or NaN) is NaN
        return x - x;
    } else { // argument reduction needed
        var y: [2]f128 = undefined;
        const n = rem_pio2.rem_pio2_128(x, &y);
        switch (n & 3) {
            0 => return kernel_cos128(y[0], y[1]),
            1 => return -sin.kernel_sin128(y[0], y[1], 1),
            2 => return -kernel_cos128(y[0], y[1]),
            else => return sin.kernel_sin128(y[0], y[1], 1),
        }
    }
}

test cos {
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x7.fffff8p-4, cos(@as(f32, 0x1.0c1524p+0)));
    try std.testing.expectEqual(0x8.00001p-4, cos(@as(f32, 0x1.0c1522p+0)));
    try std.testing.expectEqual(-0x8.00001p-4, cos(@as(f32, 0x2.182a48p+0)));
    try std.testing.expectEqual(-0x7.ffffd8p-4, cos(@as(f32, 0x2.182a44p+0)));
    try std.testing.expectEqual(-0xb.bbd2ep-28, cos(@as(f32, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0x1.4442d2p-24, cos(@as(f32, 0x1.921fb4p+0)));
    try std.testing.expectEqual(0xb.b4ff6p-4, cos(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0xf.fb702p-4, cos(@as(f32, 0x2p+64)));
    try std.testing.expectEqual(0xf.fb702p-4, cos(@as(f32, -0x2p+64)));
    try std.testing.expectEqual(0xb.201e7p-4, cos(@as(f32, 0xc.d4967p-4)));
    try std.testing.expectEqual(0xb.201e8p-4, cos(@as(f32, 0xc.d4966p-4)));
    try std.testing.expectEqual(0x2.8f3168p-20, cos(@as(f32, 0xa.217bap+12)));
    try std.testing.expectEqual(0xf.431ddp-4, cos(@as(f32, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(0xa.dd6f7p-4, cos(@as(f32, 0x2.1e19ep+72)));
    try std.testing.expectEqual(0xd.a5f96p-4, cos(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xd.a5f96p-4, cos(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xe.d0668p-4, cos(@as(f32, 0x1p+120)));
    try std.testing.expectEqual(0xc.82b8fp-4, cos(@as(f32, 0x8p+124)));
    try std.testing.expectEqual(0xf.fb2ap-4, cos(@as(f32, 0xf.ffffcp+124)));
    try std.testing.expectEqual(0xd.a5f96p-4, cos(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xd.e3b89p-4, cos(@as(f32, 0x4p+48)));
    try std.testing.expectEqual(-0x2.a62ba8p-4, cos(@as(f32, 0x1p+28)));
    try std.testing.expectEqual(0x8.a513fp-4, cos(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a514p-4, cos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513fp-4, cos(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a514p-4, cos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513fp-4, cos(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a514p-4, cos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513dp-4, cos(@as(f32, 0x1.000004p+0)));
    try std.testing.expectEqual(0x8.a513fp-4, cos(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a513bp-4, cos(@as(f32, 0x1.000006p+0)));
    try std.testing.expectEqual(0x8.a513dp-4, cos(@as(f32, 0x1.000004p+0)));
    try std.testing.expectEqual(-0xf.74fbdp-4, cos(@as(f32, 0x1.200146p+32)));
    try std.testing.expectEqual(0xf.bc96dp-4, cos(@as(f32, 0x1.200144p+32)));
    try std.testing.expectEqual(0x8.a514p-4, cos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x6.a88998p-4, cos(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(-0xf.d7026p-4, cos(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(-0xa.7553p-4, cos(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(0x4.89e16p-4, cos(@as(f32, 0x5p+0)));
    try std.testing.expectEqual(0xf.5cdb8p-4, cos(@as(f32, 0x6p+0)));
    try std.testing.expectEqual(0xc.0ffbdp-4, cos(@as(f32, 0x7p+0)));
    try std.testing.expectEqual(-0x2.53f7d8p-4, cos(@as(f32, 0x8p+0)));
    try std.testing.expectEqual(-0xe.93fd5p-4, cos(@as(f32, 0x9p+0)));
    try std.testing.expectEqual(-0xd.6cd64p-4, cos(@as(f32, 0xap+0)));
    try std.testing.expectEqual(0xf.fe001p-4, cos(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0xf.ffff8p-4, cos(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0xd.a5f96p-4, cos(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xd.a5f96p-4, cos(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0xd.a5f96p-4, cos(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xa.07bd4p-4, cos(@as(f32, 0xe.9f1e6p+112)));
    try std.testing.expectEqual(0xe.26f8bp-4, cos(@as(f32, 0xe.9f1e5p+112)));
    try std.testing.expectEqual(-0xf.dfe9p-4, cos(@as(f32, 0x4.7857dp+68)));
    try std.testing.expectEqual(0x8.7e0eap-4, cos(@as(f32, -0x1.02e34cp+0)));
    try std.testing.expectEqual(-0xd.b7f53p-4, cos(@as(f32, 0xf.f0274p+4)));
    try std.testing.expectEqual(-0xf.dfe6fp-4, cos(@as(f32, 0x3.042d88p+0)));
    try std.testing.expectEqual(0xd.a8263p-8, cos(@as(f32, 0x1.8475e6p+0)));
    try std.testing.expectEqual(0xd.a8283p-8, cos(@as(f32, 0x1.8475e4p+0)));
    try std.testing.expectEqual(-0xc.bbbd3p-24, cos(@as(f32, 0x1.921fc2p+0)));
    try std.testing.expectEqual(-0xa.bbbd3p-24, cos(@as(f32, 0x1.921fcp+0)));
    try std.testing.expectEqual(-0xa.bbbd3p-24, cos(@as(f32, -0x1.921fcp+0)));
    try std.testing.expectEqual(-0xc.bbbd3p-24, cos(@as(f32, -0x1.921fc2p+0)));
    // try std.testing.expectEqual(0xf.ffffdp-4, cos(@as(f32, 0x2.3c6ef4p-12)));
    try std.testing.expectEqual(-0xc.04e3dp-8, cos(@as(f32, 0xe.6672ep+40)));
    // try std.testing.expectEqual(0x4.92b518p-4, cos(@as(f32, 0xe.6672dp+40)));

    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x7.fffff939bdd18p-4, cos(@as(f64, 0x1.0c1524p+0)));
    try std.testing.expectEqual(0x8.000014f038b18p-4, cos(@as(f64, 0x1.0c1522p+0)));
    try std.testing.expectEqual(0x7.ffffffffffff8p-4, cos(@as(f64, 0x1.0c152382d7366p+0)));
    try std.testing.expectEqual(0x8.0000000000008p-4, cos(@as(f64, 0x1.0c152382d7365p+0)));
    try std.testing.expectEqual(-0x8.00000d8c84578p-4, cos(@as(f64, 0x2.182a48p+0)));
    try std.testing.expectEqual(-0x7.ffffd61f8e65cp-4, cos(@as(f64, 0x2.182a44p+0)));
    try std.testing.expectEqual(-0x8.000000000001p-4, cos(@as(f64, 0x2.182a4705ae6ccp+0)));
    try std.testing.expectEqual(-0x7.ffffffffffffp-4, cos(@as(f64, 0x2.182a4705ae6cap+0)));
    try std.testing.expectEqual(-0xb.bbd2e7b96766p-28, cos(@as(f64, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0x1.4442d18469893p-24, cos(@as(f64, 0x1.921fb4p+0)));
    // try std.testing.expectEqual(-0xb.9676733ae8fe8p-56, cos(@as(f64, 0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(0x4.69898cc51701cp-56, cos(@as(f64, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(0xb.b4ff632a908f8p-4, cos(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0xf.fb701e22987f8p-4, cos(@as(f64, 0x2p+64)));
    try std.testing.expectEqual(0xf.fb701e22987f8p-4, cos(@as(f64, -0x2p+64)));
    try std.testing.expectEqual(0xb.201e77869a468p-4, cos(@as(f64, 0xc.d4967p-4)));
    try std.testing.expectEqual(0xb.201e830650418p-4, cos(@as(f64, 0xc.d4966p-4)));
    try std.testing.expectEqual(0xb.201e794508848p-4, cos(@as(f64, 0xc.d4966d92d171p-4)));
    try std.testing.expectEqual(0xb.201e79450885p-4, cos(@as(f64, 0xc.d4966d92d1708p-4)));
    try std.testing.expectEqual(0x2.8f31660ce5e42p-20, cos(@as(f64, 0xa.217bap+12)));
    try std.testing.expectEqual(0xf.431dd7a36cf38p-4, cos(@as(f64, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(0xa.dd6f6bacd2068p-4, cos(@as(f64, 0x2.1e19ep+72)));
    try std.testing.expectEqual(0x8.5f167780e47ap-4, cos(@as(f64, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, cos(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xd.38cf9361195f8p-4, cos(@as(f64, 0x8p+1020)));
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, cos(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba8p-4, cos(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xe.d06685b36c67p-4, cos(@as(f64, 0x1p+120)));
    try std.testing.expectEqual(0xc.82b8ec98b5e6p-4, cos(@as(f64, 0x8p+124)));
    try std.testing.expectEqual(0xf.fb2a030c5ae2p-4, cos(@as(f64, 0xf.ffffcp+124)));
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, cos(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xd.e3b88804f0058p-4, cos(@as(f64, 0x4p+48)));
    try std.testing.expectEqual(-0x2.a62ba8824e5bcp-4, cos(@as(f64, 0x1p+28)));
    try std.testing.expectEqual(0x8.a513eced2ea58p-4, cos(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8346p-4, cos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513fcf3a90fp-4, cos(@as(f64, 0x1.000000cf4a2a2p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea58p-4, cos(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8346p-4, cos(@as(f64, 0x1p+0)));
    // try std.testing.expectEqual(0x8.a513f9cde04e8p-4, cos(@as(f64, 0x1.0000010b239a9p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea58p-4, cos(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8346p-4, cos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513f53385c6p-4, cos(@as(f64, 0x1.00000162a932bp+0)));
    try std.testing.expectEqual(0x8.a513d1ffd9e28p-4, cos(@as(f64, 0x1.000004p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea58p-4, cos(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a513e1c29117p-4, cos(@as(f64, 0x1.000002d452a1p+0)));
    try std.testing.expectEqual(0x8.a513b71284fdp-4, cos(@as(f64, 0x1.000006p+0)));
    try std.testing.expectEqual(0x8.a513d1ffd9e28p-4, cos(@as(f64, 0x1.000004p+0)));
    try std.testing.expectEqual(0x8.a513ba9f703dp-4, cos(@as(f64, 0x1.000005bc7d86dp+0)));
    try std.testing.expectEqual(-0xf.74fbd5498fe5p-4, cos(@as(f64, 0x1.200146p+32)));
    try std.testing.expectEqual(0xf.bc96ca2c658a8p-4, cos(@as(f64, 0x1.200144p+32)));
    try std.testing.expectEqual(-0x6.568e7ed3dffccp-4, cos(@as(f64, 0x1.200145a975ce6p+32)));
    try std.testing.expectEqual(0x8.a51407da8346p-4, cos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x6.a88995d4dc814p-4, cos(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(-0xf.d7025f42f2e9p-4, cos(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(-0xa.7553036d92608p-4, cos(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(0x4.89e15c1ad2b64p-4, cos(@as(f64, 0x5p+0)));
    try std.testing.expectEqual(0xf.5cdb84bc117a8p-4, cos(@as(f64, 0x6p+0)));
    try std.testing.expectEqual(0xc.0ffbcf6c900b8p-4, cos(@as(f64, 0x7p+0)));
    try std.testing.expectEqual(-0x2.53f7d7ec65f28p-4, cos(@as(f64, 0x8p+0)));
    try std.testing.expectEqual(-0xe.93fd53530cb58p-4, cos(@as(f64, 0x9p+0)));
    try std.testing.expectEqual(-0xd.6cd64486358f8p-4, cos(@as(f64, 0xap+0)));
    try std.testing.expectEqual(0xf.fe000aaa93e98p-4, cos(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0xf.ffff800000aa8p-4, cos(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0xf.ffffffep-4, cos(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffff8p-4, cos(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffffffffffep-4, cos(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, cos(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba8p-4, cos(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, cos(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba8p-4, cos(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, cos(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.febbf9949ecc1p-4, cos(@as(f64, -0x3.3de320f6be87ep+1020)));
    try std.testing.expectEqual(-0xa.07bd3ab53ab98p-4, cos(@as(f64, 0xe.9f1e6p+112)));
    try std.testing.expectEqual(0xe.26f8af8333f9p-4, cos(@as(f64, 0xe.9f1e5p+112)));
    try std.testing.expectEqual(0x1.ff01000c9ae73p-4, cos(@as(f64, 0xe.9f1e5bc3bb88p+112)));
    try std.testing.expectEqual(-0xf.dfe902135fc2p-4, cos(@as(f64, 0x4.7857dp+68)));
    try std.testing.expectEqual(0x8.7e0ea4db2f488p-4, cos(@as(f64, -0x1.02e34cp+0)));
    try std.testing.expectEqual(-0xd.b7f5359babdb8p-4, cos(@as(f64, 0xf.f0274p+4)));
    try std.testing.expectEqual(-0xf.dfe6f2169e25p-4, cos(@as(f64, 0x3.042d88p+0)));
    try std.testing.expectEqual(0xd.a8263394be6dp-8, cos(@as(f64, 0x1.8475e6p+0)));
    try std.testing.expectEqual(0xd.a82832da19f98p-8, cos(@as(f64, 0x1.8475e4p+0)));
    try std.testing.expectEqual(0xd.a82683a33cbe8p-8, cos(@as(f64, 0x1.8475e5afd4481p+0)));
    try std.testing.expectEqual(-0xc.bbbd2e7b951e8p-24, cos(@as(f64, 0x1.921fc2p+0)));
    try std.testing.expectEqual(-0xa.bbbd2e7b95a88p-24, cos(@as(f64, 0x1.921fcp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d8b95a5p-24, cos(@as(f64, 0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7b95a5p-24, cos(@as(f64, 0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(-0xa.bbbd2e7b95a88p-24, cos(@as(f64, -0x1.921fcp+0)));
    try std.testing.expectEqual(-0xc.bbbd2e7b951e8p-24, cos(@as(f64, -0x1.921fc2p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7b95a5p-24, cos(@as(f64, -0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d8b95a5p-24, cos(@as(f64, -0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(0xf.ffffd7fffffdp-4, cos(@as(f64, 0x2.3c6ef4p-12)));
    try std.testing.expectEqual(-0xc.04e3d7b33316p-8, cos(@as(f64, 0xe.6672ep+40)));
    try std.testing.expectEqual(0x4.92b51be9ed23p-4, cos(@as(f64, 0xe.6672dp+40)));
    try std.testing.expectEqual(0x1.fd4fd52bd426ap-4, cos(@as(f64, 0xe.6672d458b05fp+40)));
    try std.testing.expectEqual(0x2.053fb048fe646p-4, cos(@as(f64, 0xe.6672d458b05e8p+40)));

    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x7.fffff939bdd18038p-4, cos(@as(f80, 0x1.0c1524p+0)));
    try std.testing.expectEqual(0x8.000014f038b1ab1p-4, cos(@as(f80, 0x1.0c1522p+0)));
    try std.testing.expectEqual(0x7.ffffffffffff94f8p-4, cos(@as(f80, 0x1.0c152382d7366p+0)));
    try std.testing.expectEqual(0x8.00000000000072bp-4, cos(@as(f80, 0x1.0c152382d7365p+0)));
    try std.testing.expectEqual(0x7.ffffffffffffffe8p-4, cos(@as(f80, 0x1.0c152382d7365848p+0)));
    try std.testing.expectEqual(0x8p-4, cos(@as(f80, 0x1.0c152382d7365846p+0)));
    try std.testing.expectEqual(-0x8.00000d8c845743p-4, cos(@as(f80, 0x2.182a48p+0)));
    try std.testing.expectEqual(-0x7.ffffd61f8e65dc98p-4, cos(@as(f80, 0x2.182a44p+0)));
    try std.testing.expectEqual(-0x8.000000000000d61p-4, cos(@as(f80, 0x2.182a4705ae6ccp+0)));
    try std.testing.expectEqual(-0x7.ffffffffffff1abp-4, cos(@as(f80, 0x2.182a4705ae6cap+0)));
    try std.testing.expectEqual(-0x8.000000000000003p-4, cos(@as(f80, 0x2.182a4705ae6cb09p+0)));
    try std.testing.expectEqual(-0x7.fffffffffffffff8p-4, cos(@as(f80, 0x2.182a4705ae6cb08cp+0)));
    try std.testing.expectEqual(-0xb.bbd2e7b96766267p-28, cos(@as(f80, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0x1.4442d1846989361p-24, cos(@as(f80, 0x1.921fb4p+0)));
    try std.testing.expectEqual(-0xb.9676733ae8fe47cp-56, cos(@as(f80, 0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(0x4.69898cc51701b838p-56, cos(@as(f80, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(-0x7.6733ae8fe47c65d8p-68, cos(@as(f80, 0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(0x1.898cc51701b839a2p-64, cos(@as(f80, 0x1.921fb54442d18468p+0)));
    try std.testing.expectEqual(0xb.b4ff632a908f73fp-4, cos(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0xf.fb701e22987fbe7p-4, cos(@as(f80, 0x2p+64)));
    try std.testing.expectEqual(0xf.fb701e22987fbe7p-4, cos(@as(f80, -0x2p+64)));
    try std.testing.expectEqual(0xb.201e77869a46ae2p-4, cos(@as(f80, 0xc.d4967p-4)));
    try std.testing.expectEqual(0xb.201e83065041457p-4, cos(@as(f80, 0xc.d4966p-4)));
    try std.testing.expectEqual(0xb.201e7945088464p-4, cos(@as(f80, 0xc.d4966d92d171p-4)));
    try std.testing.expectEqual(0xb.201e79450884cp-4, cos(@as(f80, 0xc.d4966d92d1708p-4)));
    try std.testing.expectEqual(0xb.201e79450884be2p-4, cos(@as(f80, 0xc.d4966d92d17082ap-4)));
    try std.testing.expectEqual(0xb.201e79450884be3p-4, cos(@as(f80, 0xc.d4966d92d170829p-4)));
    try std.testing.expectEqual(0x2.8f31660ce5e42c04p-20, cos(@as(f80, 0xa.217bap+12)));
    try std.testing.expectEqual(0xf.431dd7a36cf37dep-4, cos(@as(f80, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(0xa.dd6f6bacd20654cp-4, cos(@as(f80, 0x2.1e19ep+72)));
    try std.testing.expectEqual(0x8.5f167780e479c9ap-4, cos(@as(f80, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, cos(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xd.38cf9361195f50bp-4, cos(@as(f80, 0x8p+1020)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, cos(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba9ep-4, cos(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xe.bcc2fc82ae39ecp-4, cos(@as(f80, 0x8p+16380)));
    try std.testing.expectEqual(-0xe.d06685b36c66c4dp-4, cos(@as(f80, 0x1p+120)));
    try std.testing.expectEqual(0xc.82b8ec98b5e62fdp-4, cos(@as(f80, 0x8p+124)));
    try std.testing.expectEqual(0xf.fb2a030c5ae20bep-4, cos(@as(f80, 0xf.ffffcp+124)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, cos(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xd.e3b88804f00552dp-4, cos(@as(f80, 0x4p+48)));
    try std.testing.expectEqual(-0x2.a62ba8824e5bcb08p-4, cos(@as(f80, 0x1p+28)));
    try std.testing.expectEqual(0x8.a513eced2ea575ep-4, cos(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8345c92p-4, cos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513fcf3a90ecp-4, cos(@as(f80, 0x1.000000cf4a2a2p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea575ep-4, cos(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8345c92p-4, cos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513f9cde04e4p-4, cos(@as(f80, 0x1.0000010b239a9p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea575ep-4, cos(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8345c92p-4, cos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513f53385c5cp-4, cos(@as(f80, 0x1.00000162a932bp+0)));
    try std.testing.expectEqual(0x8.a513d1ffd9e28e6p-4, cos(@as(f80, 0x1.000004p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea575ep-4, cos(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a513e1c29116cp-4, cos(@as(f80, 0x1.000002d452a1p+0)));
    try std.testing.expectEqual(0x8.a513b71284fd129p-4, cos(@as(f80, 0x1.000006p+0)));
    try std.testing.expectEqual(0x8.a513d1ffd9e28e6p-4, cos(@as(f80, 0x1.000004p+0)));
    try std.testing.expectEqual(0x8.a513ba9f703d4p-4, cos(@as(f80, 0x1.000005bc7d86dp+0)));
    try std.testing.expectEqual(-0xf.74fbd5498fe4c0cp-4, cos(@as(f80, 0x1.200146p+32)));
    try std.testing.expectEqual(0xf.bc96ca2c658abf6p-4, cos(@as(f80, 0x1.200144p+32)));
    try std.testing.expectEqual(-0x6.568e7ed3dffcdfep-4, cos(@as(f80, 0x1.200145a975ce6p+32)));
    try std.testing.expectEqual(0x8.a51407da8345c92p-4, cos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x6.a88995d4dc81291p-4, cos(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(-0xf.d7025f42f2e9308p-4, cos(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(-0xa.7553036d9260623p-4, cos(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(0x4.89e15c1ad2b654f8p-4, cos(@as(f80, 0x5p+0)));
    try std.testing.expectEqual(0xf.5cdb84bc117abd7p-4, cos(@as(f80, 0x6p+0)));
    try std.testing.expectEqual(0xc.0ffbcf6c900babp-4, cos(@as(f80, 0x7p+0)));
    try std.testing.expectEqual(-0x2.53f7d7ec65f271ecp-4, cos(@as(f80, 0x8p+0)));
    try std.testing.expectEqual(-0xe.93fd53530cb5b82p-4, cos(@as(f80, 0x9p+0)));
    try std.testing.expectEqual(-0xd.6cd64486358f905p-4, cos(@as(f80, 0xap+0)));
    try std.testing.expectEqual(0xf.fe000aaa93e9589p-4, cos(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0xf.ffff800000aaaabp-4, cos(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0xf.ffffffe00000001p-4, cos(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffff8p-4, cos(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffffffffffep-4, cos(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0xf.ffffffffffffff8p-4, cos(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, cos(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba9ep-4, cos(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.002ef4018753d50cp-4, cos(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, cos(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba9ep-4, cos(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.002ef4018753d50cp-4, cos(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, cos(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.febbf9949ecc1336p-4, cos(@as(f80, -0x3.3de320f6be87ep+1020)));
    try std.testing.expectEqual(-0xa.07bd3ab53ab9711p-4, cos(@as(f80, 0xe.9f1e6p+112)));
    try std.testing.expectEqual(0xe.26f8af8333f9271p-4, cos(@as(f80, 0xe.9f1e5p+112)));
    try std.testing.expectEqual(0x1.ff01000c9ae7363p-4, cos(@as(f80, 0xe.9f1e5bc3bb88p+112)));
    try std.testing.expectEqual(-0xf.dfe902135fc1c18p-4, cos(@as(f80, 0x4.7857dp+68)));
    try std.testing.expectEqual(0x8.7e0ea4db2f48867p-4, cos(@as(f80, -0x1.02e34cp+0)));
    try std.testing.expectEqual(-0xd.b7f5359babdb66cp-4, cos(@as(f80, 0xf.f0274p+4)));
    try std.testing.expectEqual(-0xf.dfe6f2169e24f27p-4, cos(@as(f80, 0x3.042d88p+0)));
    try std.testing.expectEqual(0xd.a8263394be6d0e6p-8, cos(@as(f80, 0x1.8475e6p+0)));
    try std.testing.expectEqual(0xd.a82832da19f9892p-8, cos(@as(f80, 0x1.8475e4p+0)));
    try std.testing.expectEqual(0xd.a82683a33cbecp-8, cos(@as(f80, 0x1.8475e5afd4481p+0)));
    try std.testing.expectEqual(-0xc.bbbd2e7b951e5b2p-24, cos(@as(f80, 0x1.921fc2p+0)));
    try std.testing.expectEqual(-0xa.bbbd2e7b95a85c6p-24, cos(@as(f80, 0x1.921fcp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d8b95a502fp-24, cos(@as(f80, 0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7b95a502fp-24, cos(@as(f80, 0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e89a502fp-24, cos(@as(f80, 0x1.921fc00ece4f02f4p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e87a502fp-24, cos(@as(f80, 0x1.921fc00ece4f02f2p+0)));
    try std.testing.expectEqual(-0xa.bbbd2e7b95a85c6p-24, cos(@as(f80, -0x1.921fcp+0)));
    try std.testing.expectEqual(-0xc.bbbd2e7b951e5b2p-24, cos(@as(f80, -0x1.921fc2p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7b95a502fp-24, cos(@as(f80, -0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d8b95a502fp-24, cos(@as(f80, -0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e87a502fp-24, cos(@as(f80, -0x1.921fc00ece4f02f2p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e89a502fp-24, cos(@as(f80, -0x1.921fc00ece4f02f4p+0)));
    try std.testing.expectEqual(0xf.ffffd7fffffcf5ep-4, cos(@as(f80, 0x2.3c6ef4p-12)));
    try std.testing.expectEqual(-0xc.04e3d7b33315e57p-8, cos(@as(f80, 0xe.6672ep+40)));
    try std.testing.expectEqual(0x4.92b51be9ed22fb58p-4, cos(@as(f80, 0xe.6672dp+40)));
    try std.testing.expectEqual(0x1.fd4fd52bd4269f46p-4, cos(@as(f80, 0xe.6672d458b05fp+40)));
    try std.testing.expectEqual(0x2.053fb048fe6462f8p-4, cos(@as(f80, 0xe.6672d458b05e8p+40)));
    try std.testing.expectEqual(0x1.ff55c3f07675f506p-4, cos(@as(f80, 0xe.6672d458b05edf6p+40)));
    try std.testing.expectEqual(0x1.ff56c1efc85b5ec6p-4, cos(@as(f80, 0xe.6672d458b05edf5p+40)));

    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x7.fffff939bdd18035537d20fef1b4p-4, cos(@as(f128, 0x1.0c1524p+0)));
    try std.testing.expectEqual(0x8.000014f038b1ab0e902f6811916p-4, cos(@as(f128, 0x1.0c1522p+0)));
    try std.testing.expectEqual(0x7.ffffffffffff94f4fdce055d4ed4p-4, cos(@as(f128, 0x1.0c152382d7366p+0)));
    try std.testing.expectEqual(0x8.00000000000072a8d510c7c2a25p-4, cos(@as(f128, 0x1.0c152382d7365p+0)));
    try std.testing.expectEqual(0x7.ffffffffffffffe94026ba25319cp-4, cos(@as(f128, 0x1.0c152382d7365848p+0)));
    try std.testing.expectEqual(0x8.0000000000000004f6a1a27d7e48p-4, cos(@as(f128, 0x1.0c152382d7365846p+0)));
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffffffcp-4, cos(@as(f128, 0x1.0c152382d73658465bb32e0f567bp+0)));
    try std.testing.expectEqual(0x8.0000000000000000000000000008p-4, cos(@as(f128, 0x1.0c152382d73658465bb32e0f567ap+0)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffffffb8p-4, cos(@as(f128, 0x1.0c152382d73658465bb32e0f568p+0)));
    try std.testing.expectEqual(0x8.00000000000000000000000006a8p-4, cos(@as(f128, 0x1.0c152382d73658465bb32e0f56p+0)));
    try std.testing.expectEqual(-0x8.00000d8c84574300c524d8042748p-4, cos(@as(f128, 0x2.182a48p+0)));
    try std.testing.expectEqual(-0x7.ffffd61f8e65dc9a1c1408dd99bp-4, cos(@as(f128, 0x2.182a44p+0)));
    try std.testing.expectEqual(-0x8.000000000000d6160463f5455ccp-4, cos(@as(f128, 0x2.182a4705ae6ccp+0)));
    try std.testing.expectEqual(-0x7.ffffffffffff1aae55de707ab4f4p-4, cos(@as(f128, 0x2.182a4705ae6cap+0)));
    try std.testing.expectEqual(-0x8.000000000000002d7fb28bb59cc8p-4, cos(@as(f128, 0x2.182a4705ae6cb09p+0)));
    try std.testing.expectEqual(-0x7.fffffffffffffff612bcbb050374p-4, cos(@as(f128, 0x2.182a4705ae6cb08cp+0)));
    try std.testing.expectEqual(-0x8.0000000000000000000000000008p-4, cos(@as(f128, 0x2.182a4705ae6cb08cb7665c1eacf6p+0)));
    try std.testing.expectEqual(-0x7.ffffffffffffffffffffffffffe8p-4, cos(@as(f128, 0x2.182a4705ae6cb08cb7665c1eacf4p+0)));
    try std.testing.expectEqual(-0x8.000000000000000000000000009p-4, cos(@as(f128, 0x2.182a4705ae6cb08cb7665c1eadp+0)));
    try std.testing.expectEqual(-0x7.fffffffffffffffffffffffff2b4p-4, cos(@as(f128, 0x2.182a4705ae6cb08cb7665c1eacp+0)));
    try std.testing.expectEqual(-0xb.bbd2e7b96766266f1d18f3ead01p-28, cos(@as(f128, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0x1.4442d18469893610281a0f9b0e8dp-24, cos(@as(f128, 0x1.921fb4p+0)));
    try std.testing.expectEqual(-0xb.9676733ae8fe47c65dadfb63ede8p-56, cos(@as(f128, 0x1.921fb54442d19p+0)));
    // try std.testing.expectEqual(0x4.69898cc51701b839a252049c1108p-56, cos(@as(f128, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(-0x7.6733ae8fe47c65dadfb63eeeb308p-68, cos(@as(f128, 0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(0x1.898cc51701b839a252049c1114dp-64, cos(@as(f128, 0x1.921fb54442d18468p+0)));
    try std.testing.expectEqual(-0xc.65dadfb63eeeb306717fbe882b38p-116, cos(@as(f128, 0x1.921fb54442d18469898cc51701b9p+0)));
    try std.testing.expectEqual(0x3.9a252049c1114cf98e804177d4c8p-116, cos(@as(f128, 0x1.921fb54442d18469898cc51701b8p+0)));
    try std.testing.expectEqual(-0x4.7c65dadfb63eeeb306717fbe882cp-108, cos(@as(f128, 0x1.921fb54442d18469898cc51702p+0)));
    try std.testing.expectEqual(0x3.839a252049c1114cf98e804177d4p-108, cos(@as(f128, 0x1.921fb54442d18469898cc517018p+0)));
    try std.testing.expectEqual(0xb.b4ff632a908f73ec151839cb9d98p-4, cos(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0xf.fb701e22987fbe68852ee2bc897p-4, cos(@as(f128, 0x2p+64)));
    try std.testing.expectEqual(0xf.fb701e22987fbe68852ee2bc897p-4, cos(@as(f128, -0x2p+64)));
    try std.testing.expectEqual(0xb.201e77869a46ae20ce545c5c67p-4, cos(@as(f128, 0xc.d4967p-4)));
    try std.testing.expectEqual(0xb.201e83065041456a084c70f5a128p-4, cos(@as(f128, 0xc.d4966p-4)));
    try std.testing.expectEqual(0xb.201e794508846402500c44b4f8e8p-4, cos(@as(f128, 0xc.d4966d92d171p-4)));
    try std.testing.expectEqual(0xb.201e79450884c00000000000c178p-4, cos(@as(f128, 0xc.d4966d92d1708p-4)));
    try std.testing.expectEqual(0xb.201e79450884be1d0c24406973ap-4, cos(@as(f128, 0xc.d4966d92d17082ap-4)));
    try std.testing.expectEqual(0xb.201e79450884be288bda3ee0dd18p-4, cos(@as(f128, 0xc.d4966d92d170829p-4)));
    try std.testing.expectEqual(0xb.201e79450884be22c53e47ed16d8p-4, cos(@as(f128, 0xc.d4966d92d17082980965c1a663c8p-4)));
    try std.testing.expectEqual(0xb.201e79450884be22c53e47ed16ep-4, cos(@as(f128, 0xc.d4966d92d17082980965c1a663cp-4)));
    try std.testing.expectEqual(0xb.201e79450884be22c53e47ed16bp-4, cos(@as(f128, 0xc.d4966d92d17082980965c1a664p-4)));
    try std.testing.expectEqual(0xb.201e79450884be22c53e47ed199p-4, cos(@as(f128, 0xc.d4966d92d17082980965c1a66p-4)));
    try std.testing.expectEqual(0x2.8f31660ce5e42c0544355e8e3d04p-20, cos(@as(f128, 0xa.217bap+12)));
    try std.testing.expectEqual(0xf.431dd7a36cf37de5c74544f6b438p-4, cos(@as(f128, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(0xa.dd6f6bacd20654c1404f52cde16p-4, cos(@as(f128, 0x2.1e19ep+72)));
    try std.testing.expectEqual(0x8.5f167780e479c9a5c86ffce7615p-4, cos(@as(f128, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, cos(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xd.38cf9361195f50b10fac29dd9038p-4, cos(@as(f128, 0x8p+1020)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, cos(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba9e038d934070f138p-4, cos(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xe.bcc2fc82ae39ebf8da5d687bf36p-4, cos(@as(f128, 0x8p+16380)));
    try std.testing.expectEqual(-0x5.b773d971a848e75c230605526978p-4, cos(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0xe.d06685b36c66c4cf35c11f6519p-4, cos(@as(f128, 0x1p+120)));
    try std.testing.expectEqual(0xc.82b8ec98b5e62fcf0b09fd10eb3p-4, cos(@as(f128, 0x8p+124)));
    try std.testing.expectEqual(0xf.fb2a030c5ae20bdfe29fda198eap-4, cos(@as(f128, 0xf.ffffcp+124)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, cos(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xd.e3b88804f00552d6baba709471d8p-4, cos(@as(f128, 0x4p+48)));
    try std.testing.expectEqual(-0x2.a62ba8824e5bcb065f5f3b8e4f58p-4, cos(@as(f128, 0x1p+28)));
    try std.testing.expectEqual(0x8.a513eced2ea575e738a147c82bd8p-4, cos(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8345c91c2466d9768718p-4, cos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513fcf3a90ec00000037aea619p-4, cos(@as(f128, 0x1.000000cf4a2a2p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea575e738a147c82bd8p-4, cos(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8345c91c2466d9768718p-4, cos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513f9cde04e4000000314b550f8p-4, cos(@as(f128, 0x1.0000010b239a9p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea575e738a147c82bd8p-4, cos(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a51407da8345c91c2466d9768718p-4, cos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x8.a513f53385c5c0000002a6dfa3ep-4, cos(@as(f128, 0x1.00000162a932bp+0)));
    try std.testing.expectEqual(0x8.a513d1ffd9e28e629926fb8f7fcp-4, cos(@as(f128, 0x1.000004p+0)));
    try std.testing.expectEqual(0x8.a513eced2ea575e738a147c82bd8p-4, cos(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x8.a513e1c29116c0000003f8dd14cp-4, cos(@as(f128, 0x1.000002d452a1p+0)));
    try std.testing.expectEqual(0x8.a513b71284fd128eb1ad47d820ep-4, cos(@as(f128, 0x1.000006p+0)));
    try std.testing.expectEqual(0x8.a513d1ffd9e28e629926fb8f7fcp-4, cos(@as(f128, 0x1.000004p+0)));
    try std.testing.expectEqual(0x8.a513ba9f703d3ffffffcb9235418p-4, cos(@as(f128, 0x1.000005bc7d86dp+0)));
    try std.testing.expectEqual(-0xf.74fbd5498fe4c0c71bd9e4ef59e8p-4, cos(@as(f128, 0x1.200146p+32)));
    try std.testing.expectEqual(0xf.bc96ca2c658abf5ace7b886a8fbp-4, cos(@as(f128, 0x1.200144p+32)));
    try std.testing.expectEqual(-0x6.568e7ed3dffcdfe227fd726840e4p-4, cos(@as(f128, 0x1.200145a975ce6p+32)));
    try std.testing.expectEqual(0x8.a51407da8345c91c2466d9768718p-4, cos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x6.a88995d4dc81290ccbe2b2edcac4p-4, cos(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(-0xf.d7025f42f2e9307dff82fdf6a708p-4, cos(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(-0xa.7553036d926062336d0e16e3dd5p-4, cos(@as(f128, 0x4p+0)));
    try std.testing.expectEqual(0x4.89e15c1ad2b654f99f75a35ee5fcp-4, cos(@as(f128, 0x5p+0)));
    try std.testing.expectEqual(0xf.5cdb84bc117abd74f1e700074a98p-4, cos(@as(f128, 0x6p+0)));
    try std.testing.expectEqual(0xc.0ffbcf6c900baafbd68c5a99d55p-4, cos(@as(f128, 0x7p+0)));
    try std.testing.expectEqual(-0x2.53f7d7ec65f271ec91f976afbdcep-4, cos(@as(f128, 0x8p+0)));
    try std.testing.expectEqual(-0xe.93fd53530cb5b8268bb2e8949aa8p-4, cos(@as(f128, 0x9p+0)));
    try std.testing.expectEqual(-0xd.6cd64486358f904f7e2a0b9994e8p-4, cos(@as(f128, 0xap+0)));
    try std.testing.expectEqual(0xf.fe000aaa93e9589576da4ec94948p-4, cos(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(0xf.ffff800000aaaaaa4fa4fa69a698p-4, cos(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0xf.ffffffe00000000aaaaaaaa93e9p-4, cos(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffff80000000000aaaaaaa8p-4, cos(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffffffffffe0000000000000aa8p-4, cos(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0xf.ffffffffffffff8p-4, cos(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0xf.ffffffffffffffffep-4, cos(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0xf.fffffffffffffffffff8p-4, cos(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffep-4, cos(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffff8p-4, cos(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffffep-4, cos(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, cos(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba9e038d934070f138p-4, cos(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.002ef4018753d50b7a7f6bc3f5bap-4, cos(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x4.e6dc95fb529bc365f472e4fbc1f8p-4, cos(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x5.b773d971a848e75c230605526978p-4, cos(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, cos(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fff31767d5ba9e038d934070f138p-4, cos(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.002ef4018753d50b7a7f6bc3f5bap-4, cos(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x4.e6dc95fb529bc365f472e4fbc1f8p-4, cos(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x5.b773d971a848e75c230605526978p-4, cos(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, cos(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, cos(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.febbf9949ecc133623bb8c8c5a27p-4, cos(@as(f128, -0x3.3de320f6be87ep+1020)));
    try std.testing.expectEqual(-0xa.07bd3ab53ab9710f3445538de8fp-4, cos(@as(f128, 0xe.9f1e6p+112)));
    try std.testing.expectEqual(0xe.26f8af8333f9270e9c3e9f64f94p-4, cos(@as(f128, 0xe.9f1e5p+112)));
    // try std.testing.expectEqual(0x1.ff01000c9ae73630add558c936b5p-4, cos(@as(f128, 0xe.9f1e5bc3bb88p+112)));
    try std.testing.expectEqual(-0xf.dfe902135fc1c18492e869a3f8a8p-4, cos(@as(f128, 0x4.7857dp+68)));
    try std.testing.expectEqual(0x8.7e0ea4db2f488671c85df7208968p-4, cos(@as(f128, -0x1.02e34cp+0)));
    try std.testing.expectEqual(-0xd.b7f5359babdb66be8d0cd3e293e8p-4, cos(@as(f128, 0xf.f0274p+4)));
    try std.testing.expectEqual(-0xf.dfe6f2169e24f276e8027d91ba9p-4, cos(@as(f128, 0x3.042d88p+0)));
    try std.testing.expectEqual(0xd.a8263394be6d0e58c1c35a8a3bap-8, cos(@as(f128, 0x1.8475e6p+0)));
    try std.testing.expectEqual(0xd.a82832da19f9891d9762fa659ff8p-8, cos(@as(f128, 0x1.8475e4p+0)));
    try std.testing.expectEqual(0xd.a82683a33cbebfffffffa2966878p-8, cos(@as(f128, 0x1.8475e5afd4481p+0)));
    try std.testing.expectEqual(-0xc.bbbd2e7b951e5b1e4cc9f460a12p-24, cos(@as(f128, 0x1.921fc2p+0)));
    // try std.testing.expectEqual(-0xa.bbbd2e7b95a85c638e746a5f4f6p-24, cos(@as(f128, 0x1.921fcp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d8b95a502ede0fd607f394p-24, cos(@as(f128, 0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7b95a502ede4a0f9b11c28p-24, cos(@as(f128, 0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e89a502ede3f507aaa7d8p-24, cos(@as(f128, 0x1.921fc00ece4f02f4p+0)));
    // try std.testing.expectEqual(-0xa.ca8b7d7e87a502ede3f57c1dce1p-24, cos(@as(f128, 0x1.921fc00ece4f02f2p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e881db0d49193eaab43cp-24, cos(@as(f128, 0x1.921fc00ece4f02f278ade6ad9e8ap+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e881db0d49193e9ab43cp-24, cos(@as(f128, 0x1.921fc00ece4f02f278ade6ad9e89p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e881db0d4919460ab43cp-24, cos(@as(f128, 0x1.921fc00ece4f02f278ade6ad9fp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e881db0d49193e0ab43cp-24, cos(@as(f128, 0x1.921fc00ece4f02f278ade6ad9e8p+0)));
    // try std.testing.expectEqual(-0xa.bbbd2e7b95a85c638e746a5f4f6p-24, cos(@as(f128, -0x1.921fcp+0)));
    try std.testing.expectEqual(-0xc.bbbd2e7b951e5b1e4cc9f460a12p-24, cos(@as(f128, -0x1.921fc2p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7b95a502ede4a0f9b11c28p-24, cos(@as(f128, -0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(-0xa.ca8b7d8b95a502ede0fd607f394p-24, cos(@as(f128, -0x1.921fc00ece4f1p+0)));
    // try std.testing.expectEqual(-0xa.ca8b7d7e87a502ede3f57c1dce1p-24, cos(@as(f128, -0x1.921fc00ece4f02f2p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e89a502ede3f507aaa7d8p-24, cos(@as(f128, -0x1.921fc00ece4f02f4p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e881db0d49193e9ab43cp-24, cos(@as(f128, -0x1.921fc00ece4f02f278ade6ad9e89p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e881db0d49193eaab43cp-24, cos(@as(f128, -0x1.921fc00ece4f02f278ade6ad9e8ap+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e881db0d49193e0ab43cp-24, cos(@as(f128, -0x1.921fc00ece4f02f278ade6ad9e8p+0)));
    try std.testing.expectEqual(-0xa.ca8b7d7e881db0d4919460ab43cp-24, cos(@as(f128, -0x1.921fc00ece4f02f278ade6ad9fp+0)));
    try std.testing.expectEqual(0xf.ffffd7fffffcf5e6384f874e68fp-4, cos(@as(f128, 0x2.3c6ef4p-12)));
    try std.testing.expectEqual(-0xc.04e3d7b33315e56d155d1ce9ee18p-8, cos(@as(f128, 0xe.6672ep+40)));
    try std.testing.expectEqual(0x4.92b51be9ed22fb55105031b7163p-4, cos(@as(f128, 0xe.6672dp+40)));
    try std.testing.expectEqual(0x1.fd4fd52bd4269f45bd8a6a988d49p-4, cos(@as(f128, 0xe.6672d458b05fp+40)));
    // try std.testing.expectEqual(0x2.053fb048fe6462f6c5d88f55f69cp-4, cos(@as(f128, 0xe.6672d458b05e8p+40)));
    // try std.testing.expectEqual(0x1.ff55c3f07675f5065a4c74830913p-4, cos(@as(f128, 0xe.6672d458b05edf6p+40)));
    // try std.testing.expectEqual(0x1.ff56c1efc85b5ec6aad454c88ddap-4, cos(@as(f128, 0xe.6672d458b05edf5p+40)));
    // try std.testing.expectEqual(0x1.ff56b710bf1d3d16bf762e86d21p-4, cos(@as(f128, 0xe.6672d458b05edf50af4fab1a42p+40)));
    // try std.testing.expectEqual(0x1.ff56b710bf1d3cf6ff8c0f12aa02p-4, cos(@as(f128, 0xe.6672d458b05edf50af4fab1a44p+40)));
    // try std.testing.expectEqual(0x1.ff56b710bf1d3d367f604dfafa1fp-4, cos(@as(f128, 0xe.6672d458b05edf50af4fab1a4p+40)));
}
