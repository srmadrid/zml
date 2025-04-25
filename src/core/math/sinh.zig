const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const roundeven = @import("roundeven.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn sinh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return sinh(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, sinh32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_sinhf.c
                    return sinh32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_sinh.c
                    return sinh64(x);
                },
                f80 => return cast(f80, sinh128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_sinhl.c
                    return sinh128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn sinh32(x: f32) f32 {
    const c: [4]f64 = .{ 1, 0x1.62e42fef4c4e7p-6, 0x1.ebfd1b232f475p-13, 0x1.c6b19384ecd93p-20 };
    const ch: [7]f64 = .{
        1,                     0x1.62e42fefa39efp-6,  0x1.ebfbdff82c58fp-13,
        0x1.c6b08d702e0edp-20, 0x1.3b2ab6fb92e5ep-27, 0x1.5d886e6d54203p-35,
        0x1.430976b8ce6efp-43,
    };
    const tb: [32]u64 = .{
        0x3fe0000000000000, 0x3fe059b0d3158574, 0x3fe0b5586cf9890f,
        0x3fe11301d0125b51, 0x3fe172b83c7d517b, 0x3fe1d4873168b9aa,
        0x3fe2387a6e756238, 0x3fe29e9df51fdee1, 0x3fe306fe0a31b715,
        0x3fe371a7373aa9cb, 0x3fe3dea64c123422, 0x3fe44e086061892d,
        0x3fe4bfdad5362a27, 0x3fe5342b569d4f82, 0x3fe5ab07dd485429,
        0x3fe6247eb03a5585, 0x3fe6a09e667f3bcd, 0x3fe71f75e8ec5f74,
        0x3fe7a11473eb0187, 0x3fe82589994cce13, 0x3fe8ace5422aa0db,
        0x3fe93737b0cdc5e5, 0x3fe9c49182a3f090, 0x3fea5503b23e255d,
        0x3feae89f995ad3ad, 0x3feb7f76f2fb5e47, 0x3fec199bdd85529c,
        0x3fecb720dcef9069, 0x3fed5818dcfba487, 0x3fedfc97337b9b5f,
        0x3feea4afa2a490da, 0x3fef50765b6e4540,
    };
    const st: struct { uarg: u32, rh: f32, rl: f32 } = .{ .uarg = 0x74250bfe, .rh = 0x1.250bfep-11, .rl = 0x1p-36 };
    const iln2: f64 = 0x1.71547652b82fep+5;
    const z: f64 = cast(f64, x, .{});
    const ux: u32 = @as(u32, @bitCast(x)) << 1;
    if (ux > 0x8565a9f8) { // |x| >~ 89.4
        @branchHint(.unlikely);
        const sgn: f32 = math.copysign(@as(f32, 2), x);
        if (ux >= 0xff000000) {
            if ((ux << 8) != 0)
                return x + x; // nan

            return math.copysign(std.math.inf(f32), x); // +-inf
        }
        return sgn * 0x1.fffffep127;
    }

    if (ux < 0x7c000000) { // |x| < 0.125
        @branchHint(.unlikely);
        if (ux <= 0x74250bfe) { // |x| <= 0x1.250bfep-11
            @branchHint(.unlikely);
            if (ux < 0x66000000) { // |x| < 0x1p-24
                @branchHint(.unlikely);
                return @mulAdd(f32, x, math.abs(x), x);
            }

            if (st.uarg == ux) {
                @branchHint(.unlikely);
                const sgn: f32 = math.copysign(@as(f32, 1), x);
                return sgn * st.rh + sgn * st.rl;
            }

            return (x * 0x1.555556p-3) * (x * x) + x;
        }

        const cp: [4]f64 = .{
            0x1.5555555555555p-3,  0x1.11111111146e1p-7,
            0x1.a01a00930dda6p-13, 0x1.71f92198aa6e9p-19,
        };
        const z2: f64 = z * z;
        const z4: f64 = z2 * z2;
        return cast(f32, z + (z2 * z) * ((cp[0] + z2 * cp[1]) + z4 * (cp[2] + z2 * (cp[3]))), .{});
    }

    const a: f64 = iln2 * z;
    const ia: f64 = roundeven.roundeven_finite(a);
    var h: f64 = a - ia;
    var h2: f64 = h * h;
    const jp: i64 = @bitCast(ia + 0x1.8p52);
    const jm: i64 = -jp;
    const sp: f64 = @bitCast(cast(i64, tb[@intCast(jp & 31)], .{}) + ((jp >> 5) << 52));
    const sm: f64 = @bitCast(cast(i64, tb[@intCast(jm & 31)], .{}) + ((jm >> 5) << 52));
    var te: f64 = c[0] + h2 * c[2];
    var to: f64 = (c[1] + h2 * c[3]);
    const rp: f64 = sp * (te + h * to);
    const rm: f64 = sm * (te - h * to);
    var r: f64 = rp - rm;
    var ub: f32 = cast(f32, r, .{});
    const lb: f32 = cast(f32, r - 1.52e-10 * r, .{});
    if (ub != lb) {
        @branchHint(.unlikely);
        const iln2h: f64 = 0x1.7154765p+5;
        const iln2l: f64 = 0x1.5c17f0bbbe88p-26;
        h = (iln2h * z - ia) + iln2l * z;
        h2 = h * h;
        te = ch[0] + h2 * ch[2] + (h2 * h2) * (ch[4] + h2 * ch[6]);
        to = ch[1] + h2 * (ch[3] + h2 * ch[5]);
        r = sp * (te + h * to) - sm * (te - h * to);
        ub = cast(f32, r, .{});
    }

    return ub;
}

fn sinh64(x: f64) f64 {
    const shuge: f64 = 1.0e307;

    // High word of |x|.
    var jx: i32 = undefined;
    dbl64.getHighWord(&jx, x);
    const ix: i32 = jx & 0x7fffffff;

    // x is INF or NaN
    if (ix >= 0x7ff00000) {
        @branchHint(.unlikely);
        return x + x;
    }

    var h: f64 = 0.5;
    if (jx < 0)
        h = -h;
    // |x| in [0,22], return sign(x)*0.5*(E+E/(E+1)))
    if (ix < 0x40360000) { // |x|<22
        if (ix < 0x3e300000) { // |x|<2**-28
            @branchHint(.unlikely);
            if (math.abs(x) < std.math.floatMin(f64)) {
                const vx: f64 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            if (shuge + x > 1)
                return x;

            // sinh(tiny) = tiny with inexact
        }

        const t: f64 = math.expm1(math.abs(x));
        if (ix < 0x3ff00000)
            return h * (2 * t - t * t / (t + 1));

        return h * (t + t / (t + 1));
    }

    // |x| in [22, log(maxdouble)] return 0.5*exp(|x|)
    if (ix < 0x40862e42)
        return h * math.exp(math.abs(x));

    // |x| in [log(maxdouble), overflowthresold]
    var lx: u32 = undefined;
    dbl64.getLowWord(&lx, x);
    if (ix < 0x408633ce or ((ix == 0x408633ce) and (lx <= 0x8fb9f87d))) {
        const w: f64 = math.exp(0.5 * math.abs(x));
        const t: f64 = h * w;
        return t * w;
    }

    // |x| > overflowthresold, sinh(x) overflow
    return x * shuge;
}

fn sinh128(x: f128) f128 {
    const shuge: f128 = 1.0e4931;
    const ovf_thresh: f128 = 1.1357216553474703894801348310092223067821E4;

    // Words of |x|.
    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const jx: u32 = u.w0;
    const ix: u32 = jx & 0x7fffffff;

    // x is INF or NaN
    if (ix >= 0x7fff0000)
        return x + x;

    var h: f128 = 0.5;
    if ((jx & 0x80000000) != 0)
        h = -h;

    // Absolute value of x.
    u.w0 = ix;

    // |x| in [0,40], return sign(x)*0.5*(E+E/(E+1)))
    if (ix <= 0x40044000) {
        if (ix < 0x3fc60000) { // |x| < 2^-57
            if (math.abs(x) < std.math.floatMin(f128)) {
                const vx: f128 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            if (shuge + x > 1)
                return x; // sinh(tiny) = tiny with inexact
        }
        const t: f128 = math.expm1(@as(f128, @bitCast(u)));
        if (ix < 0x3fff0000)
            return h * (2 * t - t * t / (t + 1));

        return h * (t + t / (t + 1));
    }

    // |x| in [40, log(maxdouble)] return 0.5*exp(|x|)
    if (ix <= 0x400c62e3) // 11356.375
        return h * math.exp(@as(f128, @bitCast(u)));

    // |x| in [log(maxdouble), overflowthreshold]
    // Overflow threshold is log(2 * maxdouble).
    if (@as(f128, @bitCast(u)) <= ovf_thresh) {
        const w: f128 = math.exp(0.5 * @as(f128, @bitCast(u)));
        const t: f128 = h * w;
        return t * w;
    }

    // |x| > overflowthreshold, sinhl(x) overflow
    return x * shuge;
}

test sinh {
    try std.testing.expectEqual(0x0p+0, sinh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0xd.28359p-4, sinh(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x8p-32, sinh(@as(f32, 0x8p-32)));
    try std.testing.expectEqual(0x8.00555p-8, sinh(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(-0x8.00555p-8, sinh(@as(f32, -0x8p-8)));
    try std.testing.expectEqual(0x4.000008p-12, sinh(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(-0x4.000008p-12, sinh(@as(f32, -0x4p-12)));
    try std.testing.expectEqual(0x1p-20, sinh(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(-0x1p-20, sinh(@as(f32, -0x1p-20)));
    try std.testing.expectEqual(0x4p-32, sinh(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(-0x4p-32, sinh(@as(f32, -0x4p-32)));
    try std.testing.expectEqual(0x1p-40, sinh(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(-0x1p-40, sinh(@as(f32, -0x1p-40)));
    try std.testing.expectEqual(0x4p-52, sinh(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(-0x4p-52, sinh(@as(f32, -0x4p-52)));
    try std.testing.expectEqual(0x1p-60, sinh(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(-0x1p-60, sinh(@as(f32, -0x1p-60)));
    try std.testing.expectEqual(0x4p-72, sinh(@as(f32, 0x4p-72)));
    try std.testing.expectEqual(-0x4p-72, sinh(@as(f32, -0x4p-72)));
    try std.testing.expectEqual(0x1p-100, sinh(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, sinh(@as(f32, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sinh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sinh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.ad6b7p+28, sinh(@as(f32, 0x1.6p+4)));
    try std.testing.expectEqual(0x1.226af4p+32, sinh(@as(f32, 0x1.7p+4)));
    try std.testing.expectEqual(0x3.156ff8p+32, sinh(@as(f32, 0x1.8p+4)));
    try std.testing.expectEqual(-0x7.985208p-4, sinh(@as(f32, -0x7.55d7f8p-4)));
    try std.testing.expectEqual(-0x3.fde378p-4, sinh(@as(f32, -0x3.f392f8p-4)));
    try std.testing.expectEqual(0x2.da7cd8p+0, sinh(@as(f32, 0x1.c56446p+0)));
    try std.testing.expectEqual(0x6.ff782p-4, sinh(@as(f32, 0x6.cac628p-4)));
    try std.testing.expectEqual(0x6.ff781p-4, sinh(@as(f32, 0x6.cac62p-4)));
    try std.testing.expectEqual(-0xa.0100dp+4, sinh(@as(f32, -0x5.c4cbp+0)));
    try std.testing.expectEqual(-0xa.01012p+4, sinh(@as(f32, -0x5.c4cb08p+0)));
    try std.testing.expectEqual(-0x1.e33aeep+0, sinh(@as(f32, -0x1.64685p+0)));
    try std.testing.expectEqual(-0x1.e33af2p+0, sinh(@as(f32, -0x1.646852p+0)));
    try std.testing.expectEqual(-0x7.f4861p-4, sinh(@as(f32, -0x7.a8c5fp-4)));
    try std.testing.expectEqual(-0x7.f48618p-4, sinh(@as(f32, -0x7.a8c5f8p-4)));
    try std.testing.expectEqual(0x3.4ff4d8p-4, sinh(@as(f32, 0x3.4a037p-4)));
    try std.testing.expectEqual(-0x3.f5b9acp-4, sinh(@as(f32, -0x3.eba6d8p-4)));
    try std.testing.expectEqual(-0x3.f5b9bp-4, sinh(@as(f32, -0x3.eba6dcp-4)));
    try std.testing.expectEqual(-0x5.1ed4bp+0, sinh(@as(f32, -0x2.55f63p+0)));
    try std.testing.expectEqual(-0x3.d3834cp-4, sinh(@as(f32, -0x3.ca68c8p-4)));
    try std.testing.expectEqual(-0x3.d3835p-4, sinh(@as(f32, -0x3.ca68ccp-4)));
    try std.testing.expectEqual(-0x3.9a7a2p-4, sinh(@as(f32, -0x3.92da04p-4)));
    try std.testing.expectEqual(-0x3.9a7a24p-4, sinh(@as(f32, -0x3.92da08p-4)));
    try std.testing.expectEqual(-0x3.4415b8p-4, sinh(@as(f32, -0x3.3e629p-4)));
    try std.testing.expectEqual(-0x3.4415bcp-4, sinh(@as(f32, -0x3.3e6294p-4)));
    try std.testing.expectEqual(0x7.b341ep-4, sinh(@as(f32, 0x7.6e25ap-4)));
    try std.testing.expectEqual(0x7.b341d8p-4, sinh(@as(f32, 0x7.6e2598p-4)));
    try std.testing.expectEqual(0x3.e05638p-4, sinh(@as(f32, 0x3.d6e088p-4)));
    try std.testing.expectEqual(-0x7.ad0e4p-4, sinh(@as(f32, -0x7.688eap-4)));
    try std.testing.expectEqual(-0xf.a9e6ep-4, sinh(@as(f32, -0xd.dce79p-4)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x8.a3127p+4)));
    try std.testing.expectEqual(0x1.c0709p-12, sinh(@as(f32, 0x1.c0709p-12)));
    try std.testing.expectEqual(0xc.835a6p-4, sinh(@as(f32, 0xb.7f67dp-4)));
    try std.testing.expectEqual(0xc.835a4p-4, sinh(@as(f32, 0xb.7f67cp-4)));
    try std.testing.expectEqual(-0x1.960d6ep+0, sinh(@as(f32, -0x1.3dda8ap+0)));
    try std.testing.expectEqual(-0x6.11995p-4, sinh(@as(f32, -0x5.ee9218p-4)));
    try std.testing.expectEqual(-0x2.c176ap+0, sinh(@as(f32, -0x1.bcfc98p+0)));
    try std.testing.expectEqual(-0x6.cc3dd8p-4, sinh(@as(f32, -0x6.9bbb68p-4)));
    try std.testing.expectEqual(-0x6.cc3dep-4, sinh(@as(f32, -0x6.9bbb7p-4)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-std.math.inf(f32), sinh(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x4p-128, sinh(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-128, sinh(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0xf.fffecp+124, sinh(@as(f32, 0x5.96a7ep+4)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x5.96a7e8p+4)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), sinh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(-std.math.inf(f32), sinh(@as(f32, -0x2.c678c4p+8)));
    try std.testing.expectEqual(-std.math.inf(f32), sinh(@as(f32, -0x2.c678c8p+8)));

    try std.testing.expectEqual(0x0p+0, sinh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0xd.283596e9e348p-4, sinh(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x8p-32, sinh(@as(f64, 0x8p-32)));
    try std.testing.expectEqual(0x8.0055566668068p-8, sinh(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(-0x8.0055566668068p-8, sinh(@as(f64, -0x8p-8)));
    // try std.testing.expectEqual(0x4.00000aaaaab34p-12, sinh(@as(f64, 0x4p-12)));
    // try std.testing.expectEqual(-0x4.00000aaaaab34p-12, sinh(@as(f64, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000000002abp-20, sinh(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(-0x1.00000000002abp-20, sinh(@as(f64, -0x1p-20)));
    try std.testing.expectEqual(0x4p-32, sinh(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(-0x4p-32, sinh(@as(f64, -0x4p-32)));
    try std.testing.expectEqual(0x1p-40, sinh(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(-0x1p-40, sinh(@as(f64, -0x1p-40)));
    try std.testing.expectEqual(0x4p-52, sinh(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(-0x4p-52, sinh(@as(f64, -0x4p-52)));
    try std.testing.expectEqual(0x1p-60, sinh(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(-0x1p-60, sinh(@as(f64, -0x1p-60)));
    try std.testing.expectEqual(0x4p-72, sinh(@as(f64, 0x4p-72)));
    try std.testing.expectEqual(-0x4p-72, sinh(@as(f64, -0x4p-72)));
    try std.testing.expectEqual(0x1p-100, sinh(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, sinh(@as(f64, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sinh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-1000, sinh(@as(f64, 0x1p-1000)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-1000, sinh(@as(f64, -0x1p-1000)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sinh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, sinh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sinh(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x6.ad6b6e710d8p+28, sinh(@as(f64, 0x1.6p+4)));
    try std.testing.expectEqual(0x1.226af33b1fdc1p+32, sinh(@as(f64, 0x1.7p+4)));
    try std.testing.expectEqual(0x3.156ff6a8ebf6ep+32, sinh(@as(f64, 0x1.8p+4)));
    // try std.testing.expectEqual(-0x7.9852071dfda98p-4, sinh(@as(f64, -0x7.55d7f8p-4)));
    // try std.testing.expectEqual(-0x3.fde378210a8f8p-4, sinh(@as(f64, -0x3.f392f8p-4)));
    try std.testing.expectEqual(0x2.da7cd9753b47cp+0, sinh(@as(f64, 0x1.c56446p+0)));
    // try std.testing.expectEqual(0x6.ff781ca6e6dccp-4, sinh(@as(f64, 0x6.cac628p-4)));
    try std.testing.expectEqual(0x6.ff7813eb9593cp-4, sinh(@as(f64, 0x6.cac62p-4)));
    // try std.testing.expectEqual(0x6.ff78170306f9cp-4, sinh(@as(f64, 0x6.cac622d51eebcp-4)));
    try std.testing.expectEqual(-0xa.0100cebf41c8p+4, sinh(@as(f64, -0x5.c4cbp+0)));
    try std.testing.expectEqual(-0xa.01011ec7afdap+4, sinh(@as(f64, -0x5.c4cb08p+0)));
    try std.testing.expectEqual(-0xa.0100e4f7b10f8p+4, sinh(@as(f64, -0x5.c4cb02389c094p+0)));
    // try std.testing.expectEqual(-0x1.e33aed09484p+0, sinh(@as(f64, -0x1.64685p+0)));
    // try std.testing.expectEqual(-0x1.e33af14efca0bp+0, sinh(@as(f64, -0x1.646852p+0)));
    // try std.testing.expectEqual(-0x1.e33aef14d1eap+0, sinh(@as(f64, -0x1.646850f515ef2p+0)));
    try std.testing.expectEqual(-0x7.f48612b1b30ecp-4, sinh(@as(f64, -0x7.a8c5fp-4)));
    // try std.testing.expectEqual(-0x7.f4861ba0df664p-4, sinh(@as(f64, -0x7.a8c5f8p-4)));
    try std.testing.expectEqual(-0x7.f4861a01ff01p-4, sinh(@as(f64, -0x7.a8c5f68c81facp-4)));
    try std.testing.expectEqual(-0x7.f4861a01ff014p-4, sinh(@as(f64, -0x7.a8c5f68c81fbp-4)));
    try std.testing.expectEqual(0x3.4ff4d6729691p-4, sinh(@as(f64, 0x3.4a037p-4)));
    // try std.testing.expectEqual(-0x3.f5b9aacdd086p-4, sinh(@as(f64, -0x3.eba6d8p-4)));
    // try std.testing.expectEqual(-0x3.f5b9aeecb5a4p-4, sinh(@as(f64, -0x3.eba6dcp-4)));
    // try std.testing.expectEqual(-0x3.f5b9aeb710594p-4, sinh(@as(f64, -0x3.eba6dbcbeceb2p-4)));
    try std.testing.expectEqual(-0x5.1ed4b3c8c4e08p+0, sinh(@as(f64, -0x2.55f63p+0)));
    // try std.testing.expectEqual(-0x3.d3834c8e189cp-4, sinh(@as(f64, -0x3.ca68c8p-4)));
    try std.testing.expectEqual(-0x3.d38350aaf8128p-4, sinh(@as(f64, -0x3.ca68ccp-4)));
    try std.testing.expectEqual(-0x3.d3834dfb540d6p-4, sinh(@as(f64, -0x3.ca68c96337692p-4)));
    try std.testing.expectEqual(-0x3.9a7a1fd80eae2p-4, sinh(@as(f64, -0x3.92da04p-4)));
    try std.testing.expectEqual(-0x3.9a7a23f1b49bap-4, sinh(@as(f64, -0x3.92da08p-4)));
    try std.testing.expectEqual(-0x3.9a7a218aff89p-4, sinh(@as(f64, -0x3.92da05a85024ap-4)));
    // try std.testing.expectEqual(-0x3.9a7a218aff892p-4, sinh(@as(f64, -0x3.92da05a85024cp-4)));
    try std.testing.expectEqual(-0x3.4415b63bf6484p-4, sinh(@as(f64, -0x3.3e629p-4)));
    try std.testing.expectEqual(-0x3.4415ba5113c8ap-4, sinh(@as(f64, -0x3.3e6294p-4)));
    try std.testing.expectEqual(-0x3.4415b938adb18p-4, sinh(@as(f64, -0x3.3e6292ed442d4p-4)));
    try std.testing.expectEqual(-0x3.4415b938adb1ap-4, sinh(@as(f64, -0x3.3e6292ed442d6p-4)));
    try std.testing.expectEqual(0x7.b341dd42fddep-4, sinh(@as(f64, 0x7.6e25ap-4)));
    try std.testing.expectEqual(0x7.b341d46228bc8p-4, sinh(@as(f64, 0x7.6e2598p-4)));
    try std.testing.expectEqual(0x7.b341da16deb5cp-4, sinh(@as(f64, 0x7.6e259d2436fc4p-4)));
    try std.testing.expectEqual(0x3.e0563601aac3ep-4, sinh(@as(f64, 0x3.d6e088p-4)));
    try std.testing.expectEqual(-0x7.ad0e3c83adf18p-4, sinh(@as(f64, -0x7.688eap-4)));
    // try std.testing.expectEqual(-0xf.a9e6db74e248p-4, sinh(@as(f64, -0xd.dce79p-4)));
    try std.testing.expectEqual(0x5.2a5fdb392d918p+196, sinh(@as(f64, 0x8.a3127p+4)));
    try std.testing.expectEqual(0x1.c07090e55732ap-12, sinh(@as(f64, 0x1.c0709p-12)));
    try std.testing.expectEqual(0xc.835a5a1df79bp-4, sinh(@as(f64, 0xb.7f67dp-4)));
    // try std.testing.expectEqual(0xc.835a45ce17f9p-4, sinh(@as(f64, 0xb.7f67cp-4)));
    // try std.testing.expectEqual(0xc.835a4a0d527d8p-4, sinh(@as(f64, 0xb.7f67c3586c24p-4)));
    try std.testing.expectEqual(-0x1.960d6e6e4b63cp+0, sinh(@as(f64, -0x1.3dda8ap+0)));
    // try std.testing.expectEqual(-0x6.119951b224aa4p-4, sinh(@as(f64, -0x5.ee9218p-4)));
    try std.testing.expectEqual(-0x2.c1769e4cedb6ap+0, sinh(@as(f64, -0x1.bcfc98p+0)));
    // try std.testing.expectEqual(-0x6.cc3dd8844c26p-4, sinh(@as(f64, -0x6.9bbb68p-4)));
    try std.testing.expectEqual(-0x6.cc3de135798dp-4, sinh(@as(f64, -0x6.9bbb7p-4)));
    // try std.testing.expectEqual(-0x6.cc3ddf003dcdcp-4, sinh(@as(f64, -0x6.9bbb6df7c5d08p-4)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d376167f406p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d376167f404p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-std.math.inf(f64), sinh(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-std.math.inf(f64), sinh(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4p-128, sinh(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, sinh(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, sinh(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, sinh(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, sinh(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, sinh(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, sinh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sinh(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0xf.fffec1f47394p+124, sinh(@as(f64, 0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48ep+128, sinh(@as(f64, 0x5.96a7e8p+4)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c679d4p+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, sinh(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, sinh(@as(f64, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c679d4p+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, sinh(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c679d4p+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, sinh(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, sinh(@as(f64, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c679d4p+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, sinh(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, sinh(@as(f64, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), sinh(@as(f64, 0x2.c5d37700c6bbp+12)));
    // try std.testing.expectEqual(-0xf.ef296e7b88b4p+1020, sinh(@as(f64, -0x2.c678c4p+8)));
    // try std.testing.expectEqual(-0xf.ef692ba0bc978p+1020, sinh(@as(f64, -0x2.c678c8p+8)));
    // try std.testing.expectEqual(-0xf.ef3a7e711d2c8p+1020, sinh(@as(f64, -0x2.c678c5121f428p+8)));

    try std.testing.expectEqual(0x0p+0, sinh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0xd.283596e9e347f2fp-4, sinh(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x8.000000000000005p-32, sinh(@as(f80, 0x8p-32)));
    try std.testing.expectEqual(0x8.00555666680681ep-8, sinh(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(-0x8.00555666680681ep-8, sinh(@as(f80, -0x8p-8)));
    try std.testing.expectEqual(0x4.00000aaaaab3333p-12, sinh(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(-0x4.00000aaaaab3333p-12, sinh(@as(f80, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000000002aaaaap-20, sinh(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(-0x1.00000000002aaaaap-20, sinh(@as(f80, -0x1p-20)));
    try std.testing.expectEqual(0x4.0000000000000008p-32, sinh(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(-0x4.0000000000000008p-32, sinh(@as(f80, -0x4p-32)));
    try std.testing.expectEqual(0x1p-40, sinh(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(-0x1p-40, sinh(@as(f80, -0x1p-40)));
    try std.testing.expectEqual(0x4p-52, sinh(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(-0x4p-52, sinh(@as(f80, -0x4p-52)));
    try std.testing.expectEqual(0x1p-60, sinh(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(-0x1p-60, sinh(@as(f80, -0x1p-60)));
    try std.testing.expectEqual(0x4p-72, sinh(@as(f80, 0x4p-72)));
    try std.testing.expectEqual(-0x4p-72, sinh(@as(f80, -0x4p-72)));
    try std.testing.expectEqual(0x1p-100, sinh(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, sinh(@as(f80, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sinh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-1000, sinh(@as(f80, 0x1p-1000)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-1000, sinh(@as(f80, -0x1p-1000)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sinh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, sinh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, sinh(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sinh(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1p-10000, sinh(@as(f80, -0x1p-10000)));
    try std.testing.expectEqual(0x6.ad6b6e710d7fe068p+28, sinh(@as(f80, 0x1.6p+4)));
    try std.testing.expectEqual(0x1.226af33b1fdc0a58p+32, sinh(@as(f80, 0x1.7p+4)));
    try std.testing.expectEqual(0x3.156ff6a8ebf6e67p+32, sinh(@as(f80, 0x1.8p+4)));
    try std.testing.expectEqual(-0x7.9852071dfda98d78p-4, sinh(@as(f80, -0x7.55d7f8p-4)));
    try std.testing.expectEqual(-0x3.fde378210a8f8b14p-4, sinh(@as(f80, -0x3.f392f8p-4)));
    try std.testing.expectEqual(0x2.da7cd9753b47bb4cp+0, sinh(@as(f80, 0x1.c56446p+0)));
    try std.testing.expectEqual(0x6.ff781ca6e6dca67p-4, sinh(@as(f80, 0x6.cac628p-4)));
    try std.testing.expectEqual(0x6.ff7813eb9593d89p-4, sinh(@as(f80, 0x6.cac62p-4)));
    try std.testing.expectEqual(0x6.ff78170306f9cc1p-4, sinh(@as(f80, 0x6.cac622d51eebcp-4)));
    try std.testing.expectEqual(-0xa.0100cebf41c7d7p+4, sinh(@as(f80, -0x5.c4cbp+0)));
    try std.testing.expectEqual(-0xa.01011ec7afd9d17p+4, sinh(@as(f80, -0x5.c4cb08p+0)));
    try std.testing.expectEqual(-0xa.0100e4f7b10f8b9p+4, sinh(@as(f80, -0x5.c4cb02389c094p+0)));
    try std.testing.expectEqual(-0x1.e33aed09484005ep+0, sinh(@as(f80, -0x1.64685p+0)));
    try std.testing.expectEqual(-0x1.e33af14efca0a8cap+0, sinh(@as(f80, -0x1.646852p+0)));
    try std.testing.expectEqual(-0x1.e33aef14d1ea06aep+0, sinh(@as(f80, -0x1.646850f515ef2p+0)));
    try std.testing.expectEqual(-0x7.f48612b1b30ec88p-4, sinh(@as(f80, -0x7.a8c5fp-4)));
    try std.testing.expectEqual(-0x7.f4861ba0df663478p-4, sinh(@as(f80, -0x7.a8c5f8p-4)));
    try std.testing.expectEqual(-0x7.f4861a01ff00e128p-4, sinh(@as(f80, -0x7.a8c5f68c81facp-4)));
    try std.testing.expectEqual(-0x7.f4861a01ff0128ap-4, sinh(@as(f80, -0x7.a8c5f68c81fbp-4)));
    try std.testing.expectEqual(-0x7.f4861a01ff010b6p-4, sinh(@as(f80, -0x7.a8c5f68c81fae5dp-4)));
    try std.testing.expectEqual(0x3.4ff4d672969101b8p-4, sinh(@as(f80, 0x3.4a037p-4)));
    try std.testing.expectEqual(-0x3.f5b9aacdd0860a0cp-4, sinh(@as(f80, -0x3.eba6d8p-4)));
    try std.testing.expectEqual(-0x3.f5b9aeecb5a3f93p-4, sinh(@as(f80, -0x3.eba6dcp-4)));
    try std.testing.expectEqual(-0x3.f5b9aeb7105930f8p-4, sinh(@as(f80, -0x3.eba6dbcbeceb2p-4)));
    try std.testing.expectEqual(-0x5.1ed4b3c8c4e07e8p+0, sinh(@as(f80, -0x2.55f63p+0)));
    try std.testing.expectEqual(-0x3.d3834c8e189bfe9p-4, sinh(@as(f80, -0x3.ca68c8p-4)));
    try std.testing.expectEqual(-0x3.d38350aaf8127c6p-4, sinh(@as(f80, -0x3.ca68ccp-4)));
    try std.testing.expectEqual(-0x3.d3834dfb540d632cp-4, sinh(@as(f80, -0x3.ca68c96337692p-4)));
    try std.testing.expectEqual(-0x3.9a7a1fd80eae25e8p-4, sinh(@as(f80, -0x3.92da04p-4)));
    try std.testing.expectEqual(-0x3.9a7a23f1b49b9544p-4, sinh(@as(f80, -0x3.92da08p-4)));
    try std.testing.expectEqual(-0x3.9a7a218aff88f068p-4, sinh(@as(f80, -0x3.92da05a85024ap-4)));
    try std.testing.expectEqual(-0x3.9a7a218aff891138p-4, sinh(@as(f80, -0x3.92da05a85024cp-4)));
    try std.testing.expectEqual(-0x3.9a7a218aff8903f8p-4, sinh(@as(f80, -0x3.92da05a85024b314p-4)));
    try std.testing.expectEqual(-0x3.4415b63bf64837f4p-4, sinh(@as(f80, -0x3.3e629p-4)));
    try std.testing.expectEqual(-0x3.4415ba5113c8a3bcp-4, sinh(@as(f80, -0x3.3e6294p-4)));
    try std.testing.expectEqual(-0x3.4415b938adb1860cp-4, sinh(@as(f80, -0x3.3e6292ed442d4p-4)));
    try std.testing.expectEqual(-0x3.4415b938adb1a6b4p-4, sinh(@as(f80, -0x3.3e6292ed442d6p-4)));
    try std.testing.expectEqual(-0x3.4415b938adb18b34p-4, sinh(@as(f80, -0x3.3e6292ed442d450cp-4)));
    try std.testing.expectEqual(0x7.b341dd42fdddea3p-4, sinh(@as(f80, 0x7.6e25ap-4)));
    try std.testing.expectEqual(0x7.b341d46228bc9ec8p-4, sinh(@as(f80, 0x7.6e2598p-4)));
    try std.testing.expectEqual(0x7.b341da16deb5db08p-4, sinh(@as(f80, 0x7.6e259d2436fc4p-4)));
    try std.testing.expectEqual(0x3.e0563601aac3ea64p-4, sinh(@as(f80, 0x3.d6e088p-4)));
    try std.testing.expectEqual(-0x7.ad0e3c83adf17bfp-4, sinh(@as(f80, -0x7.688eap-4)));
    try std.testing.expectEqual(-0xf.a9e6db74e247cefp-4, sinh(@as(f80, -0xd.dce79p-4)));
    try std.testing.expectEqual(0x5.2a5fdb392d919fcp+196, sinh(@as(f80, 0x8.a3127p+4)));
    try std.testing.expectEqual(0x1.c07090e55732a002p-12, sinh(@as(f80, 0x1.c0709p-12)));
    try std.testing.expectEqual(0xc.835a5a1df79ae5fp-4, sinh(@as(f80, 0xb.7f67dp-4)));
    try std.testing.expectEqual(0xc.835a45ce17f9353p-4, sinh(@as(f80, 0xb.7f67cp-4)));
    try std.testing.expectEqual(0xc.835a4a0d527d5p-4, sinh(@as(f80, 0xb.7f67c3586c24p-4)));
    try std.testing.expectEqual(-0x1.960d6e6e4b63c67p+0, sinh(@as(f80, -0x1.3dda8ap+0)));
    try std.testing.expectEqual(-0x6.119951b224aa2ab8p-4, sinh(@as(f80, -0x5.ee9218p-4)));
    try std.testing.expectEqual(-0x2.c1769e4cedb691dcp+0, sinh(@as(f80, -0x1.bcfc98p+0)));
    try std.testing.expectEqual(-0x6.cc3dd8844c261e98p-4, sinh(@as(f80, -0x6.9bbb68p-4)));
    try std.testing.expectEqual(-0x6.cc3de135798d1c48p-4, sinh(@as(f80, -0x6.9bbb7p-4)));
    try std.testing.expectEqual(-0x6.cc3ddf003dcda78p-4, sinh(@as(f80, -0x6.9bbb6df7c5d08p-4)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, sinh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(0xf.ff15bf3871a7576p+16380, sinh(@as(f80, 0x2.c5d376167f406p+12)));
    try std.testing.expectEqual(0xf.ff15bf3851a92bep+16380, sinh(@as(f80, 0x2.c5d376167f404p+12)));
    try std.testing.expectEqual(0xf.ff15bf38649c166p+16380, sinh(@as(f80, 0x2.c5d376167f4052f4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-std.math.inf(f80), sinh(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-std.math.inf(f80), sinh(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-std.math.inf(f80), sinh(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4p-128, sinh(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, sinh(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, sinh(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, sinh(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, sinh(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, sinh(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, sinh(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, sinh(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, sinh(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, sinh(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, sinh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, sinh(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sinh(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, sinh(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0xf.fffec1f473940d2p+124, sinh(@as(f80, 0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48e748p+128, sinh(@as(f80, 0x5.96a7e8p+4)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, sinh(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, sinh(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, sinh(@as(f80, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, sinh(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, sinh(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, sinh(@as(f80, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, sinh(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, sinh(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, sinh(@as(f80, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, sinh(@as(f80, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2dp+1020, sinh(@as(f80, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2dp+1020, sinh(@as(f80, 0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, sinh(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, sinh(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, sinh(@as(f80, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, sinh(@as(f80, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2dp+1020, sinh(@as(f80, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2dp+1020, sinh(@as(f80, 0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, sinh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, sinh(@as(f80, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, sinh(@as(f80, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, sinh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, sinh(@as(f80, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, sinh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, sinh(@as(f80, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, sinh(@as(f80, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, sinh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, sinh(@as(f80, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), sinh(@as(f80, 0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, sinh(@as(f80, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(-0xf.ef296e7b88b41f6p+1020, sinh(@as(f80, -0x2.c678c4p+8)));
    try std.testing.expectEqual(-0xf.ef692ba0bc97addp+1020, sinh(@as(f80, -0x2.c678c8p+8)));
    try std.testing.expectEqual(-0xf.ef3a7e711d2c75ap+1020, sinh(@as(f80, -0x2.c678c5121f428p+8)));

    try std.testing.expectEqual(0x0p+0, sinh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0xd.283596e9e347f2ee3cf47bf04bp-4, sinh(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x8.0000000000000055555555555558p-32, sinh(@as(f128, 0x8p-32)));
    try std.testing.expectEqual(0x8.00555666680681d9e591eff67c8p-8, sinh(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(-0x8.00555666680681d9e591eff67c8p-8, sinh(@as(f128, -0x8p-8)));
    try std.testing.expectEqual(0x4.00000aaaaab33333367367372c58p-12, sinh(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(-0x4.00000aaaaab33333367367372c58p-12, sinh(@as(f128, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000000002aaaaaaaaaaccccccdp-20, sinh(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(-0x1.00000000002aaaaaaaaaaccccccdp-20, sinh(@as(f128, -0x1p-20)));
    // try std.testing.expectEqual(0x4.000000000000000aaaaaaaaaaaacp-32, sinh(@as(f128, 0x4p-32)));
    // try std.testing.expectEqual(-0x4.000000000000000aaaaaaaaaaaacp-32, sinh(@as(f128, -0x4p-32)));
    try std.testing.expectEqual(0x1.000000000000000000002aaaaaabp-40, sinh(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(-0x1.000000000000000000002aaaaaabp-40, sinh(@as(f128, -0x1p-40)));
    try std.testing.expectEqual(0x4.0000000000000000000000000aacp-52, sinh(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(-0x4.0000000000000000000000000aacp-52, sinh(@as(f128, -0x4p-52)));
    try std.testing.expectEqual(0x1p-60, sinh(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(-0x1p-60, sinh(@as(f128, -0x1p-60)));
    try std.testing.expectEqual(0x4p-72, sinh(@as(f128, 0x4p-72)));
    try std.testing.expectEqual(-0x4p-72, sinh(@as(f128, -0x4p-72)));
    try std.testing.expectEqual(0x1p-100, sinh(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, sinh(@as(f128, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sinh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-1000, sinh(@as(f128, 0x1p-1000)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-1000, sinh(@as(f128, -0x1p-1000)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sinh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, sinh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, sinh(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(-0x0p+0, sinh(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sinh(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1p-10000, sinh(@as(f128, -0x1p-10000)));
    try std.testing.expectEqual(0x6.ad6b6e710d7fe065377669cb23p+28, sinh(@as(f128, 0x1.6p+4)));
    try std.testing.expectEqual(0x1.226af33b1fdc0a574c76ab216131p+32, sinh(@as(f128, 0x1.7p+4)));
    try std.testing.expectEqual(0x3.156ff6a8ebf6e66f1fb318fc8d6ap+32, sinh(@as(f128, 0x1.8p+4)));
    // try std.testing.expectEqual(-0x7.9852071dfda98d7a78bbfbeba1ccp-4, sinh(@as(f128, -0x7.55d7f8p-4)));
    try std.testing.expectEqual(-0x3.fde378210a8f8b139f4bf3918742p-4, sinh(@as(f128, -0x3.f392f8p-4)));
    // try std.testing.expectEqual(0x2.da7cd9753b47bb4b1c1b52331194p+0, sinh(@as(f128, 0x1.c56446p+0)));
    try std.testing.expectEqual(0x6.ff781ca6e6dca66ef1b4884e83a8p-4, sinh(@as(f128, 0x6.cac628p-4)));
    // try std.testing.expectEqual(0x6.ff7813eb9593d88f40043b863e1p-4, sinh(@as(f128, 0x6.cac62p-4)));
    // try std.testing.expectEqual(0x6.ff78170306f9cc0e1b26024a3bc4p-4, sinh(@as(f128, 0x6.cac622d51eebcp-4)));
    // try std.testing.expectEqual(-0xa.0100cebf41c7d702ca8ab889ec5p+4, sinh(@as(f128, -0x5.c4cbp+0)));
    // try std.testing.expectEqual(-0xa.01011ec7afd9d171664f8ff70a2p+4, sinh(@as(f128, -0x5.c4cb08p+0)));
    try std.testing.expectEqual(-0xa.0100e4f7b10f8b8ac75b9651a52p+4, sinh(@as(f128, -0x5.c4cb02389c094p+0)));
    try std.testing.expectEqual(-0x1.e33aed09484005e089b161278fe1p+0, sinh(@as(f128, -0x1.64685p+0)));
    try std.testing.expectEqual(-0x1.e33af14efca0a8c9d1b749c14b43p+0, sinh(@as(f128, -0x1.646852p+0)));
    try std.testing.expectEqual(-0x1.e33aef14d1ea06ad3181d20e64d9p+0, sinh(@as(f128, -0x1.646850f515ef2p+0)));
    try std.testing.expectEqual(-0x7.f48612b1b30ec87e5f34bd28f338p-4, sinh(@as(f128, -0x7.a8c5fp-4)));
    try std.testing.expectEqual(-0x7.f4861ba0df66347602112b03c204p-4, sinh(@as(f128, -0x7.a8c5f8p-4)));
    // try std.testing.expectEqual(-0x7.f4861a01ff00e1242960cce2c038p-4, sinh(@as(f128, -0x7.a8c5f68c81facp-4)));
    // try std.testing.expectEqual(-0x7.f4861a01ff01289d8c213e53b8dcp-4, sinh(@as(f128, -0x7.a8c5f68c81fbp-4)));
    try std.testing.expectEqual(-0x7.f4861a01ff010b5ea0f8ffe8bf1p-4, sinh(@as(f128, -0x7.a8c5f68c81fae5dp-4)));
    try std.testing.expectEqual(0x3.4ff4d672969101b81d84d928cd8p-4, sinh(@as(f128, 0x3.4a037p-4)));
    // try std.testing.expectEqual(-0x3.f5b9aacdd0860a0ddd86f9d6f6d6p-4, sinh(@as(f128, -0x3.eba6d8p-4)));
    // try std.testing.expectEqual(-0x3.f5b9aeecb5a3f92efdc01ca9652ap-4, sinh(@as(f128, -0x3.eba6dcp-4)));
    try std.testing.expectEqual(-0x3.f5b9aeb7105930f8f9931b42c7dcp-4, sinh(@as(f128, -0x3.eba6dbcbeceb2p-4)));
    try std.testing.expectEqual(-0x5.1ed4b3c8c4e07e8146d7a23bf618p+0, sinh(@as(f128, -0x2.55f63p+0)));
    try std.testing.expectEqual(-0x3.d3834c8e189bfe8faf608f8faaa6p-4, sinh(@as(f128, -0x3.ca68c8p-4)));
    try std.testing.expectEqual(-0x3.d38350aaf8127c5fcbc7cc89b398p-4, sinh(@as(f128, -0x3.ca68ccp-4)));
    // try std.testing.expectEqual(-0x3.d3834dfb540d632c64fc59c88ddep-4, sinh(@as(f128, -0x3.ca68c96337692p-4)));
    // try std.testing.expectEqual(-0x3.9a7a1fd80eae25e868477257719ep-4, sinh(@as(f128, -0x3.92da04p-4)));
    try std.testing.expectEqual(-0x3.9a7a23f1b49b9544a72e42eb797ep-4, sinh(@as(f128, -0x3.92da08p-4)));
    // try std.testing.expectEqual(-0x3.9a7a218aff88f069bca6fc2c119ep-4, sinh(@as(f128, -0x3.92da05a85024ap-4)));
    // try std.testing.expectEqual(-0x3.9a7a218aff891136ec124f8ee298p-4, sinh(@as(f128, -0x3.92da05a85024cp-4)));
    // try std.testing.expectEqual(-0x3.9a7a218aff8903f8110c58c2db72p-4, sinh(@as(f128, -0x3.92da05a85024b314p-4)));
    try std.testing.expectEqual(-0x3.4415b63bf64837f33b46cc49daf2p-4, sinh(@as(f128, -0x3.3e629p-4)));
    try std.testing.expectEqual(-0x3.4415ba5113c8a3baf20fb60a63b4p-4, sinh(@as(f128, -0x3.3e6294p-4)));
    try std.testing.expectEqual(-0x3.4415b938adb1860ce25000fa7e2p-4, sinh(@as(f128, -0x3.3e6292ed442d4p-4)));
    // try std.testing.expectEqual(-0x3.4415b938adb1a6b5ce53c0150d0ap-4, sinh(@as(f128, -0x3.3e6292ed442d6p-4)));
    try std.testing.expectEqual(-0x3.4415b938adb18b338689183e4e6ap-4, sinh(@as(f128, -0x3.3e6292ed442d450cp-4)));
    try std.testing.expectEqual(0x7.b341dd42fdddea2ea21889a6e4f8p-4, sinh(@as(f128, 0x7.6e25ap-4)));
    try std.testing.expectEqual(0x7.b341d46228bc9ecadcfa9ca3b23p-4, sinh(@as(f128, 0x7.6e2598p-4)));
    try std.testing.expectEqual(0x7.b341da16deb5db07a018f01fd9dcp-4, sinh(@as(f128, 0x7.6e259d2436fc4p-4)));
    try std.testing.expectEqual(0x3.e0563601aac3ea656b93e0306268p-4, sinh(@as(f128, 0x3.d6e088p-4)));
    try std.testing.expectEqual(-0x7.ad0e3c83adf17bed0e1571979c4p-4, sinh(@as(f128, -0x7.688eap-4)));
    try std.testing.expectEqual(-0xf.a9e6db74e247cef34f74103a47cp-4, sinh(@as(f128, -0xd.dce79p-4)));
    try std.testing.expectEqual(0x5.2a5fdb392d919fc3f2ab6db2987cp+196, sinh(@as(f128, 0x8.a3127p+4)));
    try std.testing.expectEqual(0x1.c07090e55732a001dde433d77237p-12, sinh(@as(f128, 0x1.c0709p-12)));
    // try std.testing.expectEqual(0xc.835a5a1df79ae5ec9c48cbb6bf08p-4, sinh(@as(f128, 0xb.7f67dp-4)));
    try std.testing.expectEqual(0xc.835a45ce17f9353505896689a24p-4, sinh(@as(f128, 0xb.7f67cp-4)));
    // try std.testing.expectEqual(0xc.835a4a0d527d4fff7b247fd0064p-4, sinh(@as(f128, 0xb.7f67c3586c24p-4)));
    try std.testing.expectEqual(-0x1.960d6e6e4b63c66ff64892c1bf37p+0, sinh(@as(f128, -0x1.3dda8ap+0)));
    try std.testing.expectEqual(-0x6.119951b224aa2ab9c11796817da8p-4, sinh(@as(f128, -0x5.ee9218p-4)));
    try std.testing.expectEqual(-0x2.c1769e4cedb691dd692f866d5f38p+0, sinh(@as(f128, -0x1.bcfc98p+0)));
    // try std.testing.expectEqual(-0x6.cc3dd8844c261e97377c1f6d1564p-4, sinh(@as(f128, -0x6.9bbb68p-4)));
    try std.testing.expectEqual(-0x6.cc3de135798d1c4ad5b34b3db2ep-4, sinh(@as(f128, -0x6.9bbb7p-4)));
    // try std.testing.expectEqual(-0x6.cc3ddf003dcda77f8f9e892e36d8p-4, sinh(@as(f128, -0x6.9bbb6df7c5d08p-4)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, sinh(@as(f128, 0x2.c5d374p+12)));
    // try std.testing.expectEqual(0xf.ff15bf3871a75761db61506a9bbp+16380, sinh(@as(f128, 0x2.c5d376167f406p+12)));
    // try std.testing.expectEqual(0xf.ff15bf3851a92be36a9dffe75668p+16380, sinh(@as(f128, 0x2.c5d376167f404p+12)));
    try std.testing.expectEqual(0xf.ff15bf38649c16662e1ff4acb94p+16380, sinh(@as(f128, 0x2.c5d376167f4052f4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-std.math.inf(f128), sinh(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-std.math.inf(f128), sinh(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-std.math.inf(f128), sinh(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-std.math.inf(f128), sinh(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-std.math.inf(f128), sinh(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4p-128, sinh(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, sinh(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, sinh(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, sinh(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, sinh(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, sinh(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, sinh(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, sinh(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, sinh(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, sinh(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, sinh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, sinh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, sinh(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, sinh(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, sinh(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-152, sinh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, sinh(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, sinh(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, sinh(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, sinh(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0xf.fffec1f473940d22f2195eac65ep+124, sinh(@as(f128, 0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48e7480e07d1c02e7cp+128, sinh(@as(f128, 0x5.96a7e8p+4)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, sinh(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, sinh(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, sinh(@as(f128, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, sinh(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, sinh(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, sinh(@as(f128, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, sinh(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, sinh(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, sinh(@as(f128, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, sinh(@as(f128, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2ca74dec5830328p+1020, sinh(@as(f128, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2ca74dec58303ep+1020, sinh(@as(f128, 0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffbffffffffff303a8p+1020, sinh(@as(f128, 0x2.c679d1f73f0fb624d358b213a7p+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, sinh(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, sinh(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, sinh(@as(f128, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, sinh(@as(f128, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2ca74dec5830328p+1020, sinh(@as(f128, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2ca74dec58303ep+1020, sinh(@as(f128, 0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffc0000000000303a8p+1020, sinh(@as(f128, 0x2.c679d1f73f0fb624d358b213a8p+8)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, sinh(@as(f128, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, sinh(@as(f128, 0x2.c5d37700c6bbp+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, sinh(@as(f128, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, sinh(@as(f128, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, sinh(@as(f128, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, sinh(@as(f128, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, sinh(@as(f128, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb03a8p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, sinh(@as(f128, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffe61p+16380, sinh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b494cp+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b4ap+12)));
    // try std.testing.expectEqual(0xf.ffffffffffffffffffffffb3e61p+16380, sinh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b49p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, sinh(@as(f128, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, sinh(@as(f128, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb03a8p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, sinh(@as(f128, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b494ep+12)));
    try std.testing.expectEqual(std.math.inf(f128), sinh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b4ap+12)));
    // try std.testing.expectEqual(0xf.ffffffffffffffffffffffb3e61p+16380, sinh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b49p+12)));
    try std.testing.expectEqual(-0xf.ef296e7b88b41f625301a966ffap+1020, sinh(@as(f128, -0x2.c678c4p+8)));
    try std.testing.expectEqual(-0xf.ef692ba0bc97adc852e8e105932p+1020, sinh(@as(f128, -0x2.c678c8p+8)));
    try std.testing.expectEqual(-0xf.ef3a7e711d2c75a66ea3ca1c417p+1020, sinh(@as(f128, -0x2.c678c5121f428p+8)));
}
