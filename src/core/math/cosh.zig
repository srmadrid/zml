const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const roundeven = @import("roundeven.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn cosh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return cosh(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, cosh32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_coshf.c
                    return cosh32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_cosh.c
                    return cosh64(x);
                },
                f80 => return cast(f80, cosh128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_coshl.c
                    return cosh128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn cosh32(x: f32) f32 {
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
    const iln2: f64 = 0x1.71547652b82fep+5;
    const z: f64 = cast(f64, x, .{});
    const ax: u32 = @as(u32, @bitCast(x)) << 1;
    if (ax > 0x8565a9f8) { // |x| >~ 89.4
        @branchHint(.unlikely);
        if (ax >= 0xff000000) {
            if ((ax << 8) != 0)
                return x + x; // nan

            return std.math.inf(f32); // +-inf
        }

        return 0x1p97 * 0x1p97;
    }
    if (ax < 0x7c000000) { // |x| < 0.125
        @branchHint(.unlikely);
        if (ax < 0x74000000) { // |x| < 0x1p-11
            @branchHint(.unlikely);
            if (ax < 0x66000000) { // |x| < 0x1p-24
                @branchHint(.unlikely);
                return @mulAdd(f32, math.abs(x), 0x1p-25, 1);
            }

            return (0.5 * x) * x + 1;
        }
        const cp: [4]f64 = .{
            0x1.fffffffffffe3p-2,  0x1.55555555723cfp-5,
            0x1.6c16bee4a5986p-10, 0x1.a0483fc0328f7p-16,
        };
        const z2: f64 = z * z;
        const z4: f64 = z2 * z2;
        return cast(f32, 1 + z2 * ((cp[0] + z2 * cp[1]) + z4 * (cp[2] + z2 * (cp[3]))), .{});
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
    var r: f64 = rp + rm;
    var ub: f32 = cast(f32, r, .{});
    const lb: f32 = cast(f32, r - 1.45e-10 * r, .{});
    if (ub != lb) {
        @branchHint(.unlikely);
        const iln2h: f64 = 0x1.7154765p+5;
        const iln2l: f64 = 0x1.5c17f0bbbe88p-26;
        h = (iln2h * z - ia) + iln2l * z;
        h2 = h * h;
        te = ch[0] + h2 * ch[2] + (h2 * h2) * (ch[4] + h2 * ch[6]);
        to = ch[1] + h2 * (ch[3] + h2 * ch[5]);
        r = sp * (te + h * to) + sm * (te - h * to);
        ub = cast(f32, r, .{});
    }

    return ub;
}

fn cosh64(x: f64) f64 {
    const huge: f64 = 1.0e300;

    // High word of |x|.
    var ix: i32 = undefined;
    dbl64.getHighWord(&ix, x);
    ix &= 0x7fffffff;

    // |x| in [0,22]
    if (ix < 0x40360000) {
        // |x| in [0,0.5*ln2], return 1+expm1(|x|)^2/(2*exp(|x|))
        if (ix < 0x3fd62e43) {
            if (ix < 0x3c800000) // cosh(tiny) = 1
                return 1;

            const t: f64 = math.expm1(math.abs(x));
            const w: f64 = 1 + t;
            return 1 + (t * t) / (w + w);
        }

        // |x| in [0.5*ln2,22], return (exp(|x|)+1/exp(|x|)/2;
        const t: f64 = math.exp(math.abs(x));
        return 0.5 * t + 0.5 / t;
    }

    // |x| in [22, log(maxdouble)] return 0.5*exp(|x|)
    if (ix < 0x40862e42)
        return 0.5 * math.exp(math.abs(x));

    // |x| in [log(maxdouble), overflowthresold]
    var fix: i64 = undefined;
    dbl64.extractWords64(&fix, x);
    fix &= 0x7fffffffffffffff;
    if (fix <= 0x408633ce8fb9f87d) {
        const w: f64 = math.exp(0.5 * math.abs(x));
        const t: f64 = 0.5 * w;
        return t * w;
    }

    // x is INF or NaN
    if (ix >= 0x7ff00000)
        return x * x;

    // |x| > overflowthresold, cosh(x) overflow
    return huge * huge;
}

fn cosh128(x: f128) f128 {
    const huge: f128 = 1.0e4900;
    const ovf_thresh: f128 = 1.1357216553474703894801348310092223067821e4;

    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const ex: i32 = @bitCast(u.w0 & 0x7fffffff);

    // Absolute value of x.
    u.w0 = @bitCast(ex);

    // x is INF or NaN
    if (ex >= 0x7fff0000)
        return x * x;

    // |x| in [0,0.5*ln2], return 1+expm1l(|x|)^2/(2*expl(|x|))
    if (ex < 0x3ffd62e4) { // 0.3465728759765625
        if (ex < 0x3fb80000) // |x| < 2^-116
            return 1; // cosh(tiny) = 1

        const t: f128 = math.expm1(@as(f128, @bitCast(u)));
        const w: f128 = 1 + t;
        return 1 + (t * t) / (w + w);
    }

    // |x| in [0.5*ln2,40], return (exp(|x|)+1/exp(|x|)/2;
    if (ex < 0x40044000) {
        const t: f128 = math.exp(@as(f128, @bitCast(u)));
        return 0.5 * t + 0.5 / t;
    }

    // |x| in [22, ln(maxdouble)] return 0.5*exp(|x|)
    if (ex <= 0x400c62e3) // 11356.375
        return 0.5 * math.exp(@as(f128, @bitCast(u)));

    // |x| in [log(maxdouble), overflowthresold]
    if (@as(f128, @bitCast(u)) <= ovf_thresh) {
        const w: f128 = math.exp(0.5 * @as(f128, @bitCast(u)));
        const t: f128 = 0.5 * w;
        return t * w;
    }

    // |x| > overflowthresold, cosh(x) overflow
    return huge * huge;
}

test cosh {
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1.4b705ep+0, cosh(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5e3bp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5e3acp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5e3acp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5e3bp+8)));
    try std.testing.expectEqual(0x6.ad6b7p+28, cosh(@as(f32, 0x1.6p+4)));
    try std.testing.expectEqual(0x1.226af4p+32, cosh(@as(f32, 0x1.7p+4)));
    try std.testing.expectEqual(0x3.156ff8p+32, cosh(@as(f32, 0x1.8p+4)));
    try std.testing.expectEqual(0x1.002p+0, cosh(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x1.000008p+0, cosh(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1.8b0756p+0, cosh(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x8.c881fp+68, cosh(@as(f32, 0x3.2p+4)));
    try std.testing.expectEqual(0xa.a717ap+12, cosh(@as(f32, -0xb.60713p+0)));
    try std.testing.expectEqual(0x1.68b8dcp+4, cosh(@as(f32, -0x3.cee48p+0)));
    try std.testing.expectEqual(0x9.ad527p+0, cosh(@as(f32, 0x2.f5d128p+0)));
    try std.testing.expectEqual(0x3.89993cp+16, cosh(@as(f32, -0xd.0c03p+0)));
    try std.testing.expectEqual(0x1.074e54p+0, cosh(@as(f32, -0x3.d04328p-4)));
    try std.testing.expectEqual(0x1.074e54p+0, cosh(@as(f32, -0x3.d0432cp-4)));
    try std.testing.expectEqual(0x7.d7161p+28, cosh(@as(f32, 0x1.629188p+4)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x1p-56)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, -0x1p-56)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, 0x1p-72)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f32, -0x1p-72)));
    try std.testing.expectEqual(0xf.fffecp+124, cosh(@as(f32, 0x5.96a7ep+4)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x5.96a7e8p+4)));
    try std.testing.expectEqual(0xf.fffecp+124, cosh(@as(f32, -0x5.96a7ep+4)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x5.96a7e8p+4)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d378p+12)));
    try std.testing.expectEqual(0x8.378d9p+124, cosh(@as(f32, 0x5.8bfe6p+4)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c6788cp+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c67888p+8)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f32), cosh(@as(f32, -0x2.c5d378p+12)));

    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1.4b705d1e5d6a8p+0, cosh(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x8.e6960966c8d2p+1020, cosh(@as(f64, 0x2.c5e3bp+8)));
    // try std.testing.expectEqual(0x8.e6726f55d7888p+1020, cosh(@as(f64, 0x2.c5e3acp+8)));
    try std.testing.expectEqual(0x8.e679c177a00cp+1020, cosh(@as(f64, 0x2.c5e3acd2922a6p+8)));
    // try std.testing.expectEqual(0x8.e6726f55d7888p+1020, cosh(@as(f64, -0x2.c5e3acp+8)));
    try std.testing.expectEqual(0x8.e6960966c8d2p+1020, cosh(@as(f64, -0x2.c5e3bp+8)));
    try std.testing.expectEqual(0x8.e679c177a00cp+1020, cosh(@as(f64, -0x2.c5e3acd2922a6p+8)));
    try std.testing.expectEqual(0x6.ad6b6e710d8p+28, cosh(@as(f64, 0x1.6p+4)));
    try std.testing.expectEqual(0x1.226af33b1fdc1p+32, cosh(@as(f64, 0x1.7p+4)));
    try std.testing.expectEqual(0x3.156ff6a8ebf6ep+32, cosh(@as(f64, 0x1.8p+4)));
    try std.testing.expectEqual(0x1.002000aaac16cp+0, cosh(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0x1.00000800000abp+0, cosh(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0x1.00000002p+0, cosh(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0x1.00000000008p+0, cosh(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0x1.0000000000002p+0, cosh(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1.8b07551d9f55p+0, cosh(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x8.c881f20405a28p+68, cosh(@as(f64, 0x3.2p+4)));
    try std.testing.expectEqual(0xa.a7179c1019ae8p+12, cosh(@as(f64, -0xb.60713p+0)));
    try std.testing.expectEqual(0x1.68b8dc5c49a89p+4, cosh(@as(f64, -0x3.cee48p+0)));
    try std.testing.expectEqual(0x9.ad526ad56446p+0, cosh(@as(f64, 0x2.f5d128p+0)));
    try std.testing.expectEqual(0x3.89993d3ed803p+16, cosh(@as(f64, -0xd.0c03p+0)));
    try std.testing.expectEqual(0x1.074e5452941d5p+0, cosh(@as(f64, -0x3.d04328p-4)));
    try std.testing.expectEqual(0x1.074e5461fa3e1p+0, cosh(@as(f64, -0x3.d0432cp-4)));
    try std.testing.expectEqual(0x1.074e54544d14dp+0, cosh(@as(f64, -0x3.d04328728b72cp-4)));
    try std.testing.expectEqual(0x7.d716115677b78p+28, cosh(@as(f64, 0x1.629188p+4)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x1p-56)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, -0x1p-56)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, 0x1p-72)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f64, -0x1p-72)));
    try std.testing.expectEqual(0xf.fffec1f47394p+124, cosh(@as(f64, 0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48ep+128, cosh(@as(f64, 0x5.96a7e8p+4)));
    try std.testing.expectEqual(0xf.fffec1f47394p+124, cosh(@as(f64, -0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48ep+128, cosh(@as(f64, -0x5.96a7e8p+4)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c679d4p+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, cosh(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, cosh(@as(f64, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c679d4p+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, cosh(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c679d1f73f0fcp+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, cosh(@as(f64, -0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, cosh(@as(f64, -0x2.c679d1f73f0fap+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, cosh(@as(f64, -0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c679d4p+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c679d4p+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, cosh(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, cosh(@as(f64, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c679d4p+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, cosh(@as(f64, 0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, cosh(@as(f64, 0x2.c679d1f73f0fap+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, cosh(@as(f64, -0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, cosh(@as(f64, -0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c679d1f73f0fcp+8)));
    // try std.testing.expectEqual(0xf.ffe08c2deedp+1020, cosh(@as(f64, -0x2.c679dp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d8p+1020, cosh(@as(f64, -0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0x8.378d97e8a9838p+124, cosh(@as(f64, 0x5.8bfe6p+4)));
    try std.testing.expectEqual(0xf.ebad7efd1e068p+1020, cosh(@as(f64, 0x2.c6788cp+8)));
    // try std.testing.expectEqual(0xf.eb6dd0c67ed4p+1020, cosh(@as(f64, 0x2.c67888p+8)));
    // try std.testing.expectEqual(0xf.eb9d774858638p+1020, cosh(@as(f64, 0x2.c6788afe3ccd6p+8)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d376167f406p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, 0x2.c5d376167f404p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d378p+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d376eefcd4ap+12)));
    try std.testing.expectEqual(std.math.inf(f64), cosh(@as(f64, -0x2.c5d376eefcd4cp+12)));

    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1.4b705d1e5d6a787ap+0, cosh(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x8.e6960966c8d230bp+1020, cosh(@as(f80, 0x2.c5e3bp+8)));
    try std.testing.expectEqual(0x8.e6726f55d788682p+1020, cosh(@as(f80, 0x2.c5e3acp+8)));
    try std.testing.expectEqual(0x8.e679c177a00bfb6p+1020, cosh(@as(f80, 0x2.c5e3acd2922a6p+8)));
    try std.testing.expectEqual(0x8.e6726f55d788682p+1020, cosh(@as(f80, -0x2.c5e3acp+8)));
    try std.testing.expectEqual(0x8.e6960966c8d230bp+1020, cosh(@as(f80, -0x2.c5e3bp+8)));
    try std.testing.expectEqual(0x8.e679c177a00bfb6p+1020, cosh(@as(f80, -0x2.c5e3acd2922a6p+8)));
    try std.testing.expectEqual(0x6.ad6b6e710d7fe078p+28, cosh(@as(f80, 0x1.6p+4)));
    try std.testing.expectEqual(0x1.226af33b1fdc0a58p+32, cosh(@as(f80, 0x1.7p+4)));
    try std.testing.expectEqual(0x3.156ff6a8ebf6e67p+32, cosh(@as(f80, 0x1.8p+4)));
    try std.testing.expectEqual(0x1.002000aaac16c30cp+0, cosh(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x1.00000800000aaaaap+0, cosh(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x1.00000002p+0, cosh(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0x1.00000000008p+0, cosh(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0x1.0000000000002p+0, cosh(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x1.0000000000000008p+0, cosh(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x1.8b07551d9f5504c2p+0, cosh(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x8.c881f20405a2b32p+68, cosh(@as(f80, 0x3.2p+4)));
    try std.testing.expectEqual(0xa.a7179c1019ae57ep+12, cosh(@as(f80, -0xb.60713p+0)));
    try std.testing.expectEqual(0x1.68b8dc5c49a88f56p+4, cosh(@as(f80, -0x3.cee48p+0)));
    try std.testing.expectEqual(0x9.ad526ad564464p+0, cosh(@as(f80, 0x2.f5d128p+0)));
    try std.testing.expectEqual(0x3.89993d3ed8030b98p+16, cosh(@as(f80, -0xd.0c03p+0)));
    try std.testing.expectEqual(0x1.074e5452941d4ccap+0, cosh(@as(f80, -0x3.d04328p-4)));
    try std.testing.expectEqual(0x1.074e5461fa3e0c5ep+0, cosh(@as(f80, -0x3.d0432cp-4)));
    try std.testing.expectEqual(0x1.074e54544d14c8p+0, cosh(@as(f80, -0x3.d04328728b72cp-4)));
    try std.testing.expectEqual(0x7.d716115677b7982p+28, cosh(@as(f80, 0x1.629188p+4)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x1p-56)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x1p-56)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, 0x1p-72)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f80, -0x1p-72)));
    try std.testing.expectEqual(0xf.fffec1f473940d2p+124, cosh(@as(f80, 0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48e748p+128, cosh(@as(f80, 0x5.96a7e8p+4)));
    try std.testing.expectEqual(0xf.fffec1f473940d2p+124, cosh(@as(f80, -0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48e748p+128, cosh(@as(f80, -0x5.96a7e8p+4)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, cosh(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, cosh(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, cosh(@as(f80, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, cosh(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, cosh(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, cosh(@as(f80, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, cosh(@as(f80, -0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, cosh(@as(f80, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, cosh(@as(f80, -0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, cosh(@as(f80, -0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, cosh(@as(f80, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, cosh(@as(f80, -0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, cosh(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, cosh(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, cosh(@as(f80, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, cosh(@as(f80, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2dp+1020, cosh(@as(f80, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2dp+1020, cosh(@as(f80, 0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, cosh(@as(f80, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, cosh(@as(f80, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, cosh(@as(f80, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, cosh(@as(f80, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2dp+1020, cosh(@as(f80, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2dp+1020, cosh(@as(f80, 0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, cosh(@as(f80, -0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, cosh(@as(f80, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, cosh(@as(f80, -0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, cosh(@as(f80, -0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2dp+1020, cosh(@as(f80, -0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2dp+1020, cosh(@as(f80, -0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b1p+1020, cosh(@as(f80, -0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1cp+1024, cosh(@as(f80, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72dp+1020, cosh(@as(f80, -0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.000000000009d72cp+1024, cosh(@as(f80, -0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2dp+1020, cosh(@as(f80, -0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2dp+1020, cosh(@as(f80, -0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, cosh(@as(f80, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, cosh(@as(f80, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, cosh(@as(f80, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, cosh(@as(f80, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, cosh(@as(f80, -0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, cosh(@as(f80, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, cosh(@as(f80, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, cosh(@as(f80, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, cosh(@as(f80, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, cosh(@as(f80, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, cosh(@as(f80, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, cosh(@as(f80, -0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fffffffffc593dbp+16380, cosh(@as(f80, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(0xf.fffffffffffd3dbp+16380, cosh(@as(f80, -0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(0x8.378d97e8a9838b8p+124, cosh(@as(f80, 0x5.8bfe6p+4)));
    try std.testing.expectEqual(0xf.ebad7efd1e065ep+1020, cosh(@as(f80, 0x2.c6788cp+8)));
    try std.testing.expectEqual(0xf.eb6dd0c67ed40c9p+1020, cosh(@as(f80, 0x2.c67888p+8)));
    try std.testing.expectEqual(0xf.eb9d7748586375dp+1020, cosh(@as(f80, 0x2.c6788afe3ccd6p+8)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, 0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, 0x2.c5d374p+12)));
    try std.testing.expectEqual(0xf.ff15bf3871a7576p+16380, cosh(@as(f80, 0x2.c5d376167f406p+12)));
    try std.testing.expectEqual(0xf.ff15bf3851a92bep+16380, cosh(@as(f80, 0x2.c5d376167f404p+12)));
    try std.testing.expectEqual(0xf.ff15bf38649c166p+16380, cosh(@as(f80, 0x2.c5d376167f4052f4p+12)));
    try std.testing.expectEqual(0xf.fcff8165c0f3207p+16380, cosh(@as(f80, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f80), cosh(@as(f80, -0x2.c5d378p+12)));
    try std.testing.expectEqual(0xf.ffee36237fd43a3p+16380, cosh(@as(f80, -0x2.c5d376eefcd4ap+12)));
    try std.testing.expectEqual(0xf.ffee36239fd4169p+16380, cosh(@as(f80, -0x2.c5d376eefcd4cp+12)));
    try std.testing.expectEqual(0xf.ffee36239bbc1b2p+16380, cosh(@as(f80, -0x2.c5d376eefcd4bbe8p+12)));
    try std.testing.expectEqual(0xf.ffee36239bc01b2p+16380, cosh(@as(f80, -0x2.c5d376eefcd4bbecp+12)));

    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1.4b705d1e5d6a787aa2de94beca32p+0, cosh(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x8.e6960966c8d230b719596be4b88p+1020, cosh(@as(f128, 0x2.c5e3bp+8)));
    try std.testing.expectEqual(0x8.e6726f55d78868187eba9eac383p+1020, cosh(@as(f128, 0x2.c5e3acp+8)));
    try std.testing.expectEqual(0x8.e679c177a00bfb5aec6fa96b5868p+1020, cosh(@as(f128, 0x2.c5e3acd2922a6p+8)));
    try std.testing.expectEqual(0x8.e6726f55d78868187eba9eac383p+1020, cosh(@as(f128, -0x2.c5e3acp+8)));
    try std.testing.expectEqual(0x8.e6960966c8d230b719596be4b88p+1020, cosh(@as(f128, -0x2.c5e3bp+8)));
    try std.testing.expectEqual(0x8.e679c177a00bfb5aec6fa96b5868p+1020, cosh(@as(f128, -0x2.c5e3acd2922a6p+8)));
    // try std.testing.expectEqual(0x6.ad6b6e710d7fe07862bf28dca0a4p+28, cosh(@as(f128, 0x1.6p+4)));
    try std.testing.expectEqual(0x1.226af33b1fdc0a57bd4b4ab2311bp+32, cosh(@as(f128, 0x1.7p+4)));
    try std.testing.expectEqual(0x3.156ff6a8ebf6e66f4935281c5fbp+32, cosh(@as(f128, 0x1.8p+4)));
    try std.testing.expectEqual(0x1.002000aaac16c30c31eaf1bbb19p+0, cosh(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(0x1.00000800000aaaaab05b05b1fb2p+0, cosh(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x1.0000000200000000aaaaaaaac16cp+0, cosh(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0x1.000000000080000000000aaaaaabp+0, cosh(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0x1.00000000000020000000000000abp+0, cosh(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x1.0000000000000008p+0, cosh(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0x1.000000000000000002p+0, cosh(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0x1.000000000000000000008p+0, cosh(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0x1.00000000000000000000002p+0, cosh(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0x1.00000000000000000000000008p+0, cosh(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x1p-10000)));
    // try std.testing.expectEqual(0x1.8b07551d9f5504c2bd28100196a5p+0, cosh(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x8.c881f20405a2b326bba067c62ec8p+68, cosh(@as(f128, 0x3.2p+4)));
    try std.testing.expectEqual(0xa.a7179c1019ae57dfcdfc8ae2c12p+12, cosh(@as(f128, -0xb.60713p+0)));
    try std.testing.expectEqual(0x1.68b8dc5c49a88f56145c6a6eb1fbp+4, cosh(@as(f128, -0x3.cee48p+0)));
    try std.testing.expectEqual(0x9.ad526ad564463ffecc391e2180a8p+0, cosh(@as(f128, 0x2.f5d128p+0)));
    try std.testing.expectEqual(0x3.89993d3ed8030b962f4a1d333f74p+16, cosh(@as(f128, -0xd.0c03p+0)));
    try std.testing.expectEqual(0x1.074e5452941d4cca93e217a9d915p+0, cosh(@as(f128, -0x3.d04328p-4)));
    try std.testing.expectEqual(0x1.074e5461fa3e0c5d7d941a2999d5p+0, cosh(@as(f128, -0x3.d0432cp-4)));
    try std.testing.expectEqual(0x1.074e54544d14c800f66940138bb9p+0, cosh(@as(f128, -0x3.d04328728b72cp-4)));
    try std.testing.expectEqual(0x7.d716115677b7981c1502cadb3d14p+28, cosh(@as(f128, 0x1.629188p+4)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x4p-16496)));
    // try std.testing.expectEqual(0x1.0000000000000000000000000001p+0, cosh(@as(f128, 0x1p-56)));
    // try std.testing.expectEqual(0x1.0000000000000000000000000001p+0, cosh(@as(f128, -0x1p-56)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, 0x1p-72)));
    try std.testing.expectEqual(0x1p+0, cosh(@as(f128, -0x1p-72)));
    try std.testing.expectEqual(0xf.fffec1f473940d22f2195eac65ep+124, cosh(@as(f128, 0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48e7480e07d1c02e7cp+128, cosh(@as(f128, 0x5.96a7e8p+4)));
    try std.testing.expectEqual(0xf.fffec1f473940d22f2195eac65ep+124, cosh(@as(f128, -0x5.96a7ep+4)));
    try std.testing.expectEqual(0x1.00006c1f5d48e7480e07d1c02e7cp+128, cosh(@as(f128, -0x5.96a7e8p+4)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, cosh(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, cosh(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, cosh(@as(f128, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, cosh(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, cosh(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, cosh(@as(f128, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, cosh(@as(f128, -0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, cosh(@as(f128, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, cosh(@as(f128, -0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, cosh(@as(f128, -0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, cosh(@as(f128, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, cosh(@as(f128, -0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, cosh(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, cosh(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, cosh(@as(f128, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, cosh(@as(f128, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2ca74dec5830328p+1020, cosh(@as(f128, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2ca74dec58303ep+1020, cosh(@as(f128, 0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffbffffffffff303a8p+1020, cosh(@as(f128, 0x2.c679d1f73f0fb624d358b213a7p+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, cosh(@as(f128, 0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, cosh(@as(f128, 0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, cosh(@as(f128, 0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, cosh(@as(f128, 0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2ca74dec5830328p+1020, cosh(@as(f128, 0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2ca74dec58303ep+1020, cosh(@as(f128, 0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffc0000000000303a8p+1020, cosh(@as(f128, 0x2.c679d1f73f0fb624d358b213a8p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, cosh(@as(f128, -0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, cosh(@as(f128, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, cosh(@as(f128, -0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, cosh(@as(f128, -0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2ca74dec58303ep+1020, cosh(@as(f128, -0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2ca74dec5830328p+1020, cosh(@as(f128, -0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffbffffffffff303a8p+1020, cosh(@as(f128, -0x2.c679d1f73f0fb624d358b213a7p+8)));
    try std.testing.expectEqual(0xf.ffe08c2deed02b0e9ba9e9c42178p+1020, cosh(@as(f128, -0x2.c679dp+8)));
    try std.testing.expectEqual(0x1.000208c301f36f1c494de034e38p+1024, cosh(@as(f128, -0x2.c679d4p+8)));
    try std.testing.expectEqual(0xf.fffffffffe9d72ca74ded4db59d8p+1020, cosh(@as(f128, -0x2.c679d1f73f0fap+8)));
    try std.testing.expectEqual(0x1.000000000009d72ca74dec889b32p+1024, cosh(@as(f128, -0x2.c679d1f73f0fcp+8)));
    try std.testing.expectEqual(0xf.ffffffffffffb2ca74dec58303ep+1020, cosh(@as(f128, -0x2.c679d1f73f0fb624p+8)));
    try std.testing.expectEqual(0xf.fffffffffffff2ca74dec5830328p+1020, cosh(@as(f128, -0x2.c679d1f73f0fb628p+8)));
    try std.testing.expectEqual(0xf.ffffffffffffc0000000000303a8p+1020, cosh(@as(f128, -0x2.c679d1f73f0fb624d358b213a8p+8)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, cosh(@as(f128, 0x2.c5d37700c6bbp+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, cosh(@as(f128, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, cosh(@as(f128, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb03a8p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, cosh(@as(f128, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, cosh(@as(f128, -0x2.c5d37700c6bb03a4p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, cosh(@as(f128, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb2p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, cosh(@as(f128, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb03a8p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, cosh(@as(f128, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffe61p+16380, cosh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b494cp+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b4ap+12)));
    // try std.testing.expectEqual(0xf.ffffffffffffffffffffffb3e61p+16380, cosh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b49p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, 0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, cosh(@as(f128, 0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb03a8p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, cosh(@as(f128, 0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b494ep+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b4ap+12)));
    // try std.testing.expectEqual(0xf.ffffffffffffffffffffffb3e61p+16380, cosh(@as(f128, 0x2.c5d37700c6bb03a6c24b6c9b49p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, cosh(@as(f128, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, cosh(@as(f128, -0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffe61p+16380, cosh(@as(f128, -0x2.c5d37700c6bb03a6c24b6c9b494cp+12)));
    // try std.testing.expectEqual(0xf.ffffffffffffffffffffffb3e61p+16380, cosh(@as(f128, -0x2.c5d37700c6bb03a6c24b6c9b49p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb03a6c24b6c9b4ap+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fffffffffc593db49365215d58ap+16380, cosh(@as(f128, -0x2.c5d37700c6bbp+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb2p+12)));
    // try std.testing.expectEqual(0xf.fffffffffffd3db49364b6b422f8p+16380, cosh(@as(f128, -0x2.c5d37700c6bb03a4p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb03a8p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb03a6c24b6c9b494ep+12)));
    // try std.testing.expectEqual(0xf.ffffffffffffffffffffffb3e61p+16380, cosh(@as(f128, -0x2.c5d37700c6bb03a6c24b6c9b49p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d37700c6bb03a6c24b6c9b4ap+12)));
    try std.testing.expectEqual(0x8.378d97e8a9838b8164de61b93a68p+124, cosh(@as(f128, 0x5.8bfe6p+4)));
    try std.testing.expectEqual(0xf.ebad7efd1e065dfa4889d66d8e5p+1020, cosh(@as(f128, 0x2.c6788cp+8)));
    try std.testing.expectEqual(0xf.eb6dd0c67ed40c8e47a528f28b6p+1020, cosh(@as(f128, 0x2.c67888p+8)));
    try std.testing.expectEqual(0xf.eb9d7748586375cf28c2e4264d88p+1020, cosh(@as(f128, 0x2.c6788afe3ccd6p+8)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, 0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, 0x2.c5d374p+12)));
    // try std.testing.expectEqual(0xf.ff15bf3871a75761db61506a9bbp+16380, cosh(@as(f128, 0x2.c5d376167f406p+12)));
    // try std.testing.expectEqual(0xf.ff15bf3851a92be36a9dffe75668p+16380, cosh(@as(f128, 0x2.c5d376167f404p+12)));
    try std.testing.expectEqual(0xf.ff15bf38649c16662e1ff4acb94p+16380, cosh(@as(f128, 0x2.c5d376167f4052f4p+12)));
    // try std.testing.expectEqual(0xf.fcff8165c0f3206f5cab3921788p+16380, cosh(@as(f128, -0x2.c5d374p+12)));
    try std.testing.expectEqual(std.math.inf(f128), cosh(@as(f128, -0x2.c5d378p+12)));
    // try std.testing.expectEqual(0xf.ffee36237fd43a2b15e5b20b6d68p+16380, cosh(@as(f128, -0x2.c5d376eefcd4ap+12)));
    // try std.testing.expectEqual(0xf.ffee36239fd416975d055a5c2fep+16380, cosh(@as(f128, -0x2.c5d376eefcd4cp+12)));
    try std.testing.expectEqual(0xf.ffee36239bbc1b2482e87ba9d31p+16380, cosh(@as(f128, -0x2.c5d376eefcd4bbe8p+12)));
    // try std.testing.expectEqual(0xf.ffee36239bc01b201071629959ep+16380, cosh(@as(f128, -0x2.c5d376eefcd4bbecp+12)));
    // try std.testing.expectEqual(0xf.ffee36239bbf1b257fe2a0ad4a38p+16380, cosh(@as(f128, -0x2.c5d376eefcd4bbeb000452d84662p+12)));
    // try std.testing.expectEqual(0xf.ffee36239bbf1b257fe2a04b4aa8p+16380, cosh(@as(f128, -0x2.c5d376eefcd4bbeb000452d846p+12)));
    // try std.testing.expectEqual(0xf.ffee36239bbf1b257fe2a14b4988p+16380, cosh(@as(f128, -0x2.c5d376eefcd4bbeb000452d847p+12)));
}
