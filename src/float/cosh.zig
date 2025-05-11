const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
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
                return @mulAdd(f32, float.abs(x), 0x1p-25, 1);
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

            const t: f64 = float.expm1(float.abs(x));
            const w: f64 = 1 + t;
            return 1 + (t * t) / (w + w);
        }

        // |x| in [0.5*ln2,22], return (exp(|x|)+1/exp(|x|)/2;
        const t: f64 = float.exp(float.abs(x));
        return 0.5 * t + 0.5 / t;
    }

    // |x| in [22, log(maxdouble)] return 0.5*exp(|x|)
    if (ix < 0x40862e42)
        return 0.5 * float.exp(float.abs(x));

    // |x| in [log(maxdouble), overflowthresold]
    var fix: i64 = undefined;
    dbl64.extractWords64(&fix, x);
    fix &= 0x7fffffffffffffff;
    if (fix <= 0x408633ce8fb9f87d) {
        const w: f64 = float.exp(0.5 * float.abs(x));
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

        const t: f128 = float.expm1(@as(f128, @bitCast(u)));
        const w: f128 = 1 + t;
        return 1 + (t * t) / (w + w);
    }

    // |x| in [0.5*ln2,40], return (exp(|x|)+1/exp(|x|)/2;
    if (ex < 0x40044000) {
        const t: f128 = float.exp(@as(f128, @bitCast(u)));
        return 0.5 * t + 0.5 / t;
    }

    // |x| in [22, ln(maxdouble)] return 0.5*exp(|x|)
    if (ex <= 0x400c62e3) // 11356.375
        return 0.5 * float.exp(@as(f128, @bitCast(u)));

    // |x| in [log(maxdouble), overflowthresold]
    if (@as(f128, @bitCast(u)) <= ovf_thresh) {
        const w: f128 = float.exp(0.5 * @as(f128, @bitCast(u)));
        const t: f128 = 0.5 * w;
        return t * w;
    }

    // |x| > overflowthresold, cosh(x) overflow
    return huge * huge;
}
