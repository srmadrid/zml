const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const roundeven = @import("roundeven.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn sinh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.sinh: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return scast(f16, sinh32(scast(f32, x))),
        f32 => {
            // glibc/sysdeps/ieee754/flt-32/e_sinhf.c
            return sinh32(scast(f32, x));
        },
        f64 => {
            // glibc/sysdeps/ieee754/dbl-64/e_sinh.c
            return sinh64(scast(f64, x));
        },
        f80 => return scast(f80, sinh128(scast(f128, x))),
        f128 => {
            // glibc/sysdeps/ieee754/ldbl-128/e_sinhl.c
            return sinh128(scast(f128, x));
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
    const z: f64 = scast(f64, x);
    const ux: u32 = @as(u32, @bitCast(x)) << 1;
    if (ux > 0x8565a9f8) { // |x| >~ 89.4
        @branchHint(.unlikely);
        const sgn: f32 = float.copysign(@as(f32, 2), x);
        if (ux >= 0xff000000) {
            if ((ux << 8) != 0)
                return x + x; // nan

            return float.copysign(std.math.inf(f32), x); // +-inf
        }
        return sgn * 0x1.fffffep127;
    }

    if (ux < 0x7c000000) { // |x| < 0.125
        @branchHint(.unlikely);
        if (ux <= 0x74250bfe) { // |x| <= 0x1.250bfep-11
            @branchHint(.unlikely);
            if (ux < 0x66000000) { // |x| < 0x1p-24
                @branchHint(.unlikely);
                return @mulAdd(f32, x, float.abs(x), x);
            }

            if (st.uarg == ux) {
                @branchHint(.unlikely);
                const sgn: f32 = float.copysign(@as(f32, 1), x);
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
        return scast(f32, z + (z2 * z) * ((cp[0] + z2 * cp[1]) + z4 * (cp[2] + z2 * (cp[3]))));
    }

    const a: f64 = iln2 * z;
    const ia: f64 = roundeven.roundeven_finite(a);
    var h: f64 = a - ia;
    var h2: f64 = h * h;
    const jp: i64 = @bitCast(ia + 0x1.8p52);
    const jm: i64 = -jp;
    const sp: f64 = @bitCast(scast(i64, tb[@intCast(jp & 31)]) + ((jp >> 5) << 52));
    const sm: f64 = @bitCast(scast(i64, tb[@intCast(jm & 31)]) + ((jm >> 5) << 52));
    var te: f64 = c[0] + h2 * c[2];
    var to: f64 = (c[1] + h2 * c[3]);
    const rp: f64 = sp * (te + h * to);
    const rm: f64 = sm * (te - h * to);
    var r: f64 = rp - rm;
    var ub: f32 = scast(f32, r);
    const lb: f32 = scast(f32, r - 1.52e-10 * r);
    if (ub != lb) {
        @branchHint(.unlikely);
        const iln2h: f64 = 0x1.7154765p+5;
        const iln2l: f64 = 0x1.5c17f0bbbe88p-26;
        h = (iln2h * z - ia) + iln2l * z;
        h2 = h * h;
        te = ch[0] + h2 * ch[2] + (h2 * h2) * (ch[4] + h2 * ch[6]);
        to = ch[1] + h2 * (ch[3] + h2 * ch[5]);
        r = sp * (te + h * to) - sm * (te - h * to);
        ub = scast(f32, r);
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
            if (float.abs(x) < std.math.floatMin(f64)) {
                const vx: f64 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            if (shuge + x > 1)
                return x;

            // sinh(tiny) = tiny with inexact
        }

        const t: f64 = float.expm1(float.abs(x));
        if (ix < 0x3ff00000)
            return h * (2 * t - t * t / (t + 1));

        return h * (t + t / (t + 1));
    }

    // |x| in [22, log(maxdouble)] return 0.5*exp(|x|)
    if (ix < 0x40862e42)
        return h * float.exp(float.abs(x));

    // |x| in [log(maxdouble), overflowthresold]
    var lx: u32 = undefined;
    dbl64.getLowWord(&lx, x);
    if (ix < 0x408633ce or ((ix == 0x408633ce) and (lx <= 0x8fb9f87d))) {
        const w: f64 = float.exp(0.5 * float.abs(x));
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
            if (float.abs(x) < std.math.floatMin(f128)) {
                const vx: f128 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            if (shuge + x > 1)
                return x; // sinh(tiny) = tiny with inexact
        }
        const t: f128 = float.expm1(@as(f128, @bitCast(u)));
        if (ix < 0x3fff0000)
            return h * (2 * t - t * t / (t + 1));

        return h * (t + t / (t + 1));
    }

    // |x| in [40, log(maxdouble)] return 0.5*exp(|x|)
    if (ix <= 0x400c62e3) // 11356.375
        return h * float.exp(@as(f128, @bitCast(u)));

    // |x| in [log(maxdouble), overflowthreshold]
    // Overflow threshold is log(2 * maxdouble).
    if (@as(f128, @bitCast(u)) <= ovf_thresh) {
        const w: f128 = float.exp(0.5 * @as(f128, @bitCast(u)));
        const t: f128 = h * w;
        return t * w;
    }

    // |x| > overflowthreshold, sinhl(x) overflow
    return x * shuge;
}
