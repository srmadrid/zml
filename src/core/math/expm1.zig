const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const roundeven = @import("roundeven.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn expm1(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return expm1(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, expm1_32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_expm1f.c
                    return expm1_32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_expm1.c
                    return expm1_64(x);
                },
                f80 => return cast(f80, expm1_128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_expm1l.c
                    return expm1_128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn expm1_32(x: f32) f32 {
    const c: [4]f64 = .{ 1, 0x1.62e42fef4c4e7p-6, 0x1.ebfd1b232f475p-13, 0x1.c6b19384ecd93p-20 };
    const ch: [6]f64 = .{
        0x1.62e42fefa39efp-6,  0x1.ebfbdff82c58fp-13, 0x1.c6b08d702e0edp-20,
        0x1.3b2ab6fb92e5ep-27, 0x1.5d886e6d54203p-35, 0x1.430976b8ce6efp-43,
    };
    const td: [32]f64 = .{
        0x1p+0,               0x1.059b0d3158574p+0, 0x1.0b5586cf9890fp+0,
        0x1.11301d0125b51p+0, 0x1.172b83c7d517bp+0, 0x1.1d4873168b9aap+0,
        0x1.2387a6e756238p+0, 0x1.29e9df51fdee1p+0, 0x1.306fe0a31b715p+0,
        0x1.371a7373aa9cbp+0, 0x1.3dea64c123422p+0, 0x1.44e086061892dp+0,
        0x1.4bfdad5362a27p+0, 0x1.5342b569d4f82p+0, 0x1.5ab07dd485429p+0,
        0x1.6247eb03a5585p+0, 0x1.6a09e667f3bcdp+0, 0x1.71f75e8ec5f74p+0,
        0x1.7a11473eb0187p+0, 0x1.82589994cce13p+0, 0x1.8ace5422aa0dbp+0,
        0x1.93737b0cdc5e5p+0, 0x1.9c49182a3f09p+0,  0x1.a5503b23e255dp+0,
        0x1.ae89f995ad3adp+0, 0x1.b7f76f2fb5e47p+0, 0x1.c199bdd85529cp+0,
        0x1.cb720dcef9069p+0, 0x1.d5818dcfba487p+0, 0x1.dfc97337b9b5fp+0,
        0x1.ea4afa2a490dap+0, 0x1.f50765b6e454p+0,
    };
    const iln2: f64 = 0x1.71547652b82fep+5;
    const big: f64 = 0x1.8p52;
    const z: f64 = cast(f64, x, .{});
    const ux: u32 = @bitCast(x);
    const ax: u32 = ux << 1;
    if (ax < 0x7c400000) { // |x| < 0.15625
        @branchHint(.likely);
        if (ax < 0x676a09e8) { // |x| < 0x1.6a09e8p-24
            @branchHint(.unlikely);
            if (ax == 0) {
                @branchHint(.unlikely);
                return x; // x = +-0
            }

            return @mulAdd(f32, math.abs(x), 0x1p-25, x);
        }

        const b: [8]f64 = .{
            0x1.fffffffffffc2p-2,  0x1.55555555555fep-3,  0x1.555555559767fp-5,
            0x1.1111111098dc1p-7,  0x1.6c16bca988aa9p-10, 0x1.a01a07658483fp-13,
            0x1.a05b04d2c3503p-16, 0x1.71de3a960b5e3p-19,
        };
        const z2: f64 = z * z;
        const z4: f64 = z2 * z2;
        const r: f64 = z + z2 * ((b[0] + z * b[1]) + z2 * (b[2] + z * b[3]) + z4 * ((b[4] + z * b[5]) + z2 * (b[6] + z * b[7])));
        return cast(f32, r, .{});
    }

    if (ax >= 0x8562e430) { // |x| > 88.72
        @branchHint(.unlikely);
        if (ax > (0xff << 24))
            return x + x; // nan
        if ((ux >> 31) != 0) { // x < 0
            @branchHint(.unlikely);
            if (ax == (0xff << 24))
                return -1;

            return -1 + 0x1p-26;
        }

        if (ax == (0xff << 24))
            return std.math.inf(f32);

        return 0x1p97 * 0x1p97;
    }

    const a: f64 = iln2 * z;
    const ia: f64 = roundeven.roundeven_finite(a);
    var h: f64 = a - ia;
    var h2: f64 = h * h;
    const u: u64 = @bitCast(ia + big);
    const c2: f64 = c[2] + h * c[3];
    const c0: f64 = c[0] + h * c[1];
    const sv: f64 = @bitCast(@as(u64, @bitCast(td[u & 0x1f])) +% ((u >> 5) << 52));
    var r: f64 = (c0 + h2 * c2) * sv - 1.0;
    var ub: f32 = cast(f32, r, .{});
    const lb: f32 = cast(f32, r - sv * 0x1.3b3p-33, .{});
    if (ub != lb) {
        @branchHint(.unlikely);
        if (ux > 0xc18aa123) { // x < -17.32
            @branchHint(.unlikely);
            return -1.0 + 0x1p-26;
        }

        const iln2h: f64 = 0x1.7154765p+5;
        const iln2l: f64 = 0x1.5c17f0bbbe88p-26;
        const s: f64 = sv;
        h = (iln2h * z - ia) + iln2l * z;
        h2 = h * h;
        const w: f64 = s * h;
        r = (s - 1) + w * ((ch[0] + h * ch[1]) + h2 * ((ch[2] + h * ch[3]) + h2 * (ch[4] + h * ch[5])));
        ub = cast(f32, r, .{});
    }
    return ub;
}

fn expm1_64(x: f64) f64 {
    const huge: f64 = 1.0e+300;
    const tiny: f64 = 1.0e-300;
    const o_threshold: f64 = 7.09782712893383973096e+02; // 0x40862e42, 0xfefa39ef
    const ln2_hi: f64 = 6.93147180369123816490e-01; // 0x3fe62e42, 0xfee00000
    const ln2_lo: f64 = 1.90821492927058770002e-10; // 0x3dea39ef, 0x35793c76
    const invln2: f64 = 1.44269504088896338700e+00; // 0x3ff71547, 0x652b82fe
    // scaled coefficients related to expm1
    const Q: [6]f64 = .{
        1.0, -3.33333333333331316428e-02, // bfa11111 111110f4
        1.58730158725481460165e-03, // 3f5a01a0 19fe5585
        -7.93650757867487942473e-05, // bf14ce19 9eaadbb7
        4.00821782732936239552e-06, // 3ed0cfca 86e65239
        -2.01099218183624371326e-07, // be8afdb7 6e09c32d
    };

    var hx: u32 = undefined;
    dbl64.getHighWord(&hx, x);
    const xsb: u32 = hx & 0x80000000; // sign bit of x
    var y = x;
    if (xsb != 0)
        y = -x; // y = |x|

    hx &= 0x7fffffff; // high word of |x|
    // filter out huge and non-finite argument
    if (hx >= 0x4043687a) { // if |x|>=56*ln2
        if (hx >= 0x40862e42) { // if |x|>=709.78...
            if (hx >= 0x7ff00000) {
                var low: u32 = undefined;
                dbl64.getLowWord(&low, x);
                if (((hx & 0xfffff) | low) != 0) {
                    return x + x; // NaN
                } else {
                    return if (xsb == 0) x else -1; // exp(+-inf)={inf,-1}
                }
            }
            if (x > o_threshold) {
                return huge * huge; // overflow
            }
        }
        if (xsb != 0) { // x < -56*ln2, return -1.0 with inexact
            std.mem.doNotOptimizeAway(x + tiny); // raise inexact
            return tiny - 1; // return -1
        }
    }

    // argument reduction
    var k: i32 = undefined;
    var xx: f64 = x;
    var c: f64 = undefined;
    if (hx > 0x3fd62e42) { // if  |x| > 0.5 ln2
        var hi: f64 = undefined;
        var lo: f64 = undefined;
        if (hx < 0x3ff0a2b2) { // and |x| < 1.5 ln2
            if (xsb == 0) {
                hi = x - ln2_hi;
                lo = ln2_lo;
                k = 1;
            } else {
                hi = x + ln2_hi;
                lo = -ln2_lo;
                k = -1;
            }
        } else {
            k = cast(i32, invln2 * x + (if (xsb == 0) @as(f64, 0.5) else @as(f64, -0.5)), .{});
            const t: f64 = cast(f64, k, .{});
            hi = x - t * ln2_hi; // t*ln2_hi is exact here
            lo = t * ln2_lo;
        }
        xx = hi - lo;
        c = (hi - xx) - lo;
    } else if (hx < 0x3c900000) { // when |x|<2**-54, return x
        if (math.abs(x) < std.math.floatMin(f64)) {
            const vx: f64 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }

        const t: f64 = huge + x; // return x with inexact flags when x!=0
        return x - (t - (huge + x));
    } else k = 0;

    // x is now in primary range
    const hfx: f64 = 0.5 * xx;
    const hxs: f64 = xx * hfx;
    const R1: f64 = 1 + hxs * Q[1];
    const h2: f64 = hxs * hxs;
    const R2: f64 = Q[2] + hxs * Q[3];
    const h4: f64 = h2 * h2;
    const R3: f64 = Q[4] + hxs * Q[5];
    const r1: f64 = R1 + h2 * R2 + h4 * R3;
    var t: f64 = 3 - r1 * hfx;
    var e: f64 = hxs * ((r1 - t) / (6 - xx * t));
    if (k == 0) {
        return xx - (xx * e - hxs); // c is 0
    } else {
        e = (xx * (e - c) - c);
        e -= hxs;
        if (k == -1)
            return 0.5 * (xx - e) - 0.5;

        if (k == 1) {
            if (xx < -0.25) {
                return -2 * (e - (xx + 0.5));
            } else {
                return 1 + 2.0 * (xx - e);
            }
        }
        if (k <= -2 or k > 56) { // suffice to return exp(x)-1
            y = 1 - (e - xx);
            var high: i32 = undefined;
            dbl64.getHighWord(&high, y);
            dbl64.setHighWord(&y, high + (k << 20)); // add k to y's exponent

            return y - 1;
        }

        t = 1;
        if (k < 20) {
            dbl64.setHighWord(&t, 0x3ff00000 - (@as(i32, 0x200000) >> @as(u5, @intCast(k)))); // t=1-2^-k
            y = t - (e - xx);
            var high: i32 = undefined;
            dbl64.getHighWord(&high, y);
            dbl64.setHighWord(&y, high + (k << 20)); // add k to y's exponent
        } else {
            dbl64.setHighWord(&t, ((0x3ff - k) << 20)); // 2^-k
            y = xx - (e + t);
            y += 1;
            var high: i32 = undefined;
            dbl64.getHighWord(&high, y);
            dbl64.setHighWord(&y, high + (k << 20)); // add k to y's exponent
        }
    }

    return y;
}

fn expm1_128(x: f128) f128 {
    // exp(x) - 1 = x + 0.5 x^2 + x^3 P(x)/Q(x)
    // -.5 ln 2  <  x  <  .5 ln 2
    // Theoretical peak relative error = 8.1e-36
    const P0: f128 = 2.943520915569954073888921213330863757240e8;
    const P1: f128 = -5.722847283900608941516165725053359168840e7;
    const P2: f128 = 8.944630806357575461578107295909719817253e6;
    const P3: f128 = -7.212432713558031519943281748462837065308e5;
    const P4: f128 = 4.578962475841642634225390068461943438441e4;
    const P5: f128 = -1.716772506388927649032068540558788106762e3;
    const P6: f128 = 4.401308817383362136048032038528753151144e1;
    const P7: f128 = -4.888737542888633647784737721812546636240e-1;
    const Q0: f128 = 1.766112549341972444333352727998584753865e9;
    const Q1: f128 = -7.848989743695296475743081255027098295771e8;
    const Q2: f128 = 1.615869009634292424463780387327037251069e8;
    const Q3: f128 = -2.019684072836541751428967854947019415698e7;
    const Q4: f128 = 1.682912729190313538934190635536631941751e6;
    const Q5: f128 = -9.615511549171441430850103489315371768998e4;
    const Q6: f128 = 3.697714952261803935521187272204485251835e3;
    const Q7: f128 = -8.802340681794263968892934703309274564037e1;
    // C1 + C2: f128 = ln 2
    const C1: f128 = 6.93145751953125e-1;
    const C2: f128 = 1.428606820309417232121458176568075500134e-6;
    // ln 2^-114
    const minarg: f128 = -7.9018778583833765273564461846232128760607e1;
    const big: f128 = 1e4932;

    // Detect infinity and NaN.
    const u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    var ix: i32 = @bitCast(u.w0);
    const sign: u32 = u.w0 & 0x80000000;
    ix &= 0x7fffffff;
    if (sign == 0 and ix >= 0x40060000) {
        // If num is positive and exp >= 6 use plain exp.
        return math.exp(x);
    }
    if (ix >= 0x7fff0000) {
        // Infinity (which must be negative infinity).
        if (((@as(u32, @bitCast(ix)) & 0xffff) | u.w1 | u.w2 | u.w3) == 0)
            return -1;

        // NaN.  Invalid exception if signaling.
        return x + x;
    }

    // expm1(+- 0) = +- 0.
    if ((ix == 0) and (u.w1 | u.w2 | u.w3) == 0)
        return x;

    // Minimum value.
    if (x < minarg)
        return (4 / big - 1);

    // Avoid internal underflow when result does not underflow, while
    // ensuring underflow (without returning a zero of the wrong sign)
    // when the result does underflow.
    if (math.abs(x) < 0x1p-113) {
        if (math.abs(x) < std.math.floatMin(f128)) {
            const vx: f128 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }

        return x;
    }

    // Express x = ln 2 (k + remainder), remainder not exceeding 1/2.
    var xx: f128 = C1 + C2; // ln 2.
    var px: f128 = math.floor(0.5 + x / xx);
    const k: i32 = cast(i32, px, .{});
    // remainder times ln 2
    var y: f128 = x - px * C1;
    y -= px * C2;

    // Approximate exp(remainder ln 2).
    px = (((((((P7 * y + P6) * y + P5) * y + P4) * y + P3) * y + P2) * y + P1) * y + P0) * y;

    var qx: f128 = (((((((y + Q7) * y + Q6) * y + Q5) * y + Q4) * y + Q3) * y + Q2) * y + Q1) * y + Q0;

    xx = y * y;
    qx = y + (0.5 * xx + xx * px / qx);

    // exp(x) = exp(k ln 2) exp(remainder ln 2) = 2^k exp(remainder ln 2).
    //
    // We have qx = exp(remainder ln 2) - 1, so
    // exp(x) - 1 = 2^k (qx + 1) - 1
    //            = 2^k qx + 2^k - 1.

    px = math.ldexp(1, k);
    y = px * qx + (px - 1.0);
    return y;
}

test expm1 {
    try std.testing.expectEqual(0x0p+0, expm1(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1.b7e152p+0, expm1(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.1df3b6p+0, expm1(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x6.63993p+0, expm1(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(0x1.315e5cp+4, expm1(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(0x3.599204p+4, expm1(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(0x9.369c5p+4, expm1(@as(f32, 0x5p+0)));
    try std.testing.expectEqual(0x5.60977p+12, expm1(@as(f32, 0xap+0)));
    try std.testing.expectEqual(0x3.1e1984p+20, expm1(@as(f32, 0xfp+0)));
    try std.testing.expectEqual(0x1.ceb088p+28, expm1(@as(f32, 0x1.4p+4)));
    try std.testing.expectEqual(0x1.0c3d3ap+36, expm1(@as(f32, 0x1.9p+4)));
    try std.testing.expectEqual(0x9.b8238p+40, expm1(@as(f32, 0x1.ep+4)));
    try std.testing.expectEqual(0x5.a27888p+48, expm1(@as(f32, 0x2.3p+4)));
    try std.testing.expectEqual(0x3.4441a8p+56, expm1(@as(f32, 0x2.8p+4)));
    try std.testing.expectEqual(0x1.19103ep+72, expm1(@as(f32, 0x3.2p+4)));
    try std.testing.expectEqual(0x5.e76f28p+84, expm1(@as(f32, 0x3.cp+4)));
    try std.testing.expectEqual(0x1.fbfd22p+100, expm1(@as(f32, 0x4.6p+4)));
    try std.testing.expectEqual(0xa.abbcep+112, expm1(@as(f32, 0x5p+4)));
    try std.testing.expectEqual(std.math.inf(f32), expm1(@as(f32, 0x5.ap+4)));
    try std.testing.expectEqual(std.math.inf(f32), expm1(@as(f32, 0x6.4p+4)));
    try std.testing.expectEqual(std.math.inf(f32), expm1(@as(f32, 0x7.fp+4)));
    try std.testing.expectEqual(std.math.inf(f32), expm1(@as(f32, 0x1.f4p+8)));
    try std.testing.expectEqual(std.math.inf(f32), expm1(@as(f32, 0x2.c5c4p+12)));
    try std.testing.expectEqual(-0xf.ffd06p-4, expm1(@as(f32, -0xap+0)));
    try std.testing.expectEqual(-0xf.ffffep-4, expm1(@as(f32, -0x1p+4)));
    try std.testing.expectEqual(-0xf.fffffp-4, expm1(@as(f32, -0x1.1p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x1.2p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x2.4p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x2.5p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x2.6p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x2.cp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x2.dp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x2.ep+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x4.9p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x4.ap+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x4.bp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x4.ep+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x4.fp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x5p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x6.4p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x3.e8p+8)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x2.71p+12)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f32), expm1(@as(f32, 0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f32), expm1(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x4.8b5e4p-4, expm1(@as(f32, 0x4p-4)));
    try std.testing.expectEqual(-0x3.8a083p-4, expm1(@as(f32, -0x4p-4)));
    try std.testing.expectEqual(0x4.008008p-12, expm1(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(-0x3.ff800cp-12, expm1(@as(f32, -0x4p-12)));
    try std.testing.expectEqual(0x1.000008p-20, expm1(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(-0xf.ffff8p-24, expm1(@as(f32, -0x1p-20)));
    try std.testing.expectEqual(0x8p-32, expm1(@as(f32, 0x8p-32)));
    try std.testing.expectEqual(-0x8p-32, expm1(@as(f32, -0x8p-32)));
    try std.testing.expectEqual(0x1p-32, expm1(@as(f32, 0x1p-32)));
    try std.testing.expectEqual(-0x1p-32, expm1(@as(f32, -0x1p-32)));
    try std.testing.expectEqual(0x4p-52, expm1(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(-0x4p-52, expm1(@as(f32, -0x4p-52)));
    try std.testing.expectEqual(0x1p-64, expm1(@as(f32, 0x1p-64)));
    try std.testing.expectEqual(-0x1p-64, expm1(@as(f32, -0x1p-64)));
    try std.testing.expectEqual(0x1p-100, expm1(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, expm1(@as(f32, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0xe.4152bp-60, expm1(@as(f32, 0xe.4152bp-60)));
    try std.testing.expectEqual(0xe.4152ap-60, expm1(@as(f32, 0xe.4152ap-60)));
    try std.testing.expectEqual(0x7.ddee38p-4, expm1(@as(f32, 0x6.660248p-4)));
    try std.testing.expectEqual(0x7.ddee3p-4, expm1(@as(f32, 0x6.66024p-4)));
    try std.testing.expectEqual(0x7.830428p-4, expm1(@as(f32, 0x6.289a78p-4)));
    try std.testing.expectEqual(0x7.6f805p-4, expm1(@as(f32, 0x6.1b4d38p-4)));
    try std.testing.expectEqual(0x7.6f804p-4, expm1(@as(f32, 0x6.1b4d3p-4)));
    try std.testing.expectEqual(0x7.412dep-4, expm1(@as(f32, 0x5.fb8dc8p-4)));
    try std.testing.expectEqual(0x7.412dd8p-4, expm1(@as(f32, 0x5.fb8dcp-4)));
    try std.testing.expectEqual(0x3.d9dcf4p-4, expm1(@as(f32, 0x3.735f4cp-4)));
    try std.testing.expectEqual(0x3.d9dcecp-4, expm1(@as(f32, 0x3.735f48p-4)));
    try std.testing.expectEqual(-0xf.fe62cp-4, expm1(@as(f32, -0x7.d6c508p+0)));
    try std.testing.expectEqual(-0xf.fe62cp-4, expm1(@as(f32, -0x7.d6c51p+0)));
    try std.testing.expectEqual(0x1.4aaa8ep+104, expm1(@as(f32, 0x4.857de8p+4)));
    try std.testing.expectEqual(0x7.19268p-4, expm1(@as(f32, 0x5.dfeb68p-4)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x6.a09efp-4, expm1(@as(f32, 0x5.8b912p-4)));
    try std.testing.expectEqual(0x6.a09eep-4, expm1(@as(f32, 0x5.8b9118p-4)));
    try std.testing.expectEqual(0x6.c23b78p-4, expm1(@as(f32, 0x5.a343ep-4)));
    try std.testing.expectEqual(0x6.c23b7p-4, expm1(@as(f32, 0x5.a343d8p-4)));
    try std.testing.expectEqual(0x4p-128, expm1(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-128, expm1(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f32, -0x8p-152)));

    try std.testing.expectEqual(0x0p+0, expm1(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f64, -0x0p+0)));
    // try std.testing.expectEqual(0x1.b7e151628aed3p+0, expm1(@as(f64, 0x1p+0)));
    // try std.testing.expectEqual(0x1.1df3b68cfb9efp+0, expm1(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x6.63992e35376b8p+0, expm1(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(0x1.315e5bf6fb106p+4, expm1(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(0x3.599205c4e74bp+4, expm1(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(0x9.369c4cb819c78p+4, expm1(@as(f64, 0x5p+0)));
    try std.testing.expectEqual(0x5.609773e54158p+12, expm1(@as(f64, 0xap+0)));
    try std.testing.expectEqual(0x3.1e1985f5a550ep+20, expm1(@as(f64, 0xfp+0)));
    try std.testing.expectEqual(0x1.ceb088a68e804p+28, expm1(@as(f64, 0x1.4p+4)));
    try std.testing.expectEqual(0x1.0c3d3920862c9p+36, expm1(@as(f64, 0x1.9p+4)));
    try std.testing.expectEqual(0x9.b823857613768p+40, expm1(@as(f64, 0x1.ep+4)));
    try std.testing.expectEqual(0x5.a278886f2355cp+48, expm1(@as(f64, 0x2.3p+4)));
    // try std.testing.expectEqual(0x3.4441a72f2e5d6p+56, expm1(@as(f64, 0x2.8p+4)));
    try std.testing.expectEqual(0x1.19103e4080b45p+72, expm1(@as(f64, 0x3.2p+4)));
    try std.testing.expectEqual(0x5.e76f27714f198p+84, expm1(@as(f64, 0x3.cp+4)));
    try std.testing.expectEqual(0x1.fbfd219c43b04p+100, expm1(@as(f64, 0x4.6p+4)));
    try std.testing.expectEqual(0xa.abbcdcc279f58p+112, expm1(@as(f64, 0x5p+4)));
    try std.testing.expectEqual(0x3.96211ff7d82c8p+128, expm1(@as(f64, 0x5.ap+4)));
    try std.testing.expectEqual(0x1.3494a9b171bf5p+144, expm1(@as(f64, 0x6.4p+4)));
    try std.testing.expectEqual(0x9.552183749161p+180, expm1(@as(f64, 0x7.fp+4)));
    // try std.testing.expectEqual(0x2.8b74553efc872p+720, expm1(@as(f64, 0x1.f4p+8)));
    try std.testing.expectEqual(std.math.inf(f64), expm1(@as(f64, 0x2.c5c4p+12)));
    try std.testing.expectEqual(-0xf.ffd0650c9537p-4, expm1(@as(f64, -0xap+0)));
    try std.testing.expectEqual(-0xf.ffffe1caa445p-4, expm1(@as(f64, -0x1p+4)));
    try std.testing.expectEqual(-0xf.fffff4e30e748p-4, expm1(@as(f64, -0x1.1p+4)));
    try std.testing.expectEqual(-0xf.fffffbe9675dp-4, expm1(@as(f64, -0x1.2p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffp-4, expm1(@as(f64, -0x2.4p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffff8p-4, expm1(@as(f64, -0x2.5p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x2.6p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x2.cp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x2.dp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x2.ep+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x4.9p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x4.ap+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x4.bp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x4.ep+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x4.fp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x5p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x6.4p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x3.e8p+8)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x2.71p+12)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f64), expm1(@as(f64, 0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f64), expm1(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), expm1(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4.8b5e3c3e81868p-4, expm1(@as(f64, 0x4p-4)));
    try std.testing.expectEqual(-0x3.8a0830a9befa8p-4, expm1(@as(f64, -0x4p-4)));
    try std.testing.expectEqual(0x4.00800aab555dcp-12, expm1(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(-0x3.ff800aaa00088p-12, expm1(@as(f64, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000800002abp-20, expm1(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(-0xf.ffff800002aa8p-24, expm1(@as(f64, -0x1p-20)));
    try std.testing.expectEqual(0x8.0000002p-32, expm1(@as(f64, 0x8p-32)));
    try std.testing.expectEqual(-0x7.ffffffep-32, expm1(@as(f64, -0x8p-32)));
    try std.testing.expectEqual(0x1.000000008p-32, expm1(@as(f64, 0x1p-32)));
    try std.testing.expectEqual(-0xf.fffffff8p-36, expm1(@as(f64, -0x1p-32)));
    try std.testing.expectEqual(0x4.0000000000008p-52, expm1(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(-0x3.ffffffffffff8p-52, expm1(@as(f64, -0x4p-52)));
    try std.testing.expectEqual(0x1p-64, expm1(@as(f64, 0x1p-64)));
    try std.testing.expectEqual(-0x1p-64, expm1(@as(f64, -0x1p-64)));
    try std.testing.expectEqual(0x1p-100, expm1(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, expm1(@as(f64, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, expm1(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-600, expm1(@as(f64, -0x1p-600)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, expm1(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0xe.4152bp-60, expm1(@as(f64, 0xe.4152bp-60)));
    try std.testing.expectEqual(0xe.4152ap-60, expm1(@as(f64, 0xe.4152ap-60)));
    try std.testing.expectEqual(0xe.4152ac57cd1fp-60, expm1(@as(f64, 0xe.4152ac57cd1fp-60)));
    try std.testing.expectEqual(0xe.4152ac57cd1e8p-60, expm1(@as(f64, 0xe.4152ac57cd1e8p-60)));
    try std.testing.expectEqual(0x7.ddee38beb9054p-4, expm1(@as(f64, 0x6.660248p-4)));
    try std.testing.expectEqual(0x7.ddee2ccfc1ecp-4, expm1(@as(f64, 0x6.66024p-4)));
    // try std.testing.expectEqual(0x7.ddee37ace0524p-4, expm1(@as(f64, 0x6.660247486aed8p-4)));
    try std.testing.expectEqual(0x7.8304264e39d2cp-4, expm1(@as(f64, 0x6.289a78p-4)));
    try std.testing.expectEqual(0x7.6f804c2bba678p-4, expm1(@as(f64, 0x6.1b4d38p-4)));
    try std.testing.expectEqual(0x7.6f804073fa444p-4, expm1(@as(f64, 0x6.1b4d3p-4)));
    // try std.testing.expectEqual(0x7.6f8042a9af784p-4, expm1(@as(f64, 0x6.1b4d318238d4cp-4)));
    try std.testing.expectEqual(0x7.6f8042a9af78p-4, expm1(@as(f64, 0x6.1b4d318238d48p-4)));
    // try std.testing.expectEqual(0x7.412de0a90d3dcp-4, expm1(@as(f64, 0x5.fb8dc8p-4)));
    // try std.testing.expectEqual(0x7.412dd50876504p-4, expm1(@as(f64, 0x5.fb8dcp-4)));
    try std.testing.expectEqual(0x7.412dde3318f34p-4, expm1(@as(f64, 0x5.fb8dc64e91a74p-4)));
    try std.testing.expectEqual(0x3.d9dcf29d1df02p-4, expm1(@as(f64, 0x3.735f4cp-4)));
    try std.testing.expectEqual(0x3.d9dceda6a6b42p-4, expm1(@as(f64, 0x3.735f48p-4)));
    try std.testing.expectEqual(0x3.d9dcef7e7e398p-4, expm1(@as(f64, 0x3.735f497c4e676p-4)));
    try std.testing.expectEqual(0x3.d9dcef7e7e396p-4, expm1(@as(f64, 0x3.735f497c4e674p-4)));
    try std.testing.expectEqual(-0xf.fe62c59d9de8p-4, expm1(@as(f64, -0x7.d6c508p+0)));
    try std.testing.expectEqual(-0xf.fe62c5aa87ba8p-4, expm1(@as(f64, -0x7.d6c51p+0)));
    // try std.testing.expectEqual(-0xf.fe62c5a2e7928p-4, expm1(@as(f64, -0x7.d6c50b469d404p+0)));
    try std.testing.expectEqual(0x1.4aaa8e05bcf71p+104, expm1(@as(f64, 0x4.857de8p+4)));
    try std.testing.expectEqual(0x7.19267f117e21p-4, expm1(@as(f64, 0x5.dfeb68p-4)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x6.a09eee9f80d9cp-4, expm1(@as(f64, 0x5.8b912p-4)));
    // try std.testing.expectEqual(0x6.a09ee34f31654p-4, expm1(@as(f64, 0x5.8b9118p-4)));
    // try std.testing.expectEqual(0x6.a09eeccd72f8cp-4, expm1(@as(f64, 0x5.8b911eb673348p-4)));
    try std.testing.expectEqual(0x6.a09eeccd72f84p-4, expm1(@as(f64, 0x5.8b911eb673344p-4)));
    try std.testing.expectEqual(0x6.c23b7ba6e78d8p-4, expm1(@as(f64, 0x5.a343ep-4)));
    // try std.testing.expectEqual(0x6.c23b7045c9d2cp-4, expm1(@as(f64, 0x5.a343d8p-4)));
    try std.testing.expectEqual(0x6.c23b7a4dd595cp-4, expm1(@as(f64, 0x5.a343df0d6800cp-4)));
    try std.testing.expectEqual(0x6.c23b7a4dd5958p-4, expm1(@as(f64, 0x5.a343df0d68008p-4)));
    try std.testing.expectEqual(0x4p-128, expm1(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, expm1(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, expm1(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, expm1(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, expm1(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, expm1(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, expm1(@as(f64, -0x4p-1076)));

    try std.testing.expectEqual(0x0p+0, expm1(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1.b7e151628aed2a6ap+0, expm1(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.1df3b68cfb9ef7aap+0, expm1(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x6.63992e35376b731p+0, expm1(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(0x1.315e5bf6fb105f2ep+4, expm1(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(0x3.599205c4e74b0cfp+4, expm1(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(0x9.369c4cb819c78fbp+4, expm1(@as(f80, 0x5p+0)));
    try std.testing.expectEqual(0x5.609773e54157e7cp+12, expm1(@as(f80, 0xap+0)));
    try std.testing.expectEqual(0x3.1e1985f5a550dde4p+20, expm1(@as(f80, 0xfp+0)));
    try std.testing.expectEqual(0x1.ceb088a68e804022p+28, expm1(@as(f80, 0x1.4p+4)));
    try std.testing.expectEqual(0x1.0c3d3920862c88aap+36, expm1(@as(f80, 0x1.9p+4)));
    try std.testing.expectEqual(0x9.b823857613764f4p+40, expm1(@as(f80, 0x1.ep+4)));
    try std.testing.expectEqual(0x5.a278886f2355ba68p+48, expm1(@as(f80, 0x2.3p+4)));
    try std.testing.expectEqual(0x3.4441a72f2e5d5068p+56, expm1(@as(f80, 0x2.8p+4)));
    try std.testing.expectEqual(0x1.19103e4080b45664p+72, expm1(@as(f80, 0x3.2p+4)));
    try std.testing.expectEqual(0x5.e76f27714f19925p+84, expm1(@as(f80, 0x3.cp+4)));
    try std.testing.expectEqual(0x1.fbfd219c43b0473p+100, expm1(@as(f80, 0x4.6p+4)));
    try std.testing.expectEqual(0xa.abbcdcc279f59e4p+112, expm1(@as(f80, 0x5p+4)));
    try std.testing.expectEqual(0x3.96211ff7d82c793p+128, expm1(@as(f80, 0x5.ap+4)));
    try std.testing.expectEqual(0x1.3494a9b171bf4accp+144, expm1(@as(f80, 0x6.4p+4)));
    try std.testing.expectEqual(0x9.552183749160e8bp+180, expm1(@as(f80, 0x7.fp+4)));
    try std.testing.expectEqual(0x2.8b74553efc87129p+720, expm1(@as(f80, 0x1.f4p+8)));
    try std.testing.expectEqual(0xc.2c2b72bac3ba40dp+16380, expm1(@as(f80, 0x2.c5c4p+12)));
    try std.testing.expectEqual(-0xf.ffd0650c953706dp-4, expm1(@as(f80, -0xap+0)));
    try std.testing.expectEqual(-0xf.ffffe1caa445118p-4, expm1(@as(f80, -0x1p+4)));
    try std.testing.expectEqual(-0xf.fffff4e30e7452dp-4, expm1(@as(f80, -0x1.1p+4)));
    try std.testing.expectEqual(-0xf.fffffbe9675ce5ap-4, expm1(@as(f80, -0x1.2p+4)));
    try std.testing.expectEqual(-0xf.fffffffffffef49p-4, expm1(@as(f80, -0x2.4p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffff9dap-4, expm1(@as(f80, -0x2.5p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffdbdp-4, expm1(@as(f80, -0x2.6p+4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffp-4, expm1(@as(f80, -0x2.cp+4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffp-4, expm1(@as(f80, -0x2.dp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x2.ep+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x4.9p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x4.ap+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x4.bp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x4.ep+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x4.fp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x5p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x6.4p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x3.e8p+8)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x2.71p+12)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f80), expm1(@as(f80, 0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f80), expm1(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f80), expm1(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), expm1(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4.8b5e3c3e81866768p-4, expm1(@as(f80, 0x4p-4)));
    try std.testing.expectEqual(-0x3.8a0830a9befa8bccp-4, expm1(@as(f80, -0x4p-4)));
    try std.testing.expectEqual(0x4.00800aab555dde38p-12, expm1(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(-0x3.ff800aaa0008882cp-12, expm1(@as(f80, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000800002aaaacp-20, expm1(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(-0xf.ffff800002aaaaap-24, expm1(@as(f80, -0x1p-20)));
    try std.testing.expectEqual(0x8.000000200000005p-32, expm1(@as(f80, 0x8p-32)));
    try std.testing.expectEqual(-0x7.ffffffe000000058p-32, expm1(@as(f80, -0x8p-32)));
    try std.testing.expectEqual(0x1.000000008p-32, expm1(@as(f80, 0x1p-32)));
    try std.testing.expectEqual(-0xf.fffffff8p-36, expm1(@as(f80, -0x1p-32)));
    try std.testing.expectEqual(0x4.0000000000008p-52, expm1(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(-0x3.ffffffffffff8p-52, expm1(@as(f80, -0x4p-52)));
    try std.testing.expectEqual(0x1p-64, expm1(@as(f80, 0x1p-64)));
    try std.testing.expectEqual(-0x1p-64, expm1(@as(f80, -0x1p-64)));
    try std.testing.expectEqual(0x1p-100, expm1(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, expm1(@as(f80, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, expm1(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-600, expm1(@as(f80, -0x1p-600)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, expm1(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, expm1(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1p-10000, expm1(@as(f80, -0x1p-10000)));
    try std.testing.expectEqual(0xe.4152b0000000066p-60, expm1(@as(f80, 0xe.4152bp-60)));
    try std.testing.expectEqual(0xe.4152a0000000066p-60, expm1(@as(f80, 0xe.4152ap-60)));
    try std.testing.expectEqual(0xe.4152ac57cd1f066p-60, expm1(@as(f80, 0xe.4152ac57cd1fp-60)));
    try std.testing.expectEqual(0xe.4152ac57cd1e866p-60, expm1(@as(f80, 0xe.4152ac57cd1e8p-60)));
    try std.testing.expectEqual(0xe.4152ac57cd1eaep-60, expm1(@as(f80, 0xe.4152ac57cd1ea7ap-60)));
    try std.testing.expectEqual(0x7.ddee38beb90553dp-4, expm1(@as(f80, 0x6.660248p-4)));
    try std.testing.expectEqual(0x7.ddee2ccfc1ebf03p-4, expm1(@as(f80, 0x6.66024p-4)));
    try std.testing.expectEqual(0x7.ddee37ace0525dep-4, expm1(@as(f80, 0x6.660247486aed8p-4)));
    try std.testing.expectEqual(0x7.8304264e39d2dp-4, expm1(@as(f80, 0x6.289a78p-4)));
    try std.testing.expectEqual(0x7.6f804c2bba6774a8p-4, expm1(@as(f80, 0x6.1b4d38p-4)));
    try std.testing.expectEqual(0x7.6f804073fa444cb8p-4, expm1(@as(f80, 0x6.1b4d3p-4)));
    try std.testing.expectEqual(0x7.6f8042a9af7859dp-4, expm1(@as(f80, 0x6.1b4d318238d4cp-4)));
    try std.testing.expectEqual(0x7.6f8042a9af77fc1p-4, expm1(@as(f80, 0x6.1b4d318238d48p-4)));
    try std.testing.expectEqual(0x7.6f8042a9af782ed8p-4, expm1(@as(f80, 0x6.1b4d318238d4a2a8p-4)));
    try std.testing.expectEqual(0x7.412de0a90d3dcc38p-4, expm1(@as(f80, 0x5.fb8dc8p-4)));
    try std.testing.expectEqual(0x7.412dd50876505fd8p-4, expm1(@as(f80, 0x5.fb8dcp-4)));
    try std.testing.expectEqual(0x7.412dde3318f344cp-4, expm1(@as(f80, 0x5.fb8dc64e91a74p-4)));
    try std.testing.expectEqual(0x3.d9dcf29d1df01bdp-4, expm1(@as(f80, 0x3.735f4cp-4)));
    try std.testing.expectEqual(0x3.d9dceda6a6b41354p-4, expm1(@as(f80, 0x3.735f48p-4)));
    try std.testing.expectEqual(0x3.d9dcef7e7e398c14p-4, expm1(@as(f80, 0x3.735f497c4e676p-4)));
    try std.testing.expectEqual(0x3.d9dcef7e7e39646p-4, expm1(@as(f80, 0x3.735f497c4e674p-4)));
    try std.testing.expectEqual(0x3.d9dcef7e7e397c64p-4, expm1(@as(f80, 0x3.735f497c4e67535cp-4)));
    try std.testing.expectEqual(-0xf.fe62c59d9de7d61p-4, expm1(@as(f80, -0x7.d6c508p+0)));
    try std.testing.expectEqual(-0xf.fe62c5aa87bab58p-4, expm1(@as(f80, -0x7.d6c51p+0)));
    try std.testing.expectEqual(-0xf.fe62c5a2e792cp-4, expm1(@as(f80, -0x7.d6c50b469d404p+0)));
    try std.testing.expectEqual(0x1.4aaa8e05bcf71098p+104, expm1(@as(f80, 0x4.857de8p+4)));
    try std.testing.expectEqual(0x7.19267f117e20e3b8p-4, expm1(@as(f80, 0x5.dfeb68p-4)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x4.0000000000000028p-16384, expm1(@as(f80, 0x4.0000000000000028p-16384)));
    try std.testing.expectEqual(0x6.a09eee9f80d9d6fp-4, expm1(@as(f80, 0x5.8b912p-4)));
    try std.testing.expectEqual(0x6.a09ee34f31655b48p-4, expm1(@as(f80, 0x5.8b9118p-4)));
    try std.testing.expectEqual(0x6.a09eeccd72f8a578p-4, expm1(@as(f80, 0x5.8b911eb673348p-4)));
    try std.testing.expectEqual(0x6.a09eeccd72f84afp-4, expm1(@as(f80, 0x5.8b911eb673344p-4)));
    try std.testing.expectEqual(0x6.a09eeccd72f88608p-4, expm1(@as(f80, 0x5.8b911eb6733469c8p-4)));
    try std.testing.expectEqual(0x6.c23b7ba6e78d9eb8p-4, expm1(@as(f80, 0x5.a343ep-4)));
    try std.testing.expectEqual(0x6.c23b7045c9d2a388p-4, expm1(@as(f80, 0x5.a343d8p-4)));
    try std.testing.expectEqual(0x6.c23b7a4dd595d6b8p-4, expm1(@as(f80, 0x5.a343df0d6800cp-4)));
    try std.testing.expectEqual(0x6.c23b7a4dd5957bbp-4, expm1(@as(f80, 0x5.a343df0d68008p-4)));
    try std.testing.expectEqual(0x6.c23b7a4dd595a028p-4, expm1(@as(f80, 0x5.a343df0d680099a8p-4)));
    try std.testing.expectEqual(0x6.c23b7a4dd595a02p-4, expm1(@as(f80, 0x5.a343df0d680099ap-4)));
    try std.testing.expectEqual(0x4p-128, expm1(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, expm1(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, expm1(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, expm1(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, expm1(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, expm1(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, expm1(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, expm1(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, expm1(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, expm1(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, expm1(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, expm1(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, expm1(@as(f80, -0x8p-16448)));

    try std.testing.expectEqual(0x0p+0, expm1(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1.b7e151628aed2a6abf7158809cf5p+0, expm1(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.1df3b68cfb9ef7a986addc7dcee2p+0, expm1(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x6.63992e35376b730ce8ee881ada2cp+0, expm1(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(0x1.315e5bf6fb105f2d4bdfc53744c4p+4, expm1(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(0x3.599205c4e74b0cf1ada77fb727b8p+4, expm1(@as(f128, 0x4p+0)));
    try std.testing.expectEqual(0x9.369c4cb819c78fb37d56c91ad5fp+4, expm1(@as(f128, 0x5p+0)));
    try std.testing.expectEqual(0x5.609773e54157e7c1faa3015b44d4p+12, expm1(@as(f128, 0xap+0)));
    // try std.testing.expectEqual(0x3.1e1985f5a550dde2e5fe372cd4b2p+20, expm1(@as(f128, 0xfp+0)));
    try std.testing.expectEqual(0x1.ceb088a68e80402189797f9599ccp+28, expm1(@as(f128, 0x1.4p+4)));
    try std.testing.expectEqual(0x1.0c3d3920862c88aafb2ae72d6857p+36, expm1(@as(f128, 0x1.9p+4)));
    try std.testing.expectEqual(0x9.b823857613764f43e201f73a543p+40, expm1(@as(f128, 0x1.ep+4)));
    try std.testing.expectEqual(0x5.a278886f2355ba66b452ea7226f4p+48, expm1(@as(f128, 0x2.3p+4)));
    try std.testing.expectEqual(0x3.4441a72f2e5d50686c20e8b55b3cp+56, expm1(@as(f128, 0x2.8p+4)));
    try std.testing.expectEqual(0x1.19103e4080b45664d6740cf8c5d9p+72, expm1(@as(f128, 0x3.2p+4)));
    try std.testing.expectEqual(0x5.e76f27714f19924caf2a55081894p+84, expm1(@as(f128, 0x3.cp+4)));
    try std.testing.expectEqual(0x1.fbfd219c43b04730797e2bfeb1cfp+100, expm1(@as(f128, 0x4.6p+4)));
    // try std.testing.expectEqual(0xa.abbcdcc279f59e45281da547124p+112, expm1(@as(f128, 0x5p+4)));
    try std.testing.expectEqual(0x3.96211ff7d82c792f823b2ba3a166p+128, expm1(@as(f128, 0x5.ap+4)));
    try std.testing.expectEqual(0x1.3494a9b171bf4acc225093322428p+144, expm1(@as(f128, 0x6.4p+4)));
    try std.testing.expectEqual(0x9.552183749160e8b702888dad951p+180, expm1(@as(f128, 0x7.fp+4)));
    try std.testing.expectEqual(0x2.8b74553efc87128fd5d1b2c1ea3ap+720, expm1(@as(f128, 0x1.f4p+8)));
    try std.testing.expectEqual(0xc.2c2b72bac3ba40c9d77771f196dp+16380, expm1(@as(f128, 0x2.c5c4p+12)));
    try std.testing.expectEqual(-0xf.ffd0650c953706cac749b7155edp-4, expm1(@as(f128, -0xap+0)));
    try std.testing.expectEqual(-0xf.ffffe1caa445117a35259a08c0dp-4, expm1(@as(f128, -0x1p+4)));
    try std.testing.expectEqual(-0xf.fffff4e30e7452cbb1a1331e22bp-4, expm1(@as(f128, -0x1.1p+4)));
    try std.testing.expectEqual(-0xf.fffffbe9675ce59817cddee3aa18p-4, expm1(@as(f128, -0x1.2p+4)));
    try std.testing.expectEqual(-0xf.fffffffffffef493c50221f9c7f8p-4, expm1(@as(f128, -0x2.4p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffff9d9ee380d67eac08p-4, expm1(@as(f128, -0x2.5p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffdbceea52a399f9e8p-4, expm1(@as(f128, -0x2.6p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffffe908be21e8b718p-4, expm1(@as(f128, -0x2.cp+4)));
    try std.testing.expectEqual(-0xf.fffffffffffffff78d246170056p-4, expm1(@as(f128, -0x2.dp+4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffce4543c89ec5p-4, expm1(@as(f128, -0x2.ep+4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffff99p-4, expm1(@as(f128, -0x4.9p+4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffffdap-4, expm1(@as(f128, -0x4.ap+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffffffffffffffff2p-4, expm1(@as(f128, -0x4.bp+4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffffff8p-4, expm1(@as(f128, -0x4.ep+4)));
    // try std.testing.expectEqual(-0xf.fffffffffffffffffffffffffff8p-4, expm1(@as(f128, -0x4.fp+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0x5p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0x6.4p+4)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0x3.e8p+8)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0x2.71p+12)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f128), expm1(@as(f128, 0x1.86ap+16)));
    try std.testing.expectEqual(std.math.inf(f128), expm1(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f128), expm1(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), expm1(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f128), expm1(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(std.math.inf(f128), expm1(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1p+0, expm1(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4.8b5e3c3e81866767bc3b69baabe4p-4, expm1(@as(f128, 0x4p-4)));
    try std.testing.expectEqual(-0x3.8a0830a9befa8bcbea343629c97p-4, expm1(@as(f128, -0x4p-4)));
    try std.testing.expectEqual(0x4.00800aab555dde38e6ce86e9277cp-12, expm1(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(-0x3.ff800aaa0008882d861847853132p-12, expm1(@as(f128, -0x4p-12)));
    try std.testing.expectEqual(0x1.00000800002aaaab55555777777dp-20, expm1(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(-0xf.ffff800002aaaaa00000222221c8p-24, expm1(@as(f128, -0x1p-20)));
    try std.testing.expectEqual(0x8.0000002000000055555556p-32, expm1(@as(f128, 0x8p-32)));
    try std.testing.expectEqual(-0x7.ffffffe000000055555554aaaaacp-32, expm1(@as(f128, -0x8p-32)));
    try std.testing.expectEqual(0x1.00000000800000002aaaaaaab555p-32, expm1(@as(f128, 0x1p-32)));
    try std.testing.expectEqual(-0xf.fffffff800000002aaaaaaaap-36, expm1(@as(f128, -0x1p-32)));
    try std.testing.expectEqual(0x4.0000000000008000000000000aacp-52, expm1(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(-0x3.ffffffffffff8000000000000aaap-52, expm1(@as(f128, -0x4p-52)));
    try std.testing.expectEqual(0x1.00000000000000008p-64, expm1(@as(f128, 0x1p-64)));
    try std.testing.expectEqual(-0xf.fffffffffffffff8p-68, expm1(@as(f128, -0x1p-64)));
    try std.testing.expectEqual(0x1.00000000000000000000000008p-100, expm1(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(-0xf.ffffffffffffffffffffffff8p-104, expm1(@as(f128, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, expm1(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-600, expm1(@as(f128, -0x1p-600)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, expm1(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(-0x0p+0, expm1(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, expm1(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1p-10000, expm1(@as(f128, -0x1p-10000)));
    try std.testing.expectEqual(0xe.4152b00000000659adb2c0a9c8p-60, expm1(@as(f128, 0xe.4152bp-60)));
    try std.testing.expectEqual(0xe.4152a00000000659ada47f572p-60, expm1(@as(f128, 0xe.4152ap-60)));
    try std.testing.expectEqual(0xe.4152ac57cd1f0659adaf7e8f0e98p-60, expm1(@as(f128, 0xe.4152ac57cd1fp-60)));
    try std.testing.expectEqual(0xe.4152ac57cd1e8659adaf7e8f0e28p-60, expm1(@as(f128, 0xe.4152ac57cd1e8p-60)));
    try std.testing.expectEqual(0xe.4152ac57cd1eadf9adaf7e8f0e5p-60, expm1(@as(f128, 0xe.4152ac57cd1ea7ap-60)));
    try std.testing.expectEqual(0x7.ddee38beb90553d11ec1beb27a08p-4, expm1(@as(f128, 0x6.660248p-4)));
    try std.testing.expectEqual(0x7.ddee2ccfc1ebf03262b062f7fa4p-4, expm1(@as(f128, 0x6.66024p-4)));
    try std.testing.expectEqual(0x7.ddee37ace0525de2e3d415373edp-4, expm1(@as(f128, 0x6.660247486aed8p-4)));
    try std.testing.expectEqual(0x7.8304264e39d2cffd272b76863fa8p-4, expm1(@as(f128, 0x6.289a78p-4)));
    try std.testing.expectEqual(0x7.6f804c2bba6774a43c5c34a96fcp-4, expm1(@as(f128, 0x6.1b4d38p-4)));
    try std.testing.expectEqual(0x7.6f804073fa444cb711aa5c815448p-4, expm1(@as(f128, 0x6.1b4d3p-4)));
    // try std.testing.expectEqual(0x7.6f8042a9af7859cfbcbdca8df28cp-4, expm1(@as(f128, 0x6.1b4d318238d4cp-4)));
    try std.testing.expectEqual(0x7.6f8042a9af77fc11bbb323d011ep-4, expm1(@as(f128, 0x6.1b4d318238d48p-4)));
    // try std.testing.expectEqual(0x7.6f8042a9af782ed4bf03885aa37cp-4, expm1(@as(f128, 0x6.1b4d318238d4a2a8p-4)));
    try std.testing.expectEqual(0x7.412de0a90d3dcc39dc4e01aef3b8p-4, expm1(@as(f128, 0x5.fb8dc8p-4)));
    try std.testing.expectEqual(0x7.412dd50876505fd8f901001e9238p-4, expm1(@as(f128, 0x5.fb8dcp-4)));
    try std.testing.expectEqual(0x7.412dde3318f344bc6f6cd00f8ed8p-4, expm1(@as(f128, 0x5.fb8dc64e91a74p-4)));
    try std.testing.expectEqual(0x3.d9dcf29d1df01bce7765d392e54ap-4, expm1(@as(f128, 0x3.735f4cp-4)));
    try std.testing.expectEqual(0x3.d9dceda6a6b41355e2e68cfbebcap-4, expm1(@as(f128, 0x3.735f48p-4)));
    // try std.testing.expectEqual(0x3.d9dcef7e7e398c13715cfcb81c9ap-4, expm1(@as(f128, 0x3.735f497c4e676p-4)));
    try std.testing.expectEqual(0x3.d9dcef7e7e39645fb77dffbba9aap-4, expm1(@as(f128, 0x3.735f497c4e674p-4)));
    try std.testing.expectEqual(0x3.d9dcef7e7e397c649290c708e3aap-4, expm1(@as(f128, 0x3.735f497c4e67535cp-4)));
    try std.testing.expectEqual(-0xf.fe62c59d9de7d6168bf8c31a716p-4, expm1(@as(f128, -0x7.d6c508p+0)));
    try std.testing.expectEqual(-0xf.fe62c5aa87bab580018589d526p-4, expm1(@as(f128, -0x7.d6c51p+0)));
    try std.testing.expectEqual(-0xf.fe62c5a2e792bffeb1e98cc705dp-4, expm1(@as(f128, -0x7.d6c50b469d404p+0)));
    try std.testing.expectEqual(0x1.4aaa8e05bcf71097ff88abf1c0adp+104, expm1(@as(f128, 0x4.857de8p+4)));
    // try std.testing.expectEqual(0x7.19267f117e20e3b9a8b8bdf57f74p-4, expm1(@as(f128, 0x5.dfeb68p-4)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, expm1(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x4.0000000000000028p-16384, expm1(@as(f128, 0x4.0000000000000028p-16384)));
    try std.testing.expectEqual(0x6.a09eee9f80d9d6f2256cdf5a8698p-4, expm1(@as(f128, 0x5.8b912p-4)));
    try std.testing.expectEqual(0x6.a09ee34f31655b4595dca868cc2p-4, expm1(@as(f128, 0x5.8b9118p-4)));
    try std.testing.expectEqual(0x6.a09eeccd72f8a5749e8d79476188p-4, expm1(@as(f128, 0x5.8b911eb673348p-4)));
    // try std.testing.expectEqual(0x6.a09eeccd72f84af222da437b7fa4p-4, expm1(@as(f128, 0x5.8b911eb673344p-4)));
    // try std.testing.expectEqual(0x6.a09eeccd72f8860891dba1fa3954p-4, expm1(@as(f128, 0x5.8b911eb6733469c8p-4)));
    // try std.testing.expectEqual(0x6.c23b7ba6e78d9eb5a8c9564c5134p-4, expm1(@as(f128, 0x5.a343ep-4)));
    // try std.testing.expectEqual(0x6.c23b7045c9d2a38950f5772c7164p-4, expm1(@as(f128, 0x5.a343d8p-4)));
    // try std.testing.expectEqual(0x6.c23b7a4dd595d6b551ae9c270edcp-4, expm1(@as(f128, 0x5.a343df0d6800cp-4)));
    try std.testing.expectEqual(0x6.c23b7a4dd5957bac63c564d0b838p-4, expm1(@as(f128, 0x5.a343df0d68008p-4)));
    try std.testing.expectEqual(0x6.c23b7a4dd595a02ad825029fb488p-4, expm1(@as(f128, 0x5.a343df0d680099a8p-4)));
    // try std.testing.expectEqual(0x6.c23b7a4dd595a01f77074578c9bcp-4, expm1(@as(f128, 0x5.a343df0d680099ap-4)));
    // try std.testing.expectEqual(0x6.c23b7a4dd595a02a51f35d6c1588p-4, expm1(@as(f128, 0x5.a343df0d680099a7a1a873a751a8p-4)));
    // try std.testing.expectEqual(0x6.c23b7a4dd595a02a51f35d6c1604p-4, expm1(@as(f128, 0x5.a343df0d680099a7a1a873a752p-4)));
    // try std.testing.expectEqual(0x6.c23b7a4dd595a02a51f35d6c132cp-4, expm1(@as(f128, 0x5.a343df0d680099a7a1a873a75p-4)));
    try std.testing.expectEqual(0x4p-128, expm1(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, expm1(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, expm1(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, expm1(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, expm1(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, expm1(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, expm1(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, expm1(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, expm1(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, expm1(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, expm1(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, expm1(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, expm1(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, expm1(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, expm1(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-152, expm1(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, expm1(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, expm1(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, expm1(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, expm1(@as(f128, -0x4p-16496)));
}
