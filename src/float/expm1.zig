const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
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

            return @mulAdd(f32, float.abs(x), 0x1p-25, x);
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
        if (float.abs(x) < std.math.floatMin(f64)) {
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
        return float.exp(x);
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
    if (float.abs(x) < 0x1p-113) {
        if (float.abs(x) < std.math.floatMin(f128)) {
            const vx: f128 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }

        return x;
    }

    // Express x = ln 2 (k + remainder), remainder not exceeding 1/2.
    var xx: f128 = C1 + C2; // ln 2.
    var px: f128 = float.floor(0.5 + x / xx);
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

    px = float.ldexp(1, k);
    y = px * qx + (px - 1.0);
    return y;
}
