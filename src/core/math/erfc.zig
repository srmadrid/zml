const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const erf_data = @import("erf_data.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn erfc(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return erfc(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, erfc32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_erfcf.c
                    return erfc32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_erf.c
                    return erfc64(x);
                },
                f80 => return cast(f80, erfc128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_erfl.c
                    return erfc128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn erfc32(xf: f32) f32 {
    const axf: f32 = math.abs(xf);
    const axd: f64 = cast(f64, axf, .{});
    const x2: f64 = axd * axd;
    const t: u32 = @bitCast(xf);
    const at: u32 = t & (~@as(u32, 0) >> 1);
    const sgn: u32 = t >> 31;
    const i: u32 = if (at > 0x40051000) 1 else 0;
    // for x < -0x1.ea8f94p+1, erfc(x) rounds to 2 (to nearest)
    if (t > 0xc07547ca) { // xf < -0x1.ea8f94p+1
        @branchHint(.unlikely);
        if (t >= 0xff800000) { // -Inf or NaN
            @branchHint(.unlikely);
            if (t == 0xff800000)
                return 2; // -Inf

            return xf + xf; // NaN
        }
        return 2 - 0x1p-25; // rounds to 2 or nextbelow(2)
    }

    // at is the absolute value of xf
    // for x >= 0x1.41bbf8p+3, erfc(x) < 2^-150, thus rounds to 0 or to 2^-149
    // depending on the rounding mode
    if (at >= 0x4120ddfc) { // |xf| >= 0x1.41bbf8p+3
        @branchHint(.unlikely);
        if (at >= 0x7f800000) { // +Inf or NaN
            @branchHint(.unlikely);
            if (at == 0x7f800000)
                return 0; // +Inf

            return xf + xf; // NaN
        }

        // 0x1p-149f * 0.25f rounds to 0 or 2^-149 depending on rounding
        return 0x1p-149 * 0.25;
    }

    if (at <= 0x3db80000) { // |x| <= 0x1.7p-4
        @branchHint(.unlikely);
        if (t == 0xb76c9f62) {
            @branchHint(.unlikely);
            return 0x1.00010ap+0 + 0x1p-25; // exceptional case
        }

        // for |x| <= 0x1.c5bf88p-26. erfc(x) rounds to 1 (to nearest)
        if (at <= 0x32e2dfc4) { // |x| <= 0x1.c5bf88p-26
            @branchHint(.unlikely);
            if (at == 0) {
                @branchHint(.unlikely);
                return 1;
            }

            const d: [2]f32 = .{ -0x1p-26, 0x1p-25 };
            return 1 + d[sgn];
        }

        // around 0, erfc(x) behaves as 1 - (odd polynomial)
        const c: [5]f64 = .{ 0x1.20dd750429b6dp+0, -0x1.812746b03610bp-2, 0x1.ce2f218831d2fp-4, -0x1.b82c609607dcbp-6, 0x1.553af09b8008ep-8 };
        const f0: f64 = xf * (c[0] + x2 * (c[1] + x2 * (c[2] + x2 * (c[3] + x2 * (c[4])))));

        return cast(f32, 1 - f0, .{});
    }

    // now -0x1.ea8f94p+1 <= x <= 0x1.41bbf8p+3, with |x| > 0x1.7p-4
    const iln2: f64 = 0x1.71547652b82fep+0;
    const ln2h: f64 = 0x1.62e42fefap-8;
    const ln2l: f64 = 0x1.cf79abd6f5dc8p-47;
    const jt: u64 = @bitCast(@mulAdd(f64, x2, iln2, -(1024 + 0x1p-8)));
    var j: i64 = undefined;
    {
        @setRuntimeSafety(false);
        j = @as(i64, @intCast(jt << 12)) >> 48;
    }
    const S: f64 = @bitCast(((j >> 7) + cast(i64, 0x3ff | sgn << 11, .{})) << 52);
    const ch: [4]f64 = .{ -0x1.ffffffffff333p-2, 0x1.5555555556a14p-3, -0x1.55556666659b4p-5, 0x1.1111074cc7b22p-7 };
    const d: f64 = (x2 + ln2h * cast(f64, j, .{})) + ln2l * cast(f64, j, .{});
    const d2: f64 = d * d;
    const e0: f64 = erf_data.E_32[@intCast(j & 127)];
    const f: f64 = d + d2 * ((ch[0] + d * ch[1]) + d2 * (ch[2] + d * ch[3]));
    const ct: [2][16]f64 = .{
        .{
            0x1.c162355429b28p-1,  0x1.d99999999999ap+1,  0x1.da951cece2b85p-2,
            -0x1.70ef6cff4bcc4p+0, 0x1.3d7f7b3d617dep+1,  -0x1.9d0aa47537c51p+1,
            0x1.9754ea9a3fcb1p+1,  -0x1.27a5453fcc015p+1, 0x1.1ef2e0531aebap+0,
            -0x1.eca090f5a1c06p-3, -0x1.7a3cd173a063cp-4, 0x1.30fa68a68fdddp-4,
            0x1.55ad9a326993ap-10, -0x1.07e7b0bb39fbfp-6, 0x1.2328706c0e95p-10,
            0x1.d6aa0b7b19cfep-9,
        },
        .{
            0x1.137c8983f8516p+2,  0x1.799999999999ap+1,   0x1.05b53aa241333p-3,
            -0x1.a3f53872bf87p-3,  0x1.de4c30742c9d5p-4,   -0x1.cb24bfa591986p-5,
            0x1.666aec059ca5fp-6,  -0x1.a61250eb26b0bp-8,  0x1.2b28b7924b34dp-10,
            0x1.41b13a9d45013p-15, -0x1.6dd5e8a273613p-14, 0x1.09ce8ea5e8da5p-16,
            0x1.33923b4102981p-18, -0x1.1dfd161e3f984p-19, -0x1.c87618fcae3b3p-23,
            0x1.e8a6ffa0ba2c7p-23,
        },
    };
    const z: f64 = (axd - ct[i][0]) / (axd + ct[i][1]);
    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    const z8: f64 = z4 * z4;
    const c: [16]f64 = ct[i];
    var s: f64 = (((c[3] + z * c[4]) + z2 * (c[5] + z * c[6])) + z4 * ((c[7] + z * c[8]) + z2 * (c[9] + z * c[10]))) + z8 * (((c[11] + z * c[12]) + z2 * (c[13] + z * c[14])) + z4 * (c[15]));
    s = ct[i][2] + z * s;
    const off: [2]f64 = .{ 0, 2 };
    const r: f64 = (S * (e0 - f * e0)) * s;
    const y: f64 = off[sgn] + r;

    return cast(f32, y, .{});
}

fn erfc64(x: f64) f64 {
    var hx: i32 = undefined;
    dbl64.getHighWord(&hx, x);
    const ix: i32 = hx & 0x7fffffff;
    if (ix >= 0x7ff00000) { // erfc(nan)=nan, erfc(+-inf)=0,2
        const ret: f64 = cast(f64, ((hx >> 31) << 1) + 1, .{}) / x;

        if (ret == 0)
            return 0;

        return ret;
    }

    if (ix < 0x3feb0000) { // |x|<0.84375
        if (ix < 0x3c700000) // |x|<2**-56
            return 1 - x;

        const z: f64 = x * x;
        const r1: f64 = erf_data.pp_64[0] + z * erf_data.pp_64[1];
        const z2: f64 = z * z;
        const r2: f64 = erf_data.pp_64[2] + z * erf_data.pp_64[3];
        const z4: f64 = z2 * z2;
        const s1: f64 = 1 + z * erf_data.qq_64[1];
        const s2: f64 = erf_data.qq_64[2] + z * erf_data.qq_64[3];
        const s3: f64 = erf_data.qq_64[4] + z * erf_data.qq_64[5];
        var r: f64 = r1 + z2 * r2 + z4 * erf_data.pp_64[4];
        const s: f64 = s1 + z2 * s2 + z4 * s3;
        const y: f64 = r / s;
        if (hx < 0x3fd00000) { // x<1/4
            return 1 - (x + x * y);
        } else {
            r = x * y;
            r += (x - 0.5);
            return 0.5 - r;
        }
    }

    if (ix < 0x3ff40000) { // 0.84375 <= |x| < 1.25
        const s: f64 = math.abs(x) - 1;
        const P1: f64 = erf_data.pa_64[0] + s * erf_data.pa_64[1];
        const s2: f64 = s * s;
        const Q1: f64 = 1 + s * erf_data.qa_64[1];
        const s4: f64 = s2 * s2;
        const P2: f64 = erf_data.pa_64[2] + s * erf_data.pa_64[3];
        const s6: f64 = s4 * s2;
        const Q2: f64 = erf_data.qa_64[2] + s * erf_data.qa_64[3];
        const P3: f64 = erf_data.pa_64[4] + s * erf_data.pa_64[5];
        const Q3: f64 = erf_data.qa_64[4] + s * erf_data.qa_64[5];
        const P4: f64 = erf_data.pa_64[6];
        const Q4: f64 = erf_data.qa_64[6];
        const P: f64 = P1 + s2 * P2 + s4 * P3 + s6 * P4;
        const Q: f64 = Q1 + s2 * Q2 + s4 * Q3 + s6 * Q4;
        if (hx >= 0) {
            const z: f64 = 1 - erf_data.erx_64;
            return z - P / Q;
        } else {
            const z: f64 = erf_data.erx_64 + P / Q;
            return 1 + z;
        }
    }

    if (ix < 0x403c0000) { // |x|<28
        const xx: f64 = math.abs(x);
        const s: f64 = 1 / (xx * xx);
        var R: f64 = undefined;
        var S: f64 = undefined;
        if (ix < 0x4006db6d) { // |x| < 1/.35 ~ 2.857143
            const R1: f64 = erf_data.ra_64[0] + s * erf_data.ra_64[1];
            const s2: f64 = s * s;
            const S1: f64 = 1 + s * erf_data.sa_64[1];
            const s4: f64 = s2 * s2;
            const R2: f64 = erf_data.ra_64[2] + s * erf_data.ra_64[3];
            const s6: f64 = s4 * s2;
            const S2: f64 = erf_data.sa_64[2] + s * erf_data.sa_64[3];
            const s8: f64 = s4 * s4;
            const R3: f64 = erf_data.ra_64[4] + s * erf_data.ra_64[5];
            const S3: f64 = erf_data.sa_64[4] + s * erf_data.sa_64[5];
            const R4: f64 = erf_data.ra_64[6] + s * erf_data.ra_64[7];
            const S4: f64 = erf_data.sa_64[6] + s * erf_data.sa_64[7];
            R = R1 + s2 * R2 + s4 * R3 + s6 * R4;
            S = S1 + s2 * S2 + s4 * S3 + s6 * S4 + s8 * erf_data.sa_64[8];
        } else { // |x| >= 1/.35 ~ 2.857143
            if (hx < 0 and ix >= 0x40180000)
                return 2 - erf_data.tiny_64; // x < -6
            const R1: f64 = erf_data.rb_64[0] + s * erf_data.rb_64[1];
            const s2: f64 = s * s;
            const S1: f64 = 1 + s * erf_data.sb_64[1];
            const s4: f64 = s2 * s2;
            const R2: f64 = erf_data.rb_64[2] + s * erf_data.rb_64[3];
            const s6: f64 = s4 * s2;
            const S2: f64 = erf_data.sb_64[2] + s * erf_data.sb_64[3];
            const R3: f64 = erf_data.rb_64[4] + s * erf_data.rb_64[5];
            const S3: f64 = erf_data.sb_64[4] + s * erf_data.sb_64[5];
            const S4: f64 = erf_data.sb_64[6] + s * erf_data.sb_64[7];
            R = R1 + s2 * R2 + s4 * R3 + s6 * erf_data.rb_64[6];
            S = S1 + s2 * S2 + s4 * S3 + s6 * S4;
        }

        var z: f64 = xx;
        dbl64.setLowWord(&z, @as(u32, 0));
        const r: f64 = math.exp(-z * z - 0.5625) * math.exp((z - xx) * (z + xx) + R / S);
        if (hx > 0) {
            const ret: f64 = r / xx;
            return ret;
        } else {
            return 2 - r / xx;
        }
    } else {
        if (hx > 0) {
            return erf_data.tiny_64 * erf_data.tiny_64;
        } else return 2 - erf_data.tiny_64;
    }
}

fn erfc128(x: f128) f128 {
    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const sign: i32 = @bitCast(u.w0);
    const ix: i32 = sign & 0x7fffffff;
    u.w0 = @bitCast(ix);

    if (ix >= 0x7fff0000) { // erfc(nan)=nan
        // erfc(+-inf)=0,2
        return cast(f128, ((@as(u32, @intCast(sign)) >> 31) << 1) + 1, .{}) / x;
    }

    if (ix < 0x3ffd0000) { // |x| <1/4
        if (ix < 0x3f8d0000) // |x|<2**-114
            return 1 - x;

        return 1 - math.erf(x);
    }

    if (ix < 0x3fff4000) { // 1.25
        const xx: f128 = @bitCast(u);
        const i: i32 = cast(i32, 8.0 * xx, .{});
        var y: f128 = undefined;
        switch (i) {
            2 => {
                const z: f128 = xx - 0.25;
                y = erf_data.C13b_128 + z * erf_data.neval(z, &erf_data.RNr13_128, 8) / erf_data.deval(z, &erf_data.RDr13_128, 7);
                y += erf_data.C13a_128;
            },
            3 => {
                const z: f128 = xx - 0.375;
                y = erf_data.C14b_128 + z * erf_data.neval(z, &erf_data.RNr14_128, 8) / erf_data.deval(z, &erf_data.RDr14_128, 7);
                y += erf_data.C14a_128;
            },
            4 => {
                const z: f128 = xx - 0.5;
                y = erf_data.C15b_128 + z * erf_data.neval(z, &erf_data.RNr15_128, 8) / erf_data.deval(z, &erf_data.RDr15_128, 7);
                y += erf_data.C15a_128;
            },
            5 => {
                const z: f128 = xx - 0.625;
                y = erf_data.C16b_128 + z * erf_data.neval(z, &erf_data.RNr16_128, 8) / erf_data.deval(z, &erf_data.RDr16_128, 7);
                y += erf_data.C16a_128;
            },
            6 => {
                const z: f128 = xx - 0.75;
                y = erf_data.C17b_128 + z * erf_data.neval(z, &erf_data.RNr17_128, 8) / erf_data.deval(z, &erf_data.RDr17_128, 7);
                y += erf_data.C17a_128;
            },
            7 => {
                const z: f128 = xx - 0.875;
                y = erf_data.C18b_128 + z * erf_data.neval(z, &erf_data.RNr18_128, 8) / erf_data.deval(z, &erf_data.RDr18_128, 7);
                y += erf_data.C18a_128;
            },
            8 => {
                const z: f128 = xx - 1;
                y = erf_data.C19b_128 + z * erf_data.neval(z, &erf_data.RNr19_128, 8) / erf_data.deval(z, &erf_data.RDr19_128, 7);
                y += erf_data.C19a_128;
            },
            else => {
                const z: f128 = xx - 1.125;
                y = erf_data.C20b_128 + z * erf_data.neval(z, &erf_data.RNr20_128, 8) / erf_data.deval(z, &erf_data.RDr20_128, 7);
                y += erf_data.C20a_128;
            },
        }

        if ((@as(u32, @bitCast(sign)) & 0x80000000) != 0)
            y = 2 - y;

        return y;
    }

    // 1.25 < |x| < 107
    if (ix < 0x4005ac00) {
        // x < -9
        if ((ix >= 0x40022000) and (@as(u32, @bitCast(sign)) & 0x80000000) != 0)
            return 2 - erf_data.tiny_128;

        const xx: f128 = math.abs(x);
        var z: f128 = 1 / (xx * xx);
        const i: i32 = cast(i32, 8.0 / xx, .{});
        var p: f128 = undefined;
        switch (i) {
            1 => {
                p = erf_data.neval(z, &erf_data.RNr2_128, 11) / erf_data.deval(z, &erf_data.RDr2_128, 10);
            },
            2 => {
                p = erf_data.neval(z, &erf_data.RNr3_128, 11) / erf_data.deval(z, &erf_data.RDr3_128, 10);
            },
            3 => {
                p = erf_data.neval(z, &erf_data.RNr4_128, 10) / erf_data.deval(z, &erf_data.RDr4_128, 10);
            },
            4 => {
                p = erf_data.neval(z, &erf_data.RNr5_128, 10) / erf_data.deval(z, &erf_data.RDr5_128, 9);
            },
            5 => {
                p = erf_data.neval(z, &erf_data.RNr6_128, 9) / erf_data.deval(z, &erf_data.RDr6_128, 9);
            },
            6 => {
                p = erf_data.neval(z, &erf_data.RNr7_128, 9) / erf_data.deval(z, &erf_data.RDr7_128, 9);
            },
            7 => {
                p = erf_data.neval(z, &erf_data.RNr8_128, 9) / erf_data.deval(z, &erf_data.RDr8_128, 8);
            },
            else => {
                p = erf_data.neval(z, &erf_data.RNr1_128, 9) / erf_data.deval(z, &erf_data.RDr1_128, 8);
            },
        }
        u = @bitCast(xx);
        u.w3 = 0;
        u.w2 &= 0xfe000000;
        z = @bitCast(u);
        const r: f128 = math.exp(-z * z - 0.5625) * math.exp((z - xx) * (z + xx) + p);

        if ((@as(u32, @bitCast(sign)) & 0x80000000) == 0) {
            const ret: f128 = r / xx;
            return ret;
        } else {
            return 2 - r / xx;
        }
    } else {
        if ((@as(u32, @bitCast(sign)) & 0x80000000) == 0) {
            return erf_data.tiny_128 * erf_data.tiny_128;
        } else {
            return 2 - erf_data.tiny_128;
        }
    }
}

test erfc {
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, -0x2p-56)));
    try std.testing.expectEqual(0xd.c143dp-4, erfc(@as(f32, 0x2p-4)));
    try std.testing.expectEqual(0x4.9f1b48p-4, erfc(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x2.844c2cp-4, erfc(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.d7bb3ep+0, erfc(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x1.3bcd14p-4, erfc(@as(f32, 0x1.4p+0)));
    try std.testing.expectEqual(0x1.328f5ep-8, erfc(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(0x1.fecd7p+0, erfc(@as(f32, -0x2p+0)));
    try std.testing.expectEqual(0x1.729df6p-16, erfc(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(0x1.fffe8ep+0, erfc(@as(f32, -0x3p+0)));
    try std.testing.expectEqual(0x7.4334a8p-28, erfc(@as(f32, 0x3.ee6078p+0)));
    try std.testing.expectEqual(0x4.237748p-28, erfc(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f32, -0x4p+0)));
    try std.testing.expectEqual(0x1.74b17ap-28, erfc(@as(f32, 0x4.2p+0)));
    try std.testing.expectEqual(0x1.b0c1a8p-40, erfc(@as(f32, 0x5p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f32, -0x5p+0)));
    try std.testing.expectEqual(0x1.8cf816p-56, erfc(@as(f32, 0x6p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f32, -0x6p+0)));
    try std.testing.expectEqual(0x3.294504p-76, erfc(@as(f32, 0x7p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f32, -0x7p+0)));
    try std.testing.expectEqual(0xe.3a7e2p-100, erfc(@as(f32, 0x8p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f32, -0x8p+0)));
    try std.testing.expectEqual(0x8.cc6a1p-124, erfc(@as(f32, 0x9p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f32, -0x9p+0)));
    try std.testing.expectEqual(0x8p-152, erfc(@as(f32, 0xap+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f32, -0xap+0)));
    try std.testing.expectEqual(0xf.a3372p-100, erfc(@as(f32, 0x7.fe8008p+0)));
    try std.testing.expectEqual(0xe.3b46ep-100, erfc(@as(f32, 0x7.ffff2p+0)));
    try std.testing.expectEqual(0x1.853f7ap+0, erfc(@as(f32, -0x7.fffff8p-4)));
    try std.testing.expectEqual(0x1.853f7ap+0, erfc(@as(f32, -0x8p-4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.ap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.bp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.cp+4)));
    try std.testing.expectEqual(0xe.3cd88p-100, erfc(@as(f32, 0x7.fffd6p+0)));
    try std.testing.expectEqual(0xe.3cdfbp-100, erfc(@as(f32, 0x7.fffd58p+0)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.4p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.ap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a8p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.aap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.bp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.cp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x3.e8p+8)));
    try std.testing.expectEqual(0xf.6f9d5p-4, erfc(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x1.09062cp+0, erfc(@as(f32, -0x8p-8)));
    try std.testing.expectEqual(0xf.fb7c9p-4, erfc(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0xf.ffdbep-4, erfc(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffeep-4, erfc(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffffp-4, erfc(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x4.000418p-128, erfc(@as(f32, 0x9.31cdfp+0)));
    try std.testing.expectEqual(0x3.ffff78p-128, erfc(@as(f32, 0x9.31cep+0)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.a8b13p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.a8b12ep+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.a8b13p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.a8b12ep+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.9d7adcp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.9d7adap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.9d7adcp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x1.9d7adap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a893p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a893p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a8a058p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a8a05p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a8a058p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a8a05p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a893p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0x6.a893p+4)));
    try std.testing.expectEqual(0x3.fff91cp-4, erfc(@as(f32, 0xd.03d06p-4)));
    try std.testing.expectEqual(0xd.cc226p-8, erfc(@as(f32, 0x1.5cf218p+0)));
    try std.testing.expectEqual(0xd.cc22cp-8, erfc(@as(f32, 0x1.5cf216p+0)));
    try std.testing.expectEqual(0xf.f53dp-8, erfc(@as(f32, 0x1.5166e2p+0)));
    try std.testing.expectEqual(0xf.f53d7p-8, erfc(@as(f32, 0x1.5166ep+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x7.8cde2p-8, erfc(@as(f32, 0x1.8a0c64p+0)));
    try std.testing.expectEqual(0x7.8cde58p-8, erfc(@as(f32, 0x1.8a0c62p+0)));
    try std.testing.expectEqual(0xc.766ccp-8, erfc(@as(f32, 0x1.64dafap+0)));
    try std.testing.expectEqual(0x7.23ff78p-68, erfc(@as(f32, 0x6.88fb08p+0)));
    try std.testing.expectEqual(0x3.e2fa6p-4, erfc(@as(f32, 0xd.361d9p-4)));
    try std.testing.expectEqual(0x1.eb9636p-116, erfc(@as(f32, 0x8.c66b5p+0)));
    try std.testing.expectEqual(0x1.eb9854p-116, erfc(@as(f32, 0x8.c66b4p+0)));
    try std.testing.expectEqual(0x3.ba3ac4p-12, erfc(@as(f32, 0x2.586f1cp+0)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0xb.acb72p+0)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0xb.2274ap+0)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f32, 0xb.22749p+0)));
    try std.testing.expectEqual(0x3.eaab98p-4, erfc(@as(f32, 0xd.28abfp-4)));
    try std.testing.expectEqual(0xf.bbc04p-8, erfc(@as(f32, 0x1.5289fep+0)));
    try std.testing.expectEqual(0x1.f57facp-36, erfc(@as(f32, 0x4.b48498p+0)));
    try std.testing.expectEqual(0x1.be98dep-16, erfc(@as(f32, 0x2.f8646cp+0)));
    try std.testing.expectEqual(0xf.fbeaep-8, erfc(@as(f32, 0x1.514548p+0)));
    try std.testing.expectEqual(0x7.22d058p-12, erfc(@as(f32, 0x2.36c504p+0)));
    try std.testing.expectEqual(0xc.4bf9ep-8, erfc(@as(f32, 0x1.65e31p+0)));
    try std.testing.expectEqual(0x3.da9f6p-4, erfc(@as(f32, 0xd.44cd3p-4)));
    try std.testing.expectEqual(0x3.d93aa4p-4, erfc(@as(f32, 0xd.47426p-4)));
    try std.testing.expectEqual(0x3.d93abp-4, erfc(@as(f32, 0xd.47425p-4)));
    try std.testing.expectEqual(0x1.7fefcp-4, erfc(@as(f32, 0x1.2f644ep+0)));
    try std.testing.expectEqual(0x3.dbca04p-12, erfc(@as(f32, 0x2.56af04p+0)));
    try std.testing.expectEqual(0x7.e8b2fp-16, erfc(@as(f32, 0x2.b7f8ccp+0)));
    try std.testing.expectEqual(0x7.e8b3a8p-16, erfc(@as(f32, 0x2.b7f8c8p+0)));
    try std.testing.expectEqual(0x3.281c2cp-16, erfc(@as(f32, 0x2.dfb9b4p+0)));
    try std.testing.expectEqual(0x1.f1cb04p-8, erfc(@as(f32, 0x1.e33c9ep+0)));
    try std.testing.expectEqual(0x1.3bd96p-4, erfc(@as(f32, 0x1.3ffccp+0)));
    try std.testing.expectEqual(0x1.3bd968p-4, erfc(@as(f32, 0x1.3ffcbep+0)));

    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, -0x2p-56)));
    try std.testing.expectEqual(0xd.c143cb94788fp-4, erfc(@as(f64, 0x2p-4)));
    try std.testing.expectEqual(0x4.9f1b453178d04p-4, erfc(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x2.844c2c5f7bbaap-4, erfc(@as(f64, 0x1p+0)));
    // try std.testing.expectEqual(0x1.d7bb3d3a08445p+0, erfc(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x1.3bcd133aa0ffcp-4, erfc(@as(f64, 0x1.4p+0)));
    // try std.testing.expectEqual(0x1.328f5ec350e67p-8, erfc(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(0x1.fecd70a13caf2p+0, erfc(@as(f64, -0x2p+0)));
    // try std.testing.expectEqual(0x1.729df6503422ap-16, erfc(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(0x1.fffe8d6209afdp+0, erfc(@as(f64, -0x3p+0)));
    // try std.testing.expectEqual(0x7.4334a54e1208cp-28, erfc(@as(f64, 0x3.ee6078p+0)));
    try std.testing.expectEqual(0x4.237744ef4d79cp-28, erfc(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(0x1.ffffffbdc88bbp+0, erfc(@as(f64, -0x4p+0)));
    // try std.testing.expectEqual(0x1.74b179d1eba81p-28, erfc(@as(f64, 0x4.2p+0)));
    // try std.testing.expectEqual(0x1.b0c1a759f7739p-40, erfc(@as(f64, 0x5p+0)));
    try std.testing.expectEqual(0x1.fffffffffe4f4p+0, erfc(@as(f64, -0x5p+0)));
    // try std.testing.expectEqual(0x1.8cf81557d20b6p-56, erfc(@as(f64, 0x6p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f64, -0x6p+0)));
    try std.testing.expectEqual(0x3.2945026df4e62p-76, erfc(@as(f64, 0x7p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f64, -0x7p+0)));
    // try std.testing.expectEqual(0xe.3a7e2090befd8p-100, erfc(@as(f64, 0x8p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f64, -0x8p+0)));
    try std.testing.expectEqual(0x8.cc6a115f1fc6p-124, erfc(@as(f64, 0x9p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f64, -0x9p+0)));
    try std.testing.expectEqual(0xb.ec53f9545168p-152, erfc(@as(f64, 0xap+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f64, -0xap+0)));
    try std.testing.expectEqual(0xf.a33725bea2f8p-100, erfc(@as(f64, 0x7.fe8008p+0)));
    try std.testing.expectEqual(0xe.3b46e15ad978p-100, erfc(@as(f64, 0x7.ffff2p+0)));
    try std.testing.expectEqual(0x1.853f7a704b7bep+0, erfc(@as(f64, -0x7.fffff8p-4)));
    try std.testing.expectEqual(0x1.853f7ae0c76e9p+0, erfc(@as(f64, -0x8p-4)));
    try std.testing.expectEqual(0x1.853f7ae0c76e9p+0, erfc(@as(f64, -0x7.ffffffffffffcp-4)));
    try std.testing.expectEqual(0x9.425ff0e6f512p-984, erfc(@as(f64, 0x1.ap+4)));
    try std.testing.expectEqual(0x6.783cp-1060, erfc(@as(f64, 0x1.bp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x1.cp+4)));
    try std.testing.expectEqual(0xe.3cd883e02b15p-100, erfc(@as(f64, 0x7.fffd6p+0)));
    try std.testing.expectEqual(0xe.3cdfb051e694p-100, erfc(@as(f64, 0x7.fffd58p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe84cp-100, erfc(@as(f64, 0x7.fffd59e26af38p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe8858p-100, erfc(@as(f64, 0x7.fffd59e26af34p+0)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.4p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.ap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.aap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.bp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.cp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x3.e8p+8)));
    // try std.testing.expectEqual(0xf.6f9d4dd116d68p-4, erfc(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0x1.09062b22ee92ap+0, erfc(@as(f64, -0x8p-8)));
    try std.testing.expectEqual(0xf.fb7c8a4401cdp-4, erfc(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0xf.ffdbe4515fabp-4, erfc(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffedf228afcp-4, erfc(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffff6f91458p-4, erfc(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0xf.ffffffb7c8a28p-4, erfc(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0xf.fffffffdbe45p-4, erfc(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0xf.ffffffffedf2p-4, erfc(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0xf.ffffffffff6f8p-4, erfc(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0xf.fffffffffffb8p-4, erfc(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x4p-1076)));
    // try std.testing.expectEqual(0x4.0004157f2239cp-128, erfc(@as(f64, 0x9.31cdfp+0)));
    try std.testing.expectEqual(0x3.ffff75b4a7f72p-128, erfc(@as(f64, 0x9.31cep+0)));
    try std.testing.expectEqual(0x3.fffd098f7c63cp-1024, erfc(@as(f64, 0x1.a8b13p+4)));
    // try std.testing.expectEqual(0x4.001799b7b63bcp-1024, erfc(@as(f64, 0x1.a8b12ep+4)));
    try std.testing.expectEqual(0x4.0000000000cc8p-1024, erfc(@as(f64, 0x1.a8b12fc6e4891p+4)));
    try std.testing.expectEqual(0x3.fffd098f7c63cp-1024, erfc(@as(f64, 0x1.a8b13p+4)));
    // try std.testing.expectEqual(0x4.001799b7b63bcp-1024, erfc(@as(f64, 0x1.a8b12ep+4)));
    // try std.testing.expectEqual(0x3.fffffffffff8p-1024, erfc(@as(f64, 0x1.a8b12fc6e4892p+4)));
    try std.testing.expectEqual(0x7.ffe0488939958p-972, erfc(@as(f64, 0x1.9d7adcp+4)));
    try std.testing.expectEqual(0x8.001401a2efa28p-972, erfc(@as(f64, 0x1.9d7adap+4)));
    try std.testing.expectEqual(0x7.ffffffffff3bp-972, erfc(@as(f64, 0x1.9d7adac608e86p+4)));
    try std.testing.expectEqual(0x8.0000000000d9p-972, erfc(@as(f64, 0x1.9d7adac608e85p+4)));
    try std.testing.expectEqual(0x7.ffe0488939958p-972, erfc(@as(f64, 0x1.9d7adcp+4)));
    try std.testing.expectEqual(0x8.001401a2efa28p-972, erfc(@as(f64, 0x1.9d7adap+4)));
    try std.testing.expectEqual(0x7.ffffffffff3bp-972, erfc(@as(f64, 0x1.9d7adac608e86p+4)));
    try std.testing.expectEqual(0x8.0000000000d9p-972, erfc(@as(f64, 0x1.9d7adac608e85p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893032db9054p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893032db905p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893032db9054p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893032db905p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8a058p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8a05p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8a0561d8bbecp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8a0561d8bbe8p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8a058p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8a05p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8a0561d8bbecp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a8a0561d8bbe8p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893032db9054p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893032db905p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893032db9054p+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0x6.a893032db905p+4)));
    try std.testing.expectEqual(0x3.fff91a7d782bp-4, erfc(@as(f64, 0xd.03d06p-4)));
    // try std.testing.expectEqual(0xd.cc22642cb5ab8p-8, erfc(@as(f64, 0x1.5cf218p+0)));
    try std.testing.expectEqual(0xd.cc22be4b9b328p-8, erfc(@as(f64, 0x1.5cf216p+0)));
    // try std.testing.expectEqual(0xd.cc22a7f1317fp-8, erfc(@as(f64, 0x1.5cf2167efe921p+0)));
    // try std.testing.expectEqual(0xd.cc22a7f131818p-8, erfc(@as(f64, 0x1.5cf2167efe92p+0)));
    try std.testing.expectEqual(0xf.f53d075aa92bp-8, erfc(@as(f64, 0x1.5166e2p+0)));
    // try std.testing.expectEqual(0xf.f53d6d0e58d08p-8, erfc(@as(f64, 0x1.5166ep+0)));
    try std.testing.expectEqual(0xf.f53d3d6dfa74p-8, erfc(@as(f64, 0x1.5166e0efc44aap+0)));
    // try std.testing.expectEqual(0xf.f53d3d6dfa778p-8, erfc(@as(f64, 0x1.5166e0efc44a9p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f64, -0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(0x7.8cde235791e7cp-8, erfc(@as(f64, 0x1.8a0c64p+0)));
    // try std.testing.expectEqual(0x7.8cde596304018p-8, erfc(@as(f64, 0x1.8a0c62p+0)));
    try std.testing.expectEqual(0xc.766cbf61fd648p-8, erfc(@as(f64, 0x1.64dafap+0)));
    // try std.testing.expectEqual(0x7.23ff79ae0f25cp-68, erfc(@as(f64, 0x6.88fb08p+0)));
    try std.testing.expectEqual(0x3.e2fa6064d5894p-4, erfc(@as(f64, 0xd.361d9p-4)));
    try std.testing.expectEqual(0x1.eb9635bc51eb8p-116, erfc(@as(f64, 0x8.c66b5p+0)));
    try std.testing.expectEqual(0x1.eb98546946cb2p-116, erfc(@as(f64, 0x8.c66b4p+0)));
    try std.testing.expectEqual(0x1.eb97b1f20867cp-116, erfc(@as(f64, 0x8.c66b44ca40038p+0)));
    // try std.testing.expectEqual(0x3.ba3ac339ed19p-12, erfc(@as(f64, 0x2.586f1cp+0)));
    try std.testing.expectEqual(0x7.ee2d2ec57315p-204, erfc(@as(f64, 0xb.acb72p+0)));
    // try std.testing.expectEqual(0x1.c646841c9021p-184, erfc(@as(f64, 0xb.2274ap+0)));
    // try std.testing.expectEqual(0x1.c648feeb672e9p-184, erfc(@as(f64, 0xb.22749p+0)));
    // try std.testing.expectEqual(0x1.c6479753ddcb5p-184, erfc(@as(f64, 0xb.227499103358p+0)));
    try std.testing.expectEqual(0x1.c6479753dddf2p-184, erfc(@as(f64, 0xb.2274991033578p+0)));
    // try std.testing.expectEqual(0x3.eaab96d5a2e2ap-4, erfc(@as(f64, 0xd.28abfp-4)));
    // try std.testing.expectEqual(0xf.bbc04428a3d3p-8, erfc(@as(f64, 0x1.5289fep+0)));
    // try std.testing.expectEqual(0x1.f57fab6c3db3dp-36, erfc(@as(f64, 0x4.b48498p+0)));
    try std.testing.expectEqual(0x1.be98de114e175p-16, erfc(@as(f64, 0x2.f8646cp+0)));
    try std.testing.expectEqual(0xf.fbeadad5a51f8p-8, erfc(@as(f64, 0x1.514548p+0)));
    try std.testing.expectEqual(0x7.22d059993f3f4p-12, erfc(@as(f64, 0x2.36c504p+0)));
    try std.testing.expectEqual(0xc.4bf9de451e6p-8, erfc(@as(f64, 0x1.65e31p+0)));
    // try std.testing.expectEqual(0x3.da9f608f1dd7ep-4, erfc(@as(f64, 0xd.44cd3p-4)));
    try std.testing.expectEqual(0x3.d93aa59c8f5acp-4, erfc(@as(f64, 0xd.47426p-4)));
    try std.testing.expectEqual(0x3.d93aaeadb64dp-4, erfc(@as(f64, 0xd.47425p-4)));
    try std.testing.expectEqual(0x3.d93aa84f87aap-4, erfc(@as(f64, 0xd.47425b3cafa48p-4)));
    // try std.testing.expectEqual(0x1.7fefc09137c95p-4, erfc(@as(f64, 0x1.2f644ep+0)));
    // try std.testing.expectEqual(0x3.dbca059c7e73ap-12, erfc(@as(f64, 0x2.56af04p+0)));
    try std.testing.expectEqual(0x7.e8b2efb67945p-16, erfc(@as(f64, 0x2.b7f8ccp+0)));
    try std.testing.expectEqual(0x7.e8b3a6276f04p-16, erfc(@as(f64, 0x2.b7f8c8p+0)));
    try std.testing.expectEqual(0x7.e8b308381dfc4p-16, erfc(@as(f64, 0x2.b7f8cb76737d4p+0)));
    try std.testing.expectEqual(0x7.e8b308381e02p-16, erfc(@as(f64, 0x2.b7f8cb76737d2p+0)));
    // try std.testing.expectEqual(0x3.281c2d7e470e6p-16, erfc(@as(f64, 0x2.dfb9b4p+0)));
    // try std.testing.expectEqual(0x1.f1cb04b622e6fp-8, erfc(@as(f64, 0x1.e33c9ep+0)));
    // try std.testing.expectEqual(0x1.3bd95ffe4e556p-4, erfc(@as(f64, 0x1.3ffccp+0)));
    // try std.testing.expectEqual(0x1.3bd9679020a68p-4, erfc(@as(f64, 0x1.3ffcbep+0)));
    // try std.testing.expectEqual(0x1.3bd962ebb7736p-4, erfc(@as(f64, 0x1.3ffcbf39febb4p+0)));

    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffdbep-4, erfc(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1.0000000000000242p+0, erfc(@as(f80, -0x2p-56)));
    try std.testing.expectEqual(0xd.c143cb94788ed18p-4, erfc(@as(f80, 0x2p-4)));
    try std.testing.expectEqual(0x4.9f1b453178d049d8p-4, erfc(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x2.844c2c5f7bba9c98p-4, erfc(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.d7bb3d3a08445636p+0, erfc(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x1.3bcd133aa0ffbf9ep-4, erfc(@as(f80, 0x1.4p+0)));
    try std.testing.expectEqual(0x1.328f5ec350e668d8p-8, erfc(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(0x1.fecd70a13caf1998p+0, erfc(@as(f80, -0x2p+0)));
    try std.testing.expectEqual(0x1.729df6503422a0bcp-16, erfc(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(0x1.fffe8d6209afcbdep+0, erfc(@as(f80, -0x3p+0)));
    try std.testing.expectEqual(0x7.4334a54e1208ae18p-28, erfc(@as(f80, 0x3.ee6078p+0)));
    try std.testing.expectEqual(0x4.237744ef4d79a9e8p-28, erfc(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(0x1.ffffffbdc88bb10cp+0, erfc(@as(f80, -0x4p+0)));
    try std.testing.expectEqual(0x1.74b179d1eba809f2p-28, erfc(@as(f80, 0x4.2p+0)));
    try std.testing.expectEqual(0x1.b0c1a759f7738936p-40, erfc(@as(f80, 0x5p+0)));
    try std.testing.expectEqual(0x1.fffffffffe4f3e58p+0, erfc(@as(f80, -0x5p+0)));
    try std.testing.expectEqual(0x1.8cf81557d20b61a8p-56, erfc(@as(f80, 0x6p+0)));
    try std.testing.expectEqual(0x1.fffffffffffffe74p+0, erfc(@as(f80, -0x6p+0)));
    try std.testing.expectEqual(0x3.2945026df4e62a48p-76, erfc(@as(f80, 0x7p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f80, -0x7p+0)));
    try std.testing.expectEqual(0xe.3a7e2090befdbb6p-100, erfc(@as(f80, 0x8p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f80, -0x8p+0)));
    try std.testing.expectEqual(0x8.cc6a115f1fc6137p-124, erfc(@as(f80, 0x9p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f80, -0x9p+0)));
    try std.testing.expectEqual(0xb.ec53f9545167ceap-152, erfc(@as(f80, 0xap+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f80, -0xap+0)));
    try std.testing.expectEqual(0xf.a33725bea2f7d7bp-100, erfc(@as(f80, 0x7.fe8008p+0)));
    try std.testing.expectEqual(0xe.3b46e15ad97825dp-100, erfc(@as(f80, 0x7.ffff2p+0)));
    try std.testing.expectEqual(0x1.853f7a704b7be2d6p+0, erfc(@as(f80, -0x7.fffff8p-4)));
    try std.testing.expectEqual(0x1.853f7ae0c76e915ep+0, erfc(@as(f80, -0x8p-4)));
    try std.testing.expectEqual(0x1.853f7ae0c76e8ddap+0, erfc(@as(f80, -0x7.ffffffffffffcp-4)));
    try std.testing.expectEqual(0x1.853f7ae0c76e8f9cp+0, erfc(@as(f80, -0x7.ffffffffffffep-4)));
    try std.testing.expectEqual(0x9.425ff0e6f511d75p-984, erfc(@as(f80, 0x1.ap+4)));
    try std.testing.expectEqual(0x6.783c337e0e9d7e88p-1060, erfc(@as(f80, 0x1.bp+4)));
    try std.testing.expectEqual(0x9.cd4b80875a8ec66p-1140, erfc(@as(f80, 0x1.cp+4)));
    try std.testing.expectEqual(0xe.3cd883e02b14dbp-100, erfc(@as(f80, 0x7.fffd6p+0)));
    try std.testing.expectEqual(0xe.3cdfb051e694315p-100, erfc(@as(f80, 0x7.fffd58p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe84be02p-100, erfc(@as(f80, 0x7.fffd59e26af38p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe88543cp-100, erfc(@as(f80, 0x7.fffd59e26af34p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe84fa89p-100, erfc(@as(f80, 0x7.fffd59e26af37bc8p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe84fafcp-100, erfc(@as(f80, 0x7.fffd59e26af37bcp+0)));
    try std.testing.expectEqual(0x2.fd514cef7750e588p-14436, erfc(@as(f80, 0x6.4p+4)));
    try std.testing.expectEqual(0x5.028a2f1656a432d8p-16220, erfc(@as(f80, 0x6.ap+4)));
    try std.testing.expectEqual(0x2.0b5b5b3bbf7d96a4p-16372, erfc(@as(f80, 0x6.a8p+4)));
    // try std.testing.expectEqual(0x6.0b6ee998p-16412, erfc(@as(f80, 0x6.aap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f80, 0x6.bp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f80, 0x6.cp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f80, 0x3.e8p+8)));
    try std.testing.expectEqual(0xf.6f9d4dd116d6403p-4, erfc(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x1.09062b22ee929bfcp+0, erfc(@as(f80, -0x8p-8)));
    try std.testing.expectEqual(0xf.fb7c8a4401cd1c3p-4, erfc(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0xf.ffdbe4515faaee1p-4, erfc(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffedf228afbd6bp-4, erfc(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffff6f91457debp-4, erfc(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0xf.ffffffb7c8a2befp-4, erfc(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0xf.fffffffdbe4515fp-4, erfc(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0xf.ffffffffedf228bp-4, erfc(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0xf.ffffffffff6f914p-4, erfc(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0xf.fffffffffffb7c9p-4, erfc(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffffffeep-4, erfc(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x4.0004157f2239d72p-128, erfc(@as(f80, 0x9.31cdfp+0)));
    try std.testing.expectEqual(0x3.ffff75b4a7f7172p-128, erfc(@as(f80, 0x9.31cep+0)));
    try std.testing.expectEqual(0x3.fffd098f7c63a42cp-1024, erfc(@as(f80, 0x1.a8b13p+4)));
    try std.testing.expectEqual(0x4.001799b7b63bbfp-1024, erfc(@as(f80, 0x1.a8b12ep+4)));
    try std.testing.expectEqual(0x4.0000000000cc9078p-1024, erfc(@as(f80, 0x1.a8b12fc6e4891p+4)));
    try std.testing.expectEqual(0x3.fffd098f7c63a42cp-1024, erfc(@as(f80, 0x1.a8b13p+4)));
    try std.testing.expectEqual(0x4.001799b7b63bbfp-1024, erfc(@as(f80, 0x1.a8b12ep+4)));
    try std.testing.expectEqual(0x3.fffffffffff8115cp-1024, erfc(@as(f80, 0x1.a8b12fc6e4892p+4)));
    try std.testing.expectEqual(0x7.ffe048893995704p-972, erfc(@as(f80, 0x1.9d7adcp+4)));
    try std.testing.expectEqual(0x8.001401a2efa2625p-972, erfc(@as(f80, 0x1.9d7adap+4)));
    try std.testing.expectEqual(0x7.ffffffffff3b1b68p-972, erfc(@as(f80, 0x1.9d7adac608e86p+4)));
    try std.testing.expectEqual(0x8.0000000000d8e56p-972, erfc(@as(f80, 0x1.9d7adac608e85p+4)));
    try std.testing.expectEqual(0x7.ffffffffffffe638p-972, erfc(@as(f80, 0x1.9d7adac608e85864p+4)));
    try std.testing.expectEqual(0x8.00000000000019fp-972, erfc(@as(f80, 0x1.9d7adac608e85862p+4)));
    try std.testing.expectEqual(0x7.ffe048893995704p-972, erfc(@as(f80, 0x1.9d7adcp+4)));
    try std.testing.expectEqual(0x8.001401a2efa2625p-972, erfc(@as(f80, 0x1.9d7adap+4)));
    try std.testing.expectEqual(0x7.ffffffffff3b1b68p-972, erfc(@as(f80, 0x1.9d7adac608e86p+4)));
    try std.testing.expectEqual(0x8.0000000000d8e56p-972, erfc(@as(f80, 0x1.9d7adac608e85p+4)));
    try std.testing.expectEqual(0x7.ffffffffffffe638p-972, erfc(@as(f80, 0x1.9d7adac608e85864p+4)));
    try std.testing.expectEqual(0x8.00000000000019fp-972, erfc(@as(f80, 0x1.9d7adac608e85862p+4)));
    // try std.testing.expectEqual(0x3.feff49e314f6ecd8p-16384, erfc(@as(f80, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x4.00a9613ff522441p-16384, erfc(@as(f80, 0x6.a893p+4)));
    // try std.testing.expectEqual(0x3.fffffffffad99d6p-16384, erfc(@as(f80, 0x6.a893032db9054p+4)));
    try std.testing.expectEqual(0x4.00000000082ae9d8p-16384, erfc(@as(f80, 0x6.a893032db905p+4)));
    try std.testing.expectEqual(0x4.0000000000000dfp-16384, erfc(@as(f80, 0x6.a893032db905274p+4)));
    // try std.testing.expectEqual(0x3.feff49e314f6ecd8p-16384, erfc(@as(f80, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x4.00a9613ff522441p-16384, erfc(@as(f80, 0x6.a893p+4)));
    // try std.testing.expectEqual(0x3.fffffffffad99d6p-16384, erfc(@as(f80, 0x6.a893032db9054p+4)));
    try std.testing.expectEqual(0x4.00000000082ae9d8p-16384, erfc(@as(f80, 0x6.a893032db905p+4)));
    // try std.testing.expectEqual(0x3.fffffffffffe63c8p-16384, erfc(@as(f80, 0x6.a893032db9052748p+4)));
    // try std.testing.expectEqual(0x1.ffcdcfd4f9515ad8p-16384, erfc(@as(f80, 0x6.a8a058p+4)));
    // try std.testing.expectEqual(0x2.00a2fdbcb5dc48ap-16384, erfc(@as(f80, 0x6.a8a05p+4)));
    // try std.testing.expectEqual(0x1.fffffffffb6f715p-16384, erfc(@as(f80, 0x6.a8a0561d8bbecp+4)));
    // try std.testing.expectEqual(0x2.00000000021824d8p-16384, erfc(@as(f80, 0x6.a8a0561d8bbe8p+4)));
    // try std.testing.expectEqual(0x2.0000000000001868p-16384, erfc(@as(f80, 0x6.a8a0561d8bbe942p+4)));
    // try std.testing.expectEqual(0x1.ffcdcfd4f9515ad8p-16384, erfc(@as(f80, 0x6.a8a058p+4)));
    // try std.testing.expectEqual(0x2.00a2fdbcb5dc48ap-16384, erfc(@as(f80, 0x6.a8a05p+4)));
    // try std.testing.expectEqual(0x1.fffffffffb6f715p-16384, erfc(@as(f80, 0x6.a8a0561d8bbecp+4)));
    // try std.testing.expectEqual(0x2.00000000021824d8p-16384, erfc(@as(f80, 0x6.a8a0561d8bbe8p+4)));
    // try std.testing.expectEqual(0x1.ffffffffffff435p-16384, erfc(@as(f80, 0x6.a8a0561d8bbe9428p+4)));
    // try std.testing.expectEqual(0x3.feff49e314f6ecd8p-16384, erfc(@as(f80, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x4.00a9613ff522441p-16384, erfc(@as(f80, 0x6.a893p+4)));
    // try std.testing.expectEqual(0x3.fffffffffad99d6p-16384, erfc(@as(f80, 0x6.a893032db9054p+4)));
    try std.testing.expectEqual(0x4.00000000082ae9d8p-16384, erfc(@as(f80, 0x6.a893032db905p+4)));
    // try std.testing.expectEqual(0x3.fffffffffffe63c8p-16384, erfc(@as(f80, 0x6.a893032db9052748p+4)));
    try std.testing.expectEqual(0x4.0000000000000dfp-16384, erfc(@as(f80, 0x6.a893032db905274p+4)));
    // try std.testing.expectEqual(0x3.feff49e314f6ecd8p-16384, erfc(@as(f80, 0x6.a89308p+4)));
    try std.testing.expectEqual(0x4.00a9613ff522441p-16384, erfc(@as(f80, 0x6.a893p+4)));
    // try std.testing.expectEqual(0x3.fffffffffad99d6p-16384, erfc(@as(f80, 0x6.a893032db9054p+4)));
    try std.testing.expectEqual(0x4.00000000082ae9d8p-16384, erfc(@as(f80, 0x6.a893032db905p+4)));
    // try std.testing.expectEqual(0x3.fffffffffffe63c8p-16384, erfc(@as(f80, 0x6.a893032db9052748p+4)));
    try std.testing.expectEqual(0x4.0000000000000dfp-16384, erfc(@as(f80, 0x6.a893032db905274p+4)));
    try std.testing.expectEqual(0x3.fff91a7d782b0064p-4, erfc(@as(f80, 0xd.03d06p-4)));
    try std.testing.expectEqual(0xd.cc22642cb5ab8dcp-8, erfc(@as(f80, 0x1.5cf218p+0)));
    try std.testing.expectEqual(0xd.cc22be4b9b325bcp-8, erfc(@as(f80, 0x1.5cf216p+0)));
    try std.testing.expectEqual(0xd.cc22a7f1317ede1p-8, erfc(@as(f80, 0x1.5cf2167efe921p+0)));
    try std.testing.expectEqual(0xd.cc22a7f13181af1p-8, erfc(@as(f80, 0x1.5cf2167efe92p+0)));
    try std.testing.expectEqual(0xd.cc22a7f131804ebp-8, erfc(@as(f80, 0x1.5cf2167efe9207d2p+0)));
    try std.testing.expectEqual(0xf.f53d075aa92b1fp-8, erfc(@as(f80, 0x1.5166e2p+0)));
    try std.testing.expectEqual(0xf.f53d6d0e58d08f8p-8, erfc(@as(f80, 0x1.5166ep+0)));
    try std.testing.expectEqual(0xf.f53d3d6dfa74177p-8, erfc(@as(f80, 0x1.5166e0efc44aap+0)));
    try std.testing.expectEqual(0xf.f53d3d6dfa77451p-8, erfc(@as(f80, 0x1.5166e0efc44a9p+0)));
    try std.testing.expectEqual(0xf.f53d3d6dfa747d8p-8, erfc(@as(f80, 0x1.5166e0efc44a9dfep+0)));
    try std.testing.expectEqual(0xf.f53d3d6dfa747dfp-8, erfc(@as(f80, 0x1.5166e0efc44a9dfcp+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x7.8cde235791e7d0ep-8, erfc(@as(f80, 0x1.8a0c64p+0)));
    try std.testing.expectEqual(0x7.8cde5963040180b8p-8, erfc(@as(f80, 0x1.8a0c62p+0)));
    try std.testing.expectEqual(0xc.766cbf61fd6480bp-8, erfc(@as(f80, 0x1.64dafap+0)));
    try std.testing.expectEqual(0x7.23ff79ae0f25a138p-68, erfc(@as(f80, 0x6.88fb08p+0)));
    try std.testing.expectEqual(0x3.e2fa6064d589347cp-4, erfc(@as(f80, 0xd.361d9p-4)));
    try std.testing.expectEqual(0x1.eb9635bc51eb7a94p-116, erfc(@as(f80, 0x8.c66b5p+0)));
    try std.testing.expectEqual(0x1.eb98546946cb2526p-116, erfc(@as(f80, 0x8.c66b4p+0)));
    try std.testing.expectEqual(0x1.eb97b1f20867c35ep-116, erfc(@as(f80, 0x8.c66b44ca40038p+0)));
    try std.testing.expectEqual(0x3.ba3ac339ed190204p-12, erfc(@as(f80, 0x2.586f1cp+0)));
    try std.testing.expectEqual(0x7.ee2d2ec5731504p-204, erfc(@as(f80, 0xb.acb72p+0)));
    try std.testing.expectEqual(0x1.c646841c902106e8p-184, erfc(@as(f80, 0xb.2274ap+0)));
    try std.testing.expectEqual(0x1.c648feeb672e8e58p-184, erfc(@as(f80, 0xb.22749p+0)));
    try std.testing.expectEqual(0x1.c6479753ddcb4d8cp-184, erfc(@as(f80, 0xb.227499103358p+0)));
    try std.testing.expectEqual(0x1.c6479753dddf2402p-184, erfc(@as(f80, 0xb.2274991033578p+0)));
    try std.testing.expectEqual(0x1.c6479753ddd176a6p-184, erfc(@as(f80, 0xb.227499103357d84p+0)));
    try std.testing.expectEqual(0x3.eaab96d5a2e294b8p-4, erfc(@as(f80, 0xd.28abfp-4)));
    try std.testing.expectEqual(0xf.bbc04428a3d30e7p-8, erfc(@as(f80, 0x1.5289fep+0)));
    try std.testing.expectEqual(0x1.f57fab6c3db3ce7ep-36, erfc(@as(f80, 0x4.b48498p+0)));
    try std.testing.expectEqual(0x1.be98de114e174b5p-16, erfc(@as(f80, 0x2.f8646cp+0)));
    try std.testing.expectEqual(0xf.fbeadad5a51f775p-8, erfc(@as(f80, 0x1.514548p+0)));
    try std.testing.expectEqual(0x7.22d059993f3f46dp-12, erfc(@as(f80, 0x2.36c504p+0)));
    try std.testing.expectEqual(0xc.4bf9de451e5fceep-8, erfc(@as(f80, 0x1.65e31p+0)));
    try std.testing.expectEqual(0x3.da9f608f1dd7ee3p-4, erfc(@as(f80, 0xd.44cd3p-4)));
    try std.testing.expectEqual(0x3.d93aa59c8f5abb84p-4, erfc(@as(f80, 0xd.47426p-4)));
    try std.testing.expectEqual(0x3.d93aaeadb64d00e8p-4, erfc(@as(f80, 0xd.47425p-4)));
    try std.testing.expectEqual(0x3.d93aa84f87a9ffap-4, erfc(@as(f80, 0xd.47425b3cafa48p-4)));
    try std.testing.expectEqual(0x1.7fefc09137c9485ep-4, erfc(@as(f80, 0x1.2f644ep+0)));
    try std.testing.expectEqual(0x3.dbca059c7e73a124p-12, erfc(@as(f80, 0x2.56af04p+0)));
    try std.testing.expectEqual(0x7.e8b2efb679451a4p-16, erfc(@as(f80, 0x2.b7f8ccp+0)));
    try std.testing.expectEqual(0x7.e8b3a6276f03f778p-16, erfc(@as(f80, 0x2.b7f8c8p+0)));
    try std.testing.expectEqual(0x7.e8b308381dfc55c8p-16, erfc(@as(f80, 0x2.b7f8cb76737d4p+0)));
    try std.testing.expectEqual(0x7.e8b308381e02095p-16, erfc(@as(f80, 0x2.b7f8cb76737d2p+0)));
    try std.testing.expectEqual(0x7.e8b308381e001448p-16, erfc(@as(f80, 0x2.b7f8cb76737d2afcp+0)));
    try std.testing.expectEqual(0x7.e8b308381e0015p-16, erfc(@as(f80, 0x2.b7f8cb76737d2af8p+0)));
    try std.testing.expectEqual(0x3.281c2d7e470e5084p-16, erfc(@as(f80, 0x2.dfb9b4p+0)));
    try std.testing.expectEqual(0x1.f1cb04b622e6f4d6p-8, erfc(@as(f80, 0x1.e33c9ep+0)));
    try std.testing.expectEqual(0x1.3bd95ffe4e5561c6p-4, erfc(@as(f80, 0x1.3ffccp+0)));
    try std.testing.expectEqual(0x1.3bd9679020a687cp-4, erfc(@as(f80, 0x1.3ffcbep+0)));
    try std.testing.expectEqual(0x1.3bd962ebb773644cp-4, erfc(@as(f80, 0x1.3ffcbf39febb4p+0)));

    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffdbe4515f7ac925ep-4, erfc(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1.0000000000000241baea08536da2p+0, erfc(@as(f128, -0x2p-56)));
    try std.testing.expectEqual(0xd.c143cb94788ed17a494db60f862p-4, erfc(@as(f128, 0x2p-4)));
    try std.testing.expectEqual(0x4.9f1b453178d049d79a1a68775594p-4, erfc(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x2.844c2c5f7bba9c97f778122796c8p-4, erfc(@as(f128, 0x1p+0)));
    // try std.testing.expectEqual(0x1.d7bb3d3a0844563680887edd8693p+0, erfc(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x1.3bcd133aa0ffbf9d895f72e9b1d3p-4, erfc(@as(f128, 0x1.4p+0)));
    try std.testing.expectEqual(0x1.328f5ec350e668d7fe6fb465cc11p-8, erfc(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(0x1.fecd70a13caf19972801904b9a34p+0, erfc(@as(f128, -0x2p+0)));
    // try std.testing.expectEqual(0x1.729df6503422a0bc26526214d0a4p-16, erfc(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(0x1.fffe8d6209afcbdd5f43d9ad9debp+0, erfc(@as(f128, -0x3p+0)));
    try std.testing.expectEqual(0x7.4334a54e1208ae1b8bfa15647bc4p-28, erfc(@as(f128, 0x3.ee6078p+0)));
    try std.testing.expectEqual(0x4.237744ef4d79a9ea24bfce6c7e8cp-28, erfc(@as(f128, 0x4p+0)));
    try std.testing.expectEqual(0x1.ffffffbdc88bb10b2865615db403p+0, erfc(@as(f128, -0x4p+0)));
    try std.testing.expectEqual(0x1.74b179d1eba809f2e32224074102p-28, erfc(@as(f128, 0x4.2p+0)));
    try std.testing.expectEqual(0x1.b0c1a759f7738935ea5dea8e17aap-40, erfc(@as(f128, 0x5p+0)));
    try std.testing.expectEqual(0x1.fffffffffe4f3e58a6088c76ca16p+0, erfc(@as(f128, -0x5p+0)));
    // try std.testing.expectEqual(0x1.8cf81557d20b61a7fff0cc732bfap-56, erfc(@as(f128, 0x6p+0)));
    try std.testing.expectEqual(0x1.fffffffffffffe7307eaa82df49ep+0, erfc(@as(f128, -0x6p+0)));
    // try std.testing.expectEqual(0x3.2945026df4e62a48fcf382c1cfc8p-76, erfc(@as(f128, 0x7p+0)));
    try std.testing.expectEqual(0x1.ffffffffffffffffffcd6bafd921p+0, erfc(@as(f128, -0x7p+0)));
    try std.testing.expectEqual(0xe.3a7e2090befdbb5c007d16c48e88p-100, erfc(@as(f128, 0x8p+0)));
    try std.testing.expectEqual(0x1.ffffffffffffffffffffffff1c58p+0, erfc(@as(f128, -0x8p+0)));
    try std.testing.expectEqual(0x8.cc6a115f1fc6136ba610a005ff2p-124, erfc(@as(f128, 0x9p+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f128, -0x9p+0)));
    try std.testing.expectEqual(0xb.ec53f9545167ce9b9c460ae3b268p-152, erfc(@as(f128, 0xap+0)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f128, -0xap+0)));
    // try std.testing.expectEqual(0xf.a33725bea2f7d7abe8b7461d621p-100, erfc(@as(f128, 0x7.fe8008p+0)));
    // try std.testing.expectEqual(0xe.3b46e15ad97825d129852878fecp-100, erfc(@as(f128, 0x7.ffff2p+0)));
    try std.testing.expectEqual(0x1.853f7a704b7be2d643b9e3ae3cbp+0, erfc(@as(f128, -0x7.fffff8p-4)));
    // try std.testing.expectEqual(0x1.853f7ae0c76e915e809f1a31a27bp+0, erfc(@as(f128, -0x8p-4)));
    try std.testing.expectEqual(0x1.853f7ae0c76e8ddaa10a86e7a049p+0, erfc(@as(f128, -0x7.ffffffffffffcp-4)));
    try std.testing.expectEqual(0x1.853f7ae0c76e8f9c90d4d08ca164p+0, erfc(@as(f128, -0x7.ffffffffffffep-4)));
    // try std.testing.expectEqual(0x9.425ff0e6f511d74db40cfbbceffp-984, erfc(@as(f128, 0x1.ap+4)));
    // try std.testing.expectEqual(0x6.783c337e0e9d7e84c2c58243308cp-1060, erfc(@as(f128, 0x1.bp+4)));
    // try std.testing.expectEqual(0x9.cd4b80875a8ec6603b9a1f1bead8p-1140, erfc(@as(f128, 0x1.cp+4)));
    try std.testing.expectEqual(0xe.3cd883e02b14daf90f0f812035cp-100, erfc(@as(f128, 0x7.fffd6p+0)));
    try std.testing.expectEqual(0xe.3cdfb051e6943150a3c2f2e70a58p-100, erfc(@as(f128, 0x7.fffd58p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe84be025e570dd9c7cap-100, erfc(@as(f128, 0x7.fffd59e26af38p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe88543bb57f5bfbbaa78p-100, erfc(@as(f128, 0x7.fffd59e26af34p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe84fa8965f5d5ffd3998p-100, erfc(@as(f128, 0x7.fffd59e26af37bc8p+0)));
    try std.testing.expectEqual(0xe.3cddffbbe84fafc2d20bb0997108p-100, erfc(@as(f128, 0x7.fffd59e26af37bcp+0)));
    try std.testing.expectEqual(0xe.3cddffbbe84faf818649c2377138p-100, erfc(@as(f128, 0x7.fffd59e26af37bc048d159e26ap+0)));
    // try std.testing.expectEqual(0x2.fd514cef7750e58906601ff35dcp-14436, erfc(@as(f128, 0x6.4p+4)));
    try std.testing.expectEqual(0x5.028a2f1656a432d79f76a6f2df48p-16220, erfc(@as(f128, 0x6.ap+4)));
    // try std.testing.expectEqual(0x2.0b5b5b3bbf7d96a5e595291fc8c4p-16372, erfc(@as(f128, 0x6.a8p+4)));
    // try std.testing.expectEqual(0x6.0b6ee997d3343b5bf2f08p-16412, erfc(@as(f128, 0x6.aap+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f128, 0x6.bp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f128, 0x6.cp+4)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f128, 0x3.e8p+8)));
    try std.testing.expectEqual(0xf.6f9d4dd116d640363e72a8f031b8p-4, erfc(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(0x1.09062b22ee929bfc9c18d570fce4p+0, erfc(@as(f128, -0x8p-8)));
    try std.testing.expectEqual(0xf.fb7c8a4401cd1c3347a83a17da5p-4, erfc(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0xf.ffdbe4515faaee0eb270b8c6b0dp-4, erfc(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffedf228afbd6a978b6c864b5b8p-4, erfc(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0xf.fffff6f91457deb24a37bc860f7p-4, erfc(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0xf.ffffffb7c8a2bef5924bbac83dp-4, erfc(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0xf.fffffffdbe4515f7ac925dca3bbp-4, erfc(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0xf.ffffffffedf228afbd6492ee51c8p-4, erfc(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0xf.ffffffffff6f91457deb2497729p-4, erfc(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0xf.fffffffffffb7c8a2bef5924bb98p-4, erfc(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0xf.fffffffffffffedf228afbd6493p-4, erfc(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffedf2p-4, erfc(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0x4.0004157f2239d721e27728e0acacp-128, erfc(@as(f128, 0x9.31cdfp+0)));
    try std.testing.expectEqual(0x3.ffff75b4a7f71721b89fe0646f56p-128, erfc(@as(f128, 0x9.31cep+0)));
    try std.testing.expectEqual(0x3.fffd098f7c63a42c4181f6fca376p-1024, erfc(@as(f128, 0x1.a8b13p+4)));
    try std.testing.expectEqual(0x4.001799b7b63bbeff7d28fedc0018p-1024, erfc(@as(f128, 0x1.a8b12ep+4)));
    try std.testing.expectEqual(0x4.0000000000cc9079f71a2f32ab1p-1024, erfc(@as(f128, 0x1.a8b12fc6e4891p+4)));
    try std.testing.expectEqual(0x3.fffd098f7c63a42c4181f6fca376p-1024, erfc(@as(f128, 0x1.a8b13p+4)));
    try std.testing.expectEqual(0x4.001799b7b63bbeff7d28fedc0018p-1024, erfc(@as(f128, 0x1.a8b12ep+4)));
    try std.testing.expectEqual(0x3.fffffffffff8115bf0b754d4fe1p-1024, erfc(@as(f128, 0x1.a8b12fc6e4892p+4)));
    try std.testing.expectEqual(0x7.ffe048893995703e6ead0de50f7cp-972, erfc(@as(f128, 0x1.9d7adcp+4)));
    // try std.testing.expectEqual(0x8.001401a2efa2624d0e762da13978p-972, erfc(@as(f128, 0x1.9d7adap+4)));
    // try std.testing.expectEqual(0x7.ffffffffff3b1b6aef1fdb453a14p-972, erfc(@as(f128, 0x1.9d7adac608e86p+4)));
    try std.testing.expectEqual(0x8.0000000000d8e567447df7350c7p-972, erfc(@as(f128, 0x1.9d7adac608e85p+4)));
    try std.testing.expectEqual(0x7.ffffffffffffe63ab0b952a3924cp-972, erfc(@as(f128, 0x1.9d7adac608e85864p+4)));
    try std.testing.expectEqual(0x8.00000000000019f3f043fe66cfa8p-972, erfc(@as(f128, 0x1.9d7adac608e85862p+4)));
    try std.testing.expectEqual(0x8.0000000000000000000000023cdp-972, erfc(@as(f128, 0x1.9d7adac608e8586300e6c8b99ep+4)));
    try std.testing.expectEqual(0x7.ffe048893995703e6ead0de50f7cp-972, erfc(@as(f128, 0x1.9d7adcp+4)));
    // try std.testing.expectEqual(0x8.001401a2efa2624d0e762da13978p-972, erfc(@as(f128, 0x1.9d7adap+4)));
    // try std.testing.expectEqual(0x7.ffffffffff3b1b6aef1fdb453a14p-972, erfc(@as(f128, 0x1.9d7adac608e86p+4)));
    try std.testing.expectEqual(0x8.0000000000d8e567447df7350c7p-972, erfc(@as(f128, 0x1.9d7adac608e85p+4)));
    try std.testing.expectEqual(0x7.ffffffffffffe63ab0b952a3924cp-972, erfc(@as(f128, 0x1.9d7adac608e85864p+4)));
    try std.testing.expectEqual(0x8.00000000000019f3f043fe66cfa8p-972, erfc(@as(f128, 0x1.9d7adac608e85862p+4)));
    // try std.testing.expectEqual(0x7.fffffffffffffffffffffff54e7cp-972, erfc(@as(f128, 0x1.9d7adac608e8586300e6c8b99e8p+4)));
    // try std.testing.expectEqual(0x3.feff49e314f6ecdbbfebc2b34f24p-16384, erfc(@as(f128, 0x6.a89308p+4)));
    // try std.testing.expectEqual(0x4.00a9613ff5224411b6349cce295cp-16384, erfc(@as(f128, 0x6.a893p+4)));
    // try std.testing.expectEqual(0x3.fffffffffad99d5daf268a859d98p-16384, erfc(@as(f128, 0x6.a893032db9054p+4)));
    // try std.testing.expectEqual(0x4.00000000082ae9d5a43888b96c74p-16384, erfc(@as(f128, 0x6.a893032db905p+4)));
    try std.testing.expectEqual(0x4.0000000000000df012e73ddae2b4p-16384, erfc(@as(f128, 0x6.a893032db905274p+4)));
    // try std.testing.expectEqual(0x3.feff49e314f6ecdbbfebc2b34f24p-16384, erfc(@as(f128, 0x6.a89308p+4)));
    // try std.testing.expectEqual(0x4.00a9613ff5224411b6349cce295cp-16384, erfc(@as(f128, 0x6.a893p+4)));
    // try std.testing.expectEqual(0x3.fffffffffad99d5daf268a859d98p-16384, erfc(@as(f128, 0x6.a893032db9054p+4)));
    // try std.testing.expectEqual(0x4.00000000082ae9d5a43888b96c74p-16384, erfc(@as(f128, 0x6.a893032db905p+4)));
    // try std.testing.expectEqual(0x3.fffffffffffe63c683e89c3c2e14p-16384, erfc(@as(f128, 0x6.a893032db9052748p+4)));
    // try std.testing.expectEqual(0x1.ffcdcfd4f9515ad5d9921562ca2p-16384, erfc(@as(f128, 0x6.a8a058p+4)));
    // try std.testing.expectEqual(0x2.00a2fdbcb5dc489cd9c18d5c3254p-16384, erfc(@as(f128, 0x6.a8a05p+4)));
    // try std.testing.expectEqual(0x1.fffffffffb6f714cead9a1d65bc4p-16384, erfc(@as(f128, 0x6.a8a0561d8bbecp+4)));
    // try std.testing.expectEqual(0x2.00000000021824dbaeba00661e58p-16384, erfc(@as(f128, 0x6.a8a0561d8bbe8p+4)));
    // try std.testing.expectEqual(0x2.00000000000018654a1f8eeb1fb8p-16384, erfc(@as(f128, 0x6.a8a0561d8bbe942p+4)));
    // try std.testing.expectEqual(0x1.ffcdcfd4f9515ad5d9921562ca2p-16384, erfc(@as(f128, 0x6.a8a058p+4)));
    // try std.testing.expectEqual(0x2.00a2fdbcb5dc489cd9c18d5c3254p-16384, erfc(@as(f128, 0x6.a8a05p+4)));
    // try std.testing.expectEqual(0x1.fffffffffb6f714cead9a1d65bc4p-16384, erfc(@as(f128, 0x6.a8a0561d8bbecp+4)));
    // try std.testing.expectEqual(0x2.00000000021824dbaeba00661e58p-16384, erfc(@as(f128, 0x6.a8a0561d8bbe8p+4)));
    // try std.testing.expectEqual(0x1.ffffffffffff434ed847125bd78cp-16384, erfc(@as(f128, 0x6.a8a0561d8bbe9428p+4)));
    // try std.testing.expectEqual(0x3.feff49e314f6ecdbbfebc2b34f24p-16384, erfc(@as(f128, 0x6.a89308p+4)));
    // try std.testing.expectEqual(0x4.00a9613ff5224411b6349cce295cp-16384, erfc(@as(f128, 0x6.a893p+4)));
    // try std.testing.expectEqual(0x3.fffffffffad99d5daf268a859d98p-16384, erfc(@as(f128, 0x6.a893032db9054p+4)));
    // try std.testing.expectEqual(0x4.00000000082ae9d5a43888b96c74p-16384, erfc(@as(f128, 0x6.a893032db905p+4)));
    // try std.testing.expectEqual(0x3.fffffffffffe63c683e89c3c2e14p-16384, erfc(@as(f128, 0x6.a893032db9052748p+4)));
    try std.testing.expectEqual(0x4.0000000000000df012e73ddae2b4p-16384, erfc(@as(f128, 0x6.a893032db905274p+4)));
    // try std.testing.expectEqual(0x4.0000000000000000000000006d58p-16384, erfc(@as(f128, 0x6.a893032db905274042fb05c665dcp+4)));
    // try std.testing.expectEqual(0x3.fffffffffffffffffffffff8ef9cp-16384, erfc(@as(f128, 0x6.a893032db905274042fb05c666p+4)));
    try std.testing.expectEqual(0x4.0000000000000000000000637ap-16384, erfc(@as(f128, 0x6.a893032db905274042fb05c664p+4)));
    // try std.testing.expectEqual(0x3.feff49e314f6ecdbbfebc2b34f24p-16384, erfc(@as(f128, 0x6.a89308p+4)));
    // try std.testing.expectEqual(0x4.00a9613ff5224411b6349cce295cp-16384, erfc(@as(f128, 0x6.a893p+4)));
    // try std.testing.expectEqual(0x3.fffffffffad99d5daf268a859d98p-16384, erfc(@as(f128, 0x6.a893032db9054p+4)));
    // try std.testing.expectEqual(0x4.00000000082ae9d5a43888b96c74p-16384, erfc(@as(f128, 0x6.a893032db905p+4)));
    // try std.testing.expectEqual(0x3.fffffffffffe63c683e89c3c2e14p-16384, erfc(@as(f128, 0x6.a893032db9052748p+4)));
    try std.testing.expectEqual(0x4.0000000000000df012e73ddae2b4p-16384, erfc(@as(f128, 0x6.a893032db905274p+4)));
    // try std.testing.expectEqual(0x3.ffffffffffffffffffffffff9844p-16384, erfc(@as(f128, 0x6.a893032db905274042fb05c665ep+4)));
    // try std.testing.expectEqual(0x3.fffffffffffffffffffffff8ef9cp-16384, erfc(@as(f128, 0x6.a893032db905274042fb05c666p+4)));
    try std.testing.expectEqual(0x4.0000000000000000000000637ap-16384, erfc(@as(f128, 0x6.a893032db905274042fb05c664p+4)));
    try std.testing.expectEqual(0x3.fff91a7d782b006458655c2be87cp-4, erfc(@as(f128, 0xd.03d06p-4)));
    try std.testing.expectEqual(0xd.cc22642cb5ab8dc55a1975e03bcp-8, erfc(@as(f128, 0x1.5cf218p+0)));
    try std.testing.expectEqual(0xd.cc22be4b9b325bc5efb8f07224p-8, erfc(@as(f128, 0x1.5cf216p+0)));
    // try std.testing.expectEqual(0xd.cc22a7f1317ede14f24f95f12c7p-8, erfc(@as(f128, 0x1.5cf2167efe921p+0)));
    // try std.testing.expectEqual(0xd.cc22a7f13181af0c2264782d18d8p-8, erfc(@as(f128, 0x1.5cf2167efe92p+0)));
    // try std.testing.expectEqual(0xd.cc22a7f131804ea9510443197f9p-8, erfc(@as(f128, 0x1.5cf2167efe9207d2p+0)));
    // try std.testing.expectEqual(0xf.f53d075aa92b1f075d1f393d667p-8, erfc(@as(f128, 0x1.5166e2p+0)));
    // try std.testing.expectEqual(0xf.f53d6d0e58d08f84e4986455328p-8, erfc(@as(f128, 0x1.5166ep+0)));
    try std.testing.expectEqual(0xf.f53d3d6dfa74176be2273d9c6948p-8, erfc(@as(f128, 0x1.5166e0efc44aap+0)));
    // try std.testing.expectEqual(0xf.f53d3d6dfa7745095fdac02562dp-8, erfc(@as(f128, 0x1.5166e0efc44a9p+0)));
    // try std.testing.expectEqual(0xf.f53d3d6dfa747d85458d645dd238p-8, erfc(@as(f128, 0x1.5166e0efc44a9dfep+0)));
    // try std.testing.expectEqual(0xf.f53d3d6dfa747deaf93d1ace2358p-8, erfc(@as(f128, 0x1.5166e0efc44a9dfcp+0)));
    // try std.testing.expectEqual(0xf.f53d3d6dfa747dd2cb91abfddff8p-8, erfc(@as(f128, 0x1.5166e0efc44a9dfc79b8c8873a99p+0)));
    // try std.testing.expectEqual(0xf.f53d3d6dfa747dd2cb91abfdcb8p-8, erfc(@as(f128, 0x1.5166e0efc44a9dfc79b8c8873bp+0)));
    // try std.testing.expectEqual(0xf.f53d3d6dfa747dd2cb91abfde4fp-8, erfc(@as(f128, 0x1.5166e0efc44a9dfc79b8c8873a8p+0)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1p+0, erfc(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, erfc(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2p+0, erfc(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    // try std.testing.expectEqual(0x7.8cde235791e7d0dfc843bc26d51p-8, erfc(@as(f128, 0x1.8a0c64p+0)));
    // try std.testing.expectEqual(0x7.8cde5963040180b50eb9ef7f8a04p-8, erfc(@as(f128, 0x1.8a0c62p+0)));
    // try std.testing.expectEqual(0xc.766cbf61fd6480afed02bf2b2068p-8, erfc(@as(f128, 0x1.64dafap+0)));
    try std.testing.expectEqual(0x7.23ff79ae0f25a135a4973efb8be4p-68, erfc(@as(f128, 0x6.88fb08p+0)));
    try std.testing.expectEqual(0x3.e2fa6064d589347b0f2f7aa6e388p-4, erfc(@as(f128, 0xd.361d9p-4)));
    // try std.testing.expectEqual(0x1.eb9635bc51eb7a94581f979ead1dp-116, erfc(@as(f128, 0x8.c66b5p+0)));
    // try std.testing.expectEqual(0x1.eb98546946cb2525a4905a3b1382p-116, erfc(@as(f128, 0x8.c66b4p+0)));
    // try std.testing.expectEqual(0x1.eb97b1f20867c35eff191bbeca3ep-116, erfc(@as(f128, 0x8.c66b44ca40038p+0)));
    // try std.testing.expectEqual(0x3.ba3ac339ed1902051ea00716755p-12, erfc(@as(f128, 0x2.586f1cp+0)));
    try std.testing.expectEqual(0x7.ee2d2ec57315040047a2a1252e8p-204, erfc(@as(f128, 0xb.acb72p+0)));
    try std.testing.expectEqual(0x1.c646841c902106e7ce3048dea085p-184, erfc(@as(f128, 0xb.2274ap+0)));
    try std.testing.expectEqual(0x1.c648feeb672e8e57298792150dd9p-184, erfc(@as(f128, 0xb.22749p+0)));
    // try std.testing.expectEqual(0x1.c6479753ddcb4d8c72ebbfb8ec29p-184, erfc(@as(f128, 0xb.227499103358p+0)));
    try std.testing.expectEqual(0x1.c6479753dddf2401559c4dbabe7p-184, erfc(@as(f128, 0xb.2274991033578p+0)));
    try std.testing.expectEqual(0x1.c6479753ddd176a5bf5193bad771p-184, erfc(@as(f128, 0xb.227499103357d84p+0)));
    // try std.testing.expectEqual(0x3.eaab96d5a2e294b81fff40fde9aap-4, erfc(@as(f128, 0xd.28abfp-4)));
    // try std.testing.expectEqual(0xf.bbc04428a3d30e77d2315d0046a8p-8, erfc(@as(f128, 0x1.5289fep+0)));
    // try std.testing.expectEqual(0x1.f57fab6c3db3ce7e0bd2fb137939p-36, erfc(@as(f128, 0x4.b48498p+0)));
    // try std.testing.expectEqual(0x1.be98de114e174b501b7acff72e8p-16, erfc(@as(f128, 0x2.f8646cp+0)));
    // try std.testing.expectEqual(0xf.fbeadad5a51f774a6aa2da69dad8p-8, erfc(@as(f128, 0x1.514548p+0)));
    try std.testing.expectEqual(0x7.22d059993f3f46d0e0daa16357p-12, erfc(@as(f128, 0x2.36c504p+0)));
    // try std.testing.expectEqual(0xc.4bf9de451e5fced9d5e2d18c20b8p-8, erfc(@as(f128, 0x1.65e31p+0)));
    try std.testing.expectEqual(0x3.da9f608f1dd7ee3168650dc2fb9ep-4, erfc(@as(f128, 0xd.44cd3p-4)));
    try std.testing.expectEqual(0x3.d93aa59c8f5abb821749e8017ae2p-4, erfc(@as(f128, 0xd.47426p-4)));
    try std.testing.expectEqual(0x3.d93aaeadb64d00e8ad67712ba71p-4, erfc(@as(f128, 0xd.47425p-4)));
    try std.testing.expectEqual(0x3.d93aa84f87a9ffa04577ca7dbb28p-4, erfc(@as(f128, 0xd.47425b3cafa48p-4)));
    // try std.testing.expectEqual(0x1.7fefc09137c9485d5871f07f9465p-4, erfc(@as(f128, 0x1.2f644ep+0)));
    try std.testing.expectEqual(0x3.dbca059c7e73a1239dd52028280cp-12, erfc(@as(f128, 0x2.56af04p+0)));
    try std.testing.expectEqual(0x7.e8b2efb679451a42955c7a94bfbcp-16, erfc(@as(f128, 0x2.b7f8ccp+0)));
    try std.testing.expectEqual(0x7.e8b3a6276f03f7798c2a7c4e6284p-16, erfc(@as(f128, 0x2.b7f8c8p+0)));
    try std.testing.expectEqual(0x7.e8b308381dfc55c4841397b30ae4p-16, erfc(@as(f128, 0x2.b7f8cb76737d4p+0)));
    // try std.testing.expectEqual(0x7.e8b308381e02094c04b2b4fe1ddp-16, erfc(@as(f128, 0x2.b7f8cb76737d2p+0)));
    // try std.testing.expectEqual(0x7.e8b308381e00144be16c16cf88ep-16, erfc(@as(f128, 0x2.b7f8cb76737d2afcp+0)));
    // try std.testing.expectEqual(0x7.e8b308381e001502525c2ab33258p-16, erfc(@as(f128, 0x2.b7f8cb76737d2af8p+0)));
    try std.testing.expectEqual(0x7.e8b308381e0014bb6d3bd6db599cp-16, erfc(@as(f128, 0x2.b7f8cb76737d2af98dead7c4c5eep+0)));
    try std.testing.expectEqual(0x7.e8b308381e0014bb6d3bd6db5664p-16, erfc(@as(f128, 0x2.b7f8cb76737d2af98dead7c4c6p+0)));
    try std.testing.expectEqual(0x7.e8b308381e0014bb6d3bd6db8404p-16, erfc(@as(f128, 0x2.b7f8cb76737d2af98dead7c4c5p+0)));
    try std.testing.expectEqual(0x3.281c2d7e470e5082e4209788692ap-16, erfc(@as(f128, 0x2.dfb9b4p+0)));
    // try std.testing.expectEqual(0x1.f1cb04b622e6f4d5035449633b46p-8, erfc(@as(f128, 0x1.e33c9ep+0)));
    // try std.testing.expectEqual(0x1.3bd95ffe4e5561c5991cb64b6573p-4, erfc(@as(f128, 0x1.3ffccp+0)));
    // try std.testing.expectEqual(0x1.3bd9679020a687bf0ac713ffaf7bp-4, erfc(@as(f128, 0x1.3ffcbep+0)));
    // try std.testing.expectEqual(0x1.3bd962ebb773644beafd5d55b35fp-4, erfc(@as(f128, 0x1.3ffcbf39febb4p+0)));
}
