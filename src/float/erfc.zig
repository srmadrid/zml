const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
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
    const axf: f32 = float.abs(xf);
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
        const s: f64 = float.abs(x) - 1;
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
        const xx: f64 = float.abs(x);
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
        const r: f64 = float.exp(-z * z - 0.5625) * float.exp((z - xx) * (z + xx) + R / S);
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

        return 1 - float.erf(x);
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

        const xx: f128 = float.abs(x);
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
        const r: f128 = float.exp(-z * z - 0.5625) * float.exp((z - xx) * (z + xx) + p);

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
