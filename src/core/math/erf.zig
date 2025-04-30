const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const erf_data = @import("erf_data.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn erf(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return erf(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, erf32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_erff.c
                    return erf32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_erf.c
                    return erf64(x);
                },
                f80 => return cast(f80, erf128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_erfl.c
                    return erf128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn erf32(x: f32) f32 {
    const ax: f32 = math.abs(x);
    const ux: u32 = @bitCast(ax);
    const s: f64 = cast(f64, x, .{});
    var z: f64 = cast(f64, ax, .{});
    // 0x407ad444 corresponds to x = 0x1.f5a888p+1 = 3.91921..., which is the
    // largest float such that erf(x) does not round to 1 (to nearest).
    if (ux > 0x407ad444) {
        @branchHint(.unlikely);
        const os: f32 = math.copysign(@as(f32, 1), x);

        if (ux > (0xff << 23))
            return x + x; // nan

        if (ux == (0xff << 23))
            return os; // +-inf

        return os - 0x1p-25 * os;
    }

    const v: f64 = math.floor(16 * z);
    const i: u32 = cast(u32, 16 * ax, .{});
    // 0x3ee00000 corresponds to x = 0.4375, for smaller x we have i < 7.
    if (ux < 0x3ee00000) {
        @branchHint(.unlikely);
        const c: [8]f64 = .{
            0x1.20dd750429b6dp+0,  -0x1.812746b0375fbp-2,
            0x1.ce2f219fd6f45p-4,  -0x1.b82ce2cbf0838p-6,
            0x1.565bb655adb85p-8,  -0x1.c025bfc879c94p-11,
            0x1.f81718f61309cp-14, -0x1.cc67bd88f5867p-17,
        };
        const z2: f64 = s * s;
        const z4: f64 = z2 * z2;
        const z8: f64 = z4 * z4;
        var c0: f64 = c[0] + z2 * c[1];
        const c2: f64 = c[2] + z2 * c[3];
        var c4: f64 = c[4] + z2 * c[5];
        const c6: f64 = c[6] + z2 * c[7];
        c0 += z4 * c2;
        c4 += z4 * c6;
        c0 += z8 * c4;
        return cast(f32, s * c0, .{});
    }

    z = (z - 0.03125) - 0.0625 * v;
    const c: [8]f64 = erf_data.C_32[i - 7];
    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    var c0: f64 = c[0] + z * c[1];
    const c2: f64 = c[2] + z * c[3];
    var c4: f64 = c[4] + z * c[5];
    const c6: f64 = c[6] + z * c[7];
    c0 += z2 * c2;
    c4 += z2 * c6;
    c0 += z4 * c4;
    return cast(f32, math.copysign(c0, s), .{});
}

fn erf64(x: f64) f64 {
    var hx: i32 = undefined;
    dbl64.getHighWord(&hx, x);
    const ix: i32 = hx & 0x7fffffff;
    if (ix >= 0x7ff00000) { // erf(nan)=nan
        @setRuntimeSafety(false);
        const i: i32 = cast(i32, (@as(u32, @intCast(hx)) >> 31) << 1, .{});
        return cast(f64, (1 - i) + 1, .{}) / x; // erf(+-inf)=+-1
    }

    if (ix < 0x3feb0000) { // |x|<0.84375
        if (ix < 0x3e300000) { // |x|<2**-28
            if (ix < 0x00800000) {
                // Avoid spurious underflow.
                const ret: f64 = 0.0625 * (16 * x + (16 * erf_data.efx_64) * x);

                if (math.abs(ret) < std.math.floatMin(f64)) {
                    const vret: f64 = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }

                return ret;
            }

            return x + erf_data.efx_64 * x;
        }

        const z: f64 = x * x;
        const r1: f64 = erf_data.pp_64[0] + z * erf_data.pp_64[1];
        const z2: f64 = z * z;
        const r2: f64 = erf_data.pp_64[2] + z * erf_data.pp_64[3];
        const z4: f64 = z2 * z2;
        const s1: f64 = 1 + z * erf_data.qq_64[1];
        const s2: f64 = erf_data.qq_64[2] + z * erf_data.qq_64[3];
        const s3: f64 = erf_data.qq_64[4] + z * erf_data.qq_64[5];
        const r: f64 = r1 + z2 * r2 + z4 * erf_data.pp_64[4];
        const s: f64 = s1 + z2 * s2 + z4 * s3;
        const y: f64 = r / s;
        return x + x * y;
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
        const P = P1 + s2 * P2 + s4 * P3 + s6 * P4;
        const Q = Q1 + s2 * Q2 + s4 * Q3 + s6 * Q4;
        if (hx >= 0) {
            return erf_data.erx_64 + P / Q;
        } else {
            return -erf_data.erx_64 - P / Q;
        }
    }

    if (ix >= 0x40180000) { // inf>|x|>=6
        if (hx >= 0) {
            return 1 - erf_data.tiny_64;
        } else {
            return erf_data.tiny_64 - 1;
        }
    }

    const xx: f64 = math.abs(x);
    const s: f64 = 1 / (xx * xx);
    var R: f64 = undefined;
    var S: f64 = undefined;
    if (ix < 0x4006DB6E) { // |x| < 1/0.35
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
    } else { // |x| >= 1/0.35
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

    if (hx >= 0) {
        return 1 - r / xx;
    } else {
        return r / xx - 1;
    }
}

fn erf128(x: f128) f128 {
    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const sign: i32 = @bitCast(u.w0);
    const ix: i32 = sign & 0x7fffffff;
    if (ix >= 0x7fff0000) { // erf(nan)=nan
        const i: i32 = @bitCast(((@as(u32, @bitCast(sign)) & 0xffff0000) >> 31) << 1);
        return cast(f128, (1 - i) + 1, .{}) / x; // erf(+-inf)=+-1
    }

    if (ix >= 0x3fff0000) { // |x| >= 1.0
        if (ix >= 0x40030000 and sign > 0)
            return 1; // x >= 16, avoid spurious underflow from erfc.

        const y: f128 = math.erfc(x);
        return 1 - y;
        // return (1 - __erfcl (x));
    }

    u.w0 = @bitCast(ix);
    var a: f128 = @bitCast(u);
    const z: f128 = x * x;
    var y: f128 = undefined;
    if (ix < 0x3ffec000) { // a < 0.875
        if (ix < 0x3fc60000) { // |x|<2**-57
            if (ix < 0x00080000) {
                // Avoid spurious underflow.
                const ret: f128 = 0.0625 * (16.0 * x + (16.0 * erf_data.efx_128) * x);

                if (math.abs(ret) < std.math.floatMin(f128)) {
                    const vret: f128 = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }

                return ret;
            }
            return x + erf_data.efx_128 * x;
        }

        y = a + a * erf_data.neval(z, &erf_data.TN1_128, 8) / erf_data.deval(z, &erf_data.TD1_128, 8);
    } else {
        a = a - 1;
        y = erf_data.erf_const_128 + erf_data.neval(a, &erf_data.TN2_128, 8) / erf_data.deval(a, &erf_data.TD2_128, 8);
    }

    if ((@as(u32, @bitCast(sign)) & 0x80000000) != 0) // x < 0
        y = -y;

    return y;
}

test erf {
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, erf(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x2.3ebc34p-4, erf(@as(f32, 0x2p-4)));
    try std.testing.expectEqual(0xb.60e4cp-4, erf(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0xd.7bb3dp-4, erf(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0xd.7bb3dp-4, erf(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0xe.c432fp-4, erf(@as(f32, 0x1.4p+0)));
    try std.testing.expectEqual(0xf.ecd71p-4, erf(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(-0xf.ecd71p-4, erf(@as(f32, -0x2p+0)));
    try std.testing.expectEqual(0xf.ffe8dp-4, erf(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(-0xf.ffe8dp-4, erf(@as(f32, -0x3p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f32, -0x4p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x4.2p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x5p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x6p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x7p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x8p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x9p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0xap+0)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f32, -0xap+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x1.bp+4)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f32, -0x1.bp+4)));
    try std.testing.expectEqual(-0x8.53f7ap-4, erf(@as(f32, -0x7.fffff8p-4)));
    try std.testing.expectEqual(-0x8.53f7bp-4, erf(@as(f32, -0x8p-4)));
    try std.testing.expectEqual(0x4.000018p-128, erf(@as(f32, 0x3.8b7f28p-128)));
    try std.testing.expectEqual(0x4.0000a8p-128, erf(@as(f32, 0x3.8b7fa8p-128)));
    try std.testing.expectEqual(-0x0p+0, erf(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, erf(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x3.ffffdp-128, erf(@as(f32, 0x3.8b7ee8p-128)));
    try std.testing.expectEqual(0x4.00003p-128, erf(@as(f32, 0x3.8b7f4p-128)));
    try std.testing.expectEqual(0x4.000028p-128, erf(@as(f32, 0x3.8b7f38p-128)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x1.ap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x1.cp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x6.4p+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x6.ap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x6.a8p+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x6.aap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x6.bp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x6.cp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0x3.e8p+8)));
    try std.testing.expectEqual(0x9.062b2p-8, erf(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(-0x9.062b2p-8, erf(@as(f32, -0x8p-8)));
    try std.testing.expectEqual(0x4.8375b8p-12, erf(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x2.41baecp-16, erf(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0x1.20dd76p-20, erf(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0x9.06ebbp-28, erf(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x4.8375d8p-32, erf(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x2.41baecp-36, erf(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1.20dd76p-40, erf(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x9.06ebbp-48, erf(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x4.8375d8p-52, erf(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x2.41baecp-56, erf(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1.20dd76p-60, erf(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1.20dd76p-100, erf(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d8p-128, erf(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4.8375d8p-128, erf(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, erf(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, erf(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.dded1p-4, erf(@as(f32, -0x1.ddaea4p+0)));
    try std.testing.expectEqual(-0xe.6cc4p-4, erf(@as(f32, -0x1.2b1f68p+0)));
    try std.testing.expectEqual(0xe.d6505p-4, erf(@as(f32, 0x1.44e722p+0)));
    try std.testing.expectEqual(-0xe.ad06ep-4, erf(@as(f32, -0x1.3a0d48p+0)));
    try std.testing.expectEqual(-0xf.d0e5ap-4, erf(@as(f32, -0x1.c975cap+0)));
    try std.testing.expectEqual(-0xf.e2945p-4, erf(@as(f32, -0x1.e6a006p+0)));
    try std.testing.expectEqual(-0x1.77f98ep-12, erf(@as(f32, -0x1.4d32f4p-12)));

    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, erf(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x2.3ebc346b87712p-4, erf(@as(f64, 0x2p-4)));
    try std.testing.expectEqual(0xb.60e4bace872f8p-4, erf(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0xd.7bb3d3a084458p-4, erf(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0xd.7bb3d3a084458p-4, erf(@as(f64, -0x1p+0)));
    // try std.testing.expectEqual(0xe.c432ecc55f008p-4, erf(@as(f64, 0x1.4p+0)));
    try std.testing.expectEqual(0xf.ecd70a13caf18p-4, erf(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(-0xf.ecd70a13caf18p-4, erf(@as(f64, -0x2p+0)));
    try std.testing.expectEqual(0xf.ffe8d6209afc8p-4, erf(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(-0xf.ffe8d6209afc8p-4, erf(@as(f64, -0x3p+0)));
    try std.testing.expectEqual(0xf.fffffbdc88bbp-4, erf(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(-0xf.fffffbdc88bbp-4, erf(@as(f64, -0x4p+0)));
    try std.testing.expectEqual(0xf.fffffe8b4e86p-4, erf(@as(f64, 0x4.2p+0)));
    try std.testing.expectEqual(0xf.ffffffffe4f4p-4, erf(@as(f64, 0x5p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x6p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x7p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x8p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x9p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0xap+0)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f64, -0xap+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x1.bp+4)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f64, -0x1.bp+4)));
    try std.testing.expectEqual(-0x8.53f7a704b7bep-4, erf(@as(f64, -0x7.fffff8p-4)));
    try std.testing.expectEqual(-0x8.53f7ae0c76e9p-4, erf(@as(f64, -0x8p-4)));
    try std.testing.expectEqual(-0x8.53f7ae0c76e9p-4, erf(@as(f64, -0x7.ffffffffffffcp-4)));
    try std.testing.expectEqual(0x4.000018956724p-128, erf(@as(f64, 0x3.8b7f28p-128)));
    try std.testing.expectEqual(0x4.0000a90421a64p-128, erf(@as(f64, 0x3.8b7fa8p-128)));
    try std.testing.expectEqual(-0x0p+0, erf(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x9.06eba8214db68p-152, erf(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4.0000000000004p-1024, erf(@as(f64, -0x3.8b7f12369ded8p-1024)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x3.ffffd05e09e3p-128, erf(@as(f64, 0x3.8b7ee8p-128)));
    try std.testing.expectEqual(0x4.000033aa2a1c8p-128, erf(@as(f64, 0x3.8b7f4p-128)));
    try std.testing.expectEqual(0x4.00002aa33e744p-128, erf(@as(f64, 0x3.8b7f38p-128)));
    try std.testing.expectEqual(0x4.00002f26b4488p-128, erf(@as(f64, 0x3.8b7f3cp-128)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x3.fffffffffffe8p-1024, erf(@as(f64, 0x3.8b7f12369decp-1024)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4.0000000000018p-1024, erf(@as(f64, 0x3.8b7f12369deecp-1024)));
    try std.testing.expectEqual(0x4.0000000000014p-1024, erf(@as(f64, 0x3.8b7f12369dee8p-1024)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x8p-972, erf(@as(f64, 0x7.16fe246d3bdacp-972)));
    try std.testing.expectEqual(0x7.ffffffffffffcp-972, erf(@as(f64, 0x7.16fe246d3bda8p-972)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x8p-972, erf(@as(f64, 0x7.16fe246d3bdacp-972)));
    try std.testing.expectEqual(0x7.ffffffffffffcp-972, erf(@as(f64, 0x7.16fe246d3bda8p-972)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x1.ap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x1.cp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x6.4p+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x6.ap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x6.a8p+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x6.aap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x6.bp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x6.cp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0x3.e8p+8)));
    try std.testing.expectEqual(0x9.062b22ee929cp-8, erf(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(-0x9.062b22ee929cp-8, erf(@as(f64, -0x8p-8)));
    try std.testing.expectEqual(0x4.8375bbfe32e3cp-12, erf(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0x2.41baea05511f2p-16, erf(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0x1.20dd750429568p-20, erf(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0x9.06eba8214db6p-28, erf(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x4.8375d410a6db4p-32, erf(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x2.41baea08536dap-36, erf(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0x1.20dd750429b6dp-40, erf(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x9.06eba8214db68p-48, erf(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0x4.8375d410a6db4p-52, erf(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x2.41baea08536dap-56, erf(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1.20dd750429b6dp-60, erf(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1.20dd750429b6dp-100, erf(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.20dd750429b6dp-600, erf(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x4.8375d410a6db4p-128, erf(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4.8375d410a6db4p-1024, erf(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x9.06eba8214db68p-972, erf(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4.8375d410a6db4p-128, erf(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4.8375d410a6db4p-1024, erf(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x9.06eba8214db68p-972, erf(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x9.06eba8214db68p-152, erf(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, erf(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x9.06eba8214db68p-152, erf(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, erf(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xf.dded081f58dc8p-4, erf(@as(f64, -0x1.ddaea4p+0)));
    try std.testing.expectEqual(-0xe.6cc3fab61fed8p-4, erf(@as(f64, -0x1.2b1f68p+0)));
    try std.testing.expectEqual(0xe.d6504b655135p-4, erf(@as(f64, 0x1.44e722p+0)));
    try std.testing.expectEqual(-0xe.ad06dfdab8f4p-4, erf(@as(f64, -0x1.3a0d48p+0)));
    try std.testing.expectEqual(-0xf.d0e59ffd2b02p-4, erf(@as(f64, -0x1.c975cap+0)));
    try std.testing.expectEqual(-0xf.e294502e0c118p-4, erf(@as(f64, -0x1.e6a006p+0)));
    try std.testing.expectEqual(-0x1.77f98ef609eb3p-12, erf(@as(f64, -0x1.4d32f4p-12)));

    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, erf(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x2.3ebc346b87712e84p-4, erf(@as(f80, 0x2p-4)));
    try std.testing.expectEqual(0xb.60e4bace872fb63p-4, erf(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0xd.7bb3d3a08445637p-4, erf(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0xd.7bb3d3a08445637p-4, erf(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0xe.c432ecc55f00406p-4, erf(@as(f80, 0x1.4p+0)));
    try std.testing.expectEqual(0xf.ecd70a13caf1997p-4, erf(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(-0xf.ecd70a13caf1997p-4, erf(@as(f80, -0x2p+0)));
    try std.testing.expectEqual(0xf.ffe8d6209afcbddp-4, erf(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(-0xf.ffe8d6209afcbddp-4, erf(@as(f80, -0x3p+0)));
    try std.testing.expectEqual(0xf.fffffbdc88bb10bp-4, erf(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(-0xf.fffffbdc88bb10bp-4, erf(@as(f80, -0x4p+0)));
    try std.testing.expectEqual(0xf.fffffe8b4e862e1p-4, erf(@as(f80, 0x4.2p+0)));
    try std.testing.expectEqual(0xf.ffffffffe4f3e59p-4, erf(@as(f80, 0x5p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffe73p-4, erf(@as(f80, 0x6p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x7p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x8p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x9p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0xap+0)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f80, -0xap+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x1.bp+4)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f80, -0x1.bp+4)));
    try std.testing.expectEqual(-0x8.53f7a704b7be2d6p-4, erf(@as(f80, -0x7.fffff8p-4)));
    try std.testing.expectEqual(-0x8.53f7ae0c76e915fp-4, erf(@as(f80, -0x8p-4)));
    try std.testing.expectEqual(-0x8.53f7ae0c76e8ddbp-4, erf(@as(f80, -0x7.ffffffffffffcp-4)));
    try std.testing.expectEqual(-0x8.53f7ae0c76e8f9dp-4, erf(@as(f80, -0x7.ffffffffffffep-4)));
    try std.testing.expectEqual(0x4.0000189567240f1p-128, erf(@as(f80, 0x3.8b7f28p-128)));
    try std.testing.expectEqual(0x4.0000a90421a623e8p-128, erf(@as(f80, 0x3.8b7fa8p-128)));
    try std.testing.expectEqual(-0x0p+0, erf(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x9.06eba8214db688dp-152, erf(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4.000000000000309p-1024, erf(@as(f80, -0x3.8b7f12369ded8p-1024)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x4.0000000000000028p-16384, erf(@as(f80, 0x3.8b7f12369ded5518p-16384)));
    try std.testing.expectEqual(0x3.ffffd05e09e304a4p-128, erf(@as(f80, 0x3.8b7ee8p-128)));
    try std.testing.expectEqual(0x4.000033aa2a1c72f8p-128, erf(@as(f80, 0x3.8b7f4p-128)));
    try std.testing.expectEqual(0x4.00002aa33e7451a8p-128, erf(@as(f80, 0x3.8b7f38p-128)));
    try std.testing.expectEqual(0x4.00002f26b448625p-128, erf(@as(f80, 0x3.8b7f3cp-128)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x3.fffffffffffe7f48p-1024, erf(@as(f80, 0x3.8b7f12369decp-1024)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.00000000000199a8p-1024, erf(@as(f80, 0x3.8b7f12369deecp-1024)));
    try std.testing.expectEqual(0x4.000000000001517p-1024, erf(@as(f80, 0x3.8b7f12369dee8p-1024)));
    try std.testing.expectEqual(0x4.000000000001759p-1024, erf(@as(f80, 0x3.8b7f12369deeap-1024)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x8.00000000000018fp-972, erf(@as(f80, 0x7.16fe246d3bdacp-972)));
    try std.testing.expectEqual(0x7.ffffffffffffd0b8p-972, erf(@as(f80, 0x7.16fe246d3bda8p-972)));
    try std.testing.expectEqual(0x8p-972, erf(@as(f80, 0x7.16fe246d3bdaa9e8p-972)));
    try std.testing.expectEqual(0x7.fffffffffffffff8p-972, erf(@as(f80, 0x7.16fe246d3bdaa9ep-972)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x8.00000000000018fp-972, erf(@as(f80, 0x7.16fe246d3bdacp-972)));
    try std.testing.expectEqual(0x7.ffffffffffffd0b8p-972, erf(@as(f80, 0x7.16fe246d3bda8p-972)));
    try std.testing.expectEqual(0x8p-972, erf(@as(f80, 0x7.16fe246d3bdaa9e8p-972)));
    try std.testing.expectEqual(0x7.fffffffffffffff8p-972, erf(@as(f80, 0x7.16fe246d3bdaa9ep-972)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x3.ffffffffffffffdp-16384, erf(@as(f80, 0x3.8b7f12369ded54c8p-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x4.000000000000003p-16384, erf(@as(f80, 0x3.8b7f12369ded552p-16384)));
    try std.testing.expectEqual(0x4.0000000000000028p-16384, erf(@as(f80, 0x3.8b7f12369ded5518p-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.ffffffffffffffe8p-16384, erf(@as(f80, 0x1.c5bf891b4ef6aa68p-16384)));
    try std.testing.expectEqual(0x1.ffffffffffffffep-16384, erf(@as(f80, 0x1.c5bf891b4ef6aa6p-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x2.0000000000000018p-16384, erf(@as(f80, 0x1.c5bf891b4ef6aa9p-16384)));
    try std.testing.expectEqual(0x2.000000000000001p-16384, erf(@as(f80, 0x1.c5bf891b4ef6aa88p-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x4.0000000000000008p-16384, erf(@as(f80, 0x3.8b7f12369ded54f8p-16384)));
    try std.testing.expectEqual(0x4p-16384, erf(@as(f80, 0x3.8b7f12369ded54fp-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x4.0000000000000008p-16384, erf(@as(f80, 0x3.8b7f12369ded54f8p-16384)));
    try std.testing.expectEqual(0x4p-16384, erf(@as(f80, 0x3.8b7f12369ded54fp-16384)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x1.ap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x1.cp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x6.4p+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x6.ap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x6.a8p+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x6.aap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x6.bp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x6.cp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0x3.e8p+8)));
    try std.testing.expectEqual(0x9.062b22ee929bfcap-8, erf(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(-0x9.062b22ee929bfcap-8, erf(@as(f80, -0x8p-8)));
    try std.testing.expectEqual(0x4.8375bbfe32e3ccb8p-12, erf(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x2.41baea05511f14d8p-16, erf(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0x1.20dd75042956874ap-20, erf(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0x9.06eba8214db5c84p-28, erf(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x4.8375d410a6db445p-32, erf(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x2.41baea08536da234p-36, erf(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ap-40, erf(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-48, erf(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-52, erf(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x2.41baea08536da234p-56, erf(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ap-60, erf(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ap-100, erf(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ap-600, erf(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ap-10000, erf(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-128, erf(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1024, erf(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-16384, erf(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2.41baea08536da238p-16384, erf(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-972, erf(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4.8375d410a6db4468p-128, erf(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4.8375d410a6db4468p-1024, erf(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4.8375d410a6db4468p-16384, erf(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2.41baea08536da238p-16384, erf(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x9.06eba8214db688dp-972, erf(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x9.06eba8214db688dp-152, erf(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4.8375d410a6db4468p-1076, erf(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, erf(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x9.06eba8214db688dp-152, erf(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4.8375d410a6db4468p-1076, erf(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, erf(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0xf.dded081f58dc8a9p-4, erf(@as(f80, -0x1.ddaea4p+0)));
    try std.testing.expectEqual(-0xe.6cc3fab61fed855p-4, erf(@as(f80, -0x1.2b1f68p+0)));
    try std.testing.expectEqual(0xe.d6504b655134fdcp-4, erf(@as(f80, 0x1.44e722p+0)));
    try std.testing.expectEqual(-0xe.ad06dfdab8f3efdp-4, erf(@as(f80, -0x1.3a0d48p+0)));
    try std.testing.expectEqual(-0xf.d0e59ffd2b01ccap-4, erf(@as(f80, -0x1.c975cap+0)));
    try std.testing.expectEqual(-0xf.e294502e0c11683p-4, erf(@as(f80, -0x1.e6a006p+0)));
    try std.testing.expectEqual(-0x1.77f98ef609eb313p-12, erf(@as(f80, -0x1.4d32f4p-12)));

    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, erf(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x2.3ebc346b87712e85b6b249f079e2p-4, erf(@as(f128, 0x2p-4)));
    try std.testing.expectEqual(0xb.60e4bace872fb62865e59788aa7p-4, erf(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0xd.7bb3d3a0844563680887edd86938p-4, erf(@as(f128, 0x1p+0)));
    // try std.testing.expectEqual(-0xd.7bb3d3a0844563680887edd86938p-4, erf(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0xe.c432ecc55f00406276a08d164e3p-4, erf(@as(f128, 0x1.4p+0)));
    try std.testing.expectEqual(0xf.ecd70a13caf19972801904b9a34p-4, erf(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(-0xf.ecd70a13caf19972801904b9a34p-4, erf(@as(f128, -0x2p+0)));
    try std.testing.expectEqual(0xf.ffe8d6209afcbdd5f43d9ad9debp-4, erf(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(-0xf.ffe8d6209afcbdd5f43d9ad9debp-4, erf(@as(f128, -0x3p+0)));
    try std.testing.expectEqual(0xf.fffffbdc88bb10b2865615db403p-4, erf(@as(f128, 0x4p+0)));
    try std.testing.expectEqual(-0xf.fffffbdc88bb10b2865615db403p-4, erf(@as(f128, -0x4p+0)));
    try std.testing.expectEqual(0xf.fffffe8b4e862e1457f60d1cddd8p-4, erf(@as(f128, 0x4.2p+0)));
    try std.testing.expectEqual(0xf.ffffffffe4f3e58a6088c76ca158p-4, erf(@as(f128, 0x5p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffe7307eaa82df49e8p-4, erf(@as(f128, 0x6p+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffffcd6bafd9208p-4, erf(@as(f128, 0x7p+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffff1c58p-4, erf(@as(f128, 0x8p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x9p+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0xap+0)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f128, -0xap+0)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x1.bp+4)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f128, -0x1.bp+4)));
    try std.testing.expectEqual(-0x8.53f7a704b7be2d643b9e3ae3cbp-4, erf(@as(f128, -0x7.fffff8p-4)));
    try std.testing.expectEqual(-0x8.53f7ae0c76e915e809f1a31a27a8p-4, erf(@as(f128, -0x8p-4)));
    // try std.testing.expectEqual(-0x8.53f7ae0c76e8ddaa10a86e7a049p-4, erf(@as(f128, -0x7.ffffffffffffcp-4)));
    try std.testing.expectEqual(-0x8.53f7ae0c76e8f9c90d4d08ca1638p-4, erf(@as(f128, -0x7.ffffffffffffep-4)));
    try std.testing.expectEqual(0x4.0000189567240f10919b31d6bcf4p-128, erf(@as(f128, 0x3.8b7f28p-128)));
    try std.testing.expectEqual(0x4.0000a90421a623ebfa28a3ab4774p-128, erf(@as(f128, 0x3.8b7fa8p-128)));
    try std.testing.expectEqual(-0x0p+0, erf(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4.00000000000030934524cf4ab6ep-1024, erf(@as(f128, -0x3.8b7f12369ded8p-1024)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x4.0000000000000029274014aceae4p-16384, erf(@as(f128, 0x3.8b7f12369ded5518p-16384)));
    try std.testing.expectEqual(0x3.ffffd05e09e304a2dd5478ec77b4p-128, erf(@as(f128, 0x3.8b7ee8p-128)));
    try std.testing.expectEqual(0x4.000033aa2a1c72f9b535b72e96ecp-128, erf(@as(f128, 0x3.8b7f4p-128)));
    try std.testing.expectEqual(0x4.00002aa33e7451abfeace0114e44p-128, erf(@as(f128, 0x3.8b7f38p-128)));
    try std.testing.expectEqual(0x4.00002f26b4486252d9f14b9ff298p-128, erf(@as(f128, 0x3.8b7f3cp-128)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x3.fffffffffffe7f47159e90b87d3ap-1024, erf(@as(f128, 0x3.8b7f12369decp-1024)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.00000000000199a8176a036f3c44p-1024, erf(@as(f128, 0x3.8b7f12369deecp-1024)));
    try std.testing.expectEqual(0x4.0000000000015170ba28f90187fcp-1024, erf(@as(f128, 0x3.8b7f12369dee8p-1024)));
    try std.testing.expectEqual(0x4.000000000001758c68c97e38622p-1024, erf(@as(f128, 0x3.8b7f12369deeap-1024)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x8.00000000000018ef2d089427b98p-972, erf(@as(f128, 0x7.16fe246d3bdacp-972)));
    try std.testing.expectEqual(0x7.ffffffffffffd0b7cfc789ba0534p-972, erf(@as(f128, 0x7.16fe246d3bda8p-972)));
    try std.testing.expectEqual(0x8.00000000000000011037402e1a6p-972, erf(@as(f128, 0x7.16fe246d3bdaa9e8p-972)));
    // try std.testing.expectEqual(0x7.fffffffffffffff8094b980cccacp-972, erf(@as(f128, 0x7.16fe246d3bdaa9ep-972)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffffe86cp-972, erf(@as(f128, 0x7.16fe246d3bdaa9e70ec1483562p-972)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffffeabp-972, erf(@as(f128, 0x7.16fe246d3bdaa9e70ec1483564p-972)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffffe62cp-972, erf(@as(f128, 0x7.16fe246d3bdaa9e70ec148356p-972)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x8.00000000000018ef2d089427b98p-972, erf(@as(f128, 0x7.16fe246d3bdacp-972)));
    try std.testing.expectEqual(0x7.ffffffffffffd0b7cfc789ba0534p-972, erf(@as(f128, 0x7.16fe246d3bda8p-972)));
    try std.testing.expectEqual(0x8.00000000000000011037402e1a6p-972, erf(@as(f128, 0x7.16fe246d3bdaa9e8p-972)));
    // try std.testing.expectEqual(0x7.fffffffffffffff8094b980cccacp-972, erf(@as(f128, 0x7.16fe246d3bdaa9ep-972)));
    try std.testing.expectEqual(0x8.00000000000000000000000017dp-972, erf(@as(f128, 0x7.16fe246d3bdaa9e70ec148358cp-972)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x3.ffffffffffffffcee20b835fe1c4p-16384, erf(@as(f128, 0x3.8b7f12369ded54c8p-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x4.00000000000000322e2bbcce389cp-16384, erf(@as(f128, 0x3.8b7f12369ded552p-16384)));
    try std.testing.expectEqual(0x4.0000000000000029274014aceae4p-16384, erf(@as(f128, 0x3.8b7f12369ded5518p-16384)));
    try std.testing.expectEqual(0x4.000000000000002daab5e8bd91cp-16384, erf(@as(f128, 0x3.8b7f12369ded551cp-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    // try std.testing.expectEqual(0x1.ffffffffffffffebf47b95c097bcp-16384, erf(@as(f128, 0x1.c5bf891b4ef6aa68p-16384)));
    try std.testing.expectEqual(0x1.ffffffffffffffe2ed8fed9f4a08p-16384, erf(@as(f128, 0x1.c5bf891b4ef6aa6p-16384)));
    // try std.testing.expectEqual(0x1.ffffffffffffffe77105c1aff0e4p-16384, erf(@as(f128, 0x1.c5bf891b4ef6aa64p-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x2.00000000000000191715de671c5p-16384, erf(@as(f128, 0x1.c5bf891b4ef6aa9p-16384)));
    try std.testing.expectEqual(0x2.0000000000000010102a3645ce98p-16384, erf(@as(f128, 0x1.c5bf891b4ef6aa88p-16384)));
    // try std.testing.expectEqual(0x2.000000000000001493a00a567574p-16384, erf(@as(f128, 0x1.c5bf891b4ef6aa8cp-16384)));
    try std.testing.expectEqual(0x2.0000000000000016d55af45ec8ep-16384, erf(@as(f128, 0x1.c5bf891b4ef6aa8ep-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x4.00000000000000050b917427b40cp-16384, erf(@as(f128, 0x3.8b7f12369ded54f8p-16384)));
    try std.testing.expectEqual(0x3.fffffffffffffffc04a5cc066654p-16384, erf(@as(f128, 0x3.8b7f12369ded54fp-16384)));
    try std.testing.expectEqual(0x4.0000000000000000881ba0170d3p-16384, erf(@as(f128, 0x3.8b7f12369ded54f4p-16384)));
    try std.testing.expectEqual(0x3.ffffffffffffffffffffffffffe8p-16384, erf(@as(f128, 0x3.8b7f12369ded54f38760a41abb5cp-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x4.00000000000000050b917427b40cp-16384, erf(@as(f128, 0x3.8b7f12369ded54f8p-16384)));
    try std.testing.expectEqual(0x3.fffffffffffffffc04a5cc066654p-16384, erf(@as(f128, 0x3.8b7f12369ded54fp-16384)));
    try std.testing.expectEqual(0x4.0000000000000000881ba0170d3p-16384, erf(@as(f128, 0x3.8b7f12369ded54f4p-16384)));
    try std.testing.expectEqual(0x4.0000000000000000000000000018p-16384, erf(@as(f128, 0x3.8b7f12369ded54f38760a41abb88p-16384)));
    try std.testing.expectEqual(0x4.0000000000000000000000000014p-16384, erf(@as(f128, 0x3.8b7f12369ded54f38760a41abb84p-16384)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x1.ap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x1.cp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x6.4p+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x6.ap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x6.a8p+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x6.aap+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x6.bp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x6.cp+4)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0x3.e8p+8)));
    try std.testing.expectEqual(0x9.062b22ee929bfc9c18d570fce48p-8, erf(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(-0x9.062b22ee929bfc9c18d570fce48p-8, erf(@as(f128, -0x8p-8)));
    try std.testing.expectEqual(0x4.8375bbfe32e3ccb857c5e825b18cp-12, erf(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x2.41baea05511f14d8f47394f310bp-16, erf(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0x1.20dd750429568749379b4a46c121p-20, erf(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0x9.06eba8214db5c84379f08c279848p-28, erf(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x4.8375d410a6db44537c2fe8f7e61p-32, erf(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0x2.41baea08536da235c44fdb704f8cp-36, erf(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ae3a8b4b50651p-40, erf(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7362c48p-48, erf(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb47e4p-52, erf(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x2.41baea08536da235c75229fdaff8p-56, erf(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ae3a914fed7fep-60, erf(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ae3a914fed7fep-100, erf(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ae3a914fed7fep-600, erf(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, erf(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.20dd750429b6d11ae3a914fed7fep-10000, erf(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-128, erf(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1024, erf(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-16384, erf(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2.41baea08536da235c75229fdaffcp-16384, erf(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-972, erf(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4.8375d410a6db446b8ea453fb5ff8p-128, erf(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4.8375d410a6db446b8ea453fb5ff8p-1024, erf(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4.8375d410a6db446b8ea453fb5ff8p-16384, erf(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2.41baea08536da235c75229fdaffcp-16384, erf(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x9.06eba8214db688d71d48a7f6bffp-972, erf(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x9.06eba8214db8p-16448, erf(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4.8375d410a6dcp-16448, erf(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, erf(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x9.06eba8214db688d71d48a7f6bffp-152, erf(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4.8375d410a6db446b8ea453fb5ff8p-1076, erf(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x9.06eba8214db8p-16448, erf(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4.8375d410a6dcp-16448, erf(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, erf(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1p+0, erf(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1p+0, erf(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    // try std.testing.expectEqual(-0xf.dded081f58dc8a9153bf0342acc8p-4, erf(@as(f128, -0x1.ddaea4p+0)));
    try std.testing.expectEqual(-0xe.6cc3fab61fed8550eefeff64e4fp-4, erf(@as(f128, -0x1.2b1f68p+0)));
    try std.testing.expectEqual(0xe.d6504b655134fdbfea37252f26dp-4, erf(@as(f128, 0x1.44e722p+0)));
    // try std.testing.expectEqual(-0xe.ad06dfdab8f3efcfd1feb6bc4a88p-4, erf(@as(f128, -0x1.3a0d48p+0)));
    // try std.testing.expectEqual(-0xf.d0e59ffd2b01cc98191c3d9fb338p-4, erf(@as(f128, -0x1.c975cap+0)));
    try std.testing.expectEqual(-0xf.e294502e0c11682ba40c5a30e87p-4, erf(@as(f128, -0x1.e6a006p+0)));
    try std.testing.expectEqual(-0x1.77f98ef609eb313046ceab3fa8c7p-12, erf(@as(f128, -0x1.4d32f4p-12)));
}
