const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
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
    const ax: f32 = float.abs(x);
    const ux: u32 = @bitCast(ax);
    const s: f64 = cast(f64, x, .{});
    var z: f64 = cast(f64, ax, .{});
    // 0x407ad444 corresponds to x = 0x1.f5a888p+1 = 3.91921..., which is the
    // largest float such that erf(x) does not round to 1 (to nearest).
    if (ux > 0x407ad444) {
        @branchHint(.unlikely);
        const os: f32 = float.copysign(@as(f32, 1), x);

        if (ux > (0xff << 23))
            return x + x; // nan

        if (ux == (0xff << 23))
            return os; // +-inf

        return os - 0x1p-25 * os;
    }

    const v: f64 = float.floor(16 * z);
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
    return cast(f32, float.copysign(c0, s), .{});
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

                if (float.abs(ret) < std.math.floatMin(f64)) {
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

    const xx: f64 = float.abs(x);
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
    const r: f64 = float.exp(-z * z - 0.5625) * float.exp((z - xx) * (z + xx) + R / S);

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

        const y: f128 = float.erfc(x);
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

                if (float.abs(ret) < std.math.floatMin(f128)) {
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
