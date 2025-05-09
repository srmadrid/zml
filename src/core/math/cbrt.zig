const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const classify = @import("classify.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn cbrt(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return cbrt(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, cbrt32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_cbrtf.c
                    return cbrt32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_cbrt.c
                    return cbrt64(x);
                },
                f80 => return cast(f80, cbrt128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_cbrtl.c
                    return cbrt128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn cbrt32(x: f32) f32 {
    const escale: [3]f64 = .{
        1.0,
        0x1.428a2f98d728bp+0, // 2^(1/3)
        0x1.965fea53d6e3dp+0, // 2^(2/3)
    };

    const u: u32 = @bitCast(x);
    var au: u32 = u << 1;
    const sgn: u32 = u >> 31;
    var e: u32 = au >> 24;
    if (au < 1 << 24 or au >= 0xff << 24) {
        @branchHint(.unlikely);
        if (au >= 0xff << 24)
            return x + x; // inf, nan

        if (au == 0)
            return x; // +-0

        const nz: i32 = @clz(au) - 7; // subnormal
        au <<= @as(u5, @intCast(nz));
        e -%= cast(u32, nz -% 1, .{});
    }

    const mant: u32 = au & 0xffffff;
    e +%= 899;
    const et: u32 = @divFloor(e, 3);
    const it: u32 = e % 3;
    var isc: u64 = @bitCast(escale[it]);
    isc +%= cast(u64, (et -% 342), .{}) << 52;
    isc |= cast(u64, sgn, .{}) << 63;
    const cvt2: f64 = @bitCast(isc);
    const c: [8]f64 = .{
        0x1.2319d352ea5d5p-1,  0x1.67ad8ee258d1ap-1,  -0x1.9342edf9cbad9p-2,
        0x1.b6388fc510a75p-3,  -0x1.6002455599e2fp-4, 0x1.7b096936192c4p-6,
        -0x1.e5577187e8bf8p-9, 0x1.169ef81d6c34ep-12,
    };
    const z: f64 = @bitCast(cast(u64, mant, .{}) << 28 | 0x3ff << 52);
    const r0: f64 = -0x1.9931c6c2d19d1p-6 / z;
    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    var f: f64 = ((c[0] + z * c[1]) + z2 * (c[2] + z * c[3])) + z4 * ((c[4] + z * c[5]) + z2 * (c[6] + z * c[7])) + r0;
    var r: f64 = f * cvt2;
    var ub: f32 = cast(f32, r, .{});
    const lb: f32 = cast(f32, r - cvt2 * 1.4182e-9, .{});
    if (ub == lb) {
        @branchHint(.likely);
        return ub;
    }

    const U0: f64 = -0x1.ab16ec65d138fp+3;
    const h: f64 = f * f * f - z;
    f -= (f * r0 * U0) * h;
    r = f * cvt2;
    var cvt1: u64 = @bitCast(r);
    ub = cast(f32, r, .{});
    const m0: i64 = @bitCast(cvt1 << 19);
    const m1: i64 = m0 >> 63;
    if ((m0 ^ m1) < (1 << 31)) {
        @branchHint(.unlikely);
        cvt1 = (cvt1 + (1 << 31)) & 0xffffffff00000000;
        ub = cast(f32, @as(f64, @bitCast(cvt1)), .{});
    }
    return ub;
}

fn cbrt64(x: f64) f64 {
    const factor: [5]f64 = .{ 1.0 / 1.5874010519681994748, 1.0 / 1.2599210498948731648, 1.0, 1.2599210498948731648, 1.5874010519681994748 };

    // Reduce X.  XM now is an range 1.0 to 0.5.
    var xe: i32 = undefined;
    const xm = math.frexp(math.abs(x), &xe);

    // If X is not finite or is null return it (with raising exceptions
    // if necessary.
    if (xe == 0 and classify.classify(x) <= classify.ZERO)
        return x + x;

    const u: f64 = (0.354895765043919860 + ((1.50819193781584896 + ((-2.11499494167371287 + ((2.44693122563534430 + ((-1.83469277483613086 + (0.784932344976639262 - 0.145263899385486377 * xm) * xm) * xm)) * xm)) * xm)) * xm));

    const t2: f64 = u * u * u;

    const ym: f64 = u * (t2 + 2.0 * xm) / (2.0 * t2 + xm) * factor[@intCast(2 + @mod(xe, 3))];

    return math.ldexp(if (x > 0) ym else -ym, @divFloor(xe, 3));
}

fn cbrt128(x: f128) f128 {
    const CBRT2: f128 = 1.259921049894873164767210607278228350570251;
    const CBRT4: f128 = 1.587401051968199474751705639272308260391493;
    const CBRT2I: f128 = 0.7937005259840997373758528196361541301957467;
    const CBRT4I: f128 = 0.6299605249474365823836053036391141752851257;

    if (!std.math.isFinite(x))
        return x + x;

    if (x == 0)
        return x;

    var sign: i32 = undefined;
    var xx: f128 = x;
    if (x > 0) {
        sign = 1;
    } else {
        sign = -1;
        xx = -x;
    }

    const z: f128 = xx;
    // extract power of 2, leaving mantissa between 0.5 and 1
    var e: i32 = undefined;
    xx = math.frexp(xx, &e);

    // Approximate cube root of number between .5 and 1,
    // peak relative error = 1.2e-6
    xx = ((((1.3584464340920900529734e-1 * xx - 6.3986917220457538402318e-1) * xx + 1.2875551670318751538055e0) * xx - 1.4897083391357284957891e0) * xx + 1.3304961236013647092521e0) * xx + 3.7568280825958912391243e-1;

    // exponent divided by 3
    if (e >= 0) {
        var rem: i32 = e;
        e = @divFloor(e, 3);
        rem -= 3 * e;
        if (rem == 1) {
            xx *= CBRT2;
        } else if (rem == 2) {
            xx *= CBRT4;
        }
    } else { // argument less than 1
        e = -e;
        var rem: i32 = e;
        e = @divFloor(e, 3);
        rem -= 3 * e;
        if (rem == 1) {
            xx *= CBRT2I;
        } else if (rem == 2) {
            xx *= CBRT4I;
        }
        e = -e;
    }

    // multiply by power of 2
    xx = math.ldexp(xx, e);

    // Newton iteration
    xx -= (xx - (z / (xx * xx))) * 0.3333333333333333333333333333333333333333;
    xx -= (xx - (z / (xx * xx))) * 0.3333333333333333333333333333333333333333;
    xx -= (xx - (z / (xx * xx))) * 0.3333333333333333333333333333333333333333;

    if (sign < 0)
        xx = -xx;

    return xx;
}

test cbrt {
    try std.testing.expectEqual(0x0p+0, cbrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x1.999998p-4, cbrt(@as(f32, -0x4.18937p-12)));
    try std.testing.expectEqual(-0x1.99999ap-4, cbrt(@as(f32, -0x4.189378p-12)));
    try std.testing.expectEqual(0x1.428a3p+0, cbrt(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(0x1.965feap+0, cbrt(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(0x2p+0, cbrt(@as(f32, 0x8p+0)));
    try std.testing.expectEqual(-0x2.278908p+0, cbrt(@as(f32, -0xap+0)));
    try std.testing.expectEqual(-0x3p+0, cbrt(@as(f32, -0x1.bp+4)));
    try std.testing.expectEqual(0xf.f54e3p-4, cbrt(@as(f32, 0xf.ep-4)));
    try std.testing.expectEqual(0xe.89768p-4, cbrt(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x6.597fa8p+40, cbrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x5.0a28cp-52, cbrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, cbrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.e6a77cp+4, cbrt(@as(f32, 0x1.86ap+16)));
    try std.testing.expectEqual(0x1.744268p+0, cbrt(@as(f32, 0x3.132634p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x5.0a28cp-52, cbrt(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(-0x5.80e51p+12, cbrt(@as(f32, -0xa.6b142p+40)));
    try std.testing.expectEqual(-0x3.25909cp-44, cbrt(@as(f32, -0x1.f28ab8p-128)));
    try std.testing.expectEqual(-0x3.2590ap-44, cbrt(@as(f32, -0x1.f28acp-128)));
    try std.testing.expectEqual(-0x1.64ebb2p-12, cbrt(@as(f32, -0x2.b5cd28p-36)));
    try std.testing.expectEqual(-0x3.1642p-8, cbrt(@as(f32, -0x1.d6a8bep-20)));
    try std.testing.expectEqual(-0x1.7eff28p-24, cbrt(@as(f32, -0x3.593ed8p-72)));
    try std.testing.expectEqual(0x3.07a108p-36, cbrt(@as(f32, 0x1.bd0098p-104)));
    try std.testing.expectEqual(-0x1.78c2ccp+0, cbrt(@as(f32, -0x3.300d34p+0)));
    try std.testing.expectEqual(0xb.a0f06p-4, cbrt(@as(f32, 0x6.247f5p-4)));
    try std.testing.expectEqual(-0x1.7c7862p+0, cbrt(@as(f32, -0x3.48648p+0)));
    try std.testing.expectEqual(-0x1.7c7864p+0, cbrt(@as(f32, -0x3.486484p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x5.0a28cp-52, cbrt(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.597fa8p+40, cbrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x6.597fa8p+40, cbrt(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x4p-44, cbrt(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-44, cbrt(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x5.0a28cp-52, cbrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x5.0a28cp-52, cbrt(@as(f32, -0x8p-152)));

    try std.testing.expectEqual(0x0p+0, cbrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x1.999998fbbbbb8p-4, cbrt(@as(f64, -0x4.18937p-12)));
    try std.testing.expectEqual(-0x1.99999a0666665p-4, cbrt(@as(f64, -0x4.189378p-12)));
    try std.testing.expectEqual(-0x1.9999999999999p-4, cbrt(@as(f64, -0x4.189374bc6a7ecp-12)));
    try std.testing.expectEqual(-0x1.999999999999ap-4, cbrt(@as(f64, -0x4.189374bc6a7fp-12)));
    // try std.testing.expectEqual(0x1.428a2f98d728bp+0, cbrt(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(0x1.965fea53d6e3dp+0, cbrt(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(0x2p+0, cbrt(@as(f64, 0x8p+0)));
    // try std.testing.expectEqual(-0x2.278908270e09ep+0, cbrt(@as(f64, -0xap+0)));
    // try std.testing.expectEqual(-0x3p+0, cbrt(@as(f64, -0x1.bp+4)));
    // try std.testing.expectEqual(0xf.f54e30f23e698p-4, cbrt(@as(f64, 0xf.ep-4)));
    // try std.testing.expectEqual(0xe.89768578d13f8p-4, cbrt(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x6.597fa7318656p+40, cbrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.85145f31ae516p+340, cbrt(@as(f64, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(0x5.0a28be635ca2cp-52, cbrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, cbrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-360, cbrt(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x2.e6a77a87274eap+4, cbrt(@as(f64, 0x1.86ap+16)));
    try std.testing.expectEqual(0x1.744267cbadff7p+0, cbrt(@as(f64, 0x3.132634p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f64, -0x0p+0)));
    // try std.testing.expectEqual(-0x5.0a28be635ca2cp-52, cbrt(@as(f64, -0x8p-152)));
    // try std.testing.expectEqual(-0x2.52ed87c91737ep-164, cbrt(@as(f64, -0xc.8d0442f2f0d1p-492)));
    // try std.testing.expectEqual(-0x5.80e513d887c78p+12, cbrt(@as(f64, -0xa.6b142p+40)));
    // try std.testing.expectEqual(-0x3.25909b23791cp-44, cbrt(@as(f64, -0x1.f28ab8p-128)));
    try std.testing.expectEqual(-0x3.25909f728def4p-44, cbrt(@as(f64, -0x1.f28acp-128)));
    try std.testing.expectEqual(-0x3.25909b56c104cp-44, cbrt(@as(f64, -0x1.f28ab85f3580ap-128)));
    try std.testing.expectEqual(-0x1.64ebb100c787bp-12, cbrt(@as(f64, -0x2.b5cd28p-36)));
    try std.testing.expectEqual(-0x3.164200fbbcb72p-8, cbrt(@as(f64, -0x1.d6a8bep-20)));
    // try std.testing.expectEqual(-0x1.7eff2881c395fp-24, cbrt(@as(f64, -0x3.593ed8p-72)));
    // try std.testing.expectEqual(0x3.07a108f56532ap-36, cbrt(@as(f64, 0x1.bd0098p-104)));
    // try std.testing.expectEqual(-0x1.78c2cb7ea3cdfp+0, cbrt(@as(f64, -0x3.300d34p+0)));
    try std.testing.expectEqual(0xb.a0f06280ab958p-4, cbrt(@as(f64, 0x6.247f5p-4)));
    try std.testing.expectEqual(-0x1.7c7862d51e30fp+0, cbrt(@as(f64, -0x3.48648p+0)));
    // try std.testing.expectEqual(-0x1.7c78636fa6457p+0, cbrt(@as(f64, -0x3.486484p+0)));
    // try std.testing.expectEqual(-0x1.7c7862db462edp+0, cbrt(@as(f64, -0x3.48648028cb464p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f64, -0x0p+0)));
    // try std.testing.expectEqual(-0x5.0a28be635ca2cp-52, cbrt(@as(f64, -0x8p-152)));
    // try std.testing.expectEqual(-0xb.81d0965bf918p-80, cbrt(@as(f64, -0x5.f3b076ad049c8p-232)));
    try std.testing.expectEqual(0x6.597fa7318656p+40, cbrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.85145f31ae516p+340, cbrt(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x6.597fa7318656p+40, cbrt(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.85145f31ae516p+340, cbrt(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4p-44, cbrt(@as(f64, 0x4p-128)));
    // try std.testing.expectEqual(0xa.14517cc6b9458p-344, cbrt(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x2p-324, cbrt(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-44, cbrt(@as(f64, -0x4p-128)));
    // try std.testing.expectEqual(-0xa.14517cc6b9458p-344, cbrt(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x2p-324, cbrt(@as(f64, -0x8p-972)));
    // try std.testing.expectEqual(0x5.0a28be635ca2cp-52, cbrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-360, cbrt(@as(f64, 0x4p-1076)));
    // try std.testing.expectEqual(-0x5.0a28be635ca2cp-52, cbrt(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-360, cbrt(@as(f64, -0x4p-1076)));

    try std.testing.expectEqual(0x0p+0, cbrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x1.999998fbbbbb7ee4p-4, cbrt(@as(f80, -0x4.18937p-12)));
    try std.testing.expectEqual(-0x1.99999a066666498p-4, cbrt(@as(f80, -0x4.189378p-12)));
    try std.testing.expectEqual(-0x1.9999999999999212p-4, cbrt(@as(f80, -0x4.189374bc6a7ecp-12)));
    try std.testing.expectEqual(-0x1.9999999999999a66p-4, cbrt(@as(f80, -0x4.189374bc6a7fp-12)));
    try std.testing.expectEqual(-0x1.999999999999999ap-4, cbrt(@as(f80, -0x4.189374bc6a7ef9d8p-12)));
    try std.testing.expectEqual(-0x1.999999999999999ap-4, cbrt(@as(f80, -0x4.189374bc6a7ef9ep-12)));
    try std.testing.expectEqual(0x1.428a2f98d728ae22p+0, cbrt(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(0x1.965fea53d6e3c82cp+0, cbrt(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(0x2p+0, cbrt(@as(f80, 0x8p+0)));
    try std.testing.expectEqual(-0x2.278908270e09d95p+0, cbrt(@as(f80, -0xap+0)));
    try std.testing.expectEqual(-0x3p+0, cbrt(@as(f80, -0x1.bp+4)));
    try std.testing.expectEqual(0xf.f54e30f23e69be4p-4, cbrt(@as(f80, 0xf.ep-4)));
    try std.testing.expectEqual(0xe.89768578d13f79fp-4, cbrt(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x6.597fa7318655fc48p+40, cbrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.85145f31ae51558cp+340, cbrt(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2p+5460, cbrt(@as(f80, 0x8p+16380)));
    try std.testing.expectEqual(0x5.0a28be635ca2b888p-52, cbrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, cbrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-360, cbrt(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-5464, cbrt(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x2.e6a77a87274eadc8p+4, cbrt(@as(f80, 0x1.86ap+16)));
    try std.testing.expectEqual(0x1.744267cbadff73aap+0, cbrt(@as(f80, 0x3.132634p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x5.0a28be635ca2b888p-52, cbrt(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x2.52ed87c91737d8c4p-164, cbrt(@as(f80, -0xc.8d0442f2f0d1p-492)));
    try std.testing.expectEqual(-0x5.80e513d887c77e3p+12, cbrt(@as(f80, -0xa.6b142p+40)));
    try std.testing.expectEqual(-0x3.25909b23791c01fp-44, cbrt(@as(f80, -0x1.f28ab8p-128)));
    try std.testing.expectEqual(-0x3.25909f728def3054p-44, cbrt(@as(f80, -0x1.f28acp-128)));
    try std.testing.expectEqual(-0x3.25909b56c104c22p-44, cbrt(@as(f80, -0x1.f28ab85f3580ap-128)));
    try std.testing.expectEqual(-0x1.64ebb100c787b01ep-12, cbrt(@as(f80, -0x2.b5cd28p-36)));
    try std.testing.expectEqual(-0x3.164200fbbcb7214cp-8, cbrt(@as(f80, -0x1.d6a8bep-20)));
    try std.testing.expectEqual(-0x1.7eff2881c395ed16p-24, cbrt(@as(f80, -0x3.593ed8p-72)));
    try std.testing.expectEqual(0x3.07a108f565329ec8p-36, cbrt(@as(f80, 0x1.bd0098p-104)));
    try std.testing.expectEqual(-0x1.78c2cb7ea3cdf6cap+0, cbrt(@as(f80, -0x3.300d34p+0)));
    try std.testing.expectEqual(0xb.a0f06280ab95b38p-4, cbrt(@as(f80, 0x6.247f5p-4)));
    try std.testing.expectEqual(-0x1.7c7862d51e30ebf8p+0, cbrt(@as(f80, -0x3.48648p+0)));
    try std.testing.expectEqual(-0x1.7c78636fa645682p+0, cbrt(@as(f80, -0x3.486484p+0)));
    try std.testing.expectEqual(-0x1.7c7862db462ecf3cp+0, cbrt(@as(f80, -0x3.48648028cb464p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x5.0a28be635ca2b888p-52, cbrt(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0xb.81d0965bf917c9ep-80, cbrt(@as(f80, -0x5.f3b076ad049c8p-232)));
    try std.testing.expectEqual(0x6.597fa7318655fc48p+40, cbrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.85145f31ae51558cp+340, cbrt(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.85145f31ae515c44p+5460, cbrt(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x6.597fa7318655fc48p+40, cbrt(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.85145f31ae51558cp+340, cbrt(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.85145f31ae515c44p+5460, cbrt(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4p-44, cbrt(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0xa.14517cc6b945711p-344, cbrt(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0xa.14517cc6b945711p-5464, cbrt(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x8p-5464, cbrt(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x2p-324, cbrt(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-44, cbrt(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0xa.14517cc6b945711p-344, cbrt(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0xa.14517cc6b945711p-5464, cbrt(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x8p-5464, cbrt(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x2p-324, cbrt(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x5.0a28be635ca2b888p-52, cbrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-360, cbrt(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x5.0a28be635ca2b888p-5484, cbrt(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x5.0a28be635ca2b888p-52, cbrt(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-360, cbrt(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x5.0a28be635ca2b888p-5484, cbrt(@as(f80, -0x8p-16448)));

    try std.testing.expectEqual(0x0p+0, cbrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x1.999998fbbbbb7ee38e11ce06340cp-4, cbrt(@as(f128, -0x4.18937p-12)));
    try std.testing.expectEqual(-0x1.99999a0666664980000ccb554e89p-4, cbrt(@as(f128, -0x4.189378p-12)));
    try std.testing.expectEqual(-0x1.99999999999992111111111110eep-4, cbrt(@as(f128, -0x4.189374bc6a7ecp-12)));
    try std.testing.expectEqual(-0x1.9999999999999a66666666666666p-4, cbrt(@as(f128, -0x4.189374bc6a7fp-12)));
    try std.testing.expectEqual(-0x1.9999999999999999311111111111p-4, cbrt(@as(f128, -0x4.189374bc6a7ef9d8p-12)));
    try std.testing.expectEqual(-0x1.999999999999999a3bbbbbbbbbbcp-4, cbrt(@as(f128, -0x4.189374bc6a7ef9ep-12)));
    try std.testing.expectEqual(-0x1.9999999999999999999999999999p-4, cbrt(@as(f128, -0x4.189374bc6a7ef9db22d0e5604188p-12)));
    try std.testing.expectEqual(-0x1.999999999999999999999999999ap-4, cbrt(@as(f128, -0x4.189374bc6a7ef9db22d0e560418cp-12)));
    try std.testing.expectEqual(-0x1.9999999999999999999999999966p-4, cbrt(@as(f128, -0x4.189374bc6a7ef9db22d0e5604p-12)));
    try std.testing.expectEqual(-0x1.99999999999999999999999999a9p-4, cbrt(@as(f128, -0x4.189374bc6a7ef9db22d0e56042p-12)));
    try std.testing.expectEqual(0x1.428a2f98d728ae223ddab715be25p+0, cbrt(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(0x1.965fea53d6e3c82b05999ab43dc5p+0, cbrt(@as(f128, 0x4p+0)));
    try std.testing.expectEqual(0x2p+0, cbrt(@as(f128, 0x8p+0)));
    try std.testing.expectEqual(-0x2.278908270e09d951445ae49bd412p+0, cbrt(@as(f128, -0xap+0)));
    try std.testing.expectEqual(-0x3p+0, cbrt(@as(f128, -0x1.bp+4)));
    try std.testing.expectEqual(0xf.f54e30f23e69be3850ca030dc7bp-4, cbrt(@as(f128, 0xf.ep-4)));
    try std.testing.expectEqual(0xe.89768578d13f79ed5d709a616d18p-4, cbrt(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x6.597fa7318655fc467e27422a246p+40, cbrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.85145f31ae51558c45623f054decp+340, cbrt(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2p+5460, cbrt(@as(f128, 0x8p+16380)));
    try std.testing.expectEqual(0x2.85145f31ae5158e8608bd69864eap+340, cbrt(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x5.0a28be635ca2b888f76adc56f894p-52, cbrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, cbrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-360, cbrt(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-5464, cbrt(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x2.e6a77a87274eadc9b39cffd8ab96p+4, cbrt(@as(f128, 0x1.86ap+16)));
    try std.testing.expectEqual(0x1.744267cbadff73aa2b2ff2839fe4p+0, cbrt(@as(f128, 0x3.132634p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x5.0a28be635ca2b888f76adc56f894p-52, cbrt(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x2.52ed87c91737d8c3efb34d732616p-164, cbrt(@as(f128, -0xc.8d0442f2f0d1p-492)));
    try std.testing.expectEqual(-0x5.80e513d887c77e32360beb1684d4p+12, cbrt(@as(f128, -0xa.6b142p+40)));
    try std.testing.expectEqual(-0x3.25909b23791c01f1d682fd88edbp-44, cbrt(@as(f128, -0x1.f28ab8p-128)));
    // try std.testing.expectEqual(-0x3.25909f728def3054e224bdd8a1c2p-44, cbrt(@as(f128, -0x1.f28acp-128)));
    // try std.testing.expectEqual(-0x3.25909b56c104c21f46cb785a46f6p-44, cbrt(@as(f128, -0x1.f28ab85f3580ap-128)));
    try std.testing.expectEqual(-0x1.64ebb100c787b01d76587f2a6af2p-12, cbrt(@as(f128, -0x2.b5cd28p-36)));
    try std.testing.expectEqual(-0x3.164200fbbcb7214abb0ade43ea94p-8, cbrt(@as(f128, -0x1.d6a8bep-20)));
    try std.testing.expectEqual(-0x1.7eff2881c395ed16ad29095a964cp-24, cbrt(@as(f128, -0x3.593ed8p-72)));
    try std.testing.expectEqual(0x3.07a108f565329ec8fffea3a6b35p-36, cbrt(@as(f128, 0x1.bd0098p-104)));
    try std.testing.expectEqual(-0x1.78c2cb7ea3cdf6c95160af9c9403p+0, cbrt(@as(f128, -0x3.300d34p+0)));
    try std.testing.expectEqual(0xb.a0f06280ab95b378e4672472e74p-4, cbrt(@as(f128, 0x6.247f5p-4)));
    try std.testing.expectEqual(-0x1.7c7862d51e30ebf738126e131417p+0, cbrt(@as(f128, -0x3.48648p+0)));
    try std.testing.expectEqual(-0x1.7c78636fa6456820b04e845d7b13p+0, cbrt(@as(f128, -0x3.486484p+0)));
    try std.testing.expectEqual(-0x1.7c7862db462ecf3cf63a9a1d7f4ap+0, cbrt(@as(f128, -0x3.48648028cb464p+0)));
    try std.testing.expectEqual(-0x0p+0, cbrt(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x5.0a28be635ca2b888f76adc56f894p-52, cbrt(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0xb.81d0965bf917c9dad6c6858ba3d8p-80, cbrt(@as(f128, -0x5.f3b076ad049c8p-232)));
    try std.testing.expectEqual(0x6.597fa7318655fc467e27422a246p+40, cbrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.85145f31ae51558c45623f054decp+340, cbrt(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.85145f31ae515c43a4aea3c59784p+5460, cbrt(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2.85145f31ae515c447bb56e2b7c4ap+5460, cbrt(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.85145f31ae5158e8608bd69864eap+340, cbrt(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x6.597fa7318655fc467e27422a246p+40, cbrt(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.85145f31ae51558c45623f054decp+340, cbrt(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.85145f31ae515c43a4aea3c59784p+5460, cbrt(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x2.85145f31ae515c447bb56e2b7c4ap+5460, cbrt(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x2.85145f31ae5158e8608bd69864eap+340, cbrt(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4p-44, cbrt(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0xa.14517cc6b9457111eed5b8adf128p-344, cbrt(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0xa.14517cc6b9457111eed5b8adf128p-5464, cbrt(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x8p-5464, cbrt(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x2p-324, cbrt(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-44, cbrt(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0xa.14517cc6b9457111eed5b8adf128p-344, cbrt(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0xa.14517cc6b9457111eed5b8adf128p-5464, cbrt(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x8p-5464, cbrt(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x2p-324, cbrt(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x5.0a28be635ca2b888f76adc56f894p-52, cbrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-360, cbrt(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x5.0a28be635ca2b888f76adc56f894p-5484, cbrt(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-5484, cbrt(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-5500, cbrt(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x5.0a28be635ca2b888f76adc56f894p-52, cbrt(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-360, cbrt(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x5.0a28be635ca2b888f76adc56f894p-5484, cbrt(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-5484, cbrt(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-5500, cbrt(@as(f128, -0x4p-16496)));
}
