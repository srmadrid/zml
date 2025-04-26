const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const roundeven = @import("roundeven.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn atanh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return atanh(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, atanh32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_atanhf.c
                    return atanh32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_atanh.c
                    return atanh64(x);
                },
                f80 => return cast(f80, atanh128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_atanhl.c
                    return atanh128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn as_special(x: f32) f32 {
    const ix: u32 = @bitCast(x);
    const ax: u32 = ix << 1;

    if (ax == 0x7f000000) // +-1
        return if ((ix >> 31) != 0) -1 / x else 1 / x;

    if (ax > 0xff000000)
        return x + x; // nan

    return (x - x) / (x - x);
}

fn atanh32(x: f32) f32 {
    // Calculate atanh(x) using the difference of two logarithms -- atanh(x) =
    // (ln(1+x) - ln(1-x))/2
    const tr: [64]f64 = .{
        0x1.fc07f02p-1, 0x1.f44659ep-1, 0x1.ecc07b3p-1, 0x1.e573ac9p-1,
        0x1.de5d6e4p-1, 0x1.d77b655p-1, 0x1.d0cb58fp-1, 0x1.ca4b305p-1,
        0x1.c3f8f02p-1, 0x1.bdd2b8ap-1, 0x1.b7d6c3ep-1, 0x1.b20364p-1,
        0x1.ac5701bp-1, 0x1.a6d01a7p-1, 0x1.a16d3f9p-1, 0x1.9c2d14fp-1,
        0x1.970e4f8p-1, 0x1.920fb4ap-1, 0x1.8d3018dp-1, 0x1.886e5f1p-1,
        0x1.83c977bp-1, 0x1.7f405fdp-1, 0x1.7ad2209p-1, 0x1.767dce4p-1,
        0x1.724287fp-1, 0x1.6e1f76bp-1, 0x1.6a13cd1p-1, 0x1.661ec6ap-1,
        0x1.623fa77p-1, 0x1.5e75bb9p-1, 0x1.5ac056bp-1, 0x1.571ed3cp-1,
        0x1.5390949p-1, 0x1.5015015p-1, 0x1.4cab887p-1, 0x1.49539e4p-1,
        0x1.460cbc8p-1, 0x1.42d6626p-1, 0x1.3fb014p-1,  0x1.3c995a4p-1,
        0x1.3991c2cp-1, 0x1.3698df4p-1, 0x1.33ae45bp-1, 0x1.30d1901p-1,
        0x1.2e025cp-1,  0x1.2b404adp-1, 0x1.288b013p-1, 0x1.25e2271p-1,
        0x1.2345679p-1, 0x1.20b470cp-1, 0x1.1e2ef3bp-1, 0x1.1bb4a4p-1,
        0x1.1945381p-1, 0x1.16e0689p-1, 0x1.1485f0ep-1, 0x1.12358e7p-1,
        0x1.0fef011p-1, 0x1.0db20a9p-1, 0x1.0b7e6ecp-1, 0x1.0953f39p-1,
        0x1.073260ap-1, 0x1.05197f8p-1, 0x1.03091b5p-1, 0x1.010101p-1,
    };
    const tl: [64]f64 = .{
        0x1.fe02a69106789p-9, 0x1.7b91b1155b11bp-7, 0x1.39e87ba1ebd6p-6,
        0x1.b42dd713971bfp-6, 0x1.16536ee637ae1p-5, 0x1.51b073c96183fp-5,
        0x1.8c345da019b21p-5, 0x1.c5e5492abc743p-5, 0x1.fec912fbbeabbp-5,
        0x1.1b72ad33f67ap-4,  0x1.371fc1f6e8f74p-4, 0x1.526e5e5a1b438p-4,
        0x1.6d60fe601d21dp-4, 0x1.87fa06438c911p-4, 0x1.a23bc223ab563p-4,
        0x1.bc28673a58cd6p-4, 0x1.d5c216b8fbb91p-4, 0x1.ef0adcaec5936p-4,
        0x1.040259530d041p-3, 0x1.1058bf8d24ad5p-3, 0x1.1c898c09d99fbp-3,
        0x1.2895a13e286a3p-3, 0x1.347dd9a447d55p-3, 0x1.404308716a7e4p-3,
        0x1.4be5f963b78a1p-3, 0x1.5767718015a6cp-3, 0x1.62c82f3a5c795p-3,
        0x1.6e08eab13a1e4p-3, 0x1.792a55fe147a2p-3, 0x1.842d1d9928b17p-3,
        0x1.8f11e873a62c7p-3, 0x1.99d958207e08bp-3, 0x1.a484090c1bb0ap-3,
        0x1.af129324b786bp-3, 0x1.b9858970710fbp-3, 0x1.c3dd7a6ddad4dp-3,
        0x1.ce1af0b65f3ebp-3, 0x1.d83e725022f3ep-3, 0x1.e2488197c6c26p-3,
        0x1.ec399d3d68ccp-3,  0x1.f6123fac028acp-3, 0x1.ffd2e07e7f498p-3,
        0x1.04bdf9e3b26d2p-2, 0x1.0986f4fa93521p-2, 0x1.0e4498651cc8cp-2,
        0x1.12f719595efbcp-2, 0x1.179eabb0a99a1p-2, 0x1.1c3b81e933c25p-2,
        0x1.20cdcd0e0ab6ep-2, 0x1.2555bcf50f7cbp-2, 0x1.29d37ff34b08bp-2,
        0x1.2e47437640268p-2, 0x1.32b1338401d71p-2, 0x1.37117b5c147b6p-2,
        0x1.3b6844a13fc23p-2, 0x1.3fb5b857f6f42p-2, 0x1.43f9fe2f7ce67p-2,
        0x1.48353d11488dfp-2, 0x1.4c679b014ee3ap-2, 0x1.50913cc03686bp-2,
        0x1.54b2468259498p-2, 0x1.58cadb57d7989p-2, 0x1.5cdb1dcaa1765p-2,
        0x1.60e32f46788d9p-2,
    };
    const ln2n: [24]f64 = .{
        0x1.62e42fedb2a44p-2, 0x1.62e42feeab21ap-1, 0x1.0a2b23f33e789p+0,
        0x1.62e42fef27604p+0, 0x1.bb9d3beb1048p+0,  0x1.0a2b23f37c97ep+1,
        0x1.3687a9f1710bcp+1, 0x1.62e42fef657fap+1, 0x1.8f40b5ed59f38p+1,
        0x1.bb9d3beb4e676p+1, 0x1.e7f9c1e942db4p+1, 0x1.0a2b23f39ba79p+2,
        0x1.205966f295e18p+2, 0x1.3687a9f1901b7p+2, 0x1.4cb5ecf08a556p+2,
        0x1.62e42fef848f5p+2, 0x1.791272ee7ec93p+2, 0x1.8f40b5ed79032p+2,
        0x1.a56ef8ec733d1p+2, 0x1.bb9d3beb6d77p+2,  0x1.d1cb7eea67b0fp+2,
        0x1.e7f9c1e961eaep+2, 0x1.fe2804e85c24dp+2, 0x1.0a2b23f3ab2f6p+3,
    };
    const b: [3]f64 = .{ 0x1.fffffffce5a6ap-2, -0x1.0001f81ec0ab8p-2, 0x1.555a0f53d79a5p-3 };
    const s: [2]f64 = .{ 1, -1 };

    const ux: u32 = @bitCast(x);
    const ax: u32 = ux << 1;
    if (ax < 0x7a300000 or ax >= 0x7f000000) {
        @branchHint(.unlikely);
        if (ax >= 0x7f000000) {
            @branchHint(.unlikely);
            return as_special(x);
        }

        if (ax < 0x73713744) {
            @branchHint(.unlikely);
            if (ax == 0)
                return x; // x = +-0
            return @mulAdd(f32, x, 0x1p-25, x); // |x| < 0.000352112(0x1.713744p-12)
        } else { // |x| < 0x1.3p-5
            const c: [4]f64 = .{
                0x1.5555555555527p-2, 0x1.9999999ba4ee8p-3,
                0x1.24922c280990ap-3, 0x1.c8236aae809c6p-4,
            };
            const z: f64 = cast(f64, x, .{});
            const z2: f64 = z * z;
            const z4: f64 = z2 * z2;
            const r: f64 = c[0] + z2 * c[1] + z4 * (c[2] + z2 * c[3]);
            return cast(f32, z + (z * z2) * r, .{});
        }
    }

    const sgn: f64 = s[ux >> 31];
    const e: i32 = @bitCast(ax >> 24);
    const md: u32 = ((ux << 8) | 1 << 31) >> @as(u5, @intCast(126 - e));
    var mn: u32 = -%md;
    const nz: u5 = @intCast(@clz(mn) + 1);
    mn <<= nz;
    const jn: u32 = mn >> 26;
    const jd: u32 = md >> 26;
    const tn: f64 = @bitCast((cast(i64, mn, .{}) << 20) | (1023 << 52));
    const td: f64 = @bitCast((cast(i64, md, .{}) << 20) | (1023 << 52));
    const zn: f64 = tn * tr[jn] - 1;
    const zd: f64 = td * tr[jd] - 1;
    const zn2: f64 = zn * zn;
    const zd2: f64 = zd * zd;
    const rn: f64 = ((tl[jn] - ln2n[nz - 1]) + zn * b[0]) + zn2 * (b[1] + zn * b[2]);
    const rd: f64 = (tl[jd] + zd * b[0]) + zd2 * (b[1] + zd * b[2]);
    var r: f64 = sgn * (rd - rn);
    var ub: f32 = cast(f32, r, .{});
    const lb: f32 = cast(f32, r + sgn * 0.226e-9, .{});
    if (ub != lb) {
        @branchHint(.unlikely);
        const c: [7]f64 = .{
            0x1p-1,                -0x1.000000000001bp-2, 0x1.55555555555bap-3,
            -0x1.fffffff26d72ep-4, 0x1.99999989035p-4,    -0x1.555c39cb9ee8p-4,
            0x1.24992d8b014a1p-4,
        };
        const zn4: f64 = zn2 * zn2;
        const zd4: f64 = zd2 * zd2;
        var ffn: f64 = zn * (((c[0] + zn * c[1]) + zn2 * (c[2] + zn * c[3])) + zn4 * ((c[4] + zn * c[5]) + zn2 * c[6]));
        ffn += 0x1.0ca86c3898dp-50 * @as(f64, @floatFromInt(nz));
        ffn += tl[jn];
        const en: f64 = @as(f64, @floatFromInt(nz)) * 0x1.62e42fefa3ap-2;
        var fd: f64 = zd * (((c[0] + zd * c[1]) + zd2 * (c[2] + zd * c[3])) + zd4 * ((c[4] + zd * c[5]) + zd2 * c[6]));
        fd += tl[jd];
        r = fd - ffn + en;
        ub = cast(f32, sgn * r, .{});
    }

    return ub;
}

fn atanh64(x: f64) f64 {
    const huge: f64 = 1e300;

    const xa: f64 = math.abs(x);
    var t: f64 = undefined;
    if (xa < 0.5) {
        if (xa < 0x1.0p-28) {
            @branchHint(.unlikely);
            std.mem.doNotOptimizeAway(huge + x);
            if (math.abs(x) < std.math.floatMin(f64)) {
                const vx: f64 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            return x;
        }

        t = xa + xa;
        t = 0.5 * math.log1p(t + t * xa / (1 - xa));
    } else if (xa < 1) {
        @branchHint(.likely);
        t = 0.5 * math.log1p((xa + xa) / (1 - xa));
    } else {
        if (xa > 1.0)
            return (x - x) / (x - x);

        return x / 0;
    }

    return math.copysign(t, x);
}

fn atanh128(x: f128) f128 {
    const huge: f128 = 1e4900;

    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const jx: u32 = u.w0;
    const ix: u32 = jx & 0x7fffffff;
    u.w0 = ix;
    if (ix >= 0x3fff0000) { // |x| >= 1.0 or infinity or NaN
        if (@as(f128, @bitCast(u)) == 1) {
            return x / 0;
        } else {
            return (x - x) / (x - x);
        }
    }

    if (ix < 0x3fc60000 and (huge + x) > 0) { // x < 2^-57
        if (math.abs(x) < std.math.floatMin(f128)) {
            const vx: f128 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }

        return x;
    }

    var t: f128 = undefined;
    if (ix < 0x3ffe0000) { // x < 0.5
        t = @as(f128, @bitCast(u)) + @as(f128, @bitCast(u));
        t = 0.5 * math.log1p(t + t * @as(f128, @bitCast(u)) / (1 - @as(f128, @bitCast(u))));
    } else {
        t = 0.5 * math.log1p((@as(f128, @bitCast(u)) + @as(f128, @bitCast(u))) / (1 - @as(f128, @bitCast(u))));
    }

    return if ((jx & 0x80000000) != 0) -t else t;
}

test atanh {
    try std.testing.expectEqual(0x0p+0, atanh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0xf.91395p-4, atanh(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(-0xf.91395p-4, atanh(@as(f32, -0xcp-4)));
    try std.testing.expectEqual(0x4.162bcp-4, atanh(@as(f32, 0x4p-4)));
    try std.testing.expectEqual(0x8.00aacp-8, atanh(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x4.000018p-12, atanh(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x1.2345p-20, atanh(@as(f32, 0x1.2345p-20)));
    try std.testing.expectEqual(0x1.000056p-8, atanh(@as(f32, 0x1p-8)));
    try std.testing.expectEqual(0x8.0000bp-12, atanh(@as(f32, 0x8p-12)));
    try std.testing.expectEqual(0x4.000018p-12, atanh(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x2.000004p-12, atanh(@as(f32, 0x2p-12)));
    try std.testing.expectEqual(0x1p-12, atanh(@as(f32, 0x1p-12)));
    try std.testing.expectEqual(0x8p-16, atanh(@as(f32, 0x8p-16)));
    try std.testing.expectEqual(0x1p-24, atanh(@as(f32, 0x1p-24)));
    try std.testing.expectEqual(0x8p-28, atanh(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x4p-28, atanh(@as(f32, 0x4p-28)));
    try std.testing.expectEqual(0x2p-28, atanh(@as(f32, 0x2p-28)));
    try std.testing.expectEqual(0x1p-28, atanh(@as(f32, 0x1p-28)));
    try std.testing.expectEqual(0x8p-32, atanh(@as(f32, 0x8p-32)));
    try std.testing.expectEqual(0x4p-32, atanh(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x2p-32, atanh(@as(f32, 0x2p-32)));
    try std.testing.expectEqual(0x1p-32, atanh(@as(f32, 0x1p-32)));
    try std.testing.expectEqual(0x8p-36, atanh(@as(f32, 0x8p-36)));
    try std.testing.expectEqual(0x1p-48, atanh(@as(f32, 0x1p-48)));
    try std.testing.expectEqual(0x8p-52, atanh(@as(f32, 0x8p-52)));
    try std.testing.expectEqual(0x4p-52, atanh(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x2p-52, atanh(@as(f32, 0x2p-52)));
    try std.testing.expectEqual(0x1p-52, atanh(@as(f32, 0x1p-52)));
    try std.testing.expectEqual(0x8p-56, atanh(@as(f32, 0x8p-56)));
    try std.testing.expectEqual(0x4p-56, atanh(@as(f32, 0x4p-56)));
    try std.testing.expectEqual(0x2p-56, atanh(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1p-56, atanh(@as(f32, 0x1p-56)));
    try std.testing.expectEqual(0x8p-60, atanh(@as(f32, 0x8p-60)));
    try std.testing.expectEqual(0x4p-60, atanh(@as(f32, 0x4p-60)));
    try std.testing.expectEqual(0x2p-60, atanh(@as(f32, 0x2p-60)));
    try std.testing.expectEqual(0x1p-100, atanh(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, atanh(@as(f32, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x8.aa123p+0, atanh(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa123p+0, atanh(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa123p+0, atanh(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa123p+0, atanh(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa123p+0, atanh(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa123p+0, atanh(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa123p+0, atanh(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa123p+0, atanh(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa123p+0, atanh(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa123p+0, atanh(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa123p+0, atanh(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa123p+0, atanh(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa123p+0, atanh(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa123p+0, atanh(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x6.e6c77p-20, atanh(@as(f32, -0x6.e6c77p-20)));
    try std.testing.expectEqual(0x3.379438p-4, atanh(@as(f32, 0x3.2ca824p-4)));
    try std.testing.expectEqual(-0x1.ce10a2p-4, atanh(@as(f32, -0x1.cc1d66p-4)));
    try std.testing.expectEqual(-0x2.89e0a4p+0, atanh(@as(f32, -0xf.cd38p-4)));
    try std.testing.expectEqual(-0x2.89e0ccp+0, atanh(@as(f32, -0xf.cd381p-4)));
    try std.testing.expectEqual(-0x1.054e2p-4, atanh(@as(f32, -0x1.04f386p-4)));
    try std.testing.expectEqual(-0x2.0b18b4p-4, atanh(@as(f32, -0x2.084568p-4)));
    try std.testing.expectEqual(-0x3.f4cbc4p-4, atanh(@as(f32, -0x3.e0a5d8p-4)));
    try std.testing.expectEqual(0x3.f3c8bp-4, atanh(@as(f32, 0x3.dfb1f8p-4)));
    try std.testing.expectEqual(0x3.f3c8acp-4, atanh(@as(f32, 0x3.dfb1f4p-4)));
    try std.testing.expectEqual(0x2.286e7cp-4, atanh(@as(f32, 0x2.251b2cp-4)));
    try std.testing.expectEqual(0x2.286e78p-4, atanh(@as(f32, 0x2.251b28p-4)));
    try std.testing.expectEqual(-0x2.eb75acp-4, atanh(@as(f32, -0x2.e3458cp-4)));
    try std.testing.expectEqual(0x3.a17be8p-4, atanh(@as(f32, 0x3.91d9f4p-4)));
    try std.testing.expectEqual(0x3.a17be4p-4, atanh(@as(f32, 0x3.91d9fp-4)));
    try std.testing.expectEqual(-0x2.7121d4p-4, atanh(@as(f32, -0x2.6c52cp-4)));
    try std.testing.expectEqual(-0x2.7121d8p-4, atanh(@as(f32, -0x2.6c52c4p-4)));
    try std.testing.expectEqual(0x3.b2f9d8p-4, atanh(@as(f32, 0x3.a274ecp-4)));
    try std.testing.expectEqual(-0x3.f1098p-8, atanh(@as(f32, -0x3.f0f518p-8)));
    try std.testing.expectEqual(-0x3.f10984p-8, atanh(@as(f32, -0x3.f0f51cp-8)));
    try std.testing.expectEqual(0x7.7e3f7p-4, atanh(@as(f32, 0x6.fd4ec8p-4)));
    try std.testing.expectEqual(-0x2.7183fcp-4, atanh(@as(f32, -0x2.6cb2a8p-4)));
    try std.testing.expectEqual(-0xf.dfc54p-4, atanh(@as(f32, -0xc.21df7p-4)));
    try std.testing.expectEqual(-0xf.dfc57p-4, atanh(@as(f32, -0xc.21df8p-4)));
    try std.testing.expectEqual(0x5.8be99p-40, atanh(@as(f32, 0x5.8be99p-40)));
    try std.testing.expectEqual(0x3.decf6cp-4, atanh(@as(f32, 0x3.cbed38p-4)));
    try std.testing.expectEqual(0x3.decf68p-4, atanh(@as(f32, 0x3.cbed34p-4)));
    try std.testing.expectEqual(-0x6.068ed8p-4, atanh(@as(f32, -0x5.c18b6p-4)));
    try std.testing.expectEqual(-0x7.c9279p-8, atanh(@as(f32, -0x7.c88a5p-8)));
    try std.testing.expectEqual(-0x2.ce72cp-4, atanh(@as(f32, -0x2.c72b7cp-4)));
    try std.testing.expectEqual(-0x3.a8ec7p-4, atanh(@as(f32, -0x3.98eaf4p-4)));
    try std.testing.expectEqual(0x2.c81f2cp-4, atanh(@as(f32, 0x2.c1085p-4)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1.f8032p-8, atanh(@as(f32, 0x1.f80094p-8)));
    try std.testing.expectEqual(0x2.c73a3cp-4, atanh(@as(f32, 0x2.c02a28p-4)));
    try std.testing.expectEqual(0x2.c73a38p-4, atanh(@as(f32, 0x2.c02a24p-4)));
    try std.testing.expectEqual(0x4p-128, atanh(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-128, atanh(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f32, -0x8p-152)));

    try std.testing.expectEqual(0x0p+0, atanh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0xf.913957192d2b8p-4, atanh(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(-0xf.913957192d2b8p-4, atanh(@as(f64, -0xcp-4)));
    try std.testing.expectEqual(0x4.162bbea045148p-4, atanh(@as(f64, 0x4p-4)));
    try std.testing.expectEqual(0x8.00aac448d771p-8, atanh(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0x4.0000155556224p-12, atanh(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0x1.23450000007dbp-20, atanh(@as(f64, 0x1.2345p-20)));
    try std.testing.expectEqual(0x1.000055558888bp-8, atanh(@as(f64, 0x1p-8)));
    // try std.testing.expectEqual(0x8.0000aaaac4448p-12, atanh(@as(f64, 0x8p-12)));
    try std.testing.expectEqual(0x4.0000155556224p-12, atanh(@as(f64, 0x4p-12)));
    // try std.testing.expectEqual(0x2.000002aaaab12p-12, atanh(@as(f64, 0x2p-12)));
    try std.testing.expectEqual(0x1.0000005555559p-12, atanh(@as(f64, 0x1p-12)));
    try std.testing.expectEqual(0x8.000000aaaaabp-16, atanh(@as(f64, 0x8p-16)));
    try std.testing.expectEqual(0x1.0000000000005p-24, atanh(@as(f64, 0x1p-24)));
    try std.testing.expectEqual(0x8.0000000000008p-28, atanh(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x4p-28, atanh(@as(f64, 0x4p-28)));
    try std.testing.expectEqual(0x2p-28, atanh(@as(f64, 0x2p-28)));
    try std.testing.expectEqual(0x1p-28, atanh(@as(f64, 0x1p-28)));
    try std.testing.expectEqual(0x8p-32, atanh(@as(f64, 0x8p-32)));
    try std.testing.expectEqual(0x4p-32, atanh(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x2p-32, atanh(@as(f64, 0x2p-32)));
    try std.testing.expectEqual(0x1p-32, atanh(@as(f64, 0x1p-32)));
    try std.testing.expectEqual(0x8p-36, atanh(@as(f64, 0x8p-36)));
    try std.testing.expectEqual(0x1p-48, atanh(@as(f64, 0x1p-48)));
    try std.testing.expectEqual(0x8p-52, atanh(@as(f64, 0x8p-52)));
    try std.testing.expectEqual(0x4p-52, atanh(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x2p-52, atanh(@as(f64, 0x2p-52)));
    try std.testing.expectEqual(0x1p-52, atanh(@as(f64, 0x1p-52)));
    try std.testing.expectEqual(0x8p-56, atanh(@as(f64, 0x8p-56)));
    try std.testing.expectEqual(0x4p-56, atanh(@as(f64, 0x4p-56)));
    try std.testing.expectEqual(0x2p-56, atanh(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1p-56, atanh(@as(f64, 0x1p-56)));
    try std.testing.expectEqual(0x8p-60, atanh(@as(f64, 0x8p-60)));
    try std.testing.expectEqual(0x4p-60, atanh(@as(f64, 0x4p-60)));
    try std.testing.expectEqual(0x2p-60, atanh(@as(f64, 0x2p-60)));
    try std.testing.expectEqual(0x1p-100, atanh(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, atanh(@as(f64, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, atanh(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-600, atanh(@as(f64, -0x1p-600)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atanh(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x8.aa122b59bea18p+0, atanh(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea18p+0, atanh(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea18p+0, atanh(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0xb.c8939774cec7p+0, atanh(@as(f64, 0xf.fffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea18p+0, atanh(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.c8939774cec7p+0, atanh(@as(f64, -0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x8.aa122b59bea18p+0, atanh(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.0fb6b4b37945bp+4, atanh(@as(f64, 0xf.fffffffffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea18p+0, atanh(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.0fb6b4b37945bp+4, atanh(@as(f64, -0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea18p+0, atanh(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e2p+4, atanh(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea18p+0, atanh(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e2p+4, atanh(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x8.aa122b59bea18p+0, atanh(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e2p+4, atanh(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea18p+0, atanh(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e2p+4, atanh(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x8.aa122b59bea18p+0, atanh(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e2p+4, atanh(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea18p+0, atanh(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e2p+4, atanh(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x8.aa122b59bea18p+0, atanh(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e2p+4, atanh(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea18p+0, atanh(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e2p+4, atanh(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x6.e6c770006d92cp-20, atanh(@as(f64, -0x6.e6c77p-20)));
    try std.testing.expectEqual(0x3.3794374a456e4p-4, atanh(@as(f64, 0x3.2ca824p-4)));
    try std.testing.expectEqual(-0x1.ce10a18e6cb9fp-4, atanh(@as(f64, -0x1.cc1d66p-4)));
    try std.testing.expectEqual(-0x2.89e0a3edefde6p+0, atanh(@as(f64, -0xf.cd38p-4)));
    // try std.testing.expectEqual(-0x2.89e0cc82cd374p+0, atanh(@as(f64, -0xf.cd381p-4)));
    try std.testing.expectEqual(-0x2.89e0bcc36f11p+0, atanh(@as(f64, -0xf.cd3809ca8fd28p-4)));
    // try std.testing.expectEqual(-0x1.054e200a4c76bp-4, atanh(@as(f64, -0x1.04f386p-4)));
    // try std.testing.expectEqual(-0x2.0b18b5a6aae2cp-4, atanh(@as(f64, -0x2.084568p-4)));
    try std.testing.expectEqual(-0x3.f4cbc2ee0371p-4, atanh(@as(f64, -0x3.e0a5d8p-4)));
    try std.testing.expectEqual(0x3.f3c8af642453p-4, atanh(@as(f64, 0x3.dfb1f8p-4)));
    try std.testing.expectEqual(0x3.f3c8ab2460ea4p-4, atanh(@as(f64, 0x3.dfb1f4p-4)));
    // try std.testing.expectEqual(0x3.f3c8ad1d0289cp-4, atanh(@as(f64, 0x3.dfb1f5db0ceccp-4)));
    try std.testing.expectEqual(0x2.286e7a7dea298p-4, atanh(@as(f64, 0x2.251b2cp-4)));
    try std.testing.expectEqual(0x2.286e766b2cbb6p-4, atanh(@as(f64, 0x2.251b28p-4)));
    try std.testing.expectEqual(0x2.286e78db2bfacp-4, atanh(@as(f64, 0x2.251b2a64c85dep-4)));
    try std.testing.expectEqual(-0x2.eb75aac832c62p-4, atanh(@as(f64, -0x2.e3458cp-4)));
    // try std.testing.expectEqual(0x3.a17be8186229ap-4, atanh(@as(f64, 0x3.91d9f4p-4)));
    try std.testing.expectEqual(0x3.a17be3e2bdc9p-4, atanh(@as(f64, 0x3.91d9fp-4)));
    try std.testing.expectEqual(0x3.a17be7dd80462p-4, atanh(@as(f64, 0x3.91d9f3c80c72ep-4)));
    try std.testing.expectEqual(0x3.a17be7dd8046p-4, atanh(@as(f64, 0x3.91d9f3c80c72cp-4)));
    try std.testing.expectEqual(-0x2.7121d517d0c0cp-4, atanh(@as(f64, -0x2.6c52cp-4)));
    try std.testing.expectEqual(-0x2.7121d92fda686p-4, atanh(@as(f64, -0x2.6c52c4p-4)));
    try std.testing.expectEqual(-0x2.7121d78b9e0d6p-4, atanh(@as(f64, -0x2.6c52c26567198p-4)));
    // try std.testing.expectEqual(0x3.b2f9d9f700e32p-4, atanh(@as(f64, 0x3.a274ecp-4)));
    // try std.testing.expectEqual(-0x3.f10980e9bef52p-8, atanh(@as(f64, -0x3.f0f518p-8)));
    try std.testing.expectEqual(-0x3.f10984e9fd1b2p-8, atanh(@as(f64, -0x3.f0f51cp-8)));
    // try std.testing.expectEqual(-0x3.f109829060504p-8, atanh(@as(f64, -0x3.f0f519a687b64p-8)));
    try std.testing.expectEqual(0x7.7e3f72addbf8cp-4, atanh(@as(f64, 0x6.fd4ec8p-4)));
    try std.testing.expectEqual(-0x2.7183fdca81ffap-4, atanh(@as(f64, -0x2.6cb2a8p-4)));
    // try std.testing.expectEqual(-0xf.dfc543031a8d8p-4, atanh(@as(f64, -0xc.21df7p-4)));
    try std.testing.expectEqual(-0xf.dfc568a8239cp-4, atanh(@as(f64, -0xc.21df8p-4)));
    try std.testing.expectEqual(-0xf.dfc5606a6e958p-4, atanh(@as(f64, -0xc.21df7c7f51508p-4)));
    try std.testing.expectEqual(0x5.8be99p-40, atanh(@as(f64, 0x5.8be99p-40)));
    // try std.testing.expectEqual(0x3.decf6cf9b1c12p-4, atanh(@as(f64, 0x3.cbed38p-4)));
    try std.testing.expectEqual(0x3.decf68bc9915ep-4, atanh(@as(f64, 0x3.cbed34p-4)));
    // try std.testing.expectEqual(0x3.decf6ad980fccp-4, atanh(@as(f64, 0x3.cbed35fe733d8p-4)));
    // try std.testing.expectEqual(-0x6.068ed86859d38p-4, atanh(@as(f64, -0x5.c18b6p-4)));
    try std.testing.expectEqual(-0x7.c92792d39745p-8, atanh(@as(f64, -0x7.c88a5p-8)));
    // try std.testing.expectEqual(-0x2.ce72bf32b10bcp-4, atanh(@as(f64, -0x2.c72b7cp-4)));
    try std.testing.expectEqual(-0x3.a8ec71c4ba57ep-4, atanh(@as(f64, -0x3.98eaf4p-4)));
    try std.testing.expectEqual(0x2.c81f2bf4a730cp-4, atanh(@as(f64, 0x2.c1085p-4)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-500, atanh(@as(f64, 0x1p-500)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1.f8031f3228154p-8, atanh(@as(f64, 0x1.f80094p-8)));
    try std.testing.expectEqual(0x2.c73a3db8f5678p-4, atanh(@as(f64, 0x2.c02a28p-4)));
    // try std.testing.expectEqual(0x2.c73a3999c5d42p-4, atanh(@as(f64, 0x2.c02a24p-4)));
    // try std.testing.expectEqual(0x2.c73a3a9475b5ap-4, atanh(@as(f64, 0x2.c02a24f3472c8p-4)));
    try std.testing.expectEqual(0x2.c73a3a9475b58p-4, atanh(@as(f64, 0x2.c02a24f3472c6p-4)));
    try std.testing.expectEqual(0x4p-128, atanh(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, atanh(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, atanh(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, atanh(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, atanh(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, atanh(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atanh(@as(f64, -0x4p-1076)));

    try std.testing.expectEqual(0x0p+0, atanh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0xf.913957192d2baa3p-4, atanh(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(-0xf.913957192d2baa3p-4, atanh(@as(f80, -0xcp-4)));
    try std.testing.expectEqual(0x4.162bbea0451469c8p-4, atanh(@as(f80, 0x4p-4)));
    try std.testing.expectEqual(0x8.00aac448d77125ap-8, atanh(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x4.0000155556222228p-12, atanh(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x1.23450000007daf66p-20, atanh(@as(f80, 0x1.2345p-20)));
    try std.testing.expectEqual(0x1.000055558888ad1ap-8, atanh(@as(f80, 0x1p-8)));
    try std.testing.expectEqual(0x8.0000aaaac44448dp-12, atanh(@as(f80, 0x8p-12)));
    try std.testing.expectEqual(0x4.0000155556222228p-12, atanh(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x2.000002aaaab1111p-12, atanh(@as(f80, 0x2p-12)));
    try std.testing.expectEqual(0x1.0000005555558888p-12, atanh(@as(f80, 0x1p-12)));
    try std.testing.expectEqual(0x8.000000aaaaaac44p-16, atanh(@as(f80, 0x8p-16)));
    try std.testing.expectEqual(0x1.0000000000005556p-24, atanh(@as(f80, 0x1p-24)));
    try std.testing.expectEqual(0x8.000000000000aabp-28, atanh(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x4.0000000000001558p-28, atanh(@as(f80, 0x4p-28)));
    try std.testing.expectEqual(0x2.00000000000002acp-28, atanh(@as(f80, 0x2p-28)));
    try std.testing.expectEqual(0x1.0000000000000056p-28, atanh(@as(f80, 0x1p-28)));
    try std.testing.expectEqual(0x8.00000000000000bp-32, atanh(@as(f80, 0x8p-32)));
    try std.testing.expectEqual(0x4.0000000000000018p-32, atanh(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x2.0000000000000004p-32, atanh(@as(f80, 0x2p-32)));
    try std.testing.expectEqual(0x1p-32, atanh(@as(f80, 0x1p-32)));
    try std.testing.expectEqual(0x8p-36, atanh(@as(f80, 0x8p-36)));
    try std.testing.expectEqual(0x1p-48, atanh(@as(f80, 0x1p-48)));
    try std.testing.expectEqual(0x8p-52, atanh(@as(f80, 0x8p-52)));
    try std.testing.expectEqual(0x4p-52, atanh(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x2p-52, atanh(@as(f80, 0x2p-52)));
    try std.testing.expectEqual(0x1p-52, atanh(@as(f80, 0x1p-52)));
    try std.testing.expectEqual(0x8p-56, atanh(@as(f80, 0x8p-56)));
    try std.testing.expectEqual(0x4p-56, atanh(@as(f80, 0x4p-56)));
    try std.testing.expectEqual(0x2p-56, atanh(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1p-56, atanh(@as(f80, 0x1p-56)));
    try std.testing.expectEqual(0x8p-60, atanh(@as(f80, 0x8p-60)));
    try std.testing.expectEqual(0x4p-60, atanh(@as(f80, 0x4p-60)));
    try std.testing.expectEqual(0x2p-60, atanh(@as(f80, 0x2p-60)));
    try std.testing.expectEqual(0x1p-100, atanh(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, atanh(@as(f80, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, atanh(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-600, atanh(@as(f80, -0x1p-600)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, atanh(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atanh(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1p-10000, atanh(@as(f80, -0x1p-10000)));
    try std.testing.expectEqual(0x8.aa122b59bea160ep+0, atanh(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160ep+0, atanh(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160ep+0, atanh(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0xb.c8939774cec7147p+0, atanh(@as(f80, 0xf.fffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160ep+0, atanh(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.c8939774cec7147p+0, atanh(@as(f80, -0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160ep+0, atanh(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.0fb6b4b37945ae5p+4, atanh(@as(f80, 0xf.fffffffffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160ep+0, atanh(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.0fb6b4b37945ae5p+4, atanh(@as(f80, -0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160ep+0, atanh(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e1d32p+4, atanh(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160ep+0, atanh(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e1d32p+4, atanh(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160ep+0, atanh(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e1d32p+4, atanh(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.686fc0af622d6f24p+4, atanh(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160ep+0, atanh(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e1d32p+4, atanh(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.686fc0af622d6f24p+4, atanh(@as(f80, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160ep+0, atanh(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e1d32p+4, atanh(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.686fc0af622d6f24p+4, atanh(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160ep+0, atanh(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e1d32p+4, atanh(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.686fc0af622d6f24p+4, atanh(@as(f80, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160ep+0, atanh(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e1d32p+4, atanh(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.686fc0af622d6f24p+4, atanh(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160ep+0, atanh(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e1d32p+4, atanh(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.686fc0af622d6f24p+4, atanh(@as(f80, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x6.e6c770006d92d19p-20, atanh(@as(f80, -0x6.e6c77p-20)));
    try std.testing.expectEqual(0x3.3794374a456e3294p-4, atanh(@as(f80, 0x3.2ca824p-4)));
    try std.testing.expectEqual(-0x1.ce10a18e6cb9ec12p-4, atanh(@as(f80, -0x1.cc1d66p-4)));
    try std.testing.expectEqual(-0x2.89e0a3edefde6854p+0, atanh(@as(f80, -0xf.cd38p-4)));
    try std.testing.expectEqual(-0x2.89e0cc82cd37312cp+0, atanh(@as(f80, -0xf.cd381p-4)));
    try std.testing.expectEqual(-0x2.89e0bcc36f110154p+0, atanh(@as(f80, -0xf.cd3809ca8fd28p-4)));
    try std.testing.expectEqual(-0x1.054e200a4c76aae6p-4, atanh(@as(f80, -0x1.04f386p-4)));
    try std.testing.expectEqual(-0x2.0b18b5a6aae2cfacp-4, atanh(@as(f80, -0x2.084568p-4)));
    try std.testing.expectEqual(-0x3.f4cbc2ee0371007cp-4, atanh(@as(f80, -0x3.e0a5d8p-4)));
    try std.testing.expectEqual(0x3.f3c8af642452fa74p-4, atanh(@as(f80, 0x3.dfb1f8p-4)));
    try std.testing.expectEqual(0x3.f3c8ab2460ea3ef8p-4, atanh(@as(f80, 0x3.dfb1f4p-4)));
    try std.testing.expectEqual(0x3.f3c8ad1d0289cac4p-4, atanh(@as(f80, 0x3.dfb1f5db0ceccp-4)));
    try std.testing.expectEqual(0x2.286e7a7dea2975bcp-4, atanh(@as(f80, 0x2.251b2cp-4)));
    try std.testing.expectEqual(0x2.286e766b2cbb6734p-4, atanh(@as(f80, 0x2.251b28p-4)));
    try std.testing.expectEqual(0x2.286e78db2bfabca4p-4, atanh(@as(f80, 0x2.251b2a64c85dep-4)));
    try std.testing.expectEqual(-0x2.eb75aac832c61fap-4, atanh(@as(f80, -0x2.e3458cp-4)));
    try std.testing.expectEqual(0x3.a17be81862299c04p-4, atanh(@as(f80, 0x3.91d9f4p-4)));
    try std.testing.expectEqual(0x3.a17be3e2bdc8f914p-4, atanh(@as(f80, 0x3.91d9fp-4)));
    try std.testing.expectEqual(0x3.a17be7dd8046218p-4, atanh(@as(f80, 0x3.91d9f3c80c72ep-4)));
    try std.testing.expectEqual(0x3.a17be7dd8045ffd4p-4, atanh(@as(f80, 0x3.91d9f3c80c72cp-4)));
    try std.testing.expectEqual(0x3.a17be7dd804618bcp-4, atanh(@as(f80, 0x3.91d9f3c80c72d7acp-4)));
    try std.testing.expectEqual(-0x2.7121d517d0c0b62cp-4, atanh(@as(f80, -0x2.6c52cp-4)));
    try std.testing.expectEqual(-0x2.7121d92fda685774p-4, atanh(@as(f80, -0x2.6c52c4p-4)));
    try std.testing.expectEqual(-0x2.7121d78b9e0d579p-4, atanh(@as(f80, -0x2.6c52c26567198p-4)));
    try std.testing.expectEqual(0x3.b2f9d9f700e32f28p-4, atanh(@as(f80, 0x3.a274ecp-4)));
    try std.testing.expectEqual(-0x3.f10980e9bef520d4p-8, atanh(@as(f80, -0x3.f0f518p-8)));
    try std.testing.expectEqual(-0x3.f10984e9fd1b129p-8, atanh(@as(f80, -0x3.f0f51cp-8)));
    try std.testing.expectEqual(-0x3.f109829060504074p-8, atanh(@as(f80, -0x3.f0f519a687b64p-8)));
    try std.testing.expectEqual(0x7.7e3f72addbf8ep-4, atanh(@as(f80, 0x6.fd4ec8p-4)));
    try std.testing.expectEqual(-0x2.7183fdca81ffbp-4, atanh(@as(f80, -0x2.6cb2a8p-4)));
    try std.testing.expectEqual(-0xf.dfc543031a8d535p-4, atanh(@as(f80, -0xc.21df7p-4)));
    try std.testing.expectEqual(-0xf.dfc568a8239bd4dp-4, atanh(@as(f80, -0xc.21df8p-4)));
    try std.testing.expectEqual(-0xf.dfc5606a6e957ffp-4, atanh(@as(f80, -0xc.21df7c7f51508p-4)));
    try std.testing.expectEqual(0x5.8be99p-40, atanh(@as(f80, 0x5.8be99p-40)));
    try std.testing.expectEqual(0x3.decf6cf9b1c11f34p-4, atanh(@as(f80, 0x3.cbed38p-4)));
    try std.testing.expectEqual(0x3.decf68bc9915ecc8p-4, atanh(@as(f80, 0x3.cbed34p-4)));
    try std.testing.expectEqual(0x3.decf6ad980fccfecp-4, atanh(@as(f80, 0x3.cbed35fe733d8p-4)));
    try std.testing.expectEqual(-0x6.068ed86859d36f48p-4, atanh(@as(f80, -0x5.c18b6p-4)));
    try std.testing.expectEqual(-0x7.c92792d39744e01p-8, atanh(@as(f80, -0x7.c88a5p-8)));
    try std.testing.expectEqual(-0x2.ce72bf32b10bb258p-4, atanh(@as(f80, -0x2.c72b7cp-4)));
    try std.testing.expectEqual(-0x3.a8ec71c4ba57d65cp-4, atanh(@as(f80, -0x3.98eaf4p-4)));
    try std.testing.expectEqual(0x2.c81f2bf4a730be3cp-4, atanh(@as(f80, 0x2.c1085p-4)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-500, atanh(@as(f80, 0x1p-500)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-5000, atanh(@as(f80, 0x1p-5000)));
    try std.testing.expectEqual(0x1.f8031f3228153f0ep-8, atanh(@as(f80, 0x1.f80094p-8)));
    try std.testing.expectEqual(0x2.c73a3db8f5677cccp-4, atanh(@as(f80, 0x2.c02a28p-4)));
    try std.testing.expectEqual(0x2.c73a3999c5d41404p-4, atanh(@as(f80, 0x2.c02a24p-4)));
    try std.testing.expectEqual(0x2.c73a3a9475b5ae14p-4, atanh(@as(f80, 0x2.c02a24f3472c8p-4)));
    try std.testing.expectEqual(0x2.c73a3a9475b58d1cp-4, atanh(@as(f80, 0x2.c02a24f3472c6p-4)));
    try std.testing.expectEqual(0x2.c73a3a9475b5a61cp-4, atanh(@as(f80, 0x2.c02a24f3472c7844p-4)));
    try std.testing.expectEqual(0x2.c73a3a9475b5a618p-4, atanh(@as(f80, 0x2.c02a24f3472c784p-4)));
    try std.testing.expectEqual(0x4p-128, atanh(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, atanh(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, atanh(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, atanh(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, atanh(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, atanh(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, atanh(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, atanh(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, atanh(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, atanh(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, atanh(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atanh(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, atanh(@as(f80, -0x8p-16448)));

    try std.testing.expectEqual(0x0p+0, atanh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0xf.913957192d2baa37b4a4b67930ep-4, atanh(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(-0xf.913957192d2baa37b4a4b67930ep-4, atanh(@as(f128, -0xcp-4)));
    // try std.testing.expectEqual(0x4.162bbea0451469c9daf0be0810ecp-4, atanh(@as(f128, 0x4p-4)));
    try std.testing.expectEqual(0x8.00aac448d77125a4ee9fee2db378p-8, atanh(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(0x4.000015555622222b46b4dd0dd6bp-12, atanh(@as(f128, 0x4p-12)));
    // try std.testing.expectEqual(0x1.23450000007daf665297209f19c6p-20, atanh(@as(f128, 0x1.2345p-20)));
    try std.testing.expectEqual(0x1.000055558888ad1aee1ef9340408p-8, atanh(@as(f128, 0x1p-8)));
    try std.testing.expectEqual(0x8.0000aaaac44448d68e4c64f4d81p-12, atanh(@as(f128, 0x8p-12)));
    try std.testing.expectEqual(0x4.000015555622222b46b4dd0dd6bp-12, atanh(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x2.000002aaaab11111235a35dc3dc4p-12, atanh(@as(f128, 0x2p-12)));
    try std.testing.expectEqual(0x1.000000555555888888ad1ad1c98dp-12, atanh(@as(f128, 0x1p-12)));
    try std.testing.expectEqual(0x8.000000aaaaaac4444448d68d69b8p-16, atanh(@as(f128, 0x8p-16)));
    try std.testing.expectEqual(0x1.0000000000005555555555558889p-24, atanh(@as(f128, 0x1p-24)));
    try std.testing.expectEqual(0x8.000000000000aaaaaaaaaaaac448p-28, atanh(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x4.0000000000001555555555555624p-28, atanh(@as(f128, 0x4p-28)));
    try std.testing.expectEqual(0x2.00000000000002aaaaaaaaaaaab2p-28, atanh(@as(f128, 0x2p-28)));
    try std.testing.expectEqual(0x1.0000000000000055555555555556p-28, atanh(@as(f128, 0x1p-28)));
    try std.testing.expectEqual(0x8.00000000000000aaaaaaaaaaaaa8p-32, atanh(@as(f128, 0x8p-32)));
    try std.testing.expectEqual(0x4.0000000000000015555555555554p-32, atanh(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0x2.0000000000000002aaaaaaaaaaaap-32, atanh(@as(f128, 0x2p-32)));
    try std.testing.expectEqual(0x1.0000000000000000555555555555p-32, atanh(@as(f128, 0x1p-32)));
    try std.testing.expectEqual(0x8.0000000000000000aaaaaaaaaaa8p-36, atanh(@as(f128, 0x8p-36)));
    try std.testing.expectEqual(0x1.0000000000000000000000005555p-48, atanh(@as(f128, 0x1p-48)));
    try std.testing.expectEqual(0x8.000000000000000000000000aaa8p-52, atanh(@as(f128, 0x8p-52)));
    try std.testing.expectEqual(0x4.0000000000000000000000001554p-52, atanh(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x2.00000000000000000000000002aap-52, atanh(@as(f128, 0x2p-52)));
    try std.testing.expectEqual(0x1.0000000000000000000000000055p-52, atanh(@as(f128, 0x1p-52)));
    try std.testing.expectEqual(0x8.00000000000000000000000000a8p-56, atanh(@as(f128, 0x8p-56)));
    try std.testing.expectEqual(0x4.0000000000000000000000000014p-56, atanh(@as(f128, 0x4p-56)));
    try std.testing.expectEqual(0x2.0000000000000000000000000002p-56, atanh(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1p-56, atanh(@as(f128, 0x1p-56)));
    try std.testing.expectEqual(0x8p-60, atanh(@as(f128, 0x8p-60)));
    try std.testing.expectEqual(0x4p-60, atanh(@as(f128, 0x4p-60)));
    try std.testing.expectEqual(0x2p-60, atanh(@as(f128, 0x2p-60)));
    try std.testing.expectEqual(0x1p-100, atanh(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(-0x1p-100, atanh(@as(f128, -0x1p-100)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, atanh(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1p-600, atanh(@as(f128, -0x1p-600)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, atanh(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(-0x0p+0, atanh(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atanh(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1p-10000, atanh(@as(f128, -0x1p-10000)));
    try std.testing.expectEqual(0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0xb.c8939774cec71468641eed184278p+0, atanh(@as(f128, 0xf.fffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.c8939774cec71468641eed184278p+0, atanh(@as(f128, -0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.0fb6b4b37945ae4f0d24ab00c50cp+4, atanh(@as(f128, 0xf.fffffffffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.0fb6b4b37945ae4f0d24ab00c50cp+4, atanh(@as(f128, -0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e1d31e4b03f1086aap+4, atanh(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e1d31e4b03f1086aap+4, atanh(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e1d31e4b03f1086aap+4, atanh(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.686fc0af622d6f24ee1684ccc806p+4, atanh(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e1d31e4b03f1086aap+4, atanh(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.686fc0af622d6f24ee1684ccc806p+4, atanh(@as(f128, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e1d31e4b03f1086aap+4, atanh(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.686fc0af622d6f24ee1684ccc806p+4, atanh(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x2.51558024a58dbed66b1160844d34p+4, atanh(@as(f128, 0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e1d31e4b03f1086aap+4, atanh(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.686fc0af622d6f24ee1684ccc806p+4, atanh(@as(f128, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x2.51558024a58dbed66b1160844d34p+4, atanh(@as(f128, -0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.2b708872320e1d31e4b03f1086aap+4, atanh(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.686fc0af622d6f24ee1684ccc806p+4, atanh(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x2.78267562db732173ff3b2fcd8e12p+4, atanh(@as(f128, 0xf.fffffffffffffffffffffffffff8p-4)));
    try std.testing.expectEqual(0x2.51558024a58dbed66b1160844d34p+4, atanh(@as(f128, 0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(-0x8.aa122b59bea160e35b98ef96da08p+0, atanh(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(-0x1.2b708872320e1d31e4b03f1086aap+4, atanh(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.686fc0af622d6f24ee1684ccc806p+4, atanh(@as(f128, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x2.78267562db732173ff3b2fcd8e12p+4, atanh(@as(f128, -0xf.fffffffffffffffffffffffffff8p-4)));
    try std.testing.expectEqual(-0x2.51558024a58dbed66b1160844d34p+4, atanh(@as(f128, -0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(-0x6.e6c770006d92d18e1687e22d9ffcp-20, atanh(@as(f128, -0x6.e6c77p-20)));
    // try std.testing.expectEqual(0x3.3794374a456e3292bf5cd3590f7ep-4, atanh(@as(f128, 0x3.2ca824p-4)));
    // try std.testing.expectEqual(-0x1.ce10a18e6cb9ec12c4eddab4daf7p-4, atanh(@as(f128, -0x1.cc1d66p-4)));
    try std.testing.expectEqual(-0x2.89e0a3edefde68544d26767da312p+0, atanh(@as(f128, -0xf.cd38p-4)));
    try std.testing.expectEqual(-0x2.89e0cc82cd37312bdc7f1a1b4f7ap+0, atanh(@as(f128, -0xf.cd381p-4)));
    try std.testing.expectEqual(-0x2.89e0bcc36f110155ec916486f0a4p+0, atanh(@as(f128, -0xf.cd3809ca8fd28p-4)));
    try std.testing.expectEqual(-0x1.054e200a4c76aae62cacc9b0185p-4, atanh(@as(f128, -0x1.04f386p-4)));
    try std.testing.expectEqual(-0x2.0b18b5a6aae2cfad5df39cb56c9ap-4, atanh(@as(f128, -0x2.084568p-4)));
    // try std.testing.expectEqual(-0x3.f4cbc2ee0371007c61ab1041e782p-4, atanh(@as(f128, -0x3.e0a5d8p-4)));
    try std.testing.expectEqual(0x3.f3c8af642452fa7265f1771b5b4cp-4, atanh(@as(f128, 0x3.dfb1f8p-4)));
    // try std.testing.expectEqual(0x3.f3c8ab2460ea3ef88475f9868be8p-4, atanh(@as(f128, 0x3.dfb1f4p-4)));
    // try std.testing.expectEqual(0x3.f3c8ad1d0289cac26a22cccd2efep-4, atanh(@as(f128, 0x3.dfb1f5db0ceccp-4)));
    try std.testing.expectEqual(0x2.286e7a7dea2975bc400c4029191cp-4, atanh(@as(f128, 0x2.251b2cp-4)));
    try std.testing.expectEqual(0x2.286e766b2cbb6735f2df6f49e1b6p-4, atanh(@as(f128, 0x2.251b28p-4)));
    try std.testing.expectEqual(0x2.286e78db2bfabca36df0cd858424p-4, atanh(@as(f128, 0x2.251b2a64c85dep-4)));
    // try std.testing.expectEqual(-0x2.eb75aac832c61fa1080a8277ed6p-4, atanh(@as(f128, -0x2.e3458cp-4)));
    // try std.testing.expectEqual(0x3.a17be81862299c04ac8cc24de7c2p-4, atanh(@as(f128, 0x3.91d9f4p-4)));
    // try std.testing.expectEqual(0x3.a17be3e2bdc8f913609ba2b9621ap-4, atanh(@as(f128, 0x3.91d9fp-4)));
    // try std.testing.expectEqual(0x3.a17be7dd80462181a104c9eaafc6p-4, atanh(@as(f128, 0x3.91d9f3c80c72ep-4)));
    // try std.testing.expectEqual(0x3.a17be7dd8045ffd47dfdefa0422p-4, atanh(@as(f128, 0x3.91d9f3c80c72cp-4)));
    // try std.testing.expectEqual(0x3.a17be7dd804618bdf1c7215b10fap-4, atanh(@as(f128, 0x3.91d9f3c80c72d7acp-4)));
    try std.testing.expectEqual(-0x2.7121d517d0c0b62a7a791d85633cp-4, atanh(@as(f128, -0x2.6c52cp-4)));
    try std.testing.expectEqual(-0x2.7121d92fda685772a132694ae4c8p-4, atanh(@as(f128, -0x2.6c52c4p-4)));
    try std.testing.expectEqual(-0x2.7121d78b9e0d578fbd61de7b57a8p-4, atanh(@as(f128, -0x2.6c52c26567198p-4)));
    // try std.testing.expectEqual(0x3.b2f9d9f700e32f28419a66aa3ee4p-4, atanh(@as(f128, 0x3.a274ecp-4)));
    try std.testing.expectEqual(-0x3.f10980e9bef520d2715b9fa8a22ep-8, atanh(@as(f128, -0x3.f0f518p-8)));
    try std.testing.expectEqual(-0x3.f10984e9fd1b128f333b6e4011eap-8, atanh(@as(f128, -0x3.f0f51cp-8)));
    // try std.testing.expectEqual(-0x3.f109829060504072b047c219061ap-8, atanh(@as(f128, -0x3.f0f519a687b64p-8)));
    try std.testing.expectEqual(0x7.7e3f72addbf8dffe933d8e6e8b78p-4, atanh(@as(f128, 0x6.fd4ec8p-4)));
    // try std.testing.expectEqual(-0x2.7183fdca81ffaffebc101b3793c6p-4, atanh(@as(f128, -0x2.6cb2a8p-4)));
    // try std.testing.expectEqual(-0xf.dfc543031a8d534fa78510d3d478p-4, atanh(@as(f128, -0xc.21df7p-4)));
    try std.testing.expectEqual(-0xf.dfc568a8239bd4d71e81b9dc9e9p-4, atanh(@as(f128, -0xc.21df8p-4)));
    try std.testing.expectEqual(-0xf.dfc5606a6e957febf5ef1a621ebp-4, atanh(@as(f128, -0xc.21df7c7f51508p-4)));
    try std.testing.expectEqual(0x5.8be99000000000000038e0bd45ep-40, atanh(@as(f128, 0x5.8be99p-40)));
    // try std.testing.expectEqual(0x3.decf6cf9b1c11f3526a27331f12cp-4, atanh(@as(f128, 0x3.cbed38p-4)));
    // try std.testing.expectEqual(0x3.decf68bc9915ecc9a2f8c785e97ep-4, atanh(@as(f128, 0x3.cbed34p-4)));
    try std.testing.expectEqual(0x3.decf6ad980fccfedf4ddd9d9005cp-4, atanh(@as(f128, 0x3.cbed35fe733d8p-4)));
    // try std.testing.expectEqual(-0x6.068ed86859d36f45107e2a5fbd28p-4, atanh(@as(f128, -0x5.c18b6p-4)));
    try std.testing.expectEqual(-0x7.c92792d39744e00eaa4f282934p-8, atanh(@as(f128, -0x7.c88a5p-8)));
    // try std.testing.expectEqual(-0x2.ce72bf32b10bb257a11f7f551fdcp-4, atanh(@as(f128, -0x2.c72b7cp-4)));
    try std.testing.expectEqual(-0x3.a8ec71c4ba57d65d8e2679a1418p-4, atanh(@as(f128, -0x3.98eaf4p-4)));
    // try std.testing.expectEqual(0x2.c81f2bf4a730be3ab43f058dee9p-4, atanh(@as(f128, 0x2.c1085p-4)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-500, atanh(@as(f128, 0x1p-500)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atanh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-5000, atanh(@as(f128, 0x1p-5000)));
    // try std.testing.expectEqual(0x1.f8031f3228153f0e56e4db72d2bp-8, atanh(@as(f128, 0x1.f80094p-8)));
    // try std.testing.expectEqual(0x2.c73a3db8f5677ccbe692a02b4a1ep-4, atanh(@as(f128, 0x2.c02a28p-4)));
    // try std.testing.expectEqual(0x2.c73a3999c5d414054dd8dece68d6p-4, atanh(@as(f128, 0x2.c02a24p-4)));
    // try std.testing.expectEqual(0x2.c73a3a9475b5ae1543ad4fac71f4p-4, atanh(@as(f128, 0x2.c02a24f3472c8p-4)));
    // try std.testing.expectEqual(0x2.c73a3a9475b58d1bc712cdaf1658p-4, atanh(@as(f128, 0x2.c02a24f3472c6p-4)));
    // try std.testing.expectEqual(0x2.c73a3a9475b5a61cf6af7781556ep-4, atanh(@as(f128, 0x2.c02a24f3472c7844p-4)));
    // try std.testing.expectEqual(0x2.c73a3a9475b5a618d77fe43115c2p-4, atanh(@as(f128, 0x2.c02a24f3472c784p-4)));
    // try std.testing.expectEqual(0x2.c73a3a9475b5a6198c97986c8f2cp-4, atanh(@as(f128, 0x2.c02a24f3472c7840afbd8cfb68bap-4)));
    // try std.testing.expectEqual(0x2.c73a3a9475b5a6198c97986c8f74p-4, atanh(@as(f128, 0x2.c02a24f3472c7840afbd8cfb69p-4)));
    try std.testing.expectEqual(0x2.c73a3a9475b5a6198c97986c8e6cp-4, atanh(@as(f128, 0x2.c02a24f3472c7840afbd8cfb68p-4)));
    try std.testing.expectEqual(0x4p-128, atanh(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, atanh(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, atanh(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, atanh(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, atanh(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, atanh(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, atanh(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, atanh(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, atanh(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, atanh(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, atanh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, atanh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, atanh(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, atanh(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, atanh(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-152, atanh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atanh(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, atanh(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, atanh(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, atanh(@as(f128, -0x4p-16496)));
}
