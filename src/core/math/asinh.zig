const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const roundeven = @import("roundeven.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn asinh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return asinh(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, asinh32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_asinhf.c
                    return asinh32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_asinh.c
                    return asinh64(x);
                },
                f80 => return cast(f80, asinh128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_asinhl.c
                    return asinh128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn asinh32(x: f32) f32 {
    const ix: [129]f64 = .{
        0x1p+0,           0x1.fc07f01fcp-1, 0x1.f81f81f82p-1, 0x1.f44659e4ap-1,
        0x1.f07c1f07cp-1, 0x1.ecc07b302p-1, 0x1.e9131abfp-1,  0x1.e573ac902p-1,
        0x1.e1e1e1e1ep-1, 0x1.de5d6e3f8p-1, 0x1.dae6076bap-1, 0x1.d77b654b8p-1,
        0x1.d41d41d42p-1, 0x1.d0cb58f6ep-1, 0x1.cd8568904p-1, 0x1.ca4b3055ep-1,
        0x1.c71c71c72p-1, 0x1.c3f8f01c4p-1, 0x1.c0e070382p-1, 0x1.bdd2b8994p-1,
        0x1.bacf914c2p-1, 0x1.b7d6c3ddap-1, 0x1.b4e81b4e8p-1, 0x1.b2036406cp-1,
        0x1.af286bca2p-1, 0x1.ac5701ac6p-1, 0x1.a98ef606ap-1, 0x1.a6d01a6dp-1,
        0x1.a41a41a42p-1, 0x1.a16d3f97ap-1, 0x1.9ec8e951p-1,  0x1.9c2d14ee4p-1,
        0x1.99999999ap-1, 0x1.970e4f80cp-1, 0x1.948b0fcd6p-1, 0x1.920fb49dp-1,
        0x1.8f9c18f9cp-1, 0x1.8d3018d3p-1,  0x1.8acb90f6cp-1, 0x1.886e5f0acp-1,
        0x1.861861862p-1, 0x1.83c977ab2p-1, 0x1.818181818p-1, 0x1.7f405fd02p-1,
        0x1.7d05f417ep-1, 0x1.7ad2208ep-1,  0x1.78a4c8178p-1, 0x1.767dce434p-1,
        0x1.745d1745ep-1, 0x1.724287f46p-1, 0x1.702e05c0cp-1, 0x1.6e1f76b44p-1,
        0x1.6c16c16c2p-1, 0x1.6a13cd154p-1, 0x1.681681682p-1, 0x1.661ec6a52p-1,
        0x1.642c8590cp-1, 0x1.623fa7702p-1, 0x1.605816058p-1, 0x1.5e75bb8dp-1,
        0x1.5c9882b94p-1, 0x1.5ac056b02p-1, 0x1.58ed23082p-1, 0x1.571ed3c5p-1,
        0x1.555555556p-1, 0x1.5390948f4p-1, 0x1.51d07eae2p-1, 0x1.501501502p-1,
        0x1.4e5e0a73p-1,  0x1.4cab88726p-1, 0x1.4afd6a052p-1, 0x1.49539e3b2p-1,
        0x1.47ae147aep-1, 0x1.460cbc7f6p-1, 0x1.446f86562p-1, 0x1.42d6625d6p-1,
        0x1.414141414p-1, 0x1.3fb013fbp-1,  0x1.3e22cbce4p-1, 0x1.3c995a47cp-1,
        0x1.3b13b13b2p-1, 0x1.3991c2c18p-1, 0x1.381381382p-1, 0x1.3698df3dep-1,
        0x1.3521cfb2cp-1, 0x1.33ae45b58p-1, 0x1.323e34a2cp-1, 0x1.30d19013p-1,
        0x1.2f684bda2p-1, 0x1.2e025c04cp-1, 0x1.2c9fb4d82p-1, 0x1.2b404ad02p-1,
        0x1.29e4129e4p-1, 0x1.288b01288p-1, 0x1.27350b882p-1, 0x1.25e22708p-1,
        0x1.24924924ap-1, 0x1.23456789ap-1, 0x1.21fb78122p-1, 0x1.20b470c68p-1,
        0x1.1f7047dc2p-1, 0x1.1e2ef3b4p-1,  0x1.1cf06ada2p-1, 0x1.1bb4a4046p-1,
        0x1.1a7b9611ap-1, 0x1.19453808cp-1, 0x1.181181182p-1, 0x1.16e068942p-1,
        0x1.15b1e5f76p-1, 0x1.1485f0e0ap-1, 0x1.135c81136p-1, 0x1.12358e75ep-1,
        0x1.111111112p-1, 0x1.0fef010fep-1, 0x1.0ecf56be6p-1, 0x1.0db20a89p-1,
        0x1.0c9714fbcp-1, 0x1.0b7e6ec26p-1, 0x1.0a6810a68p-1, 0x1.0953f3902p-1,
        0x1.084210842p-1, 0x1.073260a48p-1, 0x1.0624dd2f2p-1, 0x1.05197f7d8p-1,
        0x1.041041042p-1, 0x1.03091b52p-1,  0x1.020408102p-1, 0x1.01010101p-1,
        0x1p-1,
    };

    const lix: [129]f64 = .{
        0x0p+0,               0x1.fe02a6b146789p-8, 0x1.fc0a8b0fa03e4p-7,
        0x1.7b91b07de311bp-6, 0x1.f829b0e7c33p-6,   0x1.39e87b9fd7d6p-5,
        0x1.77458f63edcfcp-5, 0x1.b42dd7117b1bfp-5, 0x1.f0a30c01362a6p-5,
        0x1.16536eea7fae1p-4, 0x1.341d7961791d1p-4, 0x1.51b073f07983fp-4,
        0x1.6f0d28ae3eb4cp-4, 0x1.8c345d6383b21p-4, 0x1.a926d3a475563p-4,
        0x1.c5e548f63a743p-4, 0x1.e27076e28f2e6p-4, 0x1.fec9131dbaabbp-4,
        0x1.0d77e7ccf6e59p-3, 0x1.1b72ad52f87ap-3,  0x1.29552f81eb523p-3,
        0x1.371fc201f7f74p-3, 0x1.44d2b6ccbfd1ep-3, 0x1.526e5e3a41438p-3,
        0x1.5ff3070a613d4p-3, 0x1.6d60fe717221dp-3, 0x1.7ab890212b909p-3,
        0x1.87fa065214911p-3, 0x1.9525a9cf296b4p-3, 0x1.a23bc1fe42563p-3,
        0x1.af3c94e81bff3p-3, 0x1.bc2867430acd6p-3, 0x1.c8ff7c7989a22p-3,
        0x1.d5c216b535b91p-3, 0x1.e27076e2f92e6p-3, 0x1.ef0adcbe0d936p-3,
        0x1.fb9186d5ebe2bp-3, 0x1.0402594b51041p-2, 0x1.0a324e27370e3p-2,
        0x1.1058bf9ad7ad5p-2, 0x1.1675cabaa660ep-2, 0x1.1c898c16b91fbp-2,
        0x1.22941fbcfb966p-2, 0x1.2895a13dd2ea3p-2, 0x1.2e8e2bade7d31p-2,
        0x1.347dd9a9afd55p-2, 0x1.3a64c556b05eap-2, 0x1.40430868877e4p-2,
        0x1.4618bc219dec2p-2, 0x1.4be5f9579e0a1p-2, 0x1.51aad872c982dp-2,
        0x1.5767717432a6cp-2, 0x1.5d1bdbf5669cap-2, 0x1.62c82f2b83795p-2,
        0x1.686c81e9964afp-2, 0x1.6e08eaa2929e4p-2, 0x1.739d7f6b95007p-2,
        0x1.792a55fdb7fa2p-2, 0x1.7eaf83b82efc3p-2, 0x1.842d1da1ecb17p-2,
        0x1.89a3386be825bp-2, 0x1.8f11e87347ac7p-2, 0x1.947941c1f26fbp-2,
        0x1.99d958119208bp-2, 0x1.9f323ecbd984cp-2, 0x1.a484090e5eb0ap-2,
        0x1.a9cec9a9cf84ap-2, 0x1.af1293245606bp-2, 0x1.b44f77bc98f63p-2,
        0x1.b9858969218fbp-2, 0x1.beb4d9da96b7cp-2, 0x1.c3dd7a7d0354dp-2,
        0x1.c8ff7c79ada22p-2, 0x1.ce1af0b855bebp-2, 0x1.d32fe7e039bd5p-2,
        0x1.d83e72587673ep-2, 0x1.dd46a04c204a1p-2, 0x1.e24881a7cac26p-2,
        0x1.e744261d8a788p-2, 0x1.ec399d2457ccp-2,  0x1.f128f5fac86edp-2,
        0x1.f6123fa71c8acp-2, 0x1.faf588f76631fp-2, 0x1.ffd2e08580c98p-2,
        0x1.02552a5a4f0ffp-1, 0x1.04bdf9da8b6d2p-1, 0x1.0723e5c1b4f4p-1,
        0x1.0986f4f589521p-1, 0x1.0be72e423ca83p-1, 0x1.0e44985d0f48cp-1,
        0x1.109f39e2be497p-1, 0x1.12f71959283bcp-1, 0x1.154c3d2f4f5eap-1,
        0x1.179eabbd9c9a1p-1, 0x1.19ee6b466516fp-1, 0x1.1c3b81f723c25p-1,
        0x1.1e85f5e6ec0dp-1,  0x1.20cdcd193f76ep-1, 0x1.23130d7beb743p-1,
        0x1.2555bce9887cbp-1, 0x1.2795e1288211bp-1, 0x1.29d37fec2308bp-1,
        0x1.2c0e9ed45768cp-1, 0x1.2e47436e5ae68p-1, 0x1.307d7334ff0bep-1,
        0x1.32b1339134571p-1, 0x1.34e289d9b39d3p-1, 0x1.37117b5481bb6p-1,
        0x1.393e0d3549a1ap-1, 0x1.3b6844a017823p-1, 0x1.3d9026a70eefbp-1,
        0x1.3fb5b84cfeb42p-1, 0x1.41d8fe844b2aep-1, 0x1.43f9fe2fb9267p-1,
        0x1.4618bc21d86c2p-1, 0x1.48353d1e928dfp-1, 0x1.4a4f85db1debbp-1,
        0x1.4c679afcc323ap-1, 0x1.4e7d811b77bb1p-1, 0x1.50913cbff8c6bp-1,
        0x1.52a2d265be5abp-1, 0x1.54b2467998498p-1, 0x1.56bf9d5b34b99p-1,
        0x1.58cadb5cbe989p-1, 0x1.5ad404c33af2dp-1, 0x1.5cdb1dc6ad765p-1,
        0x1.5ee02a9241e75p-1, 0x1.60e32f447a8d9p-1, 0x1.62e42fefa39efp-1,
    };

    var t: u32 = @bitCast(x);
    t &= ~@as(u32, 0) >> 1;
    const xs: f64 = cast(f64, x, .{});
    if (t <= 0x3e815667) {
        @branchHint(.unlikely);
        if (t <= 0x39ddb3d7) {
            @branchHint(.unlikely);
            if (t == 0) {
                @branchHint(.unlikely);
                return x;
            }

            return @mulAdd(f32, x, -0x1p-25, x);
        }
        const c: [8]f64 = .{
            0x1.5555555555553p-3, -0x1.3333333330e9dp-4, 0x1.6db6db67cb37ap-5,
            -0x1.f1c71699375dp-6, 0x1.6e8a374c39ff9p-6,  -0x1.1c1e98f9d01e1p-6,
            0x1.c277e96d84026p-7, -0x1.329ff5faf02abp-7,
        };
        const x2: f64 = xs * xs;
        const x4: f64 = x2 * x2;
        const x8: f64 = x4 * x4;
        const f: f64 = x2 * (((c[0] + x2 * c[1]) + x4 * (c[2] + x2 * c[3])) + x8 * ((c[4] + x2 * c[5]) + x4 * (c[6] + x2 * c[7])));
        const r: f64 = xs - xs * f;
        return cast(f32, r, .{});
    } else {
        if (t >= 0x7f800000) {
            @branchHint(.unlikely);
            return x + x; // +-inf or nan
        }

        const xd: f64 = math.abs(xs);
        const x2: f64 = xd * xd;
        const tp: u64 = @bitCast(xd + math.sqrt(x2 + 1));
        const m: u64 = tp & (~@as(u64, 0) >> 12);
        const j: i32 = cast(i32, (m + (1 << (52 - 8))) >> (52 - 7), .{});
        const e: i32 = cast(i32, (tp >> 52) - 0x3ff, .{});
        const w: f64 = @bitCast(m | 0x3ff << 52);
        const z: f64 = w * ix[@intCast(j)] - 1;
        const c: [3]f64 = .{ 0x1.0000000066947p+0, -0x1.00007f053d8cbp-1, 0x1.555280111d914p-2 };
        var z2: f64 = z * z;
        var r: f64 = ((lix[128] * cast(f64, e, .{}) + lix[@intCast(j)]) + z * c[0]) + z2 * (c[1] + z * c[2]);
        if (((@as(u64, @bitCast(r)) + 259000) & 0xfffffff) < 260000) { // accurate path
            @branchHint(.unlikely);
            const cp: [6]f64 = .{ 0x1p+0, -0x1p-1, 0x1.55555555030bcp-2, -0x1.ffffffff2b4e5p-3, 0x1.999b5076a42f2p-3, -0x1.55570c45a647dp-3 };
            z2 = z * z;
            var c0: f64 = cp[0] + z * cp[1];
            const c2: f64 = cp[2] + z * cp[3];
            const c4: f64 = cp[4] + z * cp[5];
            c0 += z2 * (c2 + z2 * c4);
            const ln2l: f64 = 0x1.7f7d1cf79abcap-20;
            const ln2h: f64 = 0x1.62e4p-1;
            const Lh: f64 = ln2h * cast(f64, e, .{});
            const Ll: f64 = ln2l * cast(f64, e, .{});
            r = @mulAdd(f64, z, c0, Ll + lix[@intCast(j)]) + Lh;
            if ((@as(u64, @bitCast(r)) & 0xfffffff) == 0) {
                @branchHint(.unlikely);
                const h: f64 = @mulAdd(f64, z, c0, Ll + lix[@intCast(j)]) + (Lh - r);
                r = r + 64 * h;
            }
        }

        return cast(f32, math.copysign(r, xs), .{});
    }
}

fn asinh64(x: f64) f64 {
    const ln2: f64 = 6.93147180559945286227e-01; // 0x3fe62e42, 0xfefa39ef
    const huge: f64 = 1.00000000000000000000e+300;

    var hx: i32 = undefined;
    dbl64.getHighWord(&hx, x);
    const ix: i32 = hx & 0x7fffffff;
    if (ix < 0x3e300000) { // |x|<2**-28
        @branchHint(.unlikely);
        if (math.abs(x) < std.math.floatMin(f64)) {
            const vx: f64 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }

        if (huge + x > 1)
            return x; // return x inexact except 0
    }

    var w: f64 = undefined;
    if (ix > 0x41b00000) { // |x| > 2**28
        @branchHint(.unlikely);
        if (ix >= 0x7ff00000)
            return x + x; // x is inf or NaN

        w = math.log(math.abs(x)) + ln2;
    } else {
        const xa: f64 = math.abs(x);
        if (ix > 0x40000000) { // 2**28 > |x| > 2.0
            w = math.log(2.0 * xa + 1 / (math.sqrt(xa * xa + 1) + xa));
        } else { // 2.0 > |x| > 2**-28
            const t: f64 = xa * xa;
            w = math.log1p(xa + t / (1 + math.sqrt(1 + t)));
        }
    }

    return math.copysign(w, x);
}

fn asinh128(x: f128) f128 {
    const ln2: f128 = 6.931471805599453094172321214581765681e-1;
    const huge: f128 = 1.0e+4900;

    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const sign: u32 = u.w0;
    const ix: i32 = @bitCast(sign & 0x7fffffff);
    if (ix == 0x7fff0000)
        return x + x; // x is inf or NaN

    if (ix < 0x3fc70000) { // |x| < 2^ -56
        if (math.abs(x) < std.math.floatMin(f128)) {
            const vx: f128 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }

        if (huge + x > 1)
            return x; // return x inexact except 0
    }

    u.w0 = @bitCast(ix);
    var w: f128 = undefined;
    if (ix > 0x40350000) { // |x| > 2 ^ 54
        w = math.log(@as(f128, @bitCast(u))) + ln2;
    } else if (ix > 0x40000000) { // 2^ 54 > |x| > 2.0
        const t: f128 = @bitCast(u);
        w = math.log(2.0 * t + 1 / (math.sqrt(x * x + 1) + t));
    } else { // 2.0 > |x| > 2 ^ -56
        const t: f128 = x * x;
        w = math.log1p(@as(f128, @bitCast(u)) + t / (1 + math.sqrt(1 + t)));
    }

    if ((sign & 0x80000000) != 0) {
        return -w;
    } else {
        return w;
    }
}

test asinh {
    try std.testing.expectEqual(0x0p+0, asinh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, asinh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0xb.17218p-4, asinh(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0xe.1a1b3p-4, asinh(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2.ff8b8cp+0, asinh(@as(f32, 0xap+0)));
    try std.testing.expectEqual(0x5.4c6028p+0, asinh(@as(f32, 0x6.4p+4)));
    try std.testing.expectEqual(0xe.82376p+0, asinh(@as(f32, 0xf.424p+16)));
    try std.testing.expectEqual(0x6.3d0318p+0, asinh(@as(f32, 0x1p+8)));
    try std.testing.expectEqual(0x6.ee75p+0, asinh(@as(f32, 0x2p+8)));
    try std.testing.expectEqual(0x7.9fe708p+0, asinh(@as(f32, 0x4p+8)));
    try std.testing.expectEqual(0x8.51592p+0, asinh(@as(f32, 0x8p+8)));
    try std.testing.expectEqual(0x9.02cb3p+0, asinh(@as(f32, 0x1p+12)));
    try std.testing.expectEqual(0x9.b43d5p+0, asinh(@as(f32, 0x2p+12)));
    try std.testing.expectEqual(0x1.154246p+4, asinh(@as(f32, 0x1p+24)));
    try std.testing.expectEqual(0x1.205966p+4, asinh(@as(f32, 0x2p+24)));
    try std.testing.expectEqual(0x1.2b7088p+4, asinh(@as(f32, 0x4p+24)));
    try std.testing.expectEqual(0x1.3687aap+4, asinh(@as(f32, 0x8p+24)));
    try std.testing.expectEqual(0x1.419eccp+4, asinh(@as(f32, 0x1p+28)));
    try std.testing.expectEqual(0x1.4cb5ecp+4, asinh(@as(f32, 0x2p+28)));
    try std.testing.expectEqual(0x1.57cd0ep+4, asinh(@as(f32, 0x4p+28)));
    try std.testing.expectEqual(0x1.62e43p+4, asinh(@as(f32, 0x8p+28)));
    try std.testing.expectEqual(0x1.6dfb52p+4, asinh(@as(f32, 0x1p+32)));
    try std.testing.expectEqual(0x1.791272p+4, asinh(@as(f32, 0x2p+32)));
    try std.testing.expectEqual(0x2.1f6d68p+4, asinh(@as(f32, 0x1p+48)));
    try std.testing.expectEqual(0x2.2a848cp+4, asinh(@as(f32, 0x2p+48)));
    try std.testing.expectEqual(0x2.359bacp+4, asinh(@as(f32, 0x4p+48)));
    try std.testing.expectEqual(0x2.40b2ccp+4, asinh(@as(f32, 0x8p+48)));
    try std.testing.expectEqual(0x2.4bc9fp+4, asinh(@as(f32, 0x1p+52)));
    try std.testing.expectEqual(0x2.56e11p+4, asinh(@as(f32, 0x2p+52)));
    try std.testing.expectEqual(0x2.61f834p+4, asinh(@as(f32, 0x4p+52)));
    try std.testing.expectEqual(0x2.6d0f54p+4, asinh(@as(f32, 0x8p+52)));
    try std.testing.expectEqual(0x2.782674p+4, asinh(@as(f32, 0x1p+56)));
    try std.testing.expectEqual(0x2.833d98p+4, asinh(@as(f32, 0x2p+56)));
    try std.testing.expectEqual(0x2.8e54b8p+4, asinh(@as(f32, 0x4p+56)));
    try std.testing.expectEqual(0x2.996bd8p+4, asinh(@as(f32, 0x8p+56)));
    try std.testing.expectEqual(0x4.602038p+4, asinh(@as(f32, 0x1p+100)));
    try std.testing.expectEqual(0x5.96a7ep+4, asinh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x5.96a7ep+4, asinh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffd5p-12, asinh(@as(f32, 0x1p-8)));
    try std.testing.expectEqual(0x7.ffffa8p-12, asinh(@as(f32, 0x8p-12)));
    try std.testing.expectEqual(0x3.fffff4p-12, asinh(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffep-12, asinh(@as(f32, 0x2p-12)));
    try std.testing.expectEqual(0x1p-12, asinh(@as(f32, 0x1p-12)));
    try std.testing.expectEqual(0x8p-16, asinh(@as(f32, 0x8p-16)));
    try std.testing.expectEqual(0x1p-24, asinh(@as(f32, 0x1p-24)));
    try std.testing.expectEqual(0x8p-28, asinh(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x4p-28, asinh(@as(f32, 0x4p-28)));
    try std.testing.expectEqual(0x2p-28, asinh(@as(f32, 0x2p-28)));
    try std.testing.expectEqual(0x1p-28, asinh(@as(f32, 0x1p-28)));
    try std.testing.expectEqual(0x8p-32, asinh(@as(f32, 0x8p-32)));
    try std.testing.expectEqual(0x4p-32, asinh(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x2p-32, asinh(@as(f32, 0x2p-32)));
    try std.testing.expectEqual(0x1p-32, asinh(@as(f32, 0x1p-32)));
    try std.testing.expectEqual(0x8p-36, asinh(@as(f32, 0x8p-36)));
    try std.testing.expectEqual(0x1p-48, asinh(@as(f32, 0x1p-48)));
    try std.testing.expectEqual(0x8p-52, asinh(@as(f32, 0x8p-52)));
    try std.testing.expectEqual(0x4p-52, asinh(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x2p-52, asinh(@as(f32, 0x2p-52)));
    try std.testing.expectEqual(0x1p-52, asinh(@as(f32, 0x1p-52)));
    try std.testing.expectEqual(0x8p-56, asinh(@as(f32, 0x8p-56)));
    try std.testing.expectEqual(0x4p-56, asinh(@as(f32, 0x4p-56)));
    try std.testing.expectEqual(0x2p-56, asinh(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1p-56, asinh(@as(f32, 0x1p-56)));
    try std.testing.expectEqual(0x8p-60, asinh(@as(f32, 0x8p-60)));
    try std.testing.expectEqual(0x4p-60, asinh(@as(f32, 0x4p-60)));
    try std.testing.expectEqual(0x2p-60, asinh(@as(f32, 0x2p-60)));
    try std.testing.expectEqual(0x1p-100, asinh(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(-0x3.c958d8p-4, asinh(@as(f32, -0x3.d26bb4p-4)));
    try std.testing.expectEqual(-0x3.b568cp-4, asinh(@as(f32, -0x3.bdeef4p-4)));
    try std.testing.expectEqual(-0x7.fc2afp-8, asinh(@as(f32, -0x7.fc7fc8p-8)));
    try std.testing.expectEqual(-0x3.b0e33p-4, asinh(@as(f32, -0x3.b94a5p-4)));
    try std.testing.expectEqual(-0x3.b0e334p-4, asinh(@as(f32, -0x3.b94a54p-4)));
    try std.testing.expectEqual(0x7.900098p-4, asinh(@as(f32, 0x7.d8e5a8p-4)));
    try std.testing.expectEqual(-0x7.261f58p-4, asinh(@as(f32, -0x7.63a06p-4)));
    try std.testing.expectEqual(-0x7.261f6p-4, asinh(@as(f32, -0x7.63a068p-4)));
    try std.testing.expectEqual(0x6.c0ddep-4, asinh(@as(f32, 0x6.f4a93p-4)));
    try std.testing.expectEqual(-0x7.47c178p-4, asinh(@as(f32, -0x7.88bcc8p-4)));
    try std.testing.expectEqual(-0x3.0d0584p-4, asinh(@as(f32, -0x3.11c35p-4)));
    try std.testing.expectEqual(-0x4.2d24bp-4, asinh(@as(f32, -0x4.39534p-4)));
    try std.testing.expectEqual(-0x4.3170bp+4, asinh(@as(f32, -0xd.d62e8p+92)));
    try std.testing.expectEqual(-0x4.bde0b8p-4, asinh(@as(f32, -0x4.cfb98p-4)));
    try std.testing.expectEqual(-0x4.bde0cp-4, asinh(@as(f32, -0x4.cfb988p-4)));
    try std.testing.expectEqual(-0x5.ac1ebp-4, asinh(@as(f32, -0x5.cabaep-4)));
    try std.testing.expectEqual(-0x5.ac1eb8p-4, asinh(@as(f32, -0x5.cabae8p-4)));
    try std.testing.expectEqual(-0x6.b0186p-4, asinh(@as(f32, -0x6.e26358p-4)));
    try std.testing.expectEqual(0x6.98e81p-4, asinh(@as(f32, 0x6.c92c08p-4)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, asinh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, asinh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0xf.94504p-8, asinh(@as(f32, 0xf.96c69p-8)));
    try std.testing.expectEqual(0x3.fe4e64p-4, asinh(@as(f32, 0x4.08f4p-4)));
    try std.testing.expectEqual(0x3.fe4e5cp-4, asinh(@as(f32, 0x4.08f3f8p-4)));
    try std.testing.expectEqual(-0x5.8cae5p-4, asinh(@as(f32, -0x5.a9568p-4)));
    try std.testing.expectEqual(-0x5.8cae58p-4, asinh(@as(f32, -0x5.a95688p-4)));
    try std.testing.expectEqual(0x4p-128, asinh(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-128, asinh(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, asinh(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x5.96a7ep+4, asinh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x5.96a7ep+4, asinh(@as(f32, -0xf.fffffp+124)));

    try std.testing.expectEqual(0x0p+0, asinh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, asinh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0xb.17217f7d1cf78p-4, asinh(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0xe.1a1b30bcea138p-4, asinh(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x2.ff8b8a0da57b6p+0, asinh(@as(f64, 0xap+0)));
    try std.testing.expectEqual(0x5.4c602a4f4f0a8p+0, asinh(@as(f64, 0x6.4p+4)));
    try std.testing.expectEqual(0xe.823764bfd1e6p+0, asinh(@as(f64, 0xf.424p+16)));
    try std.testing.expectEqual(0x6.3d0317b6484b4p+0, asinh(@as(f64, 0x1p+8)));
    try std.testing.expectEqual(0x6.ee74ffae309acp+0, asinh(@as(f64, 0x2p+8)));
    try std.testing.expectEqual(0x7.9fe70ba603d24p+0, asinh(@as(f64, 0x4p+8)));
    try std.testing.expectEqual(0x8.5159209dd5b8p+0, asinh(@as(f64, 0x8p+8)));
    try std.testing.expectEqual(0x9.02cb37d5a789p+0, asinh(@as(f64, 0x1p+12)));
    try std.testing.expectEqual(0x9.b43d4f9d79588p+0, asinh(@as(f64, 0x2p+12)));
    try std.testing.expectEqual(0x1.1542457337d43p+4, asinh(@as(f64, 0x1p+24)));
    try std.testing.expectEqual(0x1.205966f2b4f12p+4, asinh(@as(f64, 0x2p+24)));
    try std.testing.expectEqual(0x1.2b708872320e2p+4, asinh(@as(f64, 0x4p+24)));
    try std.testing.expectEqual(0x1.3687a9f1af2b1p+4, asinh(@as(f64, 0x8p+24)));
    try std.testing.expectEqual(0x1.419ecb712c481p+4, asinh(@as(f64, 0x1p+28)));
    // try std.testing.expectEqual(0x1.4cb5ecf0a965p+4, asinh(@as(f64, 0x2p+28)));
    // try std.testing.expectEqual(0x1.57cd0e702682p+4, asinh(@as(f64, 0x4p+28)));
    try std.testing.expectEqual(0x1.62e42fefa39efp+4, asinh(@as(f64, 0x8p+28)));
    // try std.testing.expectEqual(0x1.6dfb516f20bbfp+4, asinh(@as(f64, 0x1p+32)));
    try std.testing.expectEqual(0x1.791272ee9dd8ep+4, asinh(@as(f64, 0x2p+32)));
    try std.testing.expectEqual(0x2.1f6d6966f28b6p+4, asinh(@as(f64, 0x1p+48)));
    try std.testing.expectEqual(0x2.2a848ae66fa86p+4, asinh(@as(f64, 0x2p+48)));
    try std.testing.expectEqual(0x2.359bac65ecc56p+4, asinh(@as(f64, 0x4p+48)));
    // try std.testing.expectEqual(0x2.40b2cde569e24p+4, asinh(@as(f64, 0x8p+48)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff4p+4, asinh(@as(f64, 0x1p+52)));
    try std.testing.expectEqual(0x2.56e110e4641c4p+4, asinh(@as(f64, 0x2p+52)));
    try std.testing.expectEqual(0x2.61f83263e1394p+4, asinh(@as(f64, 0x4p+52)));
    // try std.testing.expectEqual(0x2.6d0f53e35e562p+4, asinh(@as(f64, 0x8p+52)));
    try std.testing.expectEqual(0x2.78267562db732p+4, asinh(@as(f64, 0x1p+56)));
    try std.testing.expectEqual(0x2.833d96e258902p+4, asinh(@as(f64, 0x2p+56)));
    try std.testing.expectEqual(0x2.8e54b861d5ad2p+4, asinh(@as(f64, 0x4p+56)));
    // try std.testing.expectEqual(0x2.996bd9e152cap+4, asinh(@as(f64, 0x8p+56)));
    try std.testing.expectEqual(0x4.6020374c5c6dcp+4, asinh(@as(f64, 0x1p+100)));
    try std.testing.expectEqual(0x5.96a7e12e0b98cp+4, asinh(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.5b4448e7fd9b1p+8, asinh(@as(f64, 0x1p+500)));
    try std.testing.expectEqual(0x5.96a7e12e0b98cp+4, asinh(@as(f64, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0x2.c679d1f73f0fcp+8, asinh(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.fffd55568888p-12, asinh(@as(f64, 0x1p-8)));
    try std.testing.expectEqual(0x7.ffffaaaab4444p-12, asinh(@as(f64, 0x8p-12)));
    try std.testing.expectEqual(0x3.fffff55555a22p-12, asinh(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffeaaaaad1p-12, asinh(@as(f64, 0x2p-12)));
    try std.testing.expectEqual(0xf.fffffd5555568p-16, asinh(@as(f64, 0x1p-12)));
    try std.testing.expectEqual(0x7.ffffffaaaaaacp-16, asinh(@as(f64, 0x8p-16)));
    try std.testing.expectEqual(0xf.fffffffffffd8p-28, asinh(@as(f64, 0x1p-24)));
    try std.testing.expectEqual(0x7.ffffffffffffcp-28, asinh(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x4p-28, asinh(@as(f64, 0x4p-28)));
    try std.testing.expectEqual(0x2p-28, asinh(@as(f64, 0x2p-28)));
    try std.testing.expectEqual(0x1p-28, asinh(@as(f64, 0x1p-28)));
    try std.testing.expectEqual(0x8p-32, asinh(@as(f64, 0x8p-32)));
    try std.testing.expectEqual(0x4p-32, asinh(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x2p-32, asinh(@as(f64, 0x2p-32)));
    try std.testing.expectEqual(0x1p-32, asinh(@as(f64, 0x1p-32)));
    try std.testing.expectEqual(0x8p-36, asinh(@as(f64, 0x8p-36)));
    try std.testing.expectEqual(0x1p-48, asinh(@as(f64, 0x1p-48)));
    try std.testing.expectEqual(0x8p-52, asinh(@as(f64, 0x8p-52)));
    try std.testing.expectEqual(0x4p-52, asinh(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x2p-52, asinh(@as(f64, 0x2p-52)));
    try std.testing.expectEqual(0x1p-52, asinh(@as(f64, 0x1p-52)));
    try std.testing.expectEqual(0x8p-56, asinh(@as(f64, 0x8p-56)));
    try std.testing.expectEqual(0x4p-56, asinh(@as(f64, 0x4p-56)));
    try std.testing.expectEqual(0x2p-56, asinh(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1p-56, asinh(@as(f64, 0x1p-56)));
    try std.testing.expectEqual(0x8p-60, asinh(@as(f64, 0x8p-60)));
    try std.testing.expectEqual(0x4p-60, asinh(@as(f64, 0x4p-60)));
    try std.testing.expectEqual(0x2p-60, asinh(@as(f64, 0x2p-60)));
    try std.testing.expectEqual(0x1p-100, asinh(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(-0x3.c958d830129a2p-4, asinh(@as(f64, -0x3.d26bb4p-4)));
    try std.testing.expectEqual(-0x3.b568bf5eec676p-4, asinh(@as(f64, -0x3.bdeef4p-4)));
    // try std.testing.expectEqual(-0x7.fc2aec03ed36p-8, asinh(@as(f64, -0x7.fc7fc8p-8)));
    try std.testing.expectEqual(-0x3.b0e331596c778p-4, asinh(@as(f64, -0x3.b94a5p-4)));
    // try std.testing.expectEqual(-0x3.b0e3353ec4c12p-4, asinh(@as(f64, -0x3.b94a54p-4)));
    try std.testing.expectEqual(-0x3.b0e3342ca9652p-4, asinh(@as(f64, -0x3.b94a52e6913c2p-4)));
    try std.testing.expectEqual(0x7.90009894e809p-4, asinh(@as(f64, 0x7.d8e5a8p-4)));
    try std.testing.expectEqual(-0x7.261f5a1d1207p-4, asinh(@as(f64, -0x7.63a06p-4)));
    // try std.testing.expectEqual(-0x7.261f61605eb6cp-4, asinh(@as(f64, -0x7.63a068p-4)));
    try std.testing.expectEqual(-0x7.261f5cf40e168p-4, asinh(@as(f64, -0x7.63a06320c42e4p-4)));
    // try std.testing.expectEqual(0x6.c0dddeef5ea74p-4, asinh(@as(f64, 0x6.f4a93p-4)));
    try std.testing.expectEqual(-0x7.47c17bbd7ba6p-4, asinh(@as(f64, -0x7.88bcc8p-4)));
    // try std.testing.expectEqual(-0x3.0d05831101b46p-4, asinh(@as(f64, -0x3.11c35p-4)));
    try std.testing.expectEqual(-0x4.2d24ad5bedc88p-4, asinh(@as(f64, -0x4.39534p-4)));
    try std.testing.expectEqual(-0x4.3170acb265858p+4, asinh(@as(f64, -0xd.d62e8p+92)));
    try std.testing.expectEqual(-0x4.bde0b72ea682p-4, asinh(@as(f64, -0x4.cfb98p-4)));
    try std.testing.expectEqual(-0x4.bde0bed7e48ecp-4, asinh(@as(f64, -0x4.cfb988p-4)));
    try std.testing.expectEqual(-0x4.bde0b7852693p-4, asinh(@as(f64, -0x4.cfb9805a53a2p-4)));
    try std.testing.expectEqual(-0x4.bde0b78526934p-4, asinh(@as(f64, -0x4.cfb9805a53a24p-4)));
    try std.testing.expectEqual(-0x5.ac1eaf0870dccp-4, asinh(@as(f64, -0x5.cabaep-4)));
    try std.testing.expectEqual(-0x5.ac1eb68e26b14p-4, asinh(@as(f64, -0x5.cabae8p-4)));
    try std.testing.expectEqual(-0x5.ac1eb633f2fccp-4, asinh(@as(f64, -0x5.cabae7a011e3p-4)));
    try std.testing.expectEqual(-0x5.ac1eb633f2fdp-4, asinh(@as(f64, -0x5.cabae7a011e34p-4)));
    try std.testing.expectEqual(-0x6.b01863558de0cp-4, asinh(@as(f64, -0x6.e26358p-4)));
    try std.testing.expectEqual(0x6.98e810591e8cp-4, asinh(@as(f64, 0x6.c92c08p-4)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, asinh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-500, asinh(@as(f64, 0x1p-500)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, asinh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, asinh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0xf.94503821fcc68p-8, asinh(@as(f64, 0xf.96c69p-8)));
    try std.testing.expectEqual(0x3.fe4e62c525da6p-4, asinh(@as(f64, 0x4.08f4p-4)));
    try std.testing.expectEqual(0x3.fe4e5b035251p-4, asinh(@as(f64, 0x4.08f3f8p-4)));
    // try std.testing.expectEqual(0x3.fe4e5d9acef74p-4, asinh(@as(f64, 0x4.08f3faac4284cp-4)));
    try std.testing.expectEqual(-0x5.8cae501409e88p-4, asinh(@as(f64, -0x5.a9568p-4)));
    try std.testing.expectEqual(-0x5.8cae579ebc7c8p-4, asinh(@as(f64, -0x5.a95688p-4)));
    try std.testing.expectEqual(-0x5.8cae53be0eeb8p-4, asinh(@as(f64, -0x5.a95683e302a7p-4)));
    // try std.testing.expectEqual(-0x5.8cae53be0eebcp-4, asinh(@as(f64, -0x5.a95683e302a74p-4)));
    try std.testing.expectEqual(0x4p-128, asinh(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, asinh(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, asinh(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, asinh(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, asinh(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, asinh(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, asinh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-152, asinh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, asinh(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x5.96a7e12e0b98cp+4, asinh(@as(f64, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0x2.c679d1f73f0fcp+8, asinh(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x5.96a7e12e0b98cp+4, asinh(@as(f64, -0xf.fffffp+124)));
    // try std.testing.expectEqual(-0x2.c679d1f73f0fcp+8, asinh(@as(f64, -0xf.ffffffffffff8p+1020)));

    try std.testing.expectEqual(0x0p+0, asinh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, asinh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0xb.17217f7d1cf79acp-4, asinh(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0xe.1a1b30bcea13661p-4, asinh(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x2.ff8b8a0da57b5aa4p+0, asinh(@as(f80, 0xap+0)));
    try std.testing.expectEqual(0x5.4c602a4f4f0a7cfp+0, asinh(@as(f80, 0x6.4p+4)));
    try std.testing.expectEqual(0xe.823764bfd1e5fa3p+0, asinh(@as(f80, 0xf.424p+16)));
    try std.testing.expectEqual(0x6.3d0317b6484b546p+0, asinh(@as(f80, 0x1p+8)));
    try std.testing.expectEqual(0x6.ee74ffae309ac0e8p+0, asinh(@as(f80, 0x2p+8)));
    try std.testing.expectEqual(0x7.9fe70ba603d23a6p+0, asinh(@as(f80, 0x4p+8)));
    try std.testing.expectEqual(0x8.5159209dd5b8341p+0, asinh(@as(f80, 0x8p+8)));
    try std.testing.expectEqual(0x9.02cb37d5a78915cp+0, asinh(@as(f80, 0x1p+12)));
    try std.testing.expectEqual(0x9.b43d4f9d7958a5ep+0, asinh(@as(f80, 0x2p+12)));
    try std.testing.expectEqual(0x1.1542457337d4321cp+4, asinh(@as(f80, 0x1p+24)));
    try std.testing.expectEqual(0x1.205966f2b4f126b8p+4, asinh(@as(f80, 0x2p+24)));
    try std.testing.expectEqual(0x1.2b708872320e1d92p+4, asinh(@as(f80, 0x4p+24)));
    try std.testing.expectEqual(0x1.3687a9f1af2b14fcp+4, asinh(@as(f80, 0x8p+24)));
    try std.testing.expectEqual(0x1.419ecb712c480c8cp+4, asinh(@as(f80, 0x1p+28)));
    try std.testing.expectEqual(0x1.4cb5ecf0a9650424p+4, asinh(@as(f80, 0x2p+28)));
    try std.testing.expectEqual(0x1.57cd0e702681fbbep+4, asinh(@as(f80, 0x4p+28)));
    try std.testing.expectEqual(0x1.62e42fefa39ef358p+4, asinh(@as(f80, 0x8p+28)));
    try std.testing.expectEqual(0x1.6dfb516f20bbeaf2p+4, asinh(@as(f80, 0x1p+32)));
    try std.testing.expectEqual(0x1.791272ee9dd8e28ep+4, asinh(@as(f80, 0x2p+32)));
    try std.testing.expectEqual(0x2.1f6d6966f28b64ap+4, asinh(@as(f80, 0x1p+48)));
    try std.testing.expectEqual(0x2.2a848ae66fa85c38p+4, asinh(@as(f80, 0x2p+48)));
    try std.testing.expectEqual(0x2.359bac65ecc553d4p+4, asinh(@as(f80, 0x4p+48)));
    try std.testing.expectEqual(0x2.40b2cde569e24b7p+4, asinh(@as(f80, 0x8p+48)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff4308p+4, asinh(@as(f80, 0x1p+52)));
    try std.testing.expectEqual(0x2.56e110e4641c3aa4p+4, asinh(@as(f80, 0x2p+52)));
    try std.testing.expectEqual(0x2.61f83263e139324p+4, asinh(@as(f80, 0x4p+52)));
    try std.testing.expectEqual(0x2.6d0f53e35e5629d8p+4, asinh(@as(f80, 0x8p+52)));
    try std.testing.expectEqual(0x2.78267562db732174p+4, asinh(@as(f80, 0x1p+56)));
    try std.testing.expectEqual(0x2.833d96e25890191p+4, asinh(@as(f80, 0x2p+56)));
    try std.testing.expectEqual(0x2.8e54b861d5ad10a8p+4, asinh(@as(f80, 0x4p+56)));
    try std.testing.expectEqual(0x2.996bd9e152ca0844p+4, asinh(@as(f80, 0x8p+56)));
    try std.testing.expectEqual(0x4.6020374c5c6db01p+4, asinh(@as(f80, 0x1p+100)));
    try std.testing.expectEqual(0x5.96a7e12e0b98bcf8p+4, asinh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.5b4448e7fd9b091ep+8, asinh(@as(f80, 0x1p+500)));
    try std.testing.expectEqual(0x5.96a7e12e0b98bcf8p+4, asinh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb62p+8, asinh(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xd.8a6dd63831ae0fep+8, asinh(@as(f80, 0x1p+5000)));
    try std.testing.expectEqual(0xf.fffd55568887d1bp-12, asinh(@as(f80, 0x1p-8)));
    try std.testing.expectEqual(0x7.ffffaaaab44442d8p-12, asinh(@as(f80, 0x8p-12)));
    try std.testing.expectEqual(0x3.fffff55555a2222p-12, asinh(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffeaaaaad1112p-12, asinh(@as(f80, 0x2p-12)));
    try std.testing.expectEqual(0xf.fffffd555556889p-16, asinh(@as(f80, 0x1p-12)));
    try std.testing.expectEqual(0x7.ffffffaaaaaab448p-16, asinh(@as(f80, 0x8p-16)));
    try std.testing.expectEqual(0xf.fffffffffffd555p-28, asinh(@as(f80, 0x1p-24)));
    try std.testing.expectEqual(0x7.ffffffffffffaaa8p-28, asinh(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x3.fffffffffffff554p-28, asinh(@as(f80, 0x4p-28)));
    try std.testing.expectEqual(0x1.fffffffffffffeaap-28, asinh(@as(f80, 0x2p-28)));
    try std.testing.expectEqual(0xf.fffffffffffffd5p-32, asinh(@as(f80, 0x1p-28)));
    try std.testing.expectEqual(0x7.ffffffffffffffa8p-32, asinh(@as(f80, 0x8p-32)));
    try std.testing.expectEqual(0x3.fffffffffffffff4p-32, asinh(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x1.fffffffffffffffep-32, asinh(@as(f80, 0x2p-32)));
    try std.testing.expectEqual(0x1p-32, asinh(@as(f80, 0x1p-32)));
    try std.testing.expectEqual(0x8p-36, asinh(@as(f80, 0x8p-36)));
    try std.testing.expectEqual(0x1p-48, asinh(@as(f80, 0x1p-48)));
    try std.testing.expectEqual(0x8p-52, asinh(@as(f80, 0x8p-52)));
    try std.testing.expectEqual(0x4p-52, asinh(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x2p-52, asinh(@as(f80, 0x2p-52)));
    try std.testing.expectEqual(0x1p-52, asinh(@as(f80, 0x1p-52)));
    try std.testing.expectEqual(0x8p-56, asinh(@as(f80, 0x8p-56)));
    try std.testing.expectEqual(0x4p-56, asinh(@as(f80, 0x4p-56)));
    try std.testing.expectEqual(0x2p-56, asinh(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1p-56, asinh(@as(f80, 0x1p-56)));
    try std.testing.expectEqual(0x8p-60, asinh(@as(f80, 0x8p-60)));
    try std.testing.expectEqual(0x4p-60, asinh(@as(f80, 0x4p-60)));
    try std.testing.expectEqual(0x2p-60, asinh(@as(f80, 0x2p-60)));
    try std.testing.expectEqual(0x1p-100, asinh(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(-0x3.c958d830129a231p-4, asinh(@as(f80, -0x3.d26bb4p-4)));
    try std.testing.expectEqual(-0x3.b568bf5eec676954p-4, asinh(@as(f80, -0x3.bdeef4p-4)));
    try std.testing.expectEqual(-0x7.fc2aec03ed35ec5p-8, asinh(@as(f80, -0x7.fc7fc8p-8)));
    try std.testing.expectEqual(-0x3.b0e331596c7781ecp-4, asinh(@as(f80, -0x3.b94a5p-4)));
    try std.testing.expectEqual(-0x3.b0e3353ec4c110c4p-4, asinh(@as(f80, -0x3.b94a54p-4)));
    try std.testing.expectEqual(-0x3.b0e3342ca965242cp-4, asinh(@as(f80, -0x3.b94a52e6913c2p-4)));
    try std.testing.expectEqual(0x7.90009894e8091718p-4, asinh(@as(f80, 0x7.d8e5a8p-4)));
    try std.testing.expectEqual(-0x7.261f5a1d1206f028p-4, asinh(@as(f80, -0x7.63a06p-4)));
    try std.testing.expectEqual(-0x7.261f61605eb6bd18p-4, asinh(@as(f80, -0x7.63a068p-4)));
    try std.testing.expectEqual(-0x7.261f5cf40e169a4p-4, asinh(@as(f80, -0x7.63a06320c42e4p-4)));
    try std.testing.expectEqual(0x6.c0dddeef5ea744dp-4, asinh(@as(f80, 0x6.f4a93p-4)));
    try std.testing.expectEqual(-0x7.47c17bbd7ba60748p-4, asinh(@as(f80, -0x7.88bcc8p-4)));
    try std.testing.expectEqual(-0x3.0d05831101b45p-4, asinh(@as(f80, -0x3.11c35p-4)));
    try std.testing.expectEqual(-0x4.2d24ad5bedc89da8p-4, asinh(@as(f80, -0x4.39534p-4)));
    try std.testing.expectEqual(-0x4.3170acb265858p+4, asinh(@as(f80, -0xd.d62e8p+92)));
    try std.testing.expectEqual(-0x4.bde0b72ea681f6e8p-4, asinh(@as(f80, -0x4.cfb98p-4)));
    try std.testing.expectEqual(-0x4.bde0bed7e48ed178p-4, asinh(@as(f80, -0x4.cfb988p-4)));
    try std.testing.expectEqual(-0x4.bde0b78526931428p-4, asinh(@as(f80, -0x4.cfb9805a53a2p-4)));
    try std.testing.expectEqual(-0x4.bde0b7852693517p-4, asinh(@as(f80, -0x4.cfb9805a53a24p-4)));
    try std.testing.expectEqual(-0x4.bde0b78526931a3p-4, asinh(@as(f80, -0x4.cfb9805a53a2065p-4)));
    try std.testing.expectEqual(-0x5.ac1eaf0870dcb5p-4, asinh(@as(f80, -0x5.cabaep-4)));
    try std.testing.expectEqual(-0x5.ac1eb68e26b132a8p-4, asinh(@as(f80, -0x5.cabae8p-4)));
    try std.testing.expectEqual(-0x5.ac1eb633f2fcd81p-4, asinh(@as(f80, -0x5.cabae7a011e3p-4)));
    try std.testing.expectEqual(-0x5.ac1eb633f2fd1438p-4, asinh(@as(f80, -0x5.cabae7a011e34p-4)));
    try std.testing.expectEqual(-0x5.ac1eb633f2fd11fp-4, asinh(@as(f80, -0x5.cabae7a011e33d9p-4)));
    try std.testing.expectEqual(-0x6.b01863558de0abap-4, asinh(@as(f80, -0x6.e26358p-4)));
    try std.testing.expectEqual(0x6.98e810591e8c1c78p-4, asinh(@as(f80, 0x6.c92c08p-4)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, asinh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-500, asinh(@as(f80, 0x1p-500)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, asinh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, asinh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-5000, asinh(@as(f80, 0x1p-5000)));
    try std.testing.expectEqual(0xf.94503821fcc6aebp-8, asinh(@as(f80, 0xf.96c69p-8)));
    try std.testing.expectEqual(0x3.fe4e62c525da587cp-4, asinh(@as(f80, 0x4.08f4p-4)));
    try std.testing.expectEqual(0x3.fe4e5b035250f72p-4, asinh(@as(f80, 0x4.08f3f8p-4)));
    try std.testing.expectEqual(0x3.fe4e5d9acef73cdp-4, asinh(@as(f80, 0x4.08f3faac4284cp-4)));
    try std.testing.expectEqual(-0x5.8cae501409e88378p-4, asinh(@as(f80, -0x5.a9568p-4)));
    try std.testing.expectEqual(-0x5.8cae579ebc7c886p-4, asinh(@as(f80, -0x5.a95688p-4)));
    try std.testing.expectEqual(-0x5.8cae53be0eeb6f38p-4, asinh(@as(f80, -0x5.a95683e302a7p-4)));
    try std.testing.expectEqual(-0x5.8cae53be0eebab9p-4, asinh(@as(f80, -0x5.a95683e302a74p-4)));
    try std.testing.expectEqual(-0x5.8cae53be0eeb8988p-4, asinh(@as(f80, -0x5.a95683e302a71be8p-4)));
    try std.testing.expectEqual(-0x5.8cae53be0eeb899p-4, asinh(@as(f80, -0x5.a95683e302a71bfp-4)));
    try std.testing.expectEqual(0x4p-128, asinh(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, asinh(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, asinh(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, asinh(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, asinh(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, asinh(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, asinh(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, asinh(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, asinh(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, asinh(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, asinh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, asinh(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-152, asinh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, asinh(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, asinh(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x5.96a7e12e0b98bcf8p+4, asinh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb62p+8, asinh(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.c5d37700c6bb03a8p+12, asinh(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x5.96a7e12e0b98bcf8p+4, asinh(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.c679d1f73f0fb62p+8, asinh(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.c5d37700c6bb03a8p+12, asinh(@as(f80, -0xf.fffffffffffffffp+16380)));

    try std.testing.expectEqual(0x0p+0, asinh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, asinh(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0xb.17217f7d1cf79abc9e3b39803f3p-4, asinh(@as(f128, 0xcp-4)));
    // try std.testing.expectEqual(0xe.1a1b30bcea13660d8f99e8dd2518p-4, asinh(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x2.ff8b8a0da57b5aa38395e907170ep+0, asinh(@as(f128, 0xap+0)));
    try std.testing.expectEqual(0x5.4c602a4f4f0a7cedac9045f3d3f8p+0, asinh(@as(f128, 0x6.4p+4)));
    try std.testing.expectEqual(0xe.823764bfd1e5fa37c6bf52ed09dp+0, asinh(@as(f128, 0xf.424p+16)));
    try std.testing.expectEqual(0x6.3d0317b6484b545f6596abfa40dcp+0, asinh(@as(f128, 0x1p+8)));
    // try std.testing.expectEqual(0x6.ee74ffae309ac0eb383199471004p+0, asinh(@as(f128, 0x2p+8)));
    try std.testing.expectEqual(0x7.9fe70ba603d23a62821e041d812cp+0, asinh(@as(f128, 0x4p+8)));
    try std.testing.expectEqual(0x8.5159209dd5b8340d7a01c06cc4b8p+0, asinh(@as(f128, 0x8p+8)));
    try std.testing.expectEqual(0x9.02cb37d5a78915b9409d740d7fe8p+0, asinh(@as(f128, 0x1p+12)));
    try std.testing.expectEqual(0x9.b43d4f9d7958a5e50a7407a58c98p+0, asinh(@as(f128, 0x2p+12)));
    try std.testing.expectEqual(0x1.1542457337d4321c6b73c89d84acp+4, asinh(@as(f128, 0x1p+24)));
    try std.testing.expectEqual(0x1.205966f2b4f126b7281203d70653p+4, asinh(@as(f128, 0x2p+24)));
    try std.testing.expectEqual(0x1.2b708872320e1d91e4b03f1086a9p+4, asinh(@as(f128, 0x4p+24)));
    try std.testing.expectEqual(0x1.3687a9f1af2b14fca14e7a4a06e9p+4, asinh(@as(f128, 0x8p+24)));
    try std.testing.expectEqual(0x1.419ecb712c480c8b5decb5838728p+4, asinh(@as(f128, 0x1p+28)));
    try std.testing.expectEqual(0x1.4cb5ecf0a96504231a8af0bd0768p+4, asinh(@as(f128, 0x2p+28)));
    try std.testing.expectEqual(0x1.57cd0e702681fbbd17292bf687a7p+4, asinh(@as(f128, 0x4p+28)));
    try std.testing.expectEqual(0x1.62e42fefa39ef357a3c7673007e6p+4, asinh(@as(f128, 0x8p+28)));
    try std.testing.expectEqual(0x1.6dfb516f20bbeaf25465a2698825p+4, asinh(@as(f128, 0x1p+32)));
    try std.testing.expectEqual(0x1.791272ee9dd8e28d0e03dda30864p+4, asinh(@as(f128, 0x2p+32)));
    try std.testing.expectEqual(0x2.1f6d6966f28b649e1a4956019018p+4, asinh(@as(f128, 0x1p+48)));
    try std.testing.expectEqual(0x2.2a848ae66fa85c38d6e7913b0d58p+4, asinh(@as(f128, 0x2p+48)));
    try std.testing.expectEqual(0x2.359bac65ecc553d39385cc748cd6p+4, asinh(@as(f128, 0x4p+48)));
    try std.testing.expectEqual(0x2.40b2cde569e24b6e502407ae0ce6p+4, asinh(@as(f128, 0x8p+48)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff43090cc242e78d18p+4, asinh(@as(f128, 0x1p+52)));
    try std.testing.expectEqual(0x2.56e110e4641c3aa3c9607e210d56p+4, asinh(@as(f128, 0x2p+52)));
    try std.testing.expectEqual(0x2.61f83263e139323e85feb95a8d94p+4, asinh(@as(f128, 0x4p+52)));
    // try std.testing.expectEqual(0x2.6d0f53e35e5629d9429cf4940dd2p+4, asinh(@as(f128, 0x8p+52)));
    try std.testing.expectEqual(0x2.78267562db732173ff3b2fcd8e12p+4, asinh(@as(f128, 0x1p+56)));
    // try std.testing.expectEqual(0x2.833d96e25890190ebbd96b070e5p+4, asinh(@as(f128, 0x2p+56)));
    try std.testing.expectEqual(0x2.8e54b861d5ad10a97877a6408e9p+4, asinh(@as(f128, 0x4p+56)));
    try std.testing.expectEqual(0x2.996bd9e152ca08443515e17a0edp+4, asinh(@as(f128, 0x8p+56)));
    // try std.testing.expectEqual(0x4.6020374c5c6db00c6a6d5daf98ecp+4, asinh(@as(f128, 0x1p+100)));
    // try std.testing.expectEqual(0x5.96a7e12e0b98bcf90bb682a4468p+4, asinh(@as(f128, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0x1.5b4448e7fd9b091d321a9e787fbap+8, asinh(@as(f128, 0x1p+500)));
    // try std.testing.expectEqual(0x5.96a7e12e0b98bcf90bb682a4468p+4, asinh(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb620d358b213a7dp+8, asinh(@as(f128, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(0xd.8a6dd63831ae0fdceaf12f64a528p+8, asinh(@as(f128, 0x1p+5000)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb624d358b213a7dp+8, asinh(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    // try std.testing.expectEqual(0xf.fffd55568887d1ad97431894a1dp-12, asinh(@as(f128, 0x1p-8)));
    try std.testing.expectEqual(0x7.ffffaaaab44442d68da70f6582b4p-12, asinh(@as(f128, 0x8p-12)));
    try std.testing.expectEqual(0x3.fffff55555a2221f46b48a6324c4p-12, asinh(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffeaaaaad11110b5a35b2e86ep-12, asinh(@as(f128, 0x2p-12)));
    // try std.testing.expectEqual(0xf.fffffd555556888887d1ad1b4e2p-16, asinh(@as(f128, 0x1p-12)));
    try std.testing.expectEqual(0x7.ffffffaaaaaab4444442d68d6914p-16, asinh(@as(f128, 0x8p-16)));
    try std.testing.expectEqual(0xf.fffffffffffd5555555555568888p-28, asinh(@as(f128, 0x1p-24)));
    try std.testing.expectEqual(0x7.ffffffffffffaaaaaaaaaaaab444p-28, asinh(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x3.fffffffffffff5555555555555a2p-28, asinh(@as(f128, 0x4p-28)));
    try std.testing.expectEqual(0x1.fffffffffffffeaaaaaaaaaaaaadp-28, asinh(@as(f128, 0x2p-28)));
    try std.testing.expectEqual(0xf.fffffffffffffd55555555555558p-32, asinh(@as(f128, 0x1p-28)));
    try std.testing.expectEqual(0x7.ffffffffffffffaaaaaaaaaaaaacp-32, asinh(@as(f128, 0x8p-32)));
    try std.testing.expectEqual(0x3.fffffffffffffff5555555555556p-32, asinh(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0x1.fffffffffffffffeaaaaaaaaaaabp-32, asinh(@as(f128, 0x2p-32)));
    try std.testing.expectEqual(0xf.fffffffffffffffd555555555558p-36, asinh(@as(f128, 0x1p-32)));
    try std.testing.expectEqual(0x7.ffffffffffffffffaaaaaaaaaaacp-36, asinh(@as(f128, 0x8p-36)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffd5558p-52, asinh(@as(f128, 0x1p-48)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffffaaacp-52, asinh(@as(f128, 0x8p-52)));
    try std.testing.expectEqual(0x3.fffffffffffffffffffffffff556p-52, asinh(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x1.fffffffffffffffffffffffffeabp-52, asinh(@as(f128, 0x2p-52)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffd58p-56, asinh(@as(f128, 0x1p-52)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffffffacp-56, asinh(@as(f128, 0x8p-56)));
    try std.testing.expectEqual(0x3.fffffffffffffffffffffffffff6p-56, asinh(@as(f128, 0x4p-56)));
    try std.testing.expectEqual(0x1.ffffffffffffffffffffffffffffp-56, asinh(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1p-56, asinh(@as(f128, 0x1p-56)));
    try std.testing.expectEqual(0x8p-60, asinh(@as(f128, 0x8p-60)));
    try std.testing.expectEqual(0x4p-60, asinh(@as(f128, 0x4p-60)));
    try std.testing.expectEqual(0x2p-60, asinh(@as(f128, 0x2p-60)));
    try std.testing.expectEqual(0x1p-100, asinh(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(-0x3.c958d830129a2311e46d02ce948ap-4, asinh(@as(f128, -0x3.d26bb4p-4)));
    // try std.testing.expectEqual(-0x3.b568bf5eec676953e540354ab56ep-4, asinh(@as(f128, -0x3.bdeef4p-4)));
    // try std.testing.expectEqual(-0x7.fc2aec03ed35ec4f3b852298d344p-8, asinh(@as(f128, -0x7.fc7fc8p-8)));
    // try std.testing.expectEqual(-0x3.b0e331596c7781edc668b0947d86p-4, asinh(@as(f128, -0x3.b94a5p-4)));
    // try std.testing.expectEqual(-0x3.b0e3353ec4c110c2cb860e4335a6p-4, asinh(@as(f128, -0x3.b94a54p-4)));
    try std.testing.expectEqual(-0x3.b0e3342ca965242afb569c3a5ce6p-4, asinh(@as(f128, -0x3.b94a52e6913c2p-4)));
    try std.testing.expectEqual(0x7.90009894e809171b324a20cc7fc8p-4, asinh(@as(f128, 0x7.d8e5a8p-4)));
    try std.testing.expectEqual(-0x7.261f5a1d1206f0273eb68b1daaf4p-4, asinh(@as(f128, -0x7.63a06p-4)));
    // try std.testing.expectEqual(-0x7.261f61605eb6bd156f8f2c73939cp-4, asinh(@as(f128, -0x7.63a068p-4)));
    try std.testing.expectEqual(-0x7.261f5cf40e169a3c2c399a33c774p-4, asinh(@as(f128, -0x7.63a06320c42e4p-4)));
    // try std.testing.expectEqual(0x6.c0dddeef5ea744d14d99f9d11c78p-4, asinh(@as(f128, 0x6.f4a93p-4)));
    // try std.testing.expectEqual(-0x7.47c17bbd7ba607458f4e549f132cp-4, asinh(@as(f128, -0x7.88bcc8p-4)));
    try std.testing.expectEqual(-0x3.0d05831101b4500142e4b2901772p-4, asinh(@as(f128, -0x3.11c35p-4)));
    // try std.testing.expectEqual(-0x4.2d24ad5bedc89dab07914ab2cedcp-4, asinh(@as(f128, -0x4.39534p-4)));
    try std.testing.expectEqual(-0x4.3170acb265858000c5d391e6721p+4, asinh(@as(f128, -0xd.d62e8p+92)));
    // try std.testing.expectEqual(-0x4.bde0b72ea681f6e82ea91bcdc42p-4, asinh(@as(f128, -0x4.cfb98p-4)));
    try std.testing.expectEqual(-0x4.bde0bed7e48ed176770b2cee5404p-4, asinh(@as(f128, -0x4.cfb988p-4)));
    try std.testing.expectEqual(-0x4.bde0b785269314242ca206df8638p-4, asinh(@as(f128, -0x4.cfb9805a53a2p-4)));
    try std.testing.expectEqual(-0x4.bde0b7852693516e1d0cb5a0a548p-4, asinh(@as(f128, -0x4.cfb9805a53a24p-4)));
    try std.testing.expectEqual(-0x4.bde0b78526931a2fb7988c9c128cp-4, asinh(@as(f128, -0x4.cfb9805a53a2065p-4)));
    // try std.testing.expectEqual(-0x5.ac1eaf0870dcb4fc584cd1a4e9a4p-4, asinh(@as(f128, -0x5.cabaep-4)));
    // try std.testing.expectEqual(-0x5.ac1eb68e26b132a5e3ce931aa7b4p-4, asinh(@as(f128, -0x5.cabae8p-4)));
    // try std.testing.expectEqual(-0x5.ac1eb633f2fcd80e04e8e6253e2p-4, asinh(@as(f128, -0x5.cabae7a011e3p-4)));
    // try std.testing.expectEqual(-0x5.ac1eb633f2fd143bb3887641c888p-4, asinh(@as(f128, -0x5.cabae7a011e34p-4)));
    // try std.testing.expectEqual(-0x5.ac1eb633f2fd11f0f621e284b244p-4, asinh(@as(f128, -0x5.cabae7a011e33d9p-4)));
    try std.testing.expectEqual(-0x6.b01863558de0ab9db866832ea558p-4, asinh(@as(f128, -0x6.e26358p-4)));
    try std.testing.expectEqual(0x6.98e810591e8c1c7a088484b273b4p-4, asinh(@as(f128, 0x6.c92c08p-4)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, asinh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-500, asinh(@as(f128, 0x1p-500)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, asinh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, asinh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-5000, asinh(@as(f128, 0x1p-5000)));
    // try std.testing.expectEqual(0xf.94503821fcc6aead1cad28d4935p-8, asinh(@as(f128, 0xf.96c69p-8)));
    try std.testing.expectEqual(0x3.fe4e62c525da587ab7dc0384edeap-4, asinh(@as(f128, 0x4.08f4p-4)));
    // try std.testing.expectEqual(0x3.fe4e5b035250f72070a4b6fdf4c6p-4, asinh(@as(f128, 0x4.08f3f8p-4)));
    try std.testing.expectEqual(0x3.fe4e5d9acef73cd0062bafd4f44p-4, asinh(@as(f128, 0x4.08f3faac4284cp-4)));
    // try std.testing.expectEqual(-0x5.8cae501409e8837aba5edefb984cp-4, asinh(@as(f128, -0x5.a9568p-4)));
    try std.testing.expectEqual(-0x5.8cae579ebc7c885e6164f2f19838p-4, asinh(@as(f128, -0x5.a95688p-4)));
    // try std.testing.expectEqual(-0x5.8cae53be0eeb6f3a31576ce08edcp-4, asinh(@as(f128, -0x5.a95683e302a7p-4)));
    // try std.testing.expectEqual(-0x5.8cae53be0eebab8fc5f7b65e3d1p-4, asinh(@as(f128, -0x5.a95683e302a74p-4)));
    // try std.testing.expectEqual(-0x5.8cae53be0eeb89890245d0ebfbfcp-4, asinh(@as(f128, -0x5.a95683e302a71be8p-4)));
    // try std.testing.expectEqual(-0x5.8cae53be0eeb89908cf864f52bbp-4, asinh(@as(f128, -0x5.a95683e302a71bfp-4)));
    try std.testing.expectEqual(-0x5.8cae53be0eeb89898490af1ca97p-4, asinh(@as(f128, -0x5.a95683e302a71be88a35649b24ep-4)));
    // try std.testing.expectEqual(-0x5.8cae53be0eeb89898490af1ca89cp-4, asinh(@as(f128, -0x5.a95683e302a71be88a35649b24p-4)));
    try std.testing.expectEqual(-0x5.8cae53be0eeb89898490af1caa8p-4, asinh(@as(f128, -0x5.a95683e302a71be88a35649b26p-4)));
    try std.testing.expectEqual(0x4p-128, asinh(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, asinh(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, asinh(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, asinh(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, asinh(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, asinh(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, asinh(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, asinh(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, asinh(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, asinh(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, asinh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, asinh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, asinh(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, asinh(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, asinh(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-152, asinh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, asinh(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, asinh(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, asinh(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, asinh(@as(f128, -0x4p-16496)));
    // try std.testing.expectEqual(0x5.96a7e12e0b98bcf90bb682a4468p+4, asinh(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb620d358b213a7dp+8, asinh(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.c5d37700c6bb03a6c23b6c9b494cp+12, asinh(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2.c5d37700c6bb03a6c24b6c9b494cp+12, asinh(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb624d358b213a7dp+8, asinh(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    // try std.testing.expectEqual(-0x5.96a7e12e0b98bcf90bb682a4468p+4, asinh(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.c679d1f73f0fb620d358b213a7dp+8, asinh(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.c5d37700c6bb03a6c23b6c9b494cp+12, asinh(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x2.c5d37700c6bb03a6c24b6c9b494cp+12, asinh(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x2.c679d1f73f0fb624d358b213a7dp+8, asinh(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
}
