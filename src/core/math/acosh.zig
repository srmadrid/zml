const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const roundeven = @import("roundeven.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn acosh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return acosh(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, acosh32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_acoshf.c
                    return acosh32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_acosh.c
                    return acosh64(x);
                },
                f80 => return cast(f80, acosh128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_acoshl.c
                    return acosh128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn as_special(x: f32) f32 {
    const t: u32 = @bitCast(x);
    if (t == 0x3f800000)
        return 0;

    if ((t << 1) > 0xff000000)
        return x + x; // nan

    if (t == 0x7f800000)
        return x; // inf

    return (x - x) / (x - x);
}

fn acosh32(x: f32) f32 {
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

    const t: u32 = @bitCast(x);
    if (t <= 0x3f800000) {
        @branchHint(.unlikely);
        return as_special(x);
    } else if (t <= 0x3f99db23) { // x <= 0x1.33b646p+0
        @branchHint(.unlikely);
        const zf: f32 = x - 1;
        const z: f64 = cast(f64, zf, .{});
        const a = math.sqrt(2 * z);
        const c: [8]f64 = .{
            -0x1.555555555491ep-4,  0x1.333333319c2p-6,
            -0x1.6db6da26e5e0ep-8,  0x1.f1c698b7100cep-10,
            -0x1.6e7ca9e152dc2p-11, 0x1.1b4cf187158bdp-12,
            -0x1.b5137f11a88e5p-14, 0x1.056c6c57152b2p-15,
        };
        const z2: f64 = z * z;
        const z4: f64 = z2 * z2;
        const f: f64 = ((c[0] + z * c[1]) + z2 * (c[2] + z * c[3])) + z4 * ((c[4] + z * c[5]) + z2 * (c[6] + z * c[7]));
        const r: f64 = a + (a * z) * f;
        return cast(f32, r, .{});
    } else if (t < 0x7f800000) {
        @branchHint(.likely);
        const xd: f64 = cast(f64, x, .{});
        const x2: f64 = xd * xd;
        const tp: u64 = @bitCast(xd + math.sqrt(x2 - 1));
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
            const cp: [6]f64 = .{
                0x1p+0,               -0x1p-1,
                0x1.55555555030bcp-2, -0x1.ffffffff2b4e5p-3,
                0x1.999b5076a42f2p-3, -0x1.55570c45a647dp-3,
            };
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
                r = r + 64.0 * h;
            }
        }

        return cast(f32, r, .{});
    } else {
        return as_special(x);
    }
}

fn acosh64(x: f64) f64 {
    const ln2: f64 = 6.93147180559945286227e-01; // 0x3fe62e42, 0xfefa39ef

    var hx: i64 = undefined;
    dbl64.extractWords64(&hx, x);
    if (hx > 0x4000000000000000) {
        if (hx >= 0x41b0000000000000) {
            @branchHint(.unlikely);
            // x > 2**28
            if (hx >= 0x7ff0000000000000) {
                // x is inf of NaN
                return x + x;
            } else {
                return math.log(x) + ln2; // acosh(huge)=log(2x)
            }
        }

        // 2**28 > x > 2
        const t: f64 = x * x;
        return math.log(2 * x - 1 / (x + math.sqrt(t - 1)));
    } else if (hx > 0x3ff0000000000000) {
        @branchHint(.likely);
        // 1<x<2
        const t: f64 = x - 1;
        return math.log1p(t + math.sqrt(2.0 * t + t * t));
    } else if (hx == 0x3ff0000000000000) {
        @branchHint(.likely);
        return 0; // acosh(1) = 0
    } else { // x < 1
        return (x - x) / (x - x);
    }
}

fn acosh128(x: f128) f128 {
    const ln2: f128 = 0.6931471805599453094172321214581766;

    var lx: u64 = undefined;
    var hx: i64 = undefined;
    ldbl128.getWords(&hx, &lx, x);
    if (hx < 0x3fff000000000000) { // x < 1
        return (x - x) / (x - x);
    } else if (hx >= 0x4035000000000000) { // x > 2**54
        if (hx >= 0x7fff000000000000) { // x is inf of NaN
            return x + x;
        } else return math.log(x) + ln2; // acoshl(huge)=logl(2x)
    } else if (((hx - 0x3fff000000000000) | @as(i64, @bitCast(lx))) == 0) {
        return 0; // acosh(1) = 0
    } else if (hx > 0x4000000000000000) { // 2**28 > x > 2
        const t: f128 = x * x;
        return math.log(2 * x - 1 / (x + math.sqrt(t - 1)));
    } else { // 1<x<2
        const t: f128 = x - 1;
        return math.log1p(t + math.sqrt(2 * t + t * t));
    }
}

test acosh {
    try std.testing.expectEqual(0x0p+0, acosh(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2p-12, acosh(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x2p-12, acosh(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2p-12, acosh(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2p-12, acosh(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2p-12, acosh(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.11156p+0, acosh(@as(f32, 0x1.ap+0)));
    try std.testing.expectEqual(0x2.a2485p+0, acosh(@as(f32, 0x7p+0)));
    try std.testing.expectEqual(0x5.4c5cep+0, acosh(@as(f32, 0x6.4p+4)));
    try std.testing.expectEqual(0xc.34c13p+0, acosh(@as(f32, 0x1.86ap+16)));
    try std.testing.expectEqual(0x6.3d0298p+0, acosh(@as(f32, 0x1p+8)));
    try std.testing.expectEqual(0x6.ee74ep+0, acosh(@as(f32, 0x2p+8)));
    try std.testing.expectEqual(0x7.9fe7p+0, acosh(@as(f32, 0x4p+8)));
    try std.testing.expectEqual(0x8.51592p+0, acosh(@as(f32, 0x8p+8)));
    try std.testing.expectEqual(0x9.02cb3p+0, acosh(@as(f32, 0x1p+12)));
    try std.testing.expectEqual(0x9.b43d5p+0, acosh(@as(f32, 0x2p+12)));
    try std.testing.expectEqual(0x1.154246p+4, acosh(@as(f32, 0x1p+24)));
    try std.testing.expectEqual(0x1.205966p+4, acosh(@as(f32, 0x2p+24)));
    try std.testing.expectEqual(0x1.2b7088p+4, acosh(@as(f32, 0x4p+24)));
    try std.testing.expectEqual(0x1.3687aap+4, acosh(@as(f32, 0x8p+24)));
    try std.testing.expectEqual(0x1.419eccp+4, acosh(@as(f32, 0x1p+28)));
    try std.testing.expectEqual(0x1.4cb5ecp+4, acosh(@as(f32, 0x2p+28)));
    try std.testing.expectEqual(0x1.57cd0ep+4, acosh(@as(f32, 0x4p+28)));
    try std.testing.expectEqual(0x1.62e43p+4, acosh(@as(f32, 0x8p+28)));
    try std.testing.expectEqual(0x1.6dfb52p+4, acosh(@as(f32, 0x1p+32)));
    try std.testing.expectEqual(0x1.791272p+4, acosh(@as(f32, 0x2p+32)));
    try std.testing.expectEqual(0x2.1f6d68p+4, acosh(@as(f32, 0x1p+48)));
    try std.testing.expectEqual(0x2.2a848cp+4, acosh(@as(f32, 0x2p+48)));
    try std.testing.expectEqual(0x2.359bacp+4, acosh(@as(f32, 0x4p+48)));
    try std.testing.expectEqual(0x2.40b2ccp+4, acosh(@as(f32, 0x8p+48)));
    try std.testing.expectEqual(0x2.4bc9fp+4, acosh(@as(f32, 0x1p+52)));
    try std.testing.expectEqual(0x2.56e11p+4, acosh(@as(f32, 0x2p+52)));
    try std.testing.expectEqual(0x2.61f834p+4, acosh(@as(f32, 0x4p+52)));
    try std.testing.expectEqual(0x2.6d0f54p+4, acosh(@as(f32, 0x8p+52)));
    try std.testing.expectEqual(0x2.782674p+4, acosh(@as(f32, 0x1p+56)));
    try std.testing.expectEqual(0x2.833d98p+4, acosh(@as(f32, 0x2p+56)));
    try std.testing.expectEqual(0x2.8e54b8p+4, acosh(@as(f32, 0x4p+56)));
    try std.testing.expectEqual(0x2.996bd8p+4, acosh(@as(f32, 0x8p+56)));
    try std.testing.expectEqual(0x4.602038p+4, acosh(@as(f32, 0x1p+100)));
    try std.testing.expectEqual(0x5.96a7ep+4, acosh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x5.96a7ep+4, acosh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.6f355p-4, acosh(@as(f32, 0x1.80a368p+0)));
    try std.testing.expectEqual(0x4.cd32f8p-4, acosh(@as(f32, 0x1.0b9d4p+0)));
    try std.testing.expectEqual(0x4.cd329p-4, acosh(@as(f32, 0x1.0b9d3ep+0)));
    try std.testing.expectEqual(0x5.f3a5ep-4, acosh(@as(f32, 0x1.11eab6p+0)));
    try std.testing.expectEqual(0x5.a0a868p-4, acosh(@as(f32, 0x1.0fffaap+0)));
    try std.testing.expectEqual(0x3.9cf258p-4, acosh(@as(f32, 0x1.068e1p+0)));
    try std.testing.expectEqual(0x3.9cf1ccp-4, acosh(@as(f32, 0x1.068e0ep+0)));
    try std.testing.expectEqual(0x1.96d14cp+0, acosh(@as(f32, 0x2.8d4288p+0)));
    try std.testing.expectEqual(0x1.9f526p+4, acosh(@as(f32, 0x1.5d71a6p+36)));
    try std.testing.expectEqual(0x5.96a7ep+4, acosh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.38018cp+0, acosh(@as(f32, 0xc.84cb2p+0)));
    try std.testing.expectEqual(0x3.38018cp+0, acosh(@as(f32, 0xc.84cb1p+0)));
    try std.testing.expectEqual(0x3.c4bc24p-4, acosh(@as(f32, 0x1.072238p+0)));
    try std.testing.expectEqual(0x3.c4bb9cp-4, acosh(@as(f32, 0x1.072236p+0)));
    try std.testing.expectEqual(0xe.aaeb3p-4, acosh(@as(f32, 0x1.73515p+0)));
    try std.testing.expectEqual(0x7.5a6488p-4, acosh(@as(f32, 0x1.1b836p+0)));
    try std.testing.expectEqual(0x3.efb2c4p-4, acosh(@as(f32, 0x1.07c956p+0)));
    try std.testing.expectEqual(0x6.e14ea8p-4, acosh(@as(f32, 0x1.1808eep+0)));
    try std.testing.expectEqual(0x5.af0ccp-4, acosh(@as(f32, 0x1.1052c4p+0)));
    try std.testing.expectEqual(0x7.c5802p-4, acosh(@as(f32, 0x1.1ecbep+0)));
    try std.testing.expectEqual(0x7.c57fep-4, acosh(@as(f32, 0x1.1ecbdep+0)));
    try std.testing.expectEqual(0x5.8e60cp-4, acosh(@as(f32, 0x1.0f975ap+0)));
    try std.testing.expectEqual(0x5.8e6068p-4, acosh(@as(f32, 0x1.0f9758p+0)));
    try std.testing.expectEqual(0x5.96a7ep+4, acosh(@as(f32, 0xf.fffffp+124)));

    try std.testing.expectEqual(0x0p+0, acosh(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaadp-12, acosh(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaadp-12, acosh(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f64, 0x1p+0)));
    // try std.testing.expectEqual(0x5.a827999fcef3p-28, acosh(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaadp-12, acosh(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f64, 0x1p+0)));
    // try std.testing.expectEqual(0x5.a827999fcef3p-28, acosh(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaadp-12, acosh(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f64, 0x1p+0)));
    // try std.testing.expectEqual(0x5.a827999fcef3p-28, acosh(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaadp-12, acosh(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f64, 0x1p+0)));
    // try std.testing.expectEqual(0x5.a827999fcef3p-28, acosh(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.1115606c9966p+0, acosh(@as(f64, 0x1.ap+0)));
    try std.testing.expectEqual(0x2.a2484e330086ap+0, acosh(@as(f64, 0x7p+0)));
    // try std.testing.expectEqual(0x5.4c5ce372f189cp+0, acosh(@as(f64, 0x6.4p+4)));
    try std.testing.expectEqual(0xc.34c12d490b8p+0, acosh(@as(f64, 0x1.86ap+16)));
    try std.testing.expectEqual(0x6.3d0297b6484b4p+0, acosh(@as(f64, 0x1p+8)));
    try std.testing.expectEqual(0x6.ee74dfae309acp+0, acosh(@as(f64, 0x2p+8)));
    try std.testing.expectEqual(0x7.9fe703a603d24p+0, acosh(@as(f64, 0x4p+8)));
    try std.testing.expectEqual(0x8.51591e9dd5b8p+0, acosh(@as(f64, 0x8p+8)));
    try std.testing.expectEqual(0x9.02cb3755a789p+0, acosh(@as(f64, 0x1p+12)));
    try std.testing.expectEqual(0x9.b43d4f7d79588p+0, acosh(@as(f64, 0x2p+12)));
    try std.testing.expectEqual(0x1.1542457337d43p+4, acosh(@as(f64, 0x1p+24)));
    try std.testing.expectEqual(0x1.205966f2b4f12p+4, acosh(@as(f64, 0x2p+24)));
    try std.testing.expectEqual(0x1.2b708872320e2p+4, acosh(@as(f64, 0x4p+24)));
    try std.testing.expectEqual(0x1.3687a9f1af2b1p+4, acosh(@as(f64, 0x8p+24)));
    // try std.testing.expectEqual(0x1.419ecb712c481p+4, acosh(@as(f64, 0x1p+28)));
    try std.testing.expectEqual(0x1.4cb5ecf0a965p+4, acosh(@as(f64, 0x2p+28)));
    // try std.testing.expectEqual(0x1.57cd0e702682p+4, acosh(@as(f64, 0x4p+28)));
    try std.testing.expectEqual(0x1.62e42fefa39efp+4, acosh(@as(f64, 0x8p+28)));
    // try std.testing.expectEqual(0x1.6dfb516f20bbfp+4, acosh(@as(f64, 0x1p+32)));
    try std.testing.expectEqual(0x1.791272ee9dd8ep+4, acosh(@as(f64, 0x2p+32)));
    try std.testing.expectEqual(0x2.1f6d6966f28b6p+4, acosh(@as(f64, 0x1p+48)));
    try std.testing.expectEqual(0x2.2a848ae66fa86p+4, acosh(@as(f64, 0x2p+48)));
    try std.testing.expectEqual(0x2.359bac65ecc56p+4, acosh(@as(f64, 0x4p+48)));
    // try std.testing.expectEqual(0x2.40b2cde569e24p+4, acosh(@as(f64, 0x8p+48)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff4p+4, acosh(@as(f64, 0x1p+52)));
    try std.testing.expectEqual(0x2.56e110e4641c4p+4, acosh(@as(f64, 0x2p+52)));
    try std.testing.expectEqual(0x2.61f83263e1394p+4, acosh(@as(f64, 0x4p+52)));
    // try std.testing.expectEqual(0x2.6d0f53e35e562p+4, acosh(@as(f64, 0x8p+52)));
    try std.testing.expectEqual(0x2.78267562db732p+4, acosh(@as(f64, 0x1p+56)));
    try std.testing.expectEqual(0x2.833d96e258902p+4, acosh(@as(f64, 0x2p+56)));
    try std.testing.expectEqual(0x2.8e54b861d5ad2p+4, acosh(@as(f64, 0x4p+56)));
    // try std.testing.expectEqual(0x2.996bd9e152cap+4, acosh(@as(f64, 0x8p+56)));
    try std.testing.expectEqual(0x4.6020374c5c6dcp+4, acosh(@as(f64, 0x1p+100)));
    try std.testing.expectEqual(0x5.96a7e12e0b98cp+4, acosh(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.5b4448e7fd9b1p+8, acosh(@as(f64, 0x1p+500)));
    try std.testing.expectEqual(0x5.96a7e12e0b98cp+4, acosh(@as(f64, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0x2.c679d1f73f0fcp+8, acosh(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.6f3554d09745p-4, acosh(@as(f64, 0x1.80a368p+0)));
    try std.testing.expectEqual(0x4.cd32fbd26f204p-4, acosh(@as(f64, 0x1.0b9d4p+0)));
    try std.testing.expectEqual(0x4.cd3292c615cdp-4, acosh(@as(f64, 0x1.0b9d3ep+0)));
    try std.testing.expectEqual(0x4.cd32b38d26738p-4, acosh(@as(f64, 0x1.0b9d3e9fc19fbp+0)));
    try std.testing.expectEqual(0x5.f3a5dc5f386cp-4, acosh(@as(f64, 0x1.11eab6p+0)));
    try std.testing.expectEqual(0x5.a0a86a1d26c8p-4, acosh(@as(f64, 0x1.0fffaap+0)));
    // try std.testing.expectEqual(0x3.9cf2578fe1276p-4, acosh(@as(f64, 0x1.068e1p+0)));
    try std.testing.expectEqual(0x3.9cf1cb0d2c9b6p-4, acosh(@as(f64, 0x1.068e0ep+0)));
    // try std.testing.expectEqual(0x3.9cf202813d516p-4, acosh(@as(f64, 0x1.068e0eca105a6p+0)));
    try std.testing.expectEqual(0x1.96d14c020ca3cp+0, acosh(@as(f64, 0x2.8d4288p+0)));
    // try std.testing.expectEqual(0x1.9f5260500df79p+4, acosh(@as(f64, 0x1.5d71a6p+36)));
    try std.testing.expectEqual(0x5.96a7e12e0b98cp+4, acosh(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.890ca4134d0a8p+4, acosh(@as(f64, 0x1.476a3c43d7edep+288)));
    try std.testing.expectEqual(0x3.38018c2243818p+0, acosh(@as(f64, 0xc.84cb2p+0)));
    try std.testing.expectEqual(0x3.38018ada065fcp+0, acosh(@as(f64, 0xc.84cb1p+0)));
    // try std.testing.expectEqual(0x3.38018bf3c523ap+0, acosh(@as(f64, 0xc.84cb1dbbd1b1p+0)));
    try std.testing.expectEqual(0x3.c4bc22bfc0f1p-4, acosh(@as(f64, 0x1.072238p+0)));
    // try std.testing.expectEqual(0x3.c4bb9c223d182p-4, acosh(@as(f64, 0x1.072236p+0)));
    try std.testing.expectEqual(0x3.c4bba7bdf5dacp-4, acosh(@as(f64, 0x1.0722362c26ba5p+0)));
    try std.testing.expectEqual(0xe.aaeb3064737c8p-4, acosh(@as(f64, 0x1.73515p+0)));
    // try std.testing.expectEqual(0x7.5a648ae540c28p-4, acosh(@as(f64, 0x1.1b836p+0)));
    // try std.testing.expectEqual(0x3.efb2c425f3ffcp-4, acosh(@as(f64, 0x1.07c956p+0)));
    try std.testing.expectEqual(0x6.e14ea65707588p-4, acosh(@as(f64, 0x1.1808eep+0)));
    // try std.testing.expectEqual(0x5.af0cbc66dcd48p-4, acosh(@as(f64, 0x1.1052c4p+0)));
    try std.testing.expectEqual(0x7.c580228bed2a8p-4, acosh(@as(f64, 0x1.1ecbep+0)));
    try std.testing.expectEqual(0x7.c57fe32f85bf4p-4, acosh(@as(f64, 0x1.1ecbdep+0)));
    try std.testing.expectEqual(0x7.c58009b5885c8p-4, acosh(@as(f64, 0x1.1ecbdf374bce1p+0)));
    // try std.testing.expectEqual(0x7.c58009b5885a8p-4, acosh(@as(f64, 0x1.1ecbdf374bcep+0)));
    try std.testing.expectEqual(0x5.8e60c2aef6b84p-4, acosh(@as(f64, 0x1.0f975ap+0)));
    try std.testing.expectEqual(0x5.8e60685c1d694p-4, acosh(@as(f64, 0x1.0f9758p+0)));
    // try std.testing.expectEqual(0x5.8e607be4bc4ccp-4, acosh(@as(f64, 0x1.0f97586eba091p+0)));
    try std.testing.expectEqual(0x5.8e607be4bc4ap-4, acosh(@as(f64, 0x1.0f97586eba09p+0)));
    try std.testing.expectEqual(0x5.96a7e12e0b98cp+4, acosh(@as(f64, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0x2.c679d1f73f0fcp+8, acosh(@as(f64, 0xf.ffffffffffff8p+1020)));

    try std.testing.expectEqual(0x0p+0, acosh(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad112p-12, acosh(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad112p-12, acosh(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x5.a827999fcef31c98p-28, acosh(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad112p-12, acosh(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x5.a827999fcef31c98p-28, acosh(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2p-32, acosh(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad112p-12, acosh(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x5.a827999fcef31c98p-28, acosh(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2p-32, acosh(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad112p-12, acosh(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x5.a827999fcef31c98p-28, acosh(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2p-32, acosh(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x1.1115606c996603fcp+0, acosh(@as(f80, 0x1.ap+0)));
    try std.testing.expectEqual(0x2.a2484e330086937cp+0, acosh(@as(f80, 0x7p+0)));
    try std.testing.expectEqual(0x5.4c5ce372f189a288p+0, acosh(@as(f80, 0x6.4p+4)));
    try std.testing.expectEqual(0xc.34c12d490b800bap+0, acosh(@as(f80, 0x1.86ap+16)));
    try std.testing.expectEqual(0x6.3d0297b6484b39b8p+0, acosh(@as(f80, 0x1p+8)));
    try std.testing.expectEqual(0x6.ee74dfae309ac08p+0, acosh(@as(f80, 0x2p+8)));
    try std.testing.expectEqual(0x7.9fe703a603d23a6p+0, acosh(@as(f80, 0x4p+8)));
    try std.testing.expectEqual(0x8.51591e9dd5b8341p+0, acosh(@as(f80, 0x8p+8)));
    try std.testing.expectEqual(0x9.02cb3755a78915cp+0, acosh(@as(f80, 0x1p+12)));
    try std.testing.expectEqual(0x9.b43d4f7d7958a5ep+0, acosh(@as(f80, 0x2p+12)));
    try std.testing.expectEqual(0x1.1542457337d42a1cp+4, acosh(@as(f80, 0x1p+24)));
    try std.testing.expectEqual(0x1.205966f2b4f124b8p+4, acosh(@as(f80, 0x2p+24)));
    try std.testing.expectEqual(0x1.2b708872320e1d12p+4, acosh(@as(f80, 0x4p+24)));
    try std.testing.expectEqual(0x1.3687a9f1af2b14dcp+4, acosh(@as(f80, 0x8p+24)));
    try std.testing.expectEqual(0x1.419ecb712c480c84p+4, acosh(@as(f80, 0x1p+28)));
    try std.testing.expectEqual(0x1.4cb5ecf0a9650422p+4, acosh(@as(f80, 0x2p+28)));
    try std.testing.expectEqual(0x1.57cd0e702681fbbcp+4, acosh(@as(f80, 0x4p+28)));
    try std.testing.expectEqual(0x1.62e42fefa39ef358p+4, acosh(@as(f80, 0x8p+28)));
    try std.testing.expectEqual(0x1.6dfb516f20bbeaf2p+4, acosh(@as(f80, 0x1p+32)));
    try std.testing.expectEqual(0x1.791272ee9dd8e28ep+4, acosh(@as(f80, 0x2p+32)));
    try std.testing.expectEqual(0x2.1f6d6966f28b64ap+4, acosh(@as(f80, 0x1p+48)));
    try std.testing.expectEqual(0x2.2a848ae66fa85c38p+4, acosh(@as(f80, 0x2p+48)));
    try std.testing.expectEqual(0x2.359bac65ecc553d4p+4, acosh(@as(f80, 0x4p+48)));
    try std.testing.expectEqual(0x2.40b2cde569e24b7p+4, acosh(@as(f80, 0x8p+48)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff4308p+4, acosh(@as(f80, 0x1p+52)));
    try std.testing.expectEqual(0x2.56e110e4641c3aa4p+4, acosh(@as(f80, 0x2p+52)));
    try std.testing.expectEqual(0x2.61f83263e139324p+4, acosh(@as(f80, 0x4p+52)));
    try std.testing.expectEqual(0x2.6d0f53e35e5629d8p+4, acosh(@as(f80, 0x8p+52)));
    try std.testing.expectEqual(0x2.78267562db732174p+4, acosh(@as(f80, 0x1p+56)));
    try std.testing.expectEqual(0x2.833d96e25890191p+4, acosh(@as(f80, 0x2p+56)));
    try std.testing.expectEqual(0x2.8e54b861d5ad10a8p+4, acosh(@as(f80, 0x4p+56)));
    try std.testing.expectEqual(0x2.996bd9e152ca0844p+4, acosh(@as(f80, 0x8p+56)));
    try std.testing.expectEqual(0x4.6020374c5c6db01p+4, acosh(@as(f80, 0x1p+100)));
    try std.testing.expectEqual(0x5.96a7e12e0b98bcf8p+4, acosh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.5b4448e7fd9b091ep+8, acosh(@as(f80, 0x1p+500)));
    try std.testing.expectEqual(0x5.96a7e12e0b98bcf8p+4, acosh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb62p+8, acosh(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xd.8a6dd63831ae0fep+8, acosh(@as(f80, 0x1p+5000)));
    try std.testing.expectEqual(0xf.6f3554d097450e3p-4, acosh(@as(f80, 0x1.80a368p+0)));
    try std.testing.expectEqual(0x4.cd32fbd26f204f38p-4, acosh(@as(f80, 0x1.0b9d4p+0)));
    try std.testing.expectEqual(0x4.cd3292c615ccfa7p-4, acosh(@as(f80, 0x1.0b9d3ep+0)));
    try std.testing.expectEqual(0x4.cd32b38d267382a8p-4, acosh(@as(f80, 0x1.0b9d3e9fc19fbp+0)));
    try std.testing.expectEqual(0x5.f3a5dc5f386bf398p-4, acosh(@as(f80, 0x1.11eab6p+0)));
    try std.testing.expectEqual(0x5.a0a86a1d26c7fdf8p-4, acosh(@as(f80, 0x1.0fffaap+0)));
    try std.testing.expectEqual(0x3.9cf2578fe127545p-4, acosh(@as(f80, 0x1.068e1p+0)));
    try std.testing.expectEqual(0x3.9cf1cb0d2c9b6d24p-4, acosh(@as(f80, 0x1.068e0ep+0)));
    try std.testing.expectEqual(0x3.9cf202813d515f4p-4, acosh(@as(f80, 0x1.068e0eca105a6p+0)));
    try std.testing.expectEqual(0x1.96d14c020ca3cp+0, acosh(@as(f80, 0x2.8d4288p+0)));
    try std.testing.expectEqual(0x1.9f5260500df788p+4, acosh(@as(f80, 0x1.5d71a6p+36)));
    try std.testing.expectEqual(0x5.96a7e12e0b98bcf8p+4, acosh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.890ca4134d0a8p+4, acosh(@as(f80, 0x1.476a3c43d7edep+288)));
    try std.testing.expectEqual(0x3.38018c2243818418p+0, acosh(@as(f80, 0xc.84cb2p+0)));
    try std.testing.expectEqual(0x3.38018ada065fb6dp+0, acosh(@as(f80, 0xc.84cb1p+0)));
    try std.testing.expectEqual(0x3.38018bf3c5239p+0, acosh(@as(f80, 0xc.84cb1dbbd1b1p+0)));
    try std.testing.expectEqual(0x3.c4bc22bfc0f0f174p-4, acosh(@as(f80, 0x1.072238p+0)));
    try std.testing.expectEqual(0x3.c4bb9c223d182318p-4, acosh(@as(f80, 0x1.072236p+0)));
    try std.testing.expectEqual(0x3.c4bba7bdf5dac0d4p-4, acosh(@as(f80, 0x1.0722362c26ba5p+0)));
    try std.testing.expectEqual(0xe.aaeb3064737cacep-4, acosh(@as(f80, 0x1.73515p+0)));
    try std.testing.expectEqual(0x7.5a648ae540c2606p-4, acosh(@as(f80, 0x1.1b836p+0)));
    try std.testing.expectEqual(0x3.efb2c425f3ffc308p-4, acosh(@as(f80, 0x1.07c956p+0)));
    try std.testing.expectEqual(0x6.e14ea6570758831p-4, acosh(@as(f80, 0x1.1808eep+0)));
    try std.testing.expectEqual(0x5.af0cbc66dcd4993p-4, acosh(@as(f80, 0x1.1052c4p+0)));
    try std.testing.expectEqual(0x7.c580228bed2a829p-4, acosh(@as(f80, 0x1.1ecbep+0)));
    try std.testing.expectEqual(0x7.c57fe32f85bf5c48p-4, acosh(@as(f80, 0x1.1ecbdep+0)));
    try std.testing.expectEqual(0x7.c58009b5885c946p-4, acosh(@as(f80, 0x1.1ecbdf374bce1p+0)));
    try std.testing.expectEqual(0x7.c58009b5885a9978p-4, acosh(@as(f80, 0x1.1ecbdf374bcep+0)));
    try std.testing.expectEqual(0x7.c58009b5885ad228p-4, acosh(@as(f80, 0x1.1ecbdf374bce01cap+0)));
    try std.testing.expectEqual(0x5.8e60c2aef6b85ap-4, acosh(@as(f80, 0x1.0f975ap+0)));
    try std.testing.expectEqual(0x5.8e60685c1d693548p-4, acosh(@as(f80, 0x1.0f9758p+0)));
    try std.testing.expectEqual(0x5.8e607be4bc4cd598p-4, acosh(@as(f80, 0x1.0f97586eba091p+0)));
    try std.testing.expectEqual(0x5.8e607be4bc4a03p-4, acosh(@as(f80, 0x1.0f97586eba09p+0)));
    try std.testing.expectEqual(0x5.8e607be4bc4a5dbp-4, acosh(@as(f80, 0x1.0f97586eba090202p+0)));
    try std.testing.expectEqual(0x5.8e607be4bc4a5d58p-4, acosh(@as(f80, 0x1.0f97586eba0902p+0)));
    try std.testing.expectEqual(0x5.96a7e12e0b98bcf8p+4, acosh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb62p+8, acosh(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.c5d37700c6bb03a8p+12, acosh(@as(f80, 0xf.fffffffffffffffp+16380)));

    try std.testing.expectEqual(0x0p+0, acosh(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad11110fa35a369c4p-12, acosh(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad11110fa35a369c4p-12, acosh(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x5.a827999fcef31c97ec74cddd10ecp-28, acosh(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad11110fa35a369c4p-12, acosh(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x5.a827999fcef31c97ec74cddd10ecp-28, acosh(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.ffffffffffffffffaaaaaaaaaaabp-32, acosh(@as(f128, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad11110fa35a369c4p-12, acosh(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x5.a827999fcef31c97ec74cddd10ecp-28, acosh(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.ffffffffffffffffaaaaaaaaaaabp-32, acosh(@as(f128, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffff58p-56, acosh(@as(f128, 0x1.000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x1.ffffffaaaaaad11110fa35a369c4p-12, acosh(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, acosh(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x5.a827999fcef31c97ec74cddd10ecp-28, acosh(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.ffffffffffffffffaaaaaaaaaaabp-32, acosh(@as(f128, 0x1.0000000000000002p+0)));
    // try std.testing.expectEqual(0x1.6a09e667f3bcc908b2fb1366ea95p-56, acosh(@as(f128, 0x1.0000000000000000000000000001p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffff58p-56, acosh(@as(f128, 0x1.000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x1.1115606c996603fcd1822ec9d026p+0, acosh(@as(f128, 0x1.ap+0)));
    try std.testing.expectEqual(0x2.a2484e330086937cd097a0311aa6p+0, acosh(@as(f128, 0x7p+0)));
    try std.testing.expectEqual(0x5.4c5ce372f189a2883e416c22565cp+0, acosh(@as(f128, 0x6.4p+4)));
    try std.testing.expectEqual(0xc.34c12d490b800b9c3c8156140948p+0, acosh(@as(f128, 0x1.86ap+16)));
    try std.testing.expectEqual(0x6.3d0297b6484b39b4baebf4b5fc98p+0, acosh(@as(f128, 0x1p+8)));
    try std.testing.expectEqual(0x6.ee74dfae309ac0808d86ee993ef4p+0, acosh(@as(f128, 0x2p+8)));
    try std.testing.expectEqual(0x7.9fe703a603d23a60d7735972d5b8p+0, acosh(@as(f128, 0x4p+8)));
    try std.testing.expectEqual(0x8.51591e9dd5b8340d735715c21a1p+0, acosh(@as(f128, 0x8p+8)));
    try std.testing.expectEqual(0x9.02cb3755a78915b94082c962d54p+0, acosh(@as(f128, 0x1p+12)));
    try std.testing.expectEqual(0x9.b43d4f7d7958a5e50a739cfae1e8p+0, acosh(@as(f128, 0x2p+12)));
    try std.testing.expectEqual(0x1.1542457337d42a1c6b73c89d84acp+4, acosh(@as(f128, 0x1p+24)));
    try std.testing.expectEqual(0x1.205966f2b4f124b7281203d70653p+4, acosh(@as(f128, 0x2p+24)));
    try std.testing.expectEqual(0x1.2b708872320e1d11e4b03f1086a9p+4, acosh(@as(f128, 0x4p+24)));
    try std.testing.expectEqual(0x1.3687a9f1af2b14dca14e7a4a06e9p+4, acosh(@as(f128, 0x8p+24)));
    try std.testing.expectEqual(0x1.419ecb712c480c835decb5838728p+4, acosh(@as(f128, 0x1p+28)));
    try std.testing.expectEqual(0x1.4cb5ecf0a96504211a8af0bd0768p+4, acosh(@as(f128, 0x2p+28)));
    try std.testing.expectEqual(0x1.57cd0e702681fbbc97292bf687a7p+4, acosh(@as(f128, 0x4p+28)));
    try std.testing.expectEqual(0x1.62e42fefa39ef35783c7673007e6p+4, acosh(@as(f128, 0x8p+28)));
    try std.testing.expectEqual(0x1.6dfb516f20bbeaf24c65a2698825p+4, acosh(@as(f128, 0x1p+32)));
    try std.testing.expectEqual(0x1.791272ee9dd8e28d0c03dda30864p+4, acosh(@as(f128, 0x2p+32)));
    try std.testing.expectEqual(0x2.1f6d6966f28b649e1a4956018818p+4, acosh(@as(f128, 0x1p+48)));
    try std.testing.expectEqual(0x2.2a848ae66fa85c38d6e7913b0b58p+4, acosh(@as(f128, 0x2p+48)));
    try std.testing.expectEqual(0x2.359bac65ecc553d39385cc748c56p+4, acosh(@as(f128, 0x4p+48)));
    try std.testing.expectEqual(0x2.40b2cde569e24b6e502407ae0cc6p+4, acosh(@as(f128, 0x8p+48)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff43090cc242e78d1p+4, acosh(@as(f128, 0x1p+52)));
    try std.testing.expectEqual(0x2.56e110e4641c3aa3c9607e210d54p+4, acosh(@as(f128, 0x2p+52)));
    // try std.testing.expectEqual(0x2.61f83263e139323e85feb95a8d92p+4, acosh(@as(f128, 0x4p+52)));
    // try std.testing.expectEqual(0x2.6d0f53e35e5629d9429cf4940dd2p+4, acosh(@as(f128, 0x8p+52)));
    try std.testing.expectEqual(0x2.78267562db732173ff3b2fcd8e12p+4, acosh(@as(f128, 0x1p+56)));
    // try std.testing.expectEqual(0x2.833d96e25890190ebbd96b070e5p+4, acosh(@as(f128, 0x2p+56)));
    try std.testing.expectEqual(0x2.8e54b861d5ad10a97877a6408e9p+4, acosh(@as(f128, 0x4p+56)));
    try std.testing.expectEqual(0x2.996bd9e152ca08443515e17a0edp+4, acosh(@as(f128, 0x8p+56)));
    // try std.testing.expectEqual(0x4.6020374c5c6db00c6a6d5daf98ecp+4, acosh(@as(f128, 0x1p+100)));
    // try std.testing.expectEqual(0x5.96a7e12e0b98bcf90bb682a4468p+4, acosh(@as(f128, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0x1.5b4448e7fd9b091d321a9e787fbap+8, acosh(@as(f128, 0x1p+500)));
    // try std.testing.expectEqual(0x5.96a7e12e0b98bcf90bb682a4468p+4, acosh(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb620d358b213a7dp+8, acosh(@as(f128, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(0xd.8a6dd63831ae0fdceaf12f64a528p+8, acosh(@as(f128, 0x1p+5000)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb624d358b213a7dp+8, acosh(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    // try std.testing.expectEqual(0xf.6f3554d097450e37ef5beb7cf3d8p-4, acosh(@as(f128, 0x1.80a368p+0)));
    try std.testing.expectEqual(0x4.cd32fbd26f204f39efd8914f2758p-4, acosh(@as(f128, 0x1.0b9d4p+0)));
    try std.testing.expectEqual(0x4.cd3292c615ccfa71e8b890b497fcp-4, acosh(@as(f128, 0x1.0b9d3ep+0)));
    try std.testing.expectEqual(0x4.cd32b38d267382a8f4725060b34cp-4, acosh(@as(f128, 0x1.0b9d3e9fc19fbp+0)));
    // try std.testing.expectEqual(0x5.f3a5dc5f386bf397c2b2174968dcp-4, acosh(@as(f128, 0x1.11eab6p+0)));
    // try std.testing.expectEqual(0x5.a0a86a1d26c7fdfbf175ccccc738p-4, acosh(@as(f128, 0x1.0fffaap+0)));
    // try std.testing.expectEqual(0x3.9cf2578fe127544eaf4c16cf6496p-4, acosh(@as(f128, 0x1.068e1p+0)));
    try std.testing.expectEqual(0x3.9cf1cb0d2c9b6d25c39f58b1e97ep-4, acosh(@as(f128, 0x1.068e0ep+0)));
    try std.testing.expectEqual(0x3.9cf202813d515f3ef716fa7a96cp-4, acosh(@as(f128, 0x1.068e0eca105a6p+0)));
    try std.testing.expectEqual(0x1.96d14c020ca3c00072b5e01714f8p+0, acosh(@as(f128, 0x2.8d4288p+0)));
    try std.testing.expectEqual(0x1.9f5260500df78800bca5163ec377p+4, acosh(@as(f128, 0x1.5d71a6p+36)));
    // try std.testing.expectEqual(0x5.96a7e12e0b98bcf90bb682a4468p+4, acosh(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.890ca4134d0a7ffdd6350ab8e4ap+4, acosh(@as(f128, 0x1.476a3c43d7edep+288)));
    try std.testing.expectEqual(0x3.38018c2243818416c0f47ad794cep+0, acosh(@as(f128, 0xc.84cb2p+0)));
    try std.testing.expectEqual(0x3.38018ada065fb6ce976468972b42p+0, acosh(@as(f128, 0xc.84cb1p+0)));
    // try std.testing.expectEqual(0x3.38018bf3c5239000f070b3c63376p+0, acosh(@as(f128, 0xc.84cb1dbbd1b1p+0)));
    // try std.testing.expectEqual(0x3.c4bc22bfc0f0f17481e3028b44f6p-4, acosh(@as(f128, 0x1.072238p+0)));
    // try std.testing.expectEqual(0x3.c4bb9c223d182319078861df38bap-4, acosh(@as(f128, 0x1.072236p+0)));
    try std.testing.expectEqual(0x3.c4bba7bdf5dac0d2d260828b17e6p-4, acosh(@as(f128, 0x1.0722362c26ba5p+0)));
    try std.testing.expectEqual(0xe.aaeb3064737cace4c93072164fdp-4, acosh(@as(f128, 0x1.73515p+0)));
    try std.testing.expectEqual(0x7.5a648ae540c2605c952f79c9f67cp-4, acosh(@as(f128, 0x1.1b836p+0)));
    // try std.testing.expectEqual(0x3.efb2c425f3ffc3080db4776262p-4, acosh(@as(f128, 0x1.07c956p+0)));
    // try std.testing.expectEqual(0x6.e14ea6570758830c709fd6ddc388p-4, acosh(@as(f128, 0x1.1808eep+0)));
    try std.testing.expectEqual(0x5.af0cbc66dcd4992c060eb8ff4ap-4, acosh(@as(f128, 0x1.1052c4p+0)));
    try std.testing.expectEqual(0x7.c580228bed2a828f8042d267c3f4p-4, acosh(@as(f128, 0x1.1ecbep+0)));
    // try std.testing.expectEqual(0x7.c57fe32f85bf5c4bb5c8e26d43dcp-4, acosh(@as(f128, 0x1.1ecbdep+0)));
    // try std.testing.expectEqual(0x7.c58009b5885c945dbf38fa630eccp-4, acosh(@as(f128, 0x1.1ecbdf374bce1p+0)));
    try std.testing.expectEqual(0x7.c58009b5885a997a85c0ab84dd98p-4, acosh(@as(f128, 0x1.1ecbdf374bcep+0)));
    // try std.testing.expectEqual(0x7.c58009b5885ad2282e0dbf567ccp-4, acosh(@as(f128, 0x1.1ecbdf374bce01cap+0)));
    try std.testing.expectEqual(0x5.8e60c2aef6b859fea83073e03dd8p-4, acosh(@as(f128, 0x1.0f975ap+0)));
    try std.testing.expectEqual(0x5.8e60685c1d693545e67181d1a6ep-4, acosh(@as(f128, 0x1.0f9758p+0)));
    // try std.testing.expectEqual(0x5.8e607be4bc4cd59aa3e6945d8678p-4, acosh(@as(f128, 0x1.0f97586eba091p+0)));
    // try std.testing.expectEqual(0x5.8e607be4bc4a0303cbe388d3bc78p-4, acosh(@as(f128, 0x1.0f97586eba09p+0)));
    // try std.testing.expectEqual(0x5.8e607be4bc4a5db0f9beeaa67be4p-4, acosh(@as(f128, 0x1.0f97586eba090202p+0)));
    try std.testing.expectEqual(0x5.8e607be4bc4a5d56a6e3ea450a98p-4, acosh(@as(f128, 0x1.0f97586eba0902p+0)));
    // try std.testing.expectEqual(0x5.8e607be4bc4a5d59bfae715e3d1p-4, acosh(@as(f128, 0x1.0f97586eba090200118df0902f99p+0)));
    // try std.testing.expectEqual(0x5.8e607be4bc4a5d59bfae715e4f3cp-4, acosh(@as(f128, 0x1.0f97586eba090200118df0903p+0)));
    try std.testing.expectEqual(0x5.8e607be4bc4a5d59bfae715e38a8p-4, acosh(@as(f128, 0x1.0f97586eba090200118df0902f8p+0)));
    // try std.testing.expectEqual(0x5.96a7e12e0b98bcf90bb682a4468p+4, acosh(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb620d358b213a7dp+8, acosh(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.c5d37700c6bb03a6c23b6c9b494cp+12, acosh(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2.c5d37700c6bb03a6c24b6c9b494cp+12, acosh(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.c679d1f73f0fb624d358b213a7dp+8, acosh(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
}
