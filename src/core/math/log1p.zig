const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn log1p(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return log1p(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, log1p32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_log1pf.c
                    return log1p32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_log1p.c
                    return log1p64(x);
                },
                f80 => return cast(f80, log1p128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_log1pl.c
                    return log1p128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

noinline fn as_special(x: f32) f32 {
    const t: u32 = @bitCast(x);

    if (t == 0xbf800000)
        return -1 / @as(f32, 0);

    if (t == 0x7f800000)
        return x; // +inf

    const ax: u32 = t << 1;
    if (ax > 0xff000000)
        return x + x; // nan

    return (x - x) / (x - x);
}

fn log1p32(x: f32) f32 {
    const x0: [32]f64 =
        .{
            0x1.f81f82p-1,  0x1.e9131acp-1, 0x1.dae6077p-1, 0x1.cd85689p-1,
            0x1.c0e0704p-1, 0x1.b4e81b5p-1, 0x1.a98ef6p-1,  0x1.9ec8e95p-1,
            0x1.948b0fdp-1, 0x1.8acb90fp-1, 0x1.8181818p-1, 0x1.78a4c81p-1,
            0x1.702e05cp-1, 0x1.6816817p-1, 0x1.605816p-1,  0x1.58ed231p-1,
            0x1.51d07ebp-1, 0x1.4afd6ap-1,  0x1.446f865p-1, 0x1.3e22cbdp-1,
            0x1.3813814p-1, 0x1.323e34ap-1, 0x1.2c9fb4ep-1, 0x1.27350b9p-1,
            0x1.21fb781p-1, 0x1.1cf06aep-1, 0x1.1811812p-1, 0x1.135c811p-1,
            0x1.0ecf56cp-1, 0x1.0a6810ap-1, 0x1.0624dd3p-1, 0x1.0204081p-1,
        };
    const lixb: [32]f64 =
        .{
            0x1.fc0a8909b4218p-7, 0x1.77458f51aac89p-5, 0x1.341d793afb997p-4,
            0x1.a926d3a5ebd2ap-4, 0x1.0d77e7a8a823dp-3, 0x1.44d2b6c557102p-3,
            0x1.7ab89040accecp-3, 0x1.af3c94ecab3d6p-3, 0x1.e27076d54e6c9p-3,
            0x1.0a324e3888ad5p-2, 0x1.22941fc0c7357p-2, 0x1.3a64c56ae3fdbp-2,
            0x1.51aad874af21fp-2, 0x1.686c81d300eap-2,  0x1.7eaf83c7fa9b5p-2,
            0x1.947941aa610ecp-2, 0x1.a9cec9a3f023bp-2, 0x1.beb4d9ea4156ep-2,
            0x1.d32fe7f35e5c7p-2, 0x1.e7442617b817ap-2, 0x1.faf588dd5ed1p-2,
            0x1.0723e5c635c39p-1, 0x1.109f39d53c99p-1,  0x1.19ee6b38a4668p-1,
            0x1.23130d7f93c3bp-1, 0x1.2c0e9ec9b0b85p-1, 0x1.34e289cb35eccp-1,
            0x1.3d9026ad3d3f3p-1, 0x1.4618bc1eadbbbp-1, 0x1.4e7d8127dd8a9p-1,
            0x1.56bf9d5967092p-1, 0x1.5ee02a926936ep-1,
        };
    const lix: [32]f64 =
        .{
            0x1.fc0a890fc03e4p-7, 0x1.77458f532dcfcp-5, 0x1.341d793bbd1d1p-4,
            0x1.a926d3a6ad563p-4, 0x1.0d77e7a908e59p-3, 0x1.44d2b6c5b7d1ep-3,
            0x1.7ab890410d909p-3, 0x1.af3c94ed0bff3p-3, 0x1.e27076d5af2e6p-3,
            0x1.0a324e38b90e3p-2, 0x1.22941fc0f7966p-2, 0x1.3a64c56b145eap-2,
            0x1.51aad874df82dp-2, 0x1.686c81d3314afp-2, 0x1.7eaf83c82afc3p-2,
            0x1.947941aa916fbp-2, 0x1.a9cec9a42084ap-2, 0x1.beb4d9ea71b7cp-2,
            0x1.d32fe7f38ebd5p-2, 0x1.e7442617e8788p-2, 0x1.faf588dd8f31fp-2,
            0x1.0723e5c64df4p-1,  0x1.109f39d554c97p-1, 0x1.19ee6b38bc96fp-1,
            0x1.23130d7fabf43p-1, 0x1.2c0e9ec9c8e8cp-1, 0x1.34e289cb4e1d3p-1,
            0x1.3d9026ad556fbp-1, 0x1.4618bc1ec5ec2p-1, 0x1.4e7d8127f5bb1p-1,
            0x1.56bf9d597f399p-1, 0x1.5ee02a9281675p-1,
        };
    const b: [8]f64 =
        .{
            0x1p+0,
            -0x1p-1,
            0x1.5555555556f6bp-2,
            -0x1.00000000029b9p-2,
            0x1.9999988d176e4p-3,
            -0x1.55555418889a7p-3,
            0x1.24adeca50e2bcp-3,
            -0x1.001ba33bf57cfp-3,
        };

    var z: f64 = cast(f64, x, .{});
    const ux: u32 = @bitCast(x);
    const ax: u32 = ux & (~@as(u32, 0) >> 1);

    if (ax < 0x3c880000) {
        @branchHint(.likely);
        if (ax < 0x33000000) {
            @branchHint(.unlikely);

            if (ax == 0)
                return x;

            return @mulAdd(f32, x, -x, x);
        }
        const z2: f64 = z * z;
        const z4: f64 = z2 * z2;
        const f: f64 = z2 * ((b[1] + z * b[2]) + z2 * (b[3] + z * b[4]) + z4 * ((b[5] + z * b[6]) + z2 * b[7]));
        var r: f64 = z + f;
        if ((@as(u64, @bitCast(r)) & 0xfffffff) == 0) {
            @branchHint(.unlikely);
            r += 0x1p14 * (f + (z - r));
        }

        return cast(f32, r, .{});
    } else {
        if (ux >= 0xbf800000 or ax >= 0x7f800000) {
            @branchHint(.unlikely);
            return as_special(x);
        }

        const tp: u64 = @bitCast(z + 1);
        var e: i32 = @intCast(tp >> 52);
        const m52: u64 = tp & (~@as(u64, 0) >> 12);
        const j: u32 = @intCast((tp >> 47) & 31);
        e -= 0x3ff;
        const xd: f64 = @bitCast(m52 | (0x3ff << 52));
        z = xd * x0[j] - 1;
        const c: [5]f64 =
            .{ -0x1.3902c33434e7fp-43, 0x1.ffffffe1cbed5p-1, -0x1.ffffff7d1b014p-2, 0x1.5564e0ed3613ap-2, -0x1.0012232a00d4ap-2 };
        const ln2: f64 = 0x1.62e42fefa39efp-1;
        const z2: f64 = z * z;
        const r: f64 = (ln2 * cast(f64, e, .{}) + lixb[j]) + z * ((c[1] + z * c[2]) + z2 * (c[3] + z * c[4]));
        var ub: f32 = cast(f32, r, .{});
        const lb: f32 = cast(f32, r + 2.2e-11, .{});
        if (ub != lb) {
            @branchHint(.unlikely);
            const z4: f64 = z2 * z2;
            const f: f64 = z * ((b[0] + z * b[1]) + z2 * (b[2] + z * b[3]) + z4 * ((b[4] + z * b[5]) + z2 * (b[6] + z * b[7])));
            const ln2l: f64 = 0x1.7f7d1cf79abcap-20;
            const ln2h: f64 = 0x1.62e4p-1;
            const Lh: f64 = ln2h * cast(f64, e, .{});
            const Ll: f64 = ln2l * cast(f64, e, .{});
            const rl: f64 = f + Ll + lix[j];
            var tr: f64 = rl + Lh;
            if ((@as(u64, @bitCast(tr)) & 0xfffffff) == 0) {
                @branchHint(.unlikely);

                if (x == -0x1.247ab0p-6)
                    return -0x1.271f0ep-6 - 0x1p-31;

                if (x == -0x1.3a415ep-5)
                    return -0x1.407112p-5 + 0x1p-30;

                if (x == 0x1.fb035ap-2)
                    return 0x1.9bddc2p-2 + 0x1p-27;

                tr += 64 * (rl + (Lh - tr));
            } else if (rl + (Lh - tr) == 0.0) {
                if (x == 0x1.b7fd86p-4)
                    return 0x1.a1ece2p-4 + 0x1p-29;

                if (x == -0x1.3a415ep-5)
                    return -0x1.407112p-5 + 0x1p-30;

                if (x == 0x1.43c7e2p-6)
                    return 0x1.409f80p-6 + 0x1p-31;
            }

            ub = cast(f32, tr, .{});
        }
        return ub;
    }
}

fn log1p64(x: f64) f64 {
    const ln2_hi: f64 = 6.93147180369123816490e-01; // 3fe62e42 fee00000
    const ln2_lo: f64 = 1.90821492927058770002e-10; // 3dea39ef 35793c76
    const two54: f64 = 1.80143985094819840000e+16; // 43500000 00000000
    const Lp: [8]f64 = .{
        0.0,
        6.666666666666735130e-01, // 3FE55555 55555593
        3.999999999940941908e-01, // 3FD99999 9997FA04
        2.857142874366239149e-01, // 3FD24924 94229359
        2.222219843214978396e-01, // 3FCC71C5 1D8E78AF
        1.818357216161805012e-01, // 3FC74664 96CB03DE
        1.531383769920937332e-01, // 3FC39A09 D078C69F
        1.479819860511658591e-01, // 3FC2F112 DF3E5244
    };

    var hx: i32 = undefined;
    dbl64.getHighWord(&hx, x);
    const ax: i32 = hx & 0x7fffffff;

    var k: i32 = 1;
    var f: f64 = 0;
    var hu: i32 = 0;
    if (hx < 0x3fda827a) { // x < 0.41422
        if (ax >= 0x3ff00000) { // x <= -1.0
            @branchHint(.unlikely);

            if (x == -1.0) {
                return -two54 / 0; // log1p(-1)=-inf
            } else {
                return (x - x) / (x - x); // log1p(x<-1)=NaN
            }
        }

        if (ax < 0x3e200000) { // |x| < 2**-29
            @branchHint(.unlikely);
            std.mem.doNotOptimizeAway(two54 + x); // raise inexact
            if (ax < 0x3c900000) // |x| < 2**-54
            {
                if (@abs(x) < std.math.floatMin(f64)) {
                    const vx: f64 = x * x;
                    std.mem.doNotOptimizeAway(vx);
                }
                return x;
            } else {
                return x - x * x * 0.5;
            }
        }

        if (hx > 0 or hx <= @as(i32, @bitCast(@as(u32, 0xbfd2bec3)))) {
            k = 0;
            f = x;
            hu = 1;
        } // -0.2929<x<0.41422
    } else if (hx >= 0x7ff00000) {
        @branchHint(.unlikely);
        return x + x;
    }

    var c: f64 = undefined;
    if (k != 0) {
        var u: f64 = undefined;
        if (hx < 0x43400000) {
            u = 1 + x;
            dbl64.getHighWord(&hu, u);
            k = (hu >> 20) - 1023;
            c = if (k > 0) 1 - (u - x) else x - (u - 1); // correction term
            c /= u;
        } else {
            u = x;
            dbl64.getHighWord(&hu, u);
            k = (hu >> 20) - 1023;
            c = 0;
        }

        hu &= 0x000fffff;
        if (hu < 0x6a09e) {
            dbl64.setHighWord(&u, hu | 0x3ff00000); // normalize u
        } else {
            k += 1;
            dbl64.setHighWord(&u, hu | 0x3fe00000); // normalize u/2
            hu = (0x00100000 - hu) >> 2;
        }
        f = u - 1;
    }

    const hfsq: f64 = 0.5 * f * f;
    if (hu == 0) { // |f| < 2**-20
        if (f == 0) {
            if (k == 0) {
                return 0;
            } else {
                c += cast(f64, k, .{}) * ln2_lo;
                return cast(f64, k, .{}) * ln2_hi + c;
            }
        }
        const R: f64 = hfsq * (1 - 0.66666666666666666 * f);
        if (k == 0) {
            return f - R;
        } else {
            return cast(f64, k, .{}) * ln2_hi - ((R - (cast(f64, k, .{}) * ln2_lo + c)) - f);
        }
    }

    const s: f64 = f / (2 + f);
    const z: f64 = s * s;
    const R1: f64 = z * Lp[1];
    const z2: f64 = z * z;
    const R2: f64 = Lp[2] + z * Lp[3];
    const z4: f64 = z2 * z2;
    const R3: f64 = Lp[4] + z * Lp[5];
    const z6: f64 = z4 * z2;
    const R4: f64 = Lp[6] + z * Lp[7];
    const R: f64 = R1 + z2 * R2 + z4 * R3 + z6 * R4;

    if (k == 0) {
        return f - (hfsq - s * (hfsq + R));
    } else {
        return cast(f64, k, .{}) * ln2_hi - ((hfsq - (s * (hfsq + R) + (cast(f64, k, .{}) * ln2_lo + c))) - f);
    }
}

fn log1p128(xm1: f128) f128 {
    const P12: f128 = 1.538612243596254322971797716843006400388e-6;
    const P11: f128 = 4.998469661968096229986658302195402690910e-1;
    const P10: f128 = 2.321125933898420063925789532045674660756e1;
    const P9: f128 = 4.114517881637811823002128927449878962058e2;
    const P8: f128 = 3.824952356185897735160588078446136783779e3;
    const P7: f128 = 2.128857716871515081352991964243375186031e4;
    const P6: f128 = 7.594356839258970405033155585486712125861e4;
    const P5: f128 = 1.797628303815655343403735250238293741397e5;
    const P4: f128 = 2.854829159639697837788887080758954924001e5;
    const P3: f128 = 3.007007295140399532324943111654767187848e5;
    const P2: f128 = 2.014652742082537582487669938141683759923e5;
    const P1: f128 = 7.771154681358524243729929227226708890930e4;
    const P0: f128 = 1.313572404063446165910279910527789794488e4;
    // const Q12: f128 = 1.000000000000000000000000000000000000000e0;
    const Q11: f128 = 4.839208193348159620282142911143429644326e1;
    const Q10: f128 = 9.104928120962988414618126155557301584078e2;
    const Q9: f128 = 9.147150349299596453976674231612674085381e3;
    const Q8: f128 = 5.605842085972455027590989944010492125825e4;
    const Q7: f128 = 2.248234257620569139969141618556349415120e5;
    const Q6: f128 = 6.132189329546557743179177159925690841200e5;
    const Q5: f128 = 1.158019977462989115839826904108208787040e6;
    const Q4: f128 = 1.514882452993549494932585972882995548426e6;
    const Q3: f128 = 1.347518538384329112529391120390701166528e6;
    const Q2: f128 = 7.777690340007566932935753241556479363645e5;
    const Q1: f128 = 2.626900195321832660448791748036714883242e5;
    const Q0: f128 = 3.940717212190338497730839731583397586124e4;

    // Coefficients for log(x) = z + z^3 P(z^2)/Q(z^2),
    // where z = 2(x-1)/(x+1)
    // 1/sqrt(2) <= x < sqrt(2)
    // Theoretical peak relative error = 1.1e-35,
    // relative peak error spread 1.1e-9
    const R5: f128 = -8.828896441624934385266096344596648080902e-1;
    const R4: f128 = 8.057002716646055371965756206836056074715e1;
    const R3: f128 = -2.024301798136027039250415126250455056397e3;
    const R2: f128 = 2.048819892795278657810231591630928516206e4;
    const R1: f128 = -8.977257995689735303686582344659576526998e4;
    const R0: f128 = 1.418134209872192732479751274970992665513e5;
    // const S6: f128 = 1.000000000000000000000000000000000000000e0;
    const S5: f128 = -1.186359407982897997337150403816839480438e2;
    const S4: f128 = 3.998526750980007367835804959888064681098e3;
    const S3: f128 = -5.748542087379434595104154610899551484314e4;
    const S2: f128 = 4.001557694070773974936904547424676279307e5;
    const S1: f128 = -1.332535117259762928288745111081235577029e6;
    const S0: f128 = 1.701761051846631278975701529965589676574e6;

    // C1 + C2 = ln 2
    const C1: f128 = 6.93145751953125e-1;
    const C2: f128 = 1.428606820309417232121458176568075500134e-6;

    const sqrth: f128 = 0.7071067811865475244008443621048490392848;
    // ln (2^16384 * (1 - 2^-113))

    // Test for NaN or infinity input.
    const u: ldbl128.ieee_f128_shape32 = @bitCast(xm1);
    const hx: i32 = @bitCast(u.w0);
    if ((hx & 0x7fffffff) >= 0x7fff0000)
        return xm1 + @abs(xm1);

    // log1p(+- 0) = +- 0.
    if ((hx & 0x7fffffff) == 0 and (u.w1 | u.w2 | u.w3) == 0)
        return xm1;

    if ((hx & 0x7fffffff) < 0x3f8e0000) {
        if (@abs(xm1) < std.math.floatMin(f128)) {
            const vxm1: f128 = xm1 * xm1;
            std.mem.doNotOptimizeAway(vxm1);
        }
        if (xm1 == 0)
            return xm1;
    }

    var x: f128 = undefined;
    if (xm1 >= 0x1p113) {
        x = xm1;
    } else {
        x = xm1 + 1;
    }

    // log1p(-1) = -inf
    if (x <= 0) {
        if (x == 0) {
            return -1 / @as(f128, 0); // log1p(-1) = -inf
        } else {
            return 0 / (x - x);
        }
    }

    // Separate mantissa from exponent.

    // Use frexp used so that denormal numbers will be handled properly.
    var e: i32 = undefined;
    x = math.frexp(x, &e);

    // Logarithm using log(x) = z + z^3 P(z^2)/Q(z^2), where z = 2(x-1)/(x+1).
    if ((e > 2) or (e < -2)) {
        var z: f128 = undefined;
        var y: f128 = undefined;
        if (x < sqrth) { // 2( 2x-1 )/( 2x+1 )
            e -= 1;
            z = x - 0.5;
            y = 0.5 * z + 0.5;
        } else { // 2 (x-1)/(x+1)
            z = x - 0.5;
            z -= 0.5;
            y = 0.5 * x + 0.5;
        }
        x = z / y;
        z = x * x;
        const r: f128 = ((((R5 * z + R4) * z + R3) * z + R2) * z + R1) * z + R0;
        const s: f128 = (((((z + S5) * z + S4) * z + S3) * z + S2) * z + S1) * z + S0;
        z = x * (z * r / s);
        z = z + cast(f128, e, .{}) * C2;
        z = z + x;
        z = z + cast(f128, e, .{}) * C1;
        return z;
    }

    // Logarithm using log(1+x) = x - .5x^2 + x^3 P(x)/Q(x).
    if (x < sqrth) {
        e -= 1;
        if (e != 0) {
            x = 2 * x - 1; // 2x - 1
        } else {
            x = xm1;
        }
    } else {
        if (e != 0) {
            x = x - 1;
        } else {
            x = xm1;
        }
    }
    var z: f128 = x * x;
    const r: f128 = (((((((((((P12 * x + P11) * x + P10) * x + P9) * x + P8) * x + P7) * x + P6) * x + P5) * x + P4) * x + P3) * x + P2) * x + P1) * x + P0;
    const s: f128 = (((((((((((x + Q11) * x + Q10) * x + Q9) * x + Q8) * x + Q7) * x + Q6) * x + Q5) * x + Q4) * x + Q3) * x + Q2) * x + Q1) * x + Q0;
    var y: f128 = x * (z * r / s);
    y = y + cast(f128, e, .{}) * C2;
    z = y - 0.5 * z;
    z = z + x;
    z = z + cast(f128, e, .{}) * C1;
    return z;
}

test log1p {
    try std.testing.expectEqual(0x0p+0, log1p(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, log1p(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, log1p(@as(f32, 0x1.b7e152p+0)));
    try std.testing.expectEqual(0xf.fffffp-4, log1p(@as(f32, 0x1.b7e15p+0)));
    try std.testing.expectEqual(-0x4.9a5888p-4, log1p(@as(f32, -0x4p-4)));
    try std.testing.expectEqual(-0x2.145648p+0, log1p(@as(f32, -0xep-4)));
    try std.testing.expectEqual(0x7.e0a6cp-8, log1p(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x3.ff8014p-12, log1p(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffep-16, log1p(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0xf.ffff8p-24, log1p(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0x8p-28, log1p(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x4p-32, log1p(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, log1p(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, log1p(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, log1p(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, log1p(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, log1p(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, log1p(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, log1p(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, log1p(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, log1p(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x4p-128, log1p(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x4p-128, log1p(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(-0x8p-152, log1p(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.eeb4e8p+0, log1p(@as(f32, 0x4p+8)));
    try std.testing.expectEqual(0xd.dce9fp+0, log1p(@as(f32, 0x1p+20)));
    try std.testing.expectEqual(0x1.4cb5ecp+4, log1p(@as(f32, 0x4p+28)));
    try std.testing.expectEqual(0x2.2a848cp+4, log1p(@as(f32, 0x4p+48)));
    try std.testing.expectEqual(0x2.996bd8p+4, log1p(@as(f32, 0x1p+60)));
    try std.testing.expectEqual(0x4.550918p+4, log1p(@as(f32, 0x1p+100)));
    try std.testing.expectEqual(0x5.8b90cp+4, log1p(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x5.8bb35p-4, log1p(@as(f32, 0x6.a0cf48p-4)));
    try std.testing.expectEqual(0x5.8bb348p-4, log1p(@as(f32, 0x6.a0cf4p-4)));
    try std.testing.expectEqual(0x5.8b90cp+4, log1p(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x5.ebc1a8p-4, log1p(@as(f32, 0x7.2a4368p-4)));
    try std.testing.expectEqual(0x5.af7a38p-4, log1p(@as(f32, 0x6.d3a118p-4)));
    try std.testing.expectEqual(0x1.cb58e4p+0, log1p(@as(f32, 0x5.03f228p+0)));
    try std.testing.expectEqual(0x5.e90248p-4, log1p(@as(f32, 0x7.264968p-4)));
    try std.testing.expectEqual(0x5.e90248p-4, log1p(@as(f32, 0x7.26496p-4)));
    try std.testing.expectEqual(0x6.cc4a08p-4, log1p(@as(f32, 0x8.786bdp-4)));
    try std.testing.expectEqual(0x6.2d42p-4, log1p(@as(f32, 0x7.89dc18p-4)));
    try std.testing.expectEqual(0x6.2d42p-4, log1p(@as(f32, 0x7.89dc1p-4)));
    try std.testing.expectEqual(0x7.763bb8p-4, log1p(@as(f32, 0x9.81cdp-4)));
    try std.testing.expectEqual(0x7.763bbp-4, log1p(@as(f32, 0x9.81ccfp-4)));
    try std.testing.expectEqual(0x7.f5ac6p-4, log1p(@as(f32, 0xa.50287p-4)));
    try std.testing.expectEqual(0x7.f5ac58p-4, log1p(@as(f32, 0xa.50286p-4)));
    try std.testing.expectEqual(0x4.e92f7p-4, log1p(@as(f32, 0x5.bf7888p-4)));
    try std.testing.expectEqual(0x4.e92f68p-4, log1p(@as(f32, 0x5.bf788p-4)));
    try std.testing.expectEqual(0x6.43432p-4, log1p(@as(f32, 0x7.aa5198p-4)));
    try std.testing.expectEqual(0x1.34829cp+0, log1p(@as(f32, 0x2.564fap+0)));
    try std.testing.expectEqual(0x6.7a3638p-4, log1p(@as(f32, 0x7.fc243p-4)));
    try std.testing.expectEqual(0x6.7a3638p-4, log1p(@as(f32, 0x7.fc2428p-4)));
    try std.testing.expectEqual(-0x5.ec9648p-4, log1p(@as(f32, -0x4.f37d38p-4)));
    try std.testing.expectEqual(-0x5.ec9658p-4, log1p(@as(f32, -0x4.f37d4p-4)));
    try std.testing.expectEqual(0x5.eee1dp-4, log1p(@as(f32, 0x7.2eca58p-4)));
    try std.testing.expectEqual(0x5.eee1c8p-4, log1p(@as(f32, 0x7.2eca5p-4)));
    try std.testing.expectEqual(-0x7.ecba98p-4, log1p(@as(f32, -0x6.3fef3p-4)));
    try std.testing.expectEqual(-0x7.ecbaap-4, log1p(@as(f32, -0x6.3fef38p-4)));
    try std.testing.expectEqual(0x5.95f3fp-4, log1p(@as(f32, 0x6.af53d8p-4)));
    try std.testing.expectEqual(0x5.95f3fp-4, log1p(@as(f32, 0x6.af53dp-4)));

    try std.testing.expectEqual(0x0p+0, log1p(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, log1p(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1.00000039ece12p+0, log1p(@as(f64, 0x1.b7e152p+0)));
    try std.testing.expectEqual(0xf.fffff7d922f5p-4, log1p(@as(f64, 0x1.b7e15p+0)));
    try std.testing.expectEqual(0x1p+0, log1p(@as(f64, 0x1.b7e151628aed3p+0)));
    try std.testing.expectEqual(0x1p+0, log1p(@as(f64, 0x1.b7e151628aed2p+0)));
    try std.testing.expectEqual(-0x4.9a58844d36e48p-4, log1p(@as(f64, -0x4p-4)));
    try std.testing.expectEqual(-0x2.145647e7756e6p+0, log1p(@as(f64, -0xep-4)));
    try std.testing.expectEqual(0x7.e0a6c39e0ccp-8, log1p(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0x3.ff8015515622p-12, log1p(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffe0002aaa6bp-16, log1p(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0xf.ffff800005558p-24, log1p(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0x7.fffffe000000cp-28, log1p(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x3.fffffff8p-32, log1p(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x1.ffffffffep-36, log1p(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0xf.fffffffff8p-44, log1p(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x7.ffffffffffep-48, log1p(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0x3.ffffffffffff8p-52, log1p(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, log1p(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, log1p(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, log1p(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, log1p(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, log1p(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, log1p(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, log1p(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x4p-128, log1p(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, log1p(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, log1p(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, log1p(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x4p-128, log1p(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, log1p(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, log1p(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0x8p-152, log1p(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, log1p(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x6.eeb4e7af87304p+0, log1p(@as(f64, 0x4p+8)));
    try std.testing.expectEqual(0xd.dce9ef5c63b58p+0, log1p(@as(f64, 0x1p+20)));
    try std.testing.expectEqual(0x1.4cb5ecf0e965p+4, log1p(@as(f64, 0x4p+28)));
    try std.testing.expectEqual(0x2.2a848ae66fa86p+4, log1p(@as(f64, 0x4p+48)));
    try std.testing.expectEqual(0x2.996bd9e152cap+4, log1p(@as(f64, 0x1p+60)));
    try std.testing.expectEqual(0x4.550915ccdf50cp+4, log1p(@as(f64, 0x1p+100)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bcp+4, log1p(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.b525ada00b928p+8, log1p(@as(f64, 0x1p+1000)));
    try std.testing.expectEqual(0x5.8bb34ffb3a0b4p-4, log1p(@as(f64, 0x6.a0cf48p-4)));
    try std.testing.expectEqual(0x5.8bb34a531ea8cp-4, log1p(@as(f64, 0x6.a0cf4p-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0cp-4, log1p(@as(f64, 0x6.a0cf42befceap-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0cp-4, log1p(@as(f64, 0x6.a0cf42befce9cp-4)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bcp+4, log1p(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c5c85fdf473dep+8, log1p(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x5.ebc1a69570c14p-4, log1p(@as(f64, 0x7.2a4368p-4)));
    try std.testing.expectEqual(0x5.af7a38286eaf8p-4, log1p(@as(f64, 0x6.d3a118p-4)));
    // try std.testing.expectEqual(0x1.cb58e45e6b3a5p+0, log1p(@as(f64, 0x5.03f228p+0)));
    try std.testing.expectEqual(0x5.e90249b494e6p-4, log1p(@as(f64, 0x7.264968p-4)));
    try std.testing.expectEqual(0x5.e902442d1717cp-4, log1p(@as(f64, 0x7.26496p-4)));
    try std.testing.expectEqual(0x5.e902469e458c4p-4, log1p(@as(f64, 0x7.264963888ac9p-4)));
    try std.testing.expectEqual(0x6.cc4a0b2426578p-4, log1p(@as(f64, 0x8.786bdp-4)));
    try std.testing.expectEqual(0x6.2d4201bc5b648p-4, log1p(@as(f64, 0x7.89dc18p-4)));
    try std.testing.expectEqual(0x6.2d41fc4c416bcp-4, log1p(@as(f64, 0x7.89dc1p-4)));
    try std.testing.expectEqual(0x6.2d420160a1d24p-4, log1p(@as(f64, 0x7.89dc17790eeb4p-4)));
    try std.testing.expectEqual(0x7.763bb83b8a4b8p-4, log1p(@as(f64, 0x9.81cdp-4)));
    // try std.testing.expectEqual(0x7.763bae3235afcp-4, log1p(@as(f64, 0x9.81ccfp-4)));
    // try std.testing.expectEqual(0x7.763bb38c7d2d4p-4, log1p(@as(f64, 0x9.81ccf8887c25p-4)));
    // try std.testing.expectEqual(0x7.763bb38c7d2ccp-4, log1p(@as(f64, 0x9.81ccf8887c248p-4)));
    try std.testing.expectEqual(0x7.f5ac5ef8280ecp-4, log1p(@as(f64, 0xa.50287p-4)));
    try std.testing.expectEqual(0x7.f5ac553d89168p-4, log1p(@as(f64, 0xa.50286p-4)));
    try std.testing.expectEqual(0x7.f5ac559290b5cp-4, log1p(@as(f64, 0xa.5028608bd65f8p-4)));
    try std.testing.expectEqual(0x7.f5ac559290b58p-4, log1p(@as(f64, 0xa.5028608bd65fp-4)));
    try std.testing.expectEqual(0x4.e92f6c6a921dp-4, log1p(@as(f64, 0x5.bf7888p-4)));
    try std.testing.expectEqual(0x4.e92f6687da53p-4, log1p(@as(f64, 0x5.bf788p-4)));
    try std.testing.expectEqual(0x4.e92f6bdbf037p-4, log1p(@as(f64, 0x5.bf78873e20a3p-4)));
    try std.testing.expectEqual(0x4.e92f6bdbf037p-4, log1p(@as(f64, 0x5.bf78873e20a2cp-4)));
    try std.testing.expectEqual(0x6.43431f9d697f8p-4, log1p(@as(f64, 0x7.aa5198p-4)));
    // try std.testing.expectEqual(0x1.34829b3156a23p+0, log1p(@as(f64, 0x2.564fap+0)));
    try std.testing.expectEqual(0x6.7a363b8f2519cp-4, log1p(@as(f64, 0x7.fc243p-4)));
    // try std.testing.expectEqual(0x6.7a363638f41ccp-4, log1p(@as(f64, 0x7.fc2428p-4)));
    try std.testing.expectEqual(0x6.7a3637a5521e4p-4, log1p(@as(f64, 0x7.fc242a2235224p-4)));
    try std.testing.expectEqual(0x6.7a3637a5521ep-4, log1p(@as(f64, 0x7.fc242a223522p-4)));
    // try std.testing.expectEqual(-0x5.ec9649184d39cp-4, log1p(@as(f64, -0x4.f37d38p-4)));
    try std.testing.expectEqual(-0x5.ec9654ae08e38p-4, log1p(@as(f64, -0x4.f37d4p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a4p-4, log1p(@as(f64, -0x4.f37d3c9ce0b14p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a8p-4, log1p(@as(f64, -0x4.f37d3c9ce0b18p-4)));
    try std.testing.expectEqual(0x5.eee1d129eb634p-4, log1p(@as(f64, 0x7.2eca58p-4)));
    try std.testing.expectEqual(0x5.eee1cba474cc4p-4, log1p(@as(f64, 0x7.2eca5p-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c4p-4, log1p(@as(f64, 0x7.2eca50c4d931cp-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508cp-4, log1p(@as(f64, 0x7.2eca50c4d9318p-4)));
    try std.testing.expectEqual(-0x7.ecba94fcebfcp-4, log1p(@as(f64, -0x6.3fef3p-4)));
    // try std.testing.expectEqual(-0x7.ecbaa21da76b4p-4, log1p(@as(f64, -0x6.3fef38p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e86p-4, log1p(@as(f64, -0x6.3fef3067427e4p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e868p-4, log1p(@as(f64, -0x6.3fef3067427e8p-4)));
    try std.testing.expectEqual(0x5.95f3f1dfd7b6p-4, log1p(@as(f64, 0x6.af53d8p-4)));
    try std.testing.expectEqual(0x5.95f3ec3b5b154p-4, log1p(@as(f64, 0x6.af53dp-4)));
    try std.testing.expectEqual(0x5.95f3ec4683fa4p-4, log1p(@as(f64, 0x6.af53d00fd2848p-4)));
    try std.testing.expectEqual(0x5.95f3ec4683fap-4, log1p(@as(f64, 0x6.af53d00fd2844p-4)));

    try std.testing.expectEqual(0x0p+0, log1p(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, log1p(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1.00000039ece11db6p+0, log1p(@as(f80, 0x1.b7e152p+0)));
    try std.testing.expectEqual(0xf.fffff7d922f51a3p-4, log1p(@as(f80, 0x1.b7e15p+0)));
    try std.testing.expectEqual(0x1.000000000000020ep+0, log1p(@as(f80, 0x1.b7e151628aed3p+0)));
    try std.testing.expectEqual(0xf.ffffffffffffc2bp-4, log1p(@as(f80, 0x1.b7e151628aed2p+0)));
    try std.testing.expectEqual(0x1p+0, log1p(@as(f80, 0x1.b7e151628aed2a6cp+0)));
    try std.testing.expectEqual(0x1p+0, log1p(@as(f80, 0x1.b7e151628aed2a6ap+0)));
    try std.testing.expectEqual(-0x4.9a58844d36e49e1p-4, log1p(@as(f80, -0x4p-4)));
    try std.testing.expectEqual(-0x2.145647e7756e6d04p+0, log1p(@as(f80, -0xep-4)));
    try std.testing.expectEqual(0x7.e0a6c39e0cc0134p-8, log1p(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x3.ff8015515621f78p-12, log1p(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffe0002aaa6aab2p-16, log1p(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0xf.ffff80000555551p-24, log1p(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0x7.fffffe000000aaa8p-28, log1p(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x3.fffffff800000014p-32, log1p(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x1.ffffffffep-36, log1p(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0xf.fffffffff8p-44, log1p(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x7.ffffffffffep-48, log1p(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0x3.ffffffffffff8p-52, log1p(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x1.fffffffffffffep-56, log1p(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0xf.ffffffffffffff8p-64, log1p(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, log1p(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, log1p(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, log1p(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, log1p(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, log1p(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, log1p(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x4p-128, log1p(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, log1p(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, log1p(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, log1p(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, log1p(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, log1p(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, log1p(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x4p-128, log1p(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, log1p(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, log1p(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, log1p(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, log1p(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x8p-152, log1p(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, log1p(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, log1p(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x6.eeb4e7af873022d8p+0, log1p(@as(f80, 0x4p+8)));
    try std.testing.expectEqual(0xd.dce9ef5c63b5817p+0, log1p(@as(f80, 0x1p+20)));
    try std.testing.expectEqual(0x1.4cb5ecf0e9650422p+4, log1p(@as(f80, 0x4p+28)));
    try std.testing.expectEqual(0x2.2a848ae66fa86038p+4, log1p(@as(f80, 0x4p+48)));
    try std.testing.expectEqual(0x2.996bd9e152ca0844p+4, log1p(@as(f80, 0x1p+60)));
    try std.testing.expectEqual(0x4.550915ccdf50b87p+4, log1p(@as(f80, 0x1p+100)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bc56p+4, log1p(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.b525ada00b927348p+8, log1p(@as(f80, 0x1p+1000)));
    try std.testing.expectEqual(0x5.8bb34ffb3a0b27dp-4, log1p(@as(f80, 0x6.a0cf48p-4)));
    try std.testing.expectEqual(0x5.8bb34a531ea8cddp-4, log1p(@as(f80, 0x6.a0cf4p-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0d198p-4, log1p(@as(f80, 0x6.a0cf42befceap-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0a458p-4, log1p(@as(f80, 0x6.a0cf42befce9cp-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0c46p-4, log1p(@as(f80, 0x6.a0cf42befce9ed48p-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0c458p-4, log1p(@as(f80, 0x6.a0cf42befce9ed4p-4)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bc56p+4, log1p(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6a8p+8, log1p(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6bp+12, log1p(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x5.ebc1a69570c135d8p-4, log1p(@as(f80, 0x7.2a4368p-4)));
    try std.testing.expectEqual(0x5.af7a38286eaf6f08p-4, log1p(@as(f80, 0x6.d3a118p-4)));
    try std.testing.expectEqual(0x1.cb58e45e6b3a48e8p+0, log1p(@as(f80, 0x5.03f228p+0)));
    try std.testing.expectEqual(0x5.e90249b494e60a88p-4, log1p(@as(f80, 0x7.264968p-4)));
    try std.testing.expectEqual(0x5.e902442d1717cc6p-4, log1p(@as(f80, 0x7.26496p-4)));
    try std.testing.expectEqual(0x5.e902469e458c388p-4, log1p(@as(f80, 0x7.264963888ac9p-4)));
    try std.testing.expectEqual(0x6.cc4a0b2426577f78p-4, log1p(@as(f80, 0x8.786bdp-4)));
    try std.testing.expectEqual(0x6.2d4201bc5b6462cp-4, log1p(@as(f80, 0x7.89dc18p-4)));
    try std.testing.expectEqual(0x6.2d41fc4c416bbd6p-4, log1p(@as(f80, 0x7.89dc1p-4)));
    try std.testing.expectEqual(0x6.2d420160a1d24358p-4, log1p(@as(f80, 0x7.89dc17790eeb4p-4)));
    try std.testing.expectEqual(0x7.763bb83b8a4b8cdp-4, log1p(@as(f80, 0x9.81cdp-4)));
    try std.testing.expectEqual(0x7.763bae3235afd3c8p-4, log1p(@as(f80, 0x9.81ccfp-4)));
    try std.testing.expectEqual(0x7.763bb38c7d2d29dp-4, log1p(@as(f80, 0x9.81ccf8887c25p-4)));
    try std.testing.expectEqual(0x7.763bb38c7d2cd988p-4, log1p(@as(f80, 0x9.81ccf8887c248p-4)));
    try std.testing.expectEqual(0x7.763bb38c7d2cf268p-4, log1p(@as(f80, 0x9.81ccf8887c24a7bp-4)));
    try std.testing.expectEqual(0x7.f5ac5ef8280eaea8p-4, log1p(@as(f80, 0xa.50287p-4)));
    try std.testing.expectEqual(0x7.f5ac553d891675c8p-4, log1p(@as(f80, 0xa.50286p-4)));
    try std.testing.expectEqual(0x7.f5ac559290b5bb48p-4, log1p(@as(f80, 0xa.5028608bd65f8p-4)));
    try std.testing.expectEqual(0x7.f5ac559290b56d78p-4, log1p(@as(f80, 0xa.5028608bd65fp-4)));
    try std.testing.expectEqual(0x7.f5ac559290b59p-4, log1p(@as(f80, 0xa.5028608bd65f38dp-4)));
    try std.testing.expectEqual(0x4.e92f6c6a921d11d8p-4, log1p(@as(f80, 0x5.bf7888p-4)));
    try std.testing.expectEqual(0x4.e92f6687da531adp-4, log1p(@as(f80, 0x5.bf788p-4)));
    try std.testing.expectEqual(0x4.e92f6bdbf0371278p-4, log1p(@as(f80, 0x5.bf78873e20a3p-4)));
    try std.testing.expectEqual(0x4.e92f6bdbf036e36p-4, log1p(@as(f80, 0x5.bf78873e20a2cp-4)));
    try std.testing.expectEqual(0x4.e92f6bdbf036f268p-4, log1p(@as(f80, 0x5.bf78873e20a2d468p-4)));
    try std.testing.expectEqual(0x6.43431f9d697f9598p-4, log1p(@as(f80, 0x7.aa5198p-4)));
    try std.testing.expectEqual(0x1.34829b3156a228p+0, log1p(@as(f80, 0x2.564fap+0)));
    try std.testing.expectEqual(0x6.7a363b8f2519bbd8p-4, log1p(@as(f80, 0x7.fc243p-4)));
    try std.testing.expectEqual(0x6.7a363638f41ca3d8p-4, log1p(@as(f80, 0x7.fc2428p-4)));
    try std.testing.expectEqual(0x6.7a3637a5521e27ap-4, log1p(@as(f80, 0x7.fc242a2235224p-4)));
    try std.testing.expectEqual(0x6.7a3637a5521dfce8p-4, log1p(@as(f80, 0x7.fc242a223522p-4)));
    try std.testing.expectEqual(0x6.7a3637a5521e1c4p-4, log1p(@as(f80, 0x7.fc242a2235222ef8p-4)));
    try std.testing.expectEqual(-0x5.ec9649184d39a5d8p-4, log1p(@as(f80, -0x4.f37d38p-4)));
    try std.testing.expectEqual(-0x5.ec9654ae08e379bp-4, log1p(@as(f80, -0x4.f37d4p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a2d6p-4, log1p(@as(f80, -0x4.f37d3c9ce0b14p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a8a1p-4, log1p(@as(f80, -0x4.f37d3c9ce0b18p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a3e88p-4, log1p(@as(f80, -0x4.f37d3c9ce0b14bd8p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a3e9p-4, log1p(@as(f80, -0x4.f37d3c9ce0b14bep-4)));
    try std.testing.expectEqual(0x5.eee1d129eb6330c8p-4, log1p(@as(f80, 0x7.2eca58p-4)));
    try std.testing.expectEqual(0x5.eee1cba474cc2b48p-4, log1p(@as(f80, 0x7.2eca5p-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c3b38p-4, log1p(@as(f80, 0x7.2eca50c4d931cp-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c0f08p-4, log1p(@as(f80, 0x7.2eca50c4d9318p-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c1e6p-4, log1p(@as(f80, 0x7.2eca50c4d9319638p-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c1e58p-4, log1p(@as(f80, 0x7.2eca50c4d931963p-4)));
    try std.testing.expectEqual(-0x7.ecba94fcebfbf0ap-4, log1p(@as(f80, -0x6.3fef3p-4)));
    try std.testing.expectEqual(-0x7.ecbaa21da76b5868p-4, log1p(@as(f80, -0x6.3fef38p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e861fep-4, log1p(@as(f80, -0x6.3fef3067427e4p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e8688ep-4, log1p(@as(f80, -0x6.3fef3067427e8p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e86263p-4, log1p(@as(f80, -0x6.3fef3067427e43d8p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e862638p-4, log1p(@as(f80, -0x6.3fef3067427e43ep-4)));
    try std.testing.expectEqual(0x5.95f3f1dfd7b5e048p-4, log1p(@as(f80, 0x6.af53d8p-4)));
    try std.testing.expectEqual(0x5.95f3ec3b5b15424p-4, log1p(@as(f80, 0x6.af53dp-4)));
    try std.testing.expectEqual(0x5.95f3ec4683fa2d2p-4, log1p(@as(f80, 0x6.af53d00fd2848p-4)));
    try std.testing.expectEqual(0x5.95f3ec4683fap-4, log1p(@as(f80, 0x6.af53d00fd2844p-4)));
    try std.testing.expectEqual(0x5.95f3ec4683fa14ap-4, log1p(@as(f80, 0x6.af53d00fd2845d48p-4)));
    try std.testing.expectEqual(0x5.95f3ec4683fa14ap-4, log1p(@as(f80, 0x6.af53d00fd2845d4p-4)));

    try std.testing.expectEqual(0x0p+0, log1p(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, log1p(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1.00000039ece11db67b8f96c29c56p+0, log1p(@as(f128, 0x1.b7e152p+0)));
    try std.testing.expectEqual(0xf.fffff7d922f51a2d208d1c4e821p-4, log1p(@as(f128, 0x1.b7e15p+0)));
    try std.testing.expectEqual(0x1.000000000000020dcae0c29f344ep+0, log1p(@as(f128, 0x1.b7e151628aed3p+0)));
    // try std.testing.expectEqual(0xf.ffffffffffffc2af553376366578p-4, log1p(@as(f128, 0x1.b7e151628aed2p+0)));
    try std.testing.expectEqual(0x1.000000000000000075ed29d49ac4p+0, log1p(@as(f128, 0x1.b7e151628aed2a6cp+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffb9927823334ap-4, log1p(@as(f128, 0x1.b7e151628aed2a6ap+0)));
    // try std.testing.expectEqual(0x1p+0, log1p(@as(f128, 0x1.b7e151628aed2a6abf7158809cf5p+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffff8p-4, log1p(@as(f128, 0x1.b7e151628aed2a6abf7158809cf4p+0)));
    try std.testing.expectEqual(0x1.0000000000000000000000000004p+0, log1p(@as(f128, 0x1.b7e151628aed2a6abf7158809dp+0)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffd5p-4, log1p(@as(f128, 0x1.b7e151628aed2a6abf7158809c8p+0)));
    try std.testing.expectEqual(-0x4.9a58844d36e49e0efadd9db02aa8p-4, log1p(@as(f128, -0x4p-4)));
    try std.testing.expectEqual(-0x2.145647e7756e6d035dab1ac80bd8p+0, log1p(@as(f128, -0xep-4)));
    try std.testing.expectEqual(0x7.e0a6c39e0cc0133e3f04f1ef22ap-8, log1p(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(0x3.ff8015515621f7809a0a32499268p-12, log1p(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffe0002aaa6aab111066678af6bp-16, log1p(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0xf.ffff800005555515555888885dep-24, log1p(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0x7.fffffe000000aaaaaa6aaaaac444p-28, log1p(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x3.fffffff800000015555555155556p-32, log1p(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0x1.ffffffffe000000002aaaaaaaa6bp-36, log1p(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0xf.fffffffff8000000000555555558p-44, log1p(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0x7.ffffffffffe00000000000aaaaacp-48, log1p(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0x3.ffffffffffff8000000000001556p-52, log1p(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x1.fffffffffffffe00000000000003p-56, log1p(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0xf.ffffffffffffff8p-64, log1p(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffff8p-104, log1p(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, log1p(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, log1p(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, log1p(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, log1p(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, log1p(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0x4p-128, log1p(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, log1p(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, log1p(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, log1p(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, log1p(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x8p-152, log1p(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, log1p(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, log1p(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, log1p(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, log1p(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x4p-128, log1p(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, log1p(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, log1p(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, log1p(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, log1p(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x8p-152, log1p(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, log1p(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, log1p(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, log1p(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, log1p(@as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x6.eeb4e7af873022d55aeea4934c18p+0, log1p(@as(f128, 0x4p+8)));
    try std.testing.expectEqual(0xd.dce9ef5c63b581711b1f1d35a78p+0, log1p(@as(f128, 0x1p+20)));
    try std.testing.expectEqual(0x1.4cb5ecf0e96504219a8af0be5cbdp+4, log1p(@as(f128, 0x4p+28)));
    try std.testing.expectEqual(0x2.2a848ae66fa86038d6e7913b0bd8p+4, log1p(@as(f128, 0x4p+48)));
    try std.testing.expectEqual(0x2.996bd9e152ca08453515e17a0edp+4, log1p(@as(f128, 0x1p+60)));
    try std.testing.expectEqual(0x4.550915ccdf50b871adcf227619bp+4, log1p(@as(f128, 0x1p+100)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bc55e4f18476ac644p+4, log1p(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.b525ada00b9273470ca17589cf6ep+8, log1p(@as(f128, 0x1p+1000)));
    // try std.testing.expectEqual(0x5.8bb34ffb3a0b27cf845e6ca8643cp-4, log1p(@as(f128, 0x6.a0cf48p-4)));
    // try std.testing.expectEqual(0x5.8bb34a531ea8cdd0163aea3335dp-4, log1p(@as(f128, 0x6.a0cf4p-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0d19957c075571b6p-4, log1p(@as(f128, 0x6.a0cf42befceap-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0a4587cab2360533p-4, log1p(@as(f128, 0x6.a0cf42befce9cp-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0c45cffac98fe0dap-4, log1p(@as(f128, 0x6.a0cf42befce9ed48p-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0c45757913653cec8p-4, log1p(@as(f128, 0x6.a0cf42befce9ed4p-4)));
    // try std.testing.expectEqual(0x5.8bb34c4430e0c457b6453ad6d01cp-4, log1p(@as(f128, 0x6.a0cf42befce9ed4085ef59254b48p-4)));
    // try std.testing.expectEqual(0x5.8bb34c4430e0c457b6453ad6d0ap-4, log1p(@as(f128, 0x6.a0cf42befce9ed4085ef59254cp-4)));
    try std.testing.expectEqual(0x5.8bb34c4430e0c457b6453ad6cf34p-4, log1p(@as(f128, 0x6.a0cf42befce9ed4085ef59254ap-4)));
    try std.testing.expectEqual(0x5.8b90bfae8e7bc55e4f18476ac644p+4, log1p(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6a7278ece600fccp+8, log1p(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6af277ece600fccp+12, log1p(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6af278ece600fccp+12, log1p(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6ab278ece600fccp+8, log1p(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    // try std.testing.expectEqual(0x5.ebc1a69570c135d91b8705842a88p-4, log1p(@as(f128, 0x7.2a4368p-4)));
    try std.testing.expectEqual(0x5.af7a38286eaf6f0b3c7b1e4ccd4p-4, log1p(@as(f128, 0x6.d3a118p-4)));
    try std.testing.expectEqual(0x1.cb58e45e6b3a48e81d05dc4c71ebp+0, log1p(@as(f128, 0x5.03f228p+0)));
    try std.testing.expectEqual(0x5.e90249b494e60a8bfa1c60cc5b64p-4, log1p(@as(f128, 0x7.264968p-4)));
    // try std.testing.expectEqual(0x5.e902442d1717cc617b27353b5a34p-4, log1p(@as(f128, 0x7.26496p-4)));
    // try std.testing.expectEqual(0x5.e902469e458c38810487acf31074p-4, log1p(@as(f128, 0x7.264963888ac9p-4)));
    try std.testing.expectEqual(0x6.cc4a0b2426577f775c5d202625dcp-4, log1p(@as(f128, 0x8.786bdp-4)));
    try std.testing.expectEqual(0x6.2d4201bc5b6462c2eff2c81202a8p-4, log1p(@as(f128, 0x7.89dc18p-4)));
    // try std.testing.expectEqual(0x6.2d41fc4c416bbd61263d08eac78cp-4, log1p(@as(f128, 0x7.89dc1p-4)));
    try std.testing.expectEqual(0x6.2d420160a1d2435a5e14dad81f5p-4, log1p(@as(f128, 0x7.89dc17790eeb4p-4)));
    try std.testing.expectEqual(0x7.763bb83b8a4b8cd23140303f117cp-4, log1p(@as(f128, 0x9.81cdp-4)));
    // try std.testing.expectEqual(0x7.763bae3235afd3caf2f894c59a44p-4, log1p(@as(f128, 0x9.81ccfp-4)));
    try std.testing.expectEqual(0x7.763bb38c7d2d29ceaa259902b128p-4, log1p(@as(f128, 0x9.81ccf8887c25p-4)));
    try std.testing.expectEqual(0x7.763bb38c7d2cd98405497e6a1354p-4, log1p(@as(f128, 0x9.81ccf8887c248p-4)));
    try std.testing.expectEqual(0x7.763bb38c7d2cf2692a273d29257p-4, log1p(@as(f128, 0x9.81ccf8887c24a7bp-4)));
    // try std.testing.expectEqual(0x7.f5ac5ef8280eaea7ea32692292dcp-4, log1p(@as(f128, 0xa.50287p-4)));
    try std.testing.expectEqual(0x7.f5ac553d891675c6d6d447c2c1fp-4, log1p(@as(f128, 0xa.50286p-4)));
    try std.testing.expectEqual(0x7.f5ac559290b5bb49e382a3dd5384p-4, log1p(@as(f128, 0xa.5028608bd65f8p-4)));
    try std.testing.expectEqual(0x7.f5ac559290b56d74ebaad0a6f5dp-4, log1p(@as(f128, 0xa.5028608bd65fp-4)));
    // try std.testing.expectEqual(0x7.f5ac559290b59000922bfbc5f754p-4, log1p(@as(f128, 0xa.5028608bd65f38dp-4)));
    try std.testing.expectEqual(0x4.e92f6c6a921d11d700988efa8554p-4, log1p(@as(f128, 0x5.bf7888p-4)));
    try std.testing.expectEqual(0x4.e92f6687da531ad3bbd2a84c61dcp-4, log1p(@as(f128, 0x5.bf788p-4)));
    try std.testing.expectEqual(0x4.e92f6bdbf037127996e9f3eaf5ep-4, log1p(@as(f128, 0x5.bf78873e20a3p-4)));
    try std.testing.expectEqual(0x4.e92f6bdbf036e363d8a141112008p-4, log1p(@as(f128, 0x5.bf78873e20a2cp-4)));
    try std.testing.expectEqual(0x4.e92f6bdbf036f267276d2f17d4e4p-4, log1p(@as(f128, 0x5.bf78873e20a2d468p-4)));
    try std.testing.expectEqual(0x6.43431f9d697f9599160da109136p-4, log1p(@as(f128, 0x7.aa5198p-4)));
    try std.testing.expectEqual(0x1.34829b3156a22800c5c373313df7p+0, log1p(@as(f128, 0x2.564fap+0)));
    try std.testing.expectEqual(0x6.7a363b8f2519bbdbff874d3e2684p-4, log1p(@as(f128, 0x7.fc243p-4)));
    try std.testing.expectEqual(0x6.7a363638f41ca3d81d86d43434e8p-4, log1p(@as(f128, 0x7.fc2428p-4)));
    try std.testing.expectEqual(0x6.7a3637a5521e279ce027bbaab964p-4, log1p(@as(f128, 0x7.fc242a2235224p-4)));
    try std.testing.expectEqual(0x6.7a3637a5521dfceb583ba9119e34p-4, log1p(@as(f128, 0x7.fc242a223522p-4)));
    try std.testing.expectEqual(0x6.7a3637a5521e1c4061dc0937bb14p-4, log1p(@as(f128, 0x7.fc242a2235222ef8p-4)));
    try std.testing.expectEqual(-0x5.ec9649184d39a5d811e46ea6c808p-4, log1p(@as(f128, -0x4.f37d38p-4)));
    try std.testing.expectEqual(-0x5.ec9654ae08e379ad04f8cef1c53cp-4, log1p(@as(f128, -0x4.f37d4p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a2d5fb9ef23475bf8p-4, log1p(@as(f128, -0x4.f37d3c9ce0b14p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a8a0d9742e5d7514p-4, log1p(@as(f128, -0x4.f37d3c9ce0b18p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a3e8666c48368bfd4p-4, log1p(@as(f128, -0x4.f37d3c9ce0b14bd8p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a3e91fc802de111d4p-4, log1p(@as(f128, -0x4.f37d3c9ce0b14bep-4)));
    // try std.testing.expectEqual(-0x5.ec964fc6583a3e8e67ba42e3763p-4, log1p(@as(f128, -0x4.f37d3c9ce0b14bdd86eb157df5d4p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a3e8e67ba42e37388p-4, log1p(@as(f128, -0x4.f37d3c9ce0b14bdd86eb157df4p-4)));
    try std.testing.expectEqual(-0x5.ec964fc6583a3e8e67ba42e3767p-4, log1p(@as(f128, -0x4.f37d3c9ce0b14bdd86eb157df6p-4)));
    try std.testing.expectEqual(0x5.eee1d129eb6330c8efb7fbb2898p-4, log1p(@as(f128, 0x7.2eca58p-4)));
    try std.testing.expectEqual(0x5.eee1cba474cc2b4424c0cec64f04p-4, log1p(@as(f128, 0x7.2eca5p-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c3b3575d37686d98p-4, log1p(@as(f128, 0x7.2eca50c4d931cp-4)));
    // try std.testing.expectEqual(0x5.eee1cc2c508c0f09c115226258a4p-4, log1p(@as(f128, 0x7.2eca50c4d9318p-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c1e5f6d74b5d884ecp-4, log1p(@as(f128, 0x7.2eca50c4d9319638p-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c1e59e7fe1e0e005cp-4, log1p(@as(f128, 0x7.2eca50c4d931963p-4)));
    // try std.testing.expectEqual(0x5.eee1cc2c508c1e5e29fb19e8d748p-4, log1p(@as(f128, 0x7.2eca50c4d93196362b4f37f6e8dcp-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c1e5e29fb19e8d814p-4, log1p(@as(f128, 0x7.2eca50c4d93196362b4f37f6eap-4)));
    try std.testing.expectEqual(0x5.eee1cc2c508c1e5e29fb19e8d6bp-4, log1p(@as(f128, 0x7.2eca50c4d93196362b4f37f6e8p-4)));
    try std.testing.expectEqual(-0x7.ecba94fcebfbf0a3b365270e3358p-4, log1p(@as(f128, -0x6.3fef3p-4)));
    try std.testing.expectEqual(-0x7.ecbaa21da76b5866a6b58952da98p-4, log1p(@as(f128, -0x6.3fef38p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e861fde17bc7497e9ap-4, log1p(@as(f128, -0x6.3fef3067427e4p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e8688e3f310f51f1538p-4, log1p(@as(f128, -0x6.3fef3067427e8p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e86262cd1c8a7d007ccp-4, log1p(@as(f128, -0x6.3fef3067427e43d8p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e862639f284126018bp-4, log1p(@as(f128, -0x6.3fef3067427e43ep-4)));
    // try std.testing.expectEqual(-0x7.ecba95a65e862639a05337fb87ecp-4, log1p(@as(f128, -0x6.3fef3067427e43dfcde9e48f74bcp-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e862639a05337fb86b8p-4, log1p(@as(f128, -0x6.3fef3067427e43dfcde9e48f74p-4)));
    try std.testing.expectEqual(-0x7.ecba95a65e862639a05337fb8ap-4, log1p(@as(f128, -0x6.3fef3067427e43dfcde9e48f76p-4)));
    try std.testing.expectEqual(0x5.95f3f1dfd7b5e044101fab82d1fp-4, log1p(@as(f128, 0x6.af53d8p-4)));
    // try std.testing.expectEqual(0x5.95f3ec3b5b15423c0cde1ea58e4p-4, log1p(@as(f128, 0x6.af53dp-4)));
    try std.testing.expectEqual(0x5.95f3ec4683fa2d2090ee94cc945p-4, log1p(@as(f128, 0x6.af53d00fd2848p-4)));
    try std.testing.expectEqual(0x5.95f3ec4683f9fffcabe1ce64f9f4p-4, log1p(@as(f128, 0x6.af53d00fd2844p-4)));
    // try std.testing.expectEqual(0x5.95f3ec4683fa14a3b80d46ab2084p-4, log1p(@as(f128, 0x6.af53d00fd2845d48p-4)));
    try std.testing.expectEqual(0x5.95f3ec4683fa149e1390a512539p-4, log1p(@as(f128, 0x6.af53d00fd2845d4p-4)));
    // try std.testing.expectEqual(0x5.95f3ec4683fa14a354007a53e9f8p-4, log1p(@as(f128, 0x6.af53d00fd2845d4772260ef5adc4p-4)));
    // try std.testing.expectEqual(0x5.95f3ec4683fa14a354007a53ea2p-4, log1p(@as(f128, 0x6.af53d00fd2845d4772260ef5aep-4)));
    // try std.testing.expectEqual(0x5.95f3ec4683fa14a354007a53e8b8p-4, log1p(@as(f128, 0x6.af53d00fd2845d4772260ef5acp-4)));
}
