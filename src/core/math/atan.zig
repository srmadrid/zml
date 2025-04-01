const cast = @import("../types.zig").cast;
const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const dla = @import("dla.zig");
const atnat = @import("atnat.zig");
const EnsureFloat = types.EnsureFloat;

pub fn atan(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return atan(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, atan32(cast(f32, x))),
                f32 => return atan32(x),
                f64 => return atan64(x),
                f80 => return cast(f80, atan128(cast(f128, x, .{})), .{}),
                f128 => return atan128(x),
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn atan32(x: f32) f32 {
    const pi2: f64 = 0x1.921fb54442d18p+0;
    const t: u32 = @bitCast(x);
    const e: i32 = @bitCast((t >> 23) & 0xff);
    const gt: bool = e >= 127;
    const ta: u32 = t & 0x7fffffff;

    if (ta >= 0x4c700518) { // |x| > 0x1.e00a3p+25
        if (ta > 0x7f800000)
            return x + x; // nan

        return cast(f32, math.copysign(pi2, x), .{});
    }
    if (e < 127 - 13) {
        if (e < 127 - 25) {
            if ((t << 1) == 0)
                return x;

            return @mulAdd(f32, -x, @abs(x), x);
        }
        return @mulAdd(f32, -0x1.5555555555555p-2 * x, x * x, x);
    }

    // now |x| >= 0x1p-13
    var z: f64 = cast(f64, x, .{});
    if (gt)
        z = 1 / z; // gt is non-zero for |x| >= 1

    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    const z8: f64 = z4 * z4;

    // polynomials generated using rminimax
    // (https://gitlab.inria.fr/sfilip/rminimax) with the following command:
    // ./ratapprox --function="atan(x)" --dom=[0.000122070,1]
    // --num=[x,x^3,x^5,x^7,x^9,x^11,x^13] --den=[1,x^2,x^4,x^6,x^8,x^10,x^12]
    // --output=atanf.sollya --log (see output atanf.sollya) The coefficient
    // cd[0] was slightly reduced from the original value 0x1.51eccde075d67p-2 to
    // avoid an exceptional case for |x| = 0x1.1ad646p-4 and rounding to nearest.
    const cn: []const f64 =
        &.{
            0x1.51eccde075d67p-2, 0x1.a76bb5637f2f2p-1, 0x1.81e0eed20de88p-1,
            0x1.376c8ca67d11dp-2, 0x1.aec7b69202ac6p-5, 0x1.9561899acc73ep-9,
            0x1.bf9fa5b67e6p-16,
        };
    const cd: []const f64 =
        &.{ 0x1.51eccde075d66p-2, 0x1.dfbdd7b392d28p-1, 0x1p+0, 0x1.fd22bf0e89b54p-2, 0x1.d91ff8b576282p-4, 0x1.653ea99fc9bbp-7, 0x1.1e7fcc202340ap-12 };
    var cn0: f64 = cn[0] + z2 * cn[1];
    const cn2: f64 = cn[2] + z2 * cn[3];
    var cn4: f64 = cn[4] + z2 * cn[5];
    const cn6: f64 = cn[6];
    cn0 += z4 * cn2;
    cn4 += z4 * cn6;
    cn0 += z8 * cn4;
    cn0 *= z;
    var cd0: f64 = cd[0] + z2 * cd[1];
    const cd2: f64 = cd[2] + z2 * cd[3];
    var cd4: f64 = cd[4] + z2 * cd[5];
    const cd6: f64 = cd[6];
    cd0 += z4 * cd2;
    cd4 += z4 * cd6;
    cd0 += z8 * cd4;
    var r: f64 = cn0 / cd0;
    if (!gt)
        return cast(f32, r, .{}); // for |x| < 1, (float) r is correctly rounded

    const PI_OVER2_H: f64 = 0x1.9p0;
    const PI_OVER2_L: f64 = 0x1.0fdaa22168c23p-7;
    // now r approximates atan(1/x), we use atan(x) + atan(1/x) = sign(x)*pi/2,
    // where PI_OVER2_H + PI_OVER2_L approximates pi/2.
    // With sign(z)*L + (-r + sign(z)*H), it fails for x=0x1.98c252p+12 and
    // rounding upward.
    // With sign(z)*PI - r, where PI is a double approximation of pi to nearest,
    // it fails for x=0x1.ddf9f6p+0 and rounding upward. */
    r = math.copysign(PI_OVER2_L, z) - r + math.copysign(PI_OVER2_H, z);
    return cast(f32, r, .{});
}

fn atan64(x: f64) f64 {
    const TWO52: f64 = 0x1.0p52;

    const num: [2]u32 = @bitCast(x);
    const ux: u32 = num[atnat.HIGH_HALF];
    const dx: u32 = num[atnat.LOW_HALF];

    // x=NaN
    if (((ux & 0x7ff00000) == 0x7ff00000) and (((ux & 0x000fffff) | dx) != 0x00000000))
        return x + x;

    // Regular values of x, including denormals +-0 and +-INF
    const u: f64 = if (x < 0) -x else x;
    if (u < @as(f64, @bitCast(atnat.c))) {
        if (u < @as(f64, @bitCast(atnat.b))) {
            if (u < @as(f64, @bitCast(atnat.a))) {
                if (u < std.math.floatMin(f64)) {
                    const vu: f64 = if (u != 0) u else std.math.floatMin(f64);
                    const force_underflow = vu * vu;
                    std.mem.doNotOptimizeAway(force_underflow);
                }
                return x;
            } else { // a <= u < b
                const v: f64 = x * x;
                var yy: f64 = @as(f64, @bitCast(atnat.d11)) + v * @as(f64, @bitCast(atnat.d13));
                yy = @as(f64, @bitCast(atnat.d9)) + v * yy;
                yy = @as(f64, @bitCast(atnat.d7)) + v * yy;
                yy = @as(f64, @bitCast(atnat.d5)) + v * yy;
                yy = @as(f64, @bitCast(atnat.d3)) + v * yy;
                yy *= x * v;

                const y: f64 = x + yy;
                // Max ULP is 0.511.
                return y;
            }
        } else { // b <= u < c
            var i: i32 = cast(i32, (TWO52 + 256 * u) - TWO52, .{});
            i -= 16;
            const z: f64 = u - @as(f64, @bitCast(atnat.cij[@intCast(i)][0]));
            var yy: f64 = @as(f64, @bitCast(atnat.cij[@intCast(i)][5])) + z * @as(f64, @bitCast(atnat.cij[@intCast(i)][6]));
            yy = @as(f64, @bitCast(atnat.cij[@intCast(i)][4])) + z * yy;
            yy = @as(f64, @bitCast(atnat.cij[@intCast(i)][3])) + z * yy;
            yy = @as(f64, @bitCast(atnat.cij[@intCast(i)][2])) + z * yy;
            yy *= z;

            const t1: f64 = @as(f64, @bitCast(atnat.cij[@intCast(i)][1]));
            const y: f64 = t1 + yy;
            // Max ULP is 0.56.
            return math.copysign(y, x);
        }
    } else {
        if (u < @as(f64, @bitCast(atnat.d))) { // c <= u < d
            const w: f64 = 1 / u;
            var t1: f64 = undefined;
            var t2: f64 = undefined;
            dla.emulv(w, u, &t1, &t2);
            const ww: f64 = w * ((1 - t1) - t2);
            var i: i32 = cast(i32, (TWO52 + 256 * w) - TWO52, .{});
            i -= 16;
            const z: f64 = (w - @as(f64, @bitCast(atnat.cij[@intCast(i)][0]))) + ww;

            var yy: f64 = @as(f64, @bitCast(atnat.cij[@intCast(i)][5])) + z * @as(f64, @bitCast(atnat.cij[@intCast(i)][6]));
            yy = @as(f64, @bitCast(atnat.cij[@intCast(i)][4])) + z * yy;
            yy = @as(f64, @bitCast(atnat.cij[@intCast(i)][3])) + z * yy;
            yy = @as(f64, @bitCast(atnat.cij[@intCast(i)][2])) + z * yy;
            yy = @as(f64, @bitCast(atnat.hpi1)) - z * yy;

            t1 = @as(f64, @bitCast(atnat.hpi)) - @as(f64, @bitCast(atnat.cij[@intCast(i)][1]));
            const y: f64 = t1 + yy;
            // Max ULP is 0.503.

            return math.copysign(y, x);
        } else {
            if (u < @as(f64, @bitCast(atnat.e))) { // d <= u < e
                const w: f64 = 1 / u;
                const v: f64 = w * w;
                var t1: f64 = undefined;
                var t2: f64 = undefined;
                dla.emulv(w, u, &t1, &t2);

                var yy: f64 = @as(f64, @bitCast(atnat.d11)) + v * @as(f64, @bitCast(atnat.d13));
                yy = @as(f64, @bitCast(atnat.d9)) + v * yy;
                yy = @as(f64, @bitCast(atnat.d7)) + v * yy;
                yy = @as(f64, @bitCast(atnat.d5)) + v * yy;
                yy = @as(f64, @bitCast(atnat.d3)) + v * yy;
                yy *= w * v;

                const ww: f64 = w * ((1 - t1) - t2);
                var t3: f64 = undefined;
                var cor: f64 = undefined;
                dla.esub(@as(f64, @bitCast(atnat.hpi)), w, &t3, &cor);
                yy = ((@as(f64, @bitCast(atnat.hpi1)) + cor) - ww) - yy;
                const y: f64 = t3 + yy;
                // Max ULP is 0.5003.
                return math.copysign(y, x);
            } else {
                // u >= e
                if (x > 0) {
                    return @as(f64, @bitCast(atnat.hpi));
                } else {
                    return @as(f64, @bitCast(atnat.mhpi));
                }
            }
        }
    }
}

fn atan128(x: f128) f128 {
    const atantbl: [84]f128 = .{
        0.0000000000000000000000000000000000000000e0,
        1.2435499454676143503135484916387102557317e-1, // arctan(0.125)
        2.4497866312686415417208248121127581091414e-1,
        3.5877067027057222039592006392646049977698e-1,
        4.6364760900080611621425623146121440202854e-1,
        5.5859931534356243597150821640166127034645e-1,
        6.4350110879328438680280922871732263804151e-1,
        7.1882999962162450541701415152590465395142e-1,
        7.8539816339744830961566084581987572104929e-1,
        8.4415398611317100251784414827164750652594e-1,
        8.9605538457134395617480071802993782702458e-1,
        9.4200004037946366473793717053459358607166e-1,
        9.8279372324732906798571061101466601449688e-1,
        1.0191413442663497346383429170230636487744e0,
        1.0516502125483736674598673120862998296302e0,
        1.0808390005411683108871567292171998202703e0,
        1.1071487177940905030170654601785370400700e0,
        1.1309537439791604464709335155363278047493e0,
        1.1525719972156675180401498626127513797495e0,
        1.1722738811284763866005949441337046149712e0,
        1.1902899496825317329277337748293183376012e0,
        1.2068173702852525303955115800565576303133e0,
        1.2220253232109896370417417439225704908830e0,
        1.2360594894780819419094519711090786987027e0,
        1.2490457723982544258299170772810901230778e0,
        1.2610933822524404193139408812473357720101e0,
        1.2722973952087173412961937498224804940684e0,
        1.2827408797442707473628852511364955306249e0,
        1.2924966677897852679030914214070816845853e0,
        1.3016288340091961438047858503666855921414e0,
        1.3101939350475556342564376891719053122733e0,
        1.3182420510168370498593302023271362531155e0,
        1.3258176636680324650592392104284756311844e0,
        1.3329603993374458675538498697331558093700e0,
        1.3397056595989995393283037525895557411039e0,
        1.3460851583802539310489409282517796256512e0,
        1.3521273809209546571891479413898128509842e0,
        1.3578579772154994751124898859640585287459e0,
        1.3633001003596939542892985278250991189943e0,
        1.3684746984165928776366381936948529556191e0,
        1.3734007669450158608612719264449611486510e0,
        1.3780955681325110444536609641291551522494e0,
        1.3825748214901258580599674177685685125566e0,
        1.3868528702577214543289381097042486034883e0,
        1.3909428270024183486427686943836432060856e0,
        1.3948567013423687823948122092044222644895e0,
        1.3986055122719575950126700816114282335732e0,
        1.4021993871854670105330304794336492676944e0,
        1.4056476493802697809521934019958079881002e0,
        1.4089588955564736949699075250792569287156e0,
        1.4121410646084952153676136718584891599630e0,
        1.4152014988178669079462550975833894394929e0,
        1.4181469983996314594038603039700989523716e0,
        1.4209838702219992566633046424614466661176e0,
        1.4237179714064941189018190466107297503086e0,
        1.4263547484202526397918060597281265695725e0,
        1.4288992721907326964184700745371983590908e0,
        1.4313562697035588982240194668401779312122e0,
        1.4337301524847089866404719096698873648610e0,
        1.4360250423171655234964275337155008780675e0,
        1.4382447944982225979614042479354815855386e0,
        1.4403930189057632173997301031392126865694e0,
        1.4424730991091018200252920599377292525125e0,
        1.4444882097316563655148453598508037025938e0,
        1.4464413322481351841999668424758804165254e0,
        1.4483352693775551917970437843145232637695e0,
        1.4501726582147939000905940595923466567576e0,
        1.4519559822271314199339700039142990228105e0,
        1.4536875822280323362423034480994649820285e0,
        1.4553696664279718992423082296859928222270e0,
        1.4570043196511885530074841089245667532358e0,
        1.4585935117976422128825857356750737658039e0,
        1.4601391056210009726721818194296893361233e0,
        1.4616428638860188872060496086383008594310e0,
        1.4631064559620759326975975316301202111560e0,
        1.4645314639038178118428450961503371619177e0,
        1.4659193880646627234129855241049975398470e0,
        1.4672716522843522691530527207287398276197e0,
        1.4685896086876430842559640450619880951144e0,
        1.4698745421276027686510391411132998919794e0,
        1.4711276743037345918528755717617308518553e0,
        1.4723501675822635384916444186631899205983e0,
        1.4735431285433308455179928682541563973416e0, // arctan(10.25)
        1.5707963267948966192313216916397514420986e0, // pi/2
    };

    // arctan t = t + t^3 p(t^2) / q(t^2)
    // |t| <= 0.09375
    // peak relative error 5.3e-37

    const p0: f128 = -4.283708356338736809269381409828726405572e1;
    const p1: f128 = -8.636132499244548540964557273544599863825e1;
    const p2: f128 = -5.713554848244551350855604111031839613216e1;
    const p3: f128 = -1.371405711877433266573835355036413750118e1;
    const p4: f128 = -8.638214309119210906997318946650189640184e-1;
    const q0: f128 = 1.285112506901621042780814422948906537959e2;
    const q1: f128 = 3.361907253914337187957855834229672347089e2;
    const q2: f128 = 3.180448303864130128268191635189365331680e2;
    const q3: f128 = 1.307244136980865800160844625025280344686e2;
    const q4: f128 = 2.173623741810414221251136181221172551416e1;
    // const q5: f128 = 1.000000000000000000000000000000000000000e0;

    const huge: f128 = 1.0e4930;

    const s: u128 = @bitCast(x);
    var s32: [4]u32 = @bitCast(s);
    if (@import("builtin").cpu.arch.endian() == .little) {
        const t: [4]u32 = .{ s32[3], s32[2], s32[1], s32[0] };
        s32 = t;
    }
    const sign: bool = (s32[0] & 0x80000000) != 0;

    // Check for IEEE special cases.
    var k: i32 = @bitCast(s32[0] & 0x7fffffff);
    if (k >= 0x7fff0000) {
        // NaN.
        if (((k & 0xffff) | @as(i32, @bitCast(s32[1])) | @as(i32, @bitCast(s32[2])) | @as(i32, @bitCast(s32[3]))) != 0)
            return x + x;

        // Infinity.
        if (sign) {
            return -atantbl[83];
        } else {
            return atantbl[83];
        }
    }

    if (k <= 0x3fc50000) // |x| < 2**-58
    {
        if (@abs(x) < std.math.floatMin(f64)) {
            const vx: f128 = if (x != 0) x else std.math.floatMin(f128);
            const force_underflow = vx * vx;
            std.mem.doNotOptimizeAway(force_underflow);
        }

        // Raise inexact.
        if (huge + x > 0) {
            return x;
        }
    }

    if (k >= 0x40720000) // |x| > 2**115
    {
        // Saturate result to {-,+}pi/2.
        if (sign) {
            return -atantbl[83];
        } else {
            return atantbl[83];
        }
    }

    var xx: f128 = undefined;
    if (sign) {
        xx = -x;
    } else {
        xx = x;
    }

    var u: f128 = undefined;
    var t: f128 = undefined;
    if (k >= 0x40024800) // 10.25
    {
        k = 83;
        t = -1 / xx;
    } else {
        // Index of nearest table element.
        // Roundoff to integer is asymmetrical to avoid cancellation when t < 0
        // (cf. fdlibm).
        k = cast(i32, 8 * xx + 0.25, .{});
        u = 0.125 * cast(f128, k, .{});
        // Small arctan argument.
        t = (xx - u) / (1 + xx * u);
    }

    // Arctan of small argument t.
    u = t * t;
    const p: f128 = ((((p4 * u) + p3) * u + p2) * u + p1) * u + p0;
    const q: f128 = ((((u + q4) * u + q3) * u + q2) * u + q1) * u + q0;
    u = t * u * p / q + t;

    // arctan x = arctan u  +  arctan t
    u = atantbl[@intCast(k)] + u;
    if (sign) {
        return -u;
    } else {
        return u;
    }
}

test atan {
    try std.testing.expectEqual(0x1.921fb6p+0, atan(std.math.inf(f32)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan(-std.math.inf(f32)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atan(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0xc.90fdbp-4, atan(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0xa.4bc7dp-4, atan(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x7.ff557p-8, atan(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x3.ffffecp-12, atan(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x2p-16, atan(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0x1p-20, atan(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0x8p-28, atan(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x4p-32, atan(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, atan(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, atan(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, atan(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, atan(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, atan(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, atan(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1.30b6d8p+0, atan(@as(f32, 0x2.8p+0)));
    try std.testing.expectEqual(0x1.789bd2p+0, atan(@as(f32, 0xap+0)));
    try std.testing.expectEqual(0x1.921fa4p+0, atan(@as(f32, 0xf.424p+16)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan(@as(f32, 0x8p+28)));
    try std.testing.expectEqual(0x1p-100, atan(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x3.9ff7ep-4, atan(@as(f32, -0x3.b02d84p-4)));
    try std.testing.expectEqual(-0x3.348f08p-4, atan(@as(f32, -0x3.3fb708p-4)));
    try std.testing.expectEqual(-0x1.24c032p+0, atan(@as(f32, -0x2.3249ap+0)));
    try std.testing.expectEqual(-0xe.1832ep-4, atan(@as(f32, -0x1.363f46p+0)));
    try std.testing.expectEqual(-0x1.087838p+0, atan(@as(f32, -0x1.ad4c0ap+0)));
    try std.testing.expectEqual(-0x1.522f04p+0, atan(@as(f32, -0x3.eb8e18p+0)));
    try std.testing.expectEqual(0x1.476166p+0, atan(@as(f32, 0x3.53c188p+0)));
    try std.testing.expectEqual(-0xe.e9f01p-4, atan(@as(f32, -0x1.58c83p+0)));
    try std.testing.expectEqual(0x9.b0001p-4, atan(@as(f32, 0xb.133b9p-4)));
    try std.testing.expectEqual(0x4p-128, atan(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-128, atan(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, atan(@as(f32, -0x8p-152)));

    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan(std.math.inf(f64)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan(-std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atan(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0xa.4bc7d1934f708p-4, atan(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x7.ff556eea5d894p-8, atan(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0x3.ffffeaaaab778p-12, atan(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffffd55555p-16, atan(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffffaaa8p-24, atan(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0x7.ffffffffffff4p-28, atan(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x4p-32, atan(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, atan(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, atan(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, atan(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, atan(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, atan(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, atan(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1.30b6d796a4da8p+0, atan(@as(f64, 0x2.8p+0)));
    try std.testing.expectEqual(0x1.789bd2c160054p+0, atan(@as(f64, 0xap+0)));
    try std.testing.expectEqual(0x1.921fa47d4b30dp+0, atan(@as(f64, 0xf.424p+16)));
    try std.testing.expectEqual(0x1.921fb54242d18p+0, atan(@as(f64, 0x8p+28)));
    try std.testing.expectEqual(0x1p-100, atan(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, atan(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atan(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x3.9ff7e1f81370cp-4, atan(@as(f64, -0x3.b02d84p-4)));
    try std.testing.expectEqual(-0x3.348f092072332p-4, atan(@as(f64, -0x3.3fb708p-4)));
    try std.testing.expectEqual(-0x1.24c032fe9a703p+0, atan(@as(f64, -0x2.3249ap+0)));
    try std.testing.expectEqual(-0xe.1832df50b2398p-4, atan(@as(f64, -0x1.363f46p+0)));
    try std.testing.expectEqual(-0x1.0878377062daep+0, atan(@as(f64, -0x1.ad4c0ap+0)));
    try std.testing.expectEqual(-0x1.522f0408c0e8cp+0, atan(@as(f64, -0x3.eb8e18p+0)));
    // try std.testing.expectEqual(0x1.476165c27ab51p+0, atan(@as(f64, 0x3.53c188p+0)));
    try std.testing.expectEqual(-0xe.e9f00a57b1438p-4, atan(@as(f64, -0x1.58c83p+0)));
    try std.testing.expectEqual(0x9.b0000da23b9ep-4, atan(@as(f64, 0xb.133b9p-4)));
    try std.testing.expectEqual(0x4p-128, atan(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, atan(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, atan(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, atan(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, atan(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, atan(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, atan(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-152, atan(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atan(@as(f64, -0x4p-1076)));

    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan(std.math.inf(f80)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan(-std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atan(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0xa.4bc7d1934f70924p-4, atan(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x7.ff556eea5d892a1p-8, atan(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x3.ffffeaaaab77777p-12, atan(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffffd5555555cp-16, atan(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffffaaaaabp-24, atan(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0x7.ffffffffffff5558p-28, atan(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x3.ffffffffffffffecp-32, atan(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, atan(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, atan(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, atan(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, atan(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, atan(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, atan(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x1.30b6d796a4da858ap+0, atan(@as(f80, 0x2.8p+0)));
    try std.testing.expectEqual(0x1.789bd2c16005382ep+0, atan(@as(f80, 0xap+0)));
    try std.testing.expectEqual(0x1.921fa47d4b30ce82p+0, atan(@as(f80, 0xf.424p+16)));
    try std.testing.expectEqual(0x1.921fb54242d1846ap+0, atan(@as(f80, 0x8p+28)));
    try std.testing.expectEqual(0x1p-100, atan(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, atan(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atan(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, atan(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(-0x3.9ff7e1f81370b93p-4, atan(@as(f80, -0x3.b02d84p-4)));
    try std.testing.expectEqual(-0x3.348f092072331fd8p-4, atan(@as(f80, -0x3.3fb708p-4)));
    try std.testing.expectEqual(-0x1.24c032fe9a702f72p+0, atan(@as(f80, -0x2.3249ap+0)));
    try std.testing.expectEqual(-0xe.1832df50b2398e5p-4, atan(@as(f80, -0x1.363f46p+0)));
    try std.testing.expectEqual(-0x1.0878377062dada2ap+0, atan(@as(f80, -0x1.ad4c0ap+0)));
    try std.testing.expectEqual(-0x1.522f0408c0e8c2d8p+0, atan(@as(f80, -0x3.eb8e18p+0)));
    try std.testing.expectEqual(0x1.476165c27ab518p+0, atan(@as(f80, 0x3.53c188p+0)));
    try std.testing.expectEqual(-0xe.e9f00a57b143b32p-4, atan(@as(f80, -0x1.58c83p+0)));
    try std.testing.expectEqual(0x9.b0000da23b9dce3p-4, atan(@as(f80, 0xb.133b9p-4)));
    try std.testing.expectEqual(0x4p-128, atan(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, atan(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, atan(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, atan(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, atan(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, atan(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, atan(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, atan(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, atan(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, atan(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, atan(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, atan(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-152, atan(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atan(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, atan(@as(f80, -0x8p-16448)));

    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan(std.math.inf(f128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan(-std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atan(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0xa.4bc7d1934f7092419a87f2a457d8p-4, atan(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x7.ff556eea5d892a13bcebbb6ed464p-8, atan(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(0x3.ffffeaaaab77776e52e5a019fbcep-12, atan(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffffd5555555bbbbbbba97297p-16, atan(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffffaaaaaaaaaadddddddep-24, atan(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0x7.ffffffffffff5555555555556efp-28, atan(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x3.ffffffffffffffeaaaaaaaaaaaaap-32, atan(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0x1.fffffffffffffffffd5555555555p-36, atan(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffaaaaaaaa8p-44, atan(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffff555554p-48, atan(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0x3.ffffffffffffffffffffffffeaaap-52, atan(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x1.fffffffffffffffffffffffffffdp-56, atan(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, atan(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0x1.30b6d796a4da8589532c0eec663ep+0, atan(@as(f128, 0x2.8p+0)));
    try std.testing.expectEqual(0x1.789bd2c16005382eabf0cd4b6aaep+0, atan(@as(f128, 0xap+0)));
    try std.testing.expectEqual(0x1.921fa47d4b30ce822275563fcb9ap+0, atan(@as(f128, 0xf.424p+16)));
    try std.testing.expectEqual(0x1.921fb54242d18469898cc519ac63p+0, atan(@as(f128, 0x8p+28)));
    try std.testing.expectEqual(0x1p-100, atan(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, atan(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, atan(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, atan(@as(f128, 0x1p-10000)));
    // try std.testing.expectEqual(-0x3.9ff7e1f81370b93195e0357d0b5ap-4, atan(@as(f128, -0x3.b02d84p-4)));
    try std.testing.expectEqual(-0x3.348f092072331fd8ca0cbff348d6p-4, atan(@as(f128, -0x3.3fb708p-4)));
    try std.testing.expectEqual(-0x1.24c032fe9a702f7255968f75e01cp+0, atan(@as(f128, -0x2.3249ap+0)));
    // try std.testing.expectEqual(-0xe.1832df50b2398e4a96945ef0f7f8p-4, atan(@as(f128, -0x1.363f46p+0)));
    try std.testing.expectEqual(-0x1.0878377062dada2af4f466e46577p+0, atan(@as(f128, -0x1.ad4c0ap+0)));
    try std.testing.expectEqual(-0x1.522f0408c0e8c2d8d3fe682ee8bdp+0, atan(@as(f128, -0x3.eb8e18p+0)));
    try std.testing.expectEqual(0x1.476165c27ab517ff156a94e45559p+0, atan(@as(f128, 0x3.53c188p+0)));
    try std.testing.expectEqual(-0xe.e9f00a57b143b31a8f4be18ea31p-4, atan(@as(f128, -0x1.58c83p+0)));
    try std.testing.expectEqual(0x9.b0000da23b9dce2cdb744dda364p-4, atan(@as(f128, 0xb.133b9p-4)));
    try std.testing.expectEqual(0x4p-128, atan(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, atan(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, atan(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, atan(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, atan(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, atan(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, atan(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, atan(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, atan(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, atan(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, atan(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, atan(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, atan(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, atan(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, atan(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-152, atan(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, atan(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, atan(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, atan(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, atan(@as(f128, -0x4p-16496)));
}
