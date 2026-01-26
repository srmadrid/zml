const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Log(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.log: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

/// Returns the natural logarithm $\log(x)$ of a float, int or bool operand. The
/// result type is determined by coercing the operand type to a float, and the
/// operation is performed by casting the operand to the result type, then
/// computing its logarithm.
///
/// ## Signature
/// ```zig
/// float.log(x: X) float.Log(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the logarithm of.
///
/// ## Returns
/// `float.Log(@TypeOf(x))`: The natural logarithm of `x`.
pub inline fn log(x: anytype) float.Log(@TypeOf(x)) {
    switch (float.Log(@TypeOf(x))) {
        f16 => return types.scast(f16, log32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_logf.c
            return log32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_log.c
            return log64(types.scast(f64, x));
        },
        f80 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld80/e_logl.c
            return log80(types.scast(f80, x));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_logl.c
            return log128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_logf.c
//
// Original copyright notice:
// e_logf.c -- float version of e_log.c.
// Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn log32(x: f32) f32 {
    var ix: i32 = @bitCast(x);
    var k: i32 = 0;
    var xx: f32 = x;
    if (ix < 0x00800000) { // x < 2**-126
        if (ix & 0x7fffffff == 0)
            return -std.math.inf(f32);

        if (ix < 0)
            return std.math.nan(f32);

        // Subnormal number, scale up xx
        k -%= 25;
        xx *= 3.355443200e+07;
        ix = @bitCast(xx);
    }

    if (ix >= 0x7f800000)
        return x + x;

    k += (ix >> 23) -% 127;
    ix &= 0x007fffff;
    var i: i32 = (ix +% (0x95f64 << 3)) & 0x800000;
    xx = @bitCast(ix | (i ^ 0x3f800000)); // Normalize xx or xx/2
    k +%= i >> 23;
    const f: f32 = xx - 1.0;
    if ((0x007fffff & (0x8000 +% ix)) < 0xc000) { // -2**-9 <= f < 2**-9
        if (f == 0.0) {
            if (k == 0) {
                return 0.0;
            } else {
                const dk: f32 = types.scast(f32, k);
                return dk * 6.9313812256e-01 + dk * 9.0580006145e-06;
            }
        }

        const R: f32 = f * f * (0.5 - (1.0 / 3.0) * f);
        if (k == 0) {
            return f - R;
        } else {
            const dk: f32 = types.scast(f32, k);
            return dk * 6.9313812256e-01 - ((R - dk * 9.0580006145e-06) - f);
        }
    }

    const s: f32 = f / (2.0 + f);
    const dk: f32 = types.scast(f32, k);
    const z: f32 = s * s;
    i = ix -% (0x6147a << 3);
    const w: f32 = z * z;
    const j: i32 = (0x6b851 << 3) -% ix;
    const t1: f32 = w * (0xccce13.0p-25 + w * 0xf89e26.0p-26);
    const t2: f32 = z * (0xaaaaaa.0p-24 + w * 0x91e9ee.0p-25);
    i |= j;
    const R: f32 = t2 + t1;
    if (i > 0) {
        const hfsq: f32 = 0.5 * f * f;
        if (k == 0)
            return f - (hfsq - s * (hfsq + R))
        else
            return dk * 6.9313812256e-01 - ((hfsq - (s * (hfsq + R) + dk * 9.0580006145e-06)) - f);
    } else {
        if (k == 0)
            return f - s * (f - R)
        else
            return dk * 6.9313812256e-01 - ((s * (f - R) - dk * 9.0580006145e-06) - f);
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_log.c
//
// Original copyright notice:
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn log64(x: f64) f64 {
    var hx: i32 = @bitCast(dbl64.getHighPart(x));
    const lx: i32 = @bitCast(dbl64.getLowPart(x));
    var k: i32 = 0;
    var xx: f64 = x;
    if (hx < 0x00100000) { // x < 2**-1022
        if ((hx & 0x7fffffff) | lx == 0)
            return -std.math.inf(f64);

        if (hx < 0)
            return std.math.nan(f64);

        // Subnormal number, scale up xx
        k -%= 54;
        xx *= 1.80143985094819840000e+16;
        hx = @bitCast(dbl64.getHighPart(xx));
    }

    if (hx >= 0x7ff00000)
        return x + x;

    k += (hx >> 20) -% 1023;
    hx &= 0x000fffff;
    var i: i32 = (hx +% 0x95f64) & 0x100000;
    dbl64.setHighPart(&xx, @bitCast(hx | (i ^ 0x3ff00000))); // Normalize xx or xx/2
    k +%= i >> 20;
    const f: f64 = xx - 1.0;
    if ((0x000fffff & (2 +% hx)) < 3) { // -2**-20 <= f < 2**-20
        if (f == 0.0) {
            if (k == 0) {
                return 0.0;
            } else {
                const dk: f64 = types.scast(f64, k);
                return dk * 6.93147180369123816490e-1 + dk * 1.90821492927058770002e-10;
            }
        }

        const R: f64 = f * f * (0.5 - (1.0 / 3.0) * f);
        if (k == 0) {
            return f - R;
        } else {
            const dk: f64 = types.scast(f64, k);
            return dk * 6.93147180369123816490e-1 - ((R - dk * 1.90821492927058770002e-10) - f);
        }
    }

    const s: f64 = f / (2.0 + f);
    const dk: f64 = types.scast(f64, k);
    const z: f64 = s * s;
    i = hx -% 0x6147a;
    const w: f64 = z * z;
    const j: i32 = 0x6b851 -% hx;
    const t1: f64 = w *
        (3.999999999940941908e-1 + w *
            (2.222219843214978396e-1 + w *
                1.531383769920937332e-1));
    const t2: f64 = z *
        (6.666666666666735130e-1 + w *
            (2.857142874366239149e-1 + w *
                (1.818357216161805012e-1 + w *
                    1.479819860511658591e-1)));
    i |= j;
    const R: f64 = t2 + t1;
    if (i > 0) {
        const hfsq: f64 = 0.5 * f * f;
        if (k == 0)
            return f - (hfsq - s * (hfsq + R))
        else
            return dk * 6.93147180369123816490e-1 - ((hfsq - (s * (hfsq + R) + dk * 1.90821492927058770002e-10)) - f);
    } else {
        if (k == 0)
            return f - s * (f - R)
        else
            return dk * 6.93147180369123816490e-1 - ((s * (f - R) - dk * 1.90821492927058770002e-10) - f);
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld80/e_logl.c
//
// Original copyright notice:
// Copyright (c) 2008 Stephen L. Moshier <steve@moshier.net>
//
// Permission to use, copy, modify, and distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
fn log80(x: f80) f80 {
    if (std.math.isNan(x))
        return x;

    if (std.math.isInf(x))
        return x;

    // Test for domain
    if (x <= 0.0) {
        if (x == 0.0)
            return -std.math.inf(f80)
        else
            return std.math.nan(f80);
    }

    // Separate mantissa and exponent
    var e: i32 = undefined;
    var xx: f80 = float.frexp(x, &e);
    var z: f80 = undefined;
    var y: f80 = undefined;
    if (e > 2 or e < -2) {
        if (xx < 0.70710678118654752440) {
            e -%= 1;
            z = xx - 0.5;
            y = 0.5 * z + 0.5;
        } else {
            z = xx - 0.5;
            z -= 0.5;
            y = 0.5 * xx + 0.5;
        }

        xx = z / y;
        z = xx * xx;
        z = xx * (z *
            (-3.5717684488096787370998e1 + z *
                (1.0777257190312272158094e1 + z *
                    (-7.1990767473014147232598e-1 + z *
                        1.9757429581415468984296e-3))) /
            (-4.2861221385716144629696e2 + z *
                (1.9361891836232102174846e2 + z *
                    (-2.6201045551331104417768e1 + z *
                        1.00000000000000000000e0))));
        z += types.scast(f80, e) * 1.4286068203094172321215e-6;
        z += xx;
        z += types.scast(f80, e) * 6.93145751953125e-1;
        return z;
    }

    if (xx < 0.70710678118654752440) {
        e -%= 1;
        xx = xx + xx - 1.0;
    } else {
        xx = xx - 1.0;
    }

    z = xx * xx;
    y = xx * (z *
        (2.0039553499201281259648e1 + xx *
            (5.7112963590585538103336e1 + xx *
                (6.0949667980987787057556e1 + xx *
                    (2.9911919328553073277375e1 + xx *
                        (6.5787325942061044846969e0 + xx *
                            (4.9854102823193375972212e-1 + xx *
                                4.5270000862445199635215e-5)))))) /
        (6.0118660497603843919306e1 + xx *
            (2.1642788614495947685003e2 + xx *
                (3.0909872225312059774938e2 + xx *
                    (2.2176239823732856465394e2 + xx *
                        (8.3047565967967209469434e1 + xx *
                            (1.5062909083469192043167e1 + xx *
                                1.0000000000000000000000e0)))))));
    y += types.scast(f80, e) * 1.4286068203094172321215e-6;
    z = y - float.ldexp(z, -1);
    z += xx;
    z += types.scast(f80, e) * 6.93145751953125e-1;
    return z;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_logl.c
//
// Original copyright notice:
// Copyright (c) 2008 Stephen L. Moshier <steve@moshier.net>
//
// Permission to use, copy, modify, and distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
fn log128(x: f128) f128 {
    var u: ldbl128.Parts32 = .fromFloat(x);
    var m: u32 = u.mswhi;

    // Check for IEEE special cases
    var k: i32 = types.scast(i32, m & 0x7fffffff);
    if ((@as(u32, @bitCast(k)) | u.mswlo | u.lswhi | u.lswlo) == 0)
        return -std.math.inf(f128);

    if (m & 0x80000000 != 0)
        return std.math.nan(f128);

    if (k >= 0x7fff0000)
        return x + x;

    var e: i32 = undefined;
    u = .fromFloat(float.frexp(x, &e));
    m = u.mswhi & 0xffff;
    m |= 0x10000;

    // Find lookup table index k from high order bits of the significand
    var t: ldbl128.Parts32 = undefined;
    if (m < 0x16800) {
        k = types.scast(i32, (m -% 0xff00) >> 9);
        t.mswhi = types.scast(u32, 0x3fff0000 +% (k << 9));
        t.mswlo = 0;
        t.lswhi = 0;
        t.lswlo = 0;
        u.mswhi +%= 0x10000;
        e -%= 1;
        k +%= 64;
    } else {
        k = types.scast(i32, (m -% 0xfe00) >> 10);
        t.mswhi = types.scast(u32, 0x3ffe0000 +% (k << 10));
        t.mswlo = 0;
        t.lswhi = 0;
        t.lswlo = 0;
    }

    // On this interval the table is not used due to cancellation error
    var z: f128 = undefined;
    if (x <= 1.0078125 and x >= 0.9921875) {
        z = x - 1.0;
        k = 64;
        t = .fromFloat(1.0);
        e = 0;
    } else {
        z = (u.toFloat() - t.toFloat()) / t.toFloat();
    }

    const w: f128 = z * z;
    var y: f128 = ((((((((((((6.668057591071739754844678883223432347481e-2 *
        z + -7.144242754190814657241902218399056829264e-2) *
        z + 7.692307559897661630807048686258659316091e-2) *
        z + -8.333333211818065121250921925397567745734e-2) *
        z + 9.090909090915566247008015301349979892689e-2) *
        z + -1.000000000000532974938900317952530453248e-1) *
        z + 1.111111111111111093947834982832456459186e-1) *
        z + -1.249999999999999987884655626377588149000e-1) *
        z + 1.428571428571428571428808945895490721564e-1) *
        z + -1.666666666666666666666798448356171665678e-1) *
        z + 1.999999999999999999999999998515277861905e-1) *
        z + -2.499999999999999999999999999486853077002e-1) *
        z + 3.333333333333333333333333333333336096926e-1) *
        z * w;
    y -= 0.5 * w;
    y += types.scast(f128, e) * 1.4286068203094172321214581765680755001344e-6;
    y += z;
    y += logtbl_128[types.scast(u32, k - 26)];
    y += t.toFloat() - 1.0;
    y += types.scast(f128, e) * 6.93145751953125e-1;
    return y;
}

const logtbl_128: [92]f128 = .{
    -5.5345593589352099112142921677820359632418e-2,
    -5.2108257402767124761784665198737642086148e-2,
    -4.8991686870576856279407775480686721935120e-2,
    -4.5993270766361228596215288742353061431071e-2,
    -4.3110481649613269682442058976885699556950e-2,
    -4.0340872319076331310838085093194799765520e-2,
    -3.7682072451780927439219005993827431503510e-2,
    -3.5131785416234343803903228503274262719586e-2,
    -3.2687785249045246292687241862699949178831e-2,
    -3.0347913785027239068190798397055267411813e-2,
    -2.8110077931525797884641940838507561326298e-2,
    -2.5972247078357715036426583294246819637618e-2,
    -2.3932450635346084858612873953407168217307e-2,
    -2.1988775689981395152022535153795155900240e-2,
    -2.0139364778244501615441044267387667496733e-2,
    -1.8382413762093794819267536615342902718324e-2,
    -1.6716169807550022358923589720001638093023e-2,
    -1.5138929457710992616226033183958974965355e-2,
    -1.3649036795397472900424896523305726435029e-2,
    -1.2244881690473465543308397998034325468152e-2,
    -1.0924898127200937840689817557742469105693e-2,
    -9.6875626072830301572839422532631079809328e-3,
    -8.5313926245226231463436209313499745894157e-3,
    -7.4549452072765973384933565912143044991706e-3,
    -6.4568155251217050991200599386801665681310e-3,
    -5.5356355563671005131126851708522185605193e-3,
    -4.6900728132525199028885749289712348829878e-3,
    -3.9188291218610470766469347968659624282519e-3,
    -3.2206394539524058873423550293617843896540e-3,
    -2.5942708080877805657374888909297113032132e-3,
    -2.0385211375711716729239156839929281289086e-3,
    -1.5522183228760777967376942769773768850872e-3,
    -1.1342191863606077520036253234446621373191e-3,
    -7.8340854719967065861624024730268350459991e-4,
    -4.9869831458030115699628274852562992756174e-4,
    -2.7902661731604211834685052867305795169688e-4,
    -1.2335696813916860754951146082826952093496e-4,
    -3.0677461025892873184042490943581654591817e-5,
    0.0000000000000000000000000000000000000000e0,
    -3.0359557945051052537099938863236321874198e-5,
    -1.2081346403474584914595395755316412213151e-4,
    -2.7044071846562177120083903771008342059094e-4,
    -4.7834133324631162897179240322783590830326e-4,
    -7.4363569786340080624467487620270965403695e-4,
    -1.0654639687057968333207323853366578860679e-3,
    -1.4429854811877171341298062134712230604279e-3,
    -1.8753781835651574193938679595797367137975e-3,
    -2.3618380914922506054347222273705859653658e-3,
    -2.9015787624124743013946600163375853631299e-3,
    -3.4938307889254087318399313316921940859043e-3,
    -4.1378413103128673800485306215154712148146e-3,
    -4.8328735414488877044289435125365629849599e-3,
    -5.5782063183564351739381962360253116934243e-3,
    -6.3731336597098858051938306767880719015261e-3,
    -7.2169643436165454612058905294782949315193e-3,
    -8.1090214990427641365934846191367315083867e-3,
    -9.0486422112807274112838713105168375482480e-3,
    -1.0035177140880864314674126398350812606841e-2,
    -1.1067990155502102718064936259435676477423e-2,
    -1.2146457974158024928196575103115488672416e-2,
    -1.3269969823361415906628825374158424754308e-2,
    -1.4437927104692837124388550722759686270765e-2,
    -1.5649743073340777659901053944852735064621e-2,
    -1.6904842527181702880599758489058031645317e-2,
    -1.8202661505988007336096407340750378994209e-2,
    -1.9542647000370545390701192438691126552961e-2,
    -2.0924256670080119637427928803038530924742e-2,
    -2.2346958571309108496179613803760727786257e-2,
    -2.3810230892650362330447187267648486279460e-2,
    -2.5313561699385640380910474255652501521033e-2,
    -2.6856448685790244233704909690165496625399e-2,
    -2.8438398935154170008519274953860128449036e-2,
    -3.0058928687233090922411781058956589863039e-2,
    -3.1717563112854831855692484086486099896614e-2,
    -3.3413836095418743219397234253475252001090e-2,
    -3.5147290019036555862676702093393332533702e-2,
    -3.6917475563073933027920505457688955423688e-2,
    -3.8723951502862058660874073462456610731178e-2,
    -4.0566284516358241168330505467000838017425e-2,
    -4.2444048996543693813649967076598766917965e-2,
    -4.4356826869355401653098777649745233339196e-2,
    -4.6304207416957323121106944474331029996141e-2,
    -4.8285787106164123613318093945035804818364e-2,
    -5.0301169421838218987124461766244507342648e-2,
    -5.2349964705088137924875459464622098310997e-2,
    -5.4431789996103111613753440311680967840214e-2,
    -5.6546268881465384189752786409400404404794e-2,
    -5.8693031345788023909329239565012647817664e-2,
    -6.0871713627532018185577188079210189048340e-2,
    -6.3081958078862169742820420185833800925568e-2,
    -6.5323413029406789694910800219643791556918e-2,
    -6.7595732653791419081537811574227049288168e-2,
};
