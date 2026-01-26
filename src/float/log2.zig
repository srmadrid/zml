const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Log2(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.log2: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

/// Returns the base-2 logarithm $\log_2(x)$ of a float, int or bool operand.
/// The result type is determined by coercing the operand type to a float, and
/// the operation is performed by casting the operand to the result type, then
/// computing its base-2 logarithm.
///
/// ## Signature
/// ```zig
/// float.log2(x: X) float.Log2(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the base-2 logarithm of.
///
/// ## Returns
/// `float.Log2(@TypeOf(x))`: The base-2 logarithm of `x`.
pub inline fn log2(x: anytype) float.Log2(@TypeOf(x)) {
    switch (float.Log2(@TypeOf(x))) {
        f16 => return types.scast(f16, log2_32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_log2f.c
            return log2_32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_log2.c
            return log2_64(types.scast(f64, x));
        },
        f80 => {
            //
            // return log2_80(types.scast(f80, x));
            return types.scast(f80, log2_128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_log2l.c
            return log2_128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_log2f.c
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
fn log2_32(x: f32) f32 {
    var hx: i32 = @bitCast(x);

    var k: i32 = 0;
    var xx: f32 = x;
    if (hx < 0x00800000) { // x < 2**-126
        if ((hx & 0x7fffffff) == 0) // log(±0) = -inf
            return -std.math.inf(f32);

        if (hx < 0) // log(-#) = NaN
            return (x - x) / 0.0;

        k -= 25;
        xx *= 3.3554432000e+7; // Subnormal number, scale up x
        hx = @bitCast(xx);
    }

    if (hx >= 0x7f800000)
        return xx + xx;

    if (hx == 0x3f800000) // log(1) = +0
        return 0.0;

    k +%= (hx >> 23) -% 127;
    hx &= 0x007fffff;
    const i: i32 = (hx + (0x4afb0d)) & 0x800000;
    xx = @bitCast(hx | (i ^ 0x3f800000)); // Normalize x or x/2
    k +%= (i >> 23);
    const y: f32 = types.scast(f32, k);
    const f: f32 = xx - 1.0;
    const hfsq: f32 = 0.5 * f * f;
    const s: f32 = f / (2.0 + f);
    const z: f32 = s * s;
    const w: f32 = z * z;
    const t1: f32 = w * (0xccce13.0p-25 + w * 0xf89e26.0p-26);
    const t2: f32 = z * (0xaaaaaa.0p-24 + w * 0x91e9ee.0p-25);
    const R: f32 = t2 + t1;
    const r: f32 = s * (hfsq + R);

    var hi: f32 = f - hfsq;
    hx = @bitCast(hi);
    hi = @bitCast(@as(u32, @bitCast(hx)) & 0xfffff000);
    const lo: f32 = (f - hi) - hfsq + r;
    return (lo + hi) * -1.7605285393e-4 +
        lo * 1.4428710938e+0 +
        hi * 1.4428710938e+0 + y;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_log2.c
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
fn log2_64(x: f64) f64 {
    var hx: i32 = @bitCast(dbl64.getHighPart(x));
    const lx: i32 = @bitCast(dbl64.getLowPart(x));

    var k: i32 = 0;
    var xx: f64 = x;
    if (hx < 0x00100000) { // x < 2**-1022
        if (((hx & 0x7fffffff) | lx) == 0)
            return -std.math.inf(f64); // log(±0) = -inf

        if (hx < 0)
            return (x - x) / 0.0; // log(-#) = NaN

        k -= 54;
        xx *= 1.80143985094819840000e+16; // Subnormal number, scale up x
        hx = @bitCast(dbl64.getHighPart(xx));
    }

    if (hx >= 0x7ff00000)
        return xx + xx;

    if (hx == 0x3ff00000 and lx == 0)
        return 0.0; // log(1) = +0

    k +%= (hx >> 20) -% 1023;
    hx &= 0x000fffff;
    const i: i32 = (hx +% 0x95f64) & 0x100000;
    dbl64.setHighPart(&xx, @bitCast(hx | (i ^ 0x3ff00000))); // normalize x or x/2
    k +%= (i >> 20);
    const y: f64 = types.scast(f64, k);
    const f: f64 = xx - 1.0;
    const hfsq: f64 = 0.5 * f * f;
    const s: f64 = f / (2.0 + f);
    const z: f64 = s * s;
    var w: f64 = z * z;
    const t1: f64 = w *
        (3.999999999940941908e-1 + w *
            (2.222219843214978396e-1 + w *
                1.531383769920937332e-1));
    const t2: f64 = z *
        (6.666666666666735130e-1 + w *
            (2.857142874366239149e-1 + w *
                (1.818357216161805012e-1 + w *
                    1.479819860511658591e-1)));
    const R: f64 = t2 + t1;
    const r: f64 = s * (hfsq + R);

    var hi: f64 = f - hfsq;
    dbl64.setLowPart(&hi, 0);
    const lo: f64 = (f - hi) - hfsq + r;
    var val_hi: f64 = hi * 1.44269504072144627571e+0;
    var val_lo: f64 = (lo + hi) * 1.67517131648865118353e-10 + lo * 1.44269504072144627571e+0;

    w = y + val_hi;
    val_lo += (y - w) + val_hi;
    val_hi = w;

    return val_lo + val_hi;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_log2l.c
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
fn log2_128(x: f128) f128 {
    const hx: i64 = @bitCast(ldbl128.getHighPart(x));
    const lx: i64 = @bitCast(ldbl128.getLowPart(x));

    if (((hx & 0x7fffffffffffffff) | lx) == 0)
        return (-1.0 / (x - x));

    if (hx < 0)
        return (x - x) / (x - x);

    if (hx >= 0x7fff000000000000)
        return (x + x);

    // Note, frexp is used so that denormal numbers
    // will be handled properly
    var e: i32 = undefined;
    var xx: f128 = float.frexp(x, &e);

    // logarithm using log(xx) = z + z**3 P(z)/Q(z),
    // where z = 2 * (xx - 1)/(xx + 1)
    var z: f128 = undefined;
    var y: f128 = undefined;
    var skip: bool = false;
    if (e > 2 or e < -2) {
        if (xx < 7.071067811865475244008443621048490392848359e-1) { // 2 * (2 * xx - 1)/(2 * xx + 1)
            e -= 1;
            z = xx - 0.5;
            y = 0.5 * z + 0.5;
        } else { //  2 * (xx - 1)/(xx + 1)
            z = xx - 0.5;
            z -= 0.5;
            y = 0.5 * xx + 0.5;
        }
        xx = z / y;
        z = xx * xx;
        const num: f128 = 1.418134209872192732479751274970992665513e5 + z * (-8.977257995689735303686582344659576526998e4 + z *
            (2.048819892795278657810231591630928516206e4 + z *
                (-2.024301798136027039250415126250455056397e3 + z *
                    (8.057002716646055371965756206836056074715e1 + z *
                        (-8.828896441624934385266096344596648080902e-1)))));
        const den: f128 = 1.701761051846631278975701529965589676574e6 + z *
            (-1.332535117259762928288745111081235577029e6 + z *
                (4.001557694070773974936904547424676279307e5 + z *
                    (-5.748542087379434595104154610899551484314e4 + z *
                        (3.998526750980007367835804959888064681098e3 + z *
                            (-1.186359407982897997337150403816839480438e2 + z)))));
        y = xx * (z * num / den);
        skip = true;
    }

    // Logarithm using log(1 + xx) = xx - 0.5 * xx**2 + xx**3 * P(xx)/Q(xx)
    if (!skip) {
        if (xx < 7.071067811865475244008443621048490392848359e-1) {
            e -= 1;
            xx = 2.0 * xx - 1.0; //  2 * xx - 1
        } else {
            xx = xx - 1.0;
        }

        z = xx * xx;
        const num: f128 = 1.313572404063446165910279910527789794488e4 + xx *
            (7.771154681358524243729929227226708890930e4 + xx *
                (2.014652742082537582487669938141683759923e5 + xx *
                    (3.007007295140399532324943111654767187848e5 + xx *
                        (2.854829159639697837788887080758954924001e5 + xx *
                            (1.797628303815655343403735250238293741397e5 + xx *
                                (7.594356839258970405033155585486712125861e4 + xx *
                                    (2.128857716871515081352991964243375186031e4 + xx *
                                        (3.824952356185897735160588078446136783779e3 + xx *
                                            (4.114517881637811823002128927449878962058e2 + xx *
                                                (2.321125933898420063925789532045674660756e1 + xx *
                                                    (4.998469661968096229986658302195402690910e-1 + xx *
                                                        1.538612243596254322971797716843006400388e-6)))))))))));
        const den: f128 = 3.940717212190338497730839731583397586124e4 + xx *
            (2.626900195321832660448791748036714883242e5 + xx *
                (7.777690340007566932935753241556479363645e5 + xx *
                    (1.347518538384329112529391120390701166528e6 + xx *
                        (1.514882452993549494932585972882995548426e6 + xx *
                            (1.158019977462989115839826904108208787040e6 + xx *
                                (6.132189329546557743179177159925690841200e5 + xx *
                                    (2.248234257620569139969141618556349415120e5 + xx *
                                        (5.605842085972455027590989944010492125825e4 + xx *
                                            (9.147150349299596453976674231612674085381e3 + xx *
                                                (9.104928120962988414618126155557301584078e2 + xx *
                                                    (4.839208193348159620282142911143429644326e1 + xx)))))))))));
        y = xx * (z * num / den);
        y = y - 0.5 * z;
    }

    // Multiply log of fraction by log10(e)
    // and base 2 exxponent by log10(2)
    z = y * 4.4269504088896340735992468100189213742664595e-1;
    z += xx * 4.4269504088896340735992468100189213742664595e-1;
    z += y;
    z += xx;
    z += types.scast(f128, e);
    return z;
}
