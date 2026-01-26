const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Log1p(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.log1p: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

/// Returns the natural logarithm $\log(x + 1)$ of a float, int or bool operand
/// plus one. The result type is determined by coercing the operand type to a
/// float, and the operation is performed by casting the operand to the result
/// type, then computing its logarithm plus one.
///
/// ## Signature
/// ```zig
/// float.log1p(x: X) float.Log1p(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the logarithm plus one of.
///
/// ## Returns
/// `float.Log1p(@TypeOf(x))`: The logarithm of `x + 1`.
pub inline fn log1p(x: anytype) float.Log1p(@TypeOf(x)) {
    switch (float.Log1p(@TypeOf(x))) {
        f16 => return types.scast(f16, log1p32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_log1pf.c
            return log1p32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_log1p.c
            return log1p64(types.scast(f64, x));
        },
        f80 => {
            //
            // return log1p80(types.scast(f80, x));
            return types.scast(f80, log1p128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/s_log1pl.c
            return log1p128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_log1pf.c
//
// Original copyright notice:
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
fn log1p32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ax: i32 = hx & 0x7fffffff;

    var k: i32 = 1;
    var f: f32 = undefined;
    var hu: i32 = undefined;
    if (hx < 0x3ed413d0) { // 1 + x < sqrt(2)+
        if (ax >= 0x3f800000) { // x <= -1.0
            if (x == -1.0)
                return std.math.inf(f32) // log1p(-1) = +inf
            else
                return (x - x) / (x - x); // log1p(x < -1) = NaN
        }

        if (ax < 0x38000000) { // |x| < 2**-15
            if (ax < 0x33800000) // |x| < 2**-24
                return x
            else
                return x - x * x * 0.5;
        }

        if (hx > 0 or hx <= @as(i32, @bitCast(@as(u32, (0xbe95f619))))) { // sqrt(2)/2- <= 1 + x < sqrt(2)+
            k = 0;
            f = x;
            hu = 1;
        }
    }

    if (hx >= 0x7f800000)
        return x + x;

    var c: f32 = undefined;
    if (k != 0) {
        var u: f32 = undefined;
        if (hx < 0x5a000000) {
            u = 1.0 + x;
            hu = @bitCast(u);
            k = (hu >> 23) - 127;

            // Correction term
            c = if (k > 0) 1.0 - (u - x) else x - (u - 1.0);
            c /= u;
        } else {
            u = x;
            hu = @bitCast(u);
            k = (hu >> 23) -% 127;
            c = 0;
        }

        hu &= 0x007fffff;

        // The approximation to sqrt(2) used in thresholds is not
        // critical.  However, the ones used above must give less
        // strict bounds than the one here so that the k==0 case is
        // never reached from here, since here we have committed to
        // using the correction term but don't use it if k==0.
        if (hu < 0x3504f4) { // u < sqrt(2)
            u = @bitCast(hu | 0x3f800000); // Normalize u
        } else {
            k += 1;
            u = @bitCast(hu | 0x3f000000); // Normalize u/2
            hu = (0x00800000 - hu) >> 2;
        }

        f = u - 1.0;
    }

    const hfsq: f32 = 0.5 * f * f;
    if (hu == 0) { // |f| < 2**-20
        if (f == 0.0) {
            if (k == 0) {
                return 0.0;
            } else {
                c += types.scast(f32, k) * 9.0580006145e-6;
                return types.scast(f32, k) * 6.9313812256e-1 + c;
            }
        }

        const R: f32 = hfsq * (1.0 - 0.66666666666666666 * f);
        if (k == 0)
            return f - R
        else
            return types.scast(f32, k) * 6.9313812256e-1 - ((R - (types.scast(f32, k) * 9.0580006145e-6 + c)) - f);
    }

    const s: f32 = f / (2.0 + f);
    const z: f32 = s * s;
    const R: f32 = z *
        (6.6666668653e-1 + z *
            (4.0000000596e-1 + z *
                (2.8571429849e-1 + z *
                    (2.2222198546e-1 + z *
                        (1.8183572590e-1 + z *
                            (1.5313838422e-1 + z *
                                1.4798198640e-1))))));
    if (k == 0)
        return f - (hfsq - s * (hfsq + R))
    else
        return types.scast(f32, k) * 6.9313812256e-1 - ((hfsq - (s * (hfsq + R) + (types.scast(f32, k) * 9.0580006145e-6 + c))) - f);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_log1p.c
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
fn log1p64(x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ax: i32 = hx & 0x7fffffff;

    var k: i32 = 1;
    var f: f64 = undefined;
    var hu: i32 = undefined;
    if (hx < 0x3fda827a) { // 1 + x < sqrt(2)+
        if (ax >= 0x3ff00000) { // x <= -1.0
            if (x == -1.0)
                return std.math.inf(f64) // log1p(-1) = +inf
            else
                return (x - x) / (x - x); // log1p(x < -1) = NaN
        }

        if (ax < 0x3e200000) { // |x| < 2**-29
            if (ax < 0x3c900000) // |x| < 2**-54
                return x
            else
                return x - x * x * 0.5;
        }

        if (hx > 0 or hx <= @as(i32, @bitCast(@as(u32, (0xbfd2bec4))))) { // sqrt(2)/2- <= 1 + x < sqrt(2)+
            k = 0;
            f = x;
            hu = 1;
        }
    }

    if (hx >= 0x7ff00000)
        return x + x;

    var c: f64 = undefined;
    if (k != 0) {
        var u: f64 = undefined;
        if (hx < 0x43400000) {
            u = 1.0 + x;
            hu = @bitCast(dbl64.getHighPart(u));
            k = (hu >> 20) - 1023;

            // Correction term
            c = if (k > 0) 1.0 - (u - x) else x - (u - 1.0);
            c /= u;
        } else {
            u = x;
            hu = @bitCast(dbl64.getHighPart(u));
            k = (hu >> 20) - 1023;
            c = 0;
        }

        hu &= 0x000fffff;

        // The approximation to sqrt(2) used in thresholds is not
        // critical.  However, the ones used above must give less
        // strict bounds than the one here so that the k==0 case is
        // never reached from here, since here we have committed to
        // using the correction term but don't use it if k==0.
        if (hu < 0x6a09e) { // u ~< sqrt(2)
            dbl64.setHighPart(&u, @bitCast(hu | 0x3ff00000)); // Normalize u
        } else {
            k += 1;
            dbl64.setHighPart(&u, @bitCast(hu | 0x3fe00000)); // Normalize u/2
            hu = (0x00100000 - hu) >> 2;
        }

        f = u - 1.0;
    }

    const hfsq: f64 = 0.5 * f * f;
    if (hu == 0) { // |f| < 2**-20
        if (f == 0.0) {
            if (k == 0) {
                return 0.0;
            } else {
                c += types.scast(f64, k) * 1.90821492927058770002e-10;
                return types.scast(f64, k) * 6.93147180369123816490e-1 + c;
            }
        }

        const R: f64 = hfsq * (1.0 - 0.66666666666666666 * f);
        if (k == 0)
            return f - R
        else
            return types.scast(f64, k) * 6.93147180369123816490e-1 - ((R - (types.scast(f64, k) * 1.90821492927058770002e-10 + c)) - f);
    }

    const s: f64 = f / (2.0 + f);
    const z: f64 = s * s;
    const R: f64 = z *
        (6.666666666666735130e-1 + z *
            (3.999999999940941908e-1 + z *
                (2.857142874366239149e-1 + z *
                    (2.222219843214978396e-1 + z *
                        (1.818357216161805012e-1 + z *
                            (1.531383769920937332e-1 + z *
                                1.479819860511658591e-1))))));
    if (k == 0)
        return f - (hfsq - s * (hfsq + R))
    else
        return types.scast(f64, k) * 6.93147180369123816490e-1 - ((hfsq - (s * (hfsq + R) + (types.scast(f64, k) * 1.90821492927058770002e-10 + c))) - f);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/s_log1pl.c
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
fn log1p128(x: f128) f128 {
    const u: ldbl128.Parts32 = .fromFloat(x);
    const hx: i32 = @bitCast(u.mswhi);

    if (hx >= 0x7fff0000) // x is ±inf or NaN
        return x;

    if ((hx & 0x7fffffff) == 0 and (u.mswlo | u.lswhi | u.lswlo) == 0) // log1p(±0) = ±0
        return x;

    var xx: f128 = x + 1.0;

    // log1p(-1) = -inf
    if (xx <= 0.0) {
        if (xx == 0.0)
            return (-1.0 / (xx - xx))
        else
            return (0.0 / (xx - xx));
    }

    // Separate mantissa from exponent.
    // Use frexp used so that denormal numbers will be handled properly.
    var e: i32 = undefined;
    xx = float.frexp(xx, &e);

    // Logarithm using log(x) = z + z**3 P(z**2)/Q(z**2),
    // where z = 2 * (x - 1)/(x + 1)
    if (e > 2 or e < -2) {
        var z: f128 = undefined;
        var y: f128 = undefined;
        if (xx < 0.7071067811865475244008443621048490392848) { // 2 * (2 * x - 1)/(2 * x + 1)
            e -%= 1;
            z = xx - 0.5;
            y = 0.5 * z + 0.5;
        } else { // 2 * (xx - 1)/(xx + 1)
            z = xx - 0.5;
            z -= 0.5;
            y = 0.5 * xx + 0.5;
        }

        xx = z / y;
        z = xx * xx;
        const r: f128 = ((((-8.828896441624934385266096344596648080902e-1 * z +
            8.057002716646055371965756206836056074715e1) * z +
            -2.024301798136027039250415126250455056397e3) * z +
            2.048819892795278657810231591630928516206e4) * z +
            -8.977257995689735303686582344659576526998e4) * z +
            1.418134209872192732479751274970992665513e5;
        const s: f128 = (((((z +
            -1.186359407982897997337150403816839480438e2) * z +
            3.998526750980007367835804959888064681098e3) * z +
            -5.748542087379434595104154610899551484314e4) * z +
            4.001557694070773974936904547424676279307e5) * z +
            -1.332535117259762928288745111081235577029e6) * z +
            1.701761051846631278975701529965589676574e6;
        z = xx * (z * r / s);
        z = z + types.scast(f128, e) * 1.428606820309417232121458176568075500134e-6;
        z = z + xx;
        z = z + types.scast(f128, e) * 6.93145751953125e-1;
        return z;
    }

    // Logarithm using log(1 + x) = x - 0.5 * x**2 + x**3 P(x)/Q(x).
    if (xx < 0.7071067811865475244008443621048490392848) {
        e -%= 1;
        if (e != 0)
            xx = 2.0 * xx - 1.0 // 2 * xx - 1
        else
            xx = x;
    } else {
        if (e != 0)
            xx -= 1.0
        else
            xx = x;
    }

    var z: f128 = xx * xx;
    const r: f128 = (((((((((((1.538612243596254322971797716843006400388e-6 * xx +
        4.998469661968096229986658302195402690910e-1) * xx +
        2.321125933898420063925789532045674660756e1) * xx +
        4.114517881637811823002128927449878962058e2) * xx +
        3.824952356185897735160588078446136783779e3) * xx +
        2.128857716871515081352991964243375186031e4) * xx +
        7.594356839258970405033155585486712125861e4) * xx +
        1.797628303815655343403735250238293741397e5) * xx +
        2.854829159639697837788887080758954924001e5) * xx +
        3.007007295140399532324943111654767187848e5) * xx +
        2.014652742082537582487669938141683759923e5) * xx +
        7.771154681358524243729929227226708890930e4) * xx +
        1.313572404063446165910279910527789794488e4;
    const s: f128 = (((((((((((xx + 4.839208193348159620282142911143429644326e1) * xx +
        9.104928120962988414618126155557301584078e2) * xx +
        9.147150349299596453976674231612674085381e3) * xx +
        5.605842085972455027590989944010492125825e4) * xx +
        2.248234257620569139969141618556349415120e5) * xx +
        6.132189329546557743179177159925690841200e5) * xx +
        1.158019977462989115839826904108208787040e6) * xx +
        1.514882452993549494932585972882995548426e6) * xx +
        1.347518538384329112529391120390701166528e6) * xx +
        7.777690340007566932935753241556479363645e5) * xx +
        2.626900195321832660448791748036714883242e5) * xx +
        3.940717212190338497730839731583397586124e4;
    var y: f128 = xx * (z * r / s);
    y = y + types.scast(f128, e) * 1.428606820309417232121458176568075500134e-6;
    z = y - 0.5 * z;
    z = z + xx;
    z = z + types.scast(f128, e) * 6.93145751953125e-1;
    return z;
}
