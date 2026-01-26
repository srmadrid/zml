const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const erf_data = @import("erf_data.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Erf(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.erf: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

/// Returns the error function of a float, int or bool operand. The result type
/// is determined by coercing the operand type to a float, and the operation is
/// performed by casting the operand to the result type, then computing the
/// error function at that value.
///
/// The error function is defined as:
/// $$
/// \mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// float.erf(x: X) float.Erf(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the error function at.
///
/// ## Returns
/// `float.Erf(@TypeOf(x))`: The error function at `x`.
pub inline fn erf(x: anytype) float.Erf(@TypeOf(x)) {
    switch (float.Erf(@TypeOf(x))) {
        f16 => return types.scast(f16, erf32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_erff.c
            return erf32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_erf.c
            return erf64(types.scast(f64, x));
        },
        f80 => {
            //
            // return erf80(types.scast(f80, x));
            return types.scast(f80, erf128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/s_erfl.c
            return erf128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_erff.c
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
fn erf32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x7f800000) { // erf(nan) = nan
        const i: i32 = @bitCast((@as(u32, @bitCast(hx)) >> 31) << 1);
        return 1 - types.scast(f32, i); // erf(±inf) = ±1
    }

    if (ix < 0x3f580000) { // |x| < 0.84375
        if (ix < 0x38800000) { // |x| < 2**-14
            if (ix < 0x04000000) // |x| < 0x1p-119
                return (8 * x + 1.0270333290e+0 * x) / 8; // Avoid spurious underflow

            return x + 1.2837916613e-1 * x;
        }

        const z: f32 = x * x;
        const r: f32 = 1.28379166e-1 + z *
            (-3.36030394e-1 + z *
                -1.86260219e-3);
        const s: f32 = 1.0 + z *
            (3.12324286e-1 + z *
                (2.16070302e-2 + z *
                    -1.98859419e-3));
        const y: f32 = r / s;
        return x + x * y;
    }

    if (ix < 0x3fa00000) { // 0.84375 <= |x| < 1.25
        const s: f32 = float.abs(x) - 1.0;
        const P: f32 = 3.64939137e-6 + s *
            (4.15109694e-1 + s *
                (-1.65179938e-1 + s *
                    1.10914491e-1));
        const Q: f32 = 1.0 + s *
            (6.02074385e-1 + s *
                (5.35934687e-1 + s *
                    (1.68576106e-1 + s *
                        5.62181212e-2)));
        if (hx >= 0)
            return 8.42697144e-1 + P / Q
        else
            return -8.42697144e-1 - P / Q;
    }

    if (ix >= 0x40800000) { // inf > |x| >= 4
        if (hx >= 0)
            return 1.0
        else
            return -1.0;
    }

    const xx: f32 = float.abs(x);
    const s: f32 = 1.0 / (xx * xx);
    var R: f32 = undefined;
    var S: f32 = undefined;
    if (ix < 0x4036db6e) { // |x| < 1/0.35
        R = -9.87132732e-3 + s *
            (-5.53605914e-1 + s *
                (-2.17589188e+0 + s *
                    -1.43268085e+0));
        S = 1.0 + s *
            (5.45995426e+0 + s *
                (6.69798088e+0 + s *
                    (1.43113089e+0 + s *
                        -5.77397496e-2)));
    } else { // |x| >= 1/0.35
        R = -9.86494310e-3 + s *
            (-6.25171244e-1 + s *
                (-6.16498327e+0 + s *
                    (-1.66696873e+1 + s *
                        -9.53764343e+0)));
        S = 1.0 + s *
            (1.26884899e+1 + s *
                (4.51839523e+1 + s *
                    (4.72810211e+1 + s *
                        8.93033314e+0)));
    }

    const z: f32 = @bitCast(@as(u32, @bitCast(hx)) & 0xffffe000);
    const r: f32 = float.exp(-z * z - 0.5625) * float.exp((z - xx) * (z + xx) + R / S);
    if (hx >= 0)
        return 1.0 - r / xx
    else
        return r / xx - 1.0;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_erf.c
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
fn erf64(x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x7ff00000) { // erf(nan) = nan
        const i: i32 = @bitCast((@as(u32, @bitCast(hx)) >> 31) << 1);
        return 1 - types.scast(f64, i); // erf(±inf) = ±1
    }

    if (ix < 0x3feb0000) { // |x| < 0.84375
        if (ix < 0x3e300000) { // |x| < 2**-28
            if (ix < 0x00800000)
                return 0.125 * (8.0 * x + 1.02703333676410069053e+0 * x); // Avoid underflow

            return x + 1.28379167095512586316e-1 * x;
        }

        const z: f64 = x * x;
        const r: f64 = 1.28379167095512558561e-1 + z *
            (-3.25042107247001499370e-1 + z *
                (-2.84817495755985104766e-2 + z *
                    (-5.77027029648944159157e-3 + z *
                        -2.37630166566501626084e-5)));
        const s: f64 = 1.0 + z *
            (3.97917223959155352819e-1 + z *
                (6.50222499887672944485e-2 + z *
                    (5.08130628187576562776e-3 + z *
                        (1.32494738004321644526e-4 + z *
                            -3.96022827877536812320e-6))));
        const y: f64 = r / s;
        return x + x * y;
    }

    if (ix < 0x3ff40000) { // 0.84375 <= |x| < 1.25
        const s: f64 = float.abs(x) - 1.0;
        const P: f64 = -2.36211856075265944077e-3 + s *
            (4.14856118683748331666e-1 + s *
                (-3.72207876035701323847e-1 + s *
                    (3.18346619901161753674e-1 + s *
                        (-1.10894694282396677476e-1 + s *
                            (3.54783043256182359371e-2 + s *
                                -2.16637559486879084300e-3)))));
        const Q: f64 = 1.0 + s *
            (1.06420880400844228286e-1 + s *
                (5.40397917702171048937e-1 + s *
                    (7.18286544141962662868e-2 + s *
                        (1.26171219808761642112e-1 + s *
                            (1.36370839120290507362e-2 + s *
                                1.19844998467991074170e-2)))));
        if (hx >= 0)
            return 8.45062911510467529297e-1 + P / Q
        else
            return -8.45062911510467529297e-1 - P / Q;
    }

    if (ix >= 0x40180000) { // inf > |x| >= 6
        if (hx >= 0)
            return 1.0
        else
            return -1.0;
    }

    const xx: f64 = float.abs(x);
    const s: f64 = 1.0 / (xx * xx);
    var R: f64 = undefined;
    var S: f64 = undefined;
    if (ix < 0x4006db6e) { // |x| < 1/0.35
        R = -9.86494403484714822705e-3 + s *
            (-6.93858572707181764372e-1 + s *
                (-1.05586262253232909814e+1 + s *
                    (-6.23753324503260060396e+1 + s *
                        (-1.62396669462573470355e+2 + s *
                            (-1.84605092906711035994e+2 + s *
                                (-8.12874355063065934246e+1 + s *
                                    -9.81432934416914548592e+0))))));
        S = 1.0 + s *
            (1.96512716674392571292e+1 + s *
                (1.37657754143519042600e+2 + s *
                    (4.34565877475229228821e+2 + s *
                        (6.45387271733267880336e+2 + s *
                            (4.29008140027567833386e+2 + s *
                                (1.08635005541779435134e+2 + s *
                                    (6.57024977031928170135e+0 + s *
                                        -6.04244152148580987438e-2)))))));
    } else { // |x| >= 1/0.35
        R = -9.86494292470009928597e-3 + s *
            (-7.99283237680523006574e-1 + s *
                (-1.77579549177547519889e+1 + s *
                    (-1.60636384855821916062e+2 + s *
                        (-6.37566443368389627722e+2 + s *
                            (-1.02509513161107724954e+3 + s *
                                -4.83519191608651397019e+2)))));
        S = 1.0 + s *
            (3.03380607434824582924e+1 + s *
                (3.25792512996573918826e+2 + s *
                    (1.53672958608443695994e+3 + s *
                        (3.19985821950859553908e+3 + s *
                            (2.55305040643316442583e+3 + s *
                                (4.74528541206955367215e+2 + s *
                                    -2.24409524465858183362e+1))))));
    }

    var z: f64 = xx;
    dbl64.setLowPart(&z, 0);
    const r: f64 = float.exp(-z * z - 0.5625) * float.exp((z - xx) * (z + xx) + R / S);
    if (hx >= 0)
        return 1.0 - r / xx
    else
        return r / xx - 1.0;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/s_erfl.c
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
//
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
fn erf128(x: f128) f128 {
    var u: ldbl128.Parts32 = .fromFloat(x);
    const sign: i32 = @bitCast(u.mswhi);
    const ix: i32 = sign & 0x7fffffff;

    if (ix >= 0x7fff0000) { // erf(nan) = nan
        const i: i32 = @bitCast(((@as(u32, @bitCast(sign)) & 0xffff0000) >> 31) << 1);
        return types.scast(f128, (1 - i)); // erf(±inf) = ±1
    }

    if (ix >= 0x3fff0000) { // |x| >= 1.0
        const y: f128 = float.erfc(x);
        return 1 - y;
        // return (1 - __erfcl (x));
    }

    u.mswhi = @bitCast(ix);
    var a: f128 = u.toFloat();
    const z: f128 = x * x;
    var y: f128 = undefined;
    if (ix < 0x3ffec000) { // a < 0.875
        if (ix < 0x3fc60000) { // |x| < 2**-57
            if (ix < 0x00080000) // Avoid underflow
                return 0.125 * (8.0 * x + erf_data.efx8_128 * x);

            return x + erf_data.efx_128 * x;
        }

        y = a + a * erf_data.neval(z, &erf_data.TN1_128, 8) / erf_data.deval(z, &erf_data.TD1_128, 8);
    } else {
        a = a - 1.0;
        y = erf_data.erf_const_128 + erf_data.neval(a, &erf_data.TN2_128, 8) / erf_data.deval(a, &erf_data.TD2_128, 8);
    }

    if (sign < 0)
        return -y
    else
        return y;
}
