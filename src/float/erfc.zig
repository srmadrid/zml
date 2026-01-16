const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const erf_data = @import("erf_data.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Erfc(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.erfc: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

pub inline fn erfc(x: anytype) Erfc(@TypeOf(x)) {
    switch (Erfc(@TypeOf(x))) {
        f16 => return types.scast(f16, erfc32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_erff.c
            return erfc32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_erf.c
            return erfc64(types.scast(f64, x));
        },
        f80 => {
            //
            // return erfc80(types.scast(f80, x));
            return types.scast(f80, erfc128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/s_erfl.c
            return erfc128(types.scast(f128, x));
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
fn erfc32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x7f800000) { // erfc(nan) = nan
        // erfc(±inf) = 0,2
        return types.scast(f32, (@as(u32, @bitCast(hx)) >> 31) << 1);
    }

    if (ix < 0x3f580000) { // |x| < 0.84375
        if (ix < 0x33800000) // |x| < 2**-56
            return 1.0 - x;

        const z: f32 = x * x;
        var r: f32 = 1.28379166e-1 + z *
            (-3.36030394e-1 + z *
                -1.86260219e-3);
        const s: f32 = 1.0 + z *
            (3.12324286e-1 + z *
                (2.16070302e-2 + z *
                    -1.98859419e-3));
        const y: f32 = r / s;
        if (hx < 0x3e800000) { // x < 1/4
            return 1.0 - (x + x * y);
        } else {
            r = x * y;
            r += (x - 0.5);
            return 0.5 - r;
        }
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
        if (hx >= 0) {
            const z: f32 = 1.0 - 8.42697144e-1;
            return z - P / Q;
        } else {
            const z: f32 = 8.42697144e-1 + P / Q;
            return 1.0 + z;
        }
    }

    if (ix < 0x41300000) { // |x| < 28
        const xx: f32 = float.abs(x);
        const s: f32 = 1.0 / (xx * xx);
        var R: f32 = undefined;
        var S: f32 = undefined;
        if (ix < 0x4036DB6D) { // |x| < 1/.35 ~ 2.857143
            R = -9.87132732e-3 + s *
                (-5.53605914e-1 + s *
                    (-2.17589188e+0 + s *
                        -1.43268085e+0));
            S = 1.0 + s *
                (5.45995426e+0 + s *
                    (6.69798088e+0 + s *
                        (1.43113089e+0 + s *
                            -5.77397496e-2)));
        } else { // |x| >= 1/.35 ~ 2.857143
            if (hx < 0 and ix >= 0x40a00000)
                return 2.0; // x < -5

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
        if (hx > 0)
            return r / xx
        else
            return 2.0 - r / xx;
    } else {
        if (hx > 0)
            return 0.0
        else
            return 2.0;
    }
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
fn erfc64(x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x7ff00000) { // erfc(nan)=nan
        // erfc(±inf) = 0, 2
        return types.scast(f32, (@as(u32, @bitCast(hx)) >> 31) << 1);
    }

    if (ix < 0x3feb0000) { // |x| < 0.84375
        if (ix < 0x3c700000) // |x| < 2**-56
            return 1.0 - x;

        const z: f64 = x * x;
        var r: f64 = 1.28379167095512558561e-1 + z *
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
        if (hx < 0x3fd00000) { // x < 1/4
            return 1.0 - (x + x * y);
        } else {
            r = x * y;
            r += (x - 0.5);
            return 0.5 - r;
        }
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
        if (hx >= 0) {
            const z: f64 = 1.0 - 8.45062911510467529297e-1;
            return z - P / Q;
        } else {
            const z: f64 = 8.45062911510467529297e-1 + P / Q;
            return 1.0 + z;
        }
    }

    if (ix < 0x403c0000) { // |x|<28
        const xx: f64 = float.abs(x);
        const s: f64 = 1.0 / (xx * xx);
        var R: f64 = undefined;
        var S: f64 = undefined;
        if (ix < 0x4006db6d) { // |x| < 1/.35 ~ 2.857143
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
        } else { // |x| >= 1/.35 ~ 2.857143
            if (hx < 0 and ix >= 0x40180000)
                return 2.0; // x < -6

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
        if (hx > 0)
            return r / xx
        else
            return 2.0 - r / xx;
    } else {
        if (hx > 0)
            return 0.0
        else
            return 2.0;
    }
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
fn erfc128(x: f128) f128 {
    var u: ldbl128.Parts32 = .fromFloat(x);
    const sign: i32 = @bitCast(u.mswhi);
    const ix: i32 = sign & 0x7fffffff;
    u.mswhi = @bitCast(ix);

    if (ix >= 0x7fff0000) { // erfc(nan) = nan
        // erfc(±inf) = 0, 2
        return types.scast(f128, (@as(u32, @intCast(sign)) >> 31) << 1);
    }

    if (ix < 0x3ffd0000) { // |x| < 1/4
        if (ix < 0x3f8d0000) // |x| < 2**-114
            return 1 - x;

        return 1 - float.erf(x);
    }

    if (ix < 0x3fff4000) { // 1.25
        const xx: f128 = u.toFloat();
        const i: i32 = types.scast(i32, 8.0 * xx);
        var y: f128 = undefined;
        switch (i) {
            2 => {
                const z: f128 = xx - 0.25;
                y = erf_data.C13b_128 + z * erf_data.neval(z, &erf_data.RNr13_128, 8) / erf_data.deval(z, &erf_data.RDr13_128, 7);
                y += erf_data.C13a_128;
            },
            3 => {
                const z: f128 = xx - 0.375;
                y = erf_data.C14b_128 + z * erf_data.neval(z, &erf_data.RNr14_128, 8) / erf_data.deval(z, &erf_data.RDr14_128, 7);
                y += erf_data.C14a_128;
            },
            4 => {
                const z: f128 = xx - 0.5;
                y = erf_data.C15b_128 + z * erf_data.neval(z, &erf_data.RNr15_128, 8) / erf_data.deval(z, &erf_data.RDr15_128, 7);
                y += erf_data.C15a_128;
            },
            5 => {
                const z: f128 = xx - 0.625;
                y = erf_data.C16b_128 + z * erf_data.neval(z, &erf_data.RNr16_128, 8) / erf_data.deval(z, &erf_data.RDr16_128, 7);
                y += erf_data.C16a_128;
            },
            6 => {
                const z: f128 = xx - 0.75;
                y = erf_data.C17b_128 + z * erf_data.neval(z, &erf_data.RNr17_128, 8) / erf_data.deval(z, &erf_data.RDr17_128, 7);
                y += erf_data.C17a_128;
            },
            7 => {
                const z: f128 = xx - 0.875;
                y = erf_data.C18b_128 + z * erf_data.neval(z, &erf_data.RNr18_128, 8) / erf_data.deval(z, &erf_data.RDr18_128, 7);
                y += erf_data.C18a_128;
            },
            8 => {
                const z: f128 = xx - 1;
                y = erf_data.C19b_128 + z * erf_data.neval(z, &erf_data.RNr19_128, 8) / erf_data.deval(z, &erf_data.RDr19_128, 7);
                y += erf_data.C19a_128;
            },
            else => {
                const z: f128 = xx - 1.125;
                y = erf_data.C20b_128 + z * erf_data.neval(z, &erf_data.RNr20_128, 8) / erf_data.deval(z, &erf_data.RDr20_128, 7);
                y += erf_data.C20a_128;
            },
        }

        if (sign >= 0)
            return y
        else
            return 2.0 - y;
    }

    // 1.25 < |x| < 107
    if (ix < 0x4005ac00) {
        // x < -9
        if ((ix >= 0x40022000) and (@as(u32, @bitCast(sign)) & 0x80000000) != 0)
            return 2;

        const xx: f128 = float.abs(x);
        var z: f128 = 1 / (xx * xx);
        const i: i32 = types.scast(i32, 8.0 / xx);
        var p: f128 = undefined;
        switch (i) {
            1 => {
                p = erf_data.neval(z, &erf_data.RNr2_128, 11) / erf_data.deval(z, &erf_data.RDr2_128, 10);
            },
            2 => {
                p = erf_data.neval(z, &erf_data.RNr3_128, 11) / erf_data.deval(z, &erf_data.RDr3_128, 10);
            },
            3 => {
                p = erf_data.neval(z, &erf_data.RNr4_128, 10) / erf_data.deval(z, &erf_data.RDr4_128, 10);
            },
            4 => {
                p = erf_data.neval(z, &erf_data.RNr5_128, 10) / erf_data.deval(z, &erf_data.RDr5_128, 9);
            },
            5 => {
                p = erf_data.neval(z, &erf_data.RNr6_128, 9) / erf_data.deval(z, &erf_data.RDr6_128, 9);
            },
            6 => {
                p = erf_data.neval(z, &erf_data.RNr7_128, 9) / erf_data.deval(z, &erf_data.RDr7_128, 9);
            },
            7 => {
                p = erf_data.neval(z, &erf_data.RNr8_128, 9) / erf_data.deval(z, &erf_data.RDr8_128, 8);
            },
            else => {
                p = erf_data.neval(z, &erf_data.RNr1_128, 9) / erf_data.deval(z, &erf_data.RDr1_128, 8);
            },
        }

        u = .fromFloat(xx);
        u.lswlo = 0;
        u.lswhi &= 0xfe000000;
        z = u.toFloat();
        const r: f128 = float.exp(-z * z - 0.5625) * float.exp((z - xx) * (z + xx) + p);

        if (sign >= 0)
            return r / xx
        else
            return 2.0 - r / xx;
    } else {
        if (sign >= 0)
            return 0.0
        else
            return 2.0;
    }
}
