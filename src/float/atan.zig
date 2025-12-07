const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn atan(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.atan: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, atan32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_atanf.c
            return atan32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_atan.c
            return atan64(types.scast(f64, x));
        },
        f80 => {
            //
            // return atan80(types.scast(f80, x));
            return types.scast(f80, atan128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_atanl.c
            return atan128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_atanf.c
//
// Original copyright notice:
// s_atanf.c -- float version of s_atan.c.
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
fn atan32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;
    if (ix >= 0x4c800000) { // If |x| >= 2**26
        if (ix > 0x7f800000)
            return x + x; // NaN

        if (hx > 0)
            return 1.5707962513e+0 + 7.5497894159e-8
        else
            return -1.5707962513e+0 - 7.5497894159e-8;
    }

    var id: i32 = undefined;
    var xx: f32 = x;
    if (ix < 0x3ee00000) { // |x| < 0.4375
        if (ix < 0x39800000) { // |x| < 2**-12
            if (1.0e30 + x > 1.0)
                return x; // Raise inexact
        }

        id = -1;
    } else {
        xx = float.abs(xx);
        if (ix < 0x3f980000) { // |x| < 1.1875
            if (ix < 0x3f300000) { // 7/16 <= |x| < 11/16
                id = 0;
                xx = (2.0 * xx - 1.0) / (2.0 + xx);
            } else { // 11/16 <= |x| < 19/16
                id = 1;
                xx = (xx - 1.0) / (xx + 1.0);
            }
        } else {
            if (ix < 0x401c0000) { // |x| < 2.4375
                id = 2;
                xx = (xx - 1.5) / (1.0 + 1.5 * xx);
            } else { // 2.4375 <= |x| < 2**26
                id = 3;
                xx = -1.0 / xx;
            }
        }
    }

    // End of argument reduction
    var z: f32 = xx * xx;
    const w: f32 = z * z;

    // Break sum from i = 0 to 10, aT[i] * z**(i+1) into odd and even poly
    const s1: f32 = z *
        (3.3333328366e-1 + w *
            (1.4253635705e-1 + w *
                6.1687607318e-2));
    const s2: f32 = w *
        (-1.9999158382e-1 + w *
            -1.0648017377e-1);

    if (id < 0)
        return xx - xx * (s1 + s2)
    else {
        const atanhi: f32 = switch (id) {
            0 => 4.6364760399e-1,
            1 => 7.8539812565e-1,
            2 => 9.8279368877e-1,
            3 => 1.5707962513e+0,
            else => undefined,
        };
        const atanlo: f32 = switch (id) {
            0 => 5.0121582440e-8,
            1 => 3.7748947079e-8,
            2 => 3.4473217170e-8,
            3 => 7.5497894159e-8,
            else => undefined,
        };

        z = atanhi - ((xx * (s1 + s2) - atanlo) - xx);
        return if (hx < 0) -z else z;
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_atan.c
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
fn atan64(x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = hx & 0x7fffffff;
    if (ix >= 0x44100000) { // |x| >= 2^66
        const low: u32 = dbl64.getLowPart(x);
        if (ix > 0x7ff00000 or (ix == 0x7ff00000 and low != 0))
            return x + x; // NaN
        if (hx > 0)
            return 1.57079632679489655800e+0 + 6.12323399573676603587e-17
        else
            return -1.57079632679489655800e+0 - 6.12323399573676603587e-17;
    }

    var id: i32 = undefined;
    var xx: f64 = x;
    if (ix < 0x3fdc0000) { // |x| < 0.4375
        if (ix < 0x3e400000) { // |x| < 2^-27
            if (1.0e300 + x > 1.0)
                return x; // Raise inexact
        }

        id = -1;
    } else {
        xx = float.abs(xx);
        if (ix < 0x3ff30000) { // |x| < 1.1875
            if (ix < 0x3fe60000) { // 7/16 <=|x|<11/16
                id = 0;
                xx = (2.0 * xx - 1.0) / (2.0 + xx);
            } else { // 11/16 <= |x| < 19/16
                id = 1;
                xx = (xx - 1.0) / (xx + 1.0);
            }
        } else {
            if (ix < 0x40038000) { // |x| < 2.4375
                id = 2;
                xx = (xx - 1.5) / (1.0 + 1.5 * xx);
            } else { // 2.4375 <= |x| < 2^66
                id = 3;
                xx = -1.0 / xx;
            }
        }
    }

    // End of argument reduction
    var z: f64 = xx * xx;
    const w: f64 = z * z;

    // Break sum from i = 0 to 10, aT[i] * z**(i + 1) into odd and even poly
    const s1: f64 = z *
        (3.33333333333329318027e-1 + w *
            (1.42857142725034663711e-1 + w *
                (9.09088713343650656196e-2 + w *
                    (6.66107313738753120669e-2 + w *
                        (4.97687799461593236017e-2 + w *
                            1.62858201153657823623e-2)))));
    const s2: f64 = w *
        (-1.99999999998764832476e-1 + w *
            (-1.11111104054623557880e-1 + w *
                (-7.69187620504482999495e-2 + w *
                    (-5.83357013379057348645e-2 + w *
                        -3.65315727442169155270e-2))));

    if (id < 0)
        return xx - xx * (s1 + s2)
    else {
        const atanhi: f64 = switch (id) {
            0 => 4.63647609000806093515e-1,
            1 => 7.85398163397448309616e-1,
            2 => 9.82793723247329054082e-1,
            3 => 1.57079632679489655800e+0,
            else => undefined,
        };
        const atanlo: f64 = switch (id) {
            0 => 2.26987774529616870924e-17,
            1 => 3.06161699786838301793e-17,
            2 => 1.39033110312309984516e-17,
            3 => 6.12323399573676603587e-17,
            else => undefined,
        };

        z = atanhi - ((xx * (s1 + s2) - atanlo) - xx);
        return if (hx < 0) -z else z;
    }
}

fn atan80(x: f80) f80 {
    _ = x;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_atanl.c
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
fn atan128(x: f128) f128 {
    const u: ldbl128.ShapeSplit = .fromFloat(x);
    const sign: u1 = u.sign;
    const expt: i32 = @intCast(u.exponent);

    if (expt >= (16384 + 112)) { // |x| is large, atan(x) ~= pi/2
        if (expt == (16384 - 1) + 16384 and (u.mantissa_high | u.mantissa_low) != 0) // NaN
            return x + x;

        if (sign == 0)
            return 1.57079632679489661923132169163975140e+0 + 4.33590506506189051239852201302167613e-35
        else
            return -1.57079632679489661923132169163975140e+0 - 4.33590506506189051239852201302167613e-35;
    }

    const expman: i32 = (expt << 8) | @as(i32, @bitCast(@as(u32, @truncate((u.mantissa_high >> 39) & 0xff))));
    var id: i32 = undefined;
    var xx: f128 = x;
    if (expman < ((16384 - 3) << 8) + 0xc0) { // |x| < 0.4375
        if (expt < (16384 - 57)) { // |x| is small, atanl(x) ~= x
            if (1.0e300 + x > 1.0)
                return x; // Raise inexact
        }
        id = -1;
    } else {
        xx = float.abs(xx);
        if (expman < ((16384 - 1) << 8) + 0x30) { // |x| < 1.1875
            if (expman < ((16384 - 2) << 8) + 0x60) { // 7/16 <= |x| < 11/16
                id = 0;
                xx = (2.0 * xx - 1.0) / (2.0 + xx);
            } else { // 11/16 <= |x| < 19/16
                id = 1;
                xx = (xx - 1.0) / (xx + 1.0);
            }
        } else {
            if (expman < (16384 << 8) + 0x38) { // |x| < 2.4375
                id = 2;
                xx = (xx - 1.5) / (1.0 + 1.5 * xx);
            } else { // 2.4375 <= |x| < 2**113
                id = 3;
                xx = -1.0 / xx;
            }
        }
    }

    // End of argument reduction
    var z: f128 = xx * xx;
    const w: f128 = z * z;
    // Break sum aT[i] * z**(i+1) into odd and even poly
    const s1: f128 = z *
        (3.33333333333333333333333333333333125e-1 + w *
            (1.42857142857142857142857142125269827e-1 + w *
                (9.09090909090909090908522355708623681e-2 + w *
                    (6.66666666666666660390096773046256096e-2 + w *
                        (5.26315789473666478515847092020327506e-2 + w *
                            (4.34782608678695085948531993458097026e-2 + w *
                                (3.70370363987423702891250829918659723e-2 + w *
                                    (3.22579620681420149871973710852268528e-2 + w *
                                        (2.85641979882534783223403715930946138e-2 + w *
                                            (2.54194698498808542954187110873675769e-2 + w *
                                                (2.04832358998165364349957325067131428e-2 + w *
                                                    8.64492360989278761493037861575248038e-3)))))))))));
    const s2: f128 = w *
        (-1.99999999999999999999999999999180430e-1 + w *
            (-1.11111111111111111111110834490810169e-1 + w *
                (-7.69230769230769230696553844935357021e-2 + w *
                    (-5.88235294117646671706582985209643694e-2 + w *
                        (-4.76190476189855517021024424991436144e-2 + w *
                            (-3.99999999632663469330634215991142368e-2 + w *
                                (-3.44827496515048090726669907612335954e-2 + w *
                                    (-3.03020767654269261041647570626778067e-2 + w *
                                        (-2.69824879726738568189929461383741323e-2 + w *
                                            (-2.35083879708189059926183138130183215e-2 + w *
                                                (-1.54489555488544397858507248612362957e-2 + w *
                                                    -2.58521121597609872727919154569765469e-3)))))))))));

    if (id < 0)
        return xx - xx * (s1 + s2)
    else {
        const atanhi: f128 = switch (id) {
            0 => 4.63647609000806116214256231461214397e-1,
            1 => 7.85398163397448309615660845819875699e-1,
            2 => 9.82793723247329067985710611014666038e-1,
            3 => 1.57079632679489661923132169163975140e+0,
            else => undefined,
        };
        const atanlo: f128 = switch (id) {
            0 => 4.89509642257333492668618435220297706e-36,
            1 => 2.16795253253094525619926100651083806e-35,
            2 => -2.31288434538183565909319952098066272e-35,
            3 => 4.33590506506189051239852201302167613e-35,
            else => undefined,
        };

        z = atanhi - ((xx * (s1 + s2) - atanlo) - xx);
        return if (sign == 1) -z else z;
    }
}
