const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn asin(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.asin: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, asin32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_asinf.c
            return asin32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_asin.c
            return asin64(types.scast(f64, x));
        },
        f80 => {
            //
            // return asin80(types.scast(f80, x));
            return types.scast(f80, asin128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_asinl.c
            return asin128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_asinf.c
//
// Original copyright notice:
// e_asinf.c -- float version of e_asin.c.
// Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. buasiness.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn asin32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x3f800000) { // |x| >= 1
        if (ix == 0x3f800000) // |x| == 1
            return x * 1.570796326794896558e+0; // asin(±1) = ±pi/2 with inexact

        return (x - x) / (x - x); // asin(|x| > 1) is NaN
    } else if (ix < 0x3f000000) { // |x| < 0.5
        if (ix < 0x39800000) { // |x| < 2**-12
            if (1.000e+30 + x > 1.0)
                return x; // return x with inexact if x != 0
        }

        const t: f32 = x * x;
        const p: f32 = t * (1.6666586697e-1 + t * (-4.2743422091e-2 + t * -8.6563630030e-3));
        const q: f32 = 1.0 + t * -7.0662963390e-1;
        const w: f32 = p / q;
        return x + x * w;
    }

    // 1 > |x| >= 0.5
    var w: f32 = 1.0 - float.abs(x);
    var t: f32 = w * 0.5;
    const p: f32 = t * (1.6666586697e-1 + t * (-4.2743422091e-2 + t * -8.6563630030e-3));
    const q: f32 = 1.0 + t * -7.0662963390e-1;
    const s: f64 = float.sqrt(t);
    w = p / q;
    t = types.scast(f32, 1.570796326794896558e+0 - 2.0 * (s + s * types.scast(f64, w)));
    if (hx > 0)
        return t
    else
        return -t;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_asin.c
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
fn asin64(x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x3ff00000) { // |x| >= 1
        const lx: i32 = @bitCast(dbl64.getLowPart(x));
        if (((ix - 0x3ff00000) | lx) == 0) // asin(1) = ±pi/2 with inexact
            return x * 1.57079632679489655800e+0 + x * 6.12323399573676603587e-17;

        return (x - x) / (x - x); // asin(|x| > 1) is NaN
    } else if (ix < 0x3fe00000) { // |x| < 0.5
        if (ix < 0x3e500000) { // If |x| < 2**-26
            if (1.000e+300 + x > 1.0)
                return x; // Return x with inexact if x != 0
        }

        const t: f64 = x * x;
        const p: f64 = t *
            (1.66666666666666657415e-1 + t *
                (-3.25565818622400915405e-1 + t *
                    (2.01212532134862925881e-1 + t *
                        (-4.00555345006794114027e-2 + t *
                            (7.91534994289814532176e-4 + t *
                                3.47933107596021167570e-5)))));
        const q: f64 = 1.0 + t *
            (-2.40339491173441421878e+0 + t *
                (2.02094576023350569471e+0 + t *
                    (-6.88283971605453293030e-1 + t *
                        7.70381505559019352791e-2)));
        const w: f64 = p / q;
        return x + x * w;
    }

    // 1 > |x| >= 0.5
    var w: f64 = 1.0 - float.abs(x);
    var t: f64 = w * 0.5;
    var p: f64 = t *
        (1.66666666666666657415e-1 + t *
            (-3.25565818622400915405e-1 + t *
                (2.01212532134862925881e-1 + t *
                    (-4.00555345006794114027e-2 + t *
                        (7.91534994289814532176e-4 + t *
                            3.47933107596021167570e-5)))));
    var q: f64 = 1.0 + t *
        (-2.40339491173441421878e+0 + t *
            (2.02094576023350569471e+0 + t *
                (-6.88283971605453293030e-1 + t *
                    7.70381505559019352791e-2)));
    const s: f64 = float.sqrt(t);
    if (ix >= 0x3fef3333) { // If |x| > 0.975
        w = p / q;
        t = 1.57079632679489655800e+0 - (2.0 * (s + s * w) - 6.12323399573676603587e-17);
    } else {
        w = s;
        dbl64.setLowPart(&w, 0);
        const c: f64 = (t - w * w) / (s + w);
        const r: f64 = p / q;
        p = 2.0 * s * r - (6.12323399573676603587e-17 - 2.0 * c);
        q = 7.85398163397448278999e-1 - 2.0 * w;
        t = 7.85398163397448278999e-1 - (p - q);
    }

    if (hx > 0)
        return t
    else
        return -t;
}

fn asin80(x: f80) f80 {
    _ = x;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_asinl.c
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
fn asin128(x: f128) f128 {
    var u: ldbl128.ShapeSplit = .fromFloat(x);
    const sign: u1 = u.sign;
    const expt: i15 = @bitCast(u.exponent);

    if (expt >= (16384 - 1)) { // |x| >= 1
        if (expt == (16384 - 1) and (u.mantissa_high | u.mantissa_low) == 0) // asin(1) = ±pi/2 with inexact
            return x * 1.57079632679489661923132169163975140e+0 +
                x * 4.33590506506189051239852201302167613e-35;

        return (x - x) / (x - x); // asin(|x| > 1) is NaN
    } else if (expt < (16384 - 2)) { // |x| < 0.5
        if (expt < (16384 - 57)) { // if |x| is small, asinl(x) = x
            if (1.000e+300 + x > 1.0)
                return x; // Return x with inexact if x != 0
        }

        const t: f128 = x * x;
        const p: f128 = t *
            (1.66666666666666666666666666666700314e-1 + t *
                (-7.32816946414566252574527475428622708e-1 + t *
                    (1.34215708714992334609030036562143589e+0 + t *
                        (-1.32483151677116409805070261790752040e+0 + t *
                            (7.61206183613632558824485341162121989e-1 + t *
                                (-2.56165783329023486777386833928147375e-1 + t *
                                    (4.80718586374448793411019434585413855e-2 + t *
                                        (-4.42523267167024279410230886239774718e-3 + t *
                                            (1.44551535183911458253205638280410064e-4 + t *
                                                -2.10558957916600254061591040482706179e-7)))))))));
        const q: f128 = 1.0 + t *
            (-4.84690167848739751544716485245697428e+0 + t *
                (9.96619113536172610135016921140206980e+0 + t *
                    (-1.13177895428973036660836798461641458e+1 + t *
                        (7.74004374389488266169304117714658761e+0 + t *
                            (-3.25871986053534084709023539900339905e+0 + t *
                                (8.27830318881232209752469022352928864e-1 + t *
                                    (-1.18768052702942805423330715206348004e-1 + t *
                                        (8.32600764660522313269101537926539470e-3 + t *
                                            -1.99407384882605586705979504567947007e-4))))))));
        const w: f128 = p / q;
        return x + x * w;
    }

    // 1> |x| >= 0.5
    var w: f128 = 1.0 - float.abs(x);
    var t: f128 = w * 0.5;
    var p: f128 = t *
        (1.66666666666666666666666666666700314e-1 + t *
            (-7.32816946414566252574527475428622708e-1 + t *
                (1.34215708714992334609030036562143589e+0 + t *
                    (-1.32483151677116409805070261790752040e+0 + t *
                        (7.61206183613632558824485341162121989e-1 + t *
                            (-2.56165783329023486777386833928147375e-1 + t *
                                (4.80718586374448793411019434585413855e-2 + t *
                                    (-4.42523267167024279410230886239774718e-3 + t *
                                        (1.44551535183911458253205638280410064e-4 + t *
                                            -2.10558957916600254061591040482706179e-7)))))))));
    var q: f128 = 1.0 + t *
        (-4.84690167848739751544716485245697428e+0 + t *
            (9.96619113536172610135016921140206980e+0 + t *
                (-1.13177895428973036660836798461641458e+1 + t *
                    (7.74004374389488266169304117714658761e+0 + t *
                        (-3.25871986053534084709023539900339905e+0 + t *
                            (8.27830318881232209752469022352928864e-1 + t *
                                (-1.18768052702942805423330715206348004e-1 + t *
                                    (8.32600764660522313269101537926539470e-3 + t *
                                        -1.99407384882605586705979504567947007e-4))))))));
    const s: f128 = float.sqrt(t);
    if (u.mantissa_high >= (0xe666666666666666 >> 44)) { // If |x| is close to 1
        w = p / q;
        t = 1.57079632679489661923132169163975140e+0 -
            (2.0 * (s + s * w) - 4.33590506506189051239852201302167613e-35);
    } else {
        u = .fromFloat(s);
        u.mantissa_low = 0;
        w = u.toFloat();
        const c: f128 = (t - w * w) / (s + w);
        const r: f128 = p / q;
        p = 2.0 * s * r - (4.33590506506189051239852201302167613e-35 - 2.0 * c);
        q = 7.85398163397448309615660845819875699e-1 - 2.0 * w;
        t = 7.85398163397448309615660845819875699e-1 - (p - q);
    }

    if (sign == 0)
        return t
    else
        return -t;
}
