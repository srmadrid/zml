const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn acos(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.acos: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, acos32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_acosf.c
            return acos32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_acos.c
            return acos64(types.scast(f64, x));
        },
        f80 => {
            //
            // return acos80(types.scast(f80, x));
            return types.scast(f80, acos128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_acosl.c
            return acos128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_acosf.c
//
// Original copyright notice:
// e_acosf.c -- float version of e_acos.c.
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
fn acos32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x3f800000) { // |x| >= 1
        if (ix == 0x3f800000) { // |x| == 1
            if (hx > 0)
                return 0.0 // acos(1) = 0
            else
                return 3.1415925026e+0 + 2.0 * 7.5497894159e-8; // acos(-1) = pi
        }

        return (x - x) / (x - x); // acos(|x| > 1) is NaN
    }

    if (ix < 0x3f000000) { // |x| < 0.5
        if (ix <= 0x32800000) // If|x| < 2**-26
            return 1.5707962513e+0 + 7.5497894159e-8;

        const z: f32 = x * x;
        const p: f32 = z * (1.6666586697e-1 +
            z * (-4.2743422091e-2 + z * -8.6563630030e-3));
        const q: f32 = 1.0 + z * -7.0662963390e-1;
        const r: f32 = p / q;
        return 1.5707962513e+0 - (x - (7.5497894159e-8 - x * r));
    } else if (hx < 0) { // x < -0.5
        const z: f32 = (1.0 + x) * 0.5;
        const p: f32 = z * (1.6666586697e-1 +
            z * (-4.2743422091e-2 + z * -8.6563630030e-3));
        const q: f32 = 1.0 + z * -7.0662963390e-1;
        const s: f32 = float.sqrt(z);
        const r: f32 = p / q;
        const w: f32 = r * s - 7.5497894159e-8;
        return 3.1415925026e+0 - 2.0 * (s + w);
    } else { // x > 0.5
        const z: f32 = (1.0 - x) * 0.5;
        const s: f32 = float.sqrt(z);
        var df: f32 = s;
        const idf: u32 = @bitCast(df);
        df = @bitCast(idf & 0xfffff000);
        const c: f32 = (z - df * df) / (s + df);
        const p: f32 = z * (1.6666586697e-1 +
            z * (-4.2743422091e-2 + z * -8.6563630030e-3));
        const q: f32 = 1.0 + z * -7.0662963390e-1;
        const r: f32 = p / q;
        const w: f32 = r * s + c;
        return 2.0 * (df + w);
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_acos.c
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
fn acos64(x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x3ff00000) { // |x| >= 1
        const lx: i32 = @bitCast(dbl64.getLowPart(x));
        if (((ix - 0x3ff00000) | lx) == 0) { // |x| == 1
            if (hx > 0)
                return 0.0 // acos(1) = 0
            else
                return 3.14159265358979311600e+0 + 6.12323399573676603587e-17; // acos(-1) = pi
        }

        return (x - x) / (x - x); // acos(|x| > 1) is NaN
    }

    if (ix < 0x3fe00000) { // |x| < 0.5
        if (ix <= 0x3c600000) // |x| < 2**-57
            return 1.57079632679489655800e+0 + 6.12323399573676603587e-17;

        const z: f64 = x * x;
        const p: f64 = z *
            (1.66666666666666657415e-1 + z *
                (-3.25565818622400915405e-1 + z *
                    (2.01212532134862925881e-1 + z *
                        (-4.00555345006794114027e-2 + z *
                            (7.91534994289814532176e-4 + z *
                                3.47933107596021167570e-5)))));
        const q: f64 = 1.0 + z *
            (-2.40339491173441421878e+0 + z *
                (2.02094576023350569471e+0 + z *
                    (-6.88283971605453293030e-1 + z *
                        7.70381505559019352791e-2)));
        const r: f64 = p / q;
        return 1.57079632679489655800e+0 - (x - (6.12323399573676603587e-17 - x * r));
    } else if (hx < 0) { // x < -0.5
        const z: f64 = (1.0 + x) * 0.5;
        const p: f64 = z *
            (1.66666666666666657415e-1 + z *
                (-3.25565818622400915405e-1 + z *
                    (2.01212532134862925881e-1 + z *
                        (-4.00555345006794114027e-2 + z *
                            (7.91534994289814532176e-4 + z *
                                3.47933107596021167570e-5)))));
        const q: f64 = 1.0 + z *
            (-2.40339491173441421878e+0 + z *
                (2.02094576023350569471e+0 + z *
                    (-6.88283971605453293030e-1 + z *
                        7.70381505559019352791e-2)));
        const s: f64 = float.sqrt(z);
        const r: f64 = p / q;
        const w: f64 = r * s - 6.12323399573676603587e-17;
        return 3.14159265358979311600e+0 - 2.0 * (s + w);
    } else { // x > 0.5
        const z: f64 = (1.0 - x) * 0.5;
        const s: f64 = float.sqrt(z);
        var df: f64 = s;
        dbl64.setLowPart(&df, 0);
        const c: f64 = (z - df * df) / (s + df);
        const p: f64 = z *
            (1.66666666666666657415e-1 + z *
                (-3.25565818622400915405e-1 + z *
                    (2.01212532134862925881e-1 + z *
                        (-4.00555345006794114027e-2 + z *
                            (7.91534994289814532176e-4 + z *
                                3.47933107596021167570e-5)))));
        const q: f64 = 1.0 + z *
            (-2.40339491173441421878e+0 + z *
                (2.02094576023350569471e+0 + z *
                    (-6.88283971605453293030e-1 + z *
                        7.70381505559019352791e-2)));
        const r: f64 = p / q;
        const w: f64 = r * s + c;
        return 2.0 * (df + w);
    }
}

fn acos80(x: f80) f80 {
    _ = x;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_acosl.c
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
fn acos128(x: f128) f128 {
    var u: ldbl128.ShapeSplit = .fromFloat(x);
    const sign: u1 = u.sign;
    const expt: i15 = @bitCast(u.exponent);

    if (expt >= (16384 - 1)) { // |x| >= 1
        if (expt == (16384 - 1) and (u.mantissa_high | u.mantissa_low) == 0) { // |x| == 1
            if (sign == 0)
                return 0.0 // acos(1) = 0
            else
                return 3.14159265358979323846264338327950280e+0 + 2.0 * 8.67181013012378102479704402604335225e-35; // acos(-1) = pi
        }

        return (x - x) / (x - x); // acos(|x| > 1) is NaN
    } else if (expt < (16384 - 2)) { // |x| < 0.5
        if (expt < (16384 - 114)) // if |x| is small, acos(x) = pi/2
            return 1.57079632679489661923132169163975140e+0 + 4.33590506506189051239852201302167613e-35;

        const z: f128 = x * x;
        const p: f128 = z *
            (1.66666666666666666666666666666700314e-1 + z *
                (-7.32816946414566252574527475428622708e-1 + z *
                    (1.34215708714992334609030036562143589e+0 + z *
                        (-1.32483151677116409805070261790752040e+0 + z *
                            (7.61206183613632558824485341162121989e-1 + z *
                                (-2.56165783329023486777386833928147375e-1 + z *
                                    (4.80718586374448793411019434585413855e-2 + z *
                                        (-4.42523267167024279410230886239774718e-3 + z *
                                            (1.44551535183911458253205638280410064e-4 + z *
                                                -2.10558957916600254061591040482706179e-7)))))))));
        const q: f128 = 1.0 + z *
            (-4.84690167848739751544716485245697428e+0 + z *
                (9.96619113536172610135016921140206980e+0 + z *
                    (-1.13177895428973036660836798461641458e+1 + z *
                        (7.74004374389488266169304117714658761e+0 + z *
                            (-3.25871986053534084709023539900339905e+0 + z *
                                (8.27830318881232209752469022352928864e-1 + z *
                                    (-1.18768052702942805423330715206348004e-1 + z *
                                        (8.32600764660522313269101537926539470e-3 + z *
                                            -1.99407384882605586705979504567947007e-4))))))));
        const r: f128 = p / q;
        return 1.57079632679489661923132169163975140e+0 - (x - (4.33590506506189051239852201302167613e-35 - x * r));
    } else if (sign == 1) { // x < -0.5
        const z: f128 = (1.0 + x) * 0.5;
        const p: f128 = z *
            (1.66666666666666666666666666666700314e-1 + z *
                (-7.32816946414566252574527475428622708e-1 + z *
                    (1.34215708714992334609030036562143589e+0 + z *
                        (-1.32483151677116409805070261790752040e+0 + z *
                            (7.61206183613632558824485341162121989e-1 + z *
                                (-2.56165783329023486777386833928147375e-1 + z *
                                    (4.80718586374448793411019434585413855e-2 + z *
                                        (-4.42523267167024279410230886239774718e-3 + z *
                                            (1.44551535183911458253205638280410064e-4 + z *
                                                -2.10558957916600254061591040482706179e-7)))))))));
        const q: f128 = 1.0 + z *
            (-4.84690167848739751544716485245697428e+0 + z *
                (9.96619113536172610135016921140206980e+0 + z *
                    (-1.13177895428973036660836798461641458e+1 + z *
                        (7.74004374389488266169304117714658761e+0 + z *
                            (-3.25871986053534084709023539900339905e+0 + z *
                                (8.27830318881232209752469022352928864e-1 + z *
                                    (-1.18768052702942805423330715206348004e-1 + z *
                                        (8.32600764660522313269101537926539470e-3 + z *
                                            -1.99407384882605586705979504567947007e-4))))))));
        const s: f128 = float.sqrt(z);
        const r: f128 = p / q;
        const w: f128 = r * s - 4.33590506506189051239852201302167613e-35;
        return 3.14159265358979323846264338327950280e+0 - 2.0 * (s + w);
    } else { // x > 0.5
        const z: f128 = (1.0 - x) * 0.5;
        const s: f128 = float.sqrt(z);
        u = .fromFloat(s);
        u.mantissa_low = 0;
        const df: f128 = u.toFloat();
        const c: f128 = (z - df * df) / (s + df);
        const p: f128 = z *
            (1.66666666666666666666666666666700314e-1 + z *
                (-7.32816946414566252574527475428622708e-1 + z *
                    (1.34215708714992334609030036562143589e+0 + z *
                        (-1.32483151677116409805070261790752040e+0 + z *
                            (7.61206183613632558824485341162121989e-1 + z *
                                (-2.56165783329023486777386833928147375e-1 + z *
                                    (4.80718586374448793411019434585413855e-2 + z *
                                        (-4.42523267167024279410230886239774718e-3 + z *
                                            (1.44551535183911458253205638280410064e-4 + z *
                                                -2.10558957916600254061591040482706179e-7)))))))));
        const q: f128 = 1.0 + z *
            (-4.84690167848739751544716485245697428e+0 + z *
                (9.96619113536172610135016921140206980e+0 + z *
                    (-1.13177895428973036660836798461641458e+1 + z *
                        (7.74004374389488266169304117714658761e+0 + z *
                            (-3.25871986053534084709023539900339905e+0 + z *
                                (8.27830318881232209752469022352928864e-1 + z *
                                    (-1.18768052702942805423330715206348004e-1 + z *
                                        (8.32600764660522313269101537926539470e-3 + z *
                                            -1.99407384882605586705979504567947007e-4))))))));
        const r: f128 = p / q;
        const w: f128 = r * s + c;
        return 2.0 * (df + w);
    }
}
