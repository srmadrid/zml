const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn asinh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.asinh: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, asinh32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_asinhf.c
            return asinh32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_asinh.c
            return asinh64(types.scast(f64, x));
        },
        f80 => {
            //
            // return asinh80(types.scast(f80, x));
            return types.scast(f80, asinh128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/s_asinhl.c
            return asinh128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_asinhf.c
//
// Original copyright notice:
// s_asinhf.c -- float version of s_asinh.c.
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
fn asinh32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x7f800000) // x is inf or NaN
        return x + x;

    if (ix < 0x31800000) // |x| < 2**-28
        return x;

    var w: f32 = undefined;
    if (ix > 0x4d800000) { // |x| > 2**28
        w = float.log(float.abs(x)) + 6.9314718246e-1;
    } else if (ix > 0x40000000) { // 2**28 > |x| > 2.0
        const t: f32 = float.abs(x);
        w = float.log(2.0 * t + 1.0 / (float.sqrt(x * x + 1.0) + t));
    } else { // 2.0 > |x| > 2**-28
        const t: f32 = x * x;
        w = float.log1p(float.abs(x) + t / (1.0 + float.sqrt(1.0 + t)));
    }

    return if (hx > 0) w else -w;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_asinh.c
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
fn asinh64(x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = hx & 0x7fffffff;

    if (ix >= 0x7ff00000) // x is inf or NaN
        return x + x;

    if (ix < 0x3e300000) // |x| < 2**-28
        return x;

    var w: f64 = undefined;
    if (ix > 0x41b00000) { // |x| > 2**28
        w = float.log(float.abs(x)) + 6.93147180559945286227e-1;
    } else if (ix > 0x40000000) { // 2**28 > |x| > 2.0
        const t: f64 = float.abs(x);
        w = float.log(2.0 * t + 1.0 / (float.sqrt(x * x + 1.0) + t));
    } else { // 2.0 > |x| > 2**-28
        const t: f64 = x * x;
        w = float.log1p(float.abs(x) + t / (1.0 + float.sqrt(1.0 + t)));
    }

    return if (hx > 0) w else -w;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/s_asinhl.c
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
fn asinh128(x: f128) f128 {
    var u: ldbl128.Parts32 = .fromFloat(x);
    const sign: i32 = @bitCast(u.mswhi);
    const ix: i32 = sign & 0x7fffffff;

    if (ix == 0x7fff0000)
        return x + x; // x is inf or NaN

    if (ix < 0x3fc70000) // |x| < 2**-56
        return x;

    u.mswhi = @bitCast(ix);
    var w: f128 = undefined;
    if (ix > 0x40350000) { // |x| > 2**54
        w = float.log(u.toFloat()) + 6.931471805599453094172321214581765681e-1;
    } else if (ix > 0x40000000) { // 2**54 > |x| > 2.0
        const t: f128 = u.toFloat();
        w = float.log(2.0 * t + 1.0 / (float.sqrt(x * x + 1.0) + t));
    } else { // 2.0 > |x| > 2**-56
        const t: f128 = x * x;
        w = float.log1p(u.toFloat() + t / (1.0 + float.sqrt(1.0 + t)));
    }

    return if (sign > 0) w else -w;
}
