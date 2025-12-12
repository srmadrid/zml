const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn acosh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.acosh: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, acosh32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_acoshf.c
            return acosh32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_acosh.c
            return acosh64(types.scast(f64, x));
        },
        f80 => {
            //
            // return acosh80(types.scast(f80, x));
            return types.scast(f80, acosh128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_acoshl.c
            return acosh128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_acoshf.c
//
// Original copyright notice:
// s_acoshf.c -- float version of s_acosh.c.
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
fn acosh32(x: f32) f32 {
    const hx: i32 = @bitCast(x);

    if (hx < 0x3f800000) { // x < 1
        return (x - x) / (x - x);
    } else if (hx >= 0x4d800000) { // x > 2**28
        if (hx >= 0x7f800000) // x is inf of NaN
            return x + x
        else
            return float.log(x) + 6.9314718246e-1; // acosh(huge) = log(2 * x)
    } else if (hx == 0x3f800000) {
        return 0.0; // acosh(1) = 0
    } else if (hx > 0x40000000) { // 2**28 > x > 2
        const t: f32 = x * x;
        return float.log(2.0 * x - 1.0 / (x + float.sqrt(t - 1.0)));
    } else { // 1 < x < 2
        const t: f32 = x - 1.0;
        return float.log1p(t + float.sqrt(2.0 * t + t * t));
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_acosh.c
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
fn acosh64(x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const lx: i32 = @bitCast(dbl64.getLowPart(x));

    if (hx < 0x3ff00000) { // x < 1
        return (x - x) / (x - x);
    } else if (hx >= 0x41b00000) { // x > 2**28
        if (hx >= 0x7ff00000) // x is inf of NaN
            return x + x
        else
            return float.log(x) + 6.93147180559945286227e-1; // acosh(huge) = log(2 * x)
    } else if (((hx - 0x3ff00000) | lx) == 0) {
        return 0.0; // acosh(1) = 0
    } else if (hx > 0x40000000) { // 2**28 > x > 2
        const t: f64 = x * x;
        return float.log(2.0 * x - 1.0 / (x + float.sqrt(t - 1.0)));
    } else { // 1 < x < 2
        const t: f64 = x - 1.0;
        return float.log1p(t + float.sqrt(2.0 * t + t * t));
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_acoshl.c
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
fn acosh128(x: f128) f128 {
    const hx: i64 = @bitCast(ldbl128.getHighPart(x));
    const lx: i64 = @bitCast(ldbl128.getLowPart(x));

    if (hx < 0x3fff000000000000) { // x < 1
        return (x - x) / (x - x);
    } else if (hx >= 0x4035000000000000) { // x > 2**54
        if (hx >= 0x7fff000000000000) // x is inf of NaN
            return x + x
        else
            return float.log(x) + 0.6931471805599453094172321214581766; // acoshl(huge) = logl(2 * x)
    } else if (((hx - 0x3fff000000000000) | lx) == 0) {
        return 0.0; // acosh(1) = 0
    } else if (hx > 0x4000000000000000) { // 2**28 > x > 2
        const t: f128 = x * x;
        return float.log(2.0 * x - 1.0 / (x + float.sqrt(t - 1.0)));
    } else { // 1 < x < 2
        const t: f128 = x - 1.0;
        return float.log1p(t + float.sqrt(2.0 * t + t * t));
    }
}
