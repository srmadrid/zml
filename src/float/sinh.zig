const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Sinh(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.sinh: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

/// Returns the hyperbolic sine $\sinh(x)$ of a float, int or bool operand. The
/// result type is determined by coercing the operand type to a float, and the
/// operation is performed by casting the operand to the result type, then
/// computing its hyperbolic sine.
///
/// ## Signature
/// ```zig
/// float.sinh(x: X) float.Sinh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the hyperbolic sine of.
///
/// ## Returns
/// `float.Sinh(@TypeOf(x))`: The hyperbolic sine of `x`.
pub inline fn sinh(x: anytype) float.Sinh(@TypeOf(x)) {
    switch (float.Sinh(@TypeOf(x))) {
        f16 => return types.scast(f16, sinh32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_sinhf.c
            return sinh32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_sinh.c
            return sinh64(types.scast(f64, x));
        },
        f80 => {
            //
            // return sinh80(types.scast(f80, x));
            return types.scast(f80, sinh128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_sinhl.c
            return sinh128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_sinhf.c
//
// Original copyright notice:
// s_sinhf.c -- float version of s_sinh.c.
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
fn sinh32(x: f32) f32 {
    const jx: i32 = @bitCast(x);
    const ix: i32 = jx & 0x7fffffff;

    if (ix >= 0x7f800000) // x is Inf or NaN
        return x + x;

    var h: f32 = 0.5;
    if (jx < 0)
        h = -h;

    // |x| in [0,9], return sign(x) * 0.5 * (E + E/(E + 1)))
    if (ix < 0x41100000) { // |x| < 9
        if (ix < 0x39800000) // |x| < 2**-12
            return x;

        const t: f32 = float.expm1(float.abs(x));
        if (ix < 0x3f800000)
            return h * (2.0 * t - t * t / (t + 1.0));

        return h * (t + t / (t + 1.0));
    }

    if (ix < 0x42b17217) // |x| in [9, logf(maxfloat)] return 0.5 * exp(|x|)
        return h * float.exp(float.abs(x));

    if (ix <= 0x42b2d4fc) {
        // |x| in [logf(maxfloat), overflowthresold]
        var exp_x: f32 = float.exp(float.abs(x) - 162.88958740);
        const hx: u32 = @bitCast(exp_x);
        var expt: i32 = @bitCast((hx >> 23) -% (0x7f +% 127) +% 235);
        exp_x = @bitCast((hx & 0x7fffff) | ((0x7f +% 127) << 23));
        expt -%= 1;
        const scale: f32 = @bitCast((0x7f +% expt) << 23);
        return h * 2.0 * exp_x * scale;
    }

    // |x| > overflowthresold, sinh(x) overflow
    return if (jx < 0) -std.math.inf(f32) else std.math.inf(f32);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_sinh.c
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
fn sinh64(x: f64) f64 {
    const jx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = jx & 0x7fffffff;

    if (ix >= 0x7ff00000) // x is Inf or NaN
        return x + x;

    var h: f64 = 0.5;
    if (jx < 0)
        h = -h;

    // |x| in [0,22], return sign(x) * 0.5 * (E + E/(E + 1)))
    if (ix < 0x40360000) { // |x| < 22
        if (ix < 0x3e300000) // |x| < 2**-28
            return x;

        const t: f64 = float.expm1(float.abs(x));
        if (ix < 0x3ff00000)
            return h * (2.0 * t - t * t / (t + 1.0));

        return h * (t + t / (t + 1.0));
    }

    if (ix < 0x40862e42) // |x| in [22, log(maxdouble)] return 0.5 * exp(|x|)
        return h * float.exp(float.abs(x));

    if (ix <= 0x408633ce) { // |x| in [log(maxdouble), overflowthresold]
        var exp_x: f64 = float.exp(float.abs(x) - 1246.97177782734161156);
        const hx: u32 = @bitCast(dbl64.getHighPart(exp_x));
        var expt: i32 = @bitCast((hx >> 20) -% (0x3ff +% 1023) +% 1799);
        dbl64.setHighPart(&exp_x, (hx & 0xfffff) | ((0x3ff +% 1023) << 20));
        expt -%= 1;
        const scale: f64 = dbl64.Parts.toFloat(.{ .msw = @bitCast((0x3ff +% expt) << 20), .lsw = 0 });
        return h * 2.0 * exp_x * scale;
    }

    // |x| > overflowthresold, sinh(x) overflow
    return if (jx < 0) -std.math.inf(f64) else std.math.inf(f64);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_sinhl.c
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
fn sinh128(x: f128) f128 {
    var u: ldbl128.Parts32 = .fromFloat(x);
    const jx: u32 = u.mswhi;
    const ix: u32 = jx & 0x7fffffff;

    if (ix >= 0x7fff0000) // x is Inf or NaN
        return x + x;

    var h: f128 = 0.5;
    if (jx & 0x80000000 != 0)
        h = -h;

    u.mswhi = ix;

    if (ix <= 0x40044000) // |x| in [0,40], return sign(x) * 0.5 * (E + E/(E + 1)))
    {
        if (ix < 0x3fc60000) // |x| < 2**-57
            return x; // sinh(tiny) = tiny

        const t: f128 = float.expm1(u.toFloat());
        if (ix < 0x3fff0000)
            return h * (2.0 * t - t * t / (t + 1.0));

        return h * (t + t / (t + 1.0));
    }

    // |x| in [40, log(maxdouble)] return 0.5 * exp(|x|)
    if (ix <= 0x400c62e3) // 11356.375
        return h * float.exp(u.toFloat());

    // |x| in [log(maxdouble), overflowthreshold]
    // Overflow threshold is log(2 * maxdouble)
    if (u.toFloat() <= 1.1357216553474703894801348310092223067821e4) {
        const w: f128 = float.exp(0.5 * u.toFloat());
        const t: f128 = h * w;
        return t * w;
    }

    // |x| > overflowthreshold, sinhl(x) overflow
    return if (jx & 0x80000000 != 0) -std.math.inf(f128) else std.math.inf(f128);
}
