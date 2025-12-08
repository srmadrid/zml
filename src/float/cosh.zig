const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn cosh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.cosh: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, cosh32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_coshf.c
            return cosh32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_cosh.c
            return cosh64(types.scast(f64, x));
        },
        f80 => {
            //
            // return cosh80(types.scast(f80, x));
            return types.scast(f80, cosh128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_coshl.c
            return cosh128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_coshf.c
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
fn cosh32(x: f32) f32 {
    var ix: i32 = @bitCast(x);
    ix &= 0x7fffffff;

    if (ix >= 0x7f800000) // x is Inf or NaN
        return x * x;

    // |x| in [0, 0.5 * ln(2)], return 1 + expm1(|x|)**2/(2 * exp(|x|))
    if (ix < 0x3eb17218) {
        if (ix < 0x39800000) // cosh(tiny) = 1
            return 1.0;

        const t: f32 = float.expm1(float.abs(x));
        const w: f32 = 1.0 + t;
        return 1.0 + (t * t) / (2.0 * w);
    }

    // |x| in [0.5 * ln(2), 9], return (exp(|x|) + 1/exp(|x|))/2
    if (ix < 0x41100000) {
        const t: f32 = float.exp(float.abs(x));
        return 0.5 * t + 0.5 / t;
    }

    // |x| in [9, log(maxfloat)] return 0.5 * exp(|x|)
    if (ix < 0x42b17217)
        return 0.5 * float.exp(float.abs(x));

    if (ix <= 0x42b2d4fc) {
        // |x| in [logf(maxfloat), overflowthresold]
        var exp_x: f32 = float.exp(float.abs(x) - 162.88958740);
        const hx: u32 = @bitCast(exp_x);
        var expt: i32 = @bitCast((hx >> 23) -% (0x7f +% 127) +% 235);
        exp_x = @bitCast((hx & 0x7fffff) | ((0x7f +% 127) << 23));
        expt -%= 1;
        const scale: f32 = @bitCast((0x7f +% expt) << 23);
        return exp_x * scale;
    }

    // |x| > overflowthresold, cosh(x) overflow
    return std.math.inf(f32);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_cosh.c
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
fn cosh64(x: f64) f64 {
    var ix: i32 = @bitCast(dbl64.getHighPart(x));
    ix &= 0x7fffffff;

    if (ix >= 0x7ff00000) // x is Inf or NaN
        return x * x;

    // |x| in [0, 0.5 * ln(2)], return 1 + expm1(|x|)**2/(2 * exp(|x|))
    if (ix < 0x3fd62e43) {
        if (ix < 0x3c800000) // cosh(tiny) = 1
            return 1.0;

        const t: f64 = float.expm1(float.abs(x));
        const w: f64 = 1.0 + t;
        return 1.0 + (t * t) / (2.0 * w);
    }

    // |x| in [0.5 * ln(2), 22], return (exp(|x|) + 1/exp(|x|)/2
    if (ix < 0x40360000) {
        const t: f64 = float.exp(float.abs(x));
        return 0.5 * t + 0.5 / t;
    }

    // |x| in [22, log(maxdouble)] return 0.5 * exp(|x|)
    if (ix < 0x40862e42)
        return 0.5 * float.exp(float.abs(x));

    if (ix <= 0x408633ce) { // |x| in [log(maxdouble), overflowthresold]
        var exp_x: f64 = float.exp(float.abs(x) - 1246.97177782734161156);
        const hx: u32 = @bitCast(dbl64.getHighPart(exp_x));
        var expt: i32 = @bitCast((hx >> 20) -% (0x3ff +% 1023) +% 1799);
        dbl64.setHighPart(&exp_x, (hx & 0xfffff) | ((0x3ff +% 1023) << 20));
        expt -%= 1;
        const scale: f64 = dbl64.Parts.toFloat(.{ .msw = @bitCast((0x3ff +% expt) << 20), .lsw = 0 });
        return exp_x * scale;
    }

    // |x| > overflowthresold, cosh(x) overflow
    return std.math.inf(f64);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_coshl.c
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
fn cosh128(x: f128) f128 {
    var u: ldbl128.Parts32 = .fromFloat(x);
    const ex: i32 = @bitCast(u.mswhi & 0x7fffffff);
    u.mswhi = @bitCast(ex);

    // x is Inf or NaN
    if (ex >= 0x7fff0000)
        return x * x;

    // |x| in [0, 0.5 * ln(2)], return 1 + expm1l(|x|)**2/(2 * expl(|x|))
    if (ex < 0x3ffd62e4) // 0.3465728759765625
    {
        if (ex < 0x3fb80000) // coshl(tiny) = 1
            return 1.0;

        const t: f128 = float.expm1(u.toFloat());
        const w: f128 = 1.0 + t;
        return 1.0 + (t * t) / (w + w);
    }

    // |x| in [0.5 * ln(2), 40], return (exp(|x|) + 1/exp(|x|)/2
    if (ex < 0x40044000) {
        const t: f128 = float.exp(u.toFloat());
        return 0.5 * t + 0.5 / t;
    }

    // |x| in [22, ln(maxdouble)] return 0.5 * exp(|x|)
    if (ex <= 0x400c62e3) // 11356.375
        return 0.5 * float.exp(u.toFloat());

    // |x| in [log(maxdouble), overflowthresold]
    if (u.toFloat() <= 1.1357216553474703894801348310092223067821e4) {
        const w: f128 = float.exp(0.5 * u.toFloat());
        const t: f128 = 0.5 * w;
        return t * w;
    }

    // |x| > overflowthresold, coshl(x) overflow
    return std.math.inf(f128);
}
