const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn tanh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.tanh: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, tanh32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_tanhf.c
            return tanh32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_tanh.c
            return tanh64(types.scast(f64, x));
        },
        f80 => {
            //
            // return tanh80(types.scast(f80, x));
            return types.scast(f80, tanh128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_tanhl.c
            return tanh128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_tanhf.c
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
fn tanh32(x: f32) f32 {
    const jx: i32 = @bitCast(x);
    const ix: i32 = jx & 0x7fffffff;

    if (ix >= 0x7f800000) { // x is Inf or NaN
        if (jx >= 0)
            return 1.0 / x + 1.0 // tanh(±Inf) = ±1
        else
            return 1.0 / x - 1.0; // tanh(NaN) = NaN
    }

    var z: f32 = undefined;
    if (ix < 0x41100000) { // |x| < 9
        if (ix < 0x39800000) // |x| < 2**-12
            return x; // tanh(tiny) = tiny

        if (ix >= 0x3f800000) { // |x| >= 1
            const t: f32 = float.expm1(2.0 * float.abs(x));
            z = 1.0 - 2.0 / (t + 2.0);
        } else {
            const t: f32 = float.expm1(-2.0 * float.abs(x));
            z = -t / (t + 2.0);
        }
    } else { // |x| >= 9, return ±1
        z = 1.0;
    }

    return if (jx >= 0) z else -z;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_tanh.c
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
fn tanh64(x: f64) f64 {
    const jx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = jx & 0x7fffffff;

    if (ix >= 0x7ff00000) { // x is Inf or NaN
        if (jx >= 0)
            return 1.0 / x + 1.0 // tanh(±Inf) = ±1
        else
            return 1.0 / x - 1.0; // tanh(NaN) = NaN
    }

    var z: f64 = undefined;
    if (ix < 0x40360000) { // |x| < 22
        if (ix < 0x3e300000) // |x| < 2**-28
            return x; // tanh(tiny) = tiny

        if (ix >= 0x3ff00000) { // |x| >= 1
            const t: f64 = float.expm1(2.0 * float.abs(x));
            z = 1.0 - 2.0 / (t + 2.0);
        } else {
            const t: f64 = float.expm1(-2.0 * float.abs(x));
            z = -t / (t + 2.0);
        }
    } else { // |x| >= 22, return ±1
        z = 1.0;
    }

    return if (jx >= 0) z else -z;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_tanhl.c
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
fn tanh128(x: f128) f128 {
    var u: ldbl128.Parts32 = .fromFloat(x);
    const jx: u32 = u.mswhi;
    const ix: u32 = jx & 0x7fffffff;

    if (ix >= 0x7fff0000) {
        // for NaN it's not important which branch: tanhl(NaN) = NaN
        if (jx & 0x80000000 != 0)
            return 1.0 / x - 1.0 // tanhl(-inf) = -1
        else
            return 1.0 / x + 1.0; // tanhl(+inf) = +1
    }

    // |x| < 40
    var z: f128 = undefined;
    if (ix < 0x40044000) {
        if (u.toFloat() == 0)
            return x; // x == ±0

        if (ix < 0x3fc60000) // |x| < 2**-57
            return x; // tanh(small) = small

        u.mswhi = ix; // Absolute value of x
        if (ix >= 0x3fff0000) { // |x| >= 1
            const t: f128 = float.expm1(2.0 * u.toFloat());
            z = 1.0 - 2.0 / (t + 2.0);
        } else {
            const t: f128 = float.expm1(-2.0 * u.toFloat());
            z = -t / (t + 2.0);
        }
    } else { // |x| > 40, return ±1
        z = 1.0;
    }

    return if (jx & 0x80000000 != 0) -z else z;
}
