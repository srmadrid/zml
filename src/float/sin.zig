const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const cos = @import("cos.zig");
const k_cos32 = cos.k_cos32;
const k_cos64 = cos.k_cos64;
const k_cos128 = cos.k_cos128;
const rem_pio2 = @import("rem_pio2.zig");
const rem_pio2_32 = rem_pio2.rem_pio2_32;
const rem_pio2_64 = rem_pio2.rem_pio2_64;
const rem_pio2_128 = rem_pio2.rem_pio2_128;

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn sin(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.sin: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, sin32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_sinf.c
            return sin32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_sin.c
            return sin64(types.scast(f64, x));
        },
        f80 => {
            //
            // return sin80(types.scast(f80, x));
            return types.scast(f80, sin128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_sinl.c
            return sin128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_sinf.c
//
// Original copyright notice:
// s_sinf.c -- float version of s_sin.c.
// Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
// Optimized by Bruce D. Evans.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn sin32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix <= 0x3f490fda) { // |x| ~<= pi/4
        if (ix < 0x39800000) // |x| < 2**-12
            if (types.scast(i32, x) == 0)
                return x; // x with inexact if x != 0

        return k_sin32(types.scast(f64, x));
    }

    if (ix <= 0x407b53d1) { // |x| ~<= 5 * pi/4
        if (ix <= 0x4016cbe3) { // |x| ~<= 3 * pi/4
            if (hx > 0)
                return k_cos32(types.scast(f64, x) - 1.57079632679489661923)
            else
                return -k_cos32(types.scast(f64, x) + 1.57079632679489661923);
        } else {
            return k_sin32(@as(f64, if (hx > 0) 2.0 else -2.0) * 1.57079632679489661923 - types.scast(f64, x));
        }
    }

    if (ix <= 0x40e231d5) { // |x| ~<= 9 * pi/4
        if (ix <= 0x40afeddf) { // |x| ~<= 7 * pi/4
            if (hx > 0)
                return -k_cos32(types.scast(f64, x) - 3.0 * 1.57079632679489661923)
            else
                return k_cos32(types.scast(f64, x) + 3.0 * 1.57079632679489661923);
        } else {
            return k_sin32(types.scast(f64, x) + @as(f64, if (hx > 0) -4.0 else 4.0) * 1.57079632679489661923);
        }
    } else if (ix >= 0x7f800000) { // sin(Inf or NaN) is NaN
        return x - x;
    } else { // General argument reduction needed
        var y: f64 = undefined;
        const n: i32 = rem_pio2_32(x, &y);
        return switch (n & 3) {
            0 => k_sin32(y),
            1 => k_cos32(y),
            2 => k_sin32(-y),
            else => -k_cos32(y),
        };
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_sin.c
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
fn sin64(x: f64) f64 {
    const ix: i32 = @bitCast(dbl64.getHighPart(x) & 0x7fffffff);

    if (ix <= 0x3fe921fb) { // |x| ~< pi/4
        if (ix < 0x3e500000) { // |x| < 2**-26
            if (types.scast(i32, x) == 0) // Generate inexact
                return x;
        }

        return k_sin64(x, 0.0, 0);
    } else if (ix >= 0x7ff00000) { // sin(Inf or NaN) is NaN
        return x - x;
    } else {
        // Argument reduction needed
        var y: [2]f64 = undefined;
        const n: i32 = rem_pio2_64(x, &y);
        return switch (n & 3) {
            0 => k_sin64(y[0], y[1], 1),
            1 => k_cos64(y[0], y[1]),
            2 => -k_sin64(y[0], y[1], 1),
            else => -k_cos64(y[0], y[1]),
        };
    }
}

fn sin80(x: f80) f80 {
    _ = x;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_sinl.c
//
// Original copyright notice:
// Copyright (c) 2007 Steven G. Kargl
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice unmodified, this list of conditions, and the following
//    disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
fn sin128(x: f128) f128 {
    var z: ldbl128.ShapeSplit = .fromFloat(x);
    const s: u1 = z.sign;
    z.sign = 0;

    if (z.exponent == 0) // If x = ±0 or x is a subnormal number, then sin(x) = x
        return x;

    if (z.exponent == 32767) // If x = NaN or Inf, then sin(x) = NaN
        return ((x - x) / (x - x));

    if (z.toFloat() < 0.7853981633974483096156608458198757) { // Optimize the case where x is already within range
        const hi: f128 = k_sin128(z.toFloat(), 0, 0);
        return if (s != 0) -hi else hi;
    }

    var y: [2]f128 = undefined;
    const e0: i64 = rem_pio2_128(x, &y);
    var hi: f128 = y[0];
    const lo: f128 = y[1];

    switch (e0 & 3) {
        0 => hi = k_sin128(hi, lo, 1),
        1 => hi = k_cos128(hi, lo),
        2 => hi = -k_sin128(hi, lo, 1),
        3 => hi = -k_cos128(hi, lo),
        else => unreachable,
    }

    return hi;
}

pub fn k_sin32(x: f64) f32 {
    const z: f64 = x * x;
    const w: f64 = z * z;
    const r: f64 = -0x1a00f9e2cae774.0p-65 + z * 0x16cd878c3b46a7.0p-71;
    const s: f64 = z * x;
    return types.scast(f32, (x + s * (-0x15555554cbac77.0p-55 + z * 0x111110896efbb2.0p-59)) + s * w * r);
}

pub fn k_sin64(x: f64, y: f64, iy: i32) f64 {
    const z: f64 = x * x;
    const w: f64 = z * z;
    const r: f64 = 8.33333333332248946124e-3 +
        z * (-1.98412698298579493134e-4 + z * 2.75573137070700676789e-6) +
        z * w * (-2.50507602534068634195e-8 + z * 1.58969099521155010221e-10);
    const v: f64 = z * x;

    if (iy == 0)
        return x + v * (-1.66666666666666324348e-1 + z * r)
    else
        return x - ((z * (0.5 * y - v * r) - y) - v * -1.66666666666666324348e-1);
}

pub fn k_sin128(x: f128, y: f128, iy: i32) f128 {
    const z: f128 = x * x;
    const v: f128 = z * x;
    const r: f128 = 0.0083333333333333333333333333333331135404851288270047 + z *
        (-0.00019841269841269841269841269839935785325638310428717 + z *
            (0.27557319223985890652557316053039946268333231205686e-5 + z *
                (-0.25052108385441718775048214826384312253862930064745e-7 + z *
                    (0.16059043836821614596571832194524392581082444805729e-9 + z *
                        (-0.76471637318198151807063387954939213287488216303768e-12 + z *
                            (0.28114572543451292625024967174638477283187397621303e-14 + z *
                                (-0.82206352458348947812512122163446202498005154296863e-17 + z *
                                    (0.19572940011906109418080609928334380560135358385256e-19 + z *
                                        (-0.38680813379701966970673724299207480965452616911420e-22 + z *
                                            0.64038150078671872796678569586315881020659912139412e-25)))))))));

    if (iy == 0)
        return x + v * (-0.16666666666666666666666666666666666606732416116558 + z * r)
    else
        return x - ((z * (0.5 * y - v * r) - y) - v * -0.16666666666666666666666666666666666606732416116558);
}
