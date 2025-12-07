const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const sin = @import("sin.zig");
const k_sin32 = sin.k_sin32;
const k_sin64 = sin.k_sin64;
const k_sin128 = sin.k_sin128;
const rem_pio2 = @import("rem_pio2.zig");
const rem_pio2_32 = rem_pio2.rem_pio2_32;
const rem_pio2_64 = rem_pio2.rem_pio2_64;
const rem_pio2_128 = rem_pio2.rem_pio2_128;

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn cos(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.cos: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, cos32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_cosf.c
            return cos32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_cos.c
            return cos64(types.scast(f64, x));
        },
        f80 => {
            //
            // return cos80(types.scast(f80, x));
            return types.scast(f80, cos128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_cosl.c
            return cos128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_cosf.c
//
// Original copyright notice:
// s_cosf.c -- float version of s_cos.c.
// Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
// Optimized by Bruce D. Evans.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. bucosess.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn cos32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix <= 0x3f490fda) { // |x| ~<= pi/4
        if (ix < 0x39800000) // |x| < 2**-12
            if (types.scast(i32, x) == 0)
                return 1.0; // 1 with inexact if x != 0

        return k_cos32(types.scast(f64, x));
    }

    if (ix <= 0x407b53d1) { // |x| ~<= 5 * pi/4
        if (ix <= 0x4016cbe3) { // |x| ~<= 3 * pi/4
            if (hx > 0)
                return k_sin32(1.57079632679489661923 - types.scast(f64, x))
            else
                return k_sin32(types.scast(f64, x) + 1.57079632679489661923);
        } else {
            return -k_cos32(types.scast(f64, x) + @as(f64, if (hx > 0) -2.0 else 2.0) * 1.57079632679489661923);
        }
    }

    if (ix <= 0x40e231d5) { // |x| ~<= 9 * pi/4
        if (ix <= 0x40afeddf) { // |x| ~<= 7 * pi/4
            if (hx > 0)
                return k_sin32(types.scast(f64, x) - 3.0 * 1.57079632679489661923)
            else
                return k_sin32(-3.0 * 1.57079632679489661923 - types.scast(f64, x));
        } else {
            return k_cos32(types.scast(f64, x) + @as(f64, if (hx > 0) -4.0 else 4.0) * 1.57079632679489661923);
        }
    } else if (ix >= 0x7f800000) { // cos(Inf or NaN) is NaN
        return x - x;
    } else { // General argument reduction needed
        var y: f64 = undefined;
        const n: i32 = rem_pio2_32(x, &y);
        return switch (n & 3) {
            0 => k_cos32(y),
            1 => k_sin32(-y),
            2 => -k_cos32(y),
            else => k_sin32(y),
        };
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_cos.c
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
fn cos64(x: f64) f64 {
    const ix: i32 = @bitCast(dbl64.getHighPart(x) & 0x7fffffff);

    if (ix <= 0x3fe921fb) { // |x| ~< pi/4
        if (ix < 0x3e46a09e) { // x < 2**-27 * sqrt(2)
            if (types.scast(i32, x) == 0) // Generate inexact
                return 1.0;
        }

        return k_cos64(x, 0.0);
    } else if (ix >= 0x7ff00000) { // cos(Inf or NaN) is NaN
        return x - x;
    } else {
        // Argument reduction needed
        var y: [2]f64 = undefined;
        const n: i32 = rem_pio2_64(x, &y);
        return switch (n & 3) {
            0 => k_cos64(y[0], y[1]),
            1 => -k_sin64(y[0], y[1], 1),
            2 => -k_cos64(y[0], y[1]),
            else => k_sin64(y[0], y[1], 1),
        };
    }
}

fn cos80(x: f80) f80 {
    _ = x;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_cosl.c
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
fn cos128(x: f128) f128 {
    var z: ldbl128.ShapeSplit = .fromFloat(x);
    z.sign = 0;

    if (z.exponent == 0) // If x = ±0 or x is a subnormal number, then cos(x) = 1
        return 1.0;

    if (z.exponent == 32767) // If x = NaN or Inf, then cos(x) = NaN
        return ((x - x) / (x - x));

    if (z.toFloat() < 0.7853981633974483096156608458198757) // Optimize the case where x is already within range
        return k_cos128(z.toFloat(), 0.0);

    var y: [2]f128 = undefined;
    const e0: i64 = rem_pio2_128(x, &y);
    var hi: f128 = y[0];
    const lo: f128 = y[1];

    switch (e0 & 3) {
        0 => hi = k_cos128(hi, lo),
        1 => hi = -k_sin128(hi, lo, 1),
        2 => hi = -k_cos128(hi, lo),
        3 => hi = k_sin128(hi, lo, 1),
        else => unreachable,
    }

    return hi;
}

pub fn k_cos32(x: f64) f32 {
    const z: f64 = x * x;
    const w: f64 = z * z;
    const r: f64 = -0x16c087e80f1e27.0p-62 + z * 0x199342e0ee5069.0p-68;
    return types.scast(f32, ((1.0 + z * -0x1ffffffd0c5e81.0p-54) + w * 0x155553e1053a42.0p-57) + (w * z) * r);
}

pub fn k_cos64(x: f64, y: f64) f64 {
    const z: f64 = x * x;
    var w: f64 = z * z;
    const r: f64 = z * (4.16666666666666019037e-2 + z *
        (-1.38888888888741095749e-3 + z * 2.48015872894767294178e-5)) +
        w * w * (-2.75573143513906633035e-7 +
            z * (2.08757232129817482790e-9 + z * -1.13596475577881948265e-11));
    const hz: f64 = 0.5 * z;
    w = 1.0 - hz;
    return w + (((1.0 - w) - hz) + (z * r - x * y));
}

pub fn k_cos128(x: f128, y: f128) f128 {
    const z: f128 = x * x;
    const r: f128 = z *
        (0.04166666666666666666666666666666658424671 + z *
            (-0.001388888888888888888888888888863490893732 + z *
                (0.00002480158730158730158730158600795304914210 + z *
                    (-0.2755731922398589065255474947078934284324e-6 + z *
                        (0.2087675698786809897659225313136400793948e-8 + z *
                            (-0.1147074559772972315817149986812031204775e-10 + z *
                                (0.4779477332386808976875457937252120293400e-13 + z *
                                    (-0.1561920696721507929516718307820958119868e-15 + z *
                                        (0.4110317413744594971475941557607804508039e-18 + z *
                                            (-0.8896592467191938803288521958313920156409e-21 + z *
                                                0.1601061435794535138244346256065192782581e-23))))))))));
    const hz: f128 = 0.5 * z;
    const w: f128 = 1.0 - hz;
    return w + (((1.0 - w) - hz) + (z * r - x * y));
}
