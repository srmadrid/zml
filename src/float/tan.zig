const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const rem_pio2 = @import("rem_pio2.zig");
const rem_pio2_32 = rem_pio2.rem_pio2_32;
const rem_pio2_64 = rem_pio2.rem_pio2_64;
const rem_pio2_128 = rem_pio2.rem_pio2_128;

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn tan(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.tan: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, tan32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_tanf.c
            return tan32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_tan.c
            return tan64(types.scast(f64, x));
        },
        f80 => {
            //
            // return tan80(types.scast(f80, x));
            return types.scast(f80, tan128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_tanl.c
            return tan128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_tanf.c
//
// Original copyright notice:
// s_tanf.c -- float version of s_tan.c.
// Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
// Optimized by Bruce D. Evans.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. butaness.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn tan32(x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix <= 0x3f490fda) { // |x| ~<= pi/4
        if (ix < 0x39800000) // |x| < 2**-12
            return x;

        return k_tan32(x, 1);
    }

    if (ix <= 0x407b53d1) { // |x| ~<= 5 * pi/4
        if (ix <= 0x4016cbe3) // |x| ~<= 3 * pi/4
            return k_tan32(x + @as(f64, if (hx > 0) -1.0 else 1.0) * 1.57079632679489661923, -1)
        else
            return k_tan32(x + @as(f64, if (hx > 0) -2.0 else 2.0) * 1.57079632679489661923, 1);
    }

    if (ix <= 0x40e231d5) { // |x| ~<= 9 * pi/4
        if (ix <= 0x40afeddf) // |x| ~<= 7 * pi/4
            return k_tan32(x + @as(f64, if (hx > 0) -3.0 else 3.0) * 1.57079632679489661923, -1)
        else
            return k_tan32(x + @as(f64, if (hx > 0) -4.0 else 4.0) * 1.57079632679489661923, 1);
    } else if (ix >= 0x7f800000) { // tan(Inf or NaN) is NaN
        return x - x;
    } else {
        // General argument reduction needed
        var y: f64 = undefined;
        const n: i32 = rem_pio2_32(x, &y);

        // Integer parameter: 1 if n even, -1 if n odd
        return k_tan32(y, 1 - ((n & 1) << 1));
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_tan.c
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
fn tan64(x: f64) f64 {
    const ix: i32 = @bitCast(dbl64.getHighPart(x) & 0x7fffffff);

    if (ix <= 0x3fe921fb) { // |x| ~< pi/4
        if (ix < 0x3e400000) // x < 2**-27
            return x;

        return k_tan64(x, 0.0, 1);
    } else if (ix >= 0x7ff00000) { // tan(Inf or NaN) is NaN
        return x - x;
    } else {
        // General argument reduction needed
        var y: [2]f64 = undefined;
        const n: i32 = rem_pio2_64(x, &y);

        // Integer parameter: 1 if n even, -1 if n odd
        return k_tan64(y[0], y[1], 1 - ((n & 1) << 1));
    }
}

fn tan80(x: f80) f80 {
    _ = x;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_tanl.c
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
fn tan128(x: f128) f128 {
    var z: ldbl128.ShapeSplit = .fromFloat(x);
    const s: u1 = z.sign;
    z.sign = 0;

    if (z.exponent == 0) // If x = ±0 or x is subnormal, then tan(x) = x
        return x;

    if (z.exponent == 32767) // If x = NaN or Inf, then tan(x) = NaN
        return (x - x) / (x - x);

    if (z.toFloat() < 0.7853981633974483096156608458198757) { // Optimize the case where x is already within range
        const hi: f128 = k_tan128(z.toFloat(), 0.0, 0);
        return if (s != 0) -hi else hi;
    }

    var y: [2]f128 = undefined;
    const e0: i64 = rem_pio2_128(x, &y);
    var hi: f128 = y[0];
    const lo: f128 = y[1];

    switch (e0 & 3) {
        0, 2 => hi = k_tan128(hi, lo, 0),
        1, 3 => hi = k_tan128(hi, lo, 1),
        else => unreachable,
    }

    return hi;
}

pub fn k_tan32(x: f64, iy: i32) f32 {
    const z: f64 = types.scast(f64, x) * types.scast(f64, x);
    var r: f64 = 0x185dadfcecf44e.0p-61 + z * 0x1362b9bf971bcd.0p-59;
    const t: f64 = 0x1b54c91d865afe.0p-57 + z * 0x191df3908c33ce.0p-58;
    const w: f64 = z * z;
    const s: f64 = z * types.scast(f64, x);
    const u: f64 = 0x15554d3418c99f.0p-54 + z * 0x1112fd38999f72.0p-55;
    r = (types.scast(f64, x) + s * u) + (s * w) * (t + w * r);
    if (iy == 1)
        return types.scast(f32, r)
    else
        return types.scast(f32, -1.0 / r);
}

pub fn k_tan64(x: f64, y: f64, iy: i32) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = hx & 0x7fffffff;

    var xx: f64 = x;
    var yy: f64 = y;
    if (ix >= 0x3fe59428) { // |x| >= 0.6744
        if (hx < 0) {
            xx = -x;
            yy = -y;
        }

        const z: f64 = 7.85398163397448278999e-1 - xx;
        const w: f64 = 3.06161699786838301793e-17 - yy;
        xx = z + w;
        yy = 0.0;
    }

    var z: f64 = xx * xx;
    var w: f64 = z * z;
    var r: f64 = 1.33333333333201242699e-1 + w *
        (2.18694882948595424599e-2 + w *
            (3.59207910759131235356e-3 + w *
                (5.88041240820264096874e-4 + w *
                    (7.81794442939557092300e-5 + w *
                        -1.85586374855275456654e-5))));
    var v: f64 = z *
        (5.39682539762260521377e-2 + w *
            (8.86323982359930005737e-3 + w *
                (1.45620945432529025516e-3 + w *
                    (2.46463134818469906812e-4 + w *
                        (7.14072491382608190305e-5 + w *
                            2.59073051863633712884e-5)))));
    var s: f64 = z * xx;
    r = yy + z * (s * (r + v) + yy);
    r += 3.33333333333334091986e-1 * s;
    w = xx + r;
    if (ix >= 0x3fe59428) {
        v = types.scast(f64, iy);
        return types.scast(f64, 1 - ((hx >> 30) & 2)) *
            (v - 2.0 * (xx - (w * w / (w + v) - r)));
    }

    if (iy == 1) {
        return w;
    } else {
        // If allow error up to 2 ulp, simply return
        // -1.0 / (x+r) here
        // Compute -1.0 / (x+r) accurately
        z = w;
        dbl64.setLowPart(&z, 0);
        v = r - (z - xx); // z + v = r + x
        const a: f64 = -1.0 / w;
        var t: f64 = -1.0 / w;
        dbl64.setLowPart(&t, 0);
        s = 1.0 + t * z;
        return t + a * (s + t * v);
    }
}

pub fn k_tan128(x: f128, y: f128, iy: i32) f128 {
    const iiy: i32 = if (iy == 1) -1 else 1; // Recover original interface
    const osign: f128 = if (x >= 0) 1.0 else -1.0; // Slow, probably wrong for -0

    var xx: f128 = x;
    var yy: f128 = y;
    var i: i32 = 0;
    if (float.abs(xx) >= 0.67434) {
        if (xx < 0) {
            xx = -xx;
            yy = -yy;
        }

        const z: f128 = 0x1.921fb54442d18469898cc51701b8p-1 - xx;
        const w: f128 = 0x1.cd129024e088a67cc74020bbea60p-116 - yy;
        xx = z + w;
        yy = 0.0;
        i = 1;
    }

    var z: f128 = xx * xx;
    var w: f128 = z * z;
    var r: f128 = 0x1.1111111111111111111111111eb5p-3 + w *
        (0x1.664f4882c10f9f32d6bbe09d8bcdp-6 + w *
            (0x1.d6d3d0e157ddfb5fed8e84e27b37p-9 + w *
                (0x1.355824803674477dfcf726649efep-11 + w *
                    (0x1.967e18afcb180ed942dfdc518d6cp-14 + w *
                        (0x1.0b132d39f055c81be49eff7afd50p-16 + w *
                            (0x1.5ef2daf21d1113df38d0fbc00267p-19 + w *
                                (0x1.cd2a5a292b180e0bdd701057dfe3p-22 + w *
                                    (0x1.2f3190f4718a9a520f98f50081fcp-24 + w *
                                        (0x19baa1b1223219.0p-79 + w *
                                            (0x1dc6c702a05262.0p-81 + w *
                                                (0x194c0668da786a.0p-81 + w *
                                                    (0x1a92fc98c29554.0p-82 + w *
                                                        0x147edbdba6f43a.0p-85))))))))))));
    var v: f128 = z *
        (0x1.ba1ba1ba1ba1ba1ba1ba1b694cd6p-5 + w *
            (0x1.226e355e6c23c8f5b4f5762322eep-7 + w *
                (0x1.7da36452b75e2b5fce9ee7c2c92ep-10 + w *
                    (0x1.f57d7734d1656e0aceb716f614c2p-13 + w *
                        (0x1.497d8eea21e95bc7e2aa79b9f2cdp-15 + w *
                            (0x1.b0f72d33eff7bfa2fbc1059d90b6p-18 + w *
                                (0x1.1c77d6eac0234988cdaa04c96626p-20 + w *
                                    (0x1.75c7357d0298c01a31d0a6f7d518p-23 + w *
                                        (0x1e8a7592977938.0p-78 + w *
                                            (0x107385dfb24529.0p-80 + w *
                                                (-0x19ecef3569ebb6.0p-82 + w *
                                                    (-0x12e763b8845268.0p-81 + w *
                                                        -0x151106cbc779a9.0p-83))))))))))));
    var s: f128 = z * xx;
    r = yy + z * (s * (r + v) + yy);
    r += 0x1.5555555555555555555555555553p-2 * s;
    w = xx + r;
    if (i == 1) {
        v = types.scast(f128, iiy);
        return osign *
            (v - 2.0 * (xx - (w * w / (w + v) - r)));
    }

    if (iiy == 1) {
        return w;
    } else {
        // If allow error up to 2 ulp, simply return
        // -1.0 / (xx + r) here
        // Compute -1.0 / (xx + r) accurately
        z = w;
        v = r - (z - xx); // z + v = r + xx
        const a: f128 = -1.0 / w;
        const t: f128 = -1.0 / w;
        s = 1.0 + t * z;
        return t + a * (s + t * v);
    }
}
