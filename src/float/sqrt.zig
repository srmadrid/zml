const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn sqrt(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.sqrt: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, sqrt32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_sqrtf.c
            return sqrt32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_sqrt.c
            return sqrt64(types.scast(f64, x));
        },
        f80 => {
            //
            // return sqrt80(types.scast(f80, x));
            return types.scast(f80, sqrt128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_sqrtl.c
            return sqrt128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_sqrtf.c
//
// Original copyright notice:
// e_powf.c -- float version of e_pow.c.
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
fn sqrt32(x: f32) f32 {
    var ix: i32 = @bitCast(x);
    const sign: u32 = 0x80000000;

    // Take care of Inf and NaN
    if ((ix & 0x7f800000) == 0x7f800000) {
        return x * x + x; // sqrt(NaN) = NaN, sqrt(+inf) = +inf, sqrt(-inf) = sNaN
    }

    // Take care of zero
    if (ix <= 0) {
        if ((ix & (~sign)) == 0)
            return x // sqrt(±0) = ±0
        else if (ix < 0)
            return std.math.snan(f32); // sqrt(-ve) = sNaN

    }

    // Normalize x
    var m: i32 = (ix >> 23);
    if (m == 0) { // Subnormal x
        var i: i32 = 0;
        while ((ix & 0x00800000) == 0) : (i += 1)
            ix <<= 1;

        m -%= i -% 1;
    }

    m -%= 127; // Unbias exponent
    ix = (ix & 0x007fffff) | 0x00800000;

    if ((m & 1) != 0) // Odd m, double x to make it even
        ix +%= ix;

    m >>= 1; // m = [m/2]

    // Generate sqrt(x) bit by bit
    ix +%= ix;
    var s: i32 = 0;
    var q: i32 = 0; // q = sqrt(x)
    var r: i32 = 0x01000000; // r = moving bit from right to left
    while (r != 0) {
        const t: i32 = s +% r;
        if (t <= ix) {
            s = t +% r;
            ix -%= t;
            q +%= r;
        }

        ix +%= ix;
        r >>= 1;
    }

    ix = (q >> 1) +% 0x3f000000;
    ix +%= (m << 23);
    return @bitCast(ix);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_sqrt.c
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
fn sqrt64(x: f64) f64 {
    var ix0: i32 = @bitCast(dbl64.getHighPart(x));
    var ix1: u32 = dbl64.getLowPart(x);
    const sign: u32 = 0x80000000;

    // Take care of Inf and NaN
    if ((ix0 & 0x7ff00000) == 0x7ff00000) {
        return x * x + x; // sqrt(NaN) = NaN, sqrt(+inf) = +inf, sqrt(-inf) = sNaN
    }

    // Take care of zero
    if (ix0 <= 0) {
        if (((ix0 & (~sign)) | @as(i32, @bitCast(ix1))) == 0)
            return x // sqrt(±0) = ±0
        else if (ix0 < 0)
            return std.math.snan(f64); // sqrt(-ve) = sNaN
    }

    // Normalize x
    var m: i32 = (ix0 >> 20);
    if (m == 0) { // Subnormal x
        while (ix0 == 0) {
            m -%= 21;
            ix0 |= @bitCast(ix1 >> 11);
            ix1 <<= 21;
        }

        var i: i32 = 0;
        while ((ix0 & 0x00100000) == 0) : (i += 1)
            ix0 <<= 1;

        m -%= i - 1;
        ix0 |= @bitCast(ix1 >> @intCast(32 - i));
        ix1 <<= @intCast(i);
    }

    m -%= 1023; // Unbias exponent
    ix0 = (ix0 & 0x000fffff) | 0x00100000;
    if ((m & 1) != 0) { // Odd m, double x to make it even
        ix0 +%= ix0 +% types.scast(i32, (ix1 & sign) >> 31);
        ix1 +%= ix1;
    }

    m >>= 1; // m = [m/2]

    // Generate sqrt(x) bit by bit
    ix0 +%= ix0 +% types.scast(i32, (ix1 & sign) >> 31);
    ix1 +%= ix1;
    var q: i32 = 0;
    var s0: i32 = 0;
    var r: u32 = 0x00200000; // r = moving bit from right to left
    while (r != 0) {
        const t: i32 = s0 +% types.scast(i32, r);
        if (t <= ix0) {
            s0 = t +% types.scast(i32, r);
            ix0 -%= t;
            q +%= types.scast(i32, r);
        }

        ix0 +%= ix0 +% types.scast(i32, (ix1 & sign) >> 31);
        ix1 +%= ix1;
        r >>= 1;
    }

    var q1: u32 = 0; // [q, q1] = sqrt(x)
    var s1: u32 = 0;
    r = sign;
    while (r != 0) {
        const t1: u32 = s1 +% r;
        const t: i32 = s0;
        if ((t < ix0) or (t == ix0 and t1 <= ix1)) {
            s1 = t1 +% r;
            if (((t1 & sign) == sign) and (s1 & sign) == 0)
                s0 +%= 1;

            ix0 -%= t;
            if (ix1 < t1)
                ix0 -%= 1;

            ix1 -%= t1;
            q1 +%= r;
        }

        ix0 +%= ix0 +% types.scast(i32, (ix1 & sign) >> 31);
        ix1 +%= ix1;
        r >>= 1;
    }

    ix0 = (q >> 1) +% 0x3fe00000;
    ix1 = q1 >> 1;
    if ((q & 1) == 1)
        ix1 |= sign;

    ix0 +%= (m << 20);
    return dbl64.Parts.toFloat(.{ .msw = @bitCast(ix0), .lsw = ix1 });
}

fn sqrt80(x: f80) f80 {
    _ = x;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_sqrtl.c
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
fn sqrt128(x: f128) f128 {
    var u: ldbl128.ShapeSplit = .fromFloat(x);

    // If x = NaN, sqrt(x) = NaN
    // If x = +inf, sqrt(x) = +inf
    // If x = -inf, sqrt(x) = NaN
    if (u.exponent == 16384 * 2 - 1)
        return x * x + x;

    // If x = ±0, sqrt(x) = ±0
    if ((u.mantissa_high | u.mantissa_low | u.exponent) == 0)
        return x;

    // If x < 0, sqrt(x) = sNaN
    if (u.sign != 0)
        return std.math.snan(f128);

    var k: i32 = 0;
    if (u.exponent == 0) {
        // Adjust subnormal numbers
        u = .fromFloat(u.toFloat() * 0x1.0p514);
        k = -514;
    }

    // u is a normal number, so break it into u = e * 2^n. u = (2 * e) * 2^(2k) for odd n and u = (4 * e) * 2^(2k) for even n
    if ((u.exponent -% 0x3ffe) & 1 != 0) { // n is odd
        k +%= @as(i32, @intCast(u.exponent)) -% 0x3fff; // 2k = n - 1
        u.exponent = 0x3fff; // u in [1,2)
    } else {
        k +%= @as(i32, @intCast(u.exponent)) -% 0x4000; // 2k = n - 2
        u.exponent = 0x4000; // u in [2,4)
    }

    // Newton's iteration. Split u into a high and low part to achieve additional precision
    var xn: f128 = types.scast(f128, sqrt64(types.scast(f64, u.toFloat())));
    xn = (xn + (u.toFloat() / xn)) * 0.5;
    var lo: f128 = u.toFloat();
    u.mantissa_low = 0; // Zero out lower bits
    lo = (lo - u.toFloat()) / xn; // Low bits divided by xn
    xn += u.toFloat() / xn; // High portion estimate
    u = .fromFloat(xn + lo); // Combine everything
    if ((k >> 1) -% 1 > 0)
        u.exponent +%= @intCast((k >> 1) -% 1)
    else
        u.exponent -%= @intCast(-((k >> 1) -% 1));
    xn = x / u.toFloat();

    if (xn == u.toFloat()) {
        return u.toFloat();
    }

    u = .fromFloat(u.toFloat() + xn);
    u.exponent -%= 1;
    return u.toFloat();
}
