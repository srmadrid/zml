const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Frexp(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.frexp: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

pub inline fn frexp(x: anytype, e: *i32) Frexp(@TypeOf(x)) {
    switch (Frexp(@TypeOf(x))) {
        f16 => return types.scast(f16, frexp32(types.scast(f32, x), e)),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_frexpf.c
            return frexp32(types.scast(f32, x), e);
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_frexp.c
            return frexp64(types.scast(f64, x), e);
        },
        f80 => {
            //
            // return frexp80(types.types.scast(f80, x));
            return types.scast(f80, frexp128(types.scast(f128, x), e));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_frexpl.c
            return frexp128(types.scast(f128, x), e);
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_frexpf.c
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
fn frexp32(x: f32, e: *i32) f32 {
    var hx: i32 = @bitCast(x);
    var ix: i32 = 0x7fffffff & hx;
    e.* = 0;

    if (ix >= 0x7f800000 || (ix == 0))
        return x + x; // 0, inf, nan

    var xx: f32 = x;
    if (ix < 0x00800000) { // subnormal
        xx *= 3.3554432000e+7;
        hx = @bitCast(xx);
        ix = hx & 0x7fffffff;
        e.* = -25;
    }

    e.* += (ix >> 23) - 126;
    hx = (hx & 0x807fffff) | 0x3f000000;
    xx = @bitCast(hx);
    return xx;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_frexp.c
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
fn frexp64(x: f64, e: *i32) f64 {
    var hx: i32 = @bitCast(dbl64.getHighPart(x));
    const lx: i32 = @bitCast(dbl64.getLowPart(x));
    var ix: i32 = 0x7fffffff & hx;
    e.* = 0;

    if (ix >= 0x7ff00000 or ((ix | lx) == 0)) // 0, inf, nan
        return x;

    var xx: f64 = x;
    if (ix < 0x00100000) { // subnormal
        x *= 1.80143985094819840000e+16;
        hx = @bitCast(dbl64.getHighPart(x));
        ix = hx & 0x7fffffff;
        e.* = -54;
    }

    e.* += (ix >> 20) - 1022;
    hx = (hx & 0x800fffff) | 0x3fe00000;
    dbl64.setHighPart(&xx, @bitCast(hx));
    return xx;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_frexpl.c
//
// Original copyright notice:
// Copyright (c) 2004-2005 David Schultz <das@FreeBSD.ORG>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
fn frexp128(x: f128, e: *i32) f128 {
    var u: ldbl128.ShapeSplit = .fromFloat(x);
    switch (u.exponent) {
        0 => { // 0 or subnormal
            if ((u.mantissa_low | u.mantissa_high) == 0) {
                e.* = 0;
            } else {
                u = .fromFloat(u.toFloat() * 0x1.0p514);
                e.* = @as(i32, @intCast(u.exponent)) - 0x4200;
                u.exponent = 0x3ffe;
            }
        },
        0x7fff => {}, // infinity or NaN; value of e is unspecified
        else => { // normal
            e.* = @as(i32, @intCast(u.exponent)) - 0x3ffe;
            u.exponent = 0x3ffe;
        },
    }

    return u.toFloat();
}
