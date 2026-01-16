const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn rint(x: anytype) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!types.isNumeric(X) or types.numericType(X) != .float)
        @compileError("zml.float.rint: x must be a float, got \n\tx: " ++ @typeName(@TypeOf(x)) ++ "\n");

    switch (@TypeOf(x)) {
        f16 => return types.cast(f16, rint32(types.cast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_rintf.c
            return rint32(x);
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_rint.c
            return rint64(x);
        },
        f80 => return types.cast(f80, rint128(types.cast(f128, x, .{})), .{}),
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_rintl.c
            return rint128(x);
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_rintf.c
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
fn rint32(x: f32) f32 {
    const TWO23: [2]f32 = .{
        8.3886080000e+6,
        -8.3886080000e+6,
    };

    var I0: i32 = @bitCast(x);
    const sx: i32 = (I0 >> 31) & 1;
    const j0: i32 = ((I0 >> 23) & 0xff) -% 0x7f;
    if (j0 < 23) {
        if (j0 < 0) {
            if ((I0 & 0x7fffffff) == 0)
                return x;

            const w: f32 = TWO23[types.scast(u32, sx)] + x;
            var t = w - TWO23[types.scast(u32, sx)];
            I0 = @bitCast(t);
            t = @bitCast((I0 & 0x7fffffff) | (sx << 31));
            return t;
        }

        const w: f32 = TWO23[types.scast(u32, sx)] + x;
        return w - TWO23[types.scast(u32, sx)];
    }

    if (j0 == 0x80)
        return x + x // inf or NaN
    else
        return x; // x is integral
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_rint.c
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
fn rint64(x: f64) f64 {
    const TWO52: [2]f64 = .{
        4.50359962737049600000e+15,
        -4.50359962737049600000e+15,
    };

    var I0: i32 = @bitCast(dbl64.getHighPart(x));
    var I1: i32 = @bitCast(dbl64.getLowPart(x));
    const sx: i32 = (I0 >> 31) & 1;
    const j0: i32 = ((I0 >> 20) & 0x7ff) - 0x3ff;
    if (j0 < 20) {
        if (j0 < 0) {
            if (((I0 & 0x7fffffff) | I1) == 0)
                return x;

            I1 |= (I0 & 0x0fffff);
            I0 &= 0xfffe0000;
            I0 |= ((I1 | -I1) >> 12) & 0x80000;
            var xx: f64 = x;
            dbl64.setHighPart(&xx, I0);
            const w: f64 = TWO52[types.scast(u32, sx)] + xx;
            var t: f64 = w - TWO52[types.scast(u32, sx)];
            I0 = @bitCast(dbl64.getHighPart(t));
            dbl64.setHighPart(&t, (I0 & 0x7fffffff) | (sx << 31));
            return t;
        } else {
            var i: i32 = (0x000fffff) >> @intCast(j0);
            if (((I0 & i) | I1) == 0)
                return x; // x is integral

            i >>= 1;
            if (((I0 & i) | I1) != 0) {
                if (j0 == 19) I1 = 0x40000000 else if (j0 == 18)
                    I1 = 0x80000000
                else
                    I0 = (I0 & (~i)) | ((0x20000) >> @intCast(j0));
            }
        }
    } else if (j0 > 51) {
        if (j0 == 0x400) // inf or NaN
            return x + x
        else // x is integral
            return x;
    } else {
        var i: i32 = @bitCast(@as(u32, 0xffffffff) >> @intCast(j0 -% 20));
        if ((I1 & i) == 0)
            return x; // x is integral

        i >>= 1;
        if ((I1 & i) != 0)
            I1 = (I1 & (~i)) | ((0x40000000) >> @intCast(j0 - 20));
    }

    var xx: f64 = x;
    dbl64.setFromParts(&xx, @bitCast(I0), @bitCast(I1));
    const w: f64 = TWO52[types.scast(u32, sx)] + xx;
    return w - TWO52[types.scast(u32, sx)];
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_rint.c
//
// Original copyright notice:
// Copyright (c) 2008 David Schultz <das@FreeBSD.ORG>
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
fn rint128(x: f128) f128 {
    const u: ldbl128.ShapeSplit = .fromFloat(x);
    const ex: i32 = @intCast(u.exponent);

    if (ex >= (16384 - 1) + 113 - 1) {
        if (ex == (16384 - 1) + 16384)
            return x + x; // Inf, NaN, or unsupported format

        return x; // Finite and already an integer
    }

    var xx: f128 = x;
    xx += if (u.sign == 0) 0x1.0p112 else -0x1.0p112;
    xx -= if (u.sign == 0) 0x1.0p112 else -0x1.0p112;

    if (ex < (16384 - 1) and xx == 0.0)
        return if (u.sign == 0) 0.0 else -0.0;

    return xx;
}
