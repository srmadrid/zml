const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");

pub inline fn exp(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.exp: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return types.scast(f16, exp32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_expf.c
            return exp32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_exp.c
            return exp64(types.scast(f64, x));
        },
        f80 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld80/e_expl.c
            return exp80(types.scast(f80, x));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_expl.c
            return exp128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_expf.c
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
fn exp32(x: f32) f32 {
    var hx: u32 = @bitCast(x);
    const xsb: u32 = (hx >> 31) & 1;
    hx &= 0x7fffffff;

    // Filter out non-finite argument
    if (hx >= 0x42b17218) { // |x| >= 88.721...
        if (hx > 0x7f800000) // NaN
            return x + x;

        if (hx == 0x7f800000) // ±inf
            return if (xsb == 0) std.math.inf(f32) else 0.0;

        if (x > 8.8721679688e+1) // overflow
            return std.math.inf(f32);

        if (x < -1.0397208405e+2) // underflow
            return 0.0;
    }

    // Argument reduction
    var k: i32 = 0;
    var xx: f32 = x;
    var hi: f32 = undefined;
    var lo: f32 = undefined;
    if (hx > 0x3eb17218) { // |x| > 0.5 * ln(2)
        if (hx < 0x3f851592) { // |x| < 1.5 * ln(2)
            if (xsb == 0) {
                hi = x - 6.9314575195e-1;
                lo = 1.4286067653e-6;
                k = 1;
            } else {
                hi = x + 6.9314575195e-1;
                lo = -1.4286067653e-6;
                k = -1;
            }
        } else {
            k = types.scast(i32, 1.4426950216e+0 * x + @as(f32, if (xsb == 0) 0.5 else -0.5));
            const t: f32 = types.scast(f32, k);
            hi = x - t * 6.9314575195e-1;
            lo = t * 1.4286067653e-6;
        }

        xx = hi - lo;
    } else if (hx < 0x39000000) { // |x| < 2**-14
        return 1.0 + x;
    } else {
        k = 0;
    }

    // xx is now in primary range
    const t: f32 = xx * xx;
    var twopk: f32 = undefined;
    if (k >= -125)
        twopk = @bitCast(0x3f800000 +% (k << 23))
    else
        twopk = @bitCast(0x3f800000 +% ((k + 100) << 23));

    const c: f32 = xx -
        t * (1.6666625440e-1 +
            t * -2.7667332906e-3);

    var y: f32 = undefined;
    if (k == 0)
        return 1.0 - ((xx * c) / (c - 2.0) - xx)
    else
        y = 1.0 - ((lo - (xx * c) / (2.0 - c)) - hi);

    if (k >= -125) {
        if (k == 128)
            return y * 2.0 * 0x1p127;

        return y * twopk;
    } else {
        return y * twopk * 7.8886090522e-31;
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_exp.c
//
// Original copyright notice:
// ====================================================
// Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
//
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn exp64(x: f64) f64 {
    var hx: u32 = dbl64.getHighPart(x);
    const xsb: u32 = (hx >> 31) & 1;
    hx &= 0x7fffffff;

    // Filter out non-finite argument
    if (hx >= 0x40862e42) { // |x| >= 709.78...
        if (hx > 0x7ff00000) {
            const lx: u32 = dbl64.getLowPart(x);
            if ((hx & 0xfffff) | lx != 0) // NaN
                return x + x
            else // ±inf
                return if (xsb == 0) std.math.inf(f64) else 0.0;
        }

        if (x > 7.09782712893383973096e+2) // overflow
            return std.math.inf(f64);

        if (x < -7.45133219101941108420e+2) // underflow
            return 0.0;
    }

    if (x == 1.0)
        return float.e(f64);

    // Argument reduction
    var k: i32 = 0;
    var xx: f64 = x;
    var hi: f64 = undefined;
    var lo: f64 = undefined;
    if (hx > 0x3fd62e42) { // |x| > 0.5 * ln(2)
        if (hx < 0x3ff0a2b2) { // |x| < 1.5 * ln(2)
            if (xsb == 0) {
                hi = x - 6.93147180369123816490e-1;
                lo = 1.90821492927058770002e-10;
                k = 1;
            } else {
                hi = x + 6.93147180369123816490e-1;
                lo = -1.90821492927058770002e-10;
                k = -1;
            }
        } else {
            k = types.scast(i32, 1.44269504088896338700e+0 * x + @as(f64, if (xsb == 0) 0.5 else -0.5));
            const t: f64 = types.scast(f64, k);
            hi = x - t * 6.93147180369123816490e-1;
            lo = t * 1.90821492927058770002e-10;
        }

        xx = hi - lo;
    } else if (hx < 0x39000000) { // |x| < 2**-14
        return 1.0 + x;
    } else {
        k = 0;
    }

    // xx is now in primary range
    const t: f64 = xx * xx;
    var twopk: f64 = undefined;
    if (k >= -1021)
        dbl64.setFromParts(&twopk, @bitCast(0x3ff00000 +% (k << 20)), 0)
    else
        dbl64.setFromParts(&twopk, @bitCast(0x3ff00000 +% ((k + 1000) << 20)), 0);

    const c: f64 = xx -
        t * (1.66666666666666019037e-1 +
            t * (-2.77777777770155933842e-3 +
                t * (6.61375632143793436117e-5 +
                    t * (-1.65339022054652515390e-6 +
                        t * 4.13813679705723846039e-8))));

    var y: f64 = undefined;
    if (k == 0)
        return 1.0 - ((xx * c) / (c - 2.0) - xx)
    else
        y = 1.0 - ((lo - (xx * c) / (2.0 - c)) - hi);

    if (k >= -1021) {
        if (k == 1024)
            return y * 2.0 * 0x1p1023;

        return y * twopk;
    } else {
        return y * twopk * 9.33263618503218878990e-302;
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld80/e_expl.c
//
// Original copyright notice:
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
fn exp80(x: f80) f80 {
    if (std.math.isNan(x))
        return x;

    if (x > 1.1356523406294143949492e4)
        return std.math.inf(f80);

    if (x < -1.13994985314888605586758e4)
        return 0.0;

    var px: f80 = float.floor(1.4426950408889634073599e0 * x + 0.5);
    const n: i32 = types.scast(i32, px);
    var xx: f80 = x - px * 6.93145751953125e-1;
    xx -= px * 1.4286068203094172321215e-6;

    const xx2: f80 = xx * xx;
    px = xx *
        (9.9999999999999999991025e-1 +
            xx2 * (3.0299440770744196129956e-2 +
                xx2 * 1.2617719307481059087798e-4));
    xx = px /
        ((2.0000000000000000000897e0 +
            xx2 * (2.2726554820815502876593e-1 +
                xx2 * (2.5244834034968410419224e-3 +
                    xx2 * 3.0019850513866445504159e-6))) - px);

    xx = 1.0 + xx + xx;

    return float.ldexp(xx, n);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_expl.c
//
// Original copyright notice:
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
fn exp128(x: f128) f128 {
    if (std.math.isNan(x))
        return x;

    if (x > 1.1356523406294143949491931077970764891253e4)
        return std.math.inf(f128);

    if (x < -1.143276959615573793352782661133116431383730e4)
        return 0.0;

    var px: f128 = float.floor(1.442695040888963407359924681001892137426646 * x + 0.5);
    const n: i32 = types.scast(i32, px);
    var xx: f128 = x + px * -6.93145751953125e-1;
    xx += px * -1.428606820309417232121458176568075500134e-6;

    const xx2: f128 = xx * xx;
    px = xx *
        (9.999999999999999999999999999999999998502e-1 +
            xx2 * (3.508710990737834361215404761139478627390e-2 +
                xx2 * (2.708775201978218837374512615596512792224e-4 +
                    xx2 * (6.141506007208645008909088812338454698548e-7 +
                        xx2 * 3.279723985560247033712687707263393506266e-10))));
    xx = px /
        ((2.000000000000000000000000000000000000150e0 +
            xx2 * (2.368408864814233538909747618894558968880e-1 +
                xx2 * (3.611828913847589925056132680618007270344e-3 +
                    xx2 * (1.504792651814944826817779302637284053660e-5 +
                        xx2 * (1.771372078166251484503904874657985291164e-8 +
                            xx2 * 2.980756652081995192255342779918052538681e-12))))) - px);

    xx = 1.0 + xx + xx;

    return float.ldexp(xx, n);
}
