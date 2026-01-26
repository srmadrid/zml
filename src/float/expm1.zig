const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Expm1(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.expm1: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

/// Returns the exponential minus one $e^x - 1$ of a float, int or bool operand.
/// The result type is determined by coercing the operand type to a float, and
/// the operation is performed by casting the operand to the result type, then
/// computing its exponential minus one.
///
/// ## Signature
/// ```zig
/// float.expm1(x: X) float.Expm1(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the exponential minus one of.
///
/// ## Returns
/// `float.Expm1(@TypeOf(x))`: The exponential minus one of `x`.
pub inline fn expm1(x: anytype) float.Expm1(@TypeOf(x)) {
    switch (float.Expm1(@TypeOf(x))) {
        f16 => return types.scast(f16, expm1_32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_expm1f.c
            return expm1_32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_expm1.c
            return expm1_64(types.scast(f64, x));
        },
        f80 => {
            //
            // return expm1_80(types.scast(f80, x));
            return types.scast(f80, expm1_128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/s_expm1l.c
            return expm1_128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_expm1f.c
//
// Original copyright notice:
// s_expm1f.c -- float version of s_expm1.c.
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
fn expm1_32(x: f32) f32 {
    var hx: u32 = @bitCast(x);
    const xsb: u32 = hx & 0x80000000;
    hx &= 0x7fffffff;

    // Filter out huge and non-finite argument
    if (hx >= 0x4195b844) { // |x| >= 27 * ln(2)
        if (hx >= 0x42b17218) { // |x| >= 88.721...
            if (hx > 0x7f800000)
                return x + x; // NaN

            if (hx == 0x7f800000)
                return if (xsb == 0) x else -1.0; // exp(±inf) = {inf, -1}

            if (x > 8.8721679688e+1) // Overflow
                return std.math.inf(f32);
        }

        if (xsb != 0) // x < -27 * ln(2), return -1.0
            return -1.0;
    }

    // Argument reduction
    var k: i32 = undefined;
    var xx: f32 = x;
    var c: f32 = undefined;
    if (hx > 0x3eb17218) { // |x| > 0.5 * ln(2)
        var hi: f32 = undefined;
        var lo: f32 = undefined;
        if (hx < 0x3f851592) { // and |x| < 1.5 * ln(2)
            if (xsb == 0) {
                hi = x - 6.9313812256e-1;
                lo = 9.0580006145e-6;
                k = 1;
            } else {
                hi = x + 6.9313812256e-1;
                lo = -9.0580006145e-6;
                k = -1;
            }
        } else {
            k = types.scast(i32, 1.4426950216e+0 * x + @as(f32, if (xsb == 0) 0.5 else -0.5));
            const t: f32 = types.scast(f32, k);
            hi = x - t * 6.9313812256e-1; // t * ln2_hi is exact here
            lo = t * 9.0580006145e-6;
        }

        xx = hi - lo;
        c = (hi - xx) - lo;
    } else if (hx < 0x33000000) { // |x| < 2**-25, return x
        return x;
    } else {
        k = 0;
    }

    // x is now in primary range
    const hfx: f32 = 0.5 * xx;
    const hxs: f32 = xx * hfx;
    const r1: f32 = 1.0 + hxs * (-3.3333212137e-2 + hxs * 1.5807170421e-3);
    var t: f32 = 3.0 - r1 * hfx;
    var e: f32 = hxs * ((r1 - t) / (6.0 - xx * t));
    if (k == 0)
        return xx - (xx * e - hxs); // c is 0

    const twopk: f32 = @bitCast(0x3f800000 +% (k << 23)); // 2**k
    e = (xx * (e - c) - c);
    e -= hxs;
    if (k == -1)
        return 0.5 * (xx - e) - 0.5;

    if (k == 1) {
        if (xx < -0.25)
            return -2.0 * (e - (xx + 0.5))
        else
            return 1.0 + 2.0 * (xx - e);
    }

    if (k <= -2 or k > 56) { // Suffice to return exp(xx) - 1
        var y: f32 = 1.0 - (e - xx);
        if (k == 128)
            y = y * 2.0 * 0x1p127
        else
            y = y * twopk;

        return y - 1.0;
    }

    var y: f32 = undefined;
    if (k < 23) {
        t = @bitCast(0x3f800000 -% (@as(i32, 0x1000000) >> @intCast(k))); // t= 1 - 2**-k
        y = t - (e - xx);
        y = y * twopk;
    } else {
        t = @bitCast(((0x7f -% k) << 23)); // 2**-k
        y = xx - (e + t);
        y += 1.0;
        y = y * twopk;
    }

    return y;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_expm1.c
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
fn expm1_64(x: f64) f64 {
    var hx: u32 = dbl64.getHighPart(x);
    const xsb: u32 = hx & 0x80000000;
    hx &= 0x7fffffff;

    // Filter out huge and non-finite argument
    if (hx >= 0x4043687a) { // |x| >= 56 * ln(2)
        if (hx >= 0x40862e42) { // |x| >= 709.78...
            if (hx >= 0x7ff00000) {
                const low: u32 = dbl64.getLowPart(x);
                if (((hx & 0xfffff) | low) != 0)
                    return x + x // NaN
                else
                    return if (xsb == 0) x else -1.0; // exp(±inf) - 1 = {inf, -1}
            }

            if (x > 7.09782712893383973096e+2) // Overflow
                return std.math.inf(f64);
        }

        if (xsb != 0) // x < -56 * ln(2), return -1.0
            return -1.0;
    }

    // Argument reduction
    var k: i32 = undefined;
    var xx: f64 = x;
    var c: f64 = undefined;
    if (hx > 0x3fd62e42) { // |x| > 0.5 * ln(2)
        var hi: f64 = undefined;
        var lo: f64 = undefined;
        if (hx < 0x3ff0a2b2) { // and |x| < 1.5 * ln(2)
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
            hi = x - t * 6.93147180369123816490e-1; // t * ln2_hi is exact here
            lo = t * 1.90821492927058770002e-10;
        }

        xx = hi - lo;
        c = (hi - xx) - lo;
    } else if (hx < 0x3c900000) { // |x| < 2**-54, return x
        return x;
    } else {
        k = 0;
    }

    // x is now in primary range
    const hfx: f64 = 0.5 * xx;
    const hxs: f64 = xx * hfx;
    const r1: f64 = 1.0 + hxs *
        (-3.33333333333331316428e-2 + hxs *
            (1.58730158725481460165e-3 + hxs *
                (-7.93650757867487942473e-5 + hxs *
                    (4.00821782732936239552e-6 + hxs *
                        -2.01099218183624371326e-7))));
    var t: f64 = 3.0 - r1 * hfx;
    var e: f64 = hxs * ((r1 - t) / (6.0 - xx * t));
    if (k == 0)
        return xx - (xx * e - hxs); // c is 0

    const twopk: f64 = dbl64.Parts.toFloat(.{ .msw = @bitCast(0x3ff00000 +% (k << 20)), .lsw = 0 }); // 2**k
    e = (xx * (e - c) - c);
    e -= hxs;
    if (k == -1)
        return 0.5 * (xx - e) - 0.5;

    if (k == 1) {
        if (xx < -0.25)
            return -2.0 * (e - (xx + 0.5))
        else
            return 1.0 + 2.0 * (xx - e);
    }

    if (k <= -2 or k > 56) { // Suffice to return exp(x) - 1
        var y: f64 = 1.0 - (e - xx);

        if (k == 1024)
            y = y * 2.0 * 0x1p1023
        else
            y = y * twopk;

        return y - 1.0;
    }

    t = 1.0;
    var y: f64 = undefined;
    if (k < 20) {
        dbl64.setHighPart(&t, @bitCast(0x3ff00000 -% (@as(i32, 0x200000) >> @intCast(k)))); // t = 1 - 2**-k
        y = t - (e - xx);
        y = y * twopk;
    } else {
        dbl64.setHighPart(&t, @bitCast((0x3ff -% k) << 20)); // 2**-k
        y = xx - (e + t);
        y += 1.0;
        y = y * twopk;
    }

    return y;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/s_expm1l.c
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
fn expm1_128(x: f128) f128 {
    const u: ldbl128.Parts32 = .fromFloat(x);
    var ix: u32 = u.mswhi;
    const sign: u32 = ix & 0x80000000;
    ix &= 0x7fffffff;

    if (ix >= 0x7fff0000) { // |x| >= inf
        if (((ix & 0xffff) | u.mswlo | u.lswhi | u.lswlo) == 0) {
            if (sign != 0)
                return -1.0
            else
                return x;
        }

        // NaN
        return x;
    }

    if (ix == 0 and (u.mswlo | u.lswhi | u.lswlo) == 0) // expm1(±0) = ±0
        return x;

    if (x > 1.1356523406294143949491931077970764891253e4) // Overflow
        return std.math.inf(f128);

    if (x < -7.9018778583833765273564461846232128760607e1) // Underflow
        return -1.0;

    // Express x = ln(2) (k + remainder), remainder not exceeding 1/2
    var xx: f128 = 6.93145751953125e-1 + 1.428606820309417232121458176568075500134e-6; // ln(2)
    var px: f128 = float.floor(0.5 + x / xx);
    const k: i32 = types.scast(i32, px);

    // Remainder times ln(2)
    var xxx: f128 = x - px * 6.93145751953125e-1;
    xxx -= px * 1.428606820309417232121458176568075500134e-6;

    // Approximate exp(remainder ln(2))
    px = (((((((-4.888737542888633647784737721812546636240e-1 * xxx +
        4.401308817383362136048032038528753151144e1) * xxx +
        -1.716772506388927649032068540558788106762e3) * xxx +
        4.578962475841642634225390068461943438441e4) * xxx +
        -7.212432713558031519943281748462837065308e5) * xxx +
        8.944630806357575461578107295909719817253e6) * xxx +
        -5.722847283900608941516165725053359168840e7) * xxx +
        2.943520915569954073888921213330863757240e8) * xxx;

    var qx: f128 = (((((((xxx +
        -8.802340681794263968892934703309274564037e1) * xxx +
        3.697714952261803935521187272204485251835e3) * xxx +
        -9.615511549171441430850103489315371768998e4) * xxx +
        1.682912729190313538934190635536631941751e6) * xxx +
        -2.019684072836541751428967854947019415698e7) * xxx +
        1.615869009634292424463780387327037251069e8) * xxx +
        -7.848989743695296475743081255027098295771e8) * xxx +
        1.766112549341972444333352727998584753865e9;

    xx = xxx * xxx;
    qx = xxx + (0.5 * xx + xx * px / qx);

    // exp(x) = exp(k * ln(2)) * exp(remainder ln(2)) = 2**k * exp(remainder ln(2))
    // We have qx = exp(remainder ln(2)) - 1, so
    //  exp(x) - 1 = 2**k * (qx + 1) - 1 = 2**k * qx + 2**k - 1
    if (k == 16384) {
        // If k = 16384, 2**k is Inf, but the result is finite because qx < 0.
        // We compute: 2 * (2**(k - 1) * qx + 2**(k - 1)) - 1
        // This avoids intermediate infinity.
        px = std.math.ldexp(@as(f128, 1.0), k - 1);
        xxx = px * qx + (px - 0.5);
        xxx *= 2.0;
    } else {
        px = float.ldexp(@as(f128, 1.0), k);
        xxx = px * qx + (px - 1.0);
    }

    return xxx;
}
