const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const Coerce = types.Coerce;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn pow(x: anytype, y: anytype) EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.pow: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    comptime if (types.numericType(@TypeOf(y)) != .int and types.numericType(@TypeOf(y)) != .float)
        @compileError("float.pow: y must be an int or float, got " ++ @typeName(@TypeOf(y)));

    switch (EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y)))) {
        f16 => return types.scast(f16, pow32(types.scast(f32, x), types.scast(f32, y))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_powf.c
            return pow32(types.scast(f32, x), types.scast(f32, y));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_pow.c
            return pow64(types.scast(f64, x), types.scast(f64, y));
        },
        f80 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld80/e_powl.c
            // return pow80(types.scast(f80, x), types.scast(f80, y));
            // The f80 implementation is broken, fall back to f128
            return types.scast(f80, pow128(types.scast(f128, x), types.scast(f128, y)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_logl.c
            return pow128(types.scast(f128, x), types.scast(f128, y));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_logf.c
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
fn pow32(x: f32, y: f32) f32 {
    const hx: i32 = @bitCast(x);
    const hy: i32 = @bitCast(y);
    var ix: i32 = hx & 0x7fffffff;
    const iy: i32 = hy & 0x7fffffff;

    if (iy == 0) // x**0 = 1
        return 1.0;

    if (hx == 0x3f800000) // 1**y = 1, even if y is NaN
        return 1.0;

    if (ix > 0x7f800000 or iy > 0x7f800000) // x or y is NaN, with y != 0
        return (x + 0.0) + (y + 0.0);

    // Determine if y is an odd int when x < 0
    // yisint = 0 -> y is not an integer
    // yisint = 1 -> y is an odd integer
    // yisint = 2 -> y is an even integer
    var yisint: i32 = 0;
    if (hx < 0) {
        if (iy >= 0x4b800000) { // Even integer y
            yisint = 2;
        } else if (iy >= 0x3f800000) {
            const k: i32 = (iy >> 23) -% 0x7f;
            const j: i32 = iy >> @intCast(23 -% k);
            if ((j << @intCast(23 -% k)) == iy)
                yisint = 2 -% (j & 1);
        }
    }

    if (iy == 0x7f800000) { // y is ±inf
        if (ix == 0x3f800000) // (-1)**±inf
            return 1.0
        else if (ix > 0x3f800000) // (|x| > 1)**±inf
            return if (hy >= 0) std.math.inf(f32) else 0.0
        else // (|x| < 1)**±inf
            return if (hy >= 0) 0.0 else std.math.inf(f32);
    }

    if (iy == 0x3f800000) // y is ±1
        return if (hy >= 0) x else 1.0 / x;

    if (hy == 0x40000000) // y is 2
        return x * x;

    if (hy == 0x3f000000) // y is 0.5
        if (hx >= 0) // x >= 0
            return float.sqrt(x);

    var ax: f32 = float.abs(x);
    if (ix == 0x7f800000 or ix == 0 or ix == 0x3f800000) { // x is ±inf, ±0, or ±1
        var z: f32 = ax;

        if (hy < 0)
            z = 1.0 / z;

        if (hx < 0) {
            if (((ix -% 0x3f800000) | yisint) == 0) {
                z = std.math.nan(f32); // (-1)**non-int = NaN
            } else if (yisint == 1) {
                z = -z; // (x < 0)**odd = -(|x|**odd)
            }
        }

        return z;
    }

    var n: i32 = types.scast(i32, @as(u32, @bitCast(hx)) >> 31) - 1;

    // (x < 0)**non-int = NaN
    if ((n | yisint) == 0)
        return std.math.nan(f32);

    var sn: f32 = 1.0; // Sign of result
    if ((n | (yisint -% 1)) == 0)
        sn = -1.0;

    var t1: f32 = undefined;
    var t2: f32 = undefined;
    if (iy > 0x4d000000) { // |y| > 2**27
        // over/underflow if x is not close to 1
        if (ix < 0x3f7ffff6)
            return sn * (if (hy < 0) std.math.inf(f32) else 0.0);

        if (ix > 0x3f800007)
            return sn * (if (hy > 0) std.math.inf(f32) else 0.0);

        // Now |1 - x| <= 2**-20, compute log(x) directly
        const t: f32 = ax - 1.0;
        const w: f32 = (t * t) * (0.5 - t * (1.0 / 3.0 - t * 0.25));
        const u: f32 = 1.4426879883e+0 * t; // 1/log(2) * t
        const v: f32 = t * 7.0526075433e-6 - w * 1.4426950216e+0;
        t1 = u + v;
        const is: u32 = @bitCast(t1);
        t1 = @bitCast(is & 0xfffff000);
        t2 = v - (t1 - u);
    } else {
        n = 0;

        // Take care of subnormal numbers
        if (ix < 0x00800000) {
            ax *= 16777216.0; // 2**24
            n -%= 24;
            ix = @bitCast(ax);
        }

        n +%= (ix >> 23) -% 0x7f;
        const j: i32 = ix & 0x007fffff;

        // Determine interval for ax
        ix = j | 0x3f800000; // Normalize ix
        var k: i32 = undefined;
        if (j <= 0x1cc471) // |x| < sqrt(3/2)
            k = 0
        else if (j < 0x5db3d7) // |x| < sqrt(3)
            k = 1
        else {
            k = 0;
            n +%= 1;
            ix -%= 0x00800000;
        }
        ax = @bitCast(ix);

        // Compute s = s_h + s_l = (x - 1)/(x + 1) or (x - 1.5)/(x + 1.5)
        const bp: f32 = if (k == 0) 1.0 else 1.5;
        var u: f32 = ax - bp;
        var v: f32 = 1.0 / (ax + bp);
        const s: f32 = u * v;
        var s_h: f32 = s;
        var is: i32 = @bitCast(s_h);
        s_h = @bitCast(@as(u32, @bitCast(is)) & 0xfffff000);

        // t_h = ax + bp
        is = @bitCast(((@as(u32, @bitCast(ix)) >> 1) & 0xfffff000) | 0x20000000);
        var t_h: f32 = @bitCast(is +% 0x00400000 +% (k << 21));
        var t_l: f32 = ax - (t_h - bp);
        const s_l: f32 = v * ((u - s_h * t_h) - s_h * t_l);

        // Compute log(ax)
        var s2: f32 = s * s;
        var r: f32 = s2 * s2 *
            (6.0000002384e-1 + s2 *
                (4.2857143283e-1 + s2 *
                    (3.3333334327e-1 + s2 *
                        (2.7272811532e-1 + s2 *
                            (2.3066075146e-1 + s2 *
                                2.0697501302e-1)))));
        r += s_l * (s_h + s);
        s2 = s_h * s_h;
        t_h = 3.0 + s2 + r;
        is = @bitCast(t_h);
        t_h = @bitCast(@as(u32, @bitCast(is)) & 0xfffff000);
        t_l = r - ((t_h - 3.0) - s2);

        // u + v = s * (1 + ...)
        u = s_h * t_h;
        v = s_l * t_h + t_l * s;

        // 2/(3 log(2)) * (s + ...)
        var p_h: f32 = u + v;
        is = @bitCast(p_h);
        p_h = @bitCast(@as(u32, @bitCast(is)) & 0xfffff000);
        const p_l: f32 = v - (p_h - u);
        const z_h: f32 = 9.6191406250e-1 * p_h;
        const z_l: f32 = -1.1736857402e-4 * p_h +
            p_l * 9.6179670095e-1 + @as(f32, if (k == 0) 0.0 else 1.56322085e-6);

        // log2(ax) = (s + ...) * 2/(3 log(2))
        const t: f32 = types.scast(f32, n);
        t1 = (((z_h + z_l) + @as(f32, if (k == 0) 0.0 else 5.84960938e-1)) + t);
        is = @bitCast(t1);
        t1 = @bitCast(@as(u32, @bitCast(is)) & 0xfffff000);
        t2 = z_l - (((t1 - t) - @as(f32, if (k == 0) 0.0 else 5.84960938e-1)) - z_h);
    }

    var is: i32 = @bitCast(y);
    const y1: f32 = @bitCast(@as(u32, @bitCast(is)) & 0xfffff000);
    const p_l: f32 = (y - y1) * t1 + y * t2;
    var p_h: f32 = y1 * t1;
    var z: f32 = p_l + p_h;
    var j: i32 = @bitCast(z);
    if (j > 0x43000000) // z > 128
        return sn * std.math.inf(f32)
    else if (j == 0x43000000) { // z == 128
        if (p_l + 4.2995665694e-8 > z - p_h)
            return sn * std.math.inf(f32);
    } else if ((j & 0x7fffffff) > 0x43160000) // z <= -150
        return sn * 0.0
    else if (j == 0xc3160000) { // z == -150
        if (p_l <= z - p_h)
            return sn * 0.0;
    }

    // Compute 2**(p_h + p_l)
    const i: i32 = j & 0x7fffffff;
    var k: i32 = (i >> 23) -% 0x7f;
    n = 0;
    if (i > 0x3f000000) { // |z| > 0.5, set n = [z + 0.5]
        n = j +% (@as(i32, 0x00800000) >> @intCast(k + 1));
        k = ((n & 0x7fffffff) >> 23) -% 0x7f;
        const t: f32 = @bitCast(n & ~(@as(i32, 0x007fffff) >> @intCast(k)));
        n = ((n & 0x007fffff) | 0x00800000) >> @intCast(23 -% k);

        if (j < 0)
            n = -n;

        p_h -= t;
    }

    var t: f32 = p_l + p_h;
    is = @bitCast(t);
    t = @bitCast(@as(u32, @bitCast(is)) & 0xffff8000);
    const u: f32 = t * 6.93145752e-1;
    const v: f32 = (p_l - (t - p_h)) * 6.9314718246e-1 + t * 1.42860654e-6;
    z = u + v;
    const w: f32 = v - (z - u);
    t = z * z;
    t1 = z - t *
        (1.6666667163e-1 + t *
            (-2.7777778450e-3 + t *
                (6.6137559770e-5 + t *
                    (-1.6533901999e-6 + t *
                        4.1381369442e-8))));
    const r: f32 = (z * t1) / (t1 - 2.0) - (w + z * w);
    z = 1.0 - (r - z);
    j = @bitCast(z);
    j +%= n << 23;

    if ((j >> 23) <= 0) // Subnormal output
        z = float.scalbn(z, n)
    else
        z = @bitCast(j);

    return sn * z;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_pow.c
//
// Original copyright notice:
// ====================================================
// Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
//
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn pow64(x: f64, y: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const lx: u32 = dbl64.getLowPart(x);
    const hy: i32 = @bitCast(dbl64.getHighPart(y));
    const ly: i32 = @bitCast(dbl64.getLowPart(y));
    var ix: i32 = hx & 0x7fffffff;
    const iy: i32 = hy & 0x7fffffff;

    if ((iy | @as(i32, @bitCast(ly))) == 0) // x**0 = 1
        return 1.0;

    if (hx == 0x3ff00000 and lx == 0) // 1**y = 1, even if y is NaN
        return 1.0;

    if (ix > 0x7ff00000 or (ix == 0x7ff00000 and lx != 0) or
        iy > 0x7ff00000 or (iy == 0x7ff00000 and ly != 0)) // x or y is NaN, with y != 0
        return (x + 0.0) + (y + 0.0);

    // Determine if y is an odd int when x < 0
    // yisint = 0 -> y is not an integer
    // yisint = 1 -> y is an odd integer
    // yisint = 2 -> y is an even integer
    var yisint: i32 = 0;
    if (hx < 0) {
        if (iy >= 0x43400000) { // Even integer y
            yisint = 2;
        } else if (iy >= 0x3ff00000) {
            const k: i32 = (iy >> 20) -% 0x3ff;
            if (k > 20) {
                const j: i32 = types.scast(i32, ly >> @intCast(52 -% k));
                if ((j << @intCast(52 -% k)) == ly)
                    yisint = 2 - (j & 1);
            } else if (ly == 0) {
                const j: i32 = iy >> @intCast(20 -% k);
                if ((j << @intCast(20 -% k)) == iy)
                    yisint = 2 - (j & 1);
            }
        }
    }

    if (ly == 0) {
        if (iy == 0x7ff00000) { // y is ±inf
            if (((ix -% 0x3ff00000) | @as(i32, @bitCast(lx))) == 0) // (-1)**±inf
                return 1.0
            else if (ix >= 0x3ff00000) // (|x| > 1)**±inf
                return if (hy >= 0) std.math.inf(f64) else 0.0
            else // (|x| < 1)**±inf
                return if (hy >= 0) 0.0 else std.math.inf(f64);
        }

        if (iy == 0x3ff00000) // y is ±1
            return if (hy >= 0) x else 1.0 / x;

        if (hy == 0x40000000) // y is 2
            return x * x;

        if (hy == 0x40080000) // y is 3
            return x * x * x;

        if (hy == 0x40100000) { // y is 4
            const u: f64 = x * x;
            return u * u;
        }

        if (hy == 0x3fe00000) // y is 0.5
            if (hx >= 0) // x >= 0
                return float.sqrt(x);
    }

    var ax: f64 = float.abs(x);
    if (lx == 0) {
        if (ix == 0x7ff00000 or ix == 0 or ix == 0x3ff00000) { // x is ±inf, ±0, or ±1
            var z: f64 = ax;

            if (hy < 0)
                z = 1.0 / z;

            if (hx < 0) {
                if (((ix -% 0x3ff00000) | yisint) == 0) {
                    z = std.math.nan(f64); // (-1)**non-int = NaN
                } else if (yisint == 1) {
                    z = -z; // (x < 0)**odd = -(|x|**odd)
                }
            }

            return z;
        }
    }

    var n: i32 = types.scast(i32, @as(u32, @bitCast(hx)) >> 31) - 1;

    // (x < 0)**non-int = NaN
    if ((n | yisint) == 0)
        return std.math.nan(f64);

    var sn: f64 = 1.0; // Sign of result
    if ((n | (yisint -% 1)) == 0)
        sn = -1.0;

    var t1: f64 = undefined;
    var t2: f64 = undefined;
    if (iy > 0x41e00000) { // |y| > 2**31
        if (iy > 0x43f00000) { // |y| > 2**64, must over/underflow
            if (ix < 0x3fefffff)
                return if (hy < 0) std.math.inf(f64) else 0.0;

            if (ix > 0x3ff00000)
                return if (hy > 0) std.math.inf(f64) else 0.0;
        }

        // over/underflow if x is not close to 1
        if (ix < 0x3fefffff)
            return sn * (if (hy < 0) std.math.inf(f64) else 0.0);

        if (ix > 0x3ff00000)
            return sn * (if (hy > 0) std.math.inf(f64) else 0.0);

        // Now |1 - x| <= 2**-20, compute log(x) directly
        const t: f64 = ax - 1.0;
        const w: f64 = (t * t) * (0.5 - t * (1.0 / 3.0 - t * 0.25));
        const u: f64 = 1.44269502162933349609e+0 * t; // 1/log(2) * t
        const v: f64 = t * 1.92596299112661746887e-8 - w * 1.44269504088896338700e+0;
        t1 = u + v;
        dbl64.setLowPart(&t1, 0);
        t2 = v - (t1 - u);
    } else {
        n = 0;

        // Take care of subnormal numbers
        if (ix < 0x00100000) {
            ax *= 9007199254740992.0; // 2**53
            n -%= 53;
            ix = @bitCast(dbl64.getHighPart(ax));
        }

        n +%= (ix >> 20) -% 0x3ff;
        const j: i32 = ix & 0x000fffff;

        // Determine interval for ax
        ix = j | 0x3ff00000; // Normalize ix
        var k: i32 = undefined;
        if (j <= 0x3988e) // |x| < sqrt(3/2)
            k = 0
        else if (j < 0xbb67a) // |x| < sqrt(3)
            k = 1
        else {
            k = 0;
            n +%= 1;
            ix -%= 0x00100000;
        }
        dbl64.setHighPart(&ax, @bitCast(ix));

        // Compute s = s_h + s_l = (x - 1)/(x + 1) or (x - 1.5)/(x + 1.5)
        const bp: f64 = if (k == 0) 1.0 else 1.5;
        var u: f64 = ax - bp;
        var v: f64 = 1.0 / (ax + bp);
        const s: f64 = u * v;
        var s_h: f64 = s;
        dbl64.setLowPart(&s_h, 0);

        // t_h = ax + bp
        var t_h: f64 = 0.0;
        dbl64.setHighPart(&t_h, @bitCast(((ix >> 1) | 0x20000000) + 0x00080000 + (k << 18)));
        var t_l: f64 = ax - (t_h - bp);
        const s_l: f64 = v * ((u - s_h * t_h) - s_h * t_l);

        // Compute log(ax)
        var s2: f64 = s * s;
        var r: f64 = s2 * s2 *
            (5.99999999999994648725e-1 + s2 *
                (4.28571428578550184252e-1 + s2 *
                    (3.33333329818377432918e-1 + s2 *
                        (2.72728123808534006489e-1 + s2 *
                            (2.30660745775561754067e-1 + s2 *
                                2.06975017800338417784e-1)))));
        r += s_l * (s_h + s);
        s2 = s_h * s_h;
        t_h = 3.0 + s2 + r;
        dbl64.setLowPart(&t_h, 0);
        t_l = r - ((t_h - 3.0) - s2);

        // u + v = s * (1 + ...)
        u = s_h * t_h;
        v = s_l * t_h + t_l * s;

        // 2/(3 log(2)) * (s + ...)
        var p_h: f64 = u + v;
        dbl64.setLowPart(&p_h, 0);
        const p_l: f64 = v - (p_h - u);
        const z_h: f64 = 9.61796700954437255859e-1 * p_h;
        const z_l: f64 = -7.02846165095275826516e-9 * p_h +
            p_l * 9.61796693925975554329e-1 +
            @as(f64, if (k == 0) 0.0 else 1.35003920212974897128e-8);

        // log2(ax) = (s + ...) * 2/(3 log(2))
        const t: f64 = types.scast(f64, n);
        t1 = (((z_h + z_l) + @as(f64, if (k == 0) 0.0 else 5.84962487220764160156e-1)) + t);
        dbl64.setLowPart(&t1, 0);
        t2 = z_l - (((t1 - t) - @as(f64, if (k == 0) 0.0 else 5.84962487220764160156e-1)) - z_h);
    }

    var y1: f64 = y;
    dbl64.setLowPart(&y1, 0);
    const p_l: f64 = (y - y1) * t1 + y * t2;
    var p_h: f64 = y1 * t1;
    var z: f64 = p_l + p_h;
    var j: i32 = @bitCast(dbl64.getHighPart(z));
    var i: i32 = @bitCast(dbl64.getLowPart(z));
    if (j >= 0x40900000) { // z >= 1024
        if (((j -% 0x40900000) | i) != 0) // z > 1024
            return sn * std.math.inf(f64)
        else if (p_l + 8.0085662595372944372e-17 > z - p_h)
            return sn * std.math.inf(f64);
    } else if ((j & 0x7fffffff) >= 0x4090cc00) { // z <= -1075
        if (((j -% @as(i32, @bitCast(@as(u32, 0xc090cc00)))) | i) != 0) // z < -1075
            return sn * 0.0
        else if (p_l <= z - p_h)
            return sn * 0.0;
    }

    // Compute 2**(p_h + p_l)
    i = j & 0x7fffffff;
    var k: i32 = (i >> 20) -% 0x3ff;
    n = 0;
    if (i > 0x3fe00000) { // |z| > 0.5, set n = [z + 0.5]
        n = j +% (@as(i32, 0x00100000) >> @intCast(k + 1));
        k = ((n & 0x7fffffff) >> 20) -% 0x3ff;
        var t: f64 = 0.0;
        dbl64.setHighPart(&t, @bitCast(n & ~(@as(i32, 0x000fffff) >> @intCast(k))));
        n = ((n & 0x000fffff) | 0x00100000) >> @intCast(20 -% k);

        if (j < 0)
            n = -n;

        p_h -= t;
    }

    var t: f64 = p_l + p_h;
    dbl64.setLowPart(&t, 0);
    const u: f64 = t * 6.93147182464599609375e-1;
    const v: f64 = (p_l - (t - p_h)) * 6.93147180559945286227e-1 + t * -1.90465429995776804525e-9;
    z = u + v;
    const w: f64 = v - (z - u);
    t = z * z;
    t1 = z - t *
        (1.66666666666666019037e-1 + t *
            (-2.77777777770155933842e-3 + t *
                (6.61375632143793436117e-5 + t *
                    (-1.65339022054652515390e-6 + t *
                        4.13813679705723846039e-8))));
    const r: f64 = (z * t1) / (t1 - 2.0) - (w + z * w);
    z = 1.0 - (r - z);
    j = @bitCast(dbl64.getHighPart(z));
    j +%= n << 20;

    if ((j >> 20) <= 0) // Subnormal output
        z = float.scalbn(z, n)
    else
        dbl64.setHighPart(&z, @bitCast(j));

    return sn * z;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld80/e_powl.c
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
fn pow80(x: f80, y: f80) f80 {
    if (y == 0.0)
        return 1.0;

    if (x == 1.0)
        return 1.0;

    if (std.math.isNan(x) or std.math.isNan(y))
        return std.math.nan(f80);

    if (y == 1.0)
        return x;

    if (!std.math.isFinite(y) and x == -1.0)
        return 1.0;

    if (y >= std.math.floatMax(f80)) {
        if (x > 1.0)
            return std.math.inf(f80);

        if (x > 0.0 and x < 1.0)
            return 0.0;

        if (x < -1.0)
            return std.math.inf(f80);

        if (x > -1.0 and x < 0.0)
            return 0.0;
    }

    if (y <= -std.math.floatMax(f80)) {
        if (x > 1.0)
            return 0.0;

        if (x > 0.0 and x < 1.0)
            return std.math.inf(f80);

        if (x < -1.0)
            return 0.0;

        if (x > -1.0 and x < 0.0)
            return std.math.inf(f80);
    }

    if (x >= std.math.floatMax(f80)) {
        if (y > 0.0)
            return std.math.inf(f80);

        return 0.0;
    }

    var w: f80 = float.floor(y);

    // Set iyflg to true if y is an integer
    var iyflg: bool = false;
    if (w == y)
        iyflg = true;

    // Test for odd integer y
    var yoddint: bool = false;
    if (iyflg) {
        var ya: f80 = float.abs(y);
        ya = float.floor(0.5 * ya);
        const yb: f80 = 0.5 * float.abs(w);

        if (ya != yb)
            yoddint = true;
    }

    if (x <= -std.math.floatMax(f80)) {
        if (y > 0.0) {
            if (yoddint)
                return -std.math.inf(f80)
            else
                return std.math.inf(f80);
        }

        if (y < 0.0) {
            if (yoddint)
                return -0.0
            else
                return 0.0;
        }
    }

    var nflg: bool = false; // If x < 0 raised to integer power
    if (x <= 0.0) {
        if (x == 0.0) {
            if (y < 0.0) {
                if (std.math.signbit(x) and yoddint)
                    return -std.math.inf(f80);

                return std.math.inf(f80);
            }

            if (y > 0.0) {
                if (std.math.signbit(x) and yoddint)
                    return -0.0;

                return 0.0;
            }

            if (y == 0.0)
                return 1.0 // 0**0
            else
                return 0.0; // 0**y
        } else {
            if (!iyflg) // (x < 0)**non-int is NaN
                return std.math.nan(f80);

            nflg = true;
        }
    }

    // Integer power of an integer
    if (iyflg) {
        w = float.floor(x);
        if (w == x and (float.abs(y) < 32768.0))
            return powi80(x, types.scast(i32, y));
    }

    var xx: f80 = x;
    if (nflg)
        xx = float.abs(x);

    // Separate significand from exponent
    var i: i32 = undefined;
    xx = float.frexp(xx, &i);
    var e: i32 = i;

    // Find significand in antilogarithm table A_80
    i = 1;
    if (xx <= A_80[17])
        i = 17;
    if (xx <= A_80[types.scast(u32, i + 8)])
        i += 8;
    if (xx <= A_80[types.scast(u32, i + 4)])
        i += 4;
    if (xx <= A_80[types.scast(u32, i + 2)])
        i += 2;
    if (xx >= A_80[1])
        i = -1;
    i += 1;

    // Find (x - A_80[i])/A_80[i]
    // in order to compute log(x/A[i]):
    //
    // log(x) = log(a * x/a) = log(a) + log(x/a)
    //
    // log(x/a) = log(1 + v), v = x/a - 1 = (x - a)/a
    xx -= A_80[types.scast(u32, i)];
    xx -= B_80[types.scast(u32, @divTrunc(i, 2))];
    xx /= A_80[types.scast(u32, i)];

    // Rational approximation for log(1+v):
    //
    // log(1 + v) = v - v**2/2 + v**3 P(v)/Q(v)
    var z: f80 = xx * xx;
    w = xx * (z *
        (1.4000100839971580279335e0 + xx *
            (1.7500123722550302671919e0 + xx *
                (4.9000050881978028599627e-1 + xx *
                    8.3319510773868690346226e-4))) /
        (4.2000302519914740834728e0 + xx *
            (8.4000598057587009834666e0 + xx *
                (5.2500282295834889175431e0 + xx *
                    1.0000000000000000000000e0))));
    w -= float.ldexp(z, -1);

    // Convert to base 2 logarithm:
    // multiply by log2(e)
    z = 0.44269504088896340735992 * w;
    z += w;
    z += 0.44269504088896340735992 * xx;
    z += xx;

    // Compute exponent term of the base 2 logarithm
    w = types.scast(f80, -i);
    w = float.ldexp(w, -5);
    w += types.scast(f80, e);
    // Now base 2 log of x is w + z

    // Multiply base 2 log by y, in extended precision.
    // Separate y into large part ya
    // and small part yb less than 1/5
    const ya: f80 = reduc80(y);
    const yb: f80 = y - ya;

    // (w + z) * (ya + yb) = w * ya + w * yb + z * y
    const f: f80 = z * y + w * yb;
    const fa: f80 = reduc80(f);
    const fb: f80 = f - fa;

    const g: f80 = fa + w * ya;
    const ga: f80 = reduc80(g);
    const gb: f80 = g - ga;

    const h: f80 = fb + gb;
    const ha: f80 = reduc80(h);
    w = float.ldexp(ga + ha, 5);

    // Test the power of 2 for overflow
    if (w > 32.0 * 16384.0)
        return std.math.inf(f80);

    if (w < -32.0 * (16384.0 + 64.0))
        return 0.0;

    e = types.scast(i32, w);
    var hb: f80 = h - ha;

    if (hb > 0.0) {
        e +%= 1;
        hb -= 1.0 / 32.0;
    }

    // Now the product y * log2(x) = hb + e/5.
    //
    // Compute base 2 exponential of hb,
    // where -0.0625 <= hb <= 0.
    z = hb *
        (6.9314718055994530931447e-1 + hb *
            (2.4022650695910062854352e-1 + hb *
                (5.5504108664798463044015e-2 + hb *
                    (9.6181291046036762031786e-3 + hb *
                        (1.3333556028915671091390e-3 + hb *
                            (1.5402715328927013076125e-4 + hb *
                                1.5089970579127659901157e-5))))));

    // Express e/32 as an integer plus a negative number of (1/32)ths.
    // Find lookup table entry for the fractional power of 2
    if (e < 0)
        i = 0
    else
        i = 1;
    i +%= @divTrunc(e, 32);
    e = 32 *% i -% e;
    w = A_80[types.scast(u32, e)];
    z *= w;
    z += w;
    z = float.ldexp(z, i);

    if (nflg) {
        // For negative x,
        // find out if the integer exponent
        // is odd or even
        w = float.ldexp(y, -1);
        w = float.floor(w);
        w = float.ldexp(w, 1);
        if (w != y)
            z = -z;
    }

    return z;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_logl.c
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
fn pow128(x: f128, y: f128) f128 {
    const p: ldbl128.Parts32 = .fromFloat(x);
    const hx: i32 = @bitCast(p.mswhi);
    var ix: u32 = p.mswhi & 0x7fffffff;

    const q: ldbl128.Parts32 = .fromFloat(y);
    const hy: i32 = @bitCast(q.mswhi);
    const iy: u32 = q.mswhi & 0x7fffffff;

    if ((iy | q.mswlo | q.lswhi | q.lswlo) == 0) // x**0 = 1
        return 1.0;

    if (x == 1.0) // 1**y = 1, even if y is NaN
        return 1.0;

    if (x == -1.0 and iy == 0x7fff0000 and
        (q.mswlo | q.lswhi | q.lswlo) == 0) // (-1)**+-inf
        return 1.0;

    // NaN if x or y is NaN, with y != 0
    if (ix > 0x7fff0000 or
        (ix == 0x7fff0000 and (p.mswlo | p.lswhi | p.lswlo) != 0) or
        iy > 0x7fff0000 or
        (iy == 0x7fff0000 and (q.mswlo | q.lswhi | q.lswlo) != 0))
        return std.math.nan(f128);

    // Determine if y is an odd int when x < 0
    // yisint = 0 -> y is not an integer
    // yisint = 1 -> y is an odd integer
    // yisint = 2 -> y is an even integer
    var yisint: i32 = 0;
    if (hx < 0) {
        if (iy >= 0x40700000) // 2^113
            yisint = 2 // even integer y
        else if (iy >= 0x3fff0000) { // 1.0
            if (float.floor(y) == y) {
                const z: f128 = 0.5 * y;
                if (float.floor(z) == z)
                    yisint = 2
                else
                    yisint = 1;
            }
        }
    }

    if ((q.mswlo | q.lswhi | q.lswlo) == 0) {
        if (iy == 0x7fff0000) { // y is ±inf
            if (((ix -% 0x3fff0000) | p.mswlo | p.lswhi | p.lswlo) == 0)
                return y - y // ±1**inf is NaN
            else if (ix >= 0x3fff0000) // (|x| > 1)**±inf
                return if (hy >= 0) y else 0.0
            else // (|x| < 1)**±inf
                return if (hy < 0) -y else 0.0;
        }

        if (iy == 0x3fff0000) { // y is ±1
            if (hy < 0)
                return 1.0 / x
            else
                return x;
        }

        if (hy == 0x40000000)
            return x * x; // y is 2

        if (hy == 0x3ffe0000) { // y is 0.5
            if (hx >= 0) // x >= 0
                return float.sqrt(x);
        }
    }

    var ax: f128 = float.abs(x);

    if ((p.mswlo | p.lswhi | p.lswlo) == 0) {
        if (ix == 0x7fff0000 or ix == 0 or ix == 0x3fff0000) {
            var z: f128 = ax; // x is ±0,±inf,±1

            if (hy < 0)
                z = 1.0 / z; // z = (1/|x|)

            if (hx < 0) {
                if (((ix -% 0x3fff0000) | @as(u32, @bitCast(yisint))) == 0)
                    z = std.math.nan(f128) // (-1)**non-int = NaN
                else if (yisint == 1)
                    z = -z; // (x < 0)**odd = -(|x|**odd)

            }

            return z;
        }
    }

    if (((types.scast(i32, @as(u32, @bitCast(hx)) >> 31) - 1) | yisint) == 0)
        return std.math.nan(f128); // (x < 0)**non-int = NaN

    var sn: f128 = 1.0;
    if (((types.scast(i32, @as(u32, @bitCast(hx)) >> 31) - 1) | (yisint -% 1)) == 0)
        sn = -1.0;

    // |y| is huge.
    // 2**-16495 = 1/2 of smallest representable value.
    // If (1 - 1/131072)**y underflows, y > 1.4986e9
    if (iy > 0x401d654b) {
        // if (1 - 2**-113)**y underflows, y > 1.1873e38
        if (iy > 0x407d654b) {
            if (ix <= 0x3ffeffff)
                return if (hy < 0) std.math.inf(f128) else 0.0;

            if (ix >= 0x3fff0000)
                return if (hy > 0) std.math.inf(f128) else 0.0;
        }

        // over/underflow if x is not close to 1
        if (ix < 0x3ffeffff)
            return sn * (if (hy < 0) std.math.inf(f128) else 0.0);

        if (ix > 0x3fff0000)
            return sn * (if (hy > 0) std.math.inf(f128) else 0.0);
    }

    var n: i32 = 0;
    // Take care subnormal numbers
    if (ix < 0x00010000) {
        ax *= 1.0384593717069655257060992658440192e34; // 2**113
        n -%= 113;
        const o: ldbl128.Parts32 = .fromFloat(ax);
        ix = o.mswhi;
    }

    n +%= types.scast(i32, ix >> 16) -% 0x3fff;
    var j: i32 = @bitCast(ix & 0x0000ffff);
    // Determine interval
    ix = @bitCast(j | 0x3fff0000); // Normalize ix
    var k: i32 = undefined;
    if (j <= 0x3988)
        k = 0 // |x| < sqrt(3/2)
    else if (j < 0xbb67)
        k = 1 // |x| < sqrt(3)
    else {
        k = 0;
        n +%= 1;
        ix -%= 0x00010000;
    }

    var o: ldbl128.Parts32 = .fromFloat(ax);
    o.mswhi = ix;
    ax = o.toFloat();

    // Compute s = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5)
    const bp: f128 = if (k == 0) 1.0 else 1.5;
    var u: f128 = ax - bp;
    var v: f128 = 1.0 / (ax + bp);
    const s: f128 = u * v;
    var s_h: f128 = s;

    o = .fromFloat(s_h);
    o.lswlo = 0;
    o.lswhi &= 0xf8000000;
    s_h = o.toFloat();

    // t_h = ax + bp High
    var t_h: f128 = ax + bp;
    o = .fromFloat(t_h);
    o.lswlo = 0;
    o.lswhi &= 0xf8000000;
    t_h = o.toFloat();
    var t_l: f128 = ax - (t_h - bp);
    const s_l: f128 = v * ((u - s_h * t_h) - s_h * t_l);

    // Compute log(ax)
    var s2: f128 = s * s;
    u = -3.0779177200290054398792536829702930623200e1 + s2 *
        (6.5135778082209159921251824580292116201640e1 + s2 *
            (-4.6312921812152436921591152809994014413540e1 + s2 *
                (1.2510208195629420304615674658258363295208e1 + s2 *
                    -9.9266909031921425609179910128531667336670e-1)));
    v = -5.129862866715009066465422805058933131960e1 + s2 *
        (1.452015077564081884387441590064272782044e2 + s2 *
            (-1.524043275549860505277434040464085593165e2 + s2 *
                (7.236063513651544224319663428634139768808e1 + s2 *
                    (-1.494198912340228235853027849917095580053e1 + s2))));
    var r: f128 = s2 * s2 * u / v;
    r += s_l * (s_h + s);
    s2 = s_h * s_h;
    t_h = 3.0 + s2 + r;
    o = .fromFloat(t_h);
    o.lswlo = 0;
    o.lswhi &= 0xf8000000;
    t_h = o.toFloat();
    t_l = r - ((t_h - 3.0) - s2);

    // u + v = s * (1 + ...)
    u = s_h * t_h;
    v = s_l * t_h + t_l * s;

    // 2 / (3 log(2)) * (s + ...)
    var p_h: f128 = u + v;
    o = .fromFloat(p_h);
    o.lswlo = 0;
    o.lswhi &= 0xf8000000;
    p_h = o.toFloat();
    var p_l: f128 = v - (p_h - u);
    const z_h: f128 = 9.6179669392597555432899980587535537779331e-1 * p_h;
    const z_l: f128 = 5.0577616648125906047157785230014751039424e-17 * p_h +
        p_l * 9.6179669392597560490661645400126142495110e-1 +
        @as(f128, if (k == 0) 0.0 else 1.0579781240112554492329533686862998106046e-16);

    // log2(ax) = (s + ...) * 2 / (3 * log(2)) = n + dp_h + z_h + z_l
    var t: f128 = types.scast(f128, n);
    var t1: f128 = (((z_h + z_l) + @as(f128, if (k == 0) 0.0 else 5.8496250072115607565592654282227158546448e-1)) + t);
    o = .fromFloat(t1);
    o.lswlo = 0;
    o.lswhi &= 0xf8000000;
    t1 = o.toFloat();
    const t2: f128 = z_l - (((t1 - t) - @as(f128, if (k == 0) 0.0 else 5.8496250072115607565592654282227158546448e-1)) - z_h);

    // Split up y into y1 + y2 and compute (y1 + y2) * (t1 + t2)
    var y1: f128 = y;
    o = .fromFloat(y1);
    o.lswlo = 0;
    o.lswhi &= 0xf8000000;
    y1 = o.toFloat();
    p_l = (y - y1) * t1 + y * t2;
    p_h = y1 * t1;
    var z: f128 = p_l + p_h;
    o = .fromFloat(z);
    j = @bitCast(o.mswhi);
    if (j >= 0x400d0000) { // z >= 16384
        // if z > 16384
        if ((@as(u32, @bitCast(j -% 0x400d0000)) | o.mswlo | o.lswhi | o.lswlo) != 0)
            return sn * std.math.inf(f128)
        else if (p_l + 8.0085662595372944372e-17 > z - p_h)
            return sn * std.math.inf(f128);
    } else if ((j & 0x7fffffff) >= 0x400d01b9) { // z <= -16495
        // z < -16495
        if (((@as(u32, @bitCast(j)) -% 0xc00d01bc) | o.mswlo | o.lswhi | o.lswlo) != 0)
            return sn * 0.0
        else if (p_l <= z - p_h)
            return sn * 0.0;
    }

    // Compute 2**(p_h + p_l)
    const i: i32 = j & 0x7fffffff;
    k = (i >> 16) -% 0x3fff;
    n = 0;
    if (i > 0x3ffe0000) { // if |z| > 0.5, set n = [z + 0.5]
        n = types.scast(i32, float.floor(z + 0.5));
        t = types.scast(f128, n);
        p_h -= t;
    }

    t = p_l + p_h;
    o = .fromFloat(t);
    o.lswlo = 0;
    o.lswhi &= 0xf8000000;
    t = o.toFloat();
    u = t * 6.9314718055994528622676398299518041312695e-1;
    v = (p_l - (t - p_h)) * 6.9314718055994530941723212145817656807550e-1 +
        t * 2.3190468138462996154948554638754786504121e-17;
    z = u + v;
    const w: f128 = v - (z - u);

    // exp(z)
    t = z * z;
    u = 5.081801691915377692446852383385968225675e8 + t *
        (9.360895299872484512023336636427675327355e6 + t *
            (4.213701282274196030811629773097579432957e4 + t *
                (5.201006511142748908655720086041570288182e1 + t *
                    9.088368420359444263703202925095675982530e-3)));
    v = 3.049081015149226615468111430031590411682e9 + t *
        (1.069833887183886839966085436512368982758e8 + t *
            (8.259257717868875207333991924545445705394e5 + t *
                (1.872583833284143212651746812884298360922e3 + t)));
    t1 = z - t * u / v;
    r = (z * t1) / (t1 - 2.0) - (w + z * w);
    z = 1.0 - (r - z);
    o = .fromFloat(z);
    j = @bitCast(o.mswhi);
    j +%= (n << 16);
    if ((j >> 16) <= 0) {
        z = float.scalbn(z, n); // Subnormal output
    } else {
        o.mswhi = @bitCast(j);
        z = o.toFloat();
    }

    return sn * z;
}

fn reduc80(x: f80) f80 {
    var t: f80 = float.ldexp(x, 5);
    t = float.floor(t);
    t = float.ldexp(t, -5);
    return t;
}

fn powi80(x: f80, nn: i32) f80 {
    if (x == 0.0) {
        if (nn == 0)
            return 1.0
        else if (nn < 0)
            return std.math.floatMax(f80)
        else
            return 0.0;
    }

    if (nn == 0)
        return 1.0;

    var asign: i32 = 0;
    var xx: f80 = x;
    if (x < 0.0) {
        asign = -1;
        xx = -xx;
    }

    var sign: i32 = 1;
    var n: i32 = nn;
    if (nn < 0) {
        sign = -1;
        n = -nn;
    }

    // Calculate approximate logarithm of answer
    var s: f80 = xx;
    var lx: i32 = undefined;
    s = float.frexp(s, &lx);
    const e: i32 = (lx -% 1) * n;
    if (e == 0 or e > 64 or e < -64) {
        s = (s - 7.0710678118654752e-1) / (s + 7.0710678118654752e-1);
        s = (2.9142135623730950 * s - 0.5 + types.scast(f80, lx)) * types.scast(f80, nn) * 6.9314718055994530941723e-1;
    } else {
        s = 6.9314718055994530941723e-1 * types.scast(f80, e);
    }

    if (s > 1.1356523406294143949492e4)
        return std.math.inf(f80);

    if (s < -1.13994985314888605586758e4)
        return 0.0;

    // Handle tiny denormal answer, but with less accuracy
    // since roundoff error in 1.0/x will be amplified.
    // The precise demarcation should be the gradual underflow threshold.
    if (s < (-1.1356523406294143949492e4 + 2.0)) {
        xx = 1.0 / xx;
        sign = -sign;
    }

    // First bit of the power
    var y: f80 = undefined;
    if (n & 1 != 0)
        y = xx
    else {
        y = 1.0;
        asign = 0;
    }

    var ww: f80 = xx;
    n >>= 1;
    while (n != 0) : (n >>= 1) {
        ww = ww * ww; // arg to the 2-to-the-kth power
        if (n & 1 != 0) // if that bit is set, then include in product
            y *= ww;
    }

    if (asign != 0)
        y = -y; // odd power of negative number

    if (sign < 0)
        y = 1.0 / y;

    return y;
}

const A_80: [33]f80 = .{
    1.0000000000000000000000e0,
    9.7857206208770013448287e-1,
    9.5760328069857364691013e-1,
    9.3708381705514995065011e-1,
    9.1700404320467123175367e-1,
    8.9735453750155359320742e-1,
    8.7812608018664974155474e-1,
    8.5930964906123895780165e-1,
    8.4089641525371454301892e-1,
    8.2287773907698242225554e-1,
    8.0524516597462715409607e-1,
    7.8799042255394324325455e-1,
    7.7110541270397041179298e-1,
    7.5458221379671136985669e-1,
    7.3841307296974965571198e-1,
    7.2259040348852331001267e-1,
    7.0710678118654752438189e-1,
    6.9195494098191597746178e-1,
    6.7712777346844636413344e-1,
    6.6261832157987064729696e-1,
    6.4841977732550483296079e-1,
    6.3452547859586661129850e-1,
    6.2092890603674202431705e-1,
    6.0762367999023443907803e-1,
    5.9460355750136053334378e-1,
    5.8186242938878875689693e-1,
    5.6939431737834582684856e-1,
    5.5719337129794626814472e-1,
    5.4525386633262882960438e-1,
    5.3357020033841180906486e-1,
    5.2213689121370692017331e-1,
    5.1094857432705833910408e-1,
    5.0000000000000000000000e-1,
};

const B_80: [17]f80 = .{
    0.0000000000000000000000e0,
    2.6176170809902549338711e-20,
    -1.0126791927256478897086e-20,
    1.3438228172316276937655e-21,
    1.2207982955417546912101e-20,
    -6.3084814358060867200133e-21,
    1.3164426894366316434230e-20,
    -1.8527916071632873716786e-20,
    1.8950325588932570796551e-20,
    1.5564775779538780478155e-20,
    6.0859793637556860974380e-21,
    -2.0208749253662532228949e-20,
    1.4966292219224761844552e-20,
    3.3540909728056476875639e-21,
    -8.6987564101742849540743e-22,
    -1.2327176863327626135542e-20,
    0.0000000000000000000000e0,
};
