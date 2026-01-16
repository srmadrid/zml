const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Cbrt(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.cbrt: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

pub inline fn cbrt(x: anytype) Cbrt(@TypeOf(x)) {
    switch (Cbrt(@TypeOf(x))) {
        f16 => return types.scast(f16, cbrt32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_cbrtf.c
            return cbrt32(types.scast(f32, x));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_cbrt.c
            return cbrt64(types.scast(f64, x));
        },
        f80 => {
            //
            // return cbrt80(types.scast(f80, x));
            return types.scast(f80, cbrt128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_cbrtl.c
            return cbrt128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_cbrtf.c
//
// Original copyright notice:
// s_cbrtf.c -- float version of s_cbrt.c.
// Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
// Debugged and optimized by Bruce D. Evans.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn cbrt32(x: f32) f32 {
    var hx: i32 = @bitCast(x);
    const sign: u32 = @as(u32, @bitCast(hx)) & 0x80000000;
    hx &= 0x7fffffff;

    if (hx >= 0x7f800000) // cbrt(NaN) is NaN, cbrt(Inf) is Inf
        return x + x;

    // Rough cbrt to 5 bits
    var t: f32 = undefined;
    if (hx < 0x00800000) { // zero or subnormal
        if (hx == 0)
            return (x); // cbrt(±0) is itself

        t = @bitCast(@as(u32, 0x4b800000)); // Set t = 2**24
        t *= x;
        const high: u32 = @bitCast(t);
        t = @bitCast(sign | (@divTrunc(high & 0x7fffffff, 3) +% 642849266));
    } else {
        t = @bitCast(@as(i32, @bitCast(sign)) | (@divTrunc(hx, 3) +% 709958130));
    }

    // First step Newton iteration (solving t * t - x/t == 0) to 16 bits.  In
    // double precision so that its terms can be arranged for efficiency
    // without causing overflow or underflow.
    var T: f64 = types.scast(f64, t);
    var r: f64 = T * T * T;
    T = T * (2.0 * types.scast(f64, x) + r) / (types.scast(f64, x) + r + r);

    // Second step Newton iteration to 47 bits.  In double precision for
    // efficiency and accuracy.
    r = T * T * T;
    T = T * (2.0 * types.scast(f64, x) + r) / (types.scast(f64, x) + r + r);

    // Rounding to 24 bits is perfect in round-to-nearest mode
    return types.scast(f32, T);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_cbrt.c
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
// Optimized by Bruce D. Evans.
fn cbrt64(x: f64) f64 {
    var hx: i32 = @bitCast(dbl64.getHighPart(x));
    const low: u32 = dbl64.getLowPart(x);
    const sign: u32 = @as(u32, @bitCast(hx)) & 0x80000000;
    hx &= 0x7fffffff;

    if (hx >= 0x7ff00000) // cbrt(NaN) is NaN, cbrt(Inf) is Inf
        return x + x;

    // Rough cbrt to 5 bits:
    //    cbrt(2**e * (1 + m) ~= 2**(e/3) * (1 + (e%3 + m)/3)
    // where e is integral and >= 0, m is real and in [0, 1), and "/" and
    // "%" are integer division and modulus with rounding towards minus
    // infinity.  The RHS is always >= the LHS and has a maximum relative
    // error of about 1 in 16.  Adding a bias of -0.03306235651 to the
    // (e%3 + m)/3 term reduces the error to about 1 in 32. With the IEEE
    // floating point representation, for finite positive normal values,
    // ordinary integer divison of the value in bits magically gives
    // almost exactly the RHS of the above provided we first subtract the
    // exponent bias (1023 for doubles) and later add it back.  We do the
    // subtraction virtually to keep e >= 0 so that ordinary integer
    // division rounds towards minus infinity; this is also efficient.
    var t: f64 = 0.0;
    if (hx < 0x00100000) { // Zero or subnormal
        if ((@as(u32, @bitCast(hx)) | low) == 0)
            return x; // cbrt(±0) is itself

        dbl64.setHighPart(&t, 0x43500000); // Set t = 2**54
        t *= x;
        const high: u32 = dbl64.getHighPart(t);
        t = dbl64.Parts.toFloat(.{
            .msw = @bitCast(sign | (@divTrunc(high & 0x7fffffff, 3) +% 696219795)),
            .lsw = 0,
        });
    } else {
        t = dbl64.Parts.toFloat(.{
            .msw = @bitCast(@as(i32, @bitCast(sign)) | (@divTrunc(hx, 3) +% 715094163)),
            .lsw = 0,
        });
    }

    // New cbrt to 23 bits:
    //    cbrt(x) = t * cbrt(x/t**3) ~= t * P(t**3/x)
    // where P(r) is a polynomial of degree 4 that approximates 1/cbrt(r)
    // to within 2**-23.5 when |r - 1| < 1/10.  The rough approximation
    // has produced t such than |t/cbrt(x) - 1| ~< 1/32, and cubing this
    // gives us bounds for r = t**3/x.
    var r: f64 = (t * t) * (t / x);
    t = t * ((1.87595182427177009643 + r * (-1.88497979543377169875 + r * 1.621429720105354466140)) +
        ((r * r) * r) * (-0.758397934778766047437 + r * 0.145996192886612446982));

    // Round t away from zero to 23 bits (sloppily except for ensuring that
    // the result is larger in magnitude than cbrt(x) but not much more than
    // 2 23-bit ulps larger).  With rounding towards zero, the error bound
    // would be ~5/6 instead of ~4/6.  With a maximum error of 2 23-bit ulps
    // in the rounded t, the infinite-precision error in the Newton
    // approximation barely affects third digit in the final error
    // 0.667; the error in the rounded t can be up to about 3 23-bit ulps
    // before the final error is larger than 0.667 ulps.
    var u: u64 = @bitCast(t);
    u = (u +% 0x80000000) & 0xffffffffc0000000;
    t = @bitCast(u);

    // One step Newton iteration to 53 bits with error < 0.667 ulps
    const s: f64 = t * t; // t*t is exact
    r = x / s; // error <= 0.5 ulps; |r| < |t|
    const w: f64 = t + t; // t+t is exact
    r = (r - t) / (w + r); // r-t is exact; w+r ~= 3*t
    t = t + t * r; // error <= 0.5 + 0.5/3 + epsilon

    return t;
}

fn cbrt80(x: f80) f80 {
    _ = x;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_cbrtl.c
//
// Original copyright notice:
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
// Copyright (c) 2009-2011, Bruce D. Evans, Steven G. Kargl, David Schultz.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
//
// The argument reduction and testing for exceptional cases was
// written by Steven G. Kargl with input from Bruce D. Evans
// and David A. Schultz.
fn cbrt128(x: f128) f128 {
    var u: ldbl128.ShapeSplit = .fromFloat(x);
    var k: i32 = @bitCast(@as(u32, @intCast(u.exponent)));
    const sign: u1 = u.sign;

    // cbrt(NaN) is NaN, cbrt(Inf) is Inf
    if (k == 0x7fff)
        return x + x;

    if (k == 0) {
        if ((u.mantissa_high | u.mantissa_low) == 0)
            return x; // cbrt(±0) is itself

        // Adjust subnormal numbers
        u = .fromFloat(u.toFloat() * 0x1.0p514);
        k = @bitCast(@as(u32, @intCast(u.exponent)));
        k -%= (16384 - 1) + 514;
    } else {
        k -%= 16384 - 1;
    }

    u.exponent = 0x3fff;
    u.sign = 0;
    var v: ldbl128.ShapeSplit = .fromFloat(1.0);
    var xx: f128 = u.toFloat();
    switch (@rem(k, 3)) {
        1, -2 => {
            xx *= 2.0;
            k -%= 1;
        },
        2, -1 => {
            xx *= 4.0;
            k -%= 2;
        },
        else => {},
    }

    v.sign = sign;
    v.exponent = @bitCast(0x3fff +% @as(i15, @intCast(@divTrunc(k, 3))));

    // The following is the guts of cbrt32, with the handling of
    // special values removed and extra care for accuracy not taken,
    // but with most of the extra accuracy not discarded.

    // Rough cbrt to 5 bits
    const fx: f32 = types.scast(f32, xx);
    const hx: u32 = @bitCast(fx);
    const ft: f32 = @bitCast((hx & 0x7fffffff) / 3 +% 709958130);

    // 16 bit estimate
    const dx: f64 = types.scast(f64, xx);
    var dt: f64 = types.scast(f64, ft);
    var dr: f64 = dt * dt * dt;
    dt = dt * (2.0 * dx + dr) / (dx + 2.0 * dr);

    // 47 bit estimate
    dr = dt * dt * dt;
    dt = dt * (2.0 * dx + dr) / (dx + 2.0 * dr);

    // Final step Newton iteration to 64 or 113 bits with
    // error < 0.667 ulps
    var t: f128 = types.scast(f128, dt);
    const s: f128 = t * t; // t * t is exact
    var r: f128 = xx / s; // error <= 0.5 ulps; |r| < |t|
    const w: f128 = t + t; // t + t is exact
    r = (r - t) / (w + r); // r - t is exact; w + r ~= 3 * t
    t = t + t * r; // error <= 0.5 + 0.5/3 + epsilon

    return v.toFloat() * t;
}
