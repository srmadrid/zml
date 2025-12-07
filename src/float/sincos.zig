const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const Coerce = types.Coerce;
const float = @import("../float.zig");

const rem_pio2 = @import("rem_pio2.zig");
const rem_pio2_32 = rem_pio2.rem_pio2_32;
const rem_pio2_64 = rem_pio2.rem_pio2_64;

const dbl64 = @import("dbl64.zig");

pub inline fn sincos(x: anytype) struct { sinx: EnsureFloat(@TypeOf(x)), cosx: EnsureFloat(@TypeOf(x)) } {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.sincos: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (@TypeOf(x)) {
        f16 => {
            var s: f32 = undefined;
            var c: f32 = undefined;
            sincos32(types.scast(f32, x), &s, &c);

            return .{
                .sinx = types.scast(f16, s),
                .cosx = types.scast(f16, c),
            };
        },
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_sincosf.c
            var s: f32 = undefined;
            var c: f32 = undefined;
            sincos32(types.scast(f32, x), &s, &c);

            return .{
                .sinx = s,
                .cosx = c,
            };
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_sincos.c
            var s: f64 = undefined;
            var c: f64 = undefined;
            sincos64(types.scast(f64, x), &s, &c);

            return .{
                .sinx = s,
                .cosx = c,
            };
        },
        f80 => {
            //
            var s: f128 = undefined;
            var c: f128 = undefined;
            sincos128(types.scast(f128, x), &s, &c);

            return .{
                .sinx = types.scast(f80, s),
                .cosx = types.scast(f80, c),
            };
        },
        f128 => {
            //
            var s: f128 = undefined;
            var c: f128 = undefined;
            sincos128(types.scast(f128, x), &s, &c);

            return .{
                .sinx = s,
                .cosx = c,
            };
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_sincosf.c
//
// Original copyright notice:
// ====================================================
// This file is derived from fdlibm:
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
//
// ====================================================
// Copyright (C) 2013 Elliot Saba
// Developed at the University of Washington
//
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn sincos32(x: f32, sinx: *f32, cosx: *f32) void {
    sinx.* = x;
    cosx.* = x;

    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    if (ix <= 0x3f490fda) { // |x| ~<= pi/4
        if (ix < 0x39800000) { // |x| < 2**-12
            // Check if x is exactly zero
            if (types.scast(i32, x) == 0) {
                sinx.* = x;
                cosx.* = 1.0;
                return;
            }
        }

        return k_sincos32(x, sinx, cosx);
    }

    if (ix <= 0x407b53d1) { // |x| ~<= 5 * pi/4
        if (ix <= 0x4016cbe3) { // |x| ~<= 3 * pi/4
            if (hx > 0) {
                k_sincos32(1.57079632679489661923 - types.scast(f64, x), cosx, sinx);
            } else {
                k_sincos32(1.57079632679489661923 + types.scast(f64, x), cosx, sinx);
                sinx.* = -sinx.*;
            }
        } else {
            if (hx > 0) {
                k_sincos32(2.0 * 1.57079632679489661923 - types.scast(f64, x), sinx, cosx);
                cosx.* = -cosx.*;
            } else {
                k_sincos32(-2.0 * 1.57079632679489661923 - types.scast(f64, x), sinx, cosx);
                cosx.* = -cosx.*;
            }
        }

        return;
    }

    if (ix <= 0x40e231d5) { // |x| ~<= 9 * pi/4
        if (ix <= 0x40afeddf) { // |x| ~> 7 * pi/4
            if (hx > 0) {
                k_sincos32(types.scast(f64, x) - 3.0 * 1.57079632679489661923, cosx, sinx);
                sinx.* = -sinx.*;
            } else {
                k_sincos32(types.scast(f64, x) + 3.0 * 1.57079632679489661923, cosx, sinx);
                cosx.* = -cosx.*;
            }
        } else {
            if (hx > 0) {
                k_sincos32(types.scast(f64, x) - 4.0 * 1.57079632679489661923, sinx, cosx);
            } else {
                k_sincos32(types.scast(f64, x) + 4.0 * 1.57079632679489661923, sinx, cosx);
            }
        }

        return;
    } else if (ix >= 0x7f800000) { // cos(Inf or NaN) is NaN
        sinx.* = x - x;
        cosx.* = x - x;
        return;
    } else {
        // General argument reduction needed
        var y: f64 = undefined;
        const n: i32 = rem_pio2_32(x, &y);

        switch (n & 3) {
            0 => k_sincos32(y, sinx, cosx),
            1 => k_sincos32(-y, cosx, sinx),
            2 => {
                k_sincos32(-y, sinx, cosx);
                cosx.* = -cosx.*;
            },
            else => {
                k_sincos32(-y, cosx, sinx);
                sinx.* = -sinx.*;
                cosx.* = -cosx.*;
            },
        }

        return;
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_sincos.c
//
// Original copyright notice:
// ====================================================
// This file is derived from fdlibm:
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
//
// ====================================================
// Copyright (C) 2013 Elliot Saba. All rights reserved.
//
// Developed at the University of Washington.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn sincos64(x: f64, sinx: *f64, cosx: *f64) void {
    const ix: i32 = @bitCast(dbl64.getHighPart(x) & 0x7fffffff);

    if (ix <= 0x3fe921fb) { // |x| ~< pi/4
        // Check for small x for sin and cos
        if (ix < 0x3e46a09e) {
            // Check for exact zero
            if (types.scast(i32, x) == 0) {
                sinx.* = x;
                cosx.* = 1.0;
                return;
            }
        }

        // Call kernel function with 0 extra
        k_sincos64(x, 0.0, 0, sinx, cosx);
    } else if (ix >= 0x7ff00000) {
        // sincos(Inf or NaN) is NaN
        sinx.* = x - x;
        cosx.* = x - x;
    } else { // Argument reduction needed
        // Calculate remainer, then sub out to kernel */
        var y: [2]f64 = undefined;
        const n: i32 = rem_pio2_64(x, &y);
        k_sincos64(y[0], y[1], 1, sinx, cosx);

        // Figure out permutation of sin/cos outputs to true outputs
        switch (n & 3) {
            0 => {},
            1 => {
                const tmp: f64 = sinx.*;
                sinx.* = cosx.*;
                cosx.* = -tmp;
            },
            2 => {
                sinx.* = -sinx.*;
                cosx.* = -cosx.*;
            },
            else => {
                const tmp: f64 = sinx.*;
                sinx.* = -cosx.*;
                cosx.* = tmp;
            },
        }
    }
}

fn sincos128(x: f128, sinx: *f128, cosx: *f128) void {
    sinx.* = float.sin(x);
    cosx.* = float.cos(x);
}

fn k_sincos32(x: f64, sinx: *f32, cosx: *f32) void {
    const z: f64 = x * x;
    const w: f64 = z * z;

    // cos-specific computation; equivalent to calling
    // k_cos32(x, y) and storing in cosx
    var r: f64 = -0x16c087e80f1e27.0p-62 + z * 0x199342e0ee5069.0p-68;
    cosx.* = types.scast(f32, ((1.0 + z * -0x1ffffffd0c5e81.0p-54) + w * 0x155553e1053a42.0p-57) + (w * z) * r);

    // sin-specific computation; equivalent to calling
    // k_sin32(x, y, 1) and storing in sinx
    r = -0x1a00f9e2cae774.0p-65 + z * 0x16cd878c3b46a7.0p-71;
    const v: f64 = z * x;
    sinx.* = types.scast(f32, (x + v * (-0x15555554cbac77.0p-55 + z * 0x111110896efbb2.0p-59)) + v * w * r);
}

fn k_sincos64(x: f64, y: f64, iy: i32, sinx: *f64, cosx: *f64) void {
    const z: f64 = x * x;
    const w: f64 = z * z;

    // cos-specific computation; equivalent to calling
    // k_cos64(x, y) and storing in cosx
    var r: f64 = z * (4.16666666666666019037e-2 +
        z * (-1.38888888888741095749e-3 + z * 2.48015872894767294178e-5)) +
        w * w * (-2.75573143513906633035e-7 +
            z * (2.08757232129817482790e-9 + z * -1.13596475577881948265e-11));
    const hz: f64 = 0.5 * z;
    var v: f64 = 1.0 - hz;

    cosx.* = v + (((1.0 - v) - hz) + (z * r - x * y));

    // sin-specific computation; equivalent to calling
    // k_sin64(x, y, iy) and storing in sinx
    r = 8.33333333332248946124e-3 +
        z * (-1.98412698298579493134e-4 +
            z * 2.75573137070700676789e-6) +
        z * w * (-2.50507602534068634195e-8 +
            z * 1.58969099521155010221e-10);
    v = z * x;

    if (iy == 0)
        sinx.* = x + v * (-1.66666666666666324348e-1 + z * r)
    else
        sinx.* = x - ((z * (0.5 * y - v * r) - y) - v * -1.66666666666666324348e-1);
}
