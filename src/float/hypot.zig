const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Hypot(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float))
        @compileError("zml.float.hypot: x and y must be bools, ints or floats, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.EnsureFloat(types.Coerce(X, Y));
}

/// Calculates the hypotenuse $\sqrt{x^2 + y^2}$ of two operands of float, int
/// or bool types. The result type is determined by coercing the operand types,
/// and coercing the coerced type to float, and the operation is performed by
/// casting both operands to the result type, then calculating the hypotenuse.
///
/// ## Signature
/// ```zig
/// float.hypot(x: X, y: Y) float.Hypot(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `float.Hypot(@TypeOf(x), @TypeOf(y))`: The hypotenuse of `x` and `y`.
pub inline fn hypot(x: anytype, y: anytype) float.Hypot(@TypeOf(y), @TypeOf(x)) {
    switch (float.Hypot(@TypeOf(x), @TypeOf(y))) {
        f16 => return types.scast(f16, hypot32(types.scast(f32, x), types.scast(f32, y))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_hypotf.c
            return hypot32(types.scast(f32, x), types.scast(f32, y));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_hypot.c
            return hypot64(types.scast(f64, x), types.scast(f64, y));
        },
        f80 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld80/e_hypotl.c
            // return hypot80(types.scast(f80, x), types.scast(f80, y));
            return types.scast(f80, hypot128(types.scast(f128, x), types.scast(f128, y)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_hypotl.c
            return hypot128(types.scast(f128, x), types.scast(f128, y));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_hypotf.c
//
// Original copyright notice:
// e_hypotf.c -- float version of e_hypot.c.
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
fn hypot32(x: f32, y: f32) f32 {
    var ha: i32 = @bitCast(x);
    ha &= 0x7fffffff;
    var hb: i32 = @bitCast(y);
    hb &= 0x7fffffff;

    var a: f32 = undefined;
    var b: f32 = undefined;
    if (hb > ha) {
        a = y;
        b = x;
        const j: i32 = ha;
        ha = hb;
        hb = j;
    } else {
        a = x;
        b = y;
    }

    a = float.abs(a);
    b = float.abs(b);
    if (ha -% hb > 0xf000000) // x/y > 2**30
        return a + b;

    var k: i32 = 0;
    if (ha > 0x58800000) { // a > 2**50
        if (ha >= 0x7f800000) { // Inf or NaN
            // Use original arg order iff result is NaN; quieten sNaNs
            var w: f32 = float.abs(x + 0.0) - float.abs(y + 0.0);

            if (ha == 0x7f800000)
                w = a;

            if (hb == 0x7f800000)
                w = b;

            return w;
        }

        // Scale a and b by 2**-68
        ha -%= 0x22000000;
        hb -%= 0x22000000;
        k +%= 68;

        a = @bitCast(ha);
        b = @bitCast(hb);
    }

    if (hb < 0x26800000) { // b < 2**-50
        if (hb <= 0x007fffff) { // subnormal b or 0
            if (hb == 0)
                return a;

            const t1: f32 = @bitCast(@as(u32, 0x7e800000)); // t1 = 2**126
            b *= t1;
            a *= t1;
            k -%= 126;
        } else { // Scale a and b by 2**68
            ha +%= 0x22000000; // a *= 2**68
            hb +%= 0x22000000; // b *= 2**68
            k -%= 68;
            a = @bitCast(ha);
            b = @bitCast(hb);
        }
    }

    // Medium size a and b
    var w: f32 = a - b;
    if (w > b) {
        const t1: f32 = @bitCast(@as(u32, @bitCast(ha)) & 0xfffff000);
        const t2: f32 = a - t1;
        w = float.sqrt(t1 * t1 - (b * (-b) - t2 * (a + t1)));
    } else {
        a = a + a;
        const y1: f32 = @bitCast(@as(u32, @bitCast(hb)) & 0xfffff000);
        const y2: f32 = b - y1;
        const t1: f32 = @bitCast(@as(u32, @bitCast(ha +% 0x00800000)) & 0xfffff000);
        const t2: f32 = a - t1;
        w = float.sqrt(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)));
    }
    if (k != 0) {
        const t1: f32 = @bitCast(0x3f800000 +% (k << 23));
        return t1 * w;
    } else {
        return w;
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_hypot.c
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
fn hypot64(x: f64, y: f64) f64 {
    var ha: i32 = @bitCast(dbl64.getHighPart(x));
    ha &= 0x7fffffff;
    var hb: i32 = @bitCast(dbl64.getHighPart(y));
    hb &= 0x7fffffff;

    var a: f64 = undefined;
    var b: f64 = undefined;
    if (hb > ha) {
        a = y;
        b = x;
        const j: i32 = ha;
        ha = hb;
        hb = j;
    } else {
        a = x;
        b = y;
    }

    a = float.abs(a);
    b = float.abs(b);
    if ((ha - hb) > 0x3c00000) // x/y > 2**60
        return a + b;

    var k: i32 = 0;
    if (ha > 0x5f300000) { // a > 2**500
        if (ha >= 0x7ff00000) { // Inf or NaN
            // Use original arg order iff result is NaN; quieten sNaNs
            var w: f64 = float.abs(x + 0.0) - float.abs(y + 0.0);
            var low: u32 = dbl64.getLowPart(a);

            if (((@as(u32, @bitCast(ha)) & 0xfffff) | low) == 0)
                w = a;

            low = dbl64.getLowPart(b);

            if (((@as(u32, @bitCast(hb)) ^ 0x7ff00000) | low) == 0)
                w = b;

            return w;
        }

        // Scale a and b by 2**-600
        ha -%= 0x25800000;
        hb -%= 0x25800000;
        k +%= 600;
        dbl64.setHighPart(&a, @bitCast(ha));
        dbl64.setHighPart(&b, @bitCast(hb));
    }

    if (hb < 0x20b00000) { // b < 2**-500
        if (hb <= 0x000fffff) { // Subnormal b or 0
            const low: u32 = dbl64.getLowPart(b);
            if ((@as(u32, @bitCast(hb)) | low) == 0)
                return a;

            var t1: f64 = 0;
            dbl64.setHighPart(&t1, 0x7fd00000); // t1 = 2**1022
            b *= t1;
            a *= t1;
            k -%= 1022;
        } else { // Scale a and b by 2**600
            ha +%= 0x25800000; // a *= 2**600
            hb +%= 0x25800000; // b *= 2**600
            k -%= 600;
            dbl64.setHighPart(&a, @bitCast(ha));
            dbl64.setHighPart(&b, @bitCast(hb));
        }
    }

    // Medium size a and b
    var w: f64 = a - b;
    if (w > b) {
        var t1: f64 = 0;
        dbl64.setHighPart(&t1, @bitCast(ha));
        const t2: f64 = a - t1;
        w = float.sqrt(t1 * t1 - (b * (-b) - t2 * (a + t1)));
    } else {
        a = a + a;
        var y1: f64 = 0;
        dbl64.setHighPart(&y1, @bitCast(hb));
        const y2: f64 = b - y1;
        var t1: f64 = 0;
        dbl64.setHighPart(&t1, @bitCast(ha +% 0x00100000));
        const t2: f64 = a - t1;
        w = float.sqrt(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)));
    }

    if (k != 0) {
        var t1: f64 = 1.0;
        const high: i32 = @bitCast(dbl64.getHighPart(t1));
        dbl64.setHighPart(&t1, @bitCast(high +% (k << 20)));
        return t1 * w;
    } else {
        return w;
    }
}

fn hypot80(x: f80, y: f80) f80 {
    _ = x;
    _ = y;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_hypotl.c
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
fn hypot128(x: f128, y: f128) f128 {
    var ha: i64 = @bitCast(ldbl128.getHighPart(x));
    ha &= 0x7fffffffffffffff;
    var hb: i64 = @bitCast(ldbl128.getHighPart(y));
    hb &= 0x7fffffffffffffff;

    var a: f128 = undefined;
    var b: f128 = undefined;
    if (hb > ha) {
        a = y;
        b = x;
        const j: i64 = ha;
        ha = hb;
        hb = j;
    } else {
        a = x;
        b = y;
    }

    ldbl128.setHighPart(&a, @bitCast(ha)); // abs(a)
    ldbl128.setHighPart(&b, @bitCast(hb)); // abs(b)

    if ((ha - hb) > 0x78000000000000) // x/y > 2**120
        return a + b;

    var k: i64 = 0;
    if (ha > 0x5f3f000000000000) { // a > 2**8000
        if (ha >= 0x7fff000000000000) { // Inf or NaN
            var w: f128 = a + b; // For sNaN
            var low: u64 = ldbl128.getLowPart(a);

            if (((@as(u64, @bitCast(ha)) & 0xffffffffffff) | low) == 0)
                w = a;

            low = ldbl128.getLowPart(b);

            if (((@as(u64, @bitCast(hb)) ^ 0x7fff000000000000) | low) == 0)
                w = b;

            return w;
        }

        // Scale a and b by 2**-9600
        ha -%= 0x2580000000000000;
        hb -%= 0x2580000000000000;
        k +%= 9600;

        ldbl128.setHighPart(&a, @bitCast(ha));
        ldbl128.setHighPart(&b, @bitCast(hb));
    }

    if (hb < 0x20bf000000000000) { // b < 2**-8000
        if (hb <= 0x0000ffffffffffff) { // Subnormal b or 0
            const low: u64 = ldbl128.getLowPart(b);
            if ((@as(u64, @bitCast(hb)) | low) == 0)
                return a;

            var t1: f128 = 0;
            ldbl128.setHighPart(&t1, 0x7ffd000000000000); // t1 = 2**16382
            b *= t1;
            a *= t1;
            k -%= 16382;
        } else { // Scale a and b by 2**9600
            ha +%= 0x2580000000000000; // a *= 2**9600
            hb +%= 0x2580000000000000; // b *= 2**9600
            k -%= 9600;
            ldbl128.setHighPart(&a, @bitCast(ha));
            ldbl128.setHighPart(&b, @bitCast(hb));
        }
    }

    // Medium size a and b
    var w: f128 = a - b;
    if (w > b) {
        var t1: f128 = 0;
        ldbl128.setHighPart(&t1, @bitCast(ha));
        const t2: f128 = a - t1;
        w = float.sqrt(t1 * t1 - (b * (-b) - t2 * (a + t1)));
    } else {
        a = a + a;
        var yy1: f128 = 0;
        ldbl128.setHighPart(&yy1, @bitCast(hb));
        const y2: f128 = b - yy1;
        var t1: f128 = 0;
        ldbl128.setHighPart(&t1, @bitCast(ha +% 0x0001000000000000));
        const t2: f128 = a - t1;
        w = float.sqrt(t1 * yy1 - (w * (-w) - (t1 * y2 + t2 * b)));
    }

    if (k != 0) {
        var t1: f128 = 1.0;
        const high: i64 = @bitCast(ldbl128.getHighPart(t1));
        ldbl128.setHighPart(&t1, @bitCast(high +% (k << 48)));
        return t1 * w;
    } else {
        return w;
    }
}
