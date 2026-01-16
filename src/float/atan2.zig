const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Atan2(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float))
        @compileError("zml.float.atan2: x and y must be bools, ints or floats, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.EnsureFloat(types.Coerce(X, Y));
}

pub inline fn atan2(x: anytype, y: anytype) Atan2(@TypeOf(x), @TypeOf(y)) {
    switch (Atan2(@TypeOf(x), @TypeOf(y))) {
        f16 => return types.scast(f16, atan2_32(types.scast(f32, x), types.scast(f32, y))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_atan2f.c
            return atan2_32(types.scast(f32, x), types.scast(f32, y));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_atan2.c
            return atan2_64(types.scast(f64, x), types.scast(f64, y));
        },
        f80 => {
            //
            // return atan280(types.scast(f80, x), types.scast(f80, y));
            return types.scast(f80, atan2_128(types.scast(f128, x), types.scast(f128, y)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_atan2l.c
            return atan2_128(types.scast(f128, x), types.scast(f128, y));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_logf.c
//
// Original copyright notice:
// e_atan2f.c -- float version of e_atan2.c.
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
fn atan2_32(y: f32, x: f32) f32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;
    const hy: i32 = @bitCast(y);
    const iy: i32 = hy & 0x7fffffff;

    if (ix > 0x7f800000 or iy > 0x7f800000) // x or y is NaN
        return x + y;

    if (hx == 0x3f800000)
        return float.atan(y); // x = 1.0

    var m: i32 = ((hy >> 31) & 1) | ((hx >> 30) & 2); // 2 * sign(x) + sign(y)

    if (iy == 0) { // y = 0
        return switch (m) {
            0, 1 => y, // atan(±0, +anything) = ±0
            2 => 3.1415927410e+0, // atan(+0, -anything) = pi
            3 => -3.1415927410e+0, // atan(-0, -anything) = -pi
            else => unreachable,
        };
    }

    if (ix == 0) // x = 0
        return if (hy < 0)
            -1.5707963705e+0
        else
            1.5707963705e+0;

    if (ix == 0x7f800000) { // x is Inf
        if (iy == 0x7f800000) {
            return switch (m) {
                0 => 7.8539818525e-1, // atan(+Inf, +Inf)
                1 => -7.8539818525e-1, // atan(-Inf, +Inf)
                2 => 3.0 * 7.8539818525e-1, // atan(+Inf, -Inf)
                3 => -3.0 * 7.8539818525e-1, // atan(-Inf, -Inf)
                else => unreachable,
            };
        } else {
            return switch (m) {
                0 => 0.0, // atan(+..., +Inf)
                1 => -0.0, // atan(-..., +Inf)
                2 => 3.1415927410e+0, // atan(+..., -Inf)
                3 => -3.1415927410e+0, // atan(-..., -Inf)
                else => unreachable,
            };
        }
    }

    if (iy == 0x7f800000) // y is Inf
        return if (hy < 0)
            -1.5707963705e+0
        else
            1.5707963705e+0;

    // Compute y/x
    const k: i32 = (iy - ix) >> 23;
    var z: f32 = undefined;
    if (k > 26) { // |y/x| > 2**26
        z = 1.5707963705e+0 + 0.5 * -8.7422776573e-8;
        m &= 1;
    } else if (k < -26 and hx < 0) {
        z = 0.0; // 0 > |y|/x > -2**-26
    } else {
        z = float.atan(float.abs(y / x)); // Safe to do y/x
    }

    return switch (m) {
        0 => z, // atan(+, +)
        1 => -z, // atan(-, +)
        2 => 3.1415927410e+0 - (z - -8.7422776573e-8), // atan(+, -)
        else => (z - -8.7422776573e-8) - 3.1415927410e+0, // atan(-, -)
    };
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_atan2.c
//
// Original copyright notice:
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunSoft, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn atan2_64(y: f64, x: f64) f64 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const lx: u32 = dbl64.getLowPart(x);
    const ix: i32 = hx & 0x7fffffff;
    const hy: i32 = @bitCast(dbl64.getHighPart(y));
    const ly: u32 = dbl64.getLowPart(y);
    const iy: i32 = hy & 0x7fffffff;

    if ((@as(u32, @bitCast(ix)) | ((lx | (0 -% lx)) >> 31)) > 0x7ff00000 or (@as(u32, @bitCast(iy)) | ((ly | (0 -% ly)) >> 31)) > 0x7ff00000) // x or y is NaN
        return x + y;

    if (hx == 0x3ff00000 and lx == 0)
        return float.atan(y); // x = 1.0

    var m: i32 = ((hy >> 31) & 1) | ((hx >> 30) & 2); // 2 * sign(x) + sign(y)

    if ((@as(u32, @bitCast(iy)) | ly) == 0) { // y = 0
        return switch (m) {
            0, 1 => y, // atan(±0, +anything)= ±0
            2 => 3.1415926535897931160e+0, // atan(+0, -anything) = pi
            3 => -3.1415926535897931160e+0, // atan(-0, -anything) = -pi
            else => unreachable,
        };
    }

    if ((@as(u32, @bitCast(ix)) | lx) == 0) // x = 0
        return if (hy < 0)
            -1.5707963267948965580e+0
        else
            1.5707963267948965580e+0;

    if (ix == 0x7ff00000) { // x is Inf
        if (iy == 0x7ff00000) {
            return switch (m) {
                0 => 7.8539816339744827900e-1, // atan(+Inf, +Inf)
                1 => -7.8539816339744827900e-1, // atan(-Inf, +Inf)
                2 => 3.0 * 7.8539816339744827900e-1, //atan(+Inf, -Inf)
                3 => -3.0 * 7.8539816339744827900e-1, //atan(-Inf, -Inf)
                else => unreachable,
            };
        } else {
            return switch (m) {
                0 => 0.0, // atan(+..., +Inf)
                1 => -0.0, // atan(-..., +Inf)
                2 => 3.1415926535897931160e+0, // atan(+..., -Inf)
                3 => -3.1415926535897931160e+0, // atan(-..., -Inf)
                else => unreachable,
            };
        }
    }

    if (iy == 0x7ff00000) // y is Inf
        return if (hy < 0)
            -1.5707963267948965580e+0
        else
            1.5707963267948965580e+0;

    // Compute y/x
    const k: i32 = (iy - ix) >> 20;
    var z: f64 = undefined;
    if (k > 60) { // |y/x| > 2**60
        z = 1.5707963267948965580e+0 + 0.5 * 1.2246467991473531772e-16;
        m &= 1;
    } else if (k < -60 and hx < 0) {
        z = 0.0; // 0 > |y|/x > -2**-60
    } else {
        z = float.atan(float.abs(y / x)); // Safe to do y/x
    }

    return switch (m) {
        0 => z, // atan(+, +)
        1 => -z, // atan(-, +)
        2 => 3.1415926535897931160e+0 - (z - 1.2246467991473531772e-16), // atan(+, -)
        else => (z - 1.2246467991473531772e-16) - 3.1415926535897931160e+0, // atan(-, -)
    };
}

fn atan2_80(x: f80, y: f80) f80 {
    _ = x;
    _ = y;
    return std.math.nan(f80);
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_atan2l.c
//
// Original copyright notice:
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunSoft, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
fn atan2_128(y: f128, x: f128) f128 {
    const ux: ldbl128.ShapeSplit = .fromFloat(x);
    const signx: u1 = ux.sign;
    const exptx: i16 = @intCast(ux.exponent);
    const expsignx: i16 = exptx | @as(i16, @bitCast(@as(u16, @intCast(signx)) << 15));
    const uy: ldbl128.ShapeSplit = .fromFloat(y);
    const signy: u1 = uy.sign;
    const expty: i16 = @intCast(uy.exponent);
    const expsigny: i16 = expty | @as(i16, @bitCast(@as(u16, @intCast(signy)) << 15));

    if ((exptx == (16384 - 1) + 16384 and
        (ux.mantissa_high | ux.mantissa_low) != 0) or // x is NaN
        (expty == (16384 - 1) + 16384 and
            (uy.mantissa_high | uy.mantissa_low) != 0)) // y is NaN
        return x + y;

    if (expsignx == (16384 - 1) and (ux.mantissa_high | ux.mantissa_low) == 0) // x = 1.0
        return float.atan(y);

    var m: i16 = ((expsigny >> 15) & 1) | ((expsignx >> 14) & 2); // 2 * sign(x) + sign(y)

    if (expty == 0 and (uy.mantissa_high | uy.mantissa_low) == 0) { // y = 0
        return switch (m) {
            0, 1 => y, // atan(±0, +anything)= ±0
            2 => 3.14159265358979323846264338327950280e+0, // atan(+0, -anything) = pi
            3 => -3.14159265358979323846264338327950280e+0, // atan(-0, -anything) = -pi
            else => unreachable,
        };
    }

    if (exptx == 0 and (ux.mantissa_high | ux.mantissa_low) == 0) // x = 0
        return if (signy != 0)
            -1.57079632679489661923132169163975140e+0
        else
            1.57079632679489661923132169163975140e+0;

    if (exptx == (16384 - 1) + 16384) { // x is Inf
        if (expty == (16384 - 1) + 16384) {
            return switch (m) {
                0 => 1.57079632679489661923132169163975140e+0 * 0.5, // atan(+Inf, +Inf)
                1 => -1.57079632679489661923132169163975140e+0 * 0.5, // atan(-Inf, +Inf)
                2 => 1.5 * 1.57079632679489661923132169163975140e+0, //atan(+Inf, -Inf)
                3 => -1.5 * 1.57079632679489661923132169163975140e+0, //atan(-Inf, -Inf)
                else => unreachable,
            };
        } else {
            return switch (m) {
                0 => 0.0, // atan(+..., +Inf)
                1 => -0.0, // atan(-..., +Inf)
                2 => 3.14159265358979323846264338327950280e+0, // atan(+..., -Inf)
                3 => -3.14159265358979323846264338327950280e+0, // atan(-..., -Inf)
                else => unreachable,
            };
        }
    }

    if (expty == (16384 - 1) + 16384) // y is Inf
        return if (signy != 0)
            -1.57079632679489661923132169163975140e+0
        else
            1.57079632679489661923132169163975140e+0;

    // Compute y/x
    const k: i16 = expty - exptx;
    var z: f128 = undefined;
    if (k > 114) { // |y/x| huge
        z = 1.57079632679489661923132169163975140e+0 + 4.33590506506189051239852201302167613e-35;
        m &= 1;
    } else if (k < -114 and signx != 0) {
        z = 0.0; // |y/x| tiny, x < 0
    } else {
        z = float.atan(float.abs(y / x)); // Safe to do y/x
    }

    return switch (m) {
        0 => z, // atan(+, +)
        1 => -z, // atan(-, +)
        2 => 3.14159265358979323846264338327950280e+0 - (z - 8.67181013012378102479704402604335225e-35), // atan(+, -)
        else => (z - 8.67181013012378102479704402604335225e-35) - 3.14159265358979323846264338327950280e+0, // atan(-, -)
    };
}
