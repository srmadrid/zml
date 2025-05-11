const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn hypot(x: anytype, y: anytype) EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    comptime if (!types.isFixedPrecision(@TypeOf(y)) or types.isComplex(@TypeOf(y)))
        @compileError("y must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            switch (types.numericType(@TypeOf(y))) {
                .int => return hypot(cast(EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))), x, .{}), cast(EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))), y, .{})),
                .float => return hypot(cast(Coerce(@TypeOf(x), @TypeOf(y)), x, .{}), cast(Coerce(@TypeOf(x), @TypeOf(y)), y, .{})),
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(y))) {
                .int => return hypot(cast(Coerce(@TypeOf(x), @TypeOf(y)), x, .{}), cast(Coerce(@TypeOf(x), @TypeOf(y)), y, .{})),
                .float => switch (Coerce(@TypeOf(x), @TypeOf(y))) {
                    f16 => return cast(f16, hypot32(cast(f32, x, .{}), cast(f32, y, .{})), .{}),
                    f32 => {
                        // glibc/sysdeps/ieee754/flt-32/e_hypotf.c
                        return hypot32(cast(f32, x, .{}), cast(f32, y, .{}));
                    },
                    f64 => {
                        // glibc/sysdeps/ieee754/dbl-64/e_hypot.c
                        return hypot64(cast(f64, x, .{}), cast(f64, y, .{}));
                    },
                    f80 => return cast(f80, hypot128(cast(f128, x, .{}), cast(f128, y, .{})), .{}),
                    f128 => {
                        // glibc/sysdeps/ieee754/ldbl-128/e_hypotl.c
                        return hypot128(cast(f128, x, .{}), cast(f128, y, .{}));
                    },
                    else => unreachable,
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn hypot32(x: f32, y: f32) f32 {
    if (!std.math.isFinite(x) or !std.math.isFinite(y)) {
        if ((std.math.isInf(x) or std.math.isInf(y)) and !std.math.isSignalNan(x) and !std.math.isSignalNan(y))
            return std.math.inf(f32);

        return x + y;
    }

    return cast(f32, float.sqrt(cast(f64, x, .{}) * cast(f64, x, .{}) + cast(f64, y, .{}) * cast(f64, y, .{})), .{});
}

inline fn kernel64(ax: f64, ay: f64) f64 {
    if (false) {
        const t1: f64 = ay + ay;
        const t2: f64 = ax - ay;

        if (t1 >= ax) {
            return float.sqrt(@mulAdd(f64, t1, ax, t2 * t2));
        } else {
            return float.sqrt(@mulAdd(f64, ax, ax, ay * ay));
        }
    } else {
        var h: f64 = float.sqrt(ax * ax + ay * ay);
        var t1: f64 = undefined;
        var t2: f64 = undefined;
        if (h <= 2.0 * ay) {
            const delta: f64 = h - ay;
            t1 = ax * (2.0 * delta - ax);
            t2 = (delta - 2.0 * (ax - ay)) * delta;
        } else {
            const delta: f64 = h - ax;
            t1 = 2.0 * delta * (ax - 2.0 * ay);
            t2 = (4.0 * delta - ay) * ay + delta * delta;
        }

        h -= (t1 + t2) / (2.0 * h);
        return h;
    }
}

fn hypot64(x: f64, y: f64) f64 {
    const SCALE: f64 = 0x1p-600;
    const LARGE_VAL: f64 = 0x1p+511;
    const TINY_VAL: f64 = 0x1p-459;
    const EPS: f64 = 0x1p-54;

    if (!std.math.isFinite(x) or !std.math.isFinite(y)) {
        if ((std.math.isInf(x) or std.math.isInf(y)) and !std.math.isSignalNan(x) and !std.math.isSignalNan(y))
            return std.math.inf(f64);

        return x + y;
    }

    const xx: f64 = @abs(x);
    const yy: f64 = @abs(y);

    var ax: f64 = @max(xx, yy);
    const ay: f64 = @min(xx, yy);

    // If ax is huge, scale both inputs down.
    if (ax > LARGE_VAL) {
        if (ay <= ax * EPS)
            return ax + ay;

        return kernel64(ax * SCALE, ay * SCALE) / SCALE;
    }

    // If ay is tiny, scale both inputs up.
    if (ay < TINY_VAL) {
        if (ax >= ay / EPS)
            return ax + ay;

        ax = kernel64(ax / SCALE, ay / SCALE) * SCALE;
        if (ax < std.math.floatMin(f64)) {
            const vax: f64 = ax * ax;
            std.mem.doNotOptimizeAway(vax);
        }
        return ax;
    }

    // Common case: ax is not huge and ay is not tiny.
    if (ay <= ax * EPS)
        return ax + ay;

    return kernel64(ax, ay);
}

// hypot128 kernel. The inputs must be adjusted so that ax >= ay >= 0
// and squaring ax, ay and (ax - ay) does not overflow or underflow.
inline fn kernel128(ax: f128, ay: f128) f128 {
    var h: f128 = float.sqrt(ax * ax + ay * ay);
    var t1: f128 = undefined;
    var t2: f128 = undefined;
    if (h <= 2 * ay) {
        const delta: f128 = h - ay;
        t1 = ax * (2 * delta - ax);
        t2 = (delta - 2 * (ax - ay)) * delta;
    } else {
        const delta: f128 = h - ax;
        t1 = 2 * delta * (ax - 2 * ay);
        t2 = (4 * delta - ay) * ay + delta * delta;
    }

    h -= (t1 + t2) / (2 * h);
    return h;
}

fn hypot128(x: f128, y: f128) f128 {
    const SCALE: f128 = 0x1p-8303;
    const LARGE_VAL: f128 = 0x1.6a09e667f3bcc908b2fb1366ea95p+8191;
    const TINY_VAL: f128 = 0x1p-8191;
    const EPS: f128 = 0x1p-114;

    if (!std.math.isFinite(x) or !std.math.isFinite(y)) {
        if ((std.math.isInf(x) or std.math.isInf(y)) and !std.math.isSignalNan(x) and !std.math.isSignalNan(y))
            return std.math.inf(f128);

        return x + y;
    }

    const xx: f128 = @abs(x);
    const yy: f128 = @abs(y);

    var ax: f128 = @max(xx, yy);
    const ay: f128 = @min(xx, yy);

    // If ax is huge, scale both inputs down.
    if (ax > LARGE_VAL) {
        if (ay <= ax * EPS)
            return ax + ay;

        return kernel128(ax * SCALE, ay * SCALE) / SCALE;
    }

    // If ay is tiny, scale both inputs up.
    if (ay < TINY_VAL) {
        if (ax >= ay / EPS)
            return ax + ay;

        ax = kernel128(ax / SCALE, ay / SCALE) * SCALE;
        if (ax < std.math.floatMin(f128)) {
            const vax: f128 = ax * ax;
            std.mem.doNotOptimizeAway(vax);
        }
        return ax;
    }

    // Common case: ax is not huge and ay is not tiny.
    if (ay <= ax * EPS)
        return ax + ay;

    return kernel128(ax, ay);
}
