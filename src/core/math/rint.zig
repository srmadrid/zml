const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const flt32 = @import("flt32.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn rint(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return rint(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, rint32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_rintf.c
                    return rint32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_rint.c
                    return rint64(x);
                },
                f80 => return cast(f80, rint128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_rintl.c
                    return rint128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn rint32(x: f32) f32 {
    // Use generic implementation.
    const TWO23: [2]f32 = .{
        8.3886080000e+06, // 0x4b000000
        -8.3886080000e+06, // 0xcb000000
    };

    var I0: i32 = undefined;
    flt32.getWord(&I0, x);
    const sx: i32 = (I0 >> 31) & 1;
    const j0: i32 = ((I0 >> 23) & 0xff) - 0x7f;
    if (j0 < 23) {
        if (j0 < 0) {
            const w: f32 = TWO23[@intCast(sx)] + x;
            var t = w - TWO23[@intCast(sx)];
            flt32.getWord(&I0, t);
            flt32.setWord(&t, (I0 & 0x7fffffff) | (sx << 31));
            return t;
        }
    } else {
        if (j0 == 0x80) {
            return x + x; // inf or NaN
        } else {
            return x; // x is integral
        }
    }
    const w: f32 = TWO23[@intCast(sx)] + x;
    return w - TWO23[@intCast(sx)];
}

fn rint64(x: f64) f64 {
    const TWO52: [2]f64 = .{
        4.50359962737049600000e+15, // 0x43300000, 0x00000000
        -4.50359962737049600000e+15, // 0xC3300000, 0x00000000
    };

    var I0: i64 = undefined;
    dbl64.extractWords64(&I0, x);
    const sx: i64 = (I0 >> 63) & 1;
    const j0: i32 = @truncate(((I0 >> 52) & 0x7ff) - 0x3ff);
    if (j0 < 52) {
        if (j0 < 0) {
            const w: f64 = TWO52[@intCast(sx)] + x;
            var t: f64 = w - TWO52[@intCast(sx)];
            dbl64.extractWords64(&I0, t);
            dbl64.insertWords64(&t, (I0 & 0x7fffffffffffffff) | (sx << 63));
            return t;
        }
    } else {
        if (j0 == 0x400) {
            return x + x; // inf or NaN
        } else {
            return x; // x is integral
        }
    }

    const w: f64 = TWO52[@intCast(sx)] + x;
    return w - TWO52[@intCast(sx)];
}

fn rint128(x: f128) f128 {
    const TWO112: [2]f128 = .{
        5.19229685853482762853049632922009600e+33, // 0x406F000000000000, 0
        -5.19229685853482762853049632922009600e+33, // 0xC06F000000000000, 0
    };

    var I0: i64 = undefined;
    var I1: u64 = undefined;
    ldbl128.getWords(&I0, &I1, x);
    var sx: i64 = undefined;
    {
        @setRuntimeSafety(false);
        sx = @bitCast((@as(u64, @intCast(I0))) >> 63);
    }
    const j0: i64 = ((I0 >> 48) & 0x7fff) - 0x3fff;
    if (j0 < 112) {
        if (j0 < 0) {
            const w: f128 = TWO112[@intCast(sx)] + x;
            var t: f128 = w - TWO112[@intCast(sx)];
            ldbl128.getMsw(&I0, t);
            ldbl128.setMsw(&t, (I0 & 0x7fffffffffffffff) | (sx << 63));
            return t;
        }
    } else {
        if (j0 == 0x4000) {
            return x + x; // inf or NaN
        } else {
            return x; // x is integral
        }
    }

    const w: f128 = TWO112[@intCast(sx)] + x;
    return w - TWO112[@intCast(sx)];
}
