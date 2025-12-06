const std = @import("std");
const types = @import("../types.zig");
const flt32 = @import("flt32.zig");
const dbl64 = @import("dbl64.zig");
const ldbl80 = @import("ldbl80.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn frexp(x: anytype, e: *i32) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.frexp: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (@TypeOf(x)) {
        f16 => return scast(f16, frexp32(scast(f32, x), e)),
        f32 => {
            // glibc/sysdeps/ieee754/flt-32/s_frexpf.c
            return frexp32(scast(f32, x), e);
        },
        f64 => {
            // glibc/sysdeps/ieee754/dbl-64/s_frexp.c
            return frexp64(scast(f64, x), e);
        },
        f80 => {
            // glibc/sysdeps/ieee754/ldbl-96/s_frexpl.c
            return frexp80(scast(f80, x), e);
        },
        f128 => {
            // glibc/sysdeps/ieee754/ldbl-128/s_frexpl.c
            return frexp128(scast(f128, x), e);
        },
        else => unreachable,
    }
}

fn frexp32(x: f32, e: *i32) f32 {
    const two25: f32 = 3.3554432000e+07; // 0x4c000000

    var hx: i32 = undefined;
    flt32.getWord(&hx, x);
    var ix: i32 = 0x7fffffff & hx;
    e.* = 0;

    if (ix >= 0x7f800000 || (ix == 0))
        return x + x; // 0,inf,nan

    var xx: f32 = x;
    if (ix < 0x00800000) { // subnormal
        xx *= two25;
        flt32.getWord(&hx, xx);
        ix = hx & 0x7fffffff;
        e.* = -25;
    }

    e.* += (ix >> 23) - 126;
    hx = (hx & 0x807fffff) | 0x3f000000;
    flt32.setWord(&xx, hx);
    return xx;
}

fn frexp64(x: f64, e: *i32) f64 {
    var ix: u64 = undefined;
    dbl64.extractWords64(&ix, x);
    var ex: i32 = @bitCast(@as(u32, @truncate(0x7ff & (ix >> 52))));
    var ee: i32 = 0;

    var xx: f64 = x;
    if (ex != 0x7ff and x != 0.0) {
        @branchHint(.likely);
        // Not zero and finite.
        ee = ex - 1022;
        if (ex == 0) {
            @branchHint(.unlikely);
            // Subnormal.
            xx *= 0x1p54;
            dbl64.extractWords64(&ix, xx);
            ex = @bitCast(@as(u32, @truncate(0x7ff & (ix >> 52))));
            ee = ex - 1022 - 54;
        }

        ix = (ix & 0x800fffffffffffff) | 0x3fe0000000000000;
        dbl64.insertWords64(&xx, ix);
    } else {
        // Quiet signaling NaNs.
        xx += xx;
    }

    e.* = ee;
    return xx;
}

fn frexp80(x: f80, e: *i32) f80 {
    const two65: f80 = 3.68934881474191032320e+19;

    var se: u32 = undefined;
    var hx: u32 = undefined;
    var lx: u32 = undefined;
    ldbl80.getWords(&se, &hx, &lx, x);
    var ix: u32 = 0x7fff & se;
    e.* = 0;

    if (ix == 0x7fff or ((ix | hx | lx) == 0))
        return x + x; // 0,inf,nan

    var xx: f80 = x;
    if (ix == 0x0000) { // subnormal
        xx *= two65;
        ldbl80.getExp(&se, xx);
        ix = se & 0x7fff;
        e.* = -65;
    }

    e.* += scast(i32, ix) -% 16382;
    se = (se & 0x8000) | 0x3ffe;
    ldbl80.setExp(&xx, se);
    return xx;
}

fn frexp128(x: f128, e: *i32) f128 {
    const two114: f128 = 2.0769187434139310514121985316880384e+34; // 0x4071000000000000, 0

    var hx: u64 = @truncate(@as(u128, @bitCast(x)) >> 64);
    const lx: u64 = @truncate(@as(u128, @bitCast(x)) & 0xffffffffffffffff);
    // ldbl128.getWords(&hx, &lx, x);
    var ix: u64 = 0x7fffffffffffffff & hx;
    e.* = 0;
    if (ix >= 0x7fff000000000000 or (ix | lx) == 0)
        return x + x; // 0,inf,nan

    var xx: f128 = x;
    if (ix < 0x0001000000000000) { // subnormal
        xx *= two114;
        // ldbl128.getMsw(&hx, xx);
        hx = @truncate(@as(u128, @bitCast(xx)) >> 64);
        ix = hx & 0x7fffffffffffffff;
        e.* = -114;
    }

    e.* += scast(i32, ix >> 48) -% 16382;
    hx = (hx & 0x8000ffffffffffff) | 0x3ffe000000000000;
    // ldbl128.setMsw(&xx, hx);
    xx = @bitCast(@as(u128, (types.scast(u128, hx) << 64) | (@as(u128, @bitCast(xx)) & 0xffffffffffffffff)));

    return xx;
}
