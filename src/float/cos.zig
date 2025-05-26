const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const sincos_data = @import("sincos_data.zig");
const atnat = @import("atnat.zig");
const usncs = @import("usncs.zig");
const sin = @import("sin.zig");
const branred = @import("branred.zig");
const ldbl128 = @import("ldbl128.zig");
const rem_pio2 = @import("rem_pio2.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn cos(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.cos: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return scast(f16, cos32(scast(f32, x))),
        f32 => {
            // glibc/sysdeps/ieee754/flt-32/s_cosf.c
            return cos32(scast(f32, x));
        },
        f64 => {
            // glibc/sysdeps/ieee754/dbl-64/s_sin.c
            return cos64(scast(f64, x));
        },
        f80 => return scast(f80, cos128(scast(f128, x))),
        f128 => {
            // glibc/sysdeps/ieee754/ldbl-128/s_cosl.c
            return cos128(scast(f128, x));
        },
        else => unreachable,
    }
}

// Fast cosf implementation.  Worst-case ULP is 0.5607, maximum relative
// error is 0.5303 * 2^-23.  A single-step range reduction is used for
// small values.  Large inputs have their range reduced using fast integer
// arithmetic.
fn cos32(y: f32) f32 {
    var x: f64 = scast(f64, y);
    var p: *const sincos_data.Sincos32 = &sincos_data.__sincos32_table[0];

    if (sincos_data.abstop12(y) < sincos_data.abstop12(sincos_data.pio4)) {
        const x2: f64 = x * x;

        if (sincos_data.abstop12(y) < sincos_data.abstop12(0x1p-12)) {
            @branchHint(.unlikely);
            return 1.0;
        }

        return sincos_data.sin32_poly(x, x2, p, 1);
    } else if (sincos_data.abstop12(y) < sincos_data.abstop12(120.0)) {
        @branchHint(.likely);
        var n: i32 = undefined;
        x = sincos_data.reduce_fast(x, p, &n);

        // Setup the signs for sin and cos.
        const s: f64 = p.sign[@intCast(n & 3)];

        if ((n & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        return sincos_data.sin32_poly(x * s, x * x, p, n ^ 1);
    } else if (sincos_data.abstop12(y) < sincos_data.abstop12(std.math.inf(f32))) {
        const xi: u32 = @bitCast(y);
        const sign: i32 = @intCast(xi >> 31);

        var n: i32 = undefined;
        x = sincos_data.reduce_large(xi, &n);

        // Setup signs for sin and cos - include original sign.
        const s: f64 = p.sign[@intCast((n + sign) & 3)];

        if (((n + sign) & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        return sincos_data.sin32_poly(x * s, x * x, p, n ^ 1);
    } else return (y - y) / (y - y);
}

// Given a number partitioned into X and DX, this function computes the cosine
// of the number by combining the sin and cos of X (as computed by a variation
// of the Taylor series) with the values looked up from the sin/cos table to
// get the result.
pub inline fn do_cos64(x: f64, dx: f64) f64 {
    var dxnew: f64 = dx;
    if (x < 0)
        dxnew = -dxnew;

    const u: [2]i32 = @bitCast(usncs.big + float.abs(x));
    const xnew: f64 = float.abs(x) - (@as(f64, @bitCast(u)) - usncs.big) + dxnew;

    const xx: f64 = xnew * xnew;
    const s: f64 = xnew + xnew * xx * (sincos_data.sn3 + xx * sincos_data.sn5);
    const c: f64 = xx * (sincos_data.cs2 + xx * (sincos_data.cs4 + xx * sincos_data.cs6));
    const k: i32 = u[atnat.LOW_HALF] << 2;
    const sn: f64 = @bitCast(sincos_data.__sincostab[@intCast(k)]);
    const ssn: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 1)]);
    const cs: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 2)]);
    const ccs: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 3)]);
    const cor: f64 = (ccs - s * ssn - cs * c) - sn * s;
    return cs + cor;
}

//******************************************************************/
// An ultimate cos routine. Given an IEEE double machine number x
// it computes the rounded value of cos(x).
//******************************************************************/
fn cos64(x: f64) f64 {
    const u: [2]i32 = @bitCast(x);
    const m: i32 = u[atnat.HIGH_HALF];
    const k: i32 = 0x7fffffff & m;

    // |x|<2^-27 => cos(x)=1
    if (k < 0x3e400000) {
        return 1.0;
    } else if (k < 0x3feb6000) { // 2^-27 < |x| < 0.855469
        // Max ULP is 0.51.
        return do_cos64(x, 0);
    } else if (k < 0x400368fd) { // 0.855469  <|x|<2.426265
        const y: f64 = usncs.hp0 - float.abs(x);
        const a: f64 = y + usncs.hp1;
        const da: f64 = (y - a) + usncs.hp1;
        // Max ULP is 0.501 if xx < 0.01588 or 0.518 otherwise.
        // Range reduction uses 106 bits here which is sufficient.
        return sin.do_sin64(a, da);
    } else if (k < 0x419921fb) { // 2.426265<|x|< 105414350
        var a: f64 = undefined;
        var da: f64 = undefined;
        const n: i32 = sincos_data.reduce_sincos64(x, &a, &da);
        return sincos_data.do_sincos64(a, da, n + 1);
    } else if (k < 0x7ff00000) { // 105414350 <|x| <2^1024
        var a: f64 = undefined;
        var da: f64 = undefined;
        const n: i32 = branred.branred(x, &a, &da);
        return sincos_data.do_sincos64(a, da, n + 1);
    } else {
        return x / x; // |x| > 2^1024
    }
}

pub fn kernel_cos128(x: f128, y: f128) f128 {
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);
    var tix: u32 = @truncate(@as(u64, @bitCast(ix)) >> 32);
    tix &= ~@as(u32, 0x80000000); // tix = |x|'s high 32 bits
    if (tix < 0x3ffc3000) // |x| < 0.1484375
    {
        // Argument is small enough to approximate it by a Chebyshev
        // polynomial of degree 16.
        if (tix < 0x3fc60000) { // |x| < 2^-57
            if (scast(i32, x) == 1)
                return 1; // generate inexact
        }

        const z: f128 = x * x;
        return 1 + (z * (sincos_data.COS1 + z * (sincos_data.COS2 + z * (sincos_data.COS3 + z * (sincos_data.COS4 +
            z * (sincos_data.COS5 + z * (sincos_data.COS6 + z * (sincos_data.COS7 + z * sincos_data.COS8))))))));
    } else {
        // So that we don't have to use too large polynomial,  we find
        // l and h such that x = l + h,  where fabsl(l) <= 1.0/256 with 83
        // possible values for h.  We look up cosl(h) and sinl(h) in
        // pre-computed tables,  compute cosl(l) and sinl(l) using a
        // Chebyshev polynomial of degree 10(11) and compute
        // cosl(h+l) = cosl(h)cosl(l) - sinl(h)sinl(l).
        var index: u32 = 0x3ffe - (tix >> 16);
        const hix: u32 = (tix + (@as(u32, 0x200) << @as(u5, @intCast(index)))) & (@as(u32, 0xfffffc00) << @as(u5, @intCast(index)));
        var xx: f128 = x;
        var yy: f128 = y;
        if (std.math.signbit(x)) {
            xx = -x;
            yy = -y;
        }
        switch (index) {
            0 => index = ((45 << 10) + hix - 0x3ffe0000) >> 8,
            1 => index = ((13 << 11) + hix - 0x3ffd0000) >> 9,
            else => index = (hix - 0x3ffc3000) >> 10,
        }

        var h: f128 = undefined;
        ldbl128.setWords(&h, scast(u64, hix) << 32, @as(u64, 0));
        const l: f128 = yy - (h - xx);
        const z: f128 = l * l;
        const sin_l: f128 = l * (1 + z * (sincos_data.SSIN1 + z * (sincos_data.SSIN2 + z * (sincos_data.SSIN3 + z * (sincos_data.SSIN4 + z * sincos_data.SSIN5)))));
        const cos_l_m1: f128 = z * (sincos_data.SCOS1 + z * (sincos_data.SCOS2 + z * (sincos_data.SCOS3 + z * (sincos_data.SCOS4 + z * sincos_data.SCOS5))));
        return usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_LO] - (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] * sin_l - usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] * cos_l_m1));
    }
}

fn cos128(x: f128) f128 {
    // High word of x.
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);

    // |x| ~< pi/4
    ix &= 0x7fffffffffffffff;
    if (ix <= 0x3ffe921fb54442d1) {
        return kernel_cos128(x, 0);
    } else if (ix >= 0x7fff000000000000) { // cos(Inf or NaN) is NaN
        return x - x;
    } else { // argument reduction needed
        var y: [2]f128 = undefined;
        const n = rem_pio2.rem_pio2_128(x, &y);
        switch (n & 3) {
            0 => return kernel_cos128(y[0], y[1]),
            1 => return -sin.kernel_sin128(y[0], y[1], 1),
            2 => return -kernel_cos128(y[0], y[1]),
            else => return sin.kernel_sin128(y[0], y[1], 1),
        }
    }
}
