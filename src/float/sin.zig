const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const sincos_data = @import("sincos_data.zig");
const atnat = @import("atnat.zig");
const usncs = @import("usncs.zig");
const cos = @import("cos.zig");
const branred = @import("branred.zig");
const ldbl128 = @import("ldbl128.zig");
const rem_pio2 = @import("rem_pio2.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn sin(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return sin(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, sin32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_sinf.c
                    return sin32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_sin.c
                    return sin64(x);
                },
                f80 => return cast(f80, sin128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_sinl.c
                    return sin128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

// Fast sinf implementation.  Worst-case ULP is 0.5607, maximum relative
// error is 0.5303 * 2^-23.  A single-step range reduction is used for
// small values.  Large inputs have their range reduced using fast integer
// arithmetic.
fn sin32(y: f32) f32 {
    var x: f64 = cast(f64, y, .{});
    var p: *const sincos_data.Sincos32 = &sincos_data.__sincos32_table[0];

    if (sincos_data.abstop12(y) < sincos_data.abstop12(sincos_data.pio4)) {
        const s: f64 = x * x;

        if (sincos_data.abstop12(y) < sincos_data.abstop12(0x1p-12)) {
            @branchHint(.unlikely);
            // Force underflow for tiny y.
            if (sincos_data.abstop12(y) < sincos_data.abstop12(0x1p-126)) {
                @branchHint(.unlikely);
                std.mem.doNotOptimizeAway(cast(f32, s, .{}));
            }

            return y;
        }

        return sincos_data.sin32_poly(x, s, p, 0);
    } else if (sincos_data.abstop12(y) < sincos_data.abstop12(120.0)) {
        @branchHint(.likely);
        var n: i32 = undefined;
        x = sincos_data.reduce_fast(x, p, &n);

        // Setup the signs for sin and cos.
        const s: f64 = p.sign[@intCast(n & 3)];

        if ((n & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        return sincos_data.sin32_poly(x * s, x * x, p, n);
    } else if (sincos_data.abstop12(y) < sincos_data.abstop12(std.math.inf(f32))) {
        const xi: u32 = @bitCast(y);
        const sign: i32 = @intCast(xi >> 31);

        var n: i32 = undefined;
        x = sincos_data.reduce_large(xi, &n);

        // Setup signs for sin and cos - include original sign.
        const s: f64 = p.sign[@intCast((n + sign) & 3)];

        if (((n + sign) & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        return sincos_data.sin32_poly(x * s, x * x, p, n);
    } else {
        return (y - y) / (y - y);
    }
}

// Given a number partitioned into X and DX, this function computes the sine of
// the number by combining the sin and cos of X (as computed by a variation of
// the Taylor series) with the values looked up from the sin/cos table to get
// the result.
pub inline fn do_sin64(x: f64, dx: f64) f64 {
    // Max ULP is 0.501 if |x| < 0.126, otherwise ULP is 0.518.
    if (float.abs(x) < 0.126) {
        const t: f64 = ((((((usncs.s5 * x * x + usncs.s4) * x * x + usncs.s3) * x * x + usncs.s2) * x * x) + usncs.s1) * x - 0.5 * dx) * x * x + dx;
        return x + t;
    }

    var dxnew: f64 = dx;
    if (x <= 0)
        dxnew = -dxnew;

    const u: [2]i32 = @bitCast(usncs.big + float.abs(x));
    var xnew: f64 = x;
    xnew = float.abs(xnew) - (@as(f64, @bitCast(u)) - usncs.big);

    const xx: f64 = xnew * xnew;
    const s: f64 = xnew + (dxnew + xnew * xx * (sincos_data.sn3 + xx * sincos_data.sn5));
    const c: f64 = xnew * dxnew + xx * (sincos_data.cs2 + xx * (sincos_data.cs4 + xx * sincos_data.cs6));
    const k: i32 = u[atnat.LOW_HALF] << 2;
    const sn: f64 = @bitCast(sincos_data.__sincostab[@intCast(k)]);
    const ssn: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 1)]);
    const cs: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 2)]);
    const ccs: f64 = @bitCast(sincos_data.__sincostab[@intCast(k + 3)]);
    const cor: f64 = (ssn + s * ccs - sn * c) + cs * s;
    return float.copysign(sn + cor, x);
}

//******************************************************************
// An ultimate sin routine. Given an IEEE double machine number x
// it computes the rounded value of sin(x).
//******************************************************************
fn sin64(x: f64) f64 {
    const u: [2]i32 = @bitCast(x);
    const m: i32 = u[atnat.HIGH_HALF];
    const k: i32 = 0x7fffffff & m; // no sign
    if (k < 0x3e500000) { // if x->0 =>sin(x)=x
        if (float.abs(x) < std.math.floatMin(f64)) {
            const vx: f64 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }
        return x;
    }
    //--------------------------- 2^-26<|x|< 0.855469----------------------
    else if (k < 0x3feb6000) {
        // Max ULP is 0.548.
        return do_sin64(x, 0);
    }
    //----------------------- 0.855469  <|x|<2.426265  ----------------------
    else if (k < 0x400368fd) {
        const t: f64 = usncs.hp0 - float.abs(x);
        // Max ULP is 0.51.
        return float.copysign(cos.do_cos64(t, usncs.hp1), x);
    }
    //-------------------------- 2.426265<|x|< 105414350 ----------------------
    else if (k < 0x419921fb) {
        var a: f64 = undefined;
        var da: f64 = undefined;
        const n: i32 = sincos_data.reduce_sincos64(x, &a, &da);
        return sincos_data.do_sincos64(a, da, n);
    }
    //--------------------105414350 <|x| <2^1024------------------------------
    else if (k < 0x7ff00000) {
        var a: f64 = undefined;
        var da: f64 = undefined;
        const n: i32 = branred.branred(x, &a, &da);
        return sincos_data.do_sincos64(a, da, n);
    }
    //--------------------- |x| > 2^1024 ----------------------------------
    else {
        return x / x;
    }
}

pub fn kernel_sin128(x: f128, y: f128, iy: i32) f128 {
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);
    var tix: u32 = @truncate(@as(u64, @bitCast(ix)) >> 32);
    tix &= ~@as(u32, 0x80000000); // tix = |x|'s high 32 bits
    if (tix < 0x3ffc3000) { // |x| < 0.1484375
        // Argument is small enough to approximate it by a Chebyshev
        // polynomial of degree 17.
        if (tix < 0x3fc60000) { // |x| < 2^-57
            if (float.abs(x) < std.math.floatMin(f128)) {
                const vx: f128 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            if (cast(i32, x, .{}) == 0)
                return x; // generate inexact
        }
        const z: f128 = x * x;
        return x + (x * (z * (sincos_data.SIN1 + z * (sincos_data.SIN2 + z * (sincos_data.SIN3 + z * (sincos_data.SIN4 +
            z * (sincos_data.SIN5 + z * (sincos_data.SIN6 + z * (sincos_data.SIN7 + z * sincos_data.SIN8)))))))));
    } else {
        // So that we don't have to use too large polynomial,  we find
        // l and h such that x = l + h,  where fabsl(l) <= 1.0/256 with 83
        // possible values for h.  We look up cosl(h) and sinl(h) in
        // pre-computed tables,  compute cosl(l) and sinl(l) using a
        // Chebyshev polynomial of degree 10(11) and compute
        // sinl(h+l) = sinl(h)cosl(l) + cosl(h)sinl(l).
        var index: u32 = 0x3ffe - (tix >> 16);
        const hix: u32 = (tix + (@as(u32, 0x200) << @as(u5, @intCast(index)))) & (@as(u32, 0xfffffc00) << @as(u5, @intCast(index)));
        const xx: f128 = float.abs(x);
        switch (index) {
            0 => index = ((45 << 10) + hix - 0x3ffe0000) >> 8,
            1 => index = ((13 << 11) + hix - 0x3ffd0000) >> 9,
            else => index = (hix - 0x3ffc3000) >> 10,
        }

        var h: f128 = undefined;
        ldbl128.setWords(&h, cast(u64, hix, .{}) << 32, @as(u64, 0));
        var l: f128 = undefined;
        if (iy != 0) {
            l = (if (ix < 0) -y else y) - (h - xx);
        } else {
            l = xx - h;
        }

        var z: f128 = l * l;
        const sin_l: f128 = l * (1 + z * (sincos_data.SSIN1 + z * (sincos_data.SSIN2 + z * (sincos_data.SSIN3 + z * (sincos_data.SSIN4 + z * sincos_data.SSIN5)))));
        const cos_l_m1: f128 = z * (sincos_data.SCOS1 + z * (sincos_data.SCOS2 + z * (sincos_data.SCOS3 + z * (sincos_data.SCOS4 + z * sincos_data.SCOS5))));
        z = usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_LO] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] * cos_l_m1) + (usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] * sin_l));
        return if (ix < 0) -z else z;
    }
}

fn sin128(x: f128) f128 {
    // High word of x.
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);

    // |x| ~< pi/4
    ix &= 0x7fffffffffffffff;
    if (ix <= 0x3ffe921fb54442d1) {
        return kernel_sin128(x, 0, 0);
    } else if (ix >= 0x7fff000000000000) { // sin(Inf or NaN) is NaN
        return x - x;
    } else { // argument reduction needed
        var y: [2]f128 = undefined;
        const n: i32 = rem_pio2.rem_pio2_128(x, &y);
        switch (n & 3) {
            0 => return kernel_sin128(y[0], y[1], 1),
            1 => return cos.kernel_cos128(y[0], y[1]),
            2 => return -kernel_sin128(y[0], y[1], 1),
            else => return -cos.kernel_cos128(y[0], y[1]),
        }
    }
}
