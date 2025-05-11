const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const sincos_data = @import("sincos_data.zig");
const atnat = @import("atnat.zig");
const usncs = @import("usncs.zig");
const sin = @import("sin.zig");
const cos = @import("cos.zig");
const branred = @import("branred.zig");
const ldbl128 = @import("ldbl128.zig");
const rem_pio2 = @import("rem_pio2.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn sincos(x: anytype) struct { sinx: EnsureFloat(@TypeOf(x)), cosx: EnsureFloat(@TypeOf(x)) } {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return sincos(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => {
                    const res = sincos32(cast(f32, x));
                    return .{
                        .sinx = cast(f16, res.sinx, .{}),
                        .cosx = cast(f16, res.cosx, .{}),
                    };
                },
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_sincosf.c
                    const res = sincos32(x);
                    return .{
                        .sinx = res.sinx,
                        .cosx = res.cosx,
                    };
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_sincos.c
                    const res = sincos64(x);
                    return .{
                        .sinx = res.sinx,
                        .cosx = res.cosx,
                    };
                },
                f80 => {
                    const res = sincos128(cast(f128, x, .{}));
                    return .{
                        .sinx = cast(f80, res.sinx, .{}),
                        .cosx = cast(f80, res.cosx, .{}),
                    };
                },
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_sincosl.c
                    const res = sincos128(x);
                    return .{
                        .sinx = res.sinx,
                        .cosx = res.cosx,
                    };
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

// Fast sincosf implementation.  Worst-case ULP is 0.5607, maximum relative
// error is 0.5303 * 2^-23.  A single-step range reduction is used for
// small values.  Large inputs have their range reduced using fast integer
// arithmetic.
fn sincos32(x: f32) struct { sinx: f32, cosx: f32 } {
    var p: *const sincos_data.Sincos32 = &sincos_data.__sincos32_table[0];

    if (sincos_data.abstop12(x) < sincos_data.abstop12(sincos_data.pio4)) {
        const x2: f64 = x * x;

        if (sincos_data.abstop12(x) < sincos_data.abstop12(0x1p-12)) {
            @branchHint(.unlikely);
            // Force underflow for tiny y.
            if (sincos_data.abstop12(x) < sincos_data.abstop12(0x1p-126)) {
                @branchHint(.unlikely);
                std.mem.doNotOptimizeAway(cast(f32, x2, .{}));
            }

            return .{
                .sinx = x,
                .cosx = 1,
            };
        }

        var res: struct { sinx: f32, cosx: f32 } = undefined;
        sincos_data.sincos32_poly(x, x2, p, 0, &res.sinx, &res.cosx);
        return .{
            .sinx = res.sinx,
            .cosx = res.cosx,
        };
    } else if (sincos_data.abstop12(x) < sincos_data.abstop12(120.0)) {
        var n: i32 = undefined;
        const xx: f64 = sincos_data.reduce_fast(x, p, &n);

        // Setup the signs for sin and cos.
        const s: f64 = p.sign[@intCast(n & 3)];

        if ((n & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        var res: struct { sinx: f32, cosx: f32 } = undefined;
        sincos_data.sincos32_poly(xx * s, xx * xx, p, n, &res.sinx, &res.cosx);
        return .{
            .sinx = res.sinx,
            .cosx = res.cosx,
        };
    } else if (sincos_data.abstop12(x) < sincos_data.abstop12(std.math.inf(f32))) {
        @branchHint(.likely);
        const xi: u32 = @bitCast(x);
        const sign: i32 = @bitCast(xi >> 31);

        var n: i32 = undefined;
        const xx: f64 = sincos_data.reduce_large(xi, &n);

        // Setup signs for sin and cos - include original sign.
        const s: f64 = p.sign[@intCast((n + sign) & 3)];

        if (((n + sign) & 2) != 0)
            p = &sincos_data.__sincos32_table[1];

        var res: struct { sinx: f32, cosx: f32 } = undefined;
        sincos_data.sincos32_poly(xx * s, xx * xx, p, n, &res.sinx, &res.cosx);
        return .{
            .sinx = res.sinx,
            .cosx = res.cosx,
        };
    } else {
        // Return NaN if Inf or NaN for both sin and cos.
        return .{
            .sinx = x - x,
            .cosx = x - x,
        };
    }
}

// Reduce range of x to within PI/2 with abs (x) < 105414350.  The high part
// is written to *a, the low part to *da.  Range reduction is accurate to 136
// bits so that when x is large and *a very close to zero, all 53 bits of *a
// are correct.
inline fn reduce_sincos64(x: f64, a: *f64, da: *f64) i32 {
    const t: f64 = (x * usncs.hpinv + usncs.toint);
    const xn: f64 = t - usncs.toint;
    const v: [2]i32 = @bitCast(t);
    const y: f64 = (x - xn * usncs.mp1) - xn * usncs.mp2;
    const n: i32 = v[atnat.LOW_HALF] & 3;

    var t1: f64 = xn * usncs.pp3;
    const t2: f64 = y - t1;
    var db: f64 = (y - t2) - t1;

    t1 = xn * usncs.pp4;
    const b: f64 = t2 - t1;
    db += (t2 - b) - t1;

    a.* = b;
    da.* = db;
    return n;
}

fn sincos64(x: f64) struct { sinx: f64, cosx: f64 } {
    const u: [2]i32 = @bitCast(x);
    const k: i32 = u[atnat.HIGH_HALF] & 0x7fffffff;

    if (k < 0x400368fd) {
        // |x| < 2^-27 => cos (x) = 1, sin (x) = x.
        if (k < 0x3e400000) {
            if (k < 0x3e500000) {
                if (float.abs(x) < std.math.floatMin(f64)) {
                    const vx: f64 = x * x;
                    std.mem.doNotOptimizeAway(vx);
                }
            }

            return .{
                .sinx = x,
                .cosx = 1,
            };
        }
        // |x| < 0.855469.
        else if (k < 0x3feb6000) {
            return .{
                .sinx = sin.do_sin64(x, 0),
                .cosx = cos.do_cos64(x, 0),
            };
        }

        // |x| < 2.426265.
        const y: f64 = usncs.hp0 - float.abs(x);
        const a: f64 = y + usncs.hp1;
        const da: f64 = (y - a) + usncs.hp1;

        return .{
            .sinx = float.copysign(cos.do_cos64(a, da), x),
            .cosx = sin.do_sin64(a, da),
        };
    }
    // |x| < 2^1024.
    if (k < 0x7ff00000) {

        // If |x| < 105414350 use simple range reduction.
        var a: f64 = undefined;
        var da: f64 = undefined;
        var n: i32 = if (k < 0x419921fb) reduce_sincos64(x, &a, &da) else branred.branred(x, &a, &da);
        n = n & 3;

        if (n == 1 or n == 2) {
            a = -a;
            da = -da;
        }

        if ((n & 1) != 0) {
            const xx: f64 = cos.do_cos64(a, da);
            return .{
                .sinx = if ((n & 2) != 0) -xx else xx,
                .cosx = sin.do_sin64(a, da),
            };
        } else {
            const xx: f64 = cos.do_cos64(a, da);
            return .{
                .sinx = sin.do_sin64(a, da),
                .cosx = if ((n & 2) != 0) -xx else xx,
            };
        }
    }

    return .{
        .sinx = x / x,
        .cosx = x / x,
    };
}

fn kernel_sincos128(x: f128, y: f128, sinx: *f128, cosx: *f128, iy: i32) void {
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);
    var tix: u32 = @truncate((@as(u64, @bitCast(ix))) >> 32);
    tix &= ~@as(u32, 0x80000000); // tix = |x|'s high 32 bits
    if (tix < 0x3ffc3000) // |x| < 0.1484375
    {
        // Argument is small enough to approximate it by a Chebyshev
        // polynomial of degree 16(17).
        if (tix < 0x3fc60000) // |x| < 2^-57
        {
            if (float.abs(x) < std.math.floatMin(f128)) {
                const vx: f128 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            if (x == 0) { // generate inexact
                sinx.* = x;
                cosx.* = 1;
            }
        }
        const z: f128 = x * x;
        sinx.* = x + (x * (z * (sincos_data.SIN1 + z * (sincos_data.SIN2 + z * (sincos_data.SIN3 + z * (sincos_data.SIN4 +
            z * (sincos_data.SIN5 + z * (sincos_data.SIN6 + z * (sincos_data.SIN7 + z * sincos_data.SIN8)))))))));
        cosx.* = 1 + (z * (sincos_data.COS1 + z * (sincos_data.COS2 + z * (sincos_data.COS3 + z * (sincos_data.COS4 +
            z * (sincos_data.COS5 + z * (sincos_data.COS6 + z * (sincos_data.COS7 + z * sincos_data.COS8))))))));
    } else {
        // So that we don't have to use too large polynomial,  we find
        // l and h such that x = l + h,  where fabsl(l) <= 1.0/256 with 83
        // possible values for h.  We look up cosl(h) and sinl(h) in
        // pre-computed tables,  compute cosl(l) and sinl(l) using a
        // Chebyshev polynomial of degree 10(11) and compute
        // sinl(h+l) = sinl(h)cosl(l) + cosl(h)sinl(l) and
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
        ldbl128.setWords(&h, cast(u64, hix, .{}) << 32, @as(u64, 0));
        var l: f128 = undefined;
        if (iy != 0) {
            l = yy - (h - xx);
        } else {
            l = xx - h;
        }

        var z: f128 = l * l;
        const sin_l: f128 = l * (1 + z * (sincos_data.SSIN1 + z * (sincos_data.SSIN2 + z * (sincos_data.SSIN3 + z * (sincos_data.SSIN4 + z * sincos_data.SSIN5)))));
        const cos_l_m1: f128 = z * (sincos_data.SCOS1 + z * (sincos_data.SCOS2 + z * (sincos_data.SCOS3 + z * (sincos_data.SCOS4 + z * sincos_data.SCOS5))));
        z = usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_LO] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] * cos_l_m1) + (usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] * sin_l));
        sinx.* = if (ix < 0) -z else z;
        cosx.* = usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] + (usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_LO] - (usncs.__sincos128_table[index + sincos_data.SINCOSL_SIN_HI] * sin_l - usncs.__sincos128_table[index + sincos_data.SINCOSL_COS_HI] * cos_l_m1));
    }
}

fn sincos128(x: f128) struct { sinx: f128, cosx: f128 } {
    // High word of x.
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);

    // |x| ~< pi/4
    ix &= 0x7fffffffffffffff;
    if (ix <= 0x3ffe921fb54442d1) {
        var res: struct { sinx: f128, cosx: f128 } = undefined;
        kernel_sincos128(x, 0, &res.sinx, &res.cosx, 0);
        return .{
            .sinx = res.sinx,
            .cosx = res.cosx,
        };
    } else if (ix >= 0x7fff000000000000) {
        // sin(Inf or NaN) is NaN
        return .{
            .sinx = x - x,
            .cosx = x - x,
        };
    } else {
        // Argument reduction needed.
        var y: [2]f128 = undefined;
        const n: i32 = rem_pio2.rem_pio2_128(x, &y);
        switch (n & 3) {
            0 => {
                var res: struct { sinx: f128, cosx: f128 } = undefined;
                kernel_sincos128(y[0], y[1], &res.sinx, &res.cosx, 1);
                return .{
                    .sinx = res.sinx,
                    .cosx = res.cosx,
                };
            },
            1 => {
                var res: struct { sinx: f128, cosx: f128 } = undefined;
                kernel_sincos128(y[0], y[1], &res.cosx, &res.sinx, 1);
                return .{
                    .sinx = res.sinx,
                    .cosx = -res.cosx,
                };
            },
            2 => {
                var res: struct { sinx: f128, cosx: f128 } = undefined;
                kernel_sincos128(y[0], y[1], &res.sinx, &res.cosx, 1);
                return .{
                    .sinx = -res.sinx,
                    .cosx = -res.cosx,
                };
            },
            else => {
                var res: struct { sinx: f128, cosx: f128 } = undefined;
                kernel_sincos128(y[0], y[1], &res.cosx, &res.sinx, 1);
                return .{
                    .sinx = -res.sinx,
                    .cosx = res.cosx,
                };
            },
        }
    }
}
