const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
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
                if (math.abs(x) < std.math.floatMin(f64)) {
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
        const y: f64 = usncs.hp0 - math.abs(x);
        const a: f64 = y + usncs.hp1;
        const da: f64 = (y - a) + usncs.hp1;

        return .{
            .sinx = math.copysign(cos.do_cos64(a, da), x),
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
            if (math.abs(x) < std.math.floatMin(f128)) {
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

test sincos {
    try std.testing.expectEqual(0x0p+0, sincos(@as(f32, 0x0p+0)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f32, 0x0p+0)).cosx);
    try std.testing.expectEqual(-0x0p+0, sincos(@as(f32, -0x0p+0)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f32, -0x0p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f32, 0x1.921fb6p+0)).sinx);
    try std.testing.expectEqual(-0xb.bbd2ep-28, sincos(@as(f32, 0x1.921fb6p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f32, 0x1.921fb4p+0)).sinx);
    try std.testing.expectEqual(0x1.4442d2p-24, sincos(@as(f32, 0x1.921fb4p+0)).cosx);
    try std.testing.expectEqual(0x8p-4, sincos(@as(f32, 0x8.60a92p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d7p-4, sincos(@as(f32, 0x8.60a92p-4)).cosx);
    try std.testing.expectEqual(0x7.fffff8p-4, sincos(@as(f32, 0x8.60a91p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d8p-4, sincos(@as(f32, 0x8.60a91p-4)).cosx);
    try std.testing.expectEqual(0xd.db3d8p-4, sincos(@as(f32, 0x1.0c1524p+0)).sinx);
    try std.testing.expectEqual(0x7.fffff8p-4, sincos(@as(f32, 0x1.0c1524p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d7p-4, sincos(@as(f32, 0x1.0c1522p+0)).sinx);
    try std.testing.expectEqual(0x8.00001p-4, sincos(@as(f32, 0x1.0c1522p+0)).cosx);
    try std.testing.expectEqual(-0x1.777a5cp-24, sincos(@as(f32, 0x3.243f6cp+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f32, 0x3.243f6cp+0)).cosx);
    try std.testing.expectEqual(0x2.8885a4p-24, sincos(@as(f32, 0x3.243f68p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f32, 0x3.243f68p+0)).cosx);
    try std.testing.expectEqual(0x1.777a5cp-24, sincos(@as(f32, -0x3.243f6cp+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f32, -0x3.243f6cp+0)).cosx);
    try std.testing.expectEqual(-0x2.8885a4p-24, sincos(@as(f32, -0x3.243f68p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f32, -0x3.243f68p+0)).cosx);
    try std.testing.expectEqual(0xa.e7fe1p-4, sincos(@as(f32, 0xcp-4)).sinx);
    try std.testing.expectEqual(0xb.b4ff6p-4, sincos(@as(f32, 0xcp-4)).cosx);
    try std.testing.expectEqual(-0xc.143e1p-8, sincos(@as(f32, 0x2p+64)).sinx);
    try std.testing.expectEqual(0xf.fb702p-4, sincos(@as(f32, 0x2p+64)).cosx);
    try std.testing.expectEqual(0xc.143e1p-8, sincos(@as(f32, -0x2p+64)).sinx);
    try std.testing.expectEqual(0xf.fb702p-4, sincos(@as(f32, -0x2p+64)).cosx);
    try std.testing.expectEqual(0xb.7fb6p-4, sincos(@as(f32, 0xc.d4967p-4)).sinx);
    try std.testing.expectEqual(0xb.201e7p-4, sincos(@as(f32, 0xc.d4967p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fp-4, sincos(@as(f32, 0xc.d4966p-4)).sinx);
    try std.testing.expectEqual(0xb.201e8p-4, sincos(@as(f32, 0xc.d4966p-4)).cosx);
    try std.testing.expectEqual(-0x4.cd7e88p-4, sincos(@as(f32, 0x2.1e19e4p+72)).sinx);
    try std.testing.expectEqual(0xf.431ddp-4, sincos(@as(f32, 0x2.1e19e4p+72)).cosx);
    try std.testing.expectEqual(-0xb.becc4p-4, sincos(@as(f32, 0x2.1e19ep+72)).sinx);
    try std.testing.expectEqual(0xa.dd6f7p-4, sincos(@as(f32, 0x2.1e19ep+72)).cosx);
    try std.testing.expectEqual(-0x8.599b3p-4, sincos(@as(f32, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f96p-4, sincos(@as(f32, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(-0x8.599b3p-4, sincos(@as(f32, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f96p-4, sincos(@as(f32, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x6.0b8d18p-4, sincos(@as(f32, 0x1p+120)).sinx);
    try std.testing.expectEqual(-0xe.d0668p-4, sincos(@as(f32, 0x1p+120)).cosx);
    try std.testing.expectEqual(0x9.f9631p-4, sincos(@as(f32, 0x8p+124)).sinx);
    try std.testing.expectEqual(0xc.82b8fp-4, sincos(@as(f32, 0x8p+124)).cosx);
    try std.testing.expectEqual(0xc.6fa5cp-8, sincos(@as(f32, 0xf.ffffcp+124)).sinx);
    try std.testing.expectEqual(0xf.fb2ap-4, sincos(@as(f32, 0xf.ffffcp+124)).cosx);
    try std.testing.expectEqual(-0x8.599b3p-4, sincos(@as(f32, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f96p-4, sincos(@as(f32, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x7.f13d78p-4, sincos(@as(f32, 0x4p+48)).sinx);
    try std.testing.expectEqual(0xd.e3b89p-4, sincos(@as(f32, 0x4p+48)).cosx);
    try std.testing.expectEqual(-0xf.c777cp-4, sincos(@as(f32, 0x1p+28)).sinx);
    try std.testing.expectEqual(-0x2.a62ba8p-4, sincos(@as(f32, 0x1p+28)).cosx);
    try std.testing.expectEqual(0x8.599b3p-4, sincos(@as(f32, -0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f96p-4, sincos(@as(f32, -0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0xc.773a3p-4, sincos(@as(f32, 0xe.9f1e6p+112)).sinx);
    try std.testing.expectEqual(-0xa.07bd4p-4, sincos(@as(f32, 0xe.9f1e6p+112)).cosx);
    try std.testing.expectEqual(0x7.76d6p-4, sincos(@as(f32, 0xe.9f1e5p+112)).sinx);
    try std.testing.expectEqual(0xe.26f8bp-4, sincos(@as(f32, 0xe.9f1e5p+112)).cosx);
    try std.testing.expectEqual(-0x1.ffb67ap-4, sincos(@as(f32, 0x4.7857dp+68)).sinx);
    try std.testing.expectEqual(-0xf.dfe9p-4, sincos(@as(f32, 0x4.7857dp+68)).cosx);
    try std.testing.expectEqual(-0x1.fecbp-4, sincos(@as(f32, 0x6.287cdp+0)).sinx);
    try std.testing.expectEqual(0xf.e006ap-4, sincos(@as(f32, 0x6.287cdp+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb7ep-4, sincos(@as(f32, 0x6.287cc8p+0)).sinx);
    try std.testing.expectEqual(0xf.e0069p-4, sincos(@as(f32, 0x6.287cc8p+0)).cosx);
    try std.testing.expectEqual(-0xd.8f692p-4, sincos(@as(f32, -0x1.02e34cp+0)).sinx);
    try std.testing.expectEqual(0x8.7e0eap-4, sincos(@as(f32, -0x1.02e34cp+0)).cosx);
    try std.testing.expectEqual(-0x8.3beep-4, sincos(@as(f32, 0xf.f0274p+4)).sinx);
    try std.testing.expectEqual(-0xd.b7f53p-4, sincos(@as(f32, 0xf.f0274p+4)).cosx);
    try std.testing.expectEqual(0x1.ffc6dap-4, sincos(@as(f32, 0x3.042d88p+0)).sinx);
    try std.testing.expectEqual(-0xf.dfe6fp-4, sincos(@as(f32, 0x3.042d88p+0)).cosx);
    try std.testing.expectEqual(-0x8.599b3p-4, sincos(@as(f32, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f96p-4, sincos(@as(f32, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x8.599b3p-4, sincos(@as(f32, -0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f96p-4, sincos(@as(f32, -0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x4p-128, sincos(@as(f32, 0x4p-128)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f32, 0x4p-128)).cosx);
    try std.testing.expectEqual(-0x4p-128, sincos(@as(f32, -0x4p-128)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f32, -0x4p-128)).cosx);
    try std.testing.expectEqual(0x8p-152, sincos(@as(f32, 0x8p-152)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f32, 0x8p-152)).cosx);
    try std.testing.expectEqual(-0x8p-152, sincos(@as(f32, -0x8p-152)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f32, -0x8p-152)).cosx);
    try std.testing.expectEqual(0xf.fa2aep-4, sincos(@as(f32, 0x1.8475e6p+0)).sinx);
    try std.testing.expectEqual(0xd.a8263p-8, sincos(@as(f32, 0x1.8475e6p+0)).cosx);
    try std.testing.expectEqual(0xf.fa2aep-4, sincos(@as(f32, 0x1.8475e4p+0)).sinx);
    try std.testing.expectEqual(0xd.a8283p-8, sincos(@as(f32, 0x1.8475e4p+0)).cosx);

    try std.testing.expectEqual(0x0p+0, sincos(@as(f64, 0x0p+0)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, 0x0p+0)).cosx);
    try std.testing.expectEqual(-0x0p+0, sincos(@as(f64, -0x0p+0)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, -0x0p+0)).cosx);
    try std.testing.expectEqual(0xf.fffffffffffb8p-4, sincos(@as(f64, 0x1.921fb6p+0)).sinx);
    try std.testing.expectEqual(-0xb.bbd2e7b96766p-28, sincos(@as(f64, 0x1.921fb6p+0)).cosx);
    try std.testing.expectEqual(0xf.fffffffffff3p-4, sincos(@as(f64, 0x1.921fb4p+0)).sinx);
    try std.testing.expectEqual(0x1.4442d18469893p-24, sincos(@as(f64, 0x1.921fb4p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, 0x1.921fb54442d19p+0)).sinx);
    // try std.testing.expectEqual(-0xb.9676733ae8fe8p-56, sincos(@as(f64, 0x1.921fb54442d19p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, 0x1.921fb54442d18p+0)).sinx);
    try std.testing.expectEqual(0x4.69898cc51701cp-56, sincos(@as(f64, 0x1.921fb54442d18p+0)).cosx);
    try std.testing.expectEqual(0x8.0000036321168p-4, sincos(@as(f64, 0x8.60a92p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d7237832ep-4, sincos(@as(f64, 0x8.60a92p-4)).cosx);
    try std.testing.expectEqual(0x7.fffff587e3a04p-4, sincos(@as(f64, 0x8.60a91p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d7a37832a8p-4, sincos(@as(f64, 0x8.60a91p-4)).cosx);
    try std.testing.expectEqual(0x8p-4, sincos(@as(f64, 0x8.60a91c16b9b3p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c2655p-4, sincos(@as(f64, 0x8.60a91c16b9b3p-4)).cosx);
    try std.testing.expectEqual(0x7.ffffffffffffcp-4, sincos(@as(f64, 0x8.60a91c16b9b28p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c26558p-4, sincos(@as(f64, 0x8.60a91c16b9b28p-4)).cosx);
    try std.testing.expectEqual(0xd.db3d78156ca1p-4, sincos(@as(f64, 0x1.0c1524p+0)).sinx);
    try std.testing.expectEqual(0x7.fffff939bdd18p-4, sincos(@as(f64, 0x1.0c1524p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d68156c928p-4, sincos(@as(f64, 0x1.0c1522p+0)).sinx);
    try std.testing.expectEqual(0x8.000014f038b18p-4, sincos(@as(f64, 0x1.0c1522p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c26558p-4, sincos(@as(f64, 0x1.0c152382d7366p+0)).sinx);
    try std.testing.expectEqual(0x7.ffffffffffff8p-4, sincos(@as(f64, 0x1.0c152382d7366p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c2655p-4, sincos(@as(f64, 0x1.0c152382d7365p+0)).sinx);
    try std.testing.expectEqual(0x8.0000000000008p-4, sincos(@as(f64, 0x1.0c152382d7365p+0)).cosx);
    try std.testing.expectEqual(-0x1.777a5cf72cec6p-24, sincos(@as(f64, 0x3.243f6cp+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffefp-4, sincos(@as(f64, 0x3.243f6cp+0)).cosx);
    try std.testing.expectEqual(0x2.8885a308d3106p-24, sincos(@as(f64, 0x3.243f68p+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffcc8p-4, sincos(@as(f64, 0x3.243f68p+0)).cosx);
    try std.testing.expectEqual(-0x1.72cece675d1fdp-52, sincos(@as(f64, 0x3.243f6a8885a32p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f64, 0x3.243f6a8885a32p+0)).cosx);
    try std.testing.expectEqual(0x8.d313198a2e038p-56, sincos(@as(f64, 0x3.243f6a8885a3p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f64, 0x3.243f6a8885a3p+0)).cosx);
    try std.testing.expectEqual(0x1.777a5cf72cec6p-24, sincos(@as(f64, -0x3.243f6cp+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffefp-4, sincos(@as(f64, -0x3.243f6cp+0)).cosx);
    try std.testing.expectEqual(-0x2.8885a308d3106p-24, sincos(@as(f64, -0x3.243f68p+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffcc8p-4, sincos(@as(f64, -0x3.243f68p+0)).cosx);
    try std.testing.expectEqual(0x1.72cece675d1fdp-52, sincos(@as(f64, -0x3.243f6a8885a32p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f64, -0x3.243f6a8885a32p+0)).cosx);
    try std.testing.expectEqual(-0x8.d313198a2e038p-56, sincos(@as(f64, -0x3.243f6a8885a3p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f64, -0x3.243f6a8885a3p+0)).cosx);
    try std.testing.expectEqual(0xa.e7fe0b5fc7868p-4, sincos(@as(f64, 0xcp-4)).sinx);
    try std.testing.expectEqual(0xb.b4ff632a908f8p-4, sincos(@as(f64, 0xcp-4)).cosx);
    try std.testing.expectEqual(-0xc.143e153b0702p-8, sincos(@as(f64, 0x2p+64)).sinx);
    try std.testing.expectEqual(0xf.fb701e22987f8p-4, sincos(@as(f64, 0x2p+64)).cosx);
    try std.testing.expectEqual(0xc.143e153b0702p-8, sincos(@as(f64, -0x2p+64)).sinx);
    try std.testing.expectEqual(0xf.fb701e22987f8p-4, sincos(@as(f64, -0x2p+64)).cosx);
    try std.testing.expectEqual(0xb.7fb6002758778p-4, sincos(@as(f64, 0xc.d4967p-4)).sinx);
    try std.testing.expectEqual(0xb.201e77869a468p-4, sincos(@as(f64, 0xc.d4967p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5f50739fa8p-4, sincos(@as(f64, 0xc.d4966p-4)).sinx);
    try std.testing.expectEqual(0xb.201e830650418p-4, sincos(@as(f64, 0xc.d4966p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776979p-4, sincos(@as(f64, 0xc.d4966d92d171p-4)).sinx);
    try std.testing.expectEqual(0xb.201e794508848p-4, sincos(@as(f64, 0xc.d4966d92d171p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776979p-4, sincos(@as(f64, 0xc.d4966d92d1708p-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450885p-4, sincos(@as(f64, 0xc.d4966d92d1708p-4)).cosx);
    try std.testing.expectEqual(-0x4.cd7e86c4077cp-4, sincos(@as(f64, 0x2.1e19e4p+72)).sinx);
    try std.testing.expectEqual(0xf.431dd7a36cf38p-4, sincos(@as(f64, 0x2.1e19e4p+72)).cosx);
    try std.testing.expectEqual(-0xb.becc47ab1b8c8p-4, sincos(@as(f64, 0x2.1e19ep+72)).sinx);
    try std.testing.expectEqual(0xa.dd6f6bacd2068p-4, sincos(@as(f64, 0x2.1e19ep+72)).cosx);
    try std.testing.expectEqual(-0xd.a29d5bb5f9cb8p-4, sincos(@as(f64, 0x2.1e19e0c9bab24p+72)).sinx);
    try std.testing.expectEqual(0x8.5f167780e47ap-4, sincos(@as(f64, 0x2.1e19e0c9bab24p+72)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba8p-4, sincos(@as(f64, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, sincos(@as(f64, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x9.0292465edbbp-4, sincos(@as(f64, 0x8p+1020)).sinx);
    try std.testing.expectEqual(-0xd.38cf9361195f8p-4, sincos(@as(f64, 0x8p+1020)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba8p-4, sincos(@as(f64, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, sincos(@as(f64, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x1.452fc98b34e97p-8, sincos(@as(f64, 0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba8p-4, sincos(@as(f64, 0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(0x6.0b8d19579bf2cp-4, sincos(@as(f64, 0x1p+120)).sinx);
    try std.testing.expectEqual(-0xe.d06685b36c67p-4, sincos(@as(f64, 0x1p+120)).cosx);
    try std.testing.expectEqual(0x9.f963166f215e8p-4, sincos(@as(f64, 0x8p+124)).sinx);
    try std.testing.expectEqual(0xc.82b8ec98b5e6p-4, sincos(@as(f64, 0x8p+124)).cosx);
    try std.testing.expectEqual(0xc.6fa5c5665985p-8, sincos(@as(f64, 0xf.ffffcp+124)).sinx);
    try std.testing.expectEqual(0xf.fb2a030c5ae2p-4, sincos(@as(f64, 0xf.ffffcp+124)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba8p-4, sincos(@as(f64, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, sincos(@as(f64, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x7.f13d78eabb76cp-4, sincos(@as(f64, 0x4p+48)).sinx);
    try std.testing.expectEqual(0xd.e3b88804f0058p-4, sincos(@as(f64, 0x4p+48)).cosx);
    try std.testing.expectEqual(-0xf.c777c6b36a75p-4, sincos(@as(f64, 0x1p+28)).sinx);
    try std.testing.expectEqual(-0x2.a62ba8824e5bcp-4, sincos(@as(f64, 0x1p+28)).cosx);
    try std.testing.expectEqual(0x8.599b32844aba8p-4, sincos(@as(f64, -0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, sincos(@as(f64, -0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(-0xf.e00885042dd78p-4, sincos(@as(f64, -0x3.3de320f6be87ep+1020)).sinx);
    try std.testing.expectEqual(-0x1.febbf9949ecc1p-4, sincos(@as(f64, -0x3.3de320f6be87ep+1020)).cosx);
    try std.testing.expectEqual(0xc.773a2eac3001p-4, sincos(@as(f64, 0xe.9f1e6p+112)).sinx);
    try std.testing.expectEqual(-0xa.07bd3ab53ab98p-4, sincos(@as(f64, 0xe.9f1e6p+112)).cosx);
    try std.testing.expectEqual(0x7.76d600e03152p-4, sincos(@as(f64, 0xe.9f1e5p+112)).sinx);
    try std.testing.expectEqual(0xe.26f8af8333f9p-4, sincos(@as(f64, 0xe.9f1e5p+112)).cosx);
    try std.testing.expectEqual(0xf.dfffd7bde0fb8p-4, sincos(@as(f64, 0xe.9f1e5bc3bb88p+112)).sinx);
    try std.testing.expectEqual(0x1.ff01000c9ae73p-4, sincos(@as(f64, 0xe.9f1e5bc3bb88p+112)).cosx);
    try std.testing.expectEqual(-0x1.ffb679ba994b7p-4, sincos(@as(f64, 0x4.7857dp+68)).sinx);
    try std.testing.expectEqual(-0xf.dfe902135fc2p-4, sincos(@as(f64, 0x4.7857dp+68)).cosx);
    try std.testing.expectEqual(-0x1.fecaff6878a11p-4, sincos(@as(f64, 0x6.287cdp+0)).sinx);
    try std.testing.expectEqual(0xf.e006a1ad17db8p-4, sincos(@as(f64, 0x6.287cdp+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb7e68ad6eap-4, sincos(@as(f64, 0x6.287cc8p+0)).sinx);
    try std.testing.expectEqual(0xf.e00691b6bde4p-4, sincos(@as(f64, 0x6.287cc8p+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb772e1b83p-4, sincos(@as(f64, 0x6.287cc8749213p+0)).sinx);
    try std.testing.expectEqual(0xf.e006929f558d8p-4, sincos(@as(f64, 0x6.287cc8749213p+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb772e1b87p-4, sincos(@as(f64, 0x6.287cc8749212cp+0)).sinx);
    try std.testing.expectEqual(0xf.e006929f558dp-4, sincos(@as(f64, 0x6.287cc8749212cp+0)).cosx);
    try std.testing.expectEqual(-0xd.8f691a7a95428p-4, sincos(@as(f64, -0x1.02e34cp+0)).sinx);
    try std.testing.expectEqual(0x8.7e0ea4db2f488p-4, sincos(@as(f64, -0x1.02e34cp+0)).cosx);
    try std.testing.expectEqual(-0x8.3bee07bc90768p-4, sincos(@as(f64, 0xf.f0274p+4)).sinx);
    try std.testing.expectEqual(-0xd.b7f5359babdb8p-4, sincos(@as(f64, 0xf.f0274p+4)).cosx);
    try std.testing.expectEqual(0x1.ffc6da9f1ffeep-4, sincos(@as(f64, 0x3.042d88p+0)).sinx);
    try std.testing.expectEqual(-0xf.dfe6f2169e25p-4, sincos(@as(f64, 0x3.042d88p+0)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba8p-4, sincos(@as(f64, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, sincos(@as(f64, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x1.452fc98b34e97p-8, sincos(@as(f64, 0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba8p-4, sincos(@as(f64, 0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(0x8.599b32844aba8p-4, sincos(@as(f64, -0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe7p-4, sincos(@as(f64, -0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(-0x1.452fc98b34e97p-8, sincos(@as(f64, -0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba8p-4, sincos(@as(f64, -0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(0x4p-128, sincos(@as(f64, 0x4p-128)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, 0x4p-128)).cosx);
    try std.testing.expectEqual(0x4p-1024, sincos(@as(f64, 0x4p-1024)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, 0x4p-1024)).cosx);
    try std.testing.expectEqual(0x8p-972, sincos(@as(f64, 0x8p-972)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, 0x8p-972)).cosx);
    try std.testing.expectEqual(-0x4p-128, sincos(@as(f64, -0x4p-128)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, -0x4p-128)).cosx);
    try std.testing.expectEqual(-0x4p-1024, sincos(@as(f64, -0x4p-1024)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, -0x4p-1024)).cosx);
    try std.testing.expectEqual(-0x8p-972, sincos(@as(f64, -0x8p-972)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, -0x8p-972)).cosx);
    try std.testing.expectEqual(0x8p-152, sincos(@as(f64, 0x8p-152)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, 0x8p-152)).cosx);
    try std.testing.expectEqual(0x4p-1076, sincos(@as(f64, 0x4p-1076)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, 0x4p-1076)).cosx);
    try std.testing.expectEqual(-0x8p-152, sincos(@as(f64, -0x8p-152)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, -0x8p-152)).cosx);
    try std.testing.expectEqual(-0x4p-1076, sincos(@as(f64, -0x4p-1076)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f64, -0x4p-1076)).cosx);
    try std.testing.expectEqual(0xf.fa2add3e58948p-4, sincos(@as(f64, 0x1.8475e6p+0)).sinx);
    try std.testing.expectEqual(0xd.a8263394be6dp-8, sincos(@as(f64, 0x1.8475e6p+0)).cosx);
    try std.testing.expectEqual(0xf.fa2adb8953aep-4, sincos(@as(f64, 0x1.8475e4p+0)).sinx);
    try std.testing.expectEqual(0xd.a82832da19f98p-8, sincos(@as(f64, 0x1.8475e4p+0)).cosx);
    try std.testing.expectEqual(0xf.fa2adcf9ea84p-4, sincos(@as(f64, 0x1.8475e5afd4481p+0)).sinx);
    try std.testing.expectEqual(0xd.a82683a33cbe8p-8, sincos(@as(f64, 0x1.8475e5afd4481p+0)).cosx);

    try std.testing.expectEqual(0x0p+0, sincos(@as(f80, 0x0p+0)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x0p+0)).cosx);
    try std.testing.expectEqual(-0x0p+0, sincos(@as(f80, -0x0p+0)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x0p+0)).cosx);
    try std.testing.expectEqual(0xf.fffffffffffbb29p-4, sincos(@as(f80, 0x1.921fb6p+0)).sinx);
    try std.testing.expectEqual(-0xb.bbd2e7b96766267p-28, sincos(@as(f80, 0x1.921fb6p+0)).cosx);
    try std.testing.expectEqual(0xf.fffffffffff32a3p-4, sincos(@as(f80, 0x1.921fb4p+0)).sinx);
    try std.testing.expectEqual(0x1.4442d1846989361p-24, sincos(@as(f80, 0x1.921fb4p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x1.921fb54442d19p+0)).sinx);
    try std.testing.expectEqual(-0xb.9676733ae8fe47cp-56, sincos(@as(f80, 0x1.921fb54442d19p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x1.921fb54442d18p+0)).sinx);
    try std.testing.expectEqual(0x4.69898cc51701b838p-56, sincos(@as(f80, 0x1.921fb54442d18p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x1.921fb54442d1846ap+0)).sinx);
    try std.testing.expectEqual(-0x7.6733ae8fe47c65d8p-68, sincos(@as(f80, 0x1.921fb54442d1846ap+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x1.921fb54442d18468p+0)).sinx);
    try std.testing.expectEqual(0x1.898cc51701b839a2p-64, sincos(@as(f80, 0x1.921fb54442d18468p+0)).cosx);
    try std.testing.expectEqual(0x8.000003632116885p-4, sincos(@as(f80, 0x8.60a92p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d7237832e30fp-4, sincos(@as(f80, 0x8.60a92p-4)).cosx);
    try std.testing.expectEqual(0x7.fffff587e3a050dp-4, sincos(@as(f80, 0x8.60a91p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d7a37832aa68p-4, sincos(@as(f80, 0x8.60a91p-4)).cosx);
    try std.testing.expectEqual(0x8.000000000000358p-4, sincos(@as(f80, 0x8.60a91c16b9b3p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c26551afp-4, sincos(@as(f80, 0x8.60a91c16b9b3p-4)).cosx);
    try std.testing.expectEqual(0x7.ffffffffffffc6a8p-4, sincos(@as(f80, 0x8.60a91c16b9b28p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c26555afp-4, sincos(@as(f80, 0x8.60a91c16b9b28p-4)).cosx);
    try std.testing.expectEqual(0x8.000000000000001p-4, sincos(@as(f80, 0x8.60a91c16b9b2c24p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c265539dp-4, sincos(@as(f80, 0x8.60a91c16b9b2c24p-4)).cosx);
    try std.testing.expectEqual(0x8p-4, sincos(@as(f80, 0x8.60a91c16b9b2c23p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c265539ep-4, sincos(@as(f80, 0x8.60a91c16b9b2c23p-4)).cosx);
    try std.testing.expectEqual(0xd.db3d78156ca0cfbp-4, sincos(@as(f80, 0x1.0c1524p+0)).sinx);
    try std.testing.expectEqual(0x7.fffff939bdd18038p-4, sincos(@as(f80, 0x1.0c1524p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d68156c92a5cp-4, sincos(@as(f80, 0x1.0c1522p+0)).sinx);
    try std.testing.expectEqual(0x8.000014f038b1ab1p-4, sincos(@as(f80, 0x1.0c1522p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265577ap-4, sincos(@as(f80, 0x1.0c152382d7366p+0)).sinx);
    try std.testing.expectEqual(0x7.ffffffffffff94f8p-4, sincos(@as(f80, 0x1.0c152382d7366p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c2654f7ap-4, sincos(@as(f80, 0x1.0c152382d7365p+0)).sinx);
    try std.testing.expectEqual(0x8.00000000000072bp-4, sincos(@as(f80, 0x1.0c152382d7365p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265539ep-4, sincos(@as(f80, 0x1.0c152382d7365848p+0)).sinx);
    try std.testing.expectEqual(0x7.ffffffffffffffe8p-4, sincos(@as(f80, 0x1.0c152382d7365848p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265539dp-4, sincos(@as(f80, 0x1.0c152382d7365846p+0)).sinx);
    try std.testing.expectEqual(0x8p-4, sincos(@as(f80, 0x1.0c152382d7365846p+0)).cosx);
    try std.testing.expectEqual(-0x1.777a5cf72cec5fd6p-24, sincos(@as(f80, 0x3.243f6cp+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffeeca4p-4, sincos(@as(f80, 0x3.243f6cp+0)).cosx);
    try std.testing.expectEqual(0x2.8885a308d31063e4p-24, sincos(@as(f80, 0x3.243f68p+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffcca8ep-4, sincos(@as(f80, 0x3.243f68p+0)).cosx);
    try std.testing.expectEqual(-0x1.72cece675d1fc8f8p-52, sincos(@as(f80, 0x3.243f6a8885a32p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f80, 0x3.243f6a8885a32p+0)).cosx);
    try std.testing.expectEqual(0x8.d313198a2e03707p-56, sincos(@as(f80, 0x3.243f6a8885a3p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f80, 0x3.243f6a8885a3p+0)).cosx);
    try std.testing.expectEqual(-0xe.ce675d1fc8f8cbbp-68, sincos(@as(f80, 0x3.243f6a8885a308d4p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f80, 0x3.243f6a8885a308d4p+0)).cosx);
    try std.testing.expectEqual(0x3.13198a2e03707344p-64, sincos(@as(f80, 0x3.243f6a8885a308dp+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f80, 0x3.243f6a8885a308dp+0)).cosx);
    try std.testing.expectEqual(0x1.777a5cf72cec5fd6p-24, sincos(@as(f80, -0x3.243f6cp+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffeeca4p-4, sincos(@as(f80, -0x3.243f6cp+0)).cosx);
    try std.testing.expectEqual(-0x2.8885a308d31063e4p-24, sincos(@as(f80, -0x3.243f68p+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffcca8ep-4, sincos(@as(f80, -0x3.243f68p+0)).cosx);
    try std.testing.expectEqual(0x1.72cece675d1fc8f8p-52, sincos(@as(f80, -0x3.243f6a8885a32p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f80, -0x3.243f6a8885a32p+0)).cosx);
    try std.testing.expectEqual(-0x8.d313198a2e03707p-56, sincos(@as(f80, -0x3.243f6a8885a3p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f80, -0x3.243f6a8885a3p+0)).cosx);
    try std.testing.expectEqual(0xe.ce675d1fc8f8cbbp-68, sincos(@as(f80, -0x3.243f6a8885a308d4p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f80, -0x3.243f6a8885a308d4p+0)).cosx);
    try std.testing.expectEqual(-0x3.13198a2e03707344p-64, sincos(@as(f80, -0x3.243f6a8885a308dp+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f80, -0x3.243f6a8885a308dp+0)).cosx);
    try std.testing.expectEqual(0xa.e7fe0b5fc786b2ep-4, sincos(@as(f80, 0xcp-4)).sinx);
    try std.testing.expectEqual(0xb.b4ff632a908f73fp-4, sincos(@as(f80, 0xcp-4)).cosx);
    try std.testing.expectEqual(-0xc.143e153b0701e8p-8, sincos(@as(f80, 0x2p+64)).sinx);
    try std.testing.expectEqual(0xf.fb701e22987fbe7p-4, sincos(@as(f80, 0x2p+64)).cosx);
    try std.testing.expectEqual(0xc.143e153b0701e8p-8, sincos(@as(f80, -0x2p+64)).sinx);
    try std.testing.expectEqual(0xf.fb701e22987fbe7p-4, sincos(@as(f80, -0x2p+64)).cosx);
    try std.testing.expectEqual(0xb.7fb600275877a6p-4, sincos(@as(f80, 0xc.d4967p-4)).sinx);
    try std.testing.expectEqual(0xb.201e77869a46ae2p-4, sincos(@as(f80, 0xc.d4967p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5f50739fa5f9p-4, sincos(@as(f80, 0xc.d4966p-4)).sinx);
    try std.testing.expectEqual(0xb.201e83065041457p-4, sincos(@as(f80, 0xc.d4966p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe7769793e6p-4, sincos(@as(f80, 0xc.d4966d92d171p-4)).sinx);
    try std.testing.expectEqual(0xb.201e7945088464p-4, sincos(@as(f80, 0xc.d4966d92d171p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e56p-4, sincos(@as(f80, 0xc.d4966d92d1708p-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884cp-4, sincos(@as(f80, 0xc.d4966d92d1708p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e74p-4, sincos(@as(f80, 0xc.d4966d92d17082ap-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884be2p-4, sincos(@as(f80, 0xc.d4966d92d17082ap-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e73p-4, sincos(@as(f80, 0xc.d4966d92d170829p-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884be3p-4, sincos(@as(f80, 0xc.d4966d92d170829p-4)).cosx);
    try std.testing.expectEqual(-0x4.cd7e86c4077bf0ep-4, sincos(@as(f80, 0x2.1e19e4p+72)).sinx);
    try std.testing.expectEqual(0xf.431dd7a36cf37dep-4, sincos(@as(f80, 0x2.1e19e4p+72)).cosx);
    try std.testing.expectEqual(-0xb.becc47ab1b8c708p-4, sincos(@as(f80, 0x2.1e19ep+72)).sinx);
    try std.testing.expectEqual(0xa.dd6f6bacd20654cp-4, sincos(@as(f80, 0x2.1e19ep+72)).cosx);
    try std.testing.expectEqual(-0xd.a29d5bb5f9cb87dp-4, sincos(@as(f80, 0x2.1e19e0c9bab24p+72)).sinx);
    try std.testing.expectEqual(0x8.5f167780e479c9ap-4, sincos(@as(f80, 0x2.1e19e0c9bab24p+72)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba907p-4, sincos(@as(f80, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, sincos(@as(f80, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x9.0292465edbaff2dp-4, sincos(@as(f80, 0x8p+1020)).sinx);
    try std.testing.expectEqual(-0xd.38cf9361195f50bp-4, sincos(@as(f80, 0x8p+1020)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba907p-4, sincos(@as(f80, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, sincos(@as(f80, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x1.452fc98b34e96b62p-8, sincos(@as(f80, 0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba9ep-4, sincos(@as(f80, 0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(0x6.3ad4b2136cc6882p-4, sincos(@as(f80, 0x8p+16380)).sinx);
    try std.testing.expectEqual(0xe.bcc2fc82ae39ecp-4, sincos(@as(f80, 0x8p+16380)).cosx);
    try std.testing.expectEqual(0x6.0b8d19579bf2db6p-4, sincos(@as(f80, 0x1p+120)).sinx);
    try std.testing.expectEqual(-0xe.d06685b36c66c4dp-4, sincos(@as(f80, 0x1p+120)).cosx);
    try std.testing.expectEqual(0x9.f963166f215eb89p-4, sincos(@as(f80, 0x8p+124)).sinx);
    try std.testing.expectEqual(0xc.82b8ec98b5e62fdp-4, sincos(@as(f80, 0x8p+124)).cosx);
    try std.testing.expectEqual(0xc.6fa5c5665984d89p-8, sincos(@as(f80, 0xf.ffffcp+124)).sinx);
    try std.testing.expectEqual(0xf.fb2a030c5ae20bep-4, sincos(@as(f80, 0xf.ffffcp+124)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba907p-4, sincos(@as(f80, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, sincos(@as(f80, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x7.f13d78eabb76b8a8p-4, sincos(@as(f80, 0x4p+48)).sinx);
    try std.testing.expectEqual(0xd.e3b88804f00552dp-4, sincos(@as(f80, 0x4p+48)).cosx);
    try std.testing.expectEqual(-0xf.c777c6b36a750a6p-4, sincos(@as(f80, 0x1p+28)).sinx);
    try std.testing.expectEqual(-0x2.a62ba8824e5bcb08p-4, sincos(@as(f80, 0x1p+28)).cosx);
    try std.testing.expectEqual(0x8.599b32844aba907p-4, sincos(@as(f80, -0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, sincos(@as(f80, -0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(-0xf.e00885042dd770dp-4, sincos(@as(f80, -0x3.3de320f6be87ep+1020)).sinx);
    try std.testing.expectEqual(-0x1.febbf9949ecc1336p-4, sincos(@as(f80, -0x3.3de320f6be87ep+1020)).cosx);
    try std.testing.expectEqual(0xc.773a2eac3000ddfp-4, sincos(@as(f80, 0xe.9f1e6p+112)).sinx);
    try std.testing.expectEqual(-0xa.07bd3ab53ab9711p-4, sincos(@as(f80, 0xe.9f1e6p+112)).cosx);
    try std.testing.expectEqual(0x7.76d600e031521b8p-4, sincos(@as(f80, 0xe.9f1e5p+112)).sinx);
    try std.testing.expectEqual(0xe.26f8af8333f9271p-4, sincos(@as(f80, 0xe.9f1e5p+112)).cosx);
    try std.testing.expectEqual(0xf.dfffd7bde0fb4ecp-4, sincos(@as(f80, 0xe.9f1e5bc3bb88p+112)).sinx);
    try std.testing.expectEqual(0x1.ff01000c9ae7363p-4, sincos(@as(f80, 0xe.9f1e5bc3bb88p+112)).cosx);
    try std.testing.expectEqual(-0x1.ffb679ba994b7618p-4, sincos(@as(f80, 0x4.7857dp+68)).sinx);
    try std.testing.expectEqual(-0xf.dfe902135fc1c18p-4, sincos(@as(f80, 0x4.7857dp+68)).cosx);
    try std.testing.expectEqual(-0x1.fecaff6878a10ce6p-4, sincos(@as(f80, 0x6.287cdp+0)).sinx);
    try std.testing.expectEqual(0xf.e006a1ad17db69bp-4, sincos(@as(f80, 0x6.287cdp+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb7e68ad6e9c4p-4, sincos(@as(f80, 0x6.287cc8p+0)).sinx);
    try std.testing.expectEqual(0xf.e00691b6bde4252p-4, sincos(@as(f80, 0x6.287cc8p+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb772e1b8300e6p-4, sincos(@as(f80, 0x6.287cc8749213p+0)).sinx);
    try std.testing.expectEqual(0xf.e006929f558dbe6p-4, sincos(@as(f80, 0x6.287cc8749213p+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb772e1b86f8e8p-4, sincos(@as(f80, 0x6.287cc8749212cp+0)).sinx);
    try std.testing.expectEqual(0xf.e006929f558d3ebp-4, sincos(@as(f80, 0x6.287cc8749212cp+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb772e1b848bcap-4, sincos(@as(f80, 0x6.287cc8749212e72p+0)).sinx);
    try std.testing.expectEqual(0xf.e006929f558d8ccp-4, sincos(@as(f80, 0x6.287cc8749212e72p+0)).cosx);
    try std.testing.expectEqual(-0xd.8f691a7a95426p-4, sincos(@as(f80, -0x1.02e34cp+0)).sinx);
    try std.testing.expectEqual(0x8.7e0ea4db2f48867p-4, sincos(@as(f80, -0x1.02e34cp+0)).cosx);
    try std.testing.expectEqual(-0x8.3bee07bc9076425p-4, sincos(@as(f80, 0xf.f0274p+4)).sinx);
    try std.testing.expectEqual(-0xd.b7f5359babdb66cp-4, sincos(@as(f80, 0xf.f0274p+4)).cosx);
    try std.testing.expectEqual(0x1.ffc6da9f1ffed896p-4, sincos(@as(f80, 0x3.042d88p+0)).sinx);
    try std.testing.expectEqual(-0xf.dfe6f2169e24f27p-4, sincos(@as(f80, 0x3.042d88p+0)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba907p-4, sincos(@as(f80, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, sincos(@as(f80, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x1.452fc98b34e96b62p-8, sincos(@as(f80, 0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba9ep-4, sincos(@as(f80, 0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(0xf.dfd9d4b6d0e5f7cp-4, sincos(@as(f80, 0xf.fffffffffffffffp+16380)).sinx);
    try std.testing.expectEqual(-0x2.002ef4018753d50cp-4, sincos(@as(f80, 0xf.fffffffffffffffp+16380)).cosx);
    try std.testing.expectEqual(0x8.599b32844aba907p-4, sincos(@as(f80, -0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d53p-4, sincos(@as(f80, -0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(-0x1.452fc98b34e96b62p-8, sincos(@as(f80, -0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba9ep-4, sincos(@as(f80, -0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(-0xf.dfd9d4b6d0e5f7cp-4, sincos(@as(f80, -0xf.fffffffffffffffp+16380)).sinx);
    try std.testing.expectEqual(-0x2.002ef4018753d50cp-4, sincos(@as(f80, -0xf.fffffffffffffffp+16380)).cosx);
    try std.testing.expectEqual(0x4p-128, sincos(@as(f80, 0x4p-128)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x4p-128)).cosx);
    try std.testing.expectEqual(0x4p-1024, sincos(@as(f80, 0x4p-1024)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x4p-1024)).cosx);
    try std.testing.expectEqual(0x4p-16384, sincos(@as(f80, 0x4p-16384)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x4p-16384)).cosx);
    try std.testing.expectEqual(0x2p-16384, sincos(@as(f80, 0x2p-16384)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x2p-16384)).cosx);
    try std.testing.expectEqual(0x8p-972, sincos(@as(f80, 0x8p-972)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x8p-972)).cosx);
    try std.testing.expectEqual(-0x4p-128, sincos(@as(f80, -0x4p-128)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x4p-128)).cosx);
    try std.testing.expectEqual(-0x4p-1024, sincos(@as(f80, -0x4p-1024)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x4p-1024)).cosx);
    try std.testing.expectEqual(-0x4p-16384, sincos(@as(f80, -0x4p-16384)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x4p-16384)).cosx);
    try std.testing.expectEqual(-0x2p-16384, sincos(@as(f80, -0x2p-16384)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x2p-16384)).cosx);
    try std.testing.expectEqual(-0x8p-972, sincos(@as(f80, -0x8p-972)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x8p-972)).cosx);
    try std.testing.expectEqual(0x8p-152, sincos(@as(f80, 0x8p-152)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x8p-152)).cosx);
    try std.testing.expectEqual(0x4p-1076, sincos(@as(f80, 0x4p-1076)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x4p-1076)).cosx);
    try std.testing.expectEqual(0x8p-16448, sincos(@as(f80, 0x8p-16448)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, 0x8p-16448)).cosx);
    try std.testing.expectEqual(-0x8p-152, sincos(@as(f80, -0x8p-152)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x8p-152)).cosx);
    try std.testing.expectEqual(-0x4p-1076, sincos(@as(f80, -0x4p-1076)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x4p-1076)).cosx);
    try std.testing.expectEqual(-0x8p-16448, sincos(@as(f80, -0x8p-16448)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f80, -0x8p-16448)).cosx);
    try std.testing.expectEqual(0xf.fa2add3e58948d1p-4, sincos(@as(f80, 0x1.8475e6p+0)).sinx);
    try std.testing.expectEqual(0xd.a8263394be6d0e6p-8, sincos(@as(f80, 0x1.8475e6p+0)).cosx);
    try std.testing.expectEqual(0xf.fa2adb8953ae262p-4, sincos(@as(f80, 0x1.8475e4p+0)).sinx);
    try std.testing.expectEqual(0xd.a82832da19f9892p-8, sincos(@as(f80, 0x1.8475e4p+0)).cosx);
    try std.testing.expectEqual(0xf.fa2adcf9ea83dbep-4, sincos(@as(f80, 0x1.8475e5afd4481p+0)).sinx);
    try std.testing.expectEqual(0xd.a82683a33cbecp-8, sincos(@as(f80, 0x1.8475e5afd4481p+0)).cosx);

    try std.testing.expectEqual(0x0p+0, sincos(@as(f128, 0x0p+0)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x0p+0)).cosx);
    try std.testing.expectEqual(-0x0p+0, sincos(@as(f128, -0x0p+0)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x0p+0)).cosx);
    try std.testing.expectEqual(0xf.fffffffffffbb290924e3a114988p-4, sincos(@as(f128, 0x1.921fb6p+0)).sinx);
    try std.testing.expectEqual(-0xb.bbd2e7b96766266f1d18f3ead01p-28, sincos(@as(f128, 0x1.921fb6p+0)).cosx);
    try std.testing.expectEqual(0xf.fffffffffff32a3661c108e136d8p-4, sincos(@as(f128, 0x1.921fb4p+0)).sinx);
    try std.testing.expectEqual(0x1.4442d18469893610281a0f9b0e8dp-24, sincos(@as(f128, 0x1.921fb4p+0)).cosx);
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffbdp-4, sincos(@as(f128, 0x1.921fb54442d19p+0)).sinx);
    try std.testing.expectEqual(-0xb.9676733ae8fe47c65dadfb63ede8p-56, sincos(@as(f128, 0x1.921fb54442d19p+0)).cosx);
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffff68p-4, sincos(@as(f128, 0x1.921fb54442d18p+0)).sinx);
    // try std.testing.expectEqual(0x4.69898cc51701b839a252049c1108p-56, sincos(@as(f128, 0x1.921fb54442d18p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x1.921fb54442d1846ap+0)).sinx);
    try std.testing.expectEqual(-0x7.6733ae8fe47c65dadfb63eeeb308p-68, sincos(@as(f128, 0x1.921fb54442d1846ap+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x1.921fb54442d18468p+0)).sinx);
    try std.testing.expectEqual(0x1.898cc51701b839a252049c1114dp-64, sincos(@as(f128, 0x1.921fb54442d18468p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x1.921fb54442d18469898cc51701b9p+0)).sinx);
    try std.testing.expectEqual(-0xc.65dadfb63eeeb306717fbe882b38p-116, sincos(@as(f128, 0x1.921fb54442d18469898cc51701b9p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x1.921fb54442d18469898cc51701b8p+0)).sinx);
    try std.testing.expectEqual(0x3.9a252049c1114cf98e804177d4c8p-116, sincos(@as(f128, 0x1.921fb54442d18469898cc51701b8p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x1.921fb54442d18469898cc51702p+0)).sinx);
    try std.testing.expectEqual(-0x4.7c65dadfb63eeeb306717fbe882cp-108, sincos(@as(f128, 0x1.921fb54442d18469898cc51702p+0)).cosx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x1.921fb54442d18469898cc517018p+0)).sinx);
    try std.testing.expectEqual(0x3.839a252049c1114cf98e804177d4p-108, sincos(@as(f128, 0x1.921fb54442d18469898cc517018p+0)).cosx);
    try std.testing.expectEqual(0x8.0000036321168852c4130c64b4cp-4, sincos(@as(f128, 0x8.60a92p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d7237832e30f6465d599654a8p-4, sincos(@as(f128, 0x8.60a92p-4)).cosx);
    try std.testing.expectEqual(0x7.fffff587e3a050cf967fba7bc728p-4, sincos(@as(f128, 0x8.60a91p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d7a37832aa678a274956dfd3p-4, sincos(@as(f128, 0x8.60a91p-4)).cosx);
    try std.testing.expectEqual(0x8.00000000000035858118fd5157ep-4, sincos(@as(f128, 0x8.60a91c16b9b3p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c26551af29a6e23c11f48p-4, sincos(@as(f128, 0x8.60a91c16b9b3p-4)).cosx);
    try std.testing.expectEqual(0x7.ffffffffffffc6ab95779c1eae0cp-4, sincos(@as(f128, 0x8.60a91c16b9b28p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c26555af29a6e23c11f38p-4, sincos(@as(f128, 0x8.60a91c16b9b28p-4)).cosx);
    try std.testing.expectEqual(0x8.000000000000000b5feca2ed673p-4, sincos(@as(f128, 0x8.60a91c16b9b2c24p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c265539d29a6e23c11fbp-4, sincos(@as(f128, 0x8.60a91c16b9b2c24p-4)).cosx);
    try std.testing.expectEqual(0x7.fffffffffffffffd84af2ec140dcp-4, sincos(@as(f128, 0x8.60a91c16b9b2c23p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c265539da9a6e23c11fbp-4, sincos(@as(f128, 0x8.60a91c16b9b2c23p-4)).cosx);
    try std.testing.expectEqual(0x8p-4, sincos(@as(f128, 0x8.60a91c16b9b2c232dd99707ab3d8p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c5cp-4, sincos(@as(f128, 0x8.60a91c16b9b2c232dd99707ab3d8p-4)).cosx);
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffffffcp-4, sincos(@as(f128, 0x8.60a91c16b9b2c232dd99707ab3dp-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c5c8p-4, sincos(@as(f128, 0x8.60a91c16b9b2c232dd99707ab3dp-4)).cosx);
    try std.testing.expectEqual(0x8.000000000000000000000000002p-4, sincos(@as(f128, 0x8.60a91c16b9b2c232dd99707ab4p-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c5bp-4, sincos(@as(f128, 0x8.60a91c16b9b2c232dd99707ab4p-4)).cosx);
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffffcacp-4, sincos(@as(f128, 0x8.60a91c16b9b2c232dd99707abp-4)).sinx);
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c7bp-4, sincos(@as(f128, 0x8.60a91c16b9b2c232dd99707abp-4)).cosx);
    try std.testing.expectEqual(0xd.db3d78156ca0cfb4fd88fd27f7ep-4, sincos(@as(f128, 0x1.0c1524p+0)).sinx);
    try std.testing.expectEqual(0x7.fffff939bdd18035537d20fef1b4p-4, sincos(@as(f128, 0x1.0c1524p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d68156c92a5be750863ea3d58p-4, sincos(@as(f128, 0x1.0c1522p+0)).sinx);
    try std.testing.expectEqual(0x8.000014f038b1ab0e902f6811916p-4, sincos(@as(f128, 0x1.0c1522p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265577a64e07fb09105p-4, sincos(@as(f128, 0x1.0c152382d7366p+0)).sinx);
    try std.testing.expectEqual(0x7.ffffffffffff94f4fdce055d4ed4p-4, sincos(@as(f128, 0x1.0c152382d7366p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c2654f7a64e07fb09101p-4, sincos(@as(f128, 0x1.0c152382d7365p+0)).sinx);
    try std.testing.expectEqual(0x8.00000000000072a8d510c7c2a25p-4, sincos(@as(f128, 0x1.0c152382d7365p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265539e64e07fb0911e8p-4, sincos(@as(f128, 0x1.0c152382d7365848p+0)).sinx);
    try std.testing.expectEqual(0x7.ffffffffffffffe94026ba25319cp-4, sincos(@as(f128, 0x1.0c152382d7365848p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265539d64e07fb0911e8p-4, sincos(@as(f128, 0x1.0c152382d7365846p+0)).sinx);
    try std.testing.expectEqual(0x8.0000000000000004f6a1a27d7e48p-4, sincos(@as(f128, 0x1.0c152382d7365846p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c5cp-4, sincos(@as(f128, 0x1.0c152382d73658465bb32e0f567bp+0)).sinx);
    try std.testing.expectEqual(0x7.fffffffffffffffffffffffffffcp-4, sincos(@as(f128, 0x1.0c152382d73658465bb32e0f567bp+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c5b8p-4, sincos(@as(f128, 0x1.0c152382d73658465bb32e0f567ap+0)).sinx);
    try std.testing.expectEqual(0x8.0000000000000000000000000008p-4, sincos(@as(f128, 0x1.0c152382d73658465bb32e0f567ap+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c5e8p-4, sincos(@as(f128, 0x1.0c152382d73658465bb32e0f568p+0)).sinx);
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffffffb8p-4, sincos(@as(f128, 0x1.0c152382d73658465bb32e0f568p+0)).cosx);
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c1e8p-4, sincos(@as(f128, 0x1.0c152382d73658465bb32e0f56p+0)).sinx);
    try std.testing.expectEqual(0x8.00000000000000000000000006a8p-4, sincos(@as(f128, 0x1.0c152382d73658465bb32e0f56p+0)).cosx);
    try std.testing.expectEqual(-0x1.777a5cf72cec5fd61896cb4f40d2p-24, sincos(@as(f128, 0x3.243f6cp+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffeeca424938e8477678p-4, sincos(@as(f128, 0x3.243f6cp+0)).cosx);
    try std.testing.expectEqual(0x2.8885a308d31063e2b6c62b7f4d6cp-24, sincos(@as(f128, 0x3.243f68p+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffcca8d9870423997308p-4, sincos(@as(f128, 0x3.243f68p+0)).cosx);
    // try std.testing.expectEqual(-0x1.72cece675d1fc8f8cbb5bf6c7d5cp-52, sincos(@as(f128, 0x3.243f6a8885a32p+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffffffffffffffffef38p-4, sincos(@as(f128, 0x3.243f6a8885a32p+0)).cosx);
    try std.testing.expectEqual(0x8.d313198a2e03707344a4093821b8p-56, sincos(@as(f128, 0x3.243f6a8885a3p+0)).sinx);
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffffd9p-4, sincos(@as(f128, 0x3.243f6a8885a3p+0)).cosx);
    try std.testing.expectEqual(-0xe.ce675d1fc8f8cbb5bf6c7ddd661p-68, sincos(@as(f128, 0x3.243f6a8885a308d4p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, 0x3.243f6a8885a308d4p+0)).cosx);
    try std.testing.expectEqual(0x3.13198a2e03707344a409382229ap-64, sincos(@as(f128, 0x3.243f6a8885a308dp+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, 0x3.243f6a8885a308dp+0)).cosx);
    try std.testing.expectEqual(-0x1.8cbb5bf6c7ddd660ce2ff7d10567p-112, sincos(@as(f128, 0x3.243f6a8885a308d313198a2e0372p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, 0x3.243f6a8885a308d313198a2e0372p+0)).cosx);
    try std.testing.expectEqual(0x7.344a4093822299f31d0082efa99p-116, sincos(@as(f128, 0x3.243f6a8885a308d313198a2e037p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, 0x3.243f6a8885a308d313198a2e037p+0)).cosx);
    try std.testing.expectEqual(-0x8.f8cbb5bf6c7ddd660ce2ff7d1058p-108, sincos(@as(f128, 0x3.243f6a8885a308d313198a2e04p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, 0x3.243f6a8885a308d313198a2e04p+0)).cosx);
    try std.testing.expectEqual(0x7.07344a4093822299f31d0082efa8p-108, sincos(@as(f128, 0x3.243f6a8885a308d313198a2e03p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, 0x3.243f6a8885a308d313198a2e03p+0)).cosx);
    try std.testing.expectEqual(0x1.777a5cf72cec5fd61896cb4f40d2p-24, sincos(@as(f128, -0x3.243f6cp+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffeeca424938e8477678p-4, sincos(@as(f128, -0x3.243f6cp+0)).cosx);
    try std.testing.expectEqual(-0x2.8885a308d31063e2b6c62b7f4d6cp-24, sincos(@as(f128, -0x3.243f68p+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffcca8d9870423997308p-4, sincos(@as(f128, -0x3.243f68p+0)).cosx);
    // try std.testing.expectEqual(0x1.72cece675d1fc8f8cbb5bf6c7d5cp-52, sincos(@as(f128, -0x3.243f6a8885a32p+0)).sinx);
    try std.testing.expectEqual(-0xf.ffffffffffffffffffffffffef38p-4, sincos(@as(f128, -0x3.243f6a8885a32p+0)).cosx);
    try std.testing.expectEqual(-0x8.d313198a2e03707344a4093821b8p-56, sincos(@as(f128, -0x3.243f6a8885a3p+0)).sinx);
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffffd9p-4, sincos(@as(f128, -0x3.243f6a8885a3p+0)).cosx);
    try std.testing.expectEqual(0xe.ce675d1fc8f8cbb5bf6c7ddd661p-68, sincos(@as(f128, -0x3.243f6a8885a308d4p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, -0x3.243f6a8885a308d4p+0)).cosx);
    try std.testing.expectEqual(-0x3.13198a2e03707344a409382229ap-64, sincos(@as(f128, -0x3.243f6a8885a308dp+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, -0x3.243f6a8885a308dp+0)).cosx);
    try std.testing.expectEqual(0x1.8cbb5bf6c7ddd660ce2ff7d10567p-112, sincos(@as(f128, -0x3.243f6a8885a308d313198a2e0372p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, -0x3.243f6a8885a308d313198a2e0372p+0)).cosx);
    try std.testing.expectEqual(-0x7.344a4093822299f31d0082efa99p-116, sincos(@as(f128, -0x3.243f6a8885a308d313198a2e037p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, -0x3.243f6a8885a308d313198a2e037p+0)).cosx);
    try std.testing.expectEqual(0x8.f8cbb5bf6c7ddd660ce2ff7d1058p-108, sincos(@as(f128, -0x3.243f6a8885a308d313198a2e04p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, -0x3.243f6a8885a308d313198a2e04p+0)).cosx);
    try std.testing.expectEqual(-0x7.07344a4093822299f31d0082efa8p-108, sincos(@as(f128, -0x3.243f6a8885a308d313198a2e03p+0)).sinx);
    try std.testing.expectEqual(-0x1p+0, sincos(@as(f128, -0x3.243f6a8885a308d313198a2e03p+0)).cosx);
    try std.testing.expectEqual(0xa.e7fe0b5fc786b2d966e1d6af1408p-4, sincos(@as(f128, 0xcp-4)).sinx);
    try std.testing.expectEqual(0xb.b4ff632a908f73ec151839cb9d98p-4, sincos(@as(f128, 0xcp-4)).cosx);
    try std.testing.expectEqual(-0xc.143e153b0701e800f9b8a47b75b8p-8, sincos(@as(f128, 0x2p+64)).sinx);
    try std.testing.expectEqual(0xf.fb701e22987fbe68852ee2bc897p-4, sincos(@as(f128, 0x2p+64)).cosx);
    try std.testing.expectEqual(0xc.143e153b0701e800f9b8a47b75b8p-8, sincos(@as(f128, -0x2p+64)).sinx);
    try std.testing.expectEqual(0xf.fb701e22987fbe68852ee2bc897p-4, sincos(@as(f128, -0x2p+64)).cosx);
    try std.testing.expectEqual(0xb.7fb600275877a60011766c8a3178p-4, sincos(@as(f128, 0xc.d4967p-4)).sinx);
    try std.testing.expectEqual(0xb.201e77869a46ae20ce545c5c67p-4, sincos(@as(f128, 0xc.d4967p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5f50739fa5f8acc8f4f3f1b3p-4, sincos(@as(f128, 0xc.d4966p-4)).sinx);
    try std.testing.expectEqual(0xb.201e83065041456a084c70f5a128p-4, sincos(@as(f128, 0xc.d4966p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe7769793e65c978bd3cef98p-4, sincos(@as(f128, 0xc.d4966d92d171p-4)).sinx);
    try std.testing.expectEqual(0xb.201e794508846402500c44b4f8e8p-4, sincos(@as(f128, 0xc.d4966d92d171p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e564d5ae94f8cb08p-4, sincos(@as(f128, 0xc.d4966d92d1708p-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884c00000000000c178p-4, sincos(@as(f128, 0xc.d4966d92d1708p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e7381aae7a4c30dp-4, sincos(@as(f128, 0xc.d4966d92d17082ap-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884be1d0c24406973ap-4, sincos(@as(f128, 0xc.d4966d92d17082ap-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e72cfa9001072848p-4, sincos(@as(f128, 0xc.d4966d92d170829p-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884be288bda3ee0dd18p-4, sincos(@as(f128, 0xc.d4966d92d170829p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e732912810356318p-4, sincos(@as(f128, 0xc.d4966d92d17082980965c1a663c8p-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884be22c53e47ed16d8p-4, sincos(@as(f128, 0xc.d4966d92d17082980965c1a663c8p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e732912810356318p-4, sincos(@as(f128, 0xc.d4966d92d17082980965c1a663cp-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884be22c53e47ed16ep-4, sincos(@as(f128, 0xc.d4966d92d17082980965c1a663cp-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e73291281035634p-4, sincos(@as(f128, 0xc.d4966d92d17082980965c1a664p-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884be22c53e47ed16bp-4, sincos(@as(f128, 0xc.d4966d92d17082980965c1a664p-4)).cosx);
    try std.testing.expectEqual(0xb.7fb5fe776978e732912810356078p-4, sincos(@as(f128, 0xc.d4966d92d17082980965c1a66p-4)).sinx);
    try std.testing.expectEqual(0xb.201e79450884be22c53e47ed199p-4, sincos(@as(f128, 0xc.d4966d92d17082980965c1a66p-4)).cosx);
    try std.testing.expectEqual(-0x4.cd7e86c4077bf0debc87d70d196p-4, sincos(@as(f128, 0x2.1e19e4p+72)).sinx);
    try std.testing.expectEqual(0xf.431dd7a36cf37de5c74544f6b438p-4, sincos(@as(f128, 0x2.1e19e4p+72)).cosx);
    try std.testing.expectEqual(-0xb.becc47ab1b8c70793712c4ff2bcp-4, sincos(@as(f128, 0x2.1e19ep+72)).sinx);
    try std.testing.expectEqual(0xa.dd6f6bacd20654c1404f52cde16p-4, sincos(@as(f128, 0x2.1e19ep+72)).cosx);
    try std.testing.expectEqual(-0xd.a29d5bb5f9cb87d14de41dc991fp-4, sincos(@as(f128, 0x2.1e19e0c9bab24p+72)).sinx);
    try std.testing.expectEqual(0x8.5f167780e479c9a5c86ffce7615p-4, sincos(@as(f128, 0x2.1e19e0c9bab24p+72)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba906cee446be04998p-4, sincos(@as(f128, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, sincos(@as(f128, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x9.0292465edbaff2d2e64a2845e558p-4, sincos(@as(f128, 0x8p+1020)).sinx);
    try std.testing.expectEqual(-0xd.38cf9361195f50b10fac29dd9038p-4, sincos(@as(f128, 0x8p+1020)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba906cee446be04998p-4, sincos(@as(f128, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, sincos(@as(f128, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x1.452fc98b34e96b61139b09a7c84ap-8, sincos(@as(f128, 0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba9e038d934070f138p-4, sincos(@as(f128, 0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(0x6.3ad4b2136cc6881f0ca607c7946p-4, sincos(@as(f128, 0x8p+16380)).sinx);
    try std.testing.expectEqual(0xe.bcc2fc82ae39ebf8da5d687bf36p-4, sincos(@as(f128, 0x8p+16380)).cosx);
    try std.testing.expectEqual(-0xe.f1a3e1dc468a921dddb4e37fbe6p-4, sincos(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)).sinx);
    try std.testing.expectEqual(-0x5.b773d971a848e75c230605526978p-4, sincos(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)).cosx);
    try std.testing.expectEqual(0x6.0b8d19579bf2db5e5f1aa933f37cp-4, sincos(@as(f128, 0x1p+120)).sinx);
    try std.testing.expectEqual(-0xe.d06685b36c66c4cf35c11f6519p-4, sincos(@as(f128, 0x1p+120)).cosx);
    try std.testing.expectEqual(0x9.f963166f215eb89381836cafaa3p-4, sincos(@as(f128, 0x8p+124)).sinx);
    try std.testing.expectEqual(0xc.82b8ec98b5e62fcf0b09fd10eb3p-4, sincos(@as(f128, 0x8p+124)).cosx);
    try std.testing.expectEqual(0xc.6fa5c5665984d8892761be1537b8p-8, sincos(@as(f128, 0xf.ffffcp+124)).sinx);
    try std.testing.expectEqual(0xf.fb2a030c5ae20bdfe29fda198eap-4, sincos(@as(f128, 0xf.ffffcp+124)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba906cee446be04998p-4, sincos(@as(f128, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, sincos(@as(f128, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x7.f13d78eabb76b8a986d98d6703e8p-4, sincos(@as(f128, 0x4p+48)).sinx);
    try std.testing.expectEqual(0xd.e3b88804f00552d6baba709471d8p-4, sincos(@as(f128, 0x4p+48)).cosx);
    try std.testing.expectEqual(-0xf.c777c6b36a750a5fdeb8807a156p-4, sincos(@as(f128, 0x1p+28)).sinx);
    try std.testing.expectEqual(-0x2.a62ba8824e5bcb065f5f3b8e4f58p-4, sincos(@as(f128, 0x1p+28)).cosx);
    try std.testing.expectEqual(0x8.599b32844aba906cee446be04998p-4, sincos(@as(f128, -0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, sincos(@as(f128, -0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(-0xf.e00885042dd770c93962abdb61f8p-4, sincos(@as(f128, -0x3.3de320f6be87ep+1020)).sinx);
    try std.testing.expectEqual(-0x1.febbf9949ecc133623bb8c8c5a27p-4, sincos(@as(f128, -0x3.3de320f6be87ep+1020)).cosx);
    try std.testing.expectEqual(0xc.773a2eac3000ddec0c69e7ddef68p-4, sincos(@as(f128, 0xe.9f1e6p+112)).sinx);
    try std.testing.expectEqual(-0xa.07bd3ab53ab9710f3445538de8fp-4, sincos(@as(f128, 0xe.9f1e6p+112)).cosx);
    try std.testing.expectEqual(0x7.76d600e031521b7cc3cd579a135p-4, sincos(@as(f128, 0xe.9f1e5p+112)).sinx);
    try std.testing.expectEqual(0xe.26f8af8333f9270e9c3e9f64f94p-4, sincos(@as(f128, 0xe.9f1e5p+112)).cosx);
    try std.testing.expectEqual(0xf.dfffd7bde0fb4ec139784e3b799p-4, sincos(@as(f128, 0xe.9f1e5bc3bb88p+112)).sinx);
    // try std.testing.expectEqual(0x1.ff01000c9ae73630add558c936b5p-4, sincos(@as(f128, 0xe.9f1e5bc3bb88p+112)).cosx);
    try std.testing.expectEqual(-0x1.ffb679ba994b76173f9040637ff9p-4, sincos(@as(f128, 0x4.7857dp+68)).sinx);
    try std.testing.expectEqual(-0xf.dfe902135fc1c18492e869a3f8a8p-4, sincos(@as(f128, 0x4.7857dp+68)).cosx);
    // try std.testing.expectEqual(-0x1.fecaff6878a10ce5d42fde40e7p-4, sincos(@as(f128, 0x6.287cdp+0)).sinx);
    try std.testing.expectEqual(0xf.e006a1ad17db69b4cedfec37da98p-4, sincos(@as(f128, 0x6.287cdp+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb7e68ad6e9c3f77c1915bc919p-4, sincos(@as(f128, 0x6.287cc8p+0)).sinx);
    try std.testing.expectEqual(0xf.e00691b6bde4251c3b197736a7p-4, sincos(@as(f128, 0x6.287cc8p+0)).cosx);
    // try std.testing.expectEqual(-0x1.fecb772e1b8300e5ab16d9008ea9p-4, sincos(@as(f128, 0x6.287cc8749213p+0)).sinx);
    try std.testing.expectEqual(0xf.e006929f558dbe67de4071414d98p-4, sincos(@as(f128, 0x6.287cc8749213p+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb772e1b86f8e74fbeae63ee4cp-4, sincos(@as(f128, 0x6.287cc8749212cp+0)).sinx);
    try std.testing.expectEqual(0xf.e006929f558d3eb50074ea600e6p-4, sincos(@as(f128, 0x6.287cc8749212cp+0)).cosx);
    try std.testing.expectEqual(-0x1.fecb772e1b848bca4e961470b22p-4, sincos(@as(f128, 0x6.287cc8749212e72p+0)).sinx);
    try std.testing.expectEqual(0xf.e006929f558d8cc5d90bd654dfbp-4, sincos(@as(f128, 0x6.287cc8749212e72p+0)).cosx);
    try std.testing.expectEqual(-0xd.8f691a7a95425ffcb89dc2b97cep-4, sincos(@as(f128, -0x1.02e34cp+0)).sinx);
    try std.testing.expectEqual(0x8.7e0ea4db2f488671c85df7208968p-4, sincos(@as(f128, -0x1.02e34cp+0)).cosx);
    try std.testing.expectEqual(-0x8.3bee07bc9076424bef274717106p-4, sincos(@as(f128, 0xf.f0274p+4)).sinx);
    try std.testing.expectEqual(-0xd.b7f5359babdb66be8d0cd3e293e8p-4, sincos(@as(f128, 0xf.f0274p+4)).cosx);
    // try std.testing.expectEqual(0x1.ffc6da9f1ffed895f9fa424ba91p-4, sincos(@as(f128, 0x3.042d88p+0)).sinx);
    try std.testing.expectEqual(-0xf.dfe6f2169e24f276e8027d91ba9p-4, sincos(@as(f128, 0x3.042d88p+0)).cosx);
    try std.testing.expectEqual(-0x8.599b32844aba906cee446be04998p-4, sincos(@as(f128, 0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, sincos(@as(f128, 0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(0x1.452fc98b34e96b61139b09a7c84ap-8, sincos(@as(f128, 0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba9e038d934070f138p-4, sincos(@as(f128, 0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(0xf.dfd9d4b6d0e5f7b9650cab0f5438p-4, sincos(@as(f128, 0xf.fffffffffffffffp+16380)).sinx);
    try std.testing.expectEqual(-0x2.002ef4018753d50b7a7f6bc3f5bap-4, sincos(@as(f128, 0xf.fffffffffffffffp+16380)).cosx);
    try std.testing.expectEqual(0xf.3b0b11ed85b7fe43d110251580b8p-4, sincos(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)).sinx);
    try std.testing.expectEqual(-0x4.e6dc95fb529bc365f472e4fbc1f8p-4, sincos(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)).cosx);
    try std.testing.expectEqual(-0xe.f1a3e1dc468a921dddb4e37fbe6p-4, sincos(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)).sinx);
    try std.testing.expectEqual(-0x5.b773d971a848e75c230605526978p-4, sincos(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)).cosx);
    try std.testing.expectEqual(0x8.599b32844aba906cee446be04998p-4, sincos(@as(f128, -0xf.fffffp+124)).sinx);
    try std.testing.expectEqual(0xd.a5f963cdefe6d529f6b6009fb2fp-4, sincos(@as(f128, -0xf.fffffp+124)).cosx);
    try std.testing.expectEqual(-0x1.452fc98b34e96b61139b09a7c84ap-8, sincos(@as(f128, -0xf.ffffffffffff8p+1020)).sinx);
    try std.testing.expectEqual(-0xf.fff31767d5ba9e038d934070f138p-4, sincos(@as(f128, -0xf.ffffffffffff8p+1020)).cosx);
    try std.testing.expectEqual(-0xf.dfd9d4b6d0e5f7b9650cab0f5438p-4, sincos(@as(f128, -0xf.fffffffffffffffp+16380)).sinx);
    try std.testing.expectEqual(-0x2.002ef4018753d50b7a7f6bc3f5bap-4, sincos(@as(f128, -0xf.fffffffffffffffp+16380)).cosx);
    try std.testing.expectEqual(-0xf.3b0b11ed85b7fe43d110251580b8p-4, sincos(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)).sinx);
    try std.testing.expectEqual(-0x4.e6dc95fb529bc365f472e4fbc1f8p-4, sincos(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)).cosx);
    try std.testing.expectEqual(0xe.f1a3e1dc468a921dddb4e37fbe6p-4, sincos(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)).sinx);
    try std.testing.expectEqual(-0x5.b773d971a848e75c230605526978p-4, sincos(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)).cosx);
    try std.testing.expectEqual(0x4p-128, sincos(@as(f128, 0x4p-128)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x4p-128)).cosx);
    try std.testing.expectEqual(0x4p-1024, sincos(@as(f128, 0x4p-1024)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x4p-1024)).cosx);
    try std.testing.expectEqual(0x4p-16384, sincos(@as(f128, 0x4p-16384)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x4p-16384)).cosx);
    try std.testing.expectEqual(0x2p-16384, sincos(@as(f128, 0x2p-16384)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x2p-16384)).cosx);
    try std.testing.expectEqual(0x8p-972, sincos(@as(f128, 0x8p-972)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x8p-972)).cosx);
    try std.testing.expectEqual(-0x4p-128, sincos(@as(f128, -0x4p-128)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x4p-128)).cosx);
    try std.testing.expectEqual(-0x4p-1024, sincos(@as(f128, -0x4p-1024)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x4p-1024)).cosx);
    try std.testing.expectEqual(-0x4p-16384, sincos(@as(f128, -0x4p-16384)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x4p-16384)).cosx);
    try std.testing.expectEqual(-0x2p-16384, sincos(@as(f128, -0x2p-16384)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x2p-16384)).cosx);
    try std.testing.expectEqual(-0x8p-972, sincos(@as(f128, -0x8p-972)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x8p-972)).cosx);
    try std.testing.expectEqual(0x8p-152, sincos(@as(f128, 0x8p-152)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x8p-152)).cosx);
    try std.testing.expectEqual(0x4p-1076, sincos(@as(f128, 0x4p-1076)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x4p-1076)).cosx);
    try std.testing.expectEqual(0x8p-16448, sincos(@as(f128, 0x8p-16448)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x8p-16448)).cosx);
    try std.testing.expectEqual(0x4p-16448, sincos(@as(f128, 0x4p-16448)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x4p-16448)).cosx);
    try std.testing.expectEqual(0x4p-16496, sincos(@as(f128, 0x4p-16496)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, 0x4p-16496)).cosx);
    try std.testing.expectEqual(-0x8p-152, sincos(@as(f128, -0x8p-152)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x8p-152)).cosx);
    try std.testing.expectEqual(-0x4p-1076, sincos(@as(f128, -0x4p-1076)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x4p-1076)).cosx);
    try std.testing.expectEqual(-0x8p-16448, sincos(@as(f128, -0x8p-16448)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x8p-16448)).cosx);
    try std.testing.expectEqual(-0x4p-16448, sincos(@as(f128, -0x4p-16448)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x4p-16448)).cosx);
    try std.testing.expectEqual(-0x4p-16496, sincos(@as(f128, -0x4p-16496)).sinx);
    try std.testing.expectEqual(0x1p+0, sincos(@as(f128, -0x4p-16496)).cosx);
    try std.testing.expectEqual(0xf.fa2add3e58948d10238cc27b562p-4, sincos(@as(f128, 0x1.8475e6p+0)).sinx);
    try std.testing.expectEqual(0xd.a8263394be6d0e58c1c35a8a3bap-8, sincos(@as(f128, 0x1.8475e6p+0)).cosx);
    try std.testing.expectEqual(0xf.fa2adb8953ae26229c919ec8f6cp-4, sincos(@as(f128, 0x1.8475e4p+0)).sinx);
    try std.testing.expectEqual(0xd.a82832da19f9891d9762fa659ff8p-8, sincos(@as(f128, 0x1.8475e4p+0)).cosx);
    try std.testing.expectEqual(0xf.fa2adcf9ea83dbdd053ee455ea7p-4, sincos(@as(f128, 0x1.8475e5afd4481p+0)).sinx);
    try std.testing.expectEqual(0xd.a82683a33cbebfffffffa2966878p-8, sincos(@as(f128, 0x1.8475e5afd4481p+0)).cosx);
}
