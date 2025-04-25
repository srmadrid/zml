const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const roundeven = @import("roundeven.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn tanh(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return tanh(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, tanh32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_tanhf.c
                    return tanh32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_tanh.c
                    return tanh64(x);
                },
                f80 => return cast(f80, tanh128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_tanhl.c
                    return tanh128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn tanh32(x: f32) f32 {
    const z: f64 = cast(f64, x, .{});
    const ux: u32 = @bitCast(x);
    const e: i32 = @bitCast((ux >> 23) & 0xff);
    if (e == 0xff) {
        @branchHint(.unlikely);
        if ((ux << 9) != 0)
            return x + x; // x = nan

        const ir: [2]f32 = .{ 1, -1 };
        return ir[ux >> 31]; // x = +-inf
    }

    if (e < 115) {
        @branchHint(.unlikely);
        if (e < 102) {
            @branchHint(.unlikely);
            if ((ux << 1) == 0) {
                @branchHint(.unlikely);
                return x;
            }

            return @mulAdd(f32, -x, math.abs(x), x);
        }
        const x2: f32 = x * x;
        return @mulAdd(f32, x, -0x1.555556p-2 * x2, x);
    }
    if ((ux << 1) > (0x41102cb3 << 1))
        return math.copysign(@as(f32, 1), x) - math.copysign(@as(f32, 0x1p-25), x);
    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    const z8: f64 = z4 * z4;
    const cn: [8]f64 = .{
        0x1p+0,                0x1.30877b8b72d33p-3,  0x1.694aa09ae9e5ep-8,
        0x1.4101377abb729p-14, 0x1.e0392b1db0018p-22, 0x1.2533756e546f7p-30,
        0x1.d62e5abe6ae8ap-41, 0x1.b06be534182dep-54,
    };
    const cd: [8]f64 = .{
        0x1p+0,                0x1.ed99131b0ebeap-2,  0x1.0d27ed6c95a69p-5,
        0x1.7cbdaca0e9fccp-11, 0x1.b4e60b892578ep-18, 0x1.a6f707c5c71abp-26,
        0x1.35a8b6e2cd94cp-35, 0x1.ca8230677aa01p-47,
    };
    var n0: f64 = cn[0] + z2 * cn[1];
    const n2: f64 = cn[2] + z2 * cn[3];
    var n4: f64 = cn[4] + z2 * cn[5];
    const n6: f64 = cn[6] + z2 * cn[7];
    n0 += z4 * n2;
    n4 += z4 * n6;
    n0 += z8 * n4;
    var d0: f64 = cd[0] + z2 * cd[1];
    const d2: f64 = cd[2] + z2 * cd[3];
    var d4: f64 = cd[4] + z2 * cd[5];
    const d6: f64 = cd[6] + z2 * cd[7];
    d0 += z4 * d2;
    d4 += z4 * d6;
    d0 += z8 * d4;
    const r: f64 = z * n0 / d0;
    return cast(f32, r, .{});
}

fn tanh64(x: f64) f64 {
    const tiny: f64 = 1.0e-300;

    // High word of |x|.
    var jx: i32 = undefined;
    var lx: i32 = undefined;
    dbl64.extractWords(&jx, &lx, x);
    const ix: i32 = jx & 0x7fffffff;

    // x is INF or NaN
    if (ix >= 0x7ff00000) {
        if (jx >= 0) {
            return 1 / x + 1; // tanh(+-inf)=+-1
        } else {
            return 1 / x - 1; // tanh(NaN) = NaN
        }
    }

    // |x| < 22
    var z: f64 = undefined;
    if (ix < 0x40360000) { // |x|<22
        if ((ix | lx) == 0)
            return x; // x == +-0
        if (ix < 0x3c800000) { // |x|<2**-55
            if (math.abs(x) < std.math.floatMin(f64)) {
                const vx: f64 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            return x * (1 + x); // tanh(small) = small
        }
        if (ix >= 0x3ff00000) { // |x|>=1
            const t: f64 = math.expm1(2 * math.abs(x));
            z = 1 - 2 / (t + 2);
        } else {
            const t: f64 = math.expm1(-2 * math.abs(x));
            z = -t / (t + 2);
        }
        // |x| > 22, return +-1
    } else {
        z = 1 - tiny; // raised inexact flag
    }
    return if (jx >= 0) z else -z;
}

fn tanh128(x: f128) f128 {
    const tiny: f128 = 1.0e-4900;

    // Words of |x|.
    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const jx: u32 = u.w0;
    const ix: u32 = jx & 0x7fffffff;
    // x is INF or NaN
    if (ix >= 0x7fff0000) {
        // for NaN it's not important which branch: tanhl(NaN) = NaN
        if ((jx & 0x80000000) != 0) {
            return 1 / x - 1; // tanhl(-inf)= -1;
        } else {
            return 1 / x + 1; // tanhl(+inf)=+1
        }
    }

    // |x| < 40
    var z: f128 = undefined;
    if (ix < 0x40044000) {
        if (@as(f128, @bitCast(u)) == 0)
            return x; // x == +- 0
        if (ix < 0x3fc60000) { // |x| < 2^-57
            if (math.abs(x) < std.math.floatMin(f128)) {
                const vx: f128 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            return x * (1 + tiny); // tanh(small) = small
        }
        u.w0 = ix; // Absolute value of x.
        if (ix >= 0x3fff0000) { // |x| >= 1
            const t: f128 = math.expm1(2 * @as(f128, @bitCast(u)));
            z = 1 - 2 / (t + 2);
        } else {
            const t: f128 = math.expm1(-2 * @as(f128, @bitCast(u)));
            z = -t / (t + 2);
        }
        // |x| > 40, return +-1
    } else {
        z = 1 - tiny; // raised inexact flag
    }
    return if ((jx & 0x80000000) != 0) -z else z;
}

test tanh {
    try std.testing.expectEqual(0x0p+0, tanh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, tanh(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0xa.2991fp-4, tanh(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(-0xa.2991fp-4, tanh(@as(f32, -0xcp-4)));
    try std.testing.expectEqual(0xc.2f7d6p-4, tanh(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0xc.2f7d6p-4, tanh(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0xf.6ca83p-4, tanh(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(-0xf.6ca83p-4, tanh(@as(f32, -0x2p+0)));
    try std.testing.expectEqual(0xf.ebbe9p-4, tanh(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(-0xf.ebbe9p-4, tanh(@as(f32, -0x3p+0)));
    try std.testing.expectEqual(0xf.fd40cp-4, tanh(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(-0xf.fd40cp-4, tanh(@as(f32, -0x4p+0)));
    try std.testing.expectEqual(0xf.ffa0dp-4, tanh(@as(f32, 0x5p+0)));
    try std.testing.expectEqual(-0xf.ffa0dp-4, tanh(@as(f32, -0x5p+0)));
    try std.testing.expectEqual(0xf.fff32p-4, tanh(@as(f32, 0x6p+0)));
    try std.testing.expectEqual(-0xf.fff32p-4, tanh(@as(f32, -0x6p+0)));
    try std.testing.expectEqual(0xf.fffe4p-4, tanh(@as(f32, 0x7p+0)));
    try std.testing.expectEqual(-0xf.fffe4p-4, tanh(@as(f32, -0x7p+0)));
    try std.testing.expectEqual(0xf.ffffcp-4, tanh(@as(f32, 0x8p+0)));
    try std.testing.expectEqual(-0xf.ffffcp-4, tanh(@as(f32, -0x8p+0)));
    try std.testing.expectEqual(0xf.fffffp-4, tanh(@as(f32, 0x9p+0)));
    try std.testing.expectEqual(-0xf.fffffp-4, tanh(@as(f32, -0x9p+0)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0xap+0)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0xap+0)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0xfp+0)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0xfp+0)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0x1.4p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0x1.4p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0x1.6p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0x1.6p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0x1.9p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0x1.9p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0x1.ep+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0x1.ep+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0x2.3p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0x2.3p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0x2.8p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0x2.8p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0x2.dp+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0x2.dp+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0x3.2p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0x3.2p+4)));
    try std.testing.expectEqual(0x8p-60, tanh(@as(f32, 0x8p-60)));
    try std.testing.expectEqual(0xb.78df8p-4, tanh(@as(f32, 0xe.6c659p-4)));
    try std.testing.expectEqual(0x7.fa4a2p-4, tanh(@as(f32, 0x8.c259ep-4)));
    try std.testing.expectEqual(0x6.080bfp-4, tanh(@as(f32, 0x6.5821dp-4)));
    try std.testing.expectEqual(0x7.c5731p-4, tanh(@as(f32, 0x8.7c9e5p-4)));
    try std.testing.expectEqual(-0x3.a55fc8p-4, tanh(@as(f32, -0x3.b60d7cp-4)));
    try std.testing.expectEqual(0x7.2d063p-4, tanh(@as(f32, 0x7.b9985p-4)));
    try std.testing.expectEqual(0x7.19c548p-4, tanh(@as(f32, 0x7.a18e8p-4)));
    try std.testing.expectEqual(-0x2.5c12e8p-4, tanh(@as(f32, -0x2.6082fp-4)));
    try std.testing.expectEqual(0xe.05031p-16, tanh(@as(f32, 0xe.05031p-16)));
    try std.testing.expectEqual(0x3.b66d3cp-4, tanh(@as(f32, 0x3.c80eacp-4)));
    try std.testing.expectEqual(0x3.b66d38p-4, tanh(@as(f32, 0x3.c80ea8p-4)));
    try std.testing.expectEqual(0x1.fe4f3ep-4, tanh(@as(f32, 0x2.00f988p-4)));
    try std.testing.expectEqual(0x1.fe4f3ap-4, tanh(@as(f32, 0x2.00f984p-4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0xe.9e035p+0)));
    try std.testing.expectEqual(-0x3.af99fp-4, tanh(@as(f32, -0x3.c0d8b4p-4)));
    try std.testing.expectEqual(-0x3.af99f4p-4, tanh(@as(f32, -0x3.c0d8b8p-4)));
    try std.testing.expectEqual(-0x3.24bf1p-4, tanh(@as(f32, -0x3.2f59p-4)));
    try std.testing.expectEqual(0x2.deea8p-4, tanh(@as(f32, 0x2.e6f54cp-4)));
    try std.testing.expectEqual(0x3.2e7fbcp-4, tanh(@as(f32, 0x3.397f3p-4)));
    try std.testing.expectEqual(0x3.2e7fbcp-4, tanh(@as(f32, 0x3.397f2cp-4)));
    try std.testing.expectEqual(0x7.96e928p-4, tanh(@as(f32, 0x8.4024bp-4)));
    try std.testing.expectEqual(0x7.96e918p-4, tanh(@as(f32, 0x8.4024ap-4)));
    try std.testing.expectEqual(0x7.ff5568p-8, tanh(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x3.ffffecp-12, tanh(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x2p-16, tanh(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0x1p-20, tanh(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0x8p-28, tanh(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x4p-32, tanh(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, tanh(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, tanh(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, tanh(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, tanh(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, tanh(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, tanh(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, tanh(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tanh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tanh(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x3.a8babp-4, tanh(@as(f32, 0x3.b9979cp-4)));
    try std.testing.expectEqual(0x3.a8baacp-4, tanh(@as(f32, 0x3.b99798p-4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x4p-128, tanh(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-128, tanh(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, tanh(@as(f32, -0x8p-152)));

    try std.testing.expectEqual(0x0p+0, tanh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, tanh(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0xa.2991f2a97914p-4, tanh(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(-0xa.2991f2a97914p-4, tanh(@as(f64, -0xcp-4)));
    try std.testing.expectEqual(0xc.2f7d5a8a79cap-4, tanh(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0xc.2f7d5a8a79cap-4, tanh(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0xf.6ca82f0de1eap-4, tanh(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(-0xf.6ca82f0de1eap-4, tanh(@as(f64, -0x2p+0)));
    try std.testing.expectEqual(0xf.ebbe888d058p-4, tanh(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(-0xf.ebbe888d058p-4, tanh(@as(f64, -0x3p+0)));
    try std.testing.expectEqual(0xf.fd40b84505a1p-4, tanh(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(-0xf.fd40b84505a1p-4, tanh(@as(f64, -0x4p+0)));
    try std.testing.expectEqual(0xf.ffa0cb346f888p-4, tanh(@as(f64, 0x5p+0)));
    try std.testing.expectEqual(-0xf.ffa0cb346f888p-4, tanh(@as(f64, -0x5p+0)));
    try std.testing.expectEqual(0xf.fff31d5f129ep-4, tanh(@as(f64, 0x6p+0)));
    try std.testing.expectEqual(-0xf.fff31d5f129ep-4, tanh(@as(f64, -0x6p+0)));
    try std.testing.expectEqual(0xf.fffe4193a879p-4, tanh(@as(f64, 0x7p+0)));
    try std.testing.expectEqual(-0xf.fffe4193a879p-4, tanh(@as(f64, -0x7p+0)));
    try std.testing.expectEqual(0xf.ffffc39548fcp-4, tanh(@as(f64, 0x8p+0)));
    try std.testing.expectEqual(-0xf.ffffc39548fcp-4, tanh(@as(f64, -0x8p+0)));
    try std.testing.expectEqual(0xf.fffff7d2cebcp-4, tanh(@as(f64, 0x9p+0)));
    try std.testing.expectEqual(-0xf.fffff7d2cebcp-4, tanh(@as(f64, -0x9p+0)));
    try std.testing.expectEqual(0xf.fffffee4b79a8p-4, tanh(@as(f64, 0xap+0)));
    try std.testing.expectEqual(-0xf.fffffee4b79a8p-4, tanh(@as(f64, -0xap+0)));
    try std.testing.expectEqual(0xf.fffffffffcb5p-4, tanh(@as(f64, 0xfp+0)));
    try std.testing.expectEqual(-0xf.fffffffffcb5p-4, tanh(@as(f64, -0xfp+0)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0x1.4p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0x1.4p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0x1.6p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0x1.6p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0x1.9p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0x1.9p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0x1.ep+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0x1.ep+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0x2.3p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0x2.3p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0x2.8p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0x2.8p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0x2.dp+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0x2.dp+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0x3.2p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0x3.2p+4)));
    try std.testing.expectEqual(0x8p-60, tanh(@as(f64, 0x8p-60)));
    try std.testing.expectEqual(0xb.78df781e11d8p-4, tanh(@as(f64, 0xe.6c659p-4)));
    // try std.testing.expectEqual(0x7.fa4a1eea64fa4p-4, tanh(@as(f64, 0x8.c259ep-4)));
    // try std.testing.expectEqual(0x6.080bf03812d8p-4, tanh(@as(f64, 0x6.5821dp-4)));
    try std.testing.expectEqual(0x7.c57313d93519cp-4, tanh(@as(f64, 0x8.7c9e5p-4)));
    // try std.testing.expectEqual(-0x3.a55fc883707acp-4, tanh(@as(f64, -0x3.b60d7cp-4)));
    try std.testing.expectEqual(0x7.2d06324738d24p-4, tanh(@as(f64, 0x7.b9985p-4)));
    // try std.testing.expectEqual(0x7.19c5470dc5d6cp-4, tanh(@as(f64, 0x7.a18e8p-4)));
    try std.testing.expectEqual(-0x2.5c12e9588a796p-4, tanh(@as(f64, -0x2.6082fp-4)));
    // try std.testing.expectEqual(0xe.05030c697d9e8p-16, tanh(@as(f64, 0xe.05031p-16)));
    try std.testing.expectEqual(0x3.b66d3ac34ff94p-4, tanh(@as(f64, 0x3.c80eacp-4)));
    try std.testing.expectEqual(0x3.b66d36fa72348p-4, tanh(@as(f64, 0x3.c80ea8p-4)));
    try std.testing.expectEqual(0x3.b66d39531e604p-4, tanh(@as(f64, 0x3.c80eaa7adaa3p-4)));
    // try std.testing.expectEqual(0x1.fe4f3d0dd83fbp-4, tanh(@as(f64, 0x2.00f988p-4)));
    try std.testing.expectEqual(0x1.fe4f391dbd3edp-4, tanh(@as(f64, 0x2.00f984p-4)));
    try std.testing.expectEqual(0x1.fe4f3a8e05153p-4, tanh(@as(f64, 0x2.00f9857616524p-4)));
    try std.testing.expectEqual(-0xf.fffffffff8eb8p-4, tanh(@as(f64, -0xe.9e035p+0)));
    try std.testing.expectEqual(-0x3.af99f04902f54p-4, tanh(@as(f64, -0x3.c0d8b4p-4)));
    // try std.testing.expectEqual(-0x3.af99f412aab74p-4, tanh(@as(f64, -0x3.c0d8b8p-4)));
    // try std.testing.expectEqual(-0x3.af99f183b9d72p-4, tanh(@as(f64, -0x3.c0d8b54c5a488p-4)));
    try std.testing.expectEqual(-0x3.24bf114777f9p-4, tanh(@as(f64, -0x3.2f59p-4)));
    try std.testing.expectEqual(0x2.deea7ea48e5eep-4, tanh(@as(f64, 0x2.e6f54cp-4)));
    try std.testing.expectEqual(0x3.2e7fbdedf6f4ep-4, tanh(@as(f64, 0x3.397f3p-4)));
    try std.testing.expectEqual(0x3.2e7fba1674b72p-4, tanh(@as(f64, 0x3.397f2cp-4)));
    try std.testing.expectEqual(0x3.2e7fbd450f41ep-4, tanh(@as(f64, 0x3.397f2f50241d2p-4)));
    // try std.testing.expectEqual(0x3.2e7fbd450f41cp-4, tanh(@as(f64, 0x3.397f2f50241dp-4)));
    // try std.testing.expectEqual(0x7.96e925f6aa4fcp-4, tanh(@as(f64, 0x8.4024bp-4)));
    // try std.testing.expectEqual(0x7.96e9199045abcp-4, tanh(@as(f64, 0x8.4024ap-4)));
    // try std.testing.expectEqual(0x7.96e91a6be7d9cp-4, tanh(@as(f64, 0x8.4024a11b66108p-4)));
    try std.testing.expectEqual(0x7.96e91a6be7d94p-4, tanh(@as(f64, 0x8.4024a11b661p-4)));
    try std.testing.expectEqual(0x7.ff556664ac778p-8, tanh(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0x3.ffffeaaaab334p-12, tanh(@as(f64, 0x4p-12)));
    // try std.testing.expectEqual(0x1.fffffffd55555p-16, tanh(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffffaaa8p-24, tanh(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0x7.ffffffffffff4p-28, tanh(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x4p-32, tanh(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, tanh(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, tanh(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, tanh(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, tanh(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, tanh(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, tanh(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, tanh(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tanh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, tanh(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tanh(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, tanh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x3.a8baafcd6721cp-4, tanh(@as(f64, 0x3.b9979cp-4)));
    try std.testing.expectEqual(0x3.a8baac02f5784p-4, tanh(@as(f64, 0x3.b99798p-4)));
    // try std.testing.expectEqual(0x3.a8baae38037e4p-4, tanh(@as(f64, 0x3.b9979a543d0fcp-4)));
    // try std.testing.expectEqual(0x3.a8baae38037e2p-4, tanh(@as(f64, 0x3.b9979a543d0fap-4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4p-128, tanh(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, tanh(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, tanh(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, tanh(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, tanh(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, tanh(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, tanh(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-152, tanh(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, tanh(@as(f64, -0x4p-1076)));

    try std.testing.expectEqual(0x0p+0, tanh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, tanh(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0xa.2991f2a9791413ap-4, tanh(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(-0xa.2991f2a9791413ap-4, tanh(@as(f80, -0xcp-4)));
    try std.testing.expectEqual(0xc.2f7d5a8a79ca2acp-4, tanh(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0xc.2f7d5a8a79ca2acp-4, tanh(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0xf.6ca82f0de1e9e9ap-4, tanh(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(-0xf.6ca82f0de1e9e9ap-4, tanh(@as(f80, -0x2p+0)));
    try std.testing.expectEqual(0xf.ebbe888d057ff1p-4, tanh(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(-0xf.ebbe888d057ff1p-4, tanh(@as(f80, -0x3p+0)));
    try std.testing.expectEqual(0xf.fd40b84505a10b4p-4, tanh(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(-0xf.fd40b84505a10b4p-4, tanh(@as(f80, -0x4p+0)));
    try std.testing.expectEqual(0xf.ffa0cb346f889a8p-4, tanh(@as(f80, 0x5p+0)));
    try std.testing.expectEqual(-0xf.ffa0cb346f889a8p-4, tanh(@as(f80, -0x5p+0)));
    try std.testing.expectEqual(0xf.fff31d5f129deeep-4, tanh(@as(f80, 0x6p+0)));
    try std.testing.expectEqual(-0xf.fff31d5f129deeep-4, tanh(@as(f80, -0x6p+0)));
    try std.testing.expectEqual(0xf.fffe4193a878ed7p-4, tanh(@as(f80, 0x7p+0)));
    try std.testing.expectEqual(-0xf.fffe4193a878ed7p-4, tanh(@as(f80, -0x7p+0)));
    try std.testing.expectEqual(0xf.ffffc39548fc348p-4, tanh(@as(f80, 0x8p+0)));
    try std.testing.expectEqual(-0xf.ffffc39548fc348p-4, tanh(@as(f80, -0x8p+0)));
    try std.testing.expectEqual(0xf.fffff7d2cebbe21p-4, tanh(@as(f80, 0x9p+0)));
    try std.testing.expectEqual(-0xf.fffff7d2cebbe21p-4, tanh(@as(f80, -0x9p+0)));
    try std.testing.expectEqual(0xf.fffffee4b79aaa9p-4, tanh(@as(f80, 0xap+0)));
    try std.testing.expectEqual(-0xf.fffffee4b79aaa9p-4, tanh(@as(f80, -0xap+0)));
    try std.testing.expectEqual(0xf.fffffffffcb523ep-4, tanh(@as(f80, 0xfp+0)));
    try std.testing.expectEqual(-0xf.fffffffffcb523ep-4, tanh(@as(f80, -0xfp+0)));
    try std.testing.expectEqual(0xf.fffffffffffff63p-4, tanh(@as(f80, 0x1.4p+4)));
    try std.testing.expectEqual(-0xf.fffffffffffff63p-4, tanh(@as(f80, -0x1.4p+4)));
    try std.testing.expectEqual(0xf.ffffffffffffffdp-4, tanh(@as(f80, 0x1.6p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffffdp-4, tanh(@as(f80, -0x1.6p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0x1.9p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0x1.9p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0x1.ep+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0x1.ep+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0x2.3p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0x2.3p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0x2.8p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0x2.8p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0x2.dp+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0x2.dp+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0x3.2p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0x3.2p+4)));
    try std.testing.expectEqual(0x8p-60, tanh(@as(f80, 0x8p-60)));
    try std.testing.expectEqual(0xb.78df781e11d83e2p-4, tanh(@as(f80, 0xe.6c659p-4)));
    try std.testing.expectEqual(0x7.fa4a1eea64fa2838p-4, tanh(@as(f80, 0x8.c259ep-4)));
    try std.testing.expectEqual(0x6.080bf03812d804f8p-4, tanh(@as(f80, 0x6.5821dp-4)));
    try std.testing.expectEqual(0x7.c57313d93519a7fp-4, tanh(@as(f80, 0x8.7c9e5p-4)));
    try std.testing.expectEqual(-0x3.a55fc883707aca2p-4, tanh(@as(f80, -0x3.b60d7cp-4)));
    try std.testing.expectEqual(0x7.2d06324738d23d5p-4, tanh(@as(f80, 0x7.b9985p-4)));
    try std.testing.expectEqual(0x7.19c5470dc5d6c09p-4, tanh(@as(f80, 0x7.a18e8p-4)));
    try std.testing.expectEqual(-0x2.5c12e9588a795db8p-4, tanh(@as(f80, -0x2.6082fp-4)));
    try std.testing.expectEqual(0xe.05030c697d9e583p-16, tanh(@as(f80, 0xe.05031p-16)));
    try std.testing.expectEqual(0x3.b66d3ac34ff934dp-4, tanh(@as(f80, 0x3.c80eacp-4)));
    try std.testing.expectEqual(0x3.b66d36fa7234779p-4, tanh(@as(f80, 0x3.c80ea8p-4)));
    try std.testing.expectEqual(0x3.b66d39531e6043a8p-4, tanh(@as(f80, 0x3.c80eaa7adaa3p-4)));
    try std.testing.expectEqual(0x1.fe4f3d0dd83fadbp-4, tanh(@as(f80, 0x2.00f988p-4)));
    try std.testing.expectEqual(0x1.fe4f391dbd3ecd72p-4, tanh(@as(f80, 0x2.00f984p-4)));
    try std.testing.expectEqual(0x1.fe4f3a8e0515345p-4, tanh(@as(f80, 0x2.00f9857616524p-4)));
    try std.testing.expectEqual(-0xf.fffffffff8ebcp-4, tanh(@as(f80, -0xe.9e035p+0)));
    try std.testing.expectEqual(-0x3.af99f04902f54a6p-4, tanh(@as(f80, -0x3.c0d8b4p-4)));
    try std.testing.expectEqual(-0x3.af99f412aab73f58p-4, tanh(@as(f80, -0x3.c0d8b8p-4)));
    try std.testing.expectEqual(-0x3.af99f183b9d71e98p-4, tanh(@as(f80, -0x3.c0d8b54c5a488p-4)));
    try std.testing.expectEqual(-0x3.24bf114777f8faf8p-4, tanh(@as(f80, -0x3.2f59p-4)));
    try std.testing.expectEqual(0x2.deea7ea48e5ed334p-4, tanh(@as(f80, 0x2.e6f54cp-4)));
    try std.testing.expectEqual(0x3.2e7fbdedf6f4e468p-4, tanh(@as(f80, 0x3.397f3p-4)));
    try std.testing.expectEqual(0x3.2e7fba1674b721dp-4, tanh(@as(f80, 0x3.397f2cp-4)));
    try std.testing.expectEqual(0x3.2e7fbd450f41db44p-4, tanh(@as(f80, 0x3.397f2f50241d2p-4)));
    try std.testing.expectEqual(0x3.2e7fbd450f41bc84p-4, tanh(@as(f80, 0x3.397f2f50241dp-4)));
    try std.testing.expectEqual(0x3.2e7fbd450f41bf78p-4, tanh(@as(f80, 0x3.397f2f50241d031p-4)));
    try std.testing.expectEqual(0x7.96e925f6aa4fa0fp-4, tanh(@as(f80, 0x8.4024bp-4)));
    try std.testing.expectEqual(0x7.96e9199045abc438p-4, tanh(@as(f80, 0x8.4024ap-4)));
    try std.testing.expectEqual(0x7.96e91a6be7d9c2bp-4, tanh(@as(f80, 0x8.4024a11b66108p-4)));
    try std.testing.expectEqual(0x7.96e91a6be7d95f8p-4, tanh(@as(f80, 0x8.4024a11b661p-4)));
    try std.testing.expectEqual(0x7.96e91a6be7d9af78p-4, tanh(@as(f80, 0x8.4024a11b6610673p-4)));
    try std.testing.expectEqual(0x7.96e91a6be7d9af68p-4, tanh(@as(f80, 0x8.4024a11b6610672p-4)));
    try std.testing.expectEqual(0x7.ff556664ac778a1p-8, tanh(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x3.ffffeaaaab33333p-12, tanh(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffffd5555555ap-16, tanh(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffffaaaaabp-24, tanh(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0x7.ffffffffffff5558p-28, tanh(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x3.ffffffffffffffecp-32, tanh(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, tanh(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, tanh(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, tanh(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, tanh(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, tanh(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, tanh(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, tanh(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tanh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, tanh(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tanh(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, tanh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, tanh(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x3.a8baafcd6721c9d4p-4, tanh(@as(f80, 0x3.b9979cp-4)));
    try std.testing.expectEqual(0x3.a8baac02f578492p-4, tanh(@as(f80, 0x3.b99798p-4)));
    try std.testing.expectEqual(0x3.a8baae38037e31p-4, tanh(@as(f80, 0x3.b9979a543d0fcp-4)));
    try std.testing.expectEqual(0x3.a8baae38037e12acp-4, tanh(@as(f80, 0x3.b9979a543d0fap-4)));
    try std.testing.expectEqual(0x3.a8baae38037e30acp-4, tanh(@as(f80, 0x3.b9979a543d0fbfa8p-4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4p-128, tanh(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, tanh(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, tanh(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, tanh(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, tanh(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, tanh(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, tanh(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, tanh(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, tanh(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, tanh(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, tanh(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, tanh(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-152, tanh(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, tanh(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, tanh(@as(f80, -0x8p-16448)));

    try std.testing.expectEqual(0x0p+0, tanh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, tanh(@as(f128, -0x0p+0)));
    // try std.testing.expectEqual(0xa.2991f2a97914139d5832bf78fb1p-4, tanh(@as(f128, 0xcp-4)));
    // try std.testing.expectEqual(-0xa.2991f2a97914139d5832bf78fb1p-4, tanh(@as(f128, -0xcp-4)));
    // try std.testing.expectEqual(0xc.2f7d5a8a79ca2ac3195f149e2138p-4, tanh(@as(f128, 0x1p+0)));
    // try std.testing.expectEqual(-0xc.2f7d5a8a79ca2ac3195f149e2138p-4, tanh(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0xf.6ca82f0de1e9e99e2197e1f412bp-4, tanh(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(-0xf.6ca82f0de1e9e99e2197e1f412bp-4, tanh(@as(f128, -0x2p+0)));
    try std.testing.expectEqual(0xf.ebbe888d057ff1057854585bfdbp-4, tanh(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(-0xf.ebbe888d057ff1057854585bfdbp-4, tanh(@as(f128, -0x3p+0)));
    try std.testing.expectEqual(0xf.fd40b84505a10b42b92360cee308p-4, tanh(@as(f128, 0x4p+0)));
    try std.testing.expectEqual(-0xf.fd40b84505a10b42b92360cee308p-4, tanh(@as(f128, -0x4p+0)));
    try std.testing.expectEqual(0xf.ffa0cb346f889a800b7186cb573p-4, tanh(@as(f128, 0x5p+0)));
    try std.testing.expectEqual(-0xf.ffa0cb346f889a800b7186cb573p-4, tanh(@as(f128, -0x5p+0)));
    try std.testing.expectEqual(0xf.fff31d5f129deedd313b57265658p-4, tanh(@as(f128, 0x6p+0)));
    try std.testing.expectEqual(-0xf.fff31d5f129deedd313b57265658p-4, tanh(@as(f128, -0x6p+0)));
    try std.testing.expectEqual(0xf.fffe4193a878ed68e8057dafd2dp-4, tanh(@as(f128, 0x7p+0)));
    try std.testing.expectEqual(-0xf.fffe4193a878ed68e8057dafd2dp-4, tanh(@as(f128, -0x7p+0)));
    try std.testing.expectEqual(0xf.ffffc39548fc3487707369d6c578p-4, tanh(@as(f128, 0x8p+0)));
    try std.testing.expectEqual(-0xf.ffffc39548fc3487707369d6c578p-4, tanh(@as(f128, -0x8p+0)));
    try std.testing.expectEqual(0xf.fffff7d2cebbe208a50ed05e717p-4, tanh(@as(f128, 0x9p+0)));
    try std.testing.expectEqual(-0xf.fffff7d2cebbe208a50ed05e717p-4, tanh(@as(f128, -0x9p+0)));
    try std.testing.expectEqual(0xf.fffffee4b79aaa94a2b616896898p-4, tanh(@as(f128, 0xap+0)));
    try std.testing.expectEqual(-0xf.fffffee4b79aaa94a2b616896898p-4, tanh(@as(f128, -0xap+0)));
    try std.testing.expectEqual(0xf.fffffffffcb523e7aa70681dc27p-4, tanh(@as(f128, 0xfp+0)));
    try std.testing.expectEqual(-0xf.fffffffffcb523e7aa70681dc27p-4, tanh(@as(f128, -0xfp+0)));
    try std.testing.expectEqual(0xf.fffffffffffff63436db3272edfp-4, tanh(@as(f128, 0x1.4p+4)));
    try std.testing.expectEqual(-0xf.fffffffffffff63436db3272edfp-4, tanh(@as(f128, -0x1.4p+4)));
    try std.testing.expectEqual(0xf.ffffffffffffffd2117c43d16e28p-4, tanh(@as(f128, 0x1.6p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffffd2117c43d16e28p-4, tanh(@as(f128, -0x1.6p+4)));
    try std.testing.expectEqual(0xf.ffffffffffffffffe2da82ab81f8p-4, tanh(@as(f128, 0x1.9p+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffffffe2da82ab81f8p-4, tanh(@as(f128, -0x1.9p+4)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffa9479b98p-4, tanh(@as(f128, 0x1.ep+4)));
    try std.testing.expectEqual(-0xf.ffffffffffffffffffffa9479b98p-4, tanh(@as(f128, -0x1.ep+4)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffefdf8p-4, tanh(@as(f128, 0x2.3p+4)));
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffefdf8p-4, tanh(@as(f128, -0x2.3p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f128, 0x2.8p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f128, -0x2.8p+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f128, 0x2.dp+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f128, -0x2.dp+4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f128, 0x3.2p+4)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f128, -0x3.2p+4)));
    try std.testing.expectEqual(0x8p-60, tanh(@as(f128, 0x8p-60)));
    try std.testing.expectEqual(0xb.78df781e11d83e198e857e22169p-4, tanh(@as(f128, 0xe.6c659p-4)));
    // try std.testing.expectEqual(0x7.fa4a1eea64fa283aa32d48b14b94p-4, tanh(@as(f128, 0x8.c259ep-4)));
    try std.testing.expectEqual(0x6.080bf03812d804f456a4858a20dp-4, tanh(@as(f128, 0x6.5821dp-4)));
    // try std.testing.expectEqual(0x7.c57313d93519a7edb391a912d4e8p-4, tanh(@as(f128, 0x8.7c9e5p-4)));
    try std.testing.expectEqual(-0x3.a55fc883707aca21b3d3eb4c9496p-4, tanh(@as(f128, -0x3.b60d7cp-4)));
    // try std.testing.expectEqual(0x7.2d06324738d23d4d4328c1a80f9cp-4, tanh(@as(f128, 0x7.b9985p-4)));
    try std.testing.expectEqual(0x7.19c5470dc5d6c0913805237beb5p-4, tanh(@as(f128, 0x7.a18e8p-4)));
    // try std.testing.expectEqual(-0x2.5c12e9588a795db643b503e27bp-4, tanh(@as(f128, -0x2.6082fp-4)));
    // try std.testing.expectEqual(0xe.05030c697d9e582f4a79c88f0198p-16, tanh(@as(f128, 0xe.05031p-16)));
    try std.testing.expectEqual(0x3.b66d3ac34ff934cf70cbc132d382p-4, tanh(@as(f128, 0x3.c80eacp-4)));
    // try std.testing.expectEqual(0x3.b66d36fa7234778e14df5c18c67ap-4, tanh(@as(f128, 0x3.c80ea8p-4)));
    try std.testing.expectEqual(0x3.b66d39531e6043a85263d7aef20ep-4, tanh(@as(f128, 0x3.c80eaa7adaa3p-4)));
    try std.testing.expectEqual(0x1.fe4f3d0dd83fadafe273ab28dd29p-4, tanh(@as(f128, 0x2.00f988p-4)));
    try std.testing.expectEqual(0x1.fe4f391dbd3ecd714619cc709978p-4, tanh(@as(f128, 0x2.00f984p-4)));
    try std.testing.expectEqual(0x1.fe4f3a8e0515344ff794387d92d8p-4, tanh(@as(f128, 0x2.00f9857616524p-4)));
    try std.testing.expectEqual(-0xf.fffffffff8ebbffbf5b020cd6ab8p-4, tanh(@as(f128, -0xe.9e035p+0)));
    try std.testing.expectEqual(-0x3.af99f04902f54a5e1438d014c59p-4, tanh(@as(f128, -0x3.c0d8b4p-4)));
    // try std.testing.expectEqual(-0x3.af99f412aab73f59c1a2be2a32fp-4, tanh(@as(f128, -0x3.c0d8b8p-4)));
    // try std.testing.expectEqual(-0x3.af99f183b9d71e966538c40d38fep-4, tanh(@as(f128, -0x3.c0d8b54c5a488p-4)));
    // try std.testing.expectEqual(-0x3.24bf114777f8faf96902769a0d84p-4, tanh(@as(f128, -0x3.2f59p-4)));
    // try std.testing.expectEqual(0x2.deea7ea48e5ed334e492b456066ep-4, tanh(@as(f128, 0x2.e6f54cp-4)));
    try std.testing.expectEqual(0x3.2e7fbdedf6f4e4677fd41531d3b2p-4, tanh(@as(f128, 0x3.397f3p-4)));
    try std.testing.expectEqual(0x3.2e7fba1674b721d00a6064e53d74p-4, tanh(@as(f128, 0x3.397f2cp-4)));
    // try std.testing.expectEqual(0x3.2e7fbd450f41db420b102ed5c87cp-4, tanh(@as(f128, 0x3.397f2f50241d2p-4)));
    try std.testing.expectEqual(0x3.2e7fbd450f41bc85f9231ae637dp-4, tanh(@as(f128, 0x3.397f2f50241dp-4)));
    try std.testing.expectEqual(0x3.2e7fbd450f41bf76f8da4b4ea52ap-4, tanh(@as(f128, 0x3.397f2f50241d031p-4)));
    try std.testing.expectEqual(0x7.96e925f6aa4fa0f29663e3f79f08p-4, tanh(@as(f128, 0x8.4024bp-4)));
    // try std.testing.expectEqual(0x7.96e9199045abc439fc0595df5b8cp-4, tanh(@as(f128, 0x8.4024ap-4)));
    // try std.testing.expectEqual(0x7.96e91a6be7d9c2af9a5db822e29cp-4, tanh(@as(f128, 0x8.4024a11b66108p-4)));
    // try std.testing.expectEqual(0x7.96e91a6be7d95f7c75164741422p-4, tanh(@as(f128, 0x8.4024a11b661p-4)));
    // try std.testing.expectEqual(0x7.96e91a6be7d9af74d043bee12618p-4, tanh(@as(f128, 0x8.4024a11b6610673p-4)));
    // try std.testing.expectEqual(0x7.96e91a6be7d9af6869df15f309e4p-4, tanh(@as(f128, 0x8.4024a11b6610672p-4)));
    try std.testing.expectEqual(0x7.96e91a6be7d9af71106ffad34228p-4, tanh(@as(f128, 0x8.4024a11b6610672b2982b852e8p-4)));
    try std.testing.expectEqual(0x7.ff556664ac778a0c17f05ce08814p-8, tanh(@as(f128, 0x8p-8)));
    // try std.testing.expectEqual(0x3.ffffeaaaab33332fbefc0623efe6p-12, tanh(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x1.fffffffd5555555999999992b12bp-16, tanh(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0xf.fffffffffaaaaaaaaaacccccccdp-24, tanh(@as(f128, 0x1p-20)));
    // try std.testing.expectEqual(0x7.ffffffffffff5555555555556668p-28, tanh(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x3.ffffffffffffffeaaaaaaaaaaaaap-32, tanh(@as(f128, 0x4p-32)));
    // try std.testing.expectEqual(0x1.fffffffffffffffffd5555555555p-36, tanh(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffaaaaaaaa8p-44, tanh(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffff555554p-48, tanh(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0x3.ffffffffffffffffffffffffeaaap-52, tanh(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x1.fffffffffffffffffffffffffffdp-56, tanh(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, tanh(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, tanh(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tanh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, tanh(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tanh(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, tanh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, tanh(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0x3.a8baafcd6721c9d281b58b34a98ap-4, tanh(@as(f128, 0x3.b9979cp-4)));
    try std.testing.expectEqual(0x3.a8baac02f578491e7a245d92617ap-4, tanh(@as(f128, 0x3.b99798p-4)));
    try std.testing.expectEqual(0x3.a8baae38037e30fe8c8253c51894p-4, tanh(@as(f128, 0x3.b9979a543d0fcp-4)));
    // try std.testing.expectEqual(0x3.a8baae38037e12aaff36972c4886p-4, tanh(@as(f128, 0x3.b9979a543d0fap-4)));
    // try std.testing.expectEqual(0x3.a8baae38037e30ab26bdc37e7458p-4, tanh(@as(f128, 0x3.b9979a543d0fbfa8p-4)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1p+0, tanh(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1p+0, tanh(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4p-128, tanh(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, tanh(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, tanh(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, tanh(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, tanh(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, tanh(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, tanh(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, tanh(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, tanh(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, tanh(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, tanh(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, tanh(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, tanh(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, tanh(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, tanh(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-152, tanh(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, tanh(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, tanh(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, tanh(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, tanh(@as(f128, -0x4p-16496)));
}
