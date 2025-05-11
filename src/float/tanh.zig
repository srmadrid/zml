const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
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

            return @mulAdd(f32, -x, float.abs(x), x);
        }
        const x2: f32 = x * x;
        return @mulAdd(f32, x, -0x1.555556p-2 * x2, x);
    }
    if ((ux << 1) > (0x41102cb3 << 1))
        return float.copysign(@as(f32, 1), x) - float.copysign(@as(f32, 0x1p-25), x);
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
            if (float.abs(x) < std.math.floatMin(f64)) {
                const vx: f64 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            return x * (1 + x); // tanh(small) = small
        }
        if (ix >= 0x3ff00000) { // |x|>=1
            const t: f64 = float.expm1(2 * float.abs(x));
            z = 1 - 2 / (t + 2);
        } else {
            const t: f64 = float.expm1(-2 * float.abs(x));
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
            if (float.abs(x) < std.math.floatMin(f128)) {
                const vx: f128 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            return x * (1 + tiny); // tanh(small) = small
        }
        u.w0 = ix; // Absolute value of x.
        if (ix >= 0x3fff0000) { // |x| >= 1
            const t: f128 = float.expm1(2 * @as(f128, @bitCast(u)));
            z = 1 - 2 / (t + 2);
        } else {
            const t: f128 = float.expm1(-2 * @as(f128, @bitCast(u)));
            z = -t / (t + 2);
        }
        // |x| > 40, return +-1
    } else {
        z = 1 - tiny; // raised inexact flag
    }
    return if ((jx & 0x80000000) != 0) -z else z;
}
