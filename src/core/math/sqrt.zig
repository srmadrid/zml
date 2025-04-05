const cast = @import("../types.zig").cast;
const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const dla = @import("dla.zig");
const root = @import("root.zig");
const EnsureFloat = types.EnsureFloat;

const has_hardware_sqrt16: bool = @import("builtin").target.abi.float() == .hard;
const has_hardware_sqrt32: bool = false;
const has_hardware_sqrt64: bool = false;
const has_hardware_sqrt80: bool = false;
const has_hardware_sqrt128: bool = false;

pub fn sqrt(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return sqrt(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, sqrt32(cast(f32, x))),
                f32 => return sqrt32(x),
                f64 => return sqrt64(x),
                f80 => return cast(f80, sqrt128(cast(f128, x, .{})), .{}),
                f128 => return sqrt128(x),
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn sqrt32(x: f32) f32 {
    var ix: i32 = @bitCast(x);
    const sign: u32 = 0x80000000;

    // take care of Inf and NaN
    if ((ix & 0x7f800000) == 0x7f800000) {
        return x * x + x; // sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN
    }

    // take care of zero
    if (ix <= 0) {
        if ((ix & (~sign)) == 0) {
            return x; // sqrt(+-0) = +-0
        } else if (ix < 0) {
            return (x - x) / (x - x); // sqrt(-ve) = sNaN
        }
    }

    // normalize x
    var m: i32 = (ix >> 23);
    if (m == 0) { // subnormal x
        var i: i32 = 0;
        while ((ix & 0x00800000) == 0) {
            ix <<= 1;

            i += 1;
        }

        m -= i - 1;
    }

    m -= 127; // unbias exponent
    ix = (ix & 0x007fffff) | 0x00800000;

    if ((m & 1) != 0) // odd m, double x to make it even
        ix += ix;

    m >>= 1; // m = [m/2]

    // generate sqrt(x) bit by bit
    ix += ix;
    var s: i32 = 0;
    var q: i32 = 0; // q = sqrt(x)
    var r: i32 = 0x01000000; // r = moving bit from right to left

    while (r != 0) {
        const t: i32 = s + r;
        if (t <= ix) {
            s = t + r;
            ix -= t;
            q += r;
        }

        ix += ix;
        r >>= 1;
    }

    // use floating add to find out rounding direction
    var z: f32 = 0;
    if (ix != 0) {
        z = 0x1p0 - 0x1.4484cp-100; // trigger inexact flag.
        if (z >= 0x1p0) {
            z = 0x1p0 + 0x1.4484cp-100;
            if (z > 0x1p0) {
                q += 2;
            } else {
                q += (q & 1);
            }
        }
    }
    ix = (q >> 1) + 0x3f000000;
    ix += (m << 23);
    return @bitCast(ix);
}

fn sqrt64(x: f64) f64 {
    const rt0: f64 = 9.99999999859990725855365213134618e-01;
    const rt1: f64 = 4.99999999495955425917856814202739e-01;
    const rt2: f64 = 3.75017500867345182581453026130850e-01;
    const rt3: f64 = 3.12523626554518656309172508769531e-01;
    const big: f64 = 134217728.0;

    var a: [2]i32 = @bitCast(x);
    const k: i32 = a[root.HIGH_HALF];
    a[root.HIGH_HALF] = (k & 0x001fffff) | 0x3fe00000;
    var t: f64 = root.inroot[@intCast((k & 0x001fffff) >> 14)];
    const s: f64 = @bitCast(a);
    //----------------- 2^-1022  <= | x |< 2^1024  -----------------
    var c: [2]i32 = .{ 0, 0 };
    if (k > 0x000fffff and k < 0x7ff00000) {
        var y: f64 = 1.0 - t * (t * s);
        t = t * (rt0 + y * (rt1 + y * (rt2 + y * rt3)));
        c[root.HIGH_HALF] = 0x20000000 + ((k & 0x7fe00000) >> 1);
        y = t * s;
        const hy: f64 = (y + big) - big;
        const del: f64 = 0.5 * t * ((s - hy * hy) - (y - hy) * (y + hy));
        var res: f64 = y + del;
        var ret: f64 = undefined;
        if (res == (res + 1.002 * ((y - res) + del))) {
            ret = res * @as(f64, @bitCast(c));
        } else {
            const res1: f64 = res + 1.5 * ((y - res) + del);
            var z: f64 = undefined;
            var zz: f64 = undefined;
            dla.emulv(res, res1, &z, &zz); // (z+zz)=res*res1
            res = if (((z - s) + zz) < 0) @max(res, res1) else @min(res, res1);
            ret = res * @as(f64, @bitCast(c));
        }
        std.mem.doNotOptimizeAway(ret);
        const dret: f64 = x / ret;
        if (dret != ret) {
            const force_inexact: f64 = 1.0 / 3.0;
            std.mem.doNotOptimizeAway(force_inexact);
            // The square root is inexact, ret is the round-to-nearest
            // value which may need adjusting for other rounding
            // modes.
        }
        // Otherwise (x / ret == ret), either the square root was exact or
        // the division was inexact.
        return ret;
    } else {
        if ((k & 0x7ff00000) == 0x7ff00000)
            return x * x + x; // sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN

        if (x == 0)
            return x; // sqrt(+0)=+0, sqrt(-0)=-0

        if (k < 0)
            return (x - x) / (x - x); // sqrt(-ve)=sNaN

        return 0x1p-256 * sqrt64(x * 0x1p512);
    }
}

// Very close to perfect but fails too many tests, although by very little margins
fn sqrt128(x: f128) f128 {
    const t16382: f128 = 0x1p16382;
    const tm8191: f128 = 0x1p-8191;
    const big: f128 = 0x1p120;
    const big1: f128 = 0x1p120 + 1;

    const a: [2]u64 = @bitCast(x);
    const hi_index: u32 = if (@import("builtin").cpu.arch.endian() == .big) 0 else 1;
    const lo_index: u32 = 1 - hi_index;
    // Extract the components of the binary128 number
    const sign_exp_hi: u64 = a[hi_index];
    // Get the exponent bits (15 bits)
    const exp: i64 = @intCast((sign_exp_hi >> 48) & 0x7fff);
    // Mask for the high part mantissa bits
    const mantissa_hi: u64 = sign_exp_hi & 0x0000ffffffffffff;
    // Low part of mantissa
    const mantissa_lo: u64 = a[lo_index];
    //----------------- 2^-16382 <= |x| < 2^16384 -----------------
    if (exp > 0 and exp < 0x7fff) {
        // Normal number
        if (x < 0) return (big1 - big1) / (big - big); // Return NaN for negative input
        // Construct a normalized value with exponent 0x3ffe (bias 0x3fff + 0)
        const adjusted_exp: i64 = exp - 0x3fff;
        const is_odd: bool = (adjusted_exp & 1) != 0;
        const exponent_part: u64 = if (is_odd) 0x4000 else 0x3fff;
        var norm: [2]u64 = undefined;
        norm[hi_index] = (sign_exp_hi & 0x8000000000000000) | (exponent_part << 48) | mantissa_hi;
        norm[lo_index] = mantissa_lo;
        const s: f128 = @bitCast(norm);
        // Compute initial approximation using double precision sqrt
        const d_approx: f64 = sqrt(cast(f64, s, .{}));
        var i: f128 = cast(f128, d_approx, .{});
        // Set the exponent of the result: for sqrt, divide exponent by 2
        const new_exp: i64 = if (is_odd) ((adjusted_exp - 1) >> 1) + 0x3fff else (adjusted_exp >> 1) + 0x3fff;
        var c: [2]u64 = undefined;
        c[hi_index] = (a[hi_index] & 0x8000000000000000) | (@as(u64, @intCast(new_exp & 0x7fff)) << 48);
        c[lo_index] = 0;
        // Newton-Raphson iterations for binary128 precision
        const t: f128 = 0.5 * (i + s / i);
        i = 0.5 * (t + s / t);
        i = 0.5 * (i + s / i); // Extra iteration for quad precision
        return @as(f128, @bitCast(c)) * i;
    } else {
        // Handle special cases
        if (exp == 0x7fff) {
            // sqrt(-Inf) = NaN, sqrt(NaN) = NaN, sqrt(+Inf) = +Inf
            return x * x + x;
        }
        if (x == 0) return x; // sqrt(+0) = +0, sqrt(-0) = -0
        if (x < 0) return (big1 - big1) / (big - big); // Return NaN for negative input
        // Subnormal number: scale up and try again
        return tm8191 * sqrt(x * t16382);
    }
}

test sqrt {
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sqrt(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x2.fp+4, sqrt(@as(f32, 0x8.a1p+8)));
    try std.testing.expectEqual(0x2p+0, sqrt(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(0x1.6a09e6p+0, sqrt(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(0x8p-4, sqrt(@as(f32, 0x4p-4)));
    try std.testing.expectEqual(0x5.18p+4, sqrt(@as(f32, 0x1.9f24p+12)));
    try std.testing.expectEqual(0x7.b4p+4, sqrt(@as(f32, 0x3.b569p+12)));
    try std.testing.expectEqual(0xd.db3d7p-4, sqrt(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x9.9c264p+0, sqrt(@as(f32, 0x5.c59ef8p+4)));
    try std.testing.expectEqual(0x9.9c264p+0, sqrt(@as(f32, 0x5.c59efp+4)));
    try std.testing.expectEqual(0x9.a98c4p+0, sqrt(@as(f32, 0x5.d5c268p+4)));
    try std.testing.expectEqual(0x9.a98c3p+0, sqrt(@as(f32, 0x5.d5c26p+4)));
    try std.testing.expectEqual(0x9.c4b21p+0, sqrt(@as(f32, 0x5.f6ba6p+4)));
    try std.testing.expectEqual(0x9.c4b2p+0, sqrt(@as(f32, 0x5.f6ba58p+4)));
    try std.testing.expectEqual(0x9.cab2dp+0, sqrt(@as(f32, 0x5.fe1118p+4)));
    try std.testing.expectEqual(0x9.cab2dp+0, sqrt(@as(f32, 0x5.fe111p+4)));
    try std.testing.expectEqual(0x9.d755p+0, sqrt(@as(f32, 0x6.0d9198p+4)));
    try std.testing.expectEqual(0x9.d754fp+0, sqrt(@as(f32, 0x6.0d919p+4)));
    try std.testing.expectEqual(0x9.def91p+0, sqrt(@as(f32, 0x6.16fb78p+4)));
    try std.testing.expectEqual(0x9.def9p+0, sqrt(@as(f32, 0x6.16fb7p+4)));
    try std.testing.expectEqual(0x9.dfebfp+0, sqrt(@as(f32, 0x6.18274p+4)));
    try std.testing.expectEqual(0x9.dfebfp+0, sqrt(@as(f32, 0x6.182738p+4)));
    try std.testing.expectEqual(0x9.e3bf7p+0, sqrt(@as(f32, 0x6.1ce128p+4)));
    try std.testing.expectEqual(0x9.e3bf6p+0, sqrt(@as(f32, 0x6.1ce12p+4)));
    try std.testing.expectEqual(0x9.e9d39p+0, sqrt(@as(f32, 0x6.246728p+4)));
    try std.testing.expectEqual(0x9.e9d38p+0, sqrt(@as(f32, 0x6.24672p+4)));
    try std.testing.expectEqual(0x9.f93eap+0, sqrt(@as(f32, 0x6.379128p+4)));
    try std.testing.expectEqual(0x9.f93eap+0, sqrt(@as(f32, 0x6.37912p+4)));
    try std.testing.expectEqual(0xa.074abp+0, sqrt(@as(f32, 0x6.4920a8p+4)));
    try std.testing.expectEqual(0xa.074aap+0, sqrt(@as(f32, 0x6.4920ap+4)));
    try std.testing.expectEqual(0xa.07ca5p+0, sqrt(@as(f32, 0x6.49c0b8p+4)));
    try std.testing.expectEqual(0xa.07ca5p+0, sqrt(@as(f32, 0x6.49c0bp+4)));
    try std.testing.expectEqual(0xa.08ad8p+0, sqrt(@as(f32, 0x6.4add9p+4)));
    try std.testing.expectEqual(0xa.08ad7p+0, sqrt(@as(f32, 0x6.4add88p+4)));
    try std.testing.expectEqual(0xa.0e549p+0, sqrt(@as(f32, 0x6.51f688p+4)));
    try std.testing.expectEqual(0xa.0e549p+0, sqrt(@as(f32, 0x6.51f68p+4)));
    try std.testing.expectEqual(0xa.109f2p+0, sqrt(@as(f32, 0x6.54d828p+4)));
    try std.testing.expectEqual(0xa.109f1p+0, sqrt(@as(f32, 0x6.54d82p+4)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0xf.fffffp-4, sqrt(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0xf.fffffp+60, sqrt(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2p-64, sqrt(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x2.d413ccp-76, sqrt(@as(f32, 0x8p-152)));

    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sqrt(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x2.fp+4, sqrt(@as(f64, 0x8.a1p+8)));
    try std.testing.expectEqual(0x2p+0, sqrt(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(0x1.6a09e667f3bcdp+0, sqrt(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(0x8p-4, sqrt(@as(f64, 0x4p-4)));
    try std.testing.expectEqual(0x5.18p+4, sqrt(@as(f64, 0x1.9f24p+12)));
    try std.testing.expectEqual(0x7.b4p+4, sqrt(@as(f64, 0x3.b569p+12)));
    try std.testing.expectEqual(0xd.db3d742c2655p-4, sqrt(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.ffffffffffff8p+508, sqrt(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffe8p+508, sqrt(@as(f64, 0xf.fffffffffffd8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffd8p+508, sqrt(@as(f64, 0xf.fffffffffffb8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffc8p+508, sqrt(@as(f64, 0xf.fffffffffff98p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffb8p+508, sqrt(@as(f64, 0xf.fffffffffff78p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffa8p+508, sqrt(@as(f64, 0xf.fffffffffff58p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff98p+508, sqrt(@as(f64, 0xf.fffffffffff38p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff88p+508, sqrt(@as(f64, 0xf.fffffffffff18p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff78p+508, sqrt(@as(f64, 0xf.ffffffffffef8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff68p+508, sqrt(@as(f64, 0xf.ffffffffffed8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff58p+508, sqrt(@as(f64, 0xf.ffffffffffeb8p+1020)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000002p-512, sqrt(@as(f64, 0x4.000000000000cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000006p-512, sqrt(@as(f64, 0x4.000000000001cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000000ap-512, sqrt(@as(f64, 0x4.000000000002cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000000ep-512, sqrt(@as(f64, 0x4.000000000003cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000012p-512, sqrt(@as(f64, 0x4.000000000004cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000016p-512, sqrt(@as(f64, 0x4.000000000005cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000001ap-512, sqrt(@as(f64, 0x4.000000000006cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000001ep-512, sqrt(@as(f64, 0x4.000000000007cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000022p-512, sqrt(@as(f64, 0x4.000000000008cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000026p-512, sqrt(@as(f64, 0x4.000000000009cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000002ap-512, sqrt(@as(f64, 0x4.00000000000acp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000002ep-512, sqrt(@as(f64, 0x4.00000000000bcp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000032p-512, sqrt(@as(f64, 0x4.00000000000ccp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000036p-512, sqrt(@as(f64, 0x4.00000000000dcp-1024)));
    try std.testing.expectEqual(0x9.9c2644cd30bb8p+0, sqrt(@as(f64, 0x5.c59ef8p+4)));
    try std.testing.expectEqual(0x9.9c263e244ad48p+0, sqrt(@as(f64, 0x5.c59efp+4)));
    try std.testing.expectEqual(0x9.9c263edb28668p+0, sqrt(@as(f64, 0x5.c59ef0dbaa8ecp+4)));
    try std.testing.expectEqual(0x9.a98c39b89d67p+0, sqrt(@as(f64, 0x5.d5c268p+4)));
    try std.testing.expectEqual(0x9.a98c3318f39a8p+0, sqrt(@as(f64, 0x5.d5c26p+4)));
    try std.testing.expectEqual(0x9.a98c3744dc7f8p+0, sqrt(@as(f64, 0x5.d5c26509ceb5cp+4)));
    try std.testing.expectEqual(0x9.c4b20c8056ad8p+0, sqrt(@as(f64, 0x5.f6ba6p+4)));
    try std.testing.expectEqual(0x9.c4b205f315648p+0, sqrt(@as(f64, 0x5.f6ba58p+4)));
    try std.testing.expectEqual(0x9.c4b207d8c5ba8p+0, sqrt(@as(f64, 0x5.f6ba5a510bf98p+4)));
    try std.testing.expectEqual(0x9.cab2d298bebd8p+0, sqrt(@as(f64, 0x5.fe1118p+4)));
    try std.testing.expectEqual(0x9.cab2cc0f81b98p+0, sqrt(@as(f64, 0x5.fe111p+4)));
    try std.testing.expectEqual(0x9.cab2cf4a334f8p+0, sqrt(@as(f64, 0x5.fe1113f3d9f94p+4)));
    try std.testing.expectEqual(0x9.d754fb02747a8p+0, sqrt(@as(f64, 0x6.0d9198p+4)));
    try std.testing.expectEqual(0x9.d754f4819b76p+0, sqrt(@as(f64, 0x6.0d919p+4)));
    try std.testing.expectEqual(0x9.d754f7f0d1eb8p+0, sqrt(@as(f64, 0x6.0d9194398e95p+4)));
    try std.testing.expectEqual(0x9.def90901b2498p+0, sqrt(@as(f64, 0x6.16fb78p+4)));
    try std.testing.expectEqual(0x9.def90285e1fa8p+0, sqrt(@as(f64, 0x6.16fb7p+4)));
    try std.testing.expectEqual(0x9.def90643382b8p+0, sqrt(@as(f64, 0x6.16fb749d3b76p+4)));
    try std.testing.expectEqual(0x9.dfebf2f55d738p+0, sqrt(@as(f64, 0x6.18274p+4)));
    try std.testing.expectEqual(0x9.dfebec7a2ca38p+0, sqrt(@as(f64, 0x6.182738p+4)));
    try std.testing.expectEqual(0x9.dfebf0a5af558p+0, sqrt(@as(f64, 0x6.18273d25aaddcp+4)));
    try std.testing.expectEqual(0x9.e3bf6a5937aa8p+0, sqrt(@as(f64, 0x6.1ce128p+4)));
    try std.testing.expectEqual(0x9.e3bf63e088bf8p+0, sqrt(@as(f64, 0x6.1ce12p+4)));
    try std.testing.expectEqual(0x9.e3bf69a0e93f8p+0, sqrt(@as(f64, 0x6.1ce1271c28dd4p+4)));
    try std.testing.expectEqual(0x9.e9d38a9ae728p+0, sqrt(@as(f64, 0x6.246728p+4)));
    try std.testing.expectEqual(0x9.e9d384263013p+0, sqrt(@as(f64, 0x6.24672p+4)));
    try std.testing.expectEqual(0x9.e9d3889f74178p+0, sqrt(@as(f64, 0x6.2467258b2eab8p+4)));
    try std.testing.expectEqual(0x9.f93ea4af11dp+0, sqrt(@as(f64, 0x6.379128p+4)));
    try std.testing.expectEqual(0x9.f93e9e4455bp+0, sqrt(@as(f64, 0x6.37912p+4)));
    try std.testing.expectEqual(0x9.f93ea24110618p+0, sqrt(@as(f64, 0x6.379124f88b718p+4)));
    try std.testing.expectEqual(0xa.074aaaa4fe728p+0, sqrt(@as(f64, 0x6.4920a8p+4)));
    try std.testing.expectEqual(0xa.074aa4433f5p+0, sqrt(@as(f64, 0x6.4920ap+4)));
    try std.testing.expectEqual(0xa.074aa97761478p+0, sqrt(@as(f64, 0x6.4920a685e8a2p+4)));
    try std.testing.expectEqual(0xa.07ca572a4cf8p+0, sqrt(@as(f64, 0x6.49c0b8p+4)));
    try std.testing.expectEqual(0xa.07ca50c8df108p+0, sqrt(@as(f64, 0x6.49c0bp+4)));
    try std.testing.expectEqual(0xa.07ca537efeef8p+0, sqrt(@as(f64, 0x6.49c0b3664bc48p+4)));
    try std.testing.expectEqual(0xa.08ad7c223e15p+0, sqrt(@as(f64, 0x6.4add9p+4)));
    try std.testing.expectEqual(0xa.08ad75c1609fp+0, sqrt(@as(f64, 0x6.4add88p+4)));
    try std.testing.expectEqual(0xa.08ad7b0a34af8p+0, sqrt(@as(f64, 0x6.4add8ea0c47f4p+4)));
    try std.testing.expectEqual(0xa.0e548e9fa1b48p+0, sqrt(@as(f64, 0x6.51f688p+4)));
    try std.testing.expectEqual(0xa.0e5488425a1a8p+0, sqrt(@as(f64, 0x6.51f68p+4)));
    try std.testing.expectEqual(0xa.0e5488814a078p+0, sqrt(@as(f64, 0x6.51f6804f1ca4cp+4)));
    try std.testing.expectEqual(0xa.109f1c7a3781p+0, sqrt(@as(f64, 0x6.54d828p+4)));
    try std.testing.expectEqual(0xa.109f161e62cc8p+0, sqrt(@as(f64, 0x6.54d82p+4)));
    try std.testing.expectEqual(0xa.109f19a63bd68p+0, sqrt(@as(f64, 0x6.54d8247125348p+4)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x6.a62e23c62d1b4p-512, sqrt(@as(f64, 0x2.c36098cp-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xa.0b15721eac108p-512, sqrt(@as(f64, 0x6.4de27c4p-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xa.37b39b75a25d8p-512, sqrt(@as(f64, 0x6.86626dp-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xb.3d1b76201dd78p-512, sqrt(@as(f64, 0x7.e4ef24p-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xc.51155b6e7f708p-512, sqrt(@as(f64, 0x9.7b3af18p-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xe.720c54b67ac48p-512, sqrt(@as(f64, 0xd.0ac284p-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xf.2f78e32ee6758p-512, sqrt(@as(f64, 0xe.698f83cp-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.18a9f607e1701p-508, sqrt(@as(f64, 0x1.33b43b08p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.324402a00b45fp-508, sqrt(@as(f64, 0x1.6e66a858p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.3c212046bfdffp-508, sqrt(@as(f64, 0x1.8661cbf8p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.510681b939931p-508, sqrt(@as(f64, 0x1.bbb221b4p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.5461e59227ab5p-508, sqrt(@as(f64, 0x1.c4942f3cp-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.5cf7b0f78d3afp-508, sqrt(@as(f64, 0x1.dbb258c8p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.a31ab946d340bp-508, sqrt(@as(f64, 0x2.ae207d48p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.cad197e28e85bp-508, sqrt(@as(f64, 0x3.36529f1p-1016)));
    try std.testing.expectEqual(0x1.000000ffffff8p+0, sqrt(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xf.fffff7fffffep-4, sqrt(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-4, sqrt(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.ffffffffffff8p+508, sqrt(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2p-64, sqrt(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x2p-512, sqrt(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0xb.504f333f9de68p-488, sqrt(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x2.d413cccfe779ap-76, sqrt(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x8p-540, sqrt(@as(f64, 0x4p-1076)));

    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sqrt(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x2.fp+4, sqrt(@as(f80, 0x8.a1p+8)));
    try std.testing.expectEqual(0x2p+0, sqrt(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908p+0, sqrt(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(0x8p-4, sqrt(@as(f80, 0x4p-4)));
    try std.testing.expectEqual(0x5.18p+4, sqrt(@as(f80, 0x1.9f24p+12)));
    try std.testing.expectEqual(0x7.b4p+4, sqrt(@as(f80, 0x3.b569p+12)));
    try std.testing.expectEqual(0xd.db3d742c265539ep-4, sqrt(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.ffffffffffffcp+508, sqrt(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffecp+508, sqrt(@as(f80, 0xf.fffffffffffd8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffdcp+508, sqrt(@as(f80, 0xf.fffffffffffb8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffccp+508, sqrt(@as(f80, 0xf.fffffffffff98p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffbcp+508, sqrt(@as(f80, 0xf.fffffffffff78p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffacp+508, sqrt(@as(f80, 0xf.fffffffffff58p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff9cp+508, sqrt(@as(f80, 0xf.fffffffffff38p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff8cp+508, sqrt(@as(f80, 0xf.fffffffffff18p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff7cp+508, sqrt(@as(f80, 0xf.ffffffffffef8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff6cp+508, sqrt(@as(f80, 0xf.ffffffffffed8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff5cp+508, sqrt(@as(f80, 0xf.ffffffffffeb8p+1020)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000003p-512, sqrt(@as(f80, 0x4.000000000000cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000007p-512, sqrt(@as(f80, 0x4.000000000001cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000000bp-512, sqrt(@as(f80, 0x4.000000000002cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000000fp-512, sqrt(@as(f80, 0x4.000000000003cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000013p-512, sqrt(@as(f80, 0x4.000000000004cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000017p-512, sqrt(@as(f80, 0x4.000000000005cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000001bp-512, sqrt(@as(f80, 0x4.000000000006cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000001fp-512, sqrt(@as(f80, 0x4.000000000007cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000023p-512, sqrt(@as(f80, 0x4.000000000008cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000027p-512, sqrt(@as(f80, 0x4.000000000009cp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000002bp-512, sqrt(@as(f80, 0x4.00000000000acp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000002fp-512, sqrt(@as(f80, 0x4.00000000000bcp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000033p-512, sqrt(@as(f80, 0x4.00000000000ccp-1024)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000037p-512, sqrt(@as(f80, 0x4.00000000000dcp-1024)));
    try std.testing.expectEqual(0x9.9c2644cd30bbba3p+0, sqrt(@as(f80, 0x5.c59ef8p+4)));
    try std.testing.expectEqual(0x9.9c263e244ad482bp+0, sqrt(@as(f80, 0x5.c59efp+4)));
    try std.testing.expectEqual(0x9.9c263edb28664p+0, sqrt(@as(f80, 0x5.c59ef0dbaa8ecp+4)));
    try std.testing.expectEqual(0x9.a98c39b89d67019p+0, sqrt(@as(f80, 0x5.d5c268p+4)));
    try std.testing.expectEqual(0x9.a98c3318f39aaffp+0, sqrt(@as(f80, 0x5.d5c26p+4)));
    try std.testing.expectEqual(0x9.a98c3744dc7fcp+0, sqrt(@as(f80, 0x5.d5c26509ceb5cp+4)));
    try std.testing.expectEqual(0x9.c4b20c8056ad6c7p+0, sqrt(@as(f80, 0x5.f6ba6p+4)));
    try std.testing.expectEqual(0x9.c4b205f315647f8p+0, sqrt(@as(f80, 0x5.f6ba58p+4)));
    try std.testing.expectEqual(0x9.c4b207d8c5bacp+0, sqrt(@as(f80, 0x5.f6ba5a510bf98p+4)));
    try std.testing.expectEqual(0x9.cab2d298bebd96dp+0, sqrt(@as(f80, 0x5.fe1118p+4)));
    try std.testing.expectEqual(0x9.cab2cc0f81b9a91p+0, sqrt(@as(f80, 0x5.fe111p+4)));
    try std.testing.expectEqual(0x9.cab2cf4a334f4p+0, sqrt(@as(f80, 0x5.fe1113f3d9f94p+4)));
    try std.testing.expectEqual(0x9.d754fb02747aa6dp+0, sqrt(@as(f80, 0x6.0d9198p+4)));
    try std.testing.expectEqual(0x9.d754f4819b75ebfp+0, sqrt(@as(f80, 0x6.0d919p+4)));
    try std.testing.expectEqual(0x9.d754f7f0d1eb4p+0, sqrt(@as(f80, 0x6.0d9194398e95p+4)));
    try std.testing.expectEqual(0x9.def90901b2497a9p+0, sqrt(@as(f80, 0x6.16fb78p+4)));
    try std.testing.expectEqual(0x9.def90285e1fa539p+0, sqrt(@as(f80, 0x6.16fb7p+4)));
    try std.testing.expectEqual(0x9.def90643382b4p+0, sqrt(@as(f80, 0x6.16fb749d3b76p+4)));
    try std.testing.expectEqual(0x9.dfebf2f55d73acp+0, sqrt(@as(f80, 0x6.18274p+4)));
    try std.testing.expectEqual(0x9.dfebec7a2ca3556p+0, sqrt(@as(f80, 0x6.182738p+4)));
    try std.testing.expectEqual(0x9.dfebf0a5af554p+0, sqrt(@as(f80, 0x6.18273d25aaddcp+4)));
    try std.testing.expectEqual(0x9.e3bf6a5937aa469p+0, sqrt(@as(f80, 0x6.1ce128p+4)));
    try std.testing.expectEqual(0x9.e3bf63e088bf98ap+0, sqrt(@as(f80, 0x6.1ce12p+4)));
    try std.testing.expectEqual(0x9.e3bf69a0e93fcp+0, sqrt(@as(f80, 0x6.1ce1271c28dd4p+4)));
    try std.testing.expectEqual(0x9.e9d38a9ae7283d9p+0, sqrt(@as(f80, 0x6.246728p+4)));
    try std.testing.expectEqual(0x9.e9d384263012d63p+0, sqrt(@as(f80, 0x6.24672p+4)));
    try std.testing.expectEqual(0x9.e9d3889f74174p+0, sqrt(@as(f80, 0x6.2467258b2eab8p+4)));
    try std.testing.expectEqual(0x9.f93ea4af11cfcc5p+0, sqrt(@as(f80, 0x6.379128p+4)));
    try std.testing.expectEqual(0x9.f93e9e4455afe27p+0, sqrt(@as(f80, 0x6.37912p+4)));
    try std.testing.expectEqual(0x9.f93ea2411061cp+0, sqrt(@as(f80, 0x6.379124f88b718p+4)));
    try std.testing.expectEqual(0xa.074aaaa4fe728dep+0, sqrt(@as(f80, 0x6.4920a8p+4)));
    try std.testing.expectEqual(0xa.074aa4433f5023ap+0, sqrt(@as(f80, 0x6.4920ap+4)));
    try std.testing.expectEqual(0xa.074aa9776147cp+0, sqrt(@as(f80, 0x6.4920a685e8a2p+4)));
    try std.testing.expectEqual(0xa.07ca572a4cf7c2ap+0, sqrt(@as(f80, 0x6.49c0b8p+4)));
    try std.testing.expectEqual(0xa.07ca50c8df10bebp+0, sqrt(@as(f80, 0x6.49c0bp+4)));
    try std.testing.expectEqual(0xa.07ca537efeef4p+0, sqrt(@as(f80, 0x6.49c0b3664bc48p+4)));
    try std.testing.expectEqual(0xa.08ad7c223e15144p+0, sqrt(@as(f80, 0x6.4add9p+4)));
    try std.testing.expectEqual(0xa.08ad75c1609f282p+0, sqrt(@as(f80, 0x6.4add88p+4)));
    try std.testing.expectEqual(0xa.08ad7b0a34afcp+0, sqrt(@as(f80, 0x6.4add8ea0c47f4p+4)));
    try std.testing.expectEqual(0xa.0e548e9fa1b46efp+0, sqrt(@as(f80, 0x6.51f688p+4)));
    try std.testing.expectEqual(0xa.0e5488425a1a91ep+0, sqrt(@as(f80, 0x6.51f68p+4)));
    try std.testing.expectEqual(0xa.0e5488814a074p+0, sqrt(@as(f80, 0x6.51f6804f1ca4cp+4)));
    try std.testing.expectEqual(0xa.109f1c7a3780ff9p+0, sqrt(@as(f80, 0x6.54d828p+4)));
    try std.testing.expectEqual(0xa.109f161e62ccb8ep+0, sqrt(@as(f80, 0x6.54d82p+4)));
    try std.testing.expectEqual(0xa.109f19a63bd64p+0, sqrt(@as(f80, 0x6.54d8247125348p+4)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x6.a62e23c62d1b6p-512, sqrt(@as(f80, 0x2.c36098cp-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xa.0b15721eac10cp-512, sqrt(@as(f80, 0x6.4de27c4p-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xa.37b39b75a25dcp-512, sqrt(@as(f80, 0x6.86626dp-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xb.3d1b76201dd74p-512, sqrt(@as(f80, 0x7.e4ef24p-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xc.51155b6e7f70cp-512, sqrt(@as(f80, 0x9.7b3af18p-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xe.720c54b67ac4cp-512, sqrt(@as(f80, 0xd.0ac284p-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xf.2f78e32ee675cp-512, sqrt(@as(f80, 0xe.698f83cp-1020)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.18a9f607e17018p-508, sqrt(@as(f80, 0x1.33b43b08p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.324402a00b45e8p-508, sqrt(@as(f80, 0x1.6e66a858p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.3c212046bfdfe8p-508, sqrt(@as(f80, 0x1.8661cbf8p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.510681b9399308p-508, sqrt(@as(f80, 0x1.bbb221b4p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.5461e59227ab58p-508, sqrt(@as(f80, 0x1.c4942f3cp-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.5cf7b0f78d3ae8p-508, sqrt(@as(f80, 0x1.dbb258c8p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.a31ab946d340a8p-508, sqrt(@as(f80, 0x2.ae207d48p-1016)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.cad197e28e85a8p-508, sqrt(@as(f80, 0x3.36529f1p-1016)));
    try std.testing.expectEqual(0x1.000000ffffff8p+0, sqrt(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.00000000000008p+0, sqrt(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xf.fffff7fffffep-4, sqrt(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0xf.ffffffffffffcp-4, sqrt(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0xf.fffff7fffffep+60, sqrt(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.ffffffffffffcp+508, sqrt(@as(f80, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(0xf.fffffffffffffffp+8188, sqrt(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2p-64, sqrt(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x2p-512, sqrt(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x2p-8192, sqrt(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908p-8192, sqrt(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0xb.504f333f9de6484p-488, sqrt(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-76, sqrt(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x8p-540, sqrt(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x2.d413cccfe779921p-8224, sqrt(@as(f80, 0x8p-16448)));

    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, sqrt(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x2.fp+4, sqrt(@as(f128, 0x8.a1p+8)));
    try std.testing.expectEqual(0x2p+0, sqrt(@as(f128, 0x4p+0)));
    // try std.testing.expectEqual(0x1.6a09e667f3bcc908b2fb1366ea95p+0, sqrt(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(0x8p-4, sqrt(@as(f128, 0x4p-4)));
    try std.testing.expectEqual(0x5.18p+4, sqrt(@as(f128, 0x1.9f24p+12)));
    try std.testing.expectEqual(0x7.b4p+4, sqrt(@as(f128, 0x3.b569p+12)));
    try std.testing.expectEqual(0xd.db3d742c265539d92ba16b83c5cp-4, sqrt(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.ffffffffffffbfffffffffffff8p+508, sqrt(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffebffffffffffff38p+508, sqrt(@as(f128, 0xf.fffffffffffd8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffdbfffffffffffd78p+508, sqrt(@as(f128, 0xf.fffffffffffb8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffcbfffffffffffab8p+508, sqrt(@as(f128, 0xf.fffffffffff98p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffbbfffffffffff6f8p+508, sqrt(@as(f128, 0xf.fffffffffff78p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffffabfffffffffff238p+508, sqrt(@as(f128, 0xf.fffffffffff58p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff9bffffffffffec78p+508, sqrt(@as(f128, 0xf.fffffffffff38p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff8bffffffffffe5b8p+508, sqrt(@as(f128, 0xf.fffffffffff18p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff7bffffffffffddf8p+508, sqrt(@as(f128, 0xf.ffffffffffef8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff6bffffffffffd538p+508, sqrt(@as(f128, 0xf.ffffffffffed8p+1020)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff5bffffffffffcb78p+508, sqrt(@as(f128, 0xf.ffffffffffeb8p+1020)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000002ffffffffffffdcp-512, sqrt(@as(f128, 0x4.000000000000cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000006ffffffffffff3cp-512, sqrt(@as(f128, 0x4.000000000001cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000000afffffffffffe1cp-512, sqrt(@as(f128, 0x4.000000000002cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000000efffffffffffc7cp-512, sqrt(@as(f128, 0x4.000000000003cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000012fffffffffffa5cp-512, sqrt(@as(f128, 0x4.000000000004cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000016fffffffffff7bcp-512, sqrt(@as(f128, 0x4.000000000005cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000001afffffffffff49cp-512, sqrt(@as(f128, 0x4.000000000006cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000001efffffffffff0fcp-512, sqrt(@as(f128, 0x4.000000000007cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000022ffffffffffecdcp-512, sqrt(@as(f128, 0x4.000000000008cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000026ffffffffffe83cp-512, sqrt(@as(f128, 0x4.000000000009cp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000002affffffffffe31cp-512, sqrt(@as(f128, 0x4.00000000000acp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.000000000002effffffffffdd7cp-512, sqrt(@as(f128, 0x4.00000000000bcp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000032ffffffffffd75cp-512, sqrt(@as(f128, 0x4.00000000000ccp-1024)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x2.0000000000036ffffffffffd0bcp-512, sqrt(@as(f128, 0x4.00000000000dcp-1024)));
    try std.testing.expectEqual(0x9.9c2644cd30bbba2af9770d1a6adp+0, sqrt(@as(f128, 0x5.c59ef8p+4)));
    try std.testing.expectEqual(0x9.9c263e244ad482ae439b6ebb7e28p+0, sqrt(@as(f128, 0x5.c59efp+4)));
    try std.testing.expectEqual(0x9.9c263edb286640061abcbed744a8p+0, sqrt(@as(f128, 0x5.c59ef0dbaa8ecp+4)));
    // try std.testing.expectEqual(0x9.a98c39b89d67018ef271d70ac418p+0, sqrt(@as(f128, 0x5.d5c268p+4)));
    try std.testing.expectEqual(0x9.a98c3318f39aaff1fdbe5f0520ep+0, sqrt(@as(f128, 0x5.d5c26p+4)));
    // try std.testing.expectEqual(0x9.a98c3744dc7fbff970920ee38d88p+0, sqrt(@as(f128, 0x5.d5c26509ceb5cp+4)));
    try std.testing.expectEqual(0x9.c4b20c8056ad6c75de5e878cf82p+0, sqrt(@as(f128, 0x5.f6ba6p+4)));
    // try std.testing.expectEqual(0x9.c4b205f315647f84aa6427849be8p+0, sqrt(@as(f128, 0x5.f6ba58p+4)));
    try std.testing.expectEqual(0x9.c4b207d8c5babfff1af9eff60338p+0, sqrt(@as(f128, 0x5.f6ba5a510bf98p+4)));
    try std.testing.expectEqual(0x9.cab2d298bebd96cb80a18484978p+0, sqrt(@as(f128, 0x5.fe1118p+4)));
    // try std.testing.expectEqual(0x9.cab2cc0f81b9a9129e80476516c8p+0, sqrt(@as(f128, 0x5.fe111p+4)));
    try std.testing.expectEqual(0x9.cab2cf4a334f40040a75564c1b1p+0, sqrt(@as(f128, 0x5.fe1113f3d9f94p+4)));
    try std.testing.expectEqual(0x9.d754fb02747aa6d3627024158ap+0, sqrt(@as(f128, 0x6.0d9198p+4)));
    try std.testing.expectEqual(0x9.d754f4819b75ebe926d2e4b1b1ap+0, sqrt(@as(f128, 0x6.0d919p+4)));
    // try std.testing.expectEqual(0x9.d754f7f0d1eb40067ddd2f83bff8p+0, sqrt(@as(f128, 0x6.0d9194398e95p+4)));
    // try std.testing.expectEqual(0x9.def90901b2497a8da8b1b1590948p+0, sqrt(@as(f128, 0x6.16fb78p+4)));
    try std.testing.expectEqual(0x9.def90285e1fa5395a1f3d3cfbbcp+0, sqrt(@as(f128, 0x6.16fb7p+4)));
    try std.testing.expectEqual(0x9.def90643382b40078c6ec55aae08p+0, sqrt(@as(f128, 0x6.16fb749d3b76p+4)));
    try std.testing.expectEqual(0x9.dfebf2f55d73ac019895a5c2c9cp+0, sqrt(@as(f128, 0x6.18274p+4)));
    try std.testing.expectEqual(0x9.dfebec7a2ca355606b63be1f301p+0, sqrt(@as(f128, 0x6.182738p+4)));
    // try std.testing.expectEqual(0x9.dfebf0a5af5540000f31060e2b68p+0, sqrt(@as(f128, 0x6.18273d25aaddcp+4)));
    // try std.testing.expectEqual(0x9.e3bf6a5937aa46890d478fffcc08p+0, sqrt(@as(f128, 0x6.1ce128p+4)));
    // try std.testing.expectEqual(0x9.e3bf63e088bf9899430188403ea8p+0, sqrt(@as(f128, 0x6.1ce12p+4)));
    try std.testing.expectEqual(0x9.e3bf69a0e93fbffea021a294fdfp+0, sqrt(@as(f128, 0x6.1ce1271c28dd4p+4)));
    try std.testing.expectEqual(0x9.e9d38a9ae7283d96dd13217515fp+0, sqrt(@as(f128, 0x6.246728p+4)));
    // try std.testing.expectEqual(0x9.e9d384263012d635564d99c20bf8p+0, sqrt(@as(f128, 0x6.24672p+4)));
    // try std.testing.expectEqual(0x9.e9d3889f7417400693816cdbf0a8p+0, sqrt(@as(f128, 0x6.2467258b2eab8p+4)));
    try std.testing.expectEqual(0x9.f93ea4af11cfcc4c7b3ad927d8cp+0, sqrt(@as(f128, 0x6.379128p+4)));
    // try std.testing.expectEqual(0x9.f93e9e4455afe2757febddb767c8p+0, sqrt(@as(f128, 0x6.37912p+4)));
    // try std.testing.expectEqual(0x9.f93ea2411061bffe7e7f96636678p+0, sqrt(@as(f128, 0x6.379124f88b718p+4)));
    try std.testing.expectEqual(0xa.074aaaa4fe728de305512ee33bbp+0, sqrt(@as(f128, 0x6.4920a8p+4)));
    try std.testing.expectEqual(0xa.074aa4433f5023a592a623bd7a5p+0, sqrt(@as(f128, 0x6.4920ap+4)));
    try std.testing.expectEqual(0xa.074aa9776147bffed6639c1d4e8p+0, sqrt(@as(f128, 0x6.4920a685e8a2p+4)));
    try std.testing.expectEqual(0xa.07ca572a4cf7c2a718b5e0e26dd8p+0, sqrt(@as(f128, 0x6.49c0b8p+4)));
    try std.testing.expectEqual(0xa.07ca50c8df10beb4ab7d2c2955ep+0, sqrt(@as(f128, 0x6.49c0bp+4)));
    try std.testing.expectEqual(0xa.07ca537efeef4007f8bedbd94e3p+0, sqrt(@as(f128, 0x6.49c0b3664bc48p+4)));
    try std.testing.expectEqual(0xa.08ad7c223e151446a914e7db42ap+0, sqrt(@as(f128, 0x6.4add9p+4)));
    // try std.testing.expectEqual(0xa.08ad75c1609f28197eb9e111fca8p+0, sqrt(@as(f128, 0x6.4add88p+4)));
    // try std.testing.expectEqual(0xa.08ad7b0a34afbff8b29545cf6968p+0, sqrt(@as(f128, 0x6.4add8ea0c47f4p+4)));
    try std.testing.expectEqual(0xa.0e548e9fa1b46eed2440fd35521p+0, sqrt(@as(f128, 0x6.51f688p+4)));
    try std.testing.expectEqual(0xa.0e5488425a1a91e23c39fe025e1p+0, sqrt(@as(f128, 0x6.51f68p+4)));
    try std.testing.expectEqual(0xa.0e5488814a074003b5a5ffdb32dp+0, sqrt(@as(f128, 0x6.51f6804f1ca4cp+4)));
    try std.testing.expectEqual(0xa.109f1c7a3780ff90f6697001cb1p+0, sqrt(@as(f128, 0x6.54d828p+4)));
    try std.testing.expectEqual(0xa.109f161e62ccb8e65a0922adbf88p+0, sqrt(@as(f128, 0x6.54d82p+4)));
    try std.testing.expectEqual(0xa.109f19a63bd64002ee4fb2c64dd8p+0, sqrt(@as(f128, 0x6.54d8247125348p+4)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x6.a62e23c62d1b5ffe5c81a90f554p-512, sqrt(@as(f128, 0x2.c36098cp-1020)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0xa.0b15721eac10bffdd9746fa70ca8p-512, sqrt(@as(f128, 0x6.4de27c4p-1020)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0xa.37b39b75a25dbffc951409f30528p-512, sqrt(@as(f128, 0x6.86626dp-1020)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0xb.3d1b76201dd740065804ad1abea8p-512, sqrt(@as(f128, 0x7.e4ef24p-1020)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0xc.51155b6e7f70bffcf0277d2b5618p-512, sqrt(@as(f128, 0x9.7b3af18p-1020)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0xe.720c54b67ac4bfff3dde8c941bf8p-512, sqrt(@as(f128, 0xd.0ac284p-1020)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0xf.2f78e32ee675bffe6792b3ce66fp-512, sqrt(@as(f128, 0xe.698f83cp-1020)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.18a9f607e17017ff715a73e157afp-508, sqrt(@as(f128, 0x1.33b43b08p-1016)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.324402a00b45e800a761e004b92ap-508, sqrt(@as(f128, 0x1.6e66a858p-1016)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0x1.3c212046bfdfe8004a6543b0a63bp-508, sqrt(@as(f128, 0x1.8661cbf8p-1016)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.510681b93993080072e1891cced7p-508, sqrt(@as(f128, 0x1.bbb221b4p-1016)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0x1.5461e59227ab57ff0ef1d1ea7cc1p-508, sqrt(@as(f128, 0x1.c4942f3cp-1016)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.5cf7b0f78d3ae8008b2b0f38c32ep-508, sqrt(@as(f128, 0x1.dbb258c8p-1016)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.a31ab946d340a800ad52925a9b6p-508, sqrt(@as(f128, 0x2.ae207d48p-1016)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, sqrt(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.cad197e28e85a800936f11336851p-508, sqrt(@as(f128, 0x3.36529f1p-1016)));
    try std.testing.expectEqual(0x1.000000ffffff8000007fffff6p+0, sqrt(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.00000000000007ffffffffffffep+0, sqrt(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1p+0, sqrt(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p-4, sqrt(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0xf.ffffffffffffbfffffffffffff8p-4, sqrt(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0xf.fffff7fffffdfffffeffffff6p+60, sqrt(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.ffffffffffffbfffffffffffff8p+508, sqrt(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.fffffffffffffff8p+8188, sqrt(@as(f128, 0xf.fffffffffffffffp+16380)));
    // try std.testing.expectEqual(0xf.fffffffffffffffffffffffffff8p+8188, sqrt(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0xf.ffffffffffffdffffffffffffdep+508, sqrt(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x2p-64, sqrt(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x2p-512, sqrt(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x2p-8192, sqrt(@as(f128, 0x4p-16384)));
    // try std.testing.expectEqual(0x1.6a09e667f3bcc908b2fb1366ea95p-8192, sqrt(@as(f128, 0x2p-16384)));
    // try std.testing.expectEqual(0xb.504f333f9de6484597d89b3754a8p-488, sqrt(@as(f128, 0x8p-972)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-76, sqrt(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x8p-540, sqrt(@as(f128, 0x4p-1076)));
    // try std.testing.expectEqual(0x2.d413cccfe779921165f626cdd52ap-8224, sqrt(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x2p-8224, sqrt(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x2p-8248, sqrt(@as(f128, 0x4p-16496)));
}
