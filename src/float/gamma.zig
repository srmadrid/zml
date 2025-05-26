const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const roundeven = @import("roundeven.zig");
const mul_split = @import("mul_split.zig");
const lgamma = @import("lgamma.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn gamma(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.gamma: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    var local_signgam: i32 = undefined;
    const y: EnsureFloat(@TypeOf(x)) = switch (EnsureFloat(@TypeOf(x))) {
        f16 => scast(f16, gamma_r32(scast(f32, x), &local_signgam)),
        f32 => gamma_r32(scast(f32, x), &local_signgam),
        f64 => gamma_r64(scast(f64, x), &local_signgam),
        f80 => scast(f80, gamma_r128(scast(f128, x), &local_signgam)),
        f128 => gamma_r128(scast(f128, x), &local_signgam),
        else => unreachable,
    };

    return if (local_signgam < 0) -y else y;
}

fn gamma_r32(x: f32, signgamp: *i32) f32 {
    signgamp.* = 1;

    // List of exceptional cases. Each entry contains the 32-bit encoding u of x,
    // a binary32 approximation f of gamma(x), and a correction term df.
    const tb: [10]struct { u: u32, f: f32, df: f32 } = .{
        .{ .u = 0x27de86a9, .f = 0x1.268266p+47, .df = 0x1p22 }, // x = 0x1.bd0d52p-48
        .{ .u = 0x27e05475, .f = 0x1.242422p+47, .df = 0x1p22 }, // x = 0x1.c0a8eap-48
        .{ .u = 0xb63befb3, .f = -0x1.5cb6e4p+18, .df = 0x1p-7 }, // x = -0x1.77df66p-19
        .{ .u = 0x3c7bb570, .f = 0x1.021d9p+6, .df = 0x1p-19 }, // x = 0x1.f76aep-7
        .{ .u = 0x41e886d1, .f = 0x1.33136ap+98, .df = 0x1p73 }, // x = 0x1.d10da2p+4
        .{ .u = 0xc067d177, .f = 0x1.f6850cp-3, .df = 0x1p-28 }, // x = -0x1.cfa2eep+1
        .{ .u = 0xbd99da31, .f = -0x1.befe66p+3, .df = -0x1p-22 }, // x = -0x1.33b462p-4
        .{ .u = 0xbf54c45a, .f = -0x1.a6b4ecp+2, .df = 0x1p-23 }, // x = -0x1.a988b4p-1
        .{ .u = 0x41ee77fe, .f = 0x1.d3631cp+101, .df = -0x1p-76 }, // x = 0x1.dceffcp+4
        .{ .u = 0x3f843a64, .f = 0x1.f6c638p-1, .df = 0x1p-26 }, // x = 0x1.0874c8p+0
    };

    const t: u32 = @bitCast(x);
    const ax: u32 = t << 1;
    if (ax >= (0xff << 24)) { // x=NaN or +/-Inf
        @branchHint(.unlikely);
        if (ax == (0xff << 24)) { // x=+/-Inf
            if ((t >> 31) != 0) // x=-Inf
                return (x - x) / (x - x);

            return x; // x=+Inf
        }
        return x + x; // x=NaN
    }

    var z: f64 = scast(f64, x);
    if (ax < 0x6d000000) { // |x| < 0x1p-18
        @branchHint(.unlikely);
        const d: f64 = (0x1.fa658c23b1578p-1 - 0x1.d0a118f324b63p-1 * z) * z - 0x1.2788cfc6fb619p-1;
        const f: f64 = 1 / z + d;
        const r: f32 = scast(f32, f);
        const rt: u64 = @bitCast(f);
        if (((rt + 2) & 0xfffffff) < 4) {
            var i: u32 = 0;
            while (i < 10) {
                if (t == tb[i].u)
                    return tb[i].f + tb[i].df;

                i += 1;
            }
        }
        return r;
    }

    const fx: f32 = float.floor(x);
    if (x >= 0x1.18522p+5) {
        @branchHint(.unlikely);
        // Overflow case.  The original CORE-MATH code returns
        // 0x1p127f * 0x1p127f, but apparently some compilers replace this
        // by +Inf.
        return x * 0x1p127;
    }

    // compute k only after the overflow check, otherwise the case to integer
    // might overflow
    const k: i32 = scast(i32, fx);
    if (fx == x) { // x is integer
        @branchHint(.unlikely);
        if (x == 0)
            return 1 / x;

        if (x < 0)
            return (0 - 0) / (0 - 0);

        var t0: f64 = 1;
        var x0: f64 = 1;
        var i: i32 = 1;
        while (i < k) {
            t0 *= x0;

            i += 1;
            x0 += 1;
        }

        return scast(f32, t0);
    }

    if (x < -42.0) { // negative non-integer
        @branchHint(.unlikely);
        // For x < -42, x non-integer, |gamma(x)| < 2^-151.
        const sgn: [2]f32 = .{ 0x1p-127, -0x1p-127 };
        // Underflows always happens
        return 0x1p-127 * sgn[@intCast(k & 1)];
    }

    // The array c[] stores a degree-15 polynomial approximation for
    // gamma(x).
    const c: [16]f64 = .{ 0x1.c9a76be577123p+0, 0x1.8f2754ddcf90dp+0, 0x1.0d1191949419bp+0, 0x1.e1f42cf0ae4a1p-2, 0x1.82b358a3ab638p-3, 0x1.e1f2b30cd907bp-5, 0x1.240f6d4071bd8p-6, 0x1.1522c9f3cd012p-8, 0x1.1fd0051a0525bp-10, 0x1.9808a8b96c37ep-13, 0x1.b3f78e01152b5p-15, 0x1.49c85a7e1fd04p-18, 0x1.471ca49184475p-19, -0x1.368f0b7ed9e36p-23, 0x1.882222f9049efp-23, -0x1.a69ed2042842cp-25 };

    const m: f64 = z - 0x1.7p+1;
    const i: f64 = roundeven.roundeven_finite(m);
    const step: f64 = float.copysign(@as(f64, 1), i);
    const d: f64 = m - i;
    const d2: f64 = d * d;
    const d4: f64 = d2 * d2;
    const d8: f64 = d4 * d4;
    var f: f64 = (c[0] + d * c[1]) + d2 * (c[2] + d * c[3]) + d4 * ((c[4] + d * c[5]) + d2 * (c[6] + d * c[7])) + d8 * ((c[8] + d * c[9]) + d2 * (c[10] + d * c[11]) + d4 * ((c[12] + d * c[13]) + d2 * (c[14] + d * c[15])));
    const jm: i32 = scast(i32, float.abs(i));
    var w: f64 = 1;
    if (jm != 0) {
        z -= 0.5 + step * 0.5;
        w = z;
        var j: i32 = jm - 1;
        while (j != 0) {
            z -= step;
            w *= z;

            j -= 1;
        }
    }

    if (i <= -0.5)
        w = 1 / w;

    f *= w;
    const rt: u64 = @bitCast(f);
    const r: f32 = scast(f32, f);
    // Deal with exceptional cases.
    if (((rt + 2) & 0xfffffff) < 8) {
        @branchHint(.unlikely);
        var j: u32 = 0;
        while (j < 10) {
            if (t == tb[j].u)
                return tb[j].f + tb[j].df;

            j += 1;
        }
    }
    return r;
}

// Compute the product of X + X_EPS, X + X_EPS + 1, ..., X + X_EPS + N
// - 1, in the form R * (1 + *EPS) where the return value R is an
// approximation to the product and *EPS is set to indicate the
// approximate error in the return value.  X is such that all the
// values X + 1, ..., X + N - 1 are exactly representable, and X_EPS /
// X is small enough that factors quadratic in it can be
// neglected.
fn gamma_product64(x: f64, x_eps: f64, n: i32, eps: *f64) f64 {
    var ret: f64 = x;
    eps.* = x_eps / x;
    var i: i32 = 1;
    while (i < n) {
        eps.* += x_eps / (x + scast(f64, i));
        var lo: f64 = undefined;
        mul_split.mul_split64(&ret, &lo, ret, x + scast(f64, i));
        eps.* += lo / ret;

        i += 1;
    }

    return ret;
}

// Return gamma (X), for positive X less than 184, in the form R *
// 2^(*EXP2_ADJ), where R is the return value and *EXP2_ADJ is set to
// avoid overflow or underflow in intermediate calculations.
fn gamma_positive64(x: f64, exp2_adj: *i32) f64 {
    // Coefficients B_2k / 2k(2k-1) of x^-(2k-1) inside exp in Stirling's
    // approximation to gamma function.
    const gamma_coeff: [6]f64 = .{
        0x1.5555555555555p-4,
        -0xb.60b60b60b60b8p-12,
        0x3.4034034034034p-12,
        -0x2.7027027027028p-12,
        0x3.72a3c5631fe46p-12,
        -0x7.daac36664f1f4p-12,
    };

    var local_signgam: i32 = undefined;
    if (x < 0.5) {
        exp2_adj.* = 0;
        return float.exp(lgamma.lgamma_r64(x + 1, &local_signgam)) / x;
    } else if (x <= 1.5) {
        exp2_adj.* = 0;
        return float.exp(lgamma.lgamma_r64(x, &local_signgam));
    } else if (x < 6.5) {
        // Adjust into the range for using exp (lgamma).
        exp2_adj.* = 0;
        const n: f64 = float.ceil(x - 1.5);
        const x_adj: f64 = x - n;
        var eps: f64 = undefined;
        const prod: f64 = gamma_product64(x_adj, 0, scast(i32, n), &eps);
        return float.exp(lgamma.lgamma_r64(x_adj, &local_signgam)) * prod * (1 + eps);
    } else {
        var eps: f64 = 0;
        var x_eps: f64 = 0;
        var x_adj: f64 = x;
        var prod: f64 = 1;
        if (x < 12) {
            // Adjust into the range for applying Stirling's
            // approximation.
            const n: f64 = float.ceil(12 - x);
            x_adj = x + n;
            x_eps = (x - (x_adj - n));
            prod = gamma_product64(x_adj - n, x_eps, scast(i32, n), &eps);
        }
        // The result is now gamma (X_ADJ + X_EPS) / (PROD * (1 + EPS)).
        // Compute gamma (X_ADJ + X_EPS) using Stirling's approximation,
        // starting by computing pow (X_ADJ, X_ADJ) with a power of 2
        // factored out.
        const x_adj_int: f64 = float.round(x_adj);
        const x_adj_frac: f64 = x_adj - x_adj_int;
        var x_adj_log2: i32 = undefined;
        var x_adj_mant: f64 = float.frexp(x_adj, &x_adj_log2);
        if (x_adj_mant < std.math.sqrt1_2) {
            x_adj_log2 -= 1;
            x_adj_mant *= 2;
        }

        exp2_adj.* = x_adj_log2 * scast(i32, x_adj_int);
        var h1: f64 = undefined;
        var l1: f64 = undefined;
        mul_split.mul_split64(&h1, &l1, float.pow(x_adj_mant, x_adj), float.exp2(scast(f64, x_adj_log2) * x_adj_frac));
        var h2: f64 = undefined;
        var l2: f64 = undefined;
        mul_split.mul_split64(&h2, &l2, float.exp(-x_adj), float.sqrt(2 * std.math.pi / x_adj));
        mul_split.mul_expansion64(&h1, &l1, h1, l1, h2, l2);
        // Divide by prod * (1 + eps).
        mul_split.div_expansion64(&h1, &l1, h1, l1, prod, prod * eps);
        var exp_adj: f64 = x_eps * float.log(x_adj);
        var bsum: f64 = gamma_coeff[5];
        const x_adj2: f64 = x_adj * x_adj;
        var i: u32 = 1;
        while (i <= 5) {
            bsum = bsum / x_adj2 + gamma_coeff[5 - i];

            i += 1;
        }
        exp_adj += bsum / x_adj;
        // Now return (h1+l1) * exp(exp_adj), where exp_adj is small.
        l1 += h1 * float.expm1(exp_adj);
        return h1 + l1;
    }
}

fn gamma_r64(x: f64, signgamp: *i32) f64 {
    var hx: i32 = undefined;
    var lx: i32 = undefined;
    dbl64.extractWords(&hx, &lx, x);

    if (((hx & 0x7fffffff) | lx) == 0) {
        @branchHint(.unlikely);
        // Return value for x == 0 is Inf with divide by zero exception.
        signgamp.* = 0;
        return 1 / x;
    }

    if (hx < 0 and @as(u32, @bitCast(hx)) < 0xfff00000 and float.rint(x) == x) {
        // Return value for integer x < 0 is NaN with invalid exception.
        signgamp.* = 0;
        return (x - x) / (x - x);
    }

    if (@as(u32, @bitCast(hx)) == 0xfff00000 and lx == 0) {
        @branchHint(.unlikely);
        // x == -Inf.  According to ISO this is NaN.
        signgamp.* = 0;
        return x - x;
    }

    if ((hx & 0x7ff00000) == 0x7ff00000) {
        @branchHint(.unlikely);
        // Positive infinity (return positive infinity) or NaN (return NaN).
        signgamp.* = 0;
        return x + x;
    }

    var ret: f64 = undefined;
    if (x >= 172) {
        // Overflow.
        signgamp.* = 0;
        return std.math.floatMax(f64) * std.math.floatMax(f64);
    } else {
        if (x > 0) {
            signgamp.* = 0;
            var exp2_adj: i32 = undefined;
            const tret: f64 = gamma_positive64(x, &exp2_adj);
            ret = float.scalbn(tret, exp2_adj);
        } else if (x >= -std.math.floatEps(f64) / 4.0) {
            signgamp.* = 0;
            ret = 1 / x;
        } else {
            const tx: f64 = float.trunc(x);
            signgamp.* = if (tx == 2 * float.trunc(tx / 2)) -1 else 1;
            if (x <= -184.0) {
                // Underflow.
                ret = std.math.floatMin(f64) * std.math.floatMin(f64);
            } else {
                var frac: f64 = tx - x;
                if (frac > 0.5)
                    frac = 1 - frac;

                const sinpix: f64 = if (frac <= 0.25) float.sin(std.math.pi * frac) else float.cos(std.math.pi * (0.5 - frac));
                var exp2_adj: i32 = undefined;
                var h2: f64 = gamma_positive64(-x, &exp2_adj);
                var h1: f64 = undefined;
                var l1: f64 = undefined;
                mul_split.mul_split64(&h1, &l1, sinpix, h2);
                // sinpix*gamma_positive(.) = h1 + l1
                var l2: f64 = undefined;
                mul_split.mul_split64(&h2, &l2, h1, x);
                // h1*x = h2 + l2
                // (h1 + l1) * x = h1*x + l1*x = h2 + l2 + l1*x
                l2 += l1 * x;
                // x*sinpix*gamma_positive(.) ~ h2 + l2
                h1 = 0x3.243f6a8885a3p+0; // binary64 approximation of Pi
                l1 = 0x8.d313198a2e038p-56; // |h1+l1-Pi| < 3e-33
                // Now we divide h1 + l1 by h2 + l2.
                mul_split.div_expansion64(&h1, &l1, h1, l1, h2, l2);
                ret = float.scalbn(-h1, -exp2_adj);

                if (ret < std.math.floatMin(f64)) {
                    const vret: f64 = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }
            }
        }
    }

    if (std.math.isInf(ret) and x != 0) {
        if (signgamp.* < 0) {
            ret = -float.copysign(std.math.floatMax(f64), ret) * std.math.floatMax(f64);
            ret = -ret;
        } else {
            ret = float.copysign(std.math.floatMax(f64), ret) * std.math.floatMax(f64);
        }
        return ret;
    } else if (ret == 0) {
        if (signgamp.* < 0) {
            ret = -float.copysign(std.math.floatMin(f64), ret) * std.math.floatMin(f64);
            ret = -ret;
        } else {
            ret = float.copysign(std.math.floatMin(f64), ret) * std.math.floatMin(f64);
        }

        return ret;
    } else {
        return ret;
    }
}

// Compute the product of X + X_EPS, X + X_EPS + 1, ..., X + X_EPS + N
// - 1, in the form R * (1 + *EPS) where the return value R is an
// approximation to the product and *EPS is set to indicate the
// approximate error in the return value.  X is such that all the
// values X + 1, ..., X + N - 1 are exactly representable, and X_EPS /
// X is small enough that factors quadratic in it can be
// neglected.
fn gamma_product128(x: f128, x_eps: f128, n: i32, eps: *f128) f128 {
    var ret: f128 = x;
    eps.* = x_eps / x;
    var i: i32 = 1;
    while (i < n) {
        eps.* += x_eps / (x + scast(f128, i));
        var lo: f128 = undefined;
        mul_split.mul_split128(&ret, &lo, ret, x + scast(f128, i));
        eps.* += lo / ret;

        i += 1;
    }

    return ret;
}

// Return gamma (X), for positive X less than 1775, in the form R *
// 2^(*EXP2_ADJ), where R is the return value and *EXP2_ADJ is set to
// avoid overflow or underflow in intermediate calculations.
fn gamma_positive128(x: f128, exp2_adj: *i32) f128 {
    // Coefficients B_2k / 2k(2k-1) of x^-(2k-1) inside exp in Stirling's
    // approximation to gamma function.
    const gamma_coeff: [14]f128 = .{
        0x1.5555555555555555555555555555p-4,
        -0xb.60b60b60b60b60b60b60b60b60b8p-12,
        0x3.4034034034034034034034034034p-12,
        -0x2.7027027027027027027027027028p-12,
        0x3.72a3c5631fe46ae1d4e700dca8f2p-12,
        -0x7.daac36664f1f207daac36664f1f4p-12,
        0x1.a41a41a41a41a41a41a41a41a41ap-8,
        -0x7.90a1b2c3d4e5f708192a3b4c5d7p-8,
        0x2.dfd2c703c0cfff430edfd2c703cp-4,
        -0x1.6476701181f39edbdb9ce625987dp+0,
        0xd.672219167002d3a7a9c886459cp+0,
        -0x9.cd9292e6660d55b3f712eb9e07c8p+4,
        0x8.911a740da740da740da740da741p+8,
        -0x8.d0cc570e255bf59ff6eec24b49p+12,
    };

    var local_signgam: i32 = undefined;
    if (x < 0.5) {
        exp2_adj.* = 0;
        return float.exp(lgamma.lgamma_r128(x + 1, &local_signgam)) / x;
    } else if (x <= 1.5) {
        exp2_adj.* = 0;
        return float.exp(lgamma.lgamma_r128(x, &local_signgam));
    } else if (x < 12.5) {
        // Adjust into the range for using exp (lgamma).
        exp2_adj.* = 0;
        const n: f128 = float.ceil(x - 1.5);
        const x_adj: f128 = x - n;
        var eps: f128 = undefined;
        const prod: f128 = gamma_product128(x_adj, 0, scast(i32, n), &eps);
        return float.exp(lgamma.lgamma_r128(x_adj, &local_signgam)) * prod * (1 + eps);
    } else {
        var eps: f128 = 0;
        var x_eps: f128 = 0;
        var x_adj: f128 = x;
        var prod: f128 = 1;
        if (x < 24) {
            // Adjust into the range for applying Stirling's
            // approximation.
            const n: f128 = float.ceil(24 - x);
            x_adj = x + n;
            x_eps = (x - (x_adj - n));
            prod = gamma_product128(x_adj - n, x_eps, scast(i32, n), &eps);
        }
        // The result is now gamma (X_ADJ + X_EPS) / (PROD * (1 + EPS)).
        // Compute gamma (X_ADJ + X_EPS) using Stirling's approximation,
        // starting by computing pow (X_ADJ, X_ADJ) with a power of 2
        // factored out.
        var exp_adj: f128 = -eps;
        const x_adj_int: f128 = float.round(x_adj);
        const x_adj_frac: f128 = x_adj - x_adj_int;
        var x_adj_log2: i32 = undefined;
        var x_adj_mant: f128 = float.frexp(x_adj, &x_adj_log2);
        if (x_adj_mant < std.math.sqrt1_2) {
            x_adj_log2 -= 1;
            x_adj_mant *= 2;
        }

        exp2_adj.* = x_adj_log2 * scast(i32, x_adj_int);
        const ret: f128 = float.pow(x_adj_mant, x_adj) * float.exp2(scast(f128, x_adj_log2) * x_adj_frac) * float.exp(-x_adj) * float.sqrt(2 * std.math.pi / x_adj) / prod;
        exp_adj += x_eps * float.log(x_adj);
        var bsum: f128 = gamma_coeff[13];
        const x_adj2: f128 = x_adj * x_adj;
        var i: u32 = 1;
        while (i <= 13) {
            bsum = bsum / x_adj2 + gamma_coeff[13 - i];

            i += 1;
        }

        exp_adj += bsum / x_adj;
        return ret + ret * float.expm1(exp_adj);
    }
}

fn gamma_r128(x: f128, signgamp: *i32) f128 {
    var hx: i64 = undefined;
    var lx: i64 = undefined;
    ldbl128.getWords(&hx, &lx, x);

    if (((hx & 0x7fffffffffffffff) | lx) == 0) {
        // Return value for x == 0 is Inf with divide by zero exception.
        signgamp.* = 0;
        return 1 / x;
    }

    if (hx < 0 and @as(u64, @bitCast(hx)) < 0xffff000000000000 and float.rint(x) == x) {
        // Return value for integer x < 0 is NaN with invalid exception.
        signgamp.* = 0;
        return (x - x) / (x - x);
    }

    if (hx == 0xffff000000000000 and lx == 0) {
        // x == -Inf.  According to ISO this is NaN.
        signgamp.* = 0;
        return x - x;
    }

    if ((hx & 0x7fff000000000000) == 0x7fff000000000000) {
        // Positive infinity (return positive infinity) or NaN (return NaN).
        signgamp.* = 0;
        return x + x;
    }

    var ret: f128 = undefined;
    if (x >= 1756) {
        // Overflow.
        signgamp.* = 0;
        return std.math.floatMax(f128) * std.math.floatMax(f128);
    } else {
        if (x > 0) {
            signgamp.* = 0;
            var exp2_adj: i32 = undefined;
            ret = gamma_positive128(x, &exp2_adj);
            ret = float.scalbn(ret, exp2_adj);
        } else if (x >= -std.math.floatEps(f128) / 4) {
            signgamp.* = 0;
            ret = 1 / x;
        } else {
            const tx: f128 = float.trunc(x);
            signgamp.* = if (tx == 2 * float.trunc(tx / 2)) -1 else 1;
            if (x <= -1775) {
                // Underflow.
                ret = std.math.floatMin(f128) * std.math.floatMin(f128);
            } else {
                var frac: f128 = tx - x;
                if (frac > 0.5)
                    frac = 1 - frac;

                const sinpix: f128 = if (frac <= 0.25) float.sin(std.math.pi * frac) else float.cos(std.math.pi * (0.5 - frac));
                var exp2_adj: i32 = undefined;
                ret = std.math.pi / (-x * sinpix * gamma_positive128(-x, &exp2_adj));
                ret = float.scalbn(ret, -exp2_adj);

                if (ret < std.math.floatMin(f128)) {
                    const vret: f128 = ret * ret;
                    std.mem.doNotOptimizeAway(vret);
                }
            }
        }
    }

    if (std.math.isInf(ret) and x != 0) {
        if (signgamp.* < 0) {
            return -(-float.copysign(std.math.floatMax(f128), ret) * std.math.floatMax(f128));
        } else {
            return float.copysign(std.math.floatMax(f128), ret) * std.math.floatMax(f128);
        }
    } else if (ret == 0) {
        if (signgamp.* < 0) {
            return -(-float.copysign(std.math.floatMin(f128), ret) * std.math.floatMin(f128));
        } else {
            return float.copysign(std.math.floatMin(f128), ret) * std.math.floatMin(f128);
        }
    } else {
        return ret;
    }
}
