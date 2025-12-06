const types = @import("../types.zig");
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_rem_pio2f.c
//
// Original copyright notice:
// e_rem_pio2f.c -- float version of e_rem_pio2.c
// Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
// Debugged and optimized by Bruce D. Evans.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
pub fn rem_pio2_32(x: f32, y: *f64) i32 {
    const hx: i32 = @bitCast(x);
    const ix: i32 = hx & 0x7fffffff;

    // 33 + 53 bit pi is good enough for medium size
    if (ix < 0x4dc90fdb) { // |x| ~< 2**28 * (pi/2), medium size
        var f: f64 = types.scast(f64, x) * 6.36619772367581382433e-1 + 0x1.8p52;
        f -= 0x1.8p52;
        const n: i32 = types.scast(i32, f);
        const r: f64 = types.scast(f64, x) - f * 1.57079631090164184570e+0;
        const w: f64 = f * 1.58932547735281966916e-8;
        y.* = r - w;
        return n;
    }

    // All other (large) arguments
    if (ix >= 0x7f800000) { // x is inf or NaN
        y.* = x - x;
        return 0;
    }

    // Set z = scalbn(|x|, ilogb(|x|) - 23)
    const e0: i32 = (ix >> 23) -% 150; // e0 = ilogb(|x|) - 23
    const z: f32 = @bitCast(ix -% (e0 << 23));
    var tx: f64 = types.scast(f64, z);
    var ty: f64 = undefined;
    const n: i32 = k_rem_pio2_64(@ptrCast(&tx), @ptrCast(&ty), e0, 1, 0);

    if (hx < 0) {
        y.* = -ty;
        return -n;
    }

    y.* = ty;
    return n;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/e_rem_pio2.c
//
// Original copyright notice:
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
pub fn rem_pio2_64(x: f64, y: *[2]f64) i32 {
    const hx: i32 = @bitCast(dbl64.getHighPart(x));
    const ix: i32 = hx & 0x7fffffff;

    var goto_medium: bool = false;
    if (ix <= 0x400f6a7a) { // |x| ~<= 5 * pi/4
        if ((ix & 0xfffff) == 0x921fb) // |x| ~= pi/2 or 2 * pi/2
            goto_medium = true;

        if (!goto_medium) {
            if (ix <= 0x4002d97c) { // |x| ~<= 3 * pi/4
                if (hx > 0) {
                    const z: f64 = x - 1.57079632673412561417e+0; // One round good to 85 bits
                    y[0] = z - 6.07710050650619224932e-11;
                    y[1] = (z - y[0]) - 6.07710050650619224932e-11;
                    return 1;
                } else {
                    const z: f64 = x + 1.57079632673412561417e+0;
                    y[0] = z + 6.07710050650619224932e-11;
                    y[1] = (z - y[0]) + 6.07710050650619224932e-11;
                    return -1;
                }
            } else {
                if (hx > 0) {
                    const z: f64 = x - 2.0 * 1.57079632673412561417e+0;
                    y[0] = z - 2.0 * 6.07710050650619224932e-11;
                    y[1] = (z - y[0]) - 2.0 * 6.07710050650619224932e-11;
                    return 2;
                } else {
                    const z: f64 = x + 2.0 * 1.57079632673412561417e+0;
                    y[0] = z + 2.0 * 6.07710050650619224932e-11;
                    y[1] = (z - y[0]) + 2.0 * 6.07710050650619224932e-11;
                    return -2;
                }
            }
        }
    }

    if (!goto_medium) {
        if (ix <= 0x401c463b) { // |x| ~<= 9 * pi/4
            if (ix <= 0x4015fdbc) { // |x| ~<= 7 * pi/4
                if (ix == 0x4012d97c) // |x| ~= 3 * pi/2
                    goto_medium = true;

                if (!goto_medium) {
                    if (hx > 0) {
                        const z: f64 = x - 3.0 * 1.57079632673412561417e+0;
                        y[0] = z - 3.0 * 6.07710050650619224932e-11;
                        y[1] = (z - y[0]) - 3.0 * 6.07710050650619224932e-11;
                        return 3;
                    } else {
                        const z: f64 = x + 3.0 * 1.57079632673412561417e+0;
                        y[0] = z + 3.0 * 6.07710050650619224932e-11;
                        y[1] = (z - y[0]) + 3.0 * 6.07710050650619224932e-11;
                        return -3;
                    }
                } else {
                    if (ix == 0x401921fb) // |x| ~= 4 * pi/2
                        goto_medium = true;

                    if (!goto_medium) {
                        if (hx > 0) {
                            const z: f64 = x - 4.0 * 1.57079632673412561417e+0;
                            y[0] = z - 4.0 * 6.07710050650619224932e-11;
                            y[1] = (z - y[0]) - 4.0 * 6.07710050650619224932e-11;
                            return 4;
                        } else {
                            const z: f64 = x + 4.0 * 1.57079632673412561417e+0;
                            y[0] = z + 4 * 6.07710050650619224932e-11;
                            y[1] = (z - y[0]) + 4.0 * 6.07710050650619224932e-11;
                            return -4;
                        }
                    }
                }
            }
        }
    }

    if (ix < 0x413921fb or goto_medium) { // |x| ~< 2**20 * (pi/2), medium size
        var f: f64 = x * 6.36619772367581382433e-1 + 0x1.8p52;
        f -= 0x1.8p52;
        const n: i32 = types.scast(i32, f);
        var r: f64 = x - f * 1.57079632673412561417e+0;
        var w: f64 = f * 6.07710050650619224932e-11; // 1st round good to 85 bit
        {
            const j: i32 = ix >> 20;
            y[0] = r - w;
            var high: i32 = @bitCast(dbl64.getHighPart(y[0]));
            var i: i32 = j -% ((high >> 20) & 0x7ff);
            if (i > 16) { // 2nd iteration needed, good to 118
                var t: f64 = r;
                w = f * 6.07710050630396597660e-11;
                r = t - w;
                w = f * 2.02226624879595063154e-21 - ((t - r) - w);
                y[0] = r - w;
                high = @bitCast(dbl64.getHighPart(y[0]));
                i = j -% ((high >> 20) & 0x7ff);
                if (i > 49) { // 3rd iteration need, 151 bits accuracy
                    t = r; // Will cover all possible cases
                    w = f * 2.02226624871116645580e-21;
                    r = t - w;
                    w = f * 8.47842766036889956997e-32 - ((t - r) - w);
                    y[0] = r - w;
                }
            }
        }

        y[1] = (r - y[0]) - w;
        return n;
    }

    // All other (large) arguments
    if (ix >= 0x7ff00000) { // x is inf or NaN
        y[0] = x - x;
        y[1] = x - x;
        return 0;
    }

    // Set z = scalbn(|x|, ilogb(x) - 23)
    const low: u32 = dbl64.getLowPart(x);
    const e0: i32 = (ix >> 20) -% 1046; // e0 = ilogb(z) - 23
    var z: f64 = dbl64.Parts.toFloat(.{ .msw = @bitCast(ix -% (e0 << 20)), .lsw = low });

    var tx: [3]f64 = undefined;
    var i: u32 = 0;
    while (i < 2) : (i += 1) {
        tx[i] = types.scast(f64, types.scast(i32, z));
        z = (z - tx[i]) * 1.67772160000000000000e+7;
    }

    tx[2] = z;
    var nx: i32 = 3;
    while (nx > 0 and tx[types.scast(u32, nx - 1)] == 0.0) : (nx -= 1) {}

    var ty: [2]f64 = undefined;
    const n: i32 = k_rem_pio2_64(&tx, &ty, e0, nx, 1);

    if (hx < 0) {
        y[0] = -ty[0];
        y[1] = -ty[1];
        return -n;
    }

    y[0] = ty[0];
    y[1] = ty[1];
    return n;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_rem_pio2l.h
//
// Original copyright notice:
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
// Copyright (c) 2008 Steven G. Kargl, David Schultz, Bruce D. Evans.
//
// Developed at SunSoft, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
pub fn rem_pio2_128(x: f128, y: *[2]f128) i64 {
    const u: ldbl128.ShapeSplit = .fromFloat(x);
    const sign: u1 = u.sign;
    const ex: i32 = @bitCast(@as(u32, @intCast(u.exponent)));
    if (ex < (16384 - 1) + 45 or ex == (16384 - 1) + 45 and
        u.mantissa_high < 0x921fb54442d1)
    { // |x| ~< 2**45 * (pi/2), medium size
        var f: f128 = x * 6.3661977236758134307553505349005747e-1 + 0x1.8p112;
        f = f - 0x1.8p112;
        const n: i64 = types.scast(i64, f);
        var r: f128 = x - f * 1.5707963267948966192292994253909555e+0;
        var w: f128 = f * 2.0222662487959507323996846200947577e-21; // 1st round good to 180 bit
        {
            const j: i32 = ex;
            y[0] = r - w;
            var uu: ldbl128.ShapeSplit = .fromFloat(y[0]);
            var ex1: i32 = @bitCast(@as(u32, @intCast(uu.exponent)));
            var i: i32 = j -% ex1;
            if (i > 51) { // 2nd iteration needed, good to 248
                var t: f128 = r;
                w = f * 2.0222662487959507323994779168837751e-21;
                r = t - w;
                w = f * 2.0670321098263988236496903051604844e-43 - ((t - r) - w);
                y[0] = r - w;
                uu = .fromFloat(y[0]);
                ex1 = @bitCast(@as(u32, @intCast(uu.exponent)));
                i = j - ex1;
                if (i > 119) { // 3rd iteration need, 316 bits acc
                    t = r; // Will cover all possible cases
                    w = f * 2.0670321098263988236499468110329591e-43;
                    r = t - w;
                    w = f * -2.5650587247459238361625433492959285e-65 - ((t - r) - w);
                    y[0] = r - w;
                }
            }
        }

        y[1] = (r - y[0]) - w;
        return n;
    }

    // All other (large) arguments
    if (ex == 0x7fff) { // x is inf or NaN
        y[0] = x - x;
        y[1] = x - x;
        return 0;
    }

    // Set z = scalbn(|x|, ilogb(x) - 23)
    var uu: ldbl128.ShapeSplit = .fromFloat(x);
    const e0: i32 = ex -% (16384 - 1) -% 23; // e0 = ilogb(|x|) - 23
    const expsign: u32 = @bitCast(ex -% e0);
    uu.exponent = @truncate(expsign & 0x7fff);
    uu.sign = @bitCast(@as(u1, @truncate(@as(u32, @bitCast(expsign >> 15)))));
    var z: f128 = uu.toFloat();
    var i: u32 = 0;
    var tx: [5]f64 = undefined;
    while (i < 4) : (i += 1) {
        tx[i] = types.scast(f64, types.scast(i32, z));
        z = (z - types.scast(f128, tx[i])) * 1.67772160000000000000e+7;
    }
    tx[4] = types.scast(f64, z);
    var nx: i32 = 5;
    while (nx > 0 and tx[types.scast(u32, nx - 1)] == 0.0) : (nx -= 1) {} // Skip zero term
    var ty: [3]f64 = undefined;
    const n: i32 = k_rem_pio2_64(@ptrCast(&tx), @ptrCast(&ty), e0, nx, 3);
    const t: f128 = types.scast(f128, ty[2]) + types.scast(f128, ty[1]);
    const r: f128 = t + types.scast(f128, ty[0]);
    const w: f128 = types.scast(f128, ty[0]) - (r - t);
    if (sign != 0) {
        y[0] = -r;
        y[1] = -w;
        return -n;
    }

    y[0] = r;
    y[1] = w;
    return types.scast(i64, n);
}

fn goto_recompute64(x: [*]f64, j: *i32, i: *i32, f: *[20]f64, q: *[20]f64, jz: *i32, iq: *[20]i32, z: *f64, n: *i32, ih: *i32, jk: i32, jx: i32, jv: i32, q0: i32) void {
    // Distill q[] into iq[] reversingly
    i.* = 0;
    j.* = jz.*;
    z.* = q[types.scast(u32, jz.*)];
    while (j.* > 0) {
        const fw: f64 = types.scast(f64, types.scast(i32, 5.96046447753906250000e-8 * z.*));
        iq[types.scast(u32, i.*)] = types.scast(i32, z.* - 1.67772160000000000000e+7 * fw);
        z.* = q[types.scast(u32, j.* - 1)] + fw;

        i.* += 1;
        j.* -= 1;
    }

    // Compute n
    z.* = float.scalbn(z.*, q0); // Actual value of z
    z.* -= 8.0 * float.floor(z.* * 0.125); // Trim off integer >= 8
    n.* = types.scast(i32, z.*);
    z.* -= types.scast(f64, n.*);
    ih.* = 0;
    if (q0 > 0) { // Need iq[jz - 1] to determine n
        i.* = (iq[types.scast(u32, jz.* -% 1)] >> @as(u5, @intCast(24 -% q0)));
        n.* +%= i.*;
        iq[types.scast(u32, jz.* - 1)] -%= i.* << @as(u5, @intCast(24 -% q0));
        ih.* = iq[types.scast(u32, jz.* - 1)] >> @as(u5, @intCast(23 -% q0));
    } else if (q0 == 0) {
        ih.* = iq[types.scast(u32, jz.* - 1)] >> 23;
    } else if (z.* >= 0.5) {
        ih.* = 2;
    }

    if (ih.* > 0) { // q > 0.5
        n.* +%= 1;
        var carry: i32 = 0;
        i.* = 0;
        while (i.* < jz.*) : (i.* += 1) { // compute 1-q
            j.* = iq[types.scast(u32, i.*)];
            if (carry == 0) {
                if (j.* != 0) {
                    carry = 1;
                    iq[types.scast(u32, i.*)] = 0x1000000 -% j.*;
                }
            } else {
                iq[types.scast(u32, i.*)] = 0xffffff -% j.*;
            }
        }

        if (q0 > 0) { // Rare case: chance is 1 in 12
            switch (q0) {
                1 => iq[types.scast(u32, jz.* - 1)] &= 0x7fffff,
                2 => iq[types.scast(u32, jz.* - 1)] &= 0x3fffff,
                else => {},
            }
        }

        if (ih.* == 2) {
            z.* = 1.0 - z.*;
            if (carry != 0)
                z.* -= float.scalbn(@as(f64, 1.0), q0);
        }
    }

    // Check if recomputation is needed
    if (z.* == 0) {
        j.* = 0;
        i.* = jz.* -% 1;
        while (i.* >= jk) : (i.* -= 1) {
            j.* |= iq[types.scast(u32, i.*)];
        }

        if (j.* == 0) { // Need recomputation
            var k: i32 = 1;
            while (iq[types.scast(u32, jk -% k)] == 0) : (k += 1) {} // k = no. of terms needed

            i.* = jz.* +% 1;
            while (i.* <= jz.* +% k) : (i.* += 1) { // Add q[jz+1] to q[jz+k]
                f[types.scast(u32, jx +% i.*)] = types.scast(f64, ipio2_64[types.scast(u32, jv +% i.*)]);

                j.* = 0;
                var fw: f64 = 0;
                while (j.* <= jx) : (j.* += 1) {
                    fw += x[types.scast(u32, j.*)] * f[types.scast(u32, jx +% i.* -% j.*)];
                }

                q[types.scast(u32, i.*)] = fw;
            }

            jz.* += k;

            goto_recompute64(x, j, i, f, q, jz, iq, z, n, ih, jk, jx, jv, q0);
        }
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/k_rem_pio2.c
//
// Original copyright notice:
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
pub fn k_rem_pio2_64(x: [*]f64, y: [*]f64, e0: i32, nx: i32, prec: i32) i32 {

    // Initialize jk
    const jk: i32 = switch (prec) {
        0 => 3,
        1 => 4,
        2 => 4,
        3 => 6,
        else => unreachable,
    };
    const jp: i32 = jk;

    // Determine jx, jv, q0, note that 3 > q0
    const jx: i32 = nx - 1;
    var jv: i32 = @divTrunc(e0 -% 3, 24);
    if (jv < 0)
        jv = 0;
    var q0: i32 = e0 -% 24 * (jv +% 1);

    // Set up f[0] to f[jx + jk] where f[jx + jk] = ipio2_64[jv + jk]
    var j: i32 = jv -% jx;
    const m: i32 = jx +% jk;
    var i: i32 = 0;
    var f: [20]f64 = undefined;
    while (i <= m) {
        f[types.scast(u32, i)] = if (j < 0)
            0.0
        else
            types.scast(f64, ipio2_64[types.scast(u32, j)]);

        i += 1;
        j += 1;
    }

    // Compute q[0], q[1], ..., q[jk]
    var q: [20]f64 = undefined;
    i = 0;
    while (i <= jk) {
        j = 0;
        var fw: f64 = 0;
        while (j <= jx) {
            fw += x[types.scast(u32, j)] * f[types.scast(u32, jx + i - j)];

            j += 1;
        }

        q[@intCast(i)] = fw;

        i += 1;
    }

    var jz: i32 = jk;

    var iq: [20]i32 = undefined;
    var z: f64 = undefined;
    var n: i32 = undefined;
    var ih: i32 = undefined;
    goto_recompute64(x, &j, &i, &f, &q, &jz, &iq, &z, &n, &ih, jk, jx, jv, q0);

    // Chop off zero terms
    if (z == 0) {
        jz -%= 1;
        q0 -%= 24;

        while (iq[types.scast(u32, jz)] == 0) {
            jz -%= 1;
            q0 -%= 24;
        }
    } else { // Break z into 24-bit if necessary
        z = float.scalbn(z, -q0);
        if (z >= 1.67772160000000000000e+7) {
            const fw: f64 = types.scast(f64, types.scast(i32, 5.96046447753906250000e-8 * z));
            iq[types.scast(u32, jz)] = types.scast(i32, z - 1.67772160000000000000e+7 * fw);
            jz +%= 1;
            q0 +%= 24;
            iq[types.scast(u32, jz)] = types.scast(i32, fw);
        } else {
            iq[types.scast(u32, jz)] = types.scast(i32, z);
        }
    }

    // Convert integer "bit" chunk to floating-point value
    var fw: f64 = float.scalbn(@as(f64, 1.0), q0);
    i = jz;
    while (i >= 0) : (i -= 1) {
        q[types.scast(u32, i)] = fw * types.scast(f64, iq[types.scast(u32, i)]);
        fw *= 5.96046447753906250000e-8;
    }

    // Compute pio2_64[0, ..., jp] * q[jz, ..., 0]
    var fq: [20]f64 = undefined;
    i = jz;
    while (i >= 0) {
        fw = 0.0;
        var k: i32 = 0;
        while (k <= jp and k <= jz -% i) : (k += 1) {
            fw += pio2_64[types.scast(u32, k)] * q[types.scast(u32, i +% k)];
        }

        fq[types.scast(u32, jz -% i)] = fw;

        i -= 1;
    }

    // Compress fq[] into y[]
    switch (prec) {
        0 => {
            fw = 0.0;
            i = jz;
            while (i >= 0) : (i -= 1) {
                fw += fq[types.scast(u32, i)];
            }

            y[0] = if (ih == 0) fw else -fw;
        },
        1, 2 => {
            var fv: f64 = 0.0;
            i = jz;
            while (i >= 0) : (i -= 1) {
                fv = fv + fq[types.scast(u32, i)];
            }

            y[0] = if (ih == 0) fv else -fv;

            fv = fq[0] - fv;
            i = 1;
            while (i <= jz) {
                fv = fv + fq[types.scast(u32, i)];

                i += 1;
            }

            y[1] = if (ih == 0) fv else -fv;
        },
        3 => { // Painful
            i = jz;
            while (i > 0) {
                const fv: f64 = fq[types.scast(u32, i -% 1)] + fq[types.scast(u32, i)];
                fq[types.scast(u32, i)] += fq[types.scast(u32, i -% 1)] - fv;
                fq[types.scast(u32, i -% 1)] = fv;

                i -= 1;
            }

            i = jz;
            while (i > 1) : (i -= 1) {
                const fv: f64 = fq[types.scast(u32, i -% 1)] + fq[types.scast(u32, i)];
                fq[types.scast(u32, i)] += fq[types.scast(u32, i -% 1)] - fv;
                fq[types.scast(u32, i -% 1)] = fv;
            }

            fw = 0.0;
            i = jz;
            while (i >= 2) : (i -= 1) {
                fw += fq[types.scast(u32, i)];
            }

            if (ih == 0) {
                y[0] = fq[0];
                y[1] = fq[1];
                y[2] = fw;
            } else {
                y[0] = -fq[0];
                y[1] = -fq[1];
                y[2] = -fw;
            }
        },
        else => {},
    }

    return n & 7;
}

const ipio2_64: [690]i32 = .{
    0xa2f983, 0x6e4e44, 0x1529fc, 0x2757d1, 0xf534dd, 0xc0db62,
    0x95993c, 0x439041, 0xfe5163, 0xabdebb, 0xc561b7, 0x246e3a,
    0x424dd2, 0xe00649, 0x2eea09, 0xd1921c, 0xfe1deb, 0x1cb129,
    0xa73ee8, 0x8235f5, 0x2ebb44, 0x84e99c, 0x7026b4, 0x5f7e41,
    0x3991d6, 0x398353, 0x39f49c, 0x845f8b, 0xbdf928, 0x3b1ff8,
    0x97ffde, 0x05980f, 0xef2f11, 0x8b5a0a, 0x6d1f6d, 0x367ecf,
    0x27cb09, 0xb74f46, 0x3f669e, 0x5fea2d, 0x7527ba, 0xc7ebe5,
    0xf17b3d, 0x0739f7, 0x8a5292, 0xea6bfb, 0x5fb11f, 0x8d5d08,
    0x560330, 0x46fc7b, 0x6babf0, 0xcfbc20, 0x9af436, 0x1da9e3,
    0x91615e, 0xe61b08, 0x659985, 0x5f14a0, 0x68408d, 0xffd880,
    0x4d7327, 0x310606, 0x1556ca, 0x73a8c9, 0x60e27b, 0xc08c6b,
    0x47c419, 0xc367cd, 0xdce809, 0x2a8359, 0xc4768b, 0x961ca6,
    0xddaf44, 0xd15719, 0x053ea5, 0xff0705, 0x3f7e33, 0xe832c2,
    0xde4f98, 0x327dbb, 0xc33d26, 0xef6b1e, 0x5ef89f, 0x3a1f35,
    0xcaf27f, 0x1d87f1, 0x21907c, 0x7c246a, 0xfa6ed5, 0x772d30,
    0x433b15, 0xc614b5, 0x9d19c3, 0xc2c4ad, 0x414d2c, 0x5d000c,
    0x467d86, 0x2d71e3, 0x9ac69b, 0x006233, 0x7cd2b4, 0x97a7b4,
    0xd55537, 0xf63ed7, 0x1810a3, 0xfc764d, 0x2a9d64, 0xabd770,
    0xf87c63, 0x57b07a, 0xe71517, 0x5649c0, 0xd9d63b, 0x3884a7,
    0xcb2324, 0x778ad6, 0x23545a, 0xb91f00, 0x1b0af1, 0xdfce19,
    0xff319f, 0x6a1e66, 0x615799, 0x47fbac, 0xd87f7e, 0xb76522,
    0x89e832, 0x60bfe6, 0xcdc4ef, 0x09366c, 0xd43f5d, 0xd7de16,
    0xde3b58, 0x929bde, 0x2822d2, 0xe88628, 0x4d58e2, 0x32cac6,
    0x16e308, 0xcb7de0, 0x50c017, 0xa71df3, 0x5be018, 0x34132e,
    0x621283, 0x014883, 0x5b8ef5, 0x7fb0ad, 0xf2e91e, 0x434a48,
    0xd36710, 0xd8ddaa, 0x425fae, 0xce616a, 0xa4280a, 0xb499d3,
    0xf2a606, 0x7f775c, 0x83c2a3, 0x883c61, 0x78738a, 0x5a8caf,
    0xbdd76f, 0x63a62d, 0xcbbff4, 0xef818d, 0x67c126, 0x45ca55,
    0x36d9ca, 0xd2a828, 0x8d61c2, 0x77c912, 0x142604, 0x9b4612,
    0xc459c4, 0x44c5c8, 0x91b24d, 0xf31700, 0xad43d4, 0xe54929,
    0x10d5fd, 0xfcbe00, 0xcc941e, 0xeece70, 0xf53e13, 0x80f1ec,
    0xc3e7b3, 0x28f8c7, 0x940593, 0x3e71c1, 0xb3092e, 0xf3450b,
    0x9c1288, 0x7b20ab, 0x9fb52e, 0xc29247, 0x2f327b, 0x6d550c,
    0x90a772, 0x1fe76b, 0x96cb31, 0x4a1679, 0xe27941, 0x89dff4,
    0x9794e8, 0x84e6e2, 0x973199, 0x6bed88, 0x365f5f, 0x0efdbb,
    0xb49a48, 0x6ca467, 0x427271, 0x325d8d, 0xb8159f, 0x09e5bc,
    0x25318d, 0x3974f7, 0x1c0530, 0x010c0d, 0x68084b, 0x58ee2c,
    0x90aa47, 0x02e774, 0x24d6bd, 0xa67df7, 0x72486e, 0xef169f,
    0xa6948e, 0xf691b4, 0x5153d1, 0xf20acf, 0x339820, 0x7e4bf5,
    0x6863b2, 0x5f3edd, 0x035d40, 0x7f8985, 0x295255, 0xc06437,
    0x10d86d, 0x324832, 0x754c5b, 0xd4714e, 0x6e5445, 0xc1090b,
    0x69f52a, 0xd56614, 0x9d0727, 0x50045d, 0xdb3bb4, 0xc576ea,
    0x17f987, 0x7d6b49, 0xba271d, 0x296996, 0xacccc6, 0x5414ad,
    0x6ae290, 0x89d988, 0x50722c, 0xbea404, 0x940777, 0x7030f3,
    0x27fc00, 0xa871ea, 0x49c266, 0x3de064, 0x83dd97, 0x973fa3,
    0xfd9443, 0x8c860d, 0xde4131, 0x9d3992, 0x8c70dd, 0xe7b717,
    0x3bdf08, 0x2b3715, 0xa0805c, 0x93805a, 0x921110, 0xd8e80f,
    0xaf806c, 0x4bffdb, 0x0f9038, 0x761859, 0x15a562, 0xbbcb61,
    0xb989c7, 0xbd4010, 0x04f2d2, 0x277549, 0xf6b6eb, 0xbb22db,
    0xaa140a, 0x2f2689, 0x768364, 0x333b09, 0x1a940e, 0xaa3a51,
    0xc2a31d, 0xaeedaf, 0x12265c, 0x4dc26d, 0x9c7a2d, 0x9756c0,
    0x833f03, 0xf6f009, 0x8c402b, 0x99316d, 0x07b439, 0x15200c,
    0x5bc3d8, 0xc492f5, 0x4badc6, 0xa5ca4e, 0xcd37a7, 0x36a9e6,
    0x9492ab, 0x6842dd, 0xde6319, 0xef8c76, 0x528b68, 0x37dbfc,
    0xaba1ae, 0x3115df, 0xa1ae00, 0xdafb0c, 0x664d64, 0xb705ed,
    0x306529, 0xbf5657, 0x3aff47, 0xb9f96a, 0xf3be75, 0xdf9328,
    0x3080ab, 0xf68c66, 0x15cb04, 0x0622fa, 0x1de4d9, 0xa4b33d,
    0x8f1b57, 0x09cd36, 0xe9424e, 0xa4be13, 0xb52333, 0x1aaaf0,
    0xa8654f, 0xa5c1d2, 0x0f3f0b, 0xcd785b, 0x76f923, 0x048b7b,
    0x721789, 0x53a6c6, 0xe26e6f, 0x00ebef, 0x584a9b, 0xb7dac4,
    0xba66aa, 0xcfcf76, 0x1d02d1, 0x2df1b1, 0xc1998c, 0x77adc3,
    0xda4886, 0xa05df7, 0xf480c6, 0x2ff0ac, 0x9aecdd, 0xbc5c3f,
    0x6dded0, 0x1fc790, 0xb6db2a, 0x3a25a3, 0x9aaf00, 0x9353ad,
    0x0457b6, 0xb42d29, 0x7e804b, 0xa707da, 0x0eaa76, 0xa1597b,
    0x2a1216, 0x2db7dc, 0xfde5fa, 0xfedb89, 0xfdbe89, 0x6c76e4,
    0xfca906, 0x70803e, 0x156e85, 0xff87fd, 0x073e28, 0x336761,
    0x86182a, 0xeabd4d, 0xafe7b3, 0x6e6d8f, 0x396795, 0x5bbf31,
    0x48d784, 0x16df30, 0x432dc7, 0x356125, 0xce70c9, 0xb8cb30,
    0xfd6cbf, 0xa200a4, 0xe46c05, 0xa0dd5a, 0x476f21, 0xd21262,
    0x845cb9, 0x496170, 0xe0566b, 0x015299, 0x375550, 0xb7d51e,
    0xc4f133, 0x5f6e13, 0xe4305d, 0xa92e85, 0xc3b21d, 0x3632a1,
    0xa4b708, 0xd4b1ea, 0x21f716, 0xe4698f, 0x77ff27, 0x80030c,
    0x2d408d, 0xa0cd4f, 0x99a520, 0xd3a2b3, 0x0a5d2f, 0x42f9b4,
    0xcbda11, 0xd0be7d, 0xc1db9b, 0xbd17ab, 0x81a2ca, 0x5c6a08,
    0x17552e, 0x550027, 0xf0147f, 0x8607e1, 0x640b14, 0x8d4196,
    0xdebe87, 0x2afdda, 0xb6256b, 0x34897b, 0xfef305, 0x9ebfb9,
    0x4f6a68, 0xa82a4a, 0x5ac44f, 0xbcf82d, 0x985ad7, 0x95c7f4,
    0x8d4d0d, 0xa63a20, 0x5f57a4, 0xb13f14, 0x953880, 0x0120cc,
    0x86dd71, 0xb6dec9, 0xf560bf, 0x11654d, 0x6b0701, 0xacb08c,
    0xd0c0b2, 0x485551, 0x0efb1e, 0xc37295, 0x3b06a3, 0x3540c0,
    0x7bdc06, 0xcc45e0, 0xfa294e, 0xc8cad6, 0x41f3e8, 0xde647c,
    0xd8649b, 0x31bed9, 0xc397a4, 0xd45877, 0xc5e369, 0x13daf0,
    0x3c3aba, 0x461846, 0x5f7555, 0xf5bdd2, 0xc6926e, 0x5d2eac,
    0xed440e, 0x423e1c, 0x87c461, 0xe9fd29, 0xf3d6e7, 0xca7c22,
    0x35916f, 0xc5e008, 0x8dd7ff, 0xe26a6e, 0xc6fdb0, 0xc10893,
    0x745d7c, 0xb2ad6b, 0x9d6ecd, 0x7b723e, 0x6a11c6, 0xa9cff7,
    0xdf7329, 0xbac9b5, 0x5100b7, 0x0db2e2, 0x24ba74, 0x607de5,
    0x8ad874, 0x2c150d, 0x0c1881, 0x94667e, 0x162901, 0x767a9f,
    0xbefdfd, 0xef4556, 0x367ed9, 0x13d9ec, 0xb9ba8b, 0xfc97c4,
    0x27a831, 0xc36ef1, 0x36c594, 0x56a8d8, 0xb5a8b4, 0x0ecccf,
    0x2d8912, 0x34576f, 0x89562c, 0xe3ce99, 0xb920d6, 0xaa5e6b,
    0x9c2a3e, 0xcc5f11, 0x4a0bfd, 0xfbf4e1, 0x6d3b8e, 0x2c86e2,
    0x84d4e9, 0xa9b4fc, 0xd1eeef, 0xc9352e, 0x61392f, 0x442138,
    0xc8d91b, 0x0afc81, 0x6a4afb, 0xd81c2f, 0x84b453, 0x8c994e,
    0xcc2254, 0xdc552a, 0xd6c6c0, 0x96190b, 0xb8701a, 0x649569,
    0x605a26, 0xee523f, 0x0f117f, 0x11b5f4, 0xf5cbfc, 0x2dbc34,
    0xeebc34, 0xcc5de8, 0x605edd, 0x9b8e67, 0xef3392, 0xb817c9,
    0x9b5861, 0xbc57e1, 0xc68351, 0x103ed8, 0x4871dd, 0xdd1c2d,
    0xa118af, 0x462c21, 0xd7f359, 0x987ad9, 0xc0549e, 0xfa864f,
    0xfc0656, 0xae79e5, 0x362289, 0x22ad38, 0xdc9367, 0xaae855,
    0x382682, 0x9be7ca, 0xa40d51, 0xb13399, 0x0ed7a9, 0x480569,
    0xf0b265, 0xa7887f, 0x974c88, 0x36d1f9, 0xb39221, 0x4a827b,
    0x21cf98, 0xdc9f40, 0x5547dc, 0x3a74e1, 0x42eb67, 0xdf9dfe,
    0x5fd45e, 0xa4677b, 0x7aacba, 0xa2f655, 0x23882b, 0x55ba41,
    0x086e59, 0x862a21, 0x834739, 0xe6e389, 0xd49ee5, 0x40fb49,
    0xe956ff, 0xca0f1c, 0x8a59c5, 0x2bfa94, 0xc5c1d3, 0xcfc50f,
    0xae5adb, 0x86c547, 0x624385, 0x3b8621, 0x94792c, 0x876110,
    0x7b4c2a, 0x1a2c80, 0x12bf43, 0x902688, 0x893c78, 0xe4c4a8,
    0x7bdbe5, 0xc23ac4, 0xeaf426, 0x8a67f7, 0xbf920d, 0x2ba365,
    0xb1933d, 0x0b7cbd, 0xdc51a4, 0x63dd27, 0xdde169, 0x19949a,
    0x9529a8, 0x28ce68, 0xb4ed09, 0x209f44, 0xca984e, 0x638270,
    0x237c7e, 0x32b90f, 0x8ef5a7, 0xe75614, 0x08f121, 0x2a9db5,
    0x4d7e6f, 0x5119a5, 0xabf9b5, 0xd6df82, 0x61dd96, 0x023616,
    0x9f3ac4, 0xa1a283, 0x6ded72, 0x7a8d39, 0xa9b882, 0x5c326b,
    0x5b2746, 0xed3400, 0x7700d2, 0x55f4fc, 0x4d5901, 0x8071e0,
};

const pio2_64: [8]f64 = .{
    1.57079625129699707031e+00, // 0x3ff921fb, 0x40000000
    7.54978941586159635335e-08, // 0x3e74442d, 0x00000000
    5.39030252995776476554e-15, // 0x3cf84698, 0x80000000
    3.28200341580791294123e-22, // 0x3b78cc51, 0x60000000
    1.27065575308067607349e-29, // 0x39f01b83, 0x80000000
    1.22933308981111328932e-36, // 0x387a2520, 0x40000000
    2.73370053816464559624e-44, // 0x36e38222, 0x80000000
    2.16741683877804819444e-51, // 0x3569f31d, 0x00000000
};
