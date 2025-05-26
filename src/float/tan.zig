const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const roundeven = @import("roundeven.zig");
const utan = @import("utan.zig");
const branred = @import("branred.zig");
const dla = @import("dla.zig");
const ldbl128 = @import("ldbl128.zig");
const rem_pio2 = @import("rem_pio2.zig");
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn tan(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.tan: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    switch (EnsureFloat(@TypeOf(x))) {
        f16 => return scast(f16, tan32(scast(f32, x))),
        f32 => {
            // glibc/sysdeps/ieee754/flt-32/s_tanf.c
            return tan32(scast(f32, x));
        },
        f64 => {
            // glibc/sysdeps/ieee754/dbl-64/s_tan.c
            return tan64(scast(f64, x));
        },
        f80 => return scast(f80, tan128(scast(f128, x))),
        f128 => {
            // glibc/sysdeps/ieee754/ldbl-128/s_tanl.c
            return tan128(scast(f128, x));
        },
        else => unreachable,
    }
}

// argument reduction
// for |z| < 2^28, return r such that 2/pi*x = q + r
inline fn rltl(z: f32, q: *i32) f64 {
    const x: f64 = scast(f64, z);
    const idl: f64 = -0x1.b1bbead603d8bp-32 * x;
    const idh: f64 = 0x1.45f306ep-1 * x;
    const id: f64 = roundeven.roundeven_finite(idh);
    q.* = scast(i32, id);
    return (idh - id) + idl;
}

// argument reduction
// same as rltl, but for |x| >= 2^28
fn rbig(u: u32, q: *i32) f64 {
    const ipi: [4]u64 = .{
        0xfe5163abdebbc562, 0xdb6295993c439041,
        0xfc2757d1f534ddc0, 0xa2f9836e4e441529,
    };

    const e: i32 = scast(i32, (u >> 23) & 0xff);
    const m: u64 = scast(u64, (u & (~@as(u32, 0) >> 9)) | 1 << 23);
    const p0: u128 = scast(u128, m) * scast(u128, ipi[0]);
    var p1: u128 = scast(u128, m) * scast(u128, ipi[1]);
    p1 += p0 >> 64;
    var p2: u128 = scast(u128, m) * scast(u128, ipi[2]);
    p2 += p1 >> 64;
    var p3: u128 = scast(u128, m) * scast(u128, ipi[3]);
    p3 += p2 >> 64;
    const p3h: u64 = scast(u64, p3 >> 64);
    const p3l: u64 = scast(u64, p3 & 0xffffffffffffffff);
    const p2l: u64 = scast(u64, p2 & 0xffffffffffffffff);
    const p1l: u64 = scast(u64, p1 & 0xffffffffffffffff);
    const k: i32 = e - 127;
    const s: i32 = k - 23;
    // in ctanf(), rbig() is called in the case 127+28 <= e < 0xff
    // thus 155 <= e <= 254, which yields 28 <= k <= 127 and 5 <= s <= 104
    var i: i32 = undefined;
    var a: i64 = undefined;
    if (s < 64) {
        i = @bitCast(@as(u32, @truncate((p3h << @as(u6, @intCast(s)) | p3l >> @as(u6, @intCast(64 - s))) & 0xffffffff)));
        a = @bitCast(p3l << @as(u6, @intCast(s)) | p2l >> @as(u6, @intCast(64 - s)));
    } else if (s == 64) {
        i = @bitCast(@as(u32, @truncate(p3l & 0xffffffff)));
        a = @bitCast(p2l);
    } else { // s > 64
        i = @bitCast(@as(u32, @truncate((p3l << @as(u6, @intCast(s - 64)) | p2l >> @as(u6, @intCast(128 - s))) & 0xffffffff)));
        a = @bitCast(p2l << @as(u6, @intCast(s - 64)) | p1l >> @as(u6, @intCast(128 - s)));
    }

    var sgn: i32 = @bitCast(u);
    sgn >>= 31;
    const sm: i32 = @truncate(a >> 63);
    i -= sm;
    const z: f64 = scast(f64, a ^ sgn) * 0x1p-64;
    i = (i ^ sgn) - sgn;
    q.* = i;
    return z;
}

fn tan32(x: f32) f32 {
    const t: u32 = @bitCast(x);
    const e: i32 = scast(i32, (t >> 23) & 0xff);

    var i: i32 = undefined;
    var z: f64 = undefined;
    if (e < 127 + 28) { // |x| < 2^28
        @branchHint(.likely);
        if (e < 115) {
            @branchHint(.unlikely);
            if (e < 102) {
                @branchHint(.unlikely);
                return @mulAdd(f32, x, float.abs(x), x);
            }

            const x2: f32 = x * x;
            return @mulAdd(f32, x, 0x1.555556p-2 * x2, x);
        }

        z = rltl(x, &i);
    } else if (e < 0xff) {
        z = rbig(t, &i);
    } else {
        if ((t << 9) != 0)
            return x + x; // nan

        return (x - x) / (x - x);
    }

    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    const cn: [4]f64 = .{
        0x1.921fb54442d18p+0, -0x1.fd226e573289fp-2,
        0x1.b7a60c8dac9f6p-6, -0x1.725beb40f33e5p-13,
    };
    const cd: [4]f64 = .{
        0x1p+0,               -0x1.2395347fb829dp+0,
        0x1.2313660f29c36p-3, -0x1.9a707ab98d1c1p-9,
    };
    const s: [2]f64 = .{ 0, 1 };
    var n: f64 = cn[0] + z2 * cn[1];
    const n2: f64 = cn[2] + z2 * cn[3];
    n += z4 * n2;
    var d: f64 = cd[0] + z2 * cd[1];
    const d2: f64 = cd[2] + z2 * cd[3];
    d += z4 * d2;
    n *= z;
    const s0: f64 = s[@intCast(i & 1)];
    const s1: f64 = s[@intCast(1 - (i & 1))];
    const r1: f64 = (n * s1 - d * s0) / (n * s0 + d * s1);
    const tail: u64 = (@as(u64, @bitCast(r1)) + 7) & (~@as(u64, 0) >> 35);
    if (tail <= 14) {
        @branchHint(.unlikely);
        const st: [8]struct {
            arg: f32,
            rh: f32,
            rl: f32,
        } = .{
            .{ .arg = 0x1.143ec4p+0, .rh = 0x1.ddf9f6p+0, .rl = -0x1.891d24p-52 },
            .{ .arg = 0x1.ada6aap+27, .rh = 0x1.e80304p-3, .rl = 0x1.419f46p-58 },
            .{ .arg = 0x1.af61dap+48, .rh = 0x1.60d1c8p-2, .rl = -0x1.2d6c3ap-55 },
            .{ .arg = 0x1.0088bcp+52, .rh = 0x1.ca1edp+0, .rl = 0x1.f6053p-53 },
            .{ .arg = 0x1.f90dfcp+72, .rh = 0x1.597f9cp-1, .rl = 0x1.925978p-53 },
            .{ .arg = 0x1.cc4e22p+85, .rh = -0x1.f33584p+1, .rl = 0x1.d7254ap-51 },
            .{ .arg = 0x1.a6ce12p+86, .rh = -0x1.c5612ep-1, .rl = -0x1.26c33ep-53 },
            .{ .arg = 0x1.6a0b76p+102, .rh = -0x1.e42a1ep+0, .rl = -0x1.1dc906p-52 },
        };
        const ax: u32 = t & (~@as(u32, 0) >> 1);
        const sgn: u32 = t >> 31;
        var j: i32 = 0;
        while (j < 8) {
            if (st[@intCast(j)].arg == scast(f32, ax)) {
                @branchHint(.unlikely);

                if (sgn != 0) {
                    return -st[@intCast(j)].rh - st[@intCast(j)].rl;
                } else {
                    return st[@intCast(j)].rh + st[@intCast(j)].rl;
                }
            }

            j += 1;
        }
    }

    return scast(f32, r1);
}

// tan with max ULP of ~0.619 based on random sampling.
fn tan64(x: f64) f64 {
    // x=+-INF, x=NaN
    const num: [2]u32 = @bitCast(x);
    const ux: u32 = num[utan.HIGH_HALF];
    if ((ux & 0x7ff00000) == 0x7ff00000) {
        return x - x;
    }

    const w: f64 = if (x < 0) -x else x;

    // (I) The case abs(x) <= 1.259e-8
    if (w <= @as(f64, @bitCast(utan.g1))) {
        if (w < std.math.floatMin(f64)) {
            const vw: f64 = w * w;
            std.mem.doNotOptimizeAway(vw);
        }

        return x;
    }

    // (II) The case 1.259e-8 < abs(x) <= 0.0608
    if (w <= @as(f64, @bitCast(utan.g2))) {
        const x2: f64 = x * x;

        var t2: f64 = @as(f64, @bitCast(utan.d9)) + x2 * @as(f64, @bitCast(utan.d11));
        t2 = @as(f64, @bitCast(utan.d7)) + x2 * t2;
        t2 = @as(f64, @bitCast(utan.d5)) + x2 * t2;
        t2 = @as(f64, @bitCast(utan.d3)) + x2 * t2;
        t2 *= x * x2;

        // Max ULP is 0.504.
        return x + t2;
    }

    // (III) The case 0.0608 < abs(x) <= 0.787
    if (w <= @as(f64, @bitCast(utan.g3))) {
        const i: u32 = scast(u32, @as(f64, @bitCast(utan.mfftnhf)) + 256 * w);
        const z: f64 = w - @as(f64, @bitCast(utan.xfg[i][0]));
        const z2: f64 = z * z;
        const s: f64 = if (x < 0) -1 else 1;
        const pz: f64 = z + z * z2 * (@as(f64, @bitCast(utan.e0)) + z2 * @as(f64, @bitCast(utan.e1)));
        const fi: f64 = @as(f64, @bitCast(utan.xfg[i][1]));
        const gi: f64 = @as(f64, @bitCast(utan.xfg[i][2]));
        const t2: f64 = pz * (gi + fi) / (gi - pz);
        const y: f64 = fi + t2;

        // Max ULP is 0.60.
        return (s * y);
    }

    // (---) The case 0.787 < abs(x) <= 25
    if (w <= @as(f64, @bitCast(utan.g4))) {
        // Range reduction by algorithm i
        const t: f64 = (x * @as(f64, @bitCast(utan.hpinv)) + @as(f64, @bitCast(utan.toint)));
        const xn: f64 = t - @as(f64, @bitCast(utan.toint));
        const v: [2]u32 = @bitCast(t);
        var t1: f64 = (x - xn * @as(f64, @bitCast(utan.mp1))) - xn * @as(f64, @bitCast(utan.mp2));
        const n: u32 = v[utan.LOW_HALF] & 0x00000001;
        var da: f64 = xn * @as(f64, @bitCast(utan.mp3));
        const a: f64 = t1 - da;
        da = (t1 - a) - da;

        var ya: f64 = undefined;
        var yya: f64 = undefined;
        var sy: f64 = undefined;
        if (a < 0) {
            ya = -a;
            yya = -da;
            sy = -1;
        } else {
            ya = a;
            yya = da;
            sy = 1;
        }

        // (VI) The case 0.787 < abs(x) <= 25,    0 < abs(y) <= 0.0608
        if (ya <= @as(f64, @bitCast(utan.gy2))) {
            const a2: f64 = a * a;
            var t2: f64 = @as(f64, @bitCast(utan.d9)) + a2 * @as(f64, @bitCast(utan.d11));
            t2 = @as(f64, @bitCast(utan.d7)) + a2 * t2;
            t2 = @as(f64, @bitCast(utan.d5)) + a2 * t2;
            t2 = @as(f64, @bitCast(utan.d3)) + a2 * t2;
            t2 = da + a * a2 * t2;

            if (n != 0) {
                // -cot
                var b: f64 = undefined;
                var db: f64 = undefined;
                dla.eadd(a, t2, &b, &db);
                var c: f64 = undefined;
                var dc: f64 = undefined;
                var t3: f64 = undefined;
                var t4: f64 = undefined;
                dla.div2(1, 0, b, db, &c, &dc, &t1, &t2, &t3, &t4);
                const y: f64 = c + dc;

                // Max ULP is 0.506.
                return -y;
            } else {
                // tan
                const y: f64 = a + t2;

                // Max ULP is 0.506.
                return y;
            }
        }

        // (VII) The case 0.787 < abs(x) <= 25,    0.0608 < abs(y) <= 0.787
        const i: u32 = scast(u32, @as(f64, @bitCast(utan.mfftnhf)) + 256 * ya);
        const z: f64 = (ya - @as(f64, @bitCast(utan.xfg[i][0]))) + yya;
        const z2: f64 = z * z;
        const pz: f64 = z + z * z2 * (@as(f64, @bitCast(utan.e0)) + z2 * @as(f64, @bitCast(utan.e1)));
        const fi: f64 = @as(f64, @bitCast(utan.xfg[i][1]));
        const gi: f64 = @as(f64, @bitCast(utan.xfg[i][2]));

        if (n != 0) {
            // -cot
            const t2: f64 = pz * (fi + gi) / (fi + pz);
            const y: f64 = gi - t2;

            // Max ULP is 0.62.
            return -sy * y;
        } else {
            // tan
            const t2: f64 = pz * (gi + fi) / (gi - pz);
            const y: f64 = fi + t2;

            // Max ULP is 0.62.
            return sy * y;
        }
    }

    // (---) The case 25 < abs(x) <= 1e8
    if (w <= @as(f64, @bitCast(utan.g5))) {
        // Range reduction by algorithm ii
        var t: f64 = (x * @as(f64, @bitCast(utan.hpinv)) + @as(f64, @bitCast(utan.toint)));
        const xn: f64 = t - @as(f64, @bitCast(utan.toint));
        const v: [2]u32 = @bitCast(t);
        var t1: f64 = (x - xn * @as(f64, @bitCast(utan.mp1))) - xn * @as(f64, @bitCast(utan.mp2));
        const n: u32 = v[utan.LOW_HALF] & 0x00000001;
        var da: f64 = xn * @as(f64, @bitCast(utan.pp3));
        t = t1 - da;
        da = (t1 - t) - da;
        t1 = xn * @as(f64, @bitCast(utan.pp4));
        var a: f64 = t - t1;
        da = ((t - a) - t1) + da;
        var t2: f64 = undefined;
        dla.eadd(a, da, &t1, &t2);
        a = t1;
        da = t2;

        var ya: f64 = undefined;
        var yya: f64 = undefined;
        var sy: f64 = undefined;
        if (a < 0) {
            ya = -a;
            yya = -da;
            sy = -1;
        } else {
            ya = a;
            yya = da;
            sy = 1;
        }

        // (VIII) The case 25 < abs(x) <= 1e8,    0 < abs(y) <= 0.0608
        if (ya <= @as(f64, @bitCast(utan.gy2))) {
            const a2: f64 = a * a;
            t2 = @as(f64, @bitCast(utan.d9)) + a2 * @as(f64, @bitCast(utan.d11));
            t2 = @as(f64, @bitCast(utan.d7)) + a2 * t2;
            t2 = @as(f64, @bitCast(utan.d5)) + a2 * t2;
            t2 = @as(f64, @bitCast(utan.d3)) + a2 * t2;
            t2 = da + a * a2 * t2;

            if (n != 0) {
                // -cot
                var b: f64 = undefined;
                var db: f64 = undefined;
                dla.eadd(a, t2, &b, &db);
                var c: f64 = undefined;
                var dc: f64 = undefined;
                var t3: f64 = undefined;
                var t4: f64 = undefined;
                dla.div2(1, 0, b, db, &c, &dc, &t1, &t2, &t3, &t4);
                const y: f64 = c + dc;

                // Max ULP is 0.506.
                return -y;
            } else {
                // tan
                const y: f64 = a + t2;

                // Max ULP is 0.506.
                return y;
            }
        }

        // (IX) The case 25 < abs(x) <= 1e8,    0.0608 < abs(y) <= 0.787
        const i: u32 = scast(u32, @as(f64, @bitCast(utan.mfftnhf)) + 256 * ya);
        const z: f64 = (ya - @as(f64, @bitCast(utan.xfg[i][0]))) + yya;
        const z2: f64 = z * z;
        const pz: f64 = z + z * z2 * (@as(f64, @bitCast(utan.e0)) + z2 * @as(f64, @bitCast(utan.e1)));
        const fi: f64 = @as(f64, @bitCast(utan.xfg[i][1]));
        const gi: f64 = @as(f64, @bitCast(utan.xfg[i][2]));

        if (n != 0) {
            // -cot
            t2 = pz * (fi + gi) / (fi + pz);
            const y: f64 = gi - t2;

            // Max ULP is 0.62.
            return -sy * y;
        } else {
            // tan
            t2 = pz * (gi + fi) / (gi - pz);
            const y: f64 = fi + t2;

            // Max ULP is 0.62.
            return sy * y;
        }
    }

    // (---) The case 1e8 < abs(x) < 2**1024
    // Range reduction by algorithm iii
    var a: f64 = undefined;
    var da: f64 = undefined;
    const n: i32 = (branred.branred(x, &a, &da)) & 0x00000001;
    var t1: f64 = undefined;
    var t2: f64 = undefined;
    dla.eadd(a, da, &t1, &t2);
    a = t1;
    da = t2;

    var ya: f64 = undefined;
    var yya: f64 = undefined;
    var sy: f64 = undefined;
    if (a < 0) {
        ya = -a;
        yya = -da;
        sy = -1;
    } else {
        ya = a;
        yya = da;
        sy = 1;
    }

    // (X) The case 1e8 < abs(x) < 2**1024,    0 < abs(y) <= 0.0608
    if (ya <= @as(f64, @bitCast(utan.gy2))) {
        const a2: f64 = a * a;
        t2 = @as(f64, @bitCast(utan.d9)) + a2 * @as(f64, @bitCast(utan.d11));
        t2 = @as(f64, @bitCast(utan.d7)) + a2 * t2;
        t2 = @as(f64, @bitCast(utan.d5)) + a2 * t2;
        t2 = @as(f64, @bitCast(utan.d3)) + a2 * t2;
        t2 = da + a * a2 * t2;
        if (n != 0) {
            // -cot
            var b: f64 = undefined;
            var db: f64 = undefined;
            dla.eadd(a, t2, &b, &db);
            var c: f64 = undefined;
            var dc: f64 = undefined;
            var t3: f64 = undefined;
            var t4: f64 = undefined;
            dla.div2(1.0, 0.0, b, db, &c, &dc, &t1, &t2, &t3, &t4);
            const y: f64 = c + dc;

            // Max ULP is 0.506.
            return -y;
        } else {
            // tan
            const y: f64 = a + t2;

            // Max ULP is 0.507.
            return y;
        }
    }

    // (XI) The case 1e8 < abs(x) < 2**1024,    0.0608 < abs(y) <= 0.787
    const i: u32 = scast(u32, @as(f64, @bitCast(utan.mfftnhf)) + 256 * ya);
    const z: f64 = (ya - @as(f64, @bitCast(utan.xfg[i][0]))) + yya;
    const z2: f64 = z * z;
    const pz: f64 = z + z * z2 * (@as(f64, @bitCast(utan.e0)) + z2 * @as(f64, @bitCast(utan.e1)));
    const fi: f64 = @as(f64, @bitCast(utan.xfg[i][1]));
    const gi: f64 = @as(f64, @bitCast(utan.xfg[i][2]));

    if (n != 0) {
        // -cot
        t2 = pz * (fi + gi) / (fi + pz);
        const y: f64 = gi - t2;

        // Max ULP is 0.62.
        return -sy * y;
    } else {
        // tan
        t2 = pz * (gi + fi) / (gi - pz);
        const y: f64 = fi + t2;

        // Max ULP is 0.62.
        return sy * y;
    }
}

fn kernel_tan128(x: f128, y: f128, iy: i32) f128 {
    const pio4hi: f128 = 7.8539816339744830961566084581987569936977e-1;
    const pio4lo: f128 = 2.1679525325309452561992610065108379921906e-35;
    // tan x = x + x^3 / 3 + x^5 T(x^2)/U(x^2)
    // 0 <= x <= 0.6743316650390625
    // Peak relative error 8.0e-36
    const TH: f128 = 3.333333333333333333333333333333333333333e-1;
    const T0: f128 = -1.813014711743583437742363284336855889393e7;
    const T1: f128 = 1.320767960008972224312740075083259247618e6;
    const T2: f128 = -2.626775478255838182468651821863299023956e4;
    const T3: f128 = 1.764573356488504935415411383687150199315e2;
    const T4: f128 = -3.333267763822178690794678978979803526092e-1;
    const U0: f128 = -1.359761033807687578306772463253710042010e8;
    const U1: f128 = 6.494370630656893175666729313065113194784e7;
    const U2: f128 = -4.180787672237927475505536849168729386782e6;
    const U3: f128 = 8.031643765106170040139966622980914621521e4;
    const U4: f128 = -5.323131271912475695157127875560667378597e2;

    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const ix: i32 = @bitCast(u.w0 & 0x7fffffff);
    if (ix < 0x3fc60000) { // x < 2**-57
        if (scast(i32, x) == 0) { // generate inexact
            if ((ix | @as(i32, @bitCast(u.w1)) | @as(i32, @bitCast(u.w2)) | @as(i32, @bitCast(u.w3)) | (iy + 1)) == 0) {
                return 1 / float.abs(x);
            } else if (iy == 1) {
                if (float.abs(x) < std.math.floatMin(f128)) {
                    const vx: f128 = x * x;
                    std.mem.doNotOptimizeAway(vx);
                }

                return x;
            } else {
                return -1 / x;
            }
        }
    }

    var xx: f128 = x;
    var yy: f128 = y;
    var z: f128 = undefined;
    var w: f128 = undefined;
    var sign: i32 = 0;
    if (ix >= 0x3ffe5942) { // |x| >= 0.6743316650390625
        if ((u.w0 & 0x80000000) != 0) {
            xx = -xx;
            yy = -yy;
            sign = -1;
        } else {
            sign = 1;
        }

        z = pio4hi - xx;
        w = pio4lo - yy;
        xx = z + w;
        yy = 0;
    }
    z = xx * xx;
    var r: f128 = T0 + z * (T1 + z * (T2 + z * (T3 + z * T4)));
    var v: f128 = U0 + z * (U1 + z * (U2 + z * (U3 + z * (U4 + z))));
    r = r / v;

    var s: f128 = z * xx;
    r = yy + z * (s * r + yy);
    r += TH * s;
    w = xx + r;
    if (ix >= 0x3ffe5942) {
        v = scast(f128, iy);
        w = (v - 2.0 * (xx - (w * w / (w + v) - r)));
        if (sign < 0)
            w = -w;

        return w;
    }
    if (iy == 1) {
        return w;
    } else { // if allow error up to 2 ulp, simply return -1.0/(x+r) here
        // compute -1.0/(x+r) accurately
        var uu: ldbl128.ieee_f128_shape32 = @bitCast(w);
        uu.w2 = 0;
        uu.w3 = 0;
        v = r - (@as(f128, @bitCast(uu)) - xx); // uu+v = r+x
        z = -1.0 / w;
        u = @bitCast(z);
        u.w2 = 0;
        u.w3 = 0;
        s = 1.0 + @as(f128, @bitCast(u)) * @as(f128, @bitCast(uu));
        return @as(f128, @bitCast(u)) + z * (s + @as(f128, @bitCast(u)) * v);
    }
}

fn tan128(x: f128) f128 {
    // High word of x.
    var ix: i64 = undefined;
    ldbl128.getMsw(&ix, x);

    // |x| ~< pi/4
    ix &= 0x7fffffffffffffff;
    if (ix <= 0x3ffe921fb54442d1) {
        return kernel_tan128(x, 0, 1);
    } else if (ix >= 0x7fff000000000000) { // tanl(Inf or NaN) is NaN
        return x - x; // NaN
    } else { // argument reduction needed
        var y: [2]f128 = undefined;
        const n: i32 = rem_pio2.rem_pio2_128(x, &y);
        return kernel_tan128(y[0], y[1], 1 - ((n & 1) << 1)); //   1 -- n even, -1 -- n odd
    }
}
