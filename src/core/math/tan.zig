const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const roundeven = @import("roundeven.zig");
const utan = @import("utan.zig");
const branred = @import("branred.zig");
const dla = @import("dla.zig");
const ldbl128 = @import("ldbl128.zig");
const rem_pio2 = @import("rem_pio2.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn tan(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return tan(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, tan32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_tanf.c
                    return tan32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_tan.c
                    return tan64(x);
                },
                f80 => return cast(f80, tan128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_tanl.c
                    return tan128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

// argument reduction
// for |z| < 2^28, return r such that 2/pi*x = q + r
inline fn rltl(z: f32, q: *i32) f64 {
    const x: f64 = cast(f64, z, .{});
    const idl: f64 = -0x1.b1bbead603d8bp-32 * x;
    const idh: f64 = 0x1.45f306ep-1 * x;
    const id: f64 = roundeven.roundeven_finite(idh);
    q.* = cast(i32, id, .{});
    return (idh - id) + idl;
}

// argument reduction
// same as rltl, but for |x| >= 2^28
fn rbig(u: u32, q: *i32) f64 {
    const ipi: [4]u64 = .{
        0xfe5163abdebbc562, 0xdb6295993c439041,
        0xfc2757d1f534ddc0, 0xa2f9836e4e441529,
    };

    const e: i32 = cast(i32, (u >> 23) & 0xff, .{});
    const m: u64 = cast(u64, (u & (~@as(u32, 0) >> 9)) | 1 << 23, .{});
    const p0: u128 = cast(u128, m, .{}) * cast(u128, ipi[0], .{});
    var p1: u128 = cast(u128, m, .{}) * cast(u128, ipi[1], .{});
    p1 += p0 >> 64;
    var p2: u128 = cast(u128, m, .{}) * cast(u128, ipi[2], .{});
    p2 += p1 >> 64;
    var p3: u128 = cast(u128, m, .{}) * cast(u128, ipi[3], .{});
    p3 += p2 >> 64;
    const p3h: u64 = cast(u64, p3 >> 64, .{});
    const p3l: u64 = cast(u64, p3 & 0xffffffffffffffff, .{});
    const p2l: u64 = cast(u64, p2 & 0xffffffffffffffff, .{});
    const p1l: u64 = cast(u64, p1 & 0xffffffffffffffff, .{});
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
    const z: f64 = cast(f64, a ^ sgn, .{}) * 0x1p-64;
    i = (i ^ sgn) - sgn;
    q.* = i;
    return z;
}

fn tan32(x: f32) f32 {
    const t: u32 = @bitCast(x);
    const e: i32 = cast(i32, (t >> 23) & 0xff, .{});

    var i: i32 = undefined;
    var z: f64 = undefined;
    if (e < 127 + 28) { // |x| < 2^28
        @branchHint(.likely);
        if (e < 115) {
            @branchHint(.unlikely);
            if (e < 102) {
                @branchHint(.unlikely);
                return @mulAdd(f32, x, math.abs(x), x);
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
            if (st[@intCast(j)].arg == cast(f32, ax, .{})) {
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

    return cast(f32, r1, .{});
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
        const i: u32 = cast(u32, @as(f64, @bitCast(utan.mfftnhf)) + 256 * w, .{});
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
        const i: u32 = cast(u32, @as(f64, @bitCast(utan.mfftnhf)) + 256 * ya, .{});
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
        const i: u32 = cast(u32, @as(f64, @bitCast(utan.mfftnhf)) + 256 * ya, .{});
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
    const i: u32 = cast(u32, @as(f64, @bitCast(utan.mfftnhf)) + 256 * ya, .{});
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
        if (cast(i32, x, .{}) == 0) { // generate inexact
            if ((ix | @as(i32, @bitCast(u.w1)) | @as(i32, @bitCast(u.w2)) | @as(i32, @bitCast(u.w3)) | (iy + 1)) == 0) {
                return 1 / math.abs(x);
            } else if (iy == 1) {
                if (math.abs(x) < std.math.floatMin(f128)) {
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
        v = cast(f128, iy, .{});
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

test tan {
    try std.testing.expectEqual(0x0p+0, tan(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, tan(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1p+0, tan(@as(f32, 0xc.90fdbp-4)));
    try std.testing.expectEqual(0xf.fffffp-4, tan(@as(f32, 0xc.90fdap-4)));
    try std.testing.expectEqual(-0x1.5d1494p+24, tan(@as(f32, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0xc.a1bdap+20, tan(@as(f32, 0x1.921fb4p+0)));
    try std.testing.expectEqual(0x1.5d1494p+24, tan(@as(f32, -0x1.921fb6p+0)));
    try std.testing.expectEqual(-0xc.a1bdap+20, tan(@as(f32, -0x1.921fb4p+0)));
    try std.testing.expectEqual(0xe.e7d1bp-4, tan(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(-0xc.17b0cp-8, tan(@as(f32, 0x2p+64)));
    try std.testing.expectEqual(0xc.17b0cp-8, tan(@as(f32, -0x2p+64)));
    try std.testing.expectEqual(0x2p-28, tan(@as(f32, 0x2p-28)));
    try std.testing.expectEqual(-0x2p-28, tan(@as(f32, -0x2p-28)));
    try std.testing.expectEqual(0xf.fe04dp-4, tan(@as(f32, 0xc.9p-4)));
    try std.testing.expectEqual(0xf.ff04bp-4, tan(@as(f32, 0xc.908p-4)));
    try std.testing.expectEqual(0xf.ff84bp-4, tan(@as(f32, 0xc.90cp-4)));
    try std.testing.expectEqual(0xf.ffc4bp-4, tan(@as(f32, 0xc.90ep-4)));
    try std.testing.expectEqual(0xf.ffe4bp-4, tan(@as(f32, 0xc.90fp-4)));
    try std.testing.expectEqual(0xf.fff4bp-4, tan(@as(f32, 0xc.90f8p-4)));
    try std.testing.expectEqual(0xf.fffcbp-4, tan(@as(f32, 0xc.90fcp-4)));
    try std.testing.expectEqual(0xf.fffebp-4, tan(@as(f32, 0xc.90fdp-4)));
    try std.testing.expectEqual(0xf.ffffbp-4, tan(@as(f32, 0xc.90fd8p-4)));
    try std.testing.expectEqual(0xf.fffffp-4, tan(@as(f32, 0xc.90fdap-4)));
    try std.testing.expectEqual(0x1.01e21p+0, tan(@as(f32, 0xc.ap-4)));
    try std.testing.expectEqual(0x1.00e0aep+0, tan(@as(f32, 0xc.98p-4)));
    try std.testing.expectEqual(0x1.00605cp+0, tan(@as(f32, 0xc.94p-4)));
    try std.testing.expectEqual(0x1.00204cp+0, tan(@as(f32, 0xc.92p-4)));
    try std.testing.expectEqual(0x1.00004ap+0, tan(@as(f32, 0xc.91p-4)));
    try std.testing.expectEqual(0x1.00000ap+0, tan(@as(f32, 0xc.90fep-4)));
    try std.testing.expectEqual(0x1.000002p+0, tan(@as(f32, 0xc.90fdcp-4)));
    try std.testing.expectEqual(0x1p+0, tan(@as(f32, 0xc.90fdbp-4)));
    try std.testing.expectEqual(-0xf.fe04dp-4, tan(@as(f32, -0xc.9p-4)));
    try std.testing.expectEqual(-0xf.ff04bp-4, tan(@as(f32, -0xc.908p-4)));
    try std.testing.expectEqual(-0xf.ff84bp-4, tan(@as(f32, -0xc.90cp-4)));
    try std.testing.expectEqual(-0xf.ffc4bp-4, tan(@as(f32, -0xc.90ep-4)));
    try std.testing.expectEqual(-0xf.ffe4bp-4, tan(@as(f32, -0xc.90fp-4)));
    try std.testing.expectEqual(-0xf.fff4bp-4, tan(@as(f32, -0xc.90f8p-4)));
    try std.testing.expectEqual(-0xf.fffcbp-4, tan(@as(f32, -0xc.90fcp-4)));
    try std.testing.expectEqual(-0xf.fffebp-4, tan(@as(f32, -0xc.90fdp-4)));
    try std.testing.expectEqual(-0xf.ffffbp-4, tan(@as(f32, -0xc.90fd8p-4)));
    try std.testing.expectEqual(-0xf.fffffp-4, tan(@as(f32, -0xc.90fdap-4)));
    try std.testing.expectEqual(-0x1.01e21p+0, tan(@as(f32, -0xc.ap-4)));
    try std.testing.expectEqual(-0x1.00e0aep+0, tan(@as(f32, -0xc.98p-4)));
    try std.testing.expectEqual(-0x1.00605cp+0, tan(@as(f32, -0xc.94p-4)));
    try std.testing.expectEqual(-0x1.00204cp+0, tan(@as(f32, -0xc.92p-4)));
    try std.testing.expectEqual(-0x1.00004ap+0, tan(@as(f32, -0xc.91p-4)));
    try std.testing.expectEqual(-0x1.00000ap+0, tan(@as(f32, -0xc.90fep-4)));
    try std.testing.expectEqual(-0x1.000002p+0, tan(@as(f32, -0xc.90fdcp-4)));
    try std.testing.expectEqual(-0x1p+0, tan(@as(f32, -0xc.90fdbp-4)));
    try std.testing.expectEqual(-0x5.08eea8p-4, tan(@as(f32, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(-0x1.14bdfcp+0, tan(@as(f32, 0x2.1e19ep+72)));
    try std.testing.expectEqual(-0x9.c9ecap-4, tan(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x9.c9ecap-4, tan(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.8eb246p+0, tan(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x2.2f5ec4p+0, tan(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(-0x2.47dee4p-4, tan(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(0x1.2866fap+0, tan(@as(f32, 0x4p+0)));
    try std.testing.expectEqual(-0x3.61697p+0, tan(@as(f32, 0x5p+0)));
    try std.testing.expectEqual(-0x4.a7f618p-4, tan(@as(f32, 0x6p+0)));
    try std.testing.expectEqual(0xd.f1737p-4, tan(@as(f32, 0x7p+0)));
    try std.testing.expectEqual(-0x6.ccb9ep+0, tan(@as(f32, 0x8p+0)));
    try std.testing.expectEqual(-0x7.3caf58p-4, tan(@as(f32, 0x9p+0)));
    try std.testing.expectEqual(0xa.5fafap-4, tan(@as(f32, 0xap+0)));
    try std.testing.expectEqual(-0x1.a4a482p+0, tan(@as(f32, -0x1.062a48p+0)));
    try std.testing.expectEqual(-0x3.c00d44p+0, tan(@as(f32, -0x1.4f69cp+0)));
    try std.testing.expectEqual(0x6.c89cf8p+0, tan(@as(f32, 0x1.6ca7e8p+0)));
    try std.testing.expectEqual(0x7.35553p+0, tan(@as(f32, -0x1.b569cp+0)));
    try std.testing.expectEqual(0x1.d1fa34p+0, tan(@as(f32, -0x2.12bafcp+0)));
    try std.testing.expectEqual(-0x1.fe847p+0, tan(@as(f32, 0x2.091d68p+0)));
    try std.testing.expectEqual(0x1.f0dbcep+0, tan(@as(f32, -0x5.302ab8p+0)));
    try std.testing.expectEqual(0x1.f0dba8p+0, tan(@as(f32, -0x5.302acp+0)));
    try std.testing.expectEqual(0x1.fcfe68p+0, tan(@as(f32, 0x1.1ad374p+0)));
    try std.testing.expectEqual(-0x1.c074f8p+0, tan(@as(f32, -0x1.0d55b8p+0)));
    try std.testing.expectEqual(-0x1.41acc2p+20, tan(@as(f32, 0x1.921fc2p+0)));
    try std.testing.expectEqual(-0x1.7d9d38p+20, tan(@as(f32, 0x1.921fcp+0)));
    try std.testing.expectEqual(0x1.7d9d38p+20, tan(@as(f32, -0x1.921fcp+0)));
    try std.testing.expectEqual(0x1.41acc2p+20, tan(@as(f32, -0x1.921fc2p+0)));
    try std.testing.expectEqual(0x8.00aacp-8, tan(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x4.000018p-12, tan(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x2p-16, tan(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0x1p-20, tan(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0x8p-28, tan(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x4p-32, tan(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, tan(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, tan(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, tan(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, tan(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, tan(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, tan(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, tan(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tan(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tan(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x9.c9ecap-4, tan(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x9.c9ecap-4, tan(@as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x4p-128, tan(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x4p-128, tan(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-152, tan(@as(f32, -0x8p-152)));

    try std.testing.expectEqual(0x0p+0, tan(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, tan(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1.000000bbbd2ecp+0, tan(@as(f64, 0xc.90fdbp-4)));
    try std.testing.expectEqual(0xf.ffffebbbd2f48p-4, tan(@as(f64, 0xc.90fdap-4)));
    try std.testing.expectEqual(0x1.0000000000001p+0, tan(@as(f64, 0xc.90fdaa22168c8p-4)));
    try std.testing.expectEqual(0xf.ffffffffffff8p-4, tan(@as(f64, 0xc.90fdaa22168cp-4)));
    try std.testing.expectEqual(-0x1.5d14946dc9897p+24, tan(@as(f64, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0xc.a1bd99b5b586p+20, tan(@as(f64, 0x1.921fb4p+0)));
    try std.testing.expectEqual(-0x1.617a15494767ap+52, tan(@as(f64, 0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(0x3.a052cf8639b6ap+52, tan(@as(f64, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(0x1.5d14946dc9897p+24, tan(@as(f64, -0x1.921fb6p+0)));
    try std.testing.expectEqual(-0xc.a1bd99b5b586p+20, tan(@as(f64, -0x1.921fb4p+0)));
    try std.testing.expectEqual(0x1.617a15494767ap+52, tan(@as(f64, -0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(-0x3.a052cf8639b6ap+52, tan(@as(f64, -0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(0xe.e7d1b0887776p-4, tan(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(-0xc.17b0bfdb2b808p-8, tan(@as(f64, 0x2p+64)));
    try std.testing.expectEqual(0xc.17b0bfdb2b808p-8, tan(@as(f64, -0x2p+64)));
    try std.testing.expectEqual(0x2p-28, tan(@as(f64, 0x2p-28)));
    try std.testing.expectEqual(-0x2p-28, tan(@as(f64, -0x2p-28)));
    try std.testing.expectEqual(0xf.fe04cb247202p-4, tan(@as(f64, 0xc.9p-4)));
    try std.testing.expectEqual(0xf.ff04b37174f7p-4, tan(@as(f64, 0xc.908p-4)));
    try std.testing.expectEqual(0xf.ff84ad971a078p-4, tan(@as(f64, 0xc.90cp-4)));
    try std.testing.expectEqual(0xf.ffc4ac29d171p-4, tan(@as(f64, 0xc.90ep-4)));
    try std.testing.expectEqual(0xf.ffe4abd329dep-4, tan(@as(f64, 0xc.90fp-4)));
    try std.testing.expectEqual(0xf.fff4abbfd5b28p-4, tan(@as(f64, 0xc.90f8p-4)));
    try std.testing.expectEqual(0xf.fffcabbc2b928p-4, tan(@as(f64, 0xc.90fcp-4)));
    try std.testing.expectEqual(0xf.fffeabbbe10ap-4, tan(@as(f64, 0xc.90fdp-4)));
    try std.testing.expectEqual(0xf.ffffabbbd3c58p-4, tan(@as(f64, 0xc.90fd8p-4)));
    try std.testing.expectEqual(0xf.ffffebbbd2f48p-4, tan(@as(f64, 0xc.90fdap-4)));
    try std.testing.expectEqual(0x1.01e20f7e06e4bp+0, tan(@as(f64, 0xc.ap-4)));
    try std.testing.expectEqual(0x1.00e0ad36afd0ep+0, tan(@as(f64, 0xc.98p-4)));
    try std.testing.expectEqual(0x1.00605cdc5a1a2p+0, tan(@as(f64, 0xc.94p-4)));
    try std.testing.expectEqual(0x1.00204cc54b6a7p+0, tan(@as(f64, 0xc.92p-4)));
    try std.testing.expectEqual(0x1.00004abbc817p+0, tan(@as(f64, 0xc.91p-4)));
    try std.testing.expectEqual(0x1.00000abbbd681p+0, tan(@as(f64, 0xc.90fep-4)));
    try std.testing.expectEqual(0x1.000002bbbd323p+0, tan(@as(f64, 0xc.90fdcp-4)));
    try std.testing.expectEqual(0x1.000000bbbd2ecp+0, tan(@as(f64, 0xc.90fdbp-4)));
    try std.testing.expectEqual(-0xf.fe04cb247202p-4, tan(@as(f64, -0xc.9p-4)));
    try std.testing.expectEqual(-0xf.ff04b37174f7p-4, tan(@as(f64, -0xc.908p-4)));
    try std.testing.expectEqual(-0xf.ff84ad971a078p-4, tan(@as(f64, -0xc.90cp-4)));
    try std.testing.expectEqual(-0xf.ffc4ac29d171p-4, tan(@as(f64, -0xc.90ep-4)));
    try std.testing.expectEqual(-0xf.ffe4abd329dep-4, tan(@as(f64, -0xc.90fp-4)));
    try std.testing.expectEqual(-0xf.fff4abbfd5b28p-4, tan(@as(f64, -0xc.90f8p-4)));
    try std.testing.expectEqual(-0xf.fffcabbc2b928p-4, tan(@as(f64, -0xc.90fcp-4)));
    try std.testing.expectEqual(-0xf.fffeabbbe10ap-4, tan(@as(f64, -0xc.90fdp-4)));
    try std.testing.expectEqual(-0xf.ffffabbbd3c58p-4, tan(@as(f64, -0xc.90fd8p-4)));
    try std.testing.expectEqual(-0xf.ffffebbbd2f48p-4, tan(@as(f64, -0xc.90fdap-4)));
    try std.testing.expectEqual(-0x1.01e20f7e06e4bp+0, tan(@as(f64, -0xc.ap-4)));
    try std.testing.expectEqual(-0x1.00e0ad36afd0ep+0, tan(@as(f64, -0xc.98p-4)));
    try std.testing.expectEqual(-0x1.00605cdc5a1a2p+0, tan(@as(f64, -0xc.94p-4)));
    try std.testing.expectEqual(-0x1.00204cc54b6a7p+0, tan(@as(f64, -0xc.92p-4)));
    try std.testing.expectEqual(-0x1.00004abbc817p+0, tan(@as(f64, -0xc.91p-4)));
    try std.testing.expectEqual(-0x1.00000abbbd681p+0, tan(@as(f64, -0xc.90fep-4)));
    try std.testing.expectEqual(-0x1.000002bbbd323p+0, tan(@as(f64, -0xc.90fdcp-4)));
    try std.testing.expectEqual(-0x1.000000bbbd2ecp+0, tan(@as(f64, -0xc.90fdbp-4)));
    try std.testing.expectEqual(-0x5.08eea5bdd992cp-4, tan(@as(f64, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(-0x1.14bdfb7ac8b93p+0, tan(@as(f64, 0x2.1e19ep+72)));
    try std.testing.expectEqual(-0x1.a0f79c1b6b257p+0, tan(@as(f64, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(-0x9.c9eca5a4c461p-4, tan(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xa.e735a60681518p-4, tan(@as(f64, 0x8p+1020)));
    try std.testing.expectEqual(-0x9.c9eca5a4c461p-4, tan(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.4530cfe729484p-8, tan(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.8eb245cbee3a6p+0, tan(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x2.2f5ec5c12a1fp+0, tan(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(-0x2.47dee24a970dep-4, tan(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(0x1.2866f9be4de13p+0, tan(@as(f64, 0x4p+0)));
    try std.testing.expectEqual(-0x3.61696e737760ep+0, tan(@as(f64, 0x5p+0)));
    try std.testing.expectEqual(-0x4.a7f61baee56f8p-4, tan(@as(f64, 0x6p+0)));
    try std.testing.expectEqual(0xd.f173709f753c8p-4, tan(@as(f64, 0x7p+0)));
    try std.testing.expectEqual(-0x6.ccb9e3d26879p+0, tan(@as(f64, 0x8p+0)));
    try std.testing.expectEqual(-0x7.3caf584c5707p-4, tan(@as(f64, 0x9p+0)));
    try std.testing.expectEqual(0xa.5faf9a5f1bc1p-4, tan(@as(f64, 0xap+0)));
    try std.testing.expectEqual(-0x1.a4a482f560f6ep+0, tan(@as(f64, -0x1.062a48p+0)));
    try std.testing.expectEqual(-0x3.c00d4280aa7bep+0, tan(@as(f64, -0x1.4f69cp+0)));
    try std.testing.expectEqual(0x6.c89cf9333573p+0, tan(@as(f64, 0x1.6ca7e8p+0)));
    try std.testing.expectEqual(0x7.35552c167cbe8p+0, tan(@as(f64, -0x1.b569cp+0)));
    try std.testing.expectEqual(0x1.d1fa3375a3decp+0, tan(@as(f64, -0x2.12bafcp+0)));
    try std.testing.expectEqual(-0x1.fe84705639d38p+0, tan(@as(f64, 0x2.091d68p+0)));
    try std.testing.expectEqual(0x1.f0dbcee9873ccp+0, tan(@as(f64, -0x5.302ab8p+0)));
    try std.testing.expectEqual(0x1.f0dba8c6e598cp+0, tan(@as(f64, -0x5.302acp+0)));
    try std.testing.expectEqual(0x1.f0dbc6d6f7b0ap+0, tan(@as(f64, -0x5.302ab9b18593p+0)));
    try std.testing.expectEqual(0x1.f0dbc6d6f7af7p+0, tan(@as(f64, -0x5.302ab9b185934p+0)));
    try std.testing.expectEqual(0x1.fcfe678552d4ap+0, tan(@as(f64, 0x1.1ad374p+0)));
    try std.testing.expectEqual(-0x1.c074f83e72237p+0, tan(@as(f64, -0x1.0d55b8p+0)));
    try std.testing.expectEqual(-0x1.41acc1f2aebdcp+20, tan(@as(f64, 0x1.921fc2p+0)));
    try std.testing.expectEqual(-0x1.7d9d370b14023p+20, tan(@as(f64, 0x1.921fcp+0)));
    try std.testing.expectEqual(-0x1.7b91a083509cfp+20, tan(@as(f64, 0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(-0x1.7b91a08583657p+20, tan(@as(f64, 0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(0x1.7d9d370b14023p+20, tan(@as(f64, -0x1.921fcp+0)));
    try std.testing.expectEqual(0x1.41acc1f2aebdcp+20, tan(@as(f64, -0x1.921fc2p+0)));
    try std.testing.expectEqual(0x1.7b91a08583657p+20, tan(@as(f64, -0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(0x1.7b91a083509cfp+20, tan(@as(f64, -0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(0x8.00aabbbd7604p-8, tan(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0x4.0000155555ddcp-12, tan(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0x2.00000002aaaaap-16, tan(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0x1.0000000000555p-20, tan(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0x8.0000000000008p-28, tan(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x4p-32, tan(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, tan(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, tan(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, tan(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, tan(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, tan(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, tan(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, tan(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tan(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, tan(@as(f64, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tan(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, tan(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x9.c9eca5a4c461p-4, tan(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.4530cfe729484p-8, tan(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x9.c9eca5a4c461p-4, tan(@as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.4530cfe729484p-8, tan(@as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4p-128, tan(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, tan(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-972, tan(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, tan(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, tan(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x8p-972, tan(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, tan(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-152, tan(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, tan(@as(f64, -0x4p-1076)));

    try std.testing.expectEqual(0x0p+0, tan(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, tan(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1.000000bbbd2ec06ep+0, tan(@as(f80, 0xc.90fdbp-4)));
    try std.testing.expectEqual(0xf.ffffebbbd2f48f3p-4, tan(@as(f80, 0xc.90fdap-4)));
    try std.testing.expectEqual(0x1.0000000000000b96p+0, tan(@as(f80, 0xc.90fdaa22168c8p-4)));
    try std.testing.expectEqual(0xf.ffffffffffffb96p-4, tan(@as(f80, 0xc.90fdaa22168cp-4)));
    try std.testing.expectEqual(0x1p+0, tan(@as(f80, 0xc.90fdaa22168c235p-4)));
    try std.testing.expectEqual(0xf.ffffffffffffffep-4, tan(@as(f80, 0xc.90fdaa22168c234p-4)));
    try std.testing.expectEqual(-0x1.5d14946dc98975d6p+24, tan(@as(f80, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0xc.a1bd99b5b58623dp+20, tan(@as(f80, 0x1.921fb4p+0)));
    try std.testing.expectEqual(-0x1.617a15494767a048p+52, tan(@as(f80, 0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(0x3.a052cf8639b69c18p+52, tan(@as(f80, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(-0x2.29478136aaf68d7cp+64, tan(@as(f80, 0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(0xa.686780675d73f75p+60, tan(@as(f80, 0x1.921fb54442d18468p+0)));
    try std.testing.expectEqual(0x1.5d14946dc98975d6p+24, tan(@as(f80, -0x1.921fb6p+0)));
    try std.testing.expectEqual(-0xc.a1bd99b5b58623dp+20, tan(@as(f80, -0x1.921fb4p+0)));
    try std.testing.expectEqual(0x1.617a15494767a048p+52, tan(@as(f80, -0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(-0x3.a052cf8639b69c18p+52, tan(@as(f80, -0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(0x2.29478136aaf68d7cp+64, tan(@as(f80, -0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(-0xa.686780675d73f75p+60, tan(@as(f80, -0x1.921fb54442d18468p+0)));
    try std.testing.expectEqual(0xe.e7d1b0887775f06p-4, tan(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(-0xc.17b0bfdb2b8061ep-8, tan(@as(f80, 0x2p+64)));
    try std.testing.expectEqual(0xc.17b0bfdb2b8061ep-8, tan(@as(f80, -0x2p+64)));
    try std.testing.expectEqual(0x2.00000000000002acp-28, tan(@as(f80, 0x2p-28)));
    try std.testing.expectEqual(-0x2.00000000000002acp-28, tan(@as(f80, -0x2p-28)));
    try std.testing.expectEqual(0xf.fe04cb2472021f2p-4, tan(@as(f80, 0xc.9p-4)));
    try std.testing.expectEqual(0xf.ff04b37174f6f35p-4, tan(@as(f80, 0xc.908p-4)));
    try std.testing.expectEqual(0xf.ff84ad971a07664p-4, tan(@as(f80, 0xc.90cp-4)));
    try std.testing.expectEqual(0xf.ffc4ac29d1711ccp-4, tan(@as(f80, 0xc.90ep-4)));
    try std.testing.expectEqual(0xf.ffe4abd329de183p-4, tan(@as(f80, 0xc.90fp-4)));
    try std.testing.expectEqual(0xf.fff4abbfd5b29a3p-4, tan(@as(f80, 0xc.90f8p-4)));
    try std.testing.expectEqual(0xf.fffcabbc2b925c1p-4, tan(@as(f80, 0xc.90fcp-4)));
    try std.testing.expectEqual(0xf.fffeabbbe109e1fp-4, tan(@as(f80, 0xc.90fdp-4)));
    try std.testing.expectEqual(0xf.ffffabbbd3c59fep-4, tan(@as(f80, 0xc.90fd8p-4)));
    try std.testing.expectEqual(0xf.ffffebbbd2f48f3p-4, tan(@as(f80, 0xc.90fdap-4)));
    try std.testing.expectEqual(0x1.01e20f7e06e4b406p+0, tan(@as(f80, 0xc.ap-4)));
    try std.testing.expectEqual(0x1.00e0ad36afd0da36p+0, tan(@as(f80, 0xc.98p-4)));
    try std.testing.expectEqual(0x1.00605cdc5a1a1c1ep+0, tan(@as(f80, 0xc.94p-4)));
    try std.testing.expectEqual(0x1.00204cc54b6a75fcp+0, tan(@as(f80, 0xc.92p-4)));
    try std.testing.expectEqual(0x1.00004abbc817073cp+0, tan(@as(f80, 0xc.91p-4)));
    try std.testing.expectEqual(0x1.00000abbbd6815d2p+0, tan(@as(f80, 0xc.90fep-4)));
    try std.testing.expectEqual(0x1.000002bbbd3237e8p+0, tan(@as(f80, 0xc.90fdcp-4)));
    try std.testing.expectEqual(0x1.000000bbbd2ec06ep+0, tan(@as(f80, 0xc.90fdbp-4)));
    try std.testing.expectEqual(-0xf.fe04cb2472021f2p-4, tan(@as(f80, -0xc.9p-4)));
    try std.testing.expectEqual(-0xf.ff04b37174f6f35p-4, tan(@as(f80, -0xc.908p-4)));
    try std.testing.expectEqual(-0xf.ff84ad971a07664p-4, tan(@as(f80, -0xc.90cp-4)));
    try std.testing.expectEqual(-0xf.ffc4ac29d1711ccp-4, tan(@as(f80, -0xc.90ep-4)));
    try std.testing.expectEqual(-0xf.ffe4abd329de183p-4, tan(@as(f80, -0xc.90fp-4)));
    try std.testing.expectEqual(-0xf.fff4abbfd5b29a3p-4, tan(@as(f80, -0xc.90f8p-4)));
    try std.testing.expectEqual(-0xf.fffcabbc2b925c1p-4, tan(@as(f80, -0xc.90fcp-4)));
    try std.testing.expectEqual(-0xf.fffeabbbe109e1fp-4, tan(@as(f80, -0xc.90fdp-4)));
    try std.testing.expectEqual(-0xf.ffffabbbd3c59fep-4, tan(@as(f80, -0xc.90fd8p-4)));
    try std.testing.expectEqual(-0xf.ffffebbbd2f48f3p-4, tan(@as(f80, -0xc.90fdap-4)));
    try std.testing.expectEqual(-0x1.01e20f7e06e4b406p+0, tan(@as(f80, -0xc.ap-4)));
    try std.testing.expectEqual(-0x1.00e0ad36afd0da36p+0, tan(@as(f80, -0xc.98p-4)));
    try std.testing.expectEqual(-0x1.00605cdc5a1a1c1ep+0, tan(@as(f80, -0xc.94p-4)));
    try std.testing.expectEqual(-0x1.00204cc54b6a75fcp+0, tan(@as(f80, -0xc.92p-4)));
    try std.testing.expectEqual(-0x1.00004abbc817073cp+0, tan(@as(f80, -0xc.91p-4)));
    try std.testing.expectEqual(-0x1.00000abbbd6815d2p+0, tan(@as(f80, -0xc.90fep-4)));
    try std.testing.expectEqual(-0x1.000002bbbd3237e8p+0, tan(@as(f80, -0xc.90fdcp-4)));
    try std.testing.expectEqual(-0x1.000000bbbd2ec06ep+0, tan(@as(f80, -0xc.90fdbp-4)));
    try std.testing.expectEqual(-0x5.08eea5bdd992a0dp-4, tan(@as(f80, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(-0x1.14bdfb7ac8b928bap+0, tan(@as(f80, 0x2.1e19ep+72)));
    try std.testing.expectEqual(-0x1.a0f79c1b6b25774ap+0, tan(@as(f80, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(-0x9.c9eca5a4c460f93p-4, tan(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xa.e735a6068151a9ep-4, tan(@as(f80, 0x8p+1020)));
    try std.testing.expectEqual(-0x9.c9eca5a4c460f93p-4, tan(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.4530cfe729483b8ep-8, tan(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x6.c3788e85da9be4f8p-4, tan(@as(f80, 0x8p+16380)));
    try std.testing.expectEqual(0x1.8eb245cbee3a5b8ap+0, tan(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x2.2f5ec5c12a1f076cp+0, tan(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(-0x2.47dee24a970de198p-4, tan(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(0x1.2866f9be4de1370ep+0, tan(@as(f80, 0x4p+0)));
    try std.testing.expectEqual(-0x3.61696e737760d084p+0, tan(@as(f80, 0x5p+0)));
    try std.testing.expectEqual(-0x4.a7f61baee56f8c1p-4, tan(@as(f80, 0x6p+0)));
    try std.testing.expectEqual(0xd.f173709f753c4c1p-4, tan(@as(f80, 0x7p+0)));
    try std.testing.expectEqual(-0x6.ccb9e3d26878e9c8p+0, tan(@as(f80, 0x8p+0)));
    try std.testing.expectEqual(-0x7.3caf584c5706f808p-4, tan(@as(f80, 0x9p+0)));
    try std.testing.expectEqual(0xa.5faf9a5f1bc12fp-4, tan(@as(f80, 0xap+0)));
    try std.testing.expectEqual(-0x1.a4a482f560f6e4cep+0, tan(@as(f80, -0x1.062a48p+0)));
    try std.testing.expectEqual(-0x3.c00d4280aa7bede8p+0, tan(@as(f80, -0x1.4f69cp+0)));
    try std.testing.expectEqual(0x6.c89cf93335731da8p+0, tan(@as(f80, 0x1.6ca7e8p+0)));
    try std.testing.expectEqual(0x7.35552c167cbe769p+0, tan(@as(f80, -0x1.b569cp+0)));
    try std.testing.expectEqual(0x1.d1fa3375a3dec7e8p+0, tan(@as(f80, -0x2.12bafcp+0)));
    try std.testing.expectEqual(-0x1.fe84705639d38772p+0, tan(@as(f80, 0x2.091d68p+0)));
    try std.testing.expectEqual(0x1.f0dbcee9873cbf2ep+0, tan(@as(f80, -0x5.302ab8p+0)));
    try std.testing.expectEqual(0x1.f0dba8c6e598b932p+0, tan(@as(f80, -0x5.302acp+0)));
    try std.testing.expectEqual(0x1.f0dbc6d6f7b0a6e6p+0, tan(@as(f80, -0x5.302ab9b18593p+0)));
    try std.testing.expectEqual(0x1.f0dbc6d6f7af75dp+0, tan(@as(f80, -0x5.302ab9b185934p+0)));
    try std.testing.expectEqual(0x1.f0dbc6d6f7aff09p+0, tan(@as(f80, -0x5.302ab9b18593264p+0)));
    try std.testing.expectEqual(0x1.fcfe678552d4a7ccp+0, tan(@as(f80, 0x1.1ad374p+0)));
    try std.testing.expectEqual(-0x1.c074f83e72236f1ap+0, tan(@as(f80, -0x1.0d55b8p+0)));
    try std.testing.expectEqual(-0x1.41acc1f2aebdbbdp+20, tan(@as(f80, 0x1.921fc2p+0)));
    try std.testing.expectEqual(-0x1.7d9d370b140234d4p+20, tan(@as(f80, 0x1.921fcp+0)));
    try std.testing.expectEqual(-0x1.7b91a083509ce99ap+20, tan(@as(f80, 0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(-0x1.7b91a08583656db8p+20, tan(@as(f80, 0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(-0x1.7b91a0851b85eb56p+20, tan(@as(f80, 0x1.921fc00ece4f02f4p+0)));
    try std.testing.expectEqual(-0x1.7b91a0851bcc4466p+20, tan(@as(f80, 0x1.921fc00ece4f02f2p+0)));
    try std.testing.expectEqual(0x1.7d9d370b140234d4p+20, tan(@as(f80, -0x1.921fcp+0)));
    try std.testing.expectEqual(0x1.41acc1f2aebdbbdp+20, tan(@as(f80, -0x1.921fc2p+0)));
    try std.testing.expectEqual(0x1.7b91a08583656db8p+20, tan(@as(f80, -0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(0x1.7b91a083509ce99ap+20, tan(@as(f80, -0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(0x1.7b91a0851bcc4466p+20, tan(@as(f80, -0x1.921fc00ece4f02f2p+0)));
    try std.testing.expectEqual(0x1.7b91a0851b85eb56p+20, tan(@as(f80, -0x1.921fc00ece4f02f4p+0)));
    try std.testing.expectEqual(0x8.00aabbbd76042bfp-8, tan(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x4.0000155555ddddep-12, tan(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x2.00000002aaaaaabp-16, tan(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0x1.0000000000555556p-20, tan(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0x8.000000000000aabp-28, tan(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x4.0000000000000018p-32, tan(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x2p-36, tan(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0x1p-40, tan(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x8p-48, tan(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0x4p-52, tan(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x2p-56, tan(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, tan(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, tan(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tan(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, tan(@as(f80, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tan(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, tan(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, tan(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(-0x9.c9eca5a4c460f93p-4, tan(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.4530cfe729483b8ep-8, tan(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x7.ef32a4ca67437ed8p+0, tan(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x9.c9eca5a4c460f93p-4, tan(@as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.4530cfe729483b8ep-8, tan(@as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x7.ef32a4ca67437ed8p+0, tan(@as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4p-128, tan(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, tan(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, tan(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, tan(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, tan(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, tan(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, tan(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, tan(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, tan(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, tan(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, tan(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, tan(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-152, tan(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, tan(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, tan(@as(f80, -0x8p-16448)));

    try std.testing.expectEqual(0x0p+0, tan(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, tan(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1.000000bbbd2ec06d6d6fff3655a3p+0, tan(@as(f128, 0xc.90fdbp-4)));
    try std.testing.expectEqual(0xf.ffffebbbd2f48f30fa9c07dc0468p-4, tan(@as(f128, 0xc.90fdap-4)));
    try std.testing.expectEqual(0x1.0000000000000b9676733ae8fe8bp+0, tan(@as(f128, 0xc.90fdaa22168c8p-4)));
    try std.testing.expectEqual(0xf.ffffffffffffb9676733ae8fe518p-4, tan(@as(f128, 0xc.90fdaa22168cp-4)));
    try std.testing.expectEqual(0x1.000000000000000076733ae8fe48p+0, tan(@as(f128, 0xc.90fdaa22168c235p-4)));
    try std.testing.expectEqual(0xf.ffffffffffffffe76733ae8fe48p-4, tan(@as(f128, 0xc.90fdaa22168c234p-4)));
    try std.testing.expectEqual(0x1.0000000000000000000000000001p+0, tan(@as(f128, 0xc.90fdaa22168c234c4c6628b80dc8p-4)));
    try std.testing.expectEqual(0x1p+0, tan(@as(f128, 0xc.90fdaa22168c234c4c6628b80dcp-4)));
    try std.testing.expectEqual(0x1.0000000000000000000000000048p+0, tan(@as(f128, 0xc.90fdaa22168c234c4c6628b81p-4)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffc8p-4, tan(@as(f128, 0xc.90fdaa22168c234c4c6628b80cp-4)));
    try std.testing.expectEqual(-0x1.5d14946dc98975d6421a55284fep+24, tan(@as(f128, 0x1.921fb6p+0)));
    try std.testing.expectEqual(0xc.a1bd99b5b58623cd91404ccd8ca8p+20, tan(@as(f128, 0x1.921fb4p+0)));
    try std.testing.expectEqual(-0x1.617a15494767a04882c320317f3ep+52, tan(@as(f128, 0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(0x3.a052cf8639b69c1871a036cababcp+52, tan(@as(f128, 0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(-0x2.29478136aaf68d7b3b807fb349bap+64, tan(@as(f128, 0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(0xa.686780675d73f74c339c44a53238p+60, tan(@as(f128, 0x1.921fb54442d18468p+0)));
    // try std.testing.expectEqual(-0x1.4a611a1bc0c6c27a11d14bf714f1p+112, tan(@as(f128, 0x1.921fb54442d18469898cc51701b9p+0)));
    try std.testing.expectEqual(0x4.711af55c9de64e8bb98064f255ecp+112, tan(@as(f128, 0x1.921fb54442d18469898cc51701b8p+0)));
    try std.testing.expectEqual(-0x3.9113c85ed0a9399bff03c9d9d62ap+104, tan(@as(f128, 0x1.921fb54442d18469898cc51702p+0)));
    try std.testing.expectEqual(0x4.8d99880090a89825c6aaa8a4c2dcp+104, tan(@as(f128, 0x1.921fb54442d18469898cc517018p+0)));
    try std.testing.expectEqual(0x1.5d14946dc98975d6421a55284fep+24, tan(@as(f128, -0x1.921fb6p+0)));
    try std.testing.expectEqual(-0xc.a1bd99b5b58623cd91404ccd8ca8p+20, tan(@as(f128, -0x1.921fb4p+0)));
    try std.testing.expectEqual(0x1.617a15494767a04882c320317f3ep+52, tan(@as(f128, -0x1.921fb54442d19p+0)));
    try std.testing.expectEqual(-0x3.a052cf8639b69c1871a036cababcp+52, tan(@as(f128, -0x1.921fb54442d18p+0)));
    try std.testing.expectEqual(0x2.29478136aaf68d7b3b807fb349bap+64, tan(@as(f128, -0x1.921fb54442d1846ap+0)));
    try std.testing.expectEqual(-0xa.686780675d73f74c339c44a53238p+60, tan(@as(f128, -0x1.921fb54442d18468p+0)));
    // try std.testing.expectEqual(0x1.4a611a1bc0c6c27a11d14bf714f1p+112, tan(@as(f128, -0x1.921fb54442d18469898cc51701b9p+0)));
    try std.testing.expectEqual(-0x4.711af55c9de64e8bb98064f255ecp+112, tan(@as(f128, -0x1.921fb54442d18469898cc51701b8p+0)));
    try std.testing.expectEqual(0x3.9113c85ed0a9399bff03c9d9d62ap+104, tan(@as(f128, -0x1.921fb54442d18469898cc51702p+0)));
    try std.testing.expectEqual(-0x4.8d99880090a89825c6aaa8a4c2dcp+104, tan(@as(f128, -0x1.921fb54442d18469898cc517018p+0)));
    try std.testing.expectEqual(0xe.e7d1b0887775f06184cd76c016fp-4, tan(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(-0xc.17b0bfdb2b8061e7b11d50087308p-8, tan(@as(f128, 0x2p+64)));
    try std.testing.expectEqual(0xc.17b0bfdb2b8061e7b11d50087308p-8, tan(@as(f128, -0x2p+64)));
    try std.testing.expectEqual(0x2.00000000000002aaaaaaaaaaaaaep-28, tan(@as(f128, 0x2p-28)));
    try std.testing.expectEqual(-0x2.00000000000002aaaaaaaaaaaaaep-28, tan(@as(f128, -0x2p-28)));
    try std.testing.expectEqual(0xf.fe04cb2472021f1945ff358e1fa8p-4, tan(@as(f128, 0xc.9p-4)));
    try std.testing.expectEqual(0xf.ff04b37174f6f3528b2b1b69f988p-4, tan(@as(f128, 0xc.908p-4)));
    try std.testing.expectEqual(0xf.ff84ad971a07663e3871b57341ep-4, tan(@as(f128, 0xc.90cp-4)));
    try std.testing.expectEqual(0xf.ffc4ac29d1711cbfd2ecbfc18d6p-4, tan(@as(f128, 0xc.90ep-4)));
    try std.testing.expectEqual(0xf.ffe4abd329de1834a397d6e46a5p-4, tan(@as(f128, 0xc.90fp-4)));
    try std.testing.expectEqual(0xf.fff4abbfd5b29a33e190809aa1cp-4, tan(@as(f128, 0xc.90f8p-4)));
    try std.testing.expectEqual(0xf.fffcabbc2b925c0d4e40f079838p-4, tan(@as(f128, 0xc.90fcp-4)));
    try std.testing.expectEqual(0xf.fffeabbbe109e1ee8991b9d3ba38p-4, tan(@as(f128, 0xc.90fdp-4)));
    try std.testing.expectEqual(0xf.ffffabbbd3c59fe25b70006f7cbp-4, tan(@as(f128, 0xc.90fd8p-4)));
    try std.testing.expectEqual(0xf.ffffebbbd2f48f30fa9c07dc0468p-4, tan(@as(f128, 0xc.90fdap-4)));
    try std.testing.expectEqual(0x1.01e20f7e06e4b4069f6fd6809128p+0, tan(@as(f128, 0xc.ap-4)));
    try std.testing.expectEqual(0x1.00e0ad36afd0da359300dc8485abp+0, tan(@as(f128, 0xc.98p-4)));
    try std.testing.expectEqual(0x1.00605cdc5a1a1c1e7a2e9db9f98ap+0, tan(@as(f128, 0xc.94p-4)));
    try std.testing.expectEqual(0x1.00204cc54b6a75fbaa11ed8cf1c4p+0, tan(@as(f128, 0xc.92p-4)));
    try std.testing.expectEqual(0x1.00004abbc817073c57de4e2c7227p+0, tan(@as(f128, 0xc.91p-4)));
    try std.testing.expectEqual(0x1.00000abbbd6815d2da4ff16a8a5fp+0, tan(@as(f128, 0xc.90fep-4)));
    try std.testing.expectEqual(0x1.000002bbbd3237e7d114276ed32ap+0, tan(@as(f128, 0xc.90fdcp-4)));
    try std.testing.expectEqual(0x1.000000bbbd2ec06d6d6fff3655a3p+0, tan(@as(f128, 0xc.90fdbp-4)));
    try std.testing.expectEqual(-0xf.fe04cb2472021f1945ff358e1fa8p-4, tan(@as(f128, -0xc.9p-4)));
    try std.testing.expectEqual(-0xf.ff04b37174f6f3528b2b1b69f988p-4, tan(@as(f128, -0xc.908p-4)));
    try std.testing.expectEqual(-0xf.ff84ad971a07663e3871b57341ep-4, tan(@as(f128, -0xc.90cp-4)));
    try std.testing.expectEqual(-0xf.ffc4ac29d1711cbfd2ecbfc18d6p-4, tan(@as(f128, -0xc.90ep-4)));
    try std.testing.expectEqual(-0xf.ffe4abd329de1834a397d6e46a5p-4, tan(@as(f128, -0xc.90fp-4)));
    try std.testing.expectEqual(-0xf.fff4abbfd5b29a33e190809aa1cp-4, tan(@as(f128, -0xc.90f8p-4)));
    try std.testing.expectEqual(-0xf.fffcabbc2b925c0d4e40f079838p-4, tan(@as(f128, -0xc.90fcp-4)));
    try std.testing.expectEqual(-0xf.fffeabbbe109e1ee8991b9d3ba38p-4, tan(@as(f128, -0xc.90fdp-4)));
    try std.testing.expectEqual(-0xf.ffffabbbd3c59fe25b70006f7cbp-4, tan(@as(f128, -0xc.90fd8p-4)));
    try std.testing.expectEqual(-0xf.ffffebbbd2f48f30fa9c07dc0468p-4, tan(@as(f128, -0xc.90fdap-4)));
    try std.testing.expectEqual(-0x1.01e20f7e06e4b4069f6fd6809128p+0, tan(@as(f128, -0xc.ap-4)));
    try std.testing.expectEqual(-0x1.00e0ad36afd0da359300dc8485abp+0, tan(@as(f128, -0xc.98p-4)));
    try std.testing.expectEqual(-0x1.00605cdc5a1a1c1e7a2e9db9f98ap+0, tan(@as(f128, -0xc.94p-4)));
    try std.testing.expectEqual(-0x1.00204cc54b6a75fbaa11ed8cf1c4p+0, tan(@as(f128, -0xc.92p-4)));
    try std.testing.expectEqual(-0x1.00004abbc817073c57de4e2c7227p+0, tan(@as(f128, -0xc.91p-4)));
    try std.testing.expectEqual(-0x1.00000abbbd6815d2da4ff16a8a5fp+0, tan(@as(f128, -0xc.90fep-4)));
    try std.testing.expectEqual(-0x1.000002bbbd3237e7d114276ed32ap+0, tan(@as(f128, -0xc.90fdcp-4)));
    try std.testing.expectEqual(-0x1.000000bbbd2ec06d6d6fff3655a3p+0, tan(@as(f128, -0xc.90fdbp-4)));
    try std.testing.expectEqual(-0x5.08eea5bdd992a0d19c9356cc5168p-4, tan(@as(f128, 0x2.1e19e4p+72)));
    try std.testing.expectEqual(-0x1.14bdfb7ac8b928ba2c1adb3ceb44p+0, tan(@as(f128, 0x2.1e19ep+72)));
    try std.testing.expectEqual(-0x1.a0f79c1b6b257749e043d5cdf75p+0, tan(@as(f128, 0x2.1e19e0c9bab24p+72)));
    try std.testing.expectEqual(-0x9.c9eca5a4c460f92a1a2e4fbecf5p-4, tan(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xa.e735a6068151a9df841fc42ab17p-4, tan(@as(f128, 0x8p+1020)));
    try std.testing.expectEqual(-0x9.c9eca5a4c460f92a1a2e4fbecf5p-4, tan(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.4530cfe729483b8da1f7101e16cdp-8, tan(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x6.c3788e85da9be4fb78d39ebd3f3p-4, tan(@as(f128, 0x8p+16380)));
    try std.testing.expectEqual(0x2.9d36f38857f642f5fdd53dc00078p+0, tan(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.8eb245cbee3a5b8acc7d41323141p+0, tan(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x2.2f5ec5c12a1f076a13210f55ecf6p+0, tan(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(-0x2.47dee24a970de1996164fbff0a86p-4, tan(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(0x1.2866f9be4de1370db9078607012dp+0, tan(@as(f128, 0x4p+0)));
    try std.testing.expectEqual(-0x3.61696e737760d0851798bac59114p+0, tan(@as(f128, 0x5p+0)));
    try std.testing.expectEqual(-0x4.a7f61baee56f8c0d5cb480072ef8p-4, tan(@as(f128, 0x6p+0)));
    try std.testing.expectEqual(0xd.f173709f753c4c117c5feb1485a8p-4, tan(@as(f128, 0x7p+0)));
    try std.testing.expectEqual(-0x6.ccb9e3d26878e9c70c0fe7c54824p+0, tan(@as(f128, 0x8p+0)));
    try std.testing.expectEqual(-0x7.3caf584c5706f80670ce6ab1353cp-4, tan(@as(f128, 0x9p+0)));
    try std.testing.expectEqual(0xa.5faf9a5f1bc12efead12fa488fdp-4, tan(@as(f128, 0xap+0)));
    try std.testing.expectEqual(-0x1.a4a482f560f6e4ceb9d6e73567d4p+0, tan(@as(f128, -0x1.062a48p+0)));
    try std.testing.expectEqual(-0x3.c00d4280aa7bede62d35d88620cp+0, tan(@as(f128, -0x1.4f69cp+0)));
    try std.testing.expectEqual(0x6.c89cf93335731da4f3da65c3026p+0, tan(@as(f128, 0x1.6ca7e8p+0)));
    try std.testing.expectEqual(0x7.35552c167cbe768ef1a28179e914p+0, tan(@as(f128, -0x1.b569cp+0)));
    try std.testing.expectEqual(0x1.d1fa3375a3dec7e8e0c696c99bd8p+0, tan(@as(f128, -0x2.12bafcp+0)));
    try std.testing.expectEqual(-0x1.fe84705639d387710dae52455d5p+0, tan(@as(f128, 0x2.091d68p+0)));
    try std.testing.expectEqual(0x1.f0dbcee9873cbf2ee7067b00e8p+0, tan(@as(f128, -0x5.302ab8p+0)));
    try std.testing.expectEqual(0x1.f0dba8c6e598b93213c2ad3e1476p+0, tan(@as(f128, -0x5.302acp+0)));
    try std.testing.expectEqual(0x1.f0dbc6d6f7b0a6e58732bf53fb8p+0, tan(@as(f128, -0x5.302ab9b18593p+0)));
    try std.testing.expectEqual(0x1.f0dbc6d6f7af75d06f670bd22cep+0, tan(@as(f128, -0x5.302ab9b185934p+0)));
    try std.testing.expectEqual(0x1.f0dbc6d6f7aff08febfa010b436ap+0, tan(@as(f128, -0x5.302ab9b18593264p+0)));
    try std.testing.expectEqual(0x1.fcfe678552d4a7cc4536fb771832p+0, tan(@as(f128, 0x1.1ad374p+0)));
    // try std.testing.expectEqual(-0x1.c074f83e72236f1900dbba65f783p+0, tan(@as(f128, -0x1.0d55b8p+0)));
    try std.testing.expectEqual(-0x1.41acc1f2aebdbbd0c10fd22cc56ap+20, tan(@as(f128, 0x1.921fc2p+0)));
    try std.testing.expectEqual(-0x1.7d9d370b140234d37501f83845c2p+20, tan(@as(f128, 0x1.921fcp+0)));
    try std.testing.expectEqual(-0x1.7b91a083509ce99968d2431cf3a4p+20, tan(@as(f128, 0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(-0x1.7b91a08583656db8b65798b4bbabp+20, tan(@as(f128, 0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(-0x1.7b91a0851b85eb5571b2cd0ca16fp+20, tan(@as(f128, 0x1.921fc00ece4f02f4p+0)));
    try std.testing.expectEqual(-0x1.7b91a0851bcc4465f5de57e1e183p+20, tan(@as(f128, 0x1.921fc00ece4f02f2p+0)));
    try std.testing.expectEqual(-0x1.7b91a0851bbbafa14cf22c61e331p+20, tan(@as(f128, 0x1.921fc00ece4f02f278ade6ad9e8ap+0)));
    try std.testing.expectEqual(-0x1.7b91a0851bbbafa14cf22c850fb9p+20, tan(@as(f128, 0x1.921fc00ece4f02f278ade6ad9e89p+0)));
    try std.testing.expectEqual(-0x1.7b91a0851bbbafa14cf21c2b5c63p+20, tan(@as(f128, 0x1.921fc00ece4f02f278ade6ad9fp+0)));
    try std.testing.expectEqual(-0x1.7b91a0851bbbafa14cf22dc1a084p+20, tan(@as(f128, 0x1.921fc00ece4f02f278ade6ad9e8p+0)));
    try std.testing.expectEqual(0x1.7d9d370b140234d37501f83845c2p+20, tan(@as(f128, -0x1.921fcp+0)));
    try std.testing.expectEqual(0x1.41acc1f2aebdbbd0c10fd22cc56ap+20, tan(@as(f128, -0x1.921fc2p+0)));
    try std.testing.expectEqual(0x1.7b91a08583656db8b65798b4bbabp+20, tan(@as(f128, -0x1.921fc00ece4fp+0)));
    try std.testing.expectEqual(0x1.7b91a083509ce99968d2431cf3a4p+20, tan(@as(f128, -0x1.921fc00ece4f1p+0)));
    try std.testing.expectEqual(0x1.7b91a0851bcc4465f5de57e1e183p+20, tan(@as(f128, -0x1.921fc00ece4f02f2p+0)));
    try std.testing.expectEqual(0x1.7b91a0851b85eb5571b2cd0ca16fp+20, tan(@as(f128, -0x1.921fc00ece4f02f4p+0)));
    try std.testing.expectEqual(0x1.7b91a0851bbbafa14cf22c850fb9p+20, tan(@as(f128, -0x1.921fc00ece4f02f278ade6ad9e89p+0)));
    try std.testing.expectEqual(0x1.7b91a0851bbbafa14cf22c61e331p+20, tan(@as(f128, -0x1.921fc00ece4f02f278ade6ad9e8ap+0)));
    try std.testing.expectEqual(0x1.7b91a0851bbbafa14cf22dc1a084p+20, tan(@as(f128, -0x1.921fc00ece4f02f278ade6ad9e8p+0)));
    try std.testing.expectEqual(0x1.7b91a0851bbbafa14cf21c2b5c63p+20, tan(@as(f128, -0x1.921fc00ece4f02f278ade6ad9fp+0)));
    try std.testing.expectEqual(0x8.00aabbbd76042be9164bf404c3p-8, tan(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(0x4.0000155555dddde1521537b70a3cp-12, tan(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x2.00000002aaaaaaaeeeeeeef5d75ep-16, tan(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0x1.0000000000555555555577777777p-20, tan(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0x8.000000000000aaaaaaaaaaaabbb8p-28, tan(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x4.0000000000000015555555555554p-32, tan(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0x2.000000000000000002aaaaaaaaaap-36, tan(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0x1.0000000000000000000055555555p-40, tan(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0x8.0000000000000000000000aaaaa8p-48, tan(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0x4.0000000000000000000000001554p-52, tan(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x2.0000000000000000000000000002p-56, tan(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1p-60, tan(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0x1p-100, tan(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tan(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1p-600, tan(@as(f128, 0x1p-600)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, tan(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x4p-1076, tan(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-10000, tan(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(-0x9.c9eca5a4c460f92a1a2e4fbecf5p-4, tan(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.4530cfe729483b8da1f7101e16cdp-8, tan(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x7.ef32a4ca67437ed7bee844b6959cp+0, tan(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.1b6e2c58e228a81d9500ff5ce722p+0, tan(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.9d36f38857f642f5fdd53dc00078p+0, tan(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x9.c9eca5a4c460f92a1a2e4fbecf5p-4, tan(@as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.4530cfe729483b8da1f7101e16cdp-8, tan(@as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x7.ef32a4ca67437ed7bee844b6959cp+0, tan(@as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.1b6e2c58e228a81d9500ff5ce722p+0, tan(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x2.9d36f38857f642f5fdd53dc00078p+0, tan(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4p-128, tan(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x4p-1024, tan(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x4p-16384, tan(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-16384, tan(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-972, tan(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x4p-128, tan(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x4p-1024, tan(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x4p-16384, tan(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2p-16384, tan(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x8p-972, tan(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x8p-152, tan(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x4p-1076, tan(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-16448, tan(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16448, tan(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x4p-16496, tan(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-152, tan(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x4p-1076, tan(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x8p-16448, tan(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16448, tan(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x4p-16496, tan(@as(f128, -0x4p-16496)));
}
