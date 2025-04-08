const cast = @import("../types.zig").cast;
const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const dla = @import("dla.zig");
const atnat2 = @import("atnat2.zig");
const ldbl128 = @import("ldbl128.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;

pub fn atan2(y: anytype, x: anytype) EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))) {
    comptime if (!types.isFixedPrecision(@TypeOf(y)) or types.isComplex(@TypeOf(y)))
        @compileError("y must be an int or float");

    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(y))) {
        .int => {
            switch (types.numericType(@TypeOf(x))) {
                .int => return atan2(cast(EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))), y, .{}), cast(EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))), x, .{})),
                .float => return atan2(cast(Coerce(@TypeOf(y), @TypeOf(x)), y, .{}), cast(Coerce(@TypeOf(y), @TypeOf(x)), x, .{})),
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(x))) {
                .int => return atan2(cast(Coerce(@TypeOf(y), @TypeOf(x)), y, .{}), cast(Coerce(@TypeOf(y), @TypeOf(x)), x, .{})),
                .float => switch (Coerce(@TypeOf(y), @TypeOf(x))) {
                    f16 => return cast(f16, atan2_32(cast(f32, y, .{}), cast(f32, x, .{})), .{}),
                    f32 => {
                        // glibc/sysdeps/ieee754/flt-32/e_atan2f.c
                        return atan2_32(y, x);
                    },
                    f64 => {
                        // glibc/sysdeps/ieee754/dbl-64/e_atan2.c
                        return atan2_64(y, x);
                    },
                    f80 => return cast(f80, atan2_128(cast(f128, y, .{}), cast(f128, x, .{})), .{}),
                    f128 => {
                        // glibc/sysdeps/ieee754/ldbl-128/e_atan2l.c
                        return atan2_128(y, x);
                    },
                    else => unreachable,
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

inline fn muldd(xh: f64, xl: f64, ch: f64, cl: f64, l: *f64) f64 {
    const ahlh: f64 = ch * xl;
    const alhh: f64 = cl * xh;
    const ahhh: f64 = ch * xh;
    var ahhl: f64 = @mulAdd(f64, ch, xh, -ahhh);
    ahhl += alhh + ahlh;
    const hh: f64 = ahhh + ahhl;
    l.* = (ahhh - hh) + ahhl;
    return hh;
}

fn polydd(xh: f64, xl: f64, n: i64, c: []const [2]f64, l: *f64) f64 {
    var i: u64 = cast(u64, n - 2, .{});
    var ch: f64 = c[i][0];
    var cl: f64 = c[i][1];
    while (i >= 0) {
        ch = muldd(xh, xl, ch, cl, &cl);
        const th: f64 = ch + c[i][0];
        const tl: f64 = (c[i][0] - th) + ch;
        ch = th;
        cl += tl + c[i][1];

        i -= 1;
    }

    l.* = cl;
    return ch;
}

// for y/x tiny, use Taylor approximation z - z^3/3 where z=y/x
fn atan2_32_tiny(y: f32, x: f32) f32 {
    const dy: f64 = y;
    const dx: f64 = x;
    const z: f64 = dy / dx;
    var e: f64 = @mulAdd(f64, -z, x, y);
    // z * x + e = y thus y/x = z + e/x
    const c: f64 = -0x1.5555555555555p-2; // -1/3 rounded to nearest
    const zz: f64 = z * z;
    const cz: f64 = c * z;
    e = e / x + cz * zz;
    var t: u64 = @bitCast(z);
    if ((t & 0xfffffff) == 0) { // boundary case

        // If z and e are of same sign (resp. of different signs), we increase
        // (resp. decrease) the significant of t by 1 to avoid a double-rounding
        // issue when rounding t to binary32.
        if (z * e > 0) {
            t += 1;
        } else {
            t -= 1;
        }
    }

    return cast(f32, @as(f64, @bitCast(t)), .{});
}

fn atan2_32(y: f32, x: f32) f32 {
    const cn: []const f64 =
        &.{ 0x1p+0, 0x1.40e0698f94c35p+1, 0x1.248c5da347f0dp+1, 0x1.d873386572976p-1, 0x1.46fa40b20f1dp-3, 0x1.33f5e041eed0fp-7, 0x1.546bbf28667c5p-14 };
    const cd: []const f64 =
        &.{ 0x1p+0, 0x1.6b8b143a3f6dap+1, 0x1.8421201d18ed5p+1, 0x1.8221d086914ebp+0, 0x1.670657e3a07bap-2, 0x1.0f4951fd1e72dp-5, 0x1.b3874b8798286p-11 };
    const m: []const f64 = &.{ 0, 1 };
    const pi: f64 = 0x1.921fb54442d18p+1;
    const pi2: f64 = 0x1.921fb54442d18p+0;
    const pi2l: f64 = 0x1.1a62633145c07p-54;
    const off: []const f64 = &.{ 0, pi2, pi, pi2, -0.0, -pi2, -pi, -pi2 };
    const offl: []const f64 =
        &.{ 0, pi2l, 2 * pi2l, pi2l, -0.0, -pi2l, -2 * pi2l, -pi2l };
    const sgn: []const f64 = &.{ 1, -1 };
    const ux: u32 = @bitCast(x);
    const uy: u32 = @bitCast(y);
    const ax: u32 = ux & (~@as(u32, 0) >> 1);
    const ay: u32 = uy & (~@as(u32, 0) >> 1);

    if (ay >= (0xff << 23) or ax >= (0xff << 23)) {
        // we use x+y below so that the invalid exception is set for (x,y) = (qnan,snan) or (snan,qnan)
        if (ay > (0xff << 23))
            return x + y; // nan

        if (ax > (0xff << 23))
            return x + y; // nan

        const yinf: bool = ay == (0xff << 23);
        const xinf: bool = ax == (0xff << 23);
        if (yinf and xinf) {
            if ((ux >> 31) != 0) {
                return cast(f32, 0x1.2d97c7f3321d2p+1 * sgn[uy >> 31], .{}); // +/-3pi/4
            } else {
                return cast(f32, 0x1.921fb54442d18p-1 * sgn[uy >> 31], .{}); // +/-pi/4
            }
        }

        if (xinf) {
            if ((ux >> 31) != 0) {
                return cast(f32, pi * sgn[uy >> 31], .{});
            } else {
                return cast(f32, 0 * sgn[uy >> 31], .{});
            }
        }
        if (yinf)
            return cast(f32, pi2 * sgn[uy >> 31], .{});
    }

    if (ay == 0) {
        if (ax == 0) {
            const i: u32 = (uy >> 31) * 4 + (ux >> 31) * 2;
            if ((ux >> 31) != 0) {
                return cast(f32, off[i] + offl[i], .{});
            } else {
                return cast(f32, off[i], .{});
            }
        }

        if ((ux >> 31) == 0)
            return cast(f32, 0.0 * sgn[uy >> 31], .{});
    }

    const gt: u32 = cast(u32, ay > ax, .{});
    const i: u32 = (uy >> 31) * 4 + (ux >> 31) * 2 + gt;

    const zx: f64 = x;
    const zy: f64 = y;
    var z: f64 = (m[gt] * zx + m[1 - gt] * zy) / (m[gt] * zy + m[1 - gt] * zx);
    // z = x/y if |y| > |x|, and z = y/x otherwise
    var r: f64 = undefined;
    const d: i64 = cast(i64, ax, .{}) - cast(i64, ay, .{});
    if (d < (27 << 23) and d > (-(27 << 23))) {
        const z2: f64 = z * z;
        const z4: f64 = z2 * z2;
        const z8: f64 = z4 * z4;
        // z2 cannot underflow, since for |y|=0x1p-149 and |x|=0x1.fffffep+127
        // we get |z| > 2^-277 thus z2 > 2^-554, but z4 and z8 might underflow,
        // which might give spurious underflow exceptions.
        var cn0: f64 = cn[0] + z2 * cn[1];
        const cn2: f64 = cn[2] + z2 * cn[3];
        var cn4: f64 = cn[4] + z2 * cn[5];
        const cn6: f64 = cn[6];
        cn0 += z4 * cn2;
        cn4 += z4 * cn6;
        cn0 += z8 * cn4;
        var cd0: f64 = cd[0] + z2 * cd[1];
        const cd2: f64 = cd[2] + z2 * cd[3];
        var cd4: f64 = cd[4] + z2 * cd[5];
        const cd6: f64 = cd[6];
        cd0 += z4 * cd2;
        cd4 += z4 * cd6;
        cd0 += z8 * cd4;
        r = cn0 / cd0;
    } else {
        r = 1;
    }

    z *= sgn[gt];
    r = z * r + off[i];

    if (((@as(u64, @bitCast(r)) + 8) & 0xfffffff) <= 16) {
        // check tiny y/x
        if (ay < ax and ((ax - ay) >> 23 >= 25))
            return atan2_32_tiny(y, x);

        var zh: f64 = undefined;
        var zl: f64 = undefined;
        if (gt == 0) {
            zh = zy / zx;
            zl = @mulAdd(f64, zh, -zx, zy) / zx;
        } else {
            zh = zx / zy;
            zl = @mulAdd(f64, zh, -zy, zx) / zy;
        }

        var z2l: f64 = undefined;
        const z2h: f64 = muldd(zh, zl, zh, zl, &z2l);
        const c: [32][2]f64 =
            .{
                .{ 0x1p+0, -0x1.8c1dac5492248p-87 },
                .{ -0x1.5555555555555p-2, -0x1.55553bf3a2abep-56 },
                .{ 0x1.999999999999ap-3, -0x1.99deed1ec9071p-57 },
                .{ -0x1.2492492492492p-3, -0x1.fd99c8d18269ap-58 },
                .{ 0x1.c71c71c71c717p-4, -0x1.651eee4c4d9dp-61 },
                .{ -0x1.745d1745d1649p-4, -0x1.632683d6c44a6p-58 },
                .{ 0x1.3b13b13b11c63p-4, 0x1.bf69c1f8af41dp-58 },
                .{ -0x1.11111110e6338p-4, 0x1.3c3e431e8bb68p-61 },
                .{ 0x1.e1e1e1dc45c4ap-5, -0x1.be2db05c77bbfp-59 },
                .{ -0x1.af286b8164b4fp-5, 0x1.a4673491f0942p-61 },
                .{ 0x1.86185e9ad4846p-5, 0x1.e12e32d79fceep-59 },
                .{ -0x1.642c6d5161faep-5, 0x1.3ce76c1ca03fp-59 },
                .{ 0x1.47ad6f277e5bfp-5, -0x1.abd8d85bdb714p-60 },
                .{ -0x1.2f64a2ee8896dp-5, 0x1.ef87d4b615323p-61 },
                .{ 0x1.1a6a2b31741b5p-5, 0x1.a5d9d973547eep-62 },
                .{ -0x1.07fbdad65e0a6p-5, -0x1.65ac07f5d35f4p-61 },
                .{ 0x1.ee9932a9a5f8bp-6, 0x1.f8b9623f6f55ap-61 },
                .{ -0x1.ce8b5b9584dc6p-6, 0x1.fe5af96e8ea2dp-61 },
                .{ 0x1.ac9cb288087b7p-6, -0x1.450cdfceaf5cap-60 },
                .{ -0x1.84b025351f3e6p-6, 0x1.579561b0d73dap-61 },
                .{ 0x1.52f5b8ecdd52bp-6, 0x1.036bd2c6fba47p-60 },
                .{ -0x1.163a8c44909dcp-6, 0x1.18f735ffb9f16p-60 },
                .{ 0x1.a400dce3eea6fp-7, -0x1.c90569c0c1b5cp-61 },
                .{ -0x1.1caa78ae6db3ap-7, -0x1.4c60f8161ea09p-61 },
                .{ 0x1.52672453c0731p-8, 0x1.834efb598c338p-62 },
                .{ -0x1.5850c5be137cfp-9, -0x1.445fc150ca7f5p-63 },
                .{ 0x1.23eb98d22e1cap-10, -0x1.388fbaf1d783p-64 },
                .{ -0x1.8f4e974a40741p-12, 0x1.271198a97da34p-66 },
                .{ 0x1.a5cf2e9cf76e5p-14, -0x1.887eb4a63b665p-68 },
                .{ -0x1.420c270719e32p-16, 0x1.efd595b27888bp-71 },
                .{ 0x1.3ba2d69b51677p-19, -0x1.4fb06829cdfc7p-73 },
                .{ -0x1.29b7e6f676385p-23, -0x1.a783b6de718fbp-77 },
            };

        var pl: f64 = undefined;
        var ph: f64 = polydd(z2h, z2l, 32, &c, &pl);
        zh *= sgn[gt];
        zl *= sgn[gt];
        ph = muldd(zh, zl, ph, pl, &pl);
        const sh: f64 = ph + off[i];
        const sl: f64 = ((off[i] - sh) + ph) + pl + offl[i];
        const rf: f32 = cast(f32, sh, .{});
        const th: f64 = rf;
        const dh: f64 = sh - th;
        var tm: f64 = dh + sl;
        var tth: u64 = @bitCast(th);

        if (th + th * 0x1p-60 == th - th * 0x1p-60) {
            tth &= 0x7ff << 52;
            tth -= 24 << 52;
            if (@abs(tm) > @as(f64, @bitCast(tth))) {
                tm *= 1.25;
            } else {
                tm *= 0.75;
            }
        }
        r = th + tm;
    }

    return cast(f32, r, .{});
}

// atan2 with max ULP of ~0.524 based on random sampling.
fn atan2_64(y: f64, x: f64) f64 {
    const TWO52: f64 = 0x1.0p52;
    // const TWOM1022: f64 = 0x1.0p-1022;

    const ep: i32 = 59768832; // 57*16**5
    const em: i32 = -59768832; // -57*16**5

    // x = NaN or y = NaN
    var num: [2]u32 = @bitCast(x);
    const ux: u32 = num[atnat2.HIGH_HALF];
    const dx: u32 = num[atnat2.LOW_HALF];
    if ((ux & 0x7ff00000) == 0x7ff00000) {
        if (((ux & 0x000fffff) | dx) != 0x00000000)
            return x + y;
    }

    num = @bitCast(y);
    const uy: u32 = num[atnat2.HIGH_HALF];
    const dy: u32 = num[atnat2.LOW_HALF];
    if ((uy & 0x7ff00000) == 0x7ff00000) {
        if (((uy & 0x000fffff) | dy) != 0x00000000)
            return y + y;
    }

    // y = +-0
    if (uy == 0x00000000) {
        if (dy == 0x00000000) {
            if ((ux & 0x80000000) == 0x00000000) {
                return 0;
            } else {
                return @bitCast(atnat2.opi);
            }
        }
    } else if (uy == 0x80000000) {
        if (dy == 0x00000000) {
            if ((ux & 0x80000000) == 0x00000000) {
                return -0.0;
            } else {
                return @bitCast(atnat2.mopi);
            }
        }
    }

    // x = +-0
    if (x == 0) {
        if ((uy & 0x80000000) == 0x00000000) {
            return @bitCast(atnat2.hpi);
        } else {
            return @bitCast(atnat2.mhpi);
        }
    }

    // x = +-INF
    if (ux == 0x7ff00000) {
        if (dx == 0x00000000) {
            if (uy == 0x7ff00000) {
                if (dy == 0x00000000)
                    return @bitCast(atnat2.qpi);
            } else if (uy == 0xfff00000) {
                if (dy == 0x00000000)
                    return @bitCast(atnat2.mqpi);
            } else {
                if ((uy & 0x80000000) == 0x00000000) {
                    return 0;
                } else {
                    return -0.0;
                }
            }
        }
    } else if (ux == 0xfff00000) {
        if (dx == 0x00000000) {
            if (uy == 0x7ff00000) {
                if (dy == 0x00000000)
                    return @bitCast(atnat2.tqpi);
            } else if (uy == 0xfff00000) {
                if (dy == 0x00000000)
                    return @bitCast(atnat2.mtqpi);
            } else {
                if ((uy & 0x80000000) == 0x00000000) {
                    return @bitCast(atnat2.opi);
                } else {
                    return @bitCast(atnat2.mopi);
                }
            }
        }
    }

    // y=+-INF
    if (uy == 0x7ff00000) {
        if (dy == 0x00000000)
            return @bitCast(atnat2.hpi);
    } else if (uy == 0xfff00000) {
        if (dy == 0x00000000)
            return @bitCast(atnat2.mhpi);
    }

    // either x/y or y/x is very close to zero
    var ax: f64 = if (x < 0) -x else x;
    var ay: f64 = if (y < 0) -y else y;
    const de: i32 = cast(i32, uy & 0x7ff00000, .{}) - cast(i32, ux & 0x7ff00000, .{});
    if (de >= ep) {
        return if (y > 0) @bitCast(atnat2.hpi) else @bitCast(atnat2.mhpi);
    } else if (de <= em) {
        if (x > 0) {
            const z: f64 = ay / ax;
            const ret: f64 = math.copysign(z, y);
            if (@abs(ret) < std.math.floatMin(f64)) {
                const vret: f64 = if (ret != 0) ret else std.math.floatMin(f64);
                const force_underflow = vret * vret;
                std.mem.doNotOptimizeAway(force_underflow);
            }

            return ret;
        } else {
            return if (y > 0) @bitCast(atnat2.opi) else @bitCast(atnat2.mopi);
        }
    }

    // if either x or y is extremely close to zero, scale abs(x), abs(y).
    if (ax < @as(f64, @bitCast(atnat2.twom500)) or ay < @as(f64, @bitCast(atnat2.twom500))) {
        ax *= @bitCast(atnat2.two500);
        ay *= @bitCast(atnat2.two500);
    }

    // Likewise for large x and y.
    if (ax > @as(f64, @bitCast(atnat2.two500)) or ay > @as(f64, @bitCast(atnat2.two500))) {
        ax *= @bitCast(atnat2.twom500);
        ay *= @bitCast(atnat2.twom500);
    }

    // x, y which are neither special nor extreme
    var u: f64 = undefined;
    var v: f64 = undefined;
    var vv: f64 = undefined;
    var du: f64 = undefined;
    if (ay < ax) {
        u = ay / ax;
        dla.emulv(ax, u, &v, &vv);
        du = ((ay - v) - vv) / ax;
    } else {
        u = ax / ay;
        dla.emulv(ay, u, &v, &vv);
        du = ((ax - v) - vv) / ay;
    }

    if (x > 0) {
        // (i)   x>0, abs(y)< abs(x):  atan(ay/ax)
        if (ay < ax) {
            if (u < @as(f64, @bitCast(atnat2.inv16))) {
                v = u * u;

                const zz: f64 = du + u * v * (@as(f64, @bitCast(atnat2.d3)) + v * (@as(f64, @bitCast(atnat2.d5)) + v * (@as(f64, @bitCast(atnat2.d7)) + v * (@as(f64, @bitCast(atnat2.d9)) + v * (@as(f64, @bitCast(atnat2.d11)) + v * @as(f64, @bitCast(atnat2.d13)))))));

                const z: f64 = u + zz;
                // Max ULP is 0.504.
                return math.copysign(z, y);
            }

            var i: i32 = cast(i32, (TWO52 + 256 * u) - TWO52, .{});
            i -= 16;
            const t3: f64 = u - @as(f64, @bitCast(atnat2.cij[@intCast(i)][0]));
            var dv: f64 = undefined;
            dla.eadd(t3, du, &v, &dv);
            const t1: f64 = @as(f64, @bitCast(atnat2.cij[@intCast(i)][1]));
            const t2: f64 = @as(f64, @bitCast(atnat2.cij[@intCast(i)][2]));
            const zz: f64 = v * t2 + (dv * t2 + v * v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][3])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][4])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][5])) + v * @as(f64, @bitCast(atnat2.cij[@intCast(i)][6]))))));
            const z: f64 = t1 + zz;
            // Max ULP is 0.56.
            return math.copysign(z, y);
        }

        // (ii)  x > 0, abs(x) <= abs(y):  pi/2 - atan(ax/ay)
        if (u < @as(f64, @bitCast(atnat2.inv16))) {
            v = u * u;
            const zz: f64 = u * v * (@as(f64, @bitCast(atnat2.d3)) + v * (@as(f64, @bitCast(atnat2.d5)) + v * (@as(f64, @bitCast(atnat2.d7)) + v * (@as(f64, @bitCast(atnat2.d9)) + v * (@as(f64, @bitCast(atnat2.d11)) + v * @as(f64, @bitCast(atnat2.d13)))))));
            var t2: f64 = undefined;
            var cor: f64 = undefined;
            dla.esub(@as(f64, @bitCast(atnat2.hpi)), u, &t2, &cor);
            const t3: f64 = ((@as(f64, @bitCast(atnat2.hpi1)) + cor) - du) - zz;
            const z: f64 = t2 + t3;
            // Max ULP is 0.501.
            return math.copysign(z, y);
        }

        var i: i32 = cast(i32, (TWO52 + 256 * u) - TWO52, .{});
        i -= 16;
        v = (u - @as(f64, @bitCast(atnat2.cij[@intCast(i)][0]))) + du;

        const zz: f64 = @as(f64, @bitCast(atnat2.hpi1)) - v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][2])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][3])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][4])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][5])) + v * @as(f64, @bitCast(atnat2.cij[@intCast(i)][6]))))));
        const t1: f64 = @as(f64, @bitCast(atnat2.hpi)) - @as(f64, @bitCast(atnat2.cij[@intCast(i)][1]));
        const z: f64 = t1 + zz;
        // Max ULP is 0.503.
        return math.copysign(z, y);
    }

    // (iii) x < 0, abs(x) < abs(y):  pi/2 + atan(ax/ay)
    if (ax < ay) {
        if (u < @as(f64, @bitCast(atnat2.inv16))) {
            v = u * u;
            const zz: f64 = u * v * (@as(f64, @bitCast(atnat2.d3)) + v * (@as(f64, @bitCast(atnat2.d5)) + v * (@as(f64, @bitCast(atnat2.d7)) + v * (@as(f64, @bitCast(atnat2.d9)) + v * (@as(f64, @bitCast(atnat2.d11)) + v * @as(f64, @bitCast(atnat2.d13)))))));
            var t2: f64 = undefined;
            var cor: f64 = undefined;
            dla.eadd(@as(f64, @bitCast(atnat2.hpi)), u, &t2, &cor);
            const t3: f64 = ((@as(f64, @bitCast(atnat2.hpi1)) + cor) + du) + zz;
            const z: f64 = t2 + t3;
            // Max ULP is 0.501.
            return math.copysign(z, y);
        }

        var i: i32 = cast(i32, (TWO52 + 256 * u) - TWO52, .{});
        i -= 16;
        v = (u - @as(f64, @bitCast(atnat2.cij[@intCast(i)][0]))) + du;
        const zz: f64 = @as(f64, @bitCast(atnat2.hpi1)) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][2])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][3])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][4])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][5])) + v * @as(f64, @bitCast(atnat2.cij[@intCast(i)][6]))))));
        const t1: f64 = @as(f64, @bitCast(atnat2.hpi)) + @as(f64, @bitCast(atnat2.cij[@intCast(i)][1]));
        const z: f64 = t1 + zz;
        // Max ULP is 0.503.
        return math.copysign(z, y);
    }

    // (iv)  x < 0, abs(y) <= abs(x):  pi - atan(ax/ay)
    if (u < @as(f64, @bitCast(atnat2.inv16))) {
        v = u * u;
        const zz: f64 = u * v * (@as(f64, @bitCast(atnat2.d3)) + v * (@as(f64, @bitCast(atnat2.d5)) + v * (@as(f64, @bitCast(atnat2.d7)) + v * (@as(f64, @bitCast(atnat2.d9)) + v * (@as(f64, @bitCast(atnat2.d11)) + v * @as(f64, @bitCast(atnat2.d13)))))));
        var t2: f64 = undefined;
        var cor: f64 = undefined;
        dla.esub(@as(f64, @bitCast(atnat2.opi)), u, &t2, &cor);
        const t3: f64 = ((@as(f64, @bitCast(atnat2.opi1)) + cor) - du) - zz;
        const z: f64 = t2 + t3;
        // Max ULP is 0.501.
        return math.copysign(z, y);
    }

    var i: i32 = cast(i32, (TWO52 + 256 * u) - TWO52, .{});
    i -= 16;
    v = (u - @as(f64, @bitCast(atnat2.cij[@intCast(i)][0]))) + du;
    const zz: f64 = @as(f64, @bitCast(atnat2.opi1)) - v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][2])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][3])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][4])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][5])) + v * @as(f64, @bitCast(atnat2.cij[@intCast(i)][6]))))));
    const t1: f64 = @as(f64, @bitCast(atnat2.opi)) - @as(f64, @bitCast(atnat2.cij[@intCast(i)][1]));
    const z: f64 = t1 + zz;
    // Max ULP is 0.502.
    return math.copysign(z, y);
}

fn atan2_128(y: f128, x: f128) f128 {
    @setRuntimeSafety(false);
    const tiny: f128 = 1.0e-4900;
    const pi_o_4: f128 = 7.85398163397448309615660845819875699e-01; // 3ffe921fb54442d18469898cc51701b8
    const pi_o_2: f128 = 1.57079632679489661923132169163975140e+00; // 3fff921fb54442d18469898cc51701b8
    const pi: f128 = 3.14159265358979323846264338327950280e+00; // 4000921fb54442d18469898cc51701b8
    const pi_lo: f128 = 8.67181013012378102479704402604335225e-35; // 3f8dcd129024e088a67cc74020bbea64

    var hx: i64 = undefined;
    var lx: i64 = undefined;
    ldbl128.getWords(&hx, &lx, x);
    const ix: i64 = hx & 0x7fffffffffffffff;

    var hy: i64 = undefined;
    var ly: i64 = undefined;
    ldbl128.getWords(&hy, &ly, y);
    const iy: i64 = hy & 0x7fffffffffffffff;

    if ((ix | ((lx | -lx) >> 63)) > 0x7fff000000000000 or (iy | ((ly | -ly) >> 63)) > 0x7fff000000000000) // x or y is NaN
        return x + y;

    if (((hx - 0x3fff000000000000) | lx) == 0) return math.atan(y); // x = 1.0

    const m: i64 = @bitCast(((hy >> 63) & 1) | ((hx >> 62) & 2)); // 2*sign(x)+sign(y)

    // when y = 0
    if ((iy | ly) == 0) {
        switch (m) {
            0, 1 => return y, // atan(+-0,+anything)=+-0
            2 => return pi + tiny, // atan(+0,-anything) = pi
            3 => return -pi - tiny, // atan(-0,-anything) =-pi
            else => unreachable,
        }
    }

    // when x = 0
    if ((ix | lx) == 0) return if (hy < 0) -pi_o_2 - tiny else pi_o_2 + tiny;

    // when x is INF
    if (ix == 0x7fff000000000000) {
        if (iy == 0x7fff000000000000) {
            switch (m) {
                0 => return pi_o_4 + tiny, // atan(+INF,+INF)
                1 => return -pi_o_4 - tiny, // atan(-INF,+INF)
                2 => return 3 * pi_o_4 + tiny, // atan(+INF,-INF)
                3 => return -3 * pi_o_4 - tiny, // atan(-INF,-INF)
                else => unreachable,
            }
        } else {
            switch (m) {
                0 => return 0, // atan(+...,+INF)
                1 => return -0.0, // atan(-...,+INF)
                2 => return pi + tiny, // atan(+...,-INF)
                3 => return -pi - tiny, // atan(-...,-INF)
                else => unreachable,
            }
        }
    }

    // when y is INF
    if (iy == 0x7fff000000000000) return if (hy < 0) -pi_o_2 - tiny else pi_o_2 + tiny;

    // compute y/x
    const k: i64 = @bitCast((iy - ix) >> 48);
    var z: f128 = undefined;
    if (k > 120) { // |y/x| >  2**120
        z = pi_o_2 + 0.5 * pi_lo;
    } else if (hx < 0 and k < -120) { // |y|/x < -2**120
        z = 0;
    } else { // safe to do y/x
        z = math.atan(@abs(y / x));
    }

    switch (m) {
        0 => return z, // atan(+,+)
        1 => return -z, // atan(-,+)

        2 => return pi - (z - pi_lo), // atan(+,-)
        else => return (z - pi_lo) - pi, // case 3, atan(-,-)
    }
}

test atan2 {
    @setEvalBranchQuota(40000);
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x1p+0), std.math.inf(f32)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), std.math.inf(f32)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x4p-128), std.math.inf(f32)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x8p-152), std.math.inf(f32)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0xf.fffffp+124), std.math.inf(f32)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x1p+0), std.math.inf(f32)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x0p+0), std.math.inf(f32)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x4p-128), std.math.inf(f32)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x8p-152), std.math.inf(f32)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0xf.fffffp+124), std.math.inf(f32)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(std.math.inf(f32), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(-std.math.inf(f32), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x1p+0), -std.math.inf(f32)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x0p+0), -std.math.inf(f32)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x4p-128), -std.math.inf(f32)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x8p-152), -std.math.inf(f32)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0xf.fffffp+124), -std.math.inf(f32)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x1p+0), -std.math.inf(f32)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x0p+0), -std.math.inf(f32)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x4p-128), -std.math.inf(f32)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x8p-152), -std.math.inf(f32)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0xf.fffffp+124), -std.math.inf(f32)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(std.math.inf(f32), std.math.inf(f32)));
    try std.testing.expectEqual(-0xc.90fdbp-4, atan2(-std.math.inf(f32), std.math.inf(f32)));
    try std.testing.expectEqual(0x2.5b2f9p+0, atan2(std.math.inf(f32), -std.math.inf(f32)));
    try std.testing.expectEqual(-0x2.5b2f9p+0, atan2(-std.math.inf(f32), -std.math.inf(f32)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x0p+0), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x0p+0), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x0p+0), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x0p+0), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x0p+0), @as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x0p+0), @as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0x1p+0), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0x1p+0), @as(f32, -0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(@as(f32, -0x1p+0), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(@as(f32, -0x1p+0), @as(f32, -0x0p+0)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0xf.fffffp+124), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.5b2f9p+0, atan2(@as(f32, 0xf.fffffp+124), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0xc.90fdbp-4, atan2(@as(f32, -0xf.fffffp+124), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.5b2f9p+0, atan2(@as(f32, -0xf.fffffp+124), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0xf.fffffp+124), @as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(@as(f32, -0xf.fffffp+124), @as(f32, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(@as(f32, -0xf.fffffp+124), @as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0xf.fffffp+124), @as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0xf.fffffp+124), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(@as(f32, -0xf.fffffp+124), @as(f32, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(@as(f32, -0xf.fffffp+124), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0xf.fffffp+124), @as(f32, -0x8p-152)));
    try std.testing.expectEqual(0xa.4bc7dp-4, atan2(@as(f32, 0xcp-4), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0xa.4bc7dp-4, atan2(@as(f32, -0xcp-4), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2.7f82ecp+0, atan2(@as(f32, 0xcp-4), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(-0x2.7f82ecp+0, atan2(@as(f32, -0xcp-4), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x1.91ef0ep+0, atan2(@as(f32, 0x6.4p-4), @as(f32, 0x1.30164ap-12)));
    try std.testing.expectEqual(0x1.91ef0ep+0, atan2(@as(f32, 0x6.4p-4), @as(f32, 0x1.301648p-12)));
    try std.testing.expectEqual(0xf.b437ap-4, atan2(@as(f32, 0x1.64p+0), @as(f32, 0xe.ep-4)));
    try std.testing.expectEqual(-0x1.cdaa9ep+0, atan2(@as(f32, -0x1.effe8p-8), @as(f32, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9ep+0, atan2(@as(f32, -0x1.effe8p-8), @as(f32, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9ep+0, atan2(@as(f32, -0x1.effe82p-8), @as(f32, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9ep+0, atan2(@as(f32, -0x1.effe82p-8), @as(f32, -0x7.57d1ep-12)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0x1.000002p+0), @as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdcp-4, atan2(@as(f32, 0x1.000002p+0), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdap-4, atan2(@as(f32, 0x1p+0), @as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0x1p+0), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.9c22cep-4, atan2(@as(f32, 0x4.c3841p-4), @as(f32, 0x2.f2f308p+0)));
    try std.testing.expectEqual(-0x1.1dd4c4p-4, atan2(@as(f32, -0xe.cf143p-40), @as(f32, 0xd.3de7ap-36)));
    try std.testing.expectEqual(0x2.7c1784p-4, atan2(@as(f32, 0x5.576cf8p-4), @as(f32, 0x2.21e65p+0)));
    try std.testing.expectEqual(-0x2.1dbac4p-4, atan2(@as(f32, -0x4.29411p-4), @as(f32, 0x1.f4755cp+0)));
    try std.testing.expectEqual(-0x1.921fb6p+0, atan2(@as(f32, -0xa.b4101p+20), @as(f32, -0xf.9c4c8p-4)));
    try std.testing.expectEqual(0x9.23e97p-8, atan2(@as(f32, 0x4.251bb8p-4), @as(f32, 0x7.40ac68p+0)));
    try std.testing.expectEqual(0x1.fd0a44p-4, atan2(@as(f32, 0x1.47239ep+68), @as(f32, 0xa.3ac3cp+68)));
    try std.testing.expectEqual(-0x1.de8936p-4, atan2(@as(f32, -0x6.b0794p-4), @as(f32, 0x3.8ff10cp+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x0p+0), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x8p-152), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.e5617ap+0, atan2(@as(f32, 0x3.f16f1p+0), @as(f32, -0x1.546056p+0)));
    try std.testing.expectEqual(-0x1.921fb4p+0, atan2(@as(f32, -0x1.9e67cp-24), @as(f32, 0x7.40bb4p-52)));
    try std.testing.expectEqual(-0x3.f2b38p-4, atan2(@as(f32, -0x3.f39e9p+48), @as(f32, 0xf.b02ccp+48)));
    try std.testing.expectEqual(0xc.92e27p-4, atan2(@as(f32, 0x6.f2aca8p-56), @as(f32, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e27p-4, atan2(@as(f32, 0x6.f2aca8p-56), @as(f32, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e26p-4, atan2(@as(f32, 0x6.f2acap-56), @as(f32, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e27p-4, atan2(@as(f32, 0x6.f2acap-56), @as(f32, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0x8p-152), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0x8p-152), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0x8p-152), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0x8p-152), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0x8p-152), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0x8p-152), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0x8p-152), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb6p+0, atan2(@as(f32, 0x8p-152), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x0p+0), @as(f32, 0x0p+0)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0x4p-128), @as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x2.5b2f9p+0, atan2(@as(f32, 0x4p-128), @as(f32, -0x4p-128)));
    try std.testing.expectEqual(-0xc.90fdbp-4, atan2(@as(f32, -0x4p-128), @as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x2.5b2f9p+0, atan2(@as(f32, -0x4p-128), @as(f32, -0x4p-128)));
    try std.testing.expectEqual(0xc.90fdbp-4, atan2(@as(f32, 0x8p-152), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x2.5b2f9p+0, atan2(@as(f32, 0x8p-152), @as(f32, -0x8p-152)));
    try std.testing.expectEqual(-0xc.90fdbp-4, atan2(@as(f32, -0x8p-152), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x2.5b2f9p+0, atan2(@as(f32, -0x8p-152), @as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb4p+0, atan2(@as(f32, 0x4p-128), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb8p+0, atan2(@as(f32, 0x4p-128), @as(f32, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb4p+0, atan2(@as(f32, -0x4p-128), @as(f32, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb8p+0, atan2(@as(f32, -0x4p-128), @as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x2p-24, atan2(@as(f32, 0x8p-152), @as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x3.243f68p+0, atan2(@as(f32, 0x8p-152), @as(f32, -0x4p-128)));
    try std.testing.expectEqual(-0x2p-24, atan2(@as(f32, -0x8p-152), @as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x3.243f68p+0, atan2(@as(f32, -0x8p-152), @as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x1p+0), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x1p+0), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x4p-128), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x4p-128), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x8p-152), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x8p-152), @as(f32, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p-128, atan2(@as(f32, 0x1p+0), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p-128, atan2(@as(f32, -0x1p+0), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x4p-128), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x4p-128), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f32, 0x8p-152), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f32, -0x8p-152), @as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x4p-128, atan2(@as(f32, 0x4p-128), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-128, atan2(@as(f32, -0x4p-128), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x8p-152, atan2(@as(f32, 0x8p-152), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-152, atan2(@as(f32, -0x8p-152), @as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x4p-128), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x4p-128), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6cp+0, atan2(@as(f32, 0x8p-152), @as(f32, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6cp+0, atan2(@as(f32, -0x8p-152), @as(f32, -0x1p+0)));

    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x1p+0), std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x4p-128), std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x4p-1024), std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x8p-972), std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x8p-152), std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x4p-1076), std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0xf.fffffp+124), std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x1p+0), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x0p+0), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x4p-128), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x4p-1024), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x8p-972), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x8p-152), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x4p-1076), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0xf.fffffp+124), std.math.inf(f64)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), std.math.inf(f64)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(std.math.inf(f64), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(-std.math.inf(f64), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x1p+0), -std.math.inf(f64)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x0p+0), -std.math.inf(f64)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-128), -std.math.inf(f64)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1024), -std.math.inf(f64)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-972), -std.math.inf(f64)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-152), -std.math.inf(f64)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1076), -std.math.inf(f64)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0xf.fffffp+124), -std.math.inf(f64)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x1p+0), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x0p+0), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-128), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1024), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-972), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-152), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1076), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0xf.fffffp+124), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), -std.math.inf(f64)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(std.math.inf(f64), std.math.inf(f64)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan2(-std.math.inf(f64), std.math.inf(f64)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, atan2(std.math.inf(f64), -std.math.inf(f64)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, atan2(-std.math.inf(f64), -std.math.inf(f64)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x0p+0), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x0p+0), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x0p+0), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x0p+0), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x0p+0), @as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x0p+0), @as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x1p+0), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x1p+0), @as(f64, -0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x1p+0), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x1p+0), @as(f64, -0x0p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0xf.fffffp+124), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffff00000008p-900, atan2(@as(f64, 0xf.fffffp+124), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan2(@as(f64, -0xf.fffffp+124), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fffff00000008p-900, atan2(@as(f64, -0xf.fffffp+124), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.fffffp+124), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0xf.ffffffffffff8p+1020), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.fffffp+124), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xf.ffffffffffff8p+1020), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0xa.4bc7d1934f708p-4, atan2(@as(f64, 0xcp-4), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0xa.4bc7d1934f708p-4, atan2(@as(f64, -0xcp-4), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x2.7f82ed6f50acp+0, atan2(@as(f64, 0xcp-4), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(-0x2.7f82ed6f50acp+0, atan2(@as(f64, -0xcp-4), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x1.91ef0ddcd8c8cp+0, atan2(@as(f64, 0x6.4p-4), @as(f64, 0x1.30164ap-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2ab45p+0, atan2(@as(f64, 0x6.4p-4), @as(f64, 0x1.301648p-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052cp+0, atan2(@as(f64, 0x6.4p-4), @as(f64, 0x1.30164840e171ap-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052cp+0, atan2(@as(f64, 0x6.4p-4), @as(f64, 0x1.30164840e1719p-12)));
    try std.testing.expectEqual(0xf.b437a72087798p-4, atan2(@as(f64, 0x1.64p+0), @as(f64, 0xe.ep-4)));
    try std.testing.expectEqual(-0x1.cdaa9db2d1b45p+0, atan2(@as(f64, -0x1.effe8p-8), @as(f64, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9df15fde1p+0, atan2(@as(f64, -0x1.effe8p-8), @as(f64, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46ep+0, atan2(@as(f64, -0x1.effe8p-8), @as(f64, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46ep+0, atan2(@as(f64, -0x1.effe8p-8), @as(f64, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9d778c452p+0, atan2(@as(f64, -0x1.effe82p-8), @as(f64, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db61a6ebp+0, atan2(@as(f64, -0x1.effe82p-8), @as(f64, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da6e6d79p+0, atan2(@as(f64, -0x1.effe82p-8), @as(f64, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da6e6d79p+0, atan2(@as(f64, -0x1.effe82p-8), @as(f64, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9d786fcfap+0, atan2(@as(f64, -0x1.effe81f852716p-8), @as(f64, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db6fdf93p+0, atan2(@as(f64, -0x1.effe81f852716p-8), @as(f64, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca621p+0, atan2(@as(f64, -0x1.effe81f852716p-8), @as(f64, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca621p+0, atan2(@as(f64, -0x1.effe81f852716p-8), @as(f64, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9d786fcfap+0, atan2(@as(f64, -0x1.effe81f852717p-8), @as(f64, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db6fdf93p+0, atan2(@as(f64, -0x1.effe81f852717p-8), @as(f64, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca621p+0, atan2(@as(f64, -0x1.effe81f852717p-8), @as(f64, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca621p+0, atan2(@as(f64, -0x1.effe81f852717p-8), @as(f64, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x1.000002p+0), @as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdba22167cp-4, atan2(@as(f64, 0x1.000002p+0), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdba22167b8p-4, atan2(@as(f64, 0x1.000002p+0), @as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fd9a22169cp-4, atan2(@as(f64, 0x1p+0), @as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x1p+0), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168b8p-4, atan2(@as(f64, 0x1p+0), @as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fd9a22169c8p-4, atan2(@as(f64, 0x1.0000000000001p+0), @as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c8p-4, atan2(@as(f64, 0x1.0000000000001p+0), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x1.0000000000001p+0), @as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x1.9c22ce44a722ap-4, atan2(@as(f64, 0x4.c3841p-4), @as(f64, 0x2.f2f308p+0)));
    try std.testing.expectEqual(-0x1.1dd4c4e264577p-4, atan2(@as(f64, -0xe.cf143p-40), @as(f64, 0xd.3de7ap-36)));
    try std.testing.expectEqual(0x2.7c1782a75e16cp-4, atan2(@as(f64, 0x5.576cf8p-4), @as(f64, 0x2.21e65p+0)));
    try std.testing.expectEqual(-0x2.1dbac4fa4bfecp-4, atan2(@as(f64, -0x4.29411p-4), @as(f64, 0x1.f4755cp+0)));
    try std.testing.expectEqual(-0x1.921fb6b9a118dp+0, atan2(@as(f64, -0xa.b4101p+20), @as(f64, -0xf.9c4c8p-4)));
    try std.testing.expectEqual(0x9.23e97736442d8p-8, atan2(@as(f64, 0x4.251bb8p-4), @as(f64, 0x7.40ac68p+0)));
    try std.testing.expectEqual(0x1.fd0a44d0aba44p-4, atan2(@as(f64, 0x1.47239ep+68), @as(f64, 0xa.3ac3cp+68)));
    try std.testing.expectEqual(-0x1.de89352a0e839p-4, atan2(@as(f64, -0x6.b0794p-4), @as(f64, 0x3.8ff10cp+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x0p+0), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x8.000008000008p-280, atan2(@as(f64, -0x8p-152), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x7.15e7b61fff36cp-852, atan2(@as(f64, -0x7.15e7af0a1780cp-724), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.e5617957818bcp+0, atan2(@as(f64, 0x3.f16f1p+0), @as(f64, -0x1.546056p+0)));
    try std.testing.expectEqual(-0x1.921fb4fc92693p+0, atan2(@as(f64, -0x1.9e657cp-24), @as(f64, 0x7.40bb4p-52)));
    try std.testing.expectEqual(-0x3.f2b37e0ca7ac4p-4, atan2(@as(f64, -0x3.f39e9p+48), @as(f64, 0xf.b02ccp+48)));
    try std.testing.expectEqual(0xc.92e26b8fc63d8p-4, atan2(@as(f64, 0x6.f2aca8p-56), @as(f64, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e274c80fb8p-4, atan2(@as(f64, 0x6.f2aca8p-56), @as(f64, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e270ff3b358p-4, atan2(@as(f64, 0x6.f2aca8p-56), @as(f64, 0x6.f107d348a52ep-56)));
    try std.testing.expectEqual(0xc.92e26259ab2ep-4, atan2(@as(f64, 0x6.f2acap-56), @as(f64, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e26b91f4a8p-4, atan2(@as(f64, 0x6.f2acap-56), @as(f64, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e267c920258p-4, atan2(@as(f64, 0x6.f2acap-56), @as(f64, 0x6.f107d348a52ep-56)));
    try std.testing.expectEqual(0xc.92e26ae105628p-4, atan2(@as(f64, 0x6.f2aca7683a51cp-56), @as(f64, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e274194edc8p-4, atan2(@as(f64, 0x6.f2aca7683a51cp-56), @as(f64, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e270507a5ap-4, atan2(@as(f64, 0x6.f2aca7683a51cp-56), @as(f64, 0x6.f107d348a52ep-56)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x8p-152), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0x3.423123d8009a8p-200, atan2(@as(f64, 0x1.a11891ec004d4p-348), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x1.a11891ec004d4p-348), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xd.334f8c6d4a1bp-4, atan2(@as(f64, 0x1.a11891ec004d4p-348), @as(f64, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x8p-152), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0x3.706ddacf17c52p-440, atan2(@as(f64, 0x1.b836ed678be29p-588), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x1.b836ed678be29p-588), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xc.932e87412478p-4, atan2(@as(f64, 0x1.b836ed678be29p-588), @as(f64, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x8p-152), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0x1.a83f842ef3f73p-484, atan2(@as(f64, 0xd.41fc21779fb98p-636), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0xd.41fc21779fb98p-636), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xc.941e0644c21bp-4, atan2(@as(f64, 0xd.41fc21779fb98p-636), @as(f64, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x8p-152), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x0p+0), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-928, atan2(@as(f64, 0x4p-1076), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x4p-1076), @as(f64, 0x0p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x4p-1076), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x4p-128), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x4p-128), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x4p-128), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1p-896, atan2(@as(f64, 0x4p-1024), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x4p-1024), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-56, atan2(@as(f64, 0x4p-1024), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x2p-844, atan2(@as(f64, 0x8p-972), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-972), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x8p-972), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, atan2(@as(f64, 0x4p-128), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x4p-128), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x4p-128), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1024), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, atan2(@as(f64, 0x4p-1024), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1024), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-972), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d19p+0, atan2(@as(f64, 0x8p-972), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, atan2(@as(f64, 0x8p-972), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan2(@as(f64, -0x4p-128), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x4p-128), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x4p-128), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-896, atan2(@as(f64, -0x4p-1024), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan2(@as(f64, -0x4p-1024), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x8p-56, atan2(@as(f64, -0x4p-1024), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x2p-844, atan2(@as(f64, -0x8p-972), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-972), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan2(@as(f64, -0x8p-972), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, atan2(@as(f64, -0x4p-128), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x4p-128), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x4p-128), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1024), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, atan2(@as(f64, -0x4p-1024), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1024), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-972), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d19p+0, atan2(@as(f64, -0x8p-972), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, atan2(@as(f64, -0x8p-972), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x8p-152), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-928, atan2(@as(f64, 0x4p-1076), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0xc.90fdaa22168cp-4, atan2(@as(f64, 0x4p-1076), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, atan2(@as(f64, 0x8p-152), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1076), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a4p+0, atan2(@as(f64, 0x4p-1076), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan2(@as(f64, -0x8p-152), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-152), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-928, atan2(@as(f64, -0x4p-1076), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0xc.90fdaa22168cp-4, atan2(@as(f64, -0x4p-1076), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, atan2(@as(f64, -0x8p-152), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-152), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1076), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a4p+0, atan2(@as(f64, -0x4p-1076), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb34442d18p+0, atan2(@as(f64, 0x4p-128), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x4p-128), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x8p-876, atan2(@as(f64, 0x4p-1024), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d17p+0, atan2(@as(f64, 0x4p-1024), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1p-820, atan2(@as(f64, 0x8p-972), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-972), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb74442d18p+0, atan2(@as(f64, 0x4p-128), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x4p-128), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1024), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d19p+0, atan2(@as(f64, 0x4p-1024), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-972), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-972), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb34442d18p+0, atan2(@as(f64, -0x4p-128), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x4p-128), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x8p-876, atan2(@as(f64, -0x4p-1024), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d17p+0, atan2(@as(f64, -0x4p-1024), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x1p-820, atan2(@as(f64, -0x8p-972), @as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-972), @as(f64, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb74442d18p+0, atan2(@as(f64, -0x4p-128), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x4p-128), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1024), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d19p+0, atan2(@as(f64, -0x4p-1024), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-972), @as(f64, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-972), @as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x1.fffffffffffd5p-24, atan2(@as(f64, 0x8p-152), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1p-948, atan2(@as(f64, 0x4p-1076), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1p-52, atan2(@as(f64, 0x4p-1076), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x8p-108, atan2(@as(f64, 0x4p-1076), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x3.243f688885a3p+0, atan2(@as(f64, 0x8p-152), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, atan2(@as(f64, 0x8p-152), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1076), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1076), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1076), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0x1.fffffffffffd5p-24, atan2(@as(f64, -0x8p-152), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-152), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-152), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-948, atan2(@as(f64, -0x4p-1076), @as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x1p-52, atan2(@as(f64, -0x4p-1076), @as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x8p-108, atan2(@as(f64, -0x4p-1076), @as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x3.243f688885a3p+0, atan2(@as(f64, -0x8p-152), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-152), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18p+0, atan2(@as(f64, -0x8p-152), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1076), @as(f64, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1076), @as(f64, -0x4p-1024)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1076), @as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x1p+0), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x1p+0), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x1p+0), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x1p+0), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-128), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-128), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1024), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1024), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-972), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-972), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-128), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-128), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1024), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1024), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-972), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-972), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-152), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-152), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1076), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1076), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-152), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-152), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1076), @as(f64, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1076), @as(f64, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.000001000001p-128, atan2(@as(f64, 0x1p+0), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1p-1024, atan2(@as(f64, 0x1p+0), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.000001000001p-128, atan2(@as(f64, -0x1p+0), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1p-1024, atan2(@as(f64, -0x1p+0), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4.000004000004p-256, atan2(@as(f64, 0x4p-128), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x4p-128), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x4p-1024), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x4p-1024), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x8p-972), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x8p-972), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x4.000004000004p-256, atan2(@as(f64, -0x4p-128), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x4p-128), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x4p-1024), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x4p-1024), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x8p-972), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x8p-972), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x8.000008000008p-280, atan2(@as(f64, 0x8p-152), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x8p-152), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x4p-1076), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f64, 0x4p-1076), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x8.000008000008p-280, atan2(@as(f64, -0x8p-152), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x8p-152), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x4p-1076), @as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f64, -0x4p-1076), @as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4p-128, atan2(@as(f64, 0x4p-128), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x4p-1024, atan2(@as(f64, 0x4p-1024), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x8p-972, atan2(@as(f64, 0x8p-972), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-128, atan2(@as(f64, -0x4p-128), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-1024, atan2(@as(f64, -0x4p-1024), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-972, atan2(@as(f64, -0x8p-972), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x8p-152, atan2(@as(f64, 0x8p-152), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x4p-1076, atan2(@as(f64, 0x4p-1076), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-152, atan2(@as(f64, -0x8p-152), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-1076, atan2(@as(f64, -0x4p-1076), @as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-128), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1024), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-972), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-128), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1024), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-972), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x8p-152), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, atan2(@as(f64, 0x4p-1076), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x8p-152), @as(f64, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a3p+0, atan2(@as(f64, -0x4p-1076), @as(f64, -0x1p+0)));

    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x1p+0), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-128), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-1024), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-16384), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x2p-16384), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x8p-972), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x8p-152), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-1076), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x8p-16448), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0xf.fffffp+124), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x1p+0), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x0p+0), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-128), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-1024), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-16384), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x2p-16384), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x8p-972), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x8p-152), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-1076), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x8p-16448), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0xf.fffffp+124), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), std.math.inf(f80)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), std.math.inf(f80)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(std.math.inf(f80), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(-std.math.inf(f80), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x1p+0), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x0p+0), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-128), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1024), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-972), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-152), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1076), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0xf.fffffp+124), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), -std.math.inf(f80)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x1p+0), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x0p+0), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-128), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1024), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-972), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-152), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1076), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0xf.fffffp+124), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), -std.math.inf(f80)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(std.math.inf(f80), std.math.inf(f80)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(-std.math.inf(f80), std.math.inf(f80)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(std.math.inf(f80), -std.math.inf(f80)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(-std.math.inf(f80), -std.math.inf(f80)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x0p+0), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x0p+0), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x0p+0), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x0p+0), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x0p+0), @as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x0p+0), @as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x1p+0), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x1p+0), @as(f80, -0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x1p+0), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x1p+0), @as(f80, -0x0p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffff00000008p-900, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.fffff0000000001p-16260, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.ffffffffffff801p-15364, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fffff00000008p-900, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xf.fffff0000000001p-16260, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xf.ffffffffffff801p-15364, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffp+124), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.ffffffffffff8p+1020), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0xf.fffffffffffffffp+16380), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffp+124), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.ffffffffffff8p+1020), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xf.fffffffffffffffp+16380), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0xa.4bc7d1934f70924p-4, atan2(@as(f80, 0xcp-4), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0xa.4bc7d1934f70924p-4, atan2(@as(f80, -0xcp-4), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x2.7f82ed6f50abffbp+0, atan2(@as(f80, 0xcp-4), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x2.7f82ed6f50abffbp+0, atan2(@as(f80, -0xcp-4), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x1.91ef0ddcd8c8c7bp+0, atan2(@as(f80, 0x6.4p-4), @as(f80, 0x1.30164ap-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2ab449d8p+0, atan2(@as(f80, 0x6.4p-4), @as(f80, 0x1.301648p-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c766p+0, atan2(@as(f80, 0x6.4p-4), @as(f80, 0x1.30164840e171ap-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c76ap+0, atan2(@as(f80, 0x6.4p-4), @as(f80, 0x1.30164840e1719p-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c766p+0, atan2(@as(f80, 0x6.4p-4), @as(f80, 0x1.30164840e1719f8p-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c766p+0, atan2(@as(f80, 0x6.4p-4), @as(f80, 0x1.30164840e1719f7ep-12)));
    try std.testing.expectEqual(0xf.b437a72087797cfp-4, atan2(@as(f80, 0x1.64p+0), @as(f80, 0xe.ep-4)));
    try std.testing.expectEqual(-0x1.cdaa9db2d1b4481p+0, atan2(@as(f80, -0x1.effe8p-8), @as(f80, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9df15fde0914p+0, atan2(@as(f80, -0x1.effe8p-8), @as(f80, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46dd6ap+0, atan2(@as(f80, -0x1.effe8p-8), @as(f80, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46df5ep+0, atan2(@as(f80, -0x1.effe8p-8), @as(f80, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46de96p+0, atan2(@as(f80, -0x1.effe8p-8), @as(f80, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46de96p+0, atan2(@as(f80, -0x1.effe8p-8), @as(f80, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(-0x1.cdaa9d778c4525bep+0, atan2(@as(f80, -0x1.effe82p-8), @as(f80, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db61a6ead0cp+0, atan2(@as(f80, -0x1.effe82p-8), @as(f80, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da6e6d78f68p+0, atan2(@as(f80, -0x1.effe82p-8), @as(f80, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da6e6d7915cp+0, atan2(@as(f80, -0x1.effe82p-8), @as(f80, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da6e6d79094p+0, atan2(@as(f80, -0x1.effe82p-8), @as(f80, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da6e6d79094p+0, atan2(@as(f80, -0x1.effe82p-8), @as(f80, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(-0x1.cdaa9d786fcfa2bp+0, atan2(@as(f80, -0x1.effe81f852716p-8), @as(f80, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db6fdf92adcp+0, atan2(@as(f80, -0x1.effe81f852716p-8), @as(f80, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620d02p+0, atan2(@as(f80, -0x1.effe81f852716p-8), @as(f80, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620ef6p+0, atan2(@as(f80, -0x1.effe81f852716p-8), @as(f80, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620e2cp+0, atan2(@as(f80, -0x1.effe81f852716p-8), @as(f80, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620e2ep+0, atan2(@as(f80, -0x1.effe81f852716p-8), @as(f80, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d6p+0, atan2(@as(f80, -0x1.effe81f852717p-8), @as(f80, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db6fdf92902p+0, atan2(@as(f80, -0x1.effe81f852717p-8), @as(f80, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620b28p+0, atan2(@as(f80, -0x1.effe81f852717p-8), @as(f80, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620d1cp+0, atan2(@as(f80, -0x1.effe81f852717p-8), @as(f80, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c52p+0, atan2(@as(f80, -0x1.effe81f852717p-8), @as(f80, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c54p+0, atan2(@as(f80, -0x1.effe81f852717p-8), @as(f80, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d6p+0, atan2(@as(f80, -0x1.effe81f852716ffcp-8), @as(f80, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db6fdf92902p+0, atan2(@as(f80, -0x1.effe81f852716ffcp-8), @as(f80, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620b28p+0, atan2(@as(f80, -0x1.effe81f852716ffcp-8), @as(f80, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620d1cp+0, atan2(@as(f80, -0x1.effe81f852716ffcp-8), @as(f80, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c54p+0, atan2(@as(f80, -0x1.effe81f852716ffcp-8), @as(f80, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c54p+0, atan2(@as(f80, -0x1.effe81f852716ffcp-8), @as(f80, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d6p+0, atan2(@as(f80, -0x1.effe81f852716ffep-8), @as(f80, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db6fdf92902p+0, atan2(@as(f80, -0x1.effe81f852716ffep-8), @as(f80, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620b28p+0, atan2(@as(f80, -0x1.effe81f852716ffep-8), @as(f80, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620d1cp+0, atan2(@as(f80, -0x1.effe81f852716ffep-8), @as(f80, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c54p+0, atan2(@as(f80, -0x1.effe81f852716ffep-8), @as(f80, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c54p+0, atan2(@as(f80, -0x1.effe81f852716ffep-8), @as(f80, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x1.000002p+0), @as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdba22167c235p-4, atan2(@as(f80, 0x1.000002p+0), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdba22167ba35p-4, atan2(@as(f80, 0x1.000002p+0), @as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fdba22167c1b5p-4, atan2(@as(f80, 0x1.000002p+0), @as(f80, 0x1.00000000000001p+0)));
    try std.testing.expectEqual(0xc.90fd9a22169c235p-4, atan2(@as(f80, 0x1p+0), @as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x1p+0), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168ba35p-4, atan2(@as(f80, 0x1p+0), @as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c1b5p-4, atan2(@as(f80, 0x1p+0), @as(f80, 0x1.00000000000001p+0)));
    try std.testing.expectEqual(0xc.90fd9a22169ca35p-4, atan2(@as(f80, 0x1.0000000000001p+0), @as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168ca35p-4, atan2(@as(f80, 0x1.0000000000001p+0), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x1.0000000000001p+0), @as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c9b5p-4, atan2(@as(f80, 0x1.0000000000001p+0), @as(f80, 0x1.00000000000001p+0)));
    try std.testing.expectEqual(0xc.90fd9a22169c2b5p-4, atan2(@as(f80, 0x1.00000000000001p+0), @as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c2b5p-4, atan2(@as(f80, 0x1.00000000000001p+0), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168bab5p-4, atan2(@as(f80, 0x1.00000000000001p+0), @as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x1.00000000000001p+0), @as(f80, 0x1.00000000000001p+0)));
    try std.testing.expectEqual(0x1.9c22ce44a7229d12p-4, atan2(@as(f80, 0x4.c3841p-4), @as(f80, 0x2.f2f308p+0)));
    try std.testing.expectEqual(-0x1.1dd4c4e2645769d2p-4, atan2(@as(f80, -0xe.cf143p-40), @as(f80, 0xd.3de7ap-36)));
    try std.testing.expectEqual(0x2.7c1782a75e16b744p-4, atan2(@as(f80, 0x5.576cf8p-4), @as(f80, 0x2.21e65p+0)));
    try std.testing.expectEqual(-0x2.1dbac4fa4bfeb75p-4, atan2(@as(f80, -0x4.29411p-4), @as(f80, 0x1.f4755cp+0)));
    try std.testing.expectEqual(-0x1.921fb6b9a118c896p+0, atan2(@as(f80, -0xa.b4101p+20), @as(f80, -0xf.9c4c8p-4)));
    try std.testing.expectEqual(0x9.23e97736442d916p-8, atan2(@as(f80, 0x4.251bb8p-4), @as(f80, 0x7.40ac68p+0)));
    try std.testing.expectEqual(0x1.fd0a44d0aba440f4p-4, atan2(@as(f80, 0x1.47239ep+68), @as(f80, 0xa.3ac3cp+68)));
    try std.testing.expectEqual(-0x1.de89352a0e839634p-4, atan2(@as(f80, -0x6.b0794p-4), @as(f80, 0x3.8ff10cp+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x0p+0), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x8.000008000008p-280, atan2(@as(f80, -0x8p-152), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x7.15e7b61fff36ep-852, atan2(@as(f80, -0x7.15e7af0a1780cp-724), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.e5617957818bbb3ap+0, atan2(@as(f80, 0x3.f16f1p+0), @as(f80, -0x1.546056p+0)));
    try std.testing.expectEqual(-0x1.921fb4fc926936dep+0, atan2(@as(f80, -0x1.9e657cp-24), @as(f80, 0x7.40bb4p-52)));
    try std.testing.expectEqual(-0x3.f2b37e0ca7ac31ap-4, atan2(@as(f80, -0x3.f39e9p+48), @as(f80, 0xf.b02ccp+48)));
    try std.testing.expectEqual(0xc.92e26b8fc63dbafp-4, atan2(@as(f80, 0x6.f2aca8p-56), @as(f80, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e274c80fb7ce9p-4, atan2(@as(f80, 0x6.f2aca8p-56), @as(f80, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e270ff3b35849p-4, atan2(@as(f80, 0x6.f2aca8p-56), @as(f80, 0x6.f107d348a52ep-56)));
    try std.testing.expectEqual(0xc.92e26259ab2ddbfp-4, atan2(@as(f80, 0x6.f2acap-56), @as(f80, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e26b91f4a7f21p-4, atan2(@as(f80, 0x6.f2acap-56), @as(f80, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e267c92025a71p-4, atan2(@as(f80, 0x6.f2acap-56), @as(f80, 0x6.f107d348a52ep-56)));
    try std.testing.expectEqual(0xc.92e26ae105624fbp-4, atan2(@as(f80, 0x6.f2aca7683a51cp-56), @as(f80, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e274194edc638p-4, atan2(@as(f80, 0x6.f2aca7683a51cp-56), @as(f80, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e270507a5a197p-4, atan2(@as(f80, 0x6.f2aca7683a51cp-56), @as(f80, 0x6.f107d348a52ep-56)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x8p-152), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0x3.423123d8009a8p-200, atan2(@as(f80, 0x1.a11891ec004d4p-348), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x1.a11891ec004d4p-348), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xd.334f8c6d4a1b32p-4, atan2(@as(f80, 0x1.a11891ec004d4p-348), @as(f80, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x8p-152), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0x3.706ddacf17c52p-440, atan2(@as(f80, 0x1.b836ed678be29p-588), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x1.b836ed678be29p-588), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xc.932e874124780cbp-4, atan2(@as(f80, 0x1.b836ed678be29p-588), @as(f80, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x8p-152), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0x1.a83f842ef3f73p-484, atan2(@as(f80, 0xd.41fc21779fb98p-636), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0xd.41fc21779fb98p-636), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xc.941e0644c21b3c1p-4, atan2(@as(f80, 0xd.41fc21779fb98p-636), @as(f80, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x8p-152), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x0p+0), @as(f80, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x8p-928, atan2(@as(f80, 0x4p-1076), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x4p-1076), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x2.83beb54428c259p-13168, atan2(@as(f80, 0x1.41df5aa214612c8p-13316), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x1.41df5aa214612c8p-13316), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x5.077d6a885184b2p-12244, atan2(@as(f80, 0x1.41df5aa214612c8p-13316), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edbcp-4, atan2(@as(f80, 0x1.41df5aa214612c8p-13316), @as(f80, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edcp-4, atan2(@as(f80, 0x1.41df5aa214612c8p-13316), @as(f80, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x2.83beb54428c258fcp-13168, atan2(@as(f80, 0x1.41df5aa214612c7ep-13316), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x1.41df5aa214612c7ep-13316), @as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x5.077d6a885184b1f8p-12244, atan2(@as(f80, 0x1.41df5aa214612c7ep-13316), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edb8p-4, atan2(@as(f80, 0x1.41df5aa214612c7ep-13316), @as(f80, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edbcp-4, atan2(@as(f80, 0x1.41df5aa214612c7ep-13316), @as(f80, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x4p-128), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p-896, atan2(@as(f80, 0x4p-1024), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x4p-1024), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-56, atan2(@as(f80, 0x4p-1024), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p-16256, atan2(@as(f80, 0x4p-16384), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1p-15360, atan2(@as(f80, 0x4p-16384), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x4p-16384), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.1b6e192ebbe446c6p+0, atan2(@as(f80, 0x4p-16384), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-15416, atan2(@as(f80, 0x4p-16384), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x8p-16260, atan2(@as(f80, 0x2p-16384), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x8p-15364, atan2(@as(f80, 0x2p-16384), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x7.6b19c1586ed3da28p-4, atan2(@as(f80, 0x2p-16384), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x2p-16384), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x4p-15416, atan2(@as(f80, 0x2p-16384), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x2p-844, atan2(@as(f80, 0x8p-972), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d17c6ap+0, atan2(@as(f80, 0x8p-972), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-972), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-972), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x8p-972), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a300d4p+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x2.08d15159c9bec20cp+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x2.ad8dce72feb5cb3p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18c6ap+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0x4p-128), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-896, atan2(@as(f80, -0x4p-1024), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0x4p-1024), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x8p-56, atan2(@as(f80, -0x4p-1024), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-16256, atan2(@as(f80, -0x4p-16384), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1p-15360, atan2(@as(f80, -0x4p-16384), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0x4p-16384), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.1b6e192ebbe446c6p+0, atan2(@as(f80, -0x4p-16384), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x8p-15416, atan2(@as(f80, -0x4p-16384), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x8p-16260, atan2(@as(f80, -0x2p-16384), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x8p-15364, atan2(@as(f80, -0x2p-16384), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x7.6b19c1586ed3da28p-4, atan2(@as(f80, -0x2p-16384), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0x2p-16384), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x4p-15416, atan2(@as(f80, -0x2p-16384), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x2p-844, atan2(@as(f80, -0x8p-972), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d17c6ap+0, atan2(@as(f80, -0x8p-972), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-972), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-972), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0x8p-972), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a300d4p+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2.08d15159c9bec20cp+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x2.ad8dce72feb5cb3p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18c6ap+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x8p-152), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x8p-928, atan2(@as(f80, 0x4p-1076), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x4p-1076), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1p-16296, atan2(@as(f80, 0x8p-16448), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2p-15372, atan2(@as(f80, 0x8p-16448), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0xc.90fdaa22168c235p-4, atan2(@as(f80, 0x8p-16448), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0x8p-152), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-928, atan2(@as(f80, -0x4p-1076), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0x4p-1076), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1076), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x1p-16296, atan2(@as(f80, -0x8p-16448), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x2p-15372, atan2(@as(f80, -0x8p-16448), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0xc.90fdaa22168c235p-4, atan2(@as(f80, -0x8p-16448), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a46ap+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb34442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x8p-876, atan2(@as(f80, 0x4p-1024), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1746ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x8p-16236, atan2(@as(f80, 0x4p-16384), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1p-15308, atan2(@as(f80, 0x4p-16384), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18468p+0, atan2(@as(f80, 0x4p-16384), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x4p-16236, atan2(@as(f80, 0x2p-16384), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x8p-15312, atan2(@as(f80, 0x2p-16384), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18466p+0, atan2(@as(f80, 0x2p-16384), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1p-820, atan2(@as(f80, 0x8p-972), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-972), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-972), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb74442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1946ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846cp+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ep+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb34442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-876, atan2(@as(f80, -0x4p-1024), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1746ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x8p-16236, atan2(@as(f80, -0x4p-16384), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1p-15308, atan2(@as(f80, -0x4p-16384), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18468p+0, atan2(@as(f80, -0x4p-16384), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x4p-16236, atan2(@as(f80, -0x2p-16384), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-15312, atan2(@as(f80, -0x2p-16384), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18466p+0, atan2(@as(f80, -0x2p-16384), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x1p-820, atan2(@as(f80, -0x8p-972), @as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-972), @as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-972), @as(f80, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb74442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1946ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846cp+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ep+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x1.fffffffffffd5556p-24, atan2(@as(f80, 0x8p-152), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1p-948, atan2(@as(f80, 0x4p-1076), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1p-52, atan2(@as(f80, 0x4p-1076), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-108, atan2(@as(f80, 0x4p-1076), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x2p-16320, atan2(@as(f80, 0x8p-16448), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x2p-15424, atan2(@as(f80, 0x8p-16448), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x2p-64, atan2(@as(f80, 0x8p-16448), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x4p-64, atan2(@as(f80, 0x8p-16448), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1p-15476, atan2(@as(f80, 0x8p-16448), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x3.243f688885a308d4p+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a2f8d4p+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x3.243f6a8885a308dp+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308dp+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x1.fffffffffffd5556p-24, atan2(@as(f80, -0x8p-152), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-948, atan2(@as(f80, -0x4p-1076), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x1p-52, atan2(@as(f80, -0x4p-1076), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1076), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1076), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x8p-108, atan2(@as(f80, -0x4p-1076), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x2p-16320, atan2(@as(f80, -0x8p-16448), @as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x2p-15424, atan2(@as(f80, -0x8p-16448), @as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x2p-64, atan2(@as(f80, -0x8p-16448), @as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x4p-64, atan2(@as(f80, -0x8p-16448), @as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x1p-15476, atan2(@as(f80, -0x8p-16448), @as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x3.243f688885a308d4p+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a2f8d4p+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d1846ap+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x4p-1024)));
    try std.testing.expectEqual(-0x3.243f6a8885a308dp+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x4p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308dp+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x1p+0), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x1p+0), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x1p+0), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x1p+0), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x1p+0), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x1p+0), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-128), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-128), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-128), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-972), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-972), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-972), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-128), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-128), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-128), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-972), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-972), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-972), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-152), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-152), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-152), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-152), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-152), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-152), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.000001000001p-128, atan2(@as(f80, 0x1p+0), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.00000000000008p-1024, atan2(@as(f80, 0x1p+0), @as(f80, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(0x1p-16384, atan2(@as(f80, 0x1p+0), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.000001000001p-128, atan2(@as(f80, -0x1p+0), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.00000000000008p-1024, atan2(@as(f80, -0x1p+0), @as(f80, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(-0x1p-16384, atan2(@as(f80, -0x1p+0), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4.000004000004p-256, atan2(@as(f80, 0x4p-128), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x4.0000000000002p-1152, atan2(@as(f80, 0x4p-128), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-128), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4.000004000004p-1152, atan2(@as(f80, 0x4p-1024), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x4.0000000000002p-2048, atan2(@as(f80, 0x4p-1024), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-1024), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-16384), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-16384), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-16384), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x2p-16384), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x2p-16384), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x2p-16384), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x8.000008000008p-1100, atan2(@as(f80, 0x8p-972), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x8.0000000000004p-1996, atan2(@as(f80, 0x8p-972), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x8p-972), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x4.000004000004p-256, atan2(@as(f80, -0x4p-128), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x4.0000000000002p-1152, atan2(@as(f80, -0x4p-128), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-128), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x4.000004000004p-1152, atan2(@as(f80, -0x4p-1024), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x4.0000000000002p-2048, atan2(@as(f80, -0x4p-1024), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-1024), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-16384), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-16384), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-16384), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x2p-16384), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x2p-16384), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x2p-16384), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x8.000008000008p-1100, atan2(@as(f80, -0x8p-972), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x8.0000000000004p-1996, atan2(@as(f80, -0x8p-972), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x8p-972), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x8.000008000008p-280, atan2(@as(f80, 0x8p-152), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x8.0000000000004p-1176, atan2(@as(f80, 0x8p-152), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x8p-152), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4.000004000004p-1204, atan2(@as(f80, 0x4p-1076), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x4.0000000000002p-2100, atan2(@as(f80, 0x4p-1076), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x4p-1076), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x8p-16448), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x8p-16448), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f80, 0x8p-16448), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x8.000008000008p-280, atan2(@as(f80, -0x8p-152), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x8.0000000000004p-1176, atan2(@as(f80, -0x8p-152), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x8p-152), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x4.000004000004p-1204, atan2(@as(f80, -0x4p-1076), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x4.0000000000002p-2100, atan2(@as(f80, -0x4p-1076), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x4p-1076), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x8p-16448), @as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x8p-16448), @as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f80, -0x8p-16448), @as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4p-128, atan2(@as(f80, 0x4p-128), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x4p-1024, atan2(@as(f80, 0x4p-1024), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x4p-16384, atan2(@as(f80, 0x4p-16384), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x2p-16384, atan2(@as(f80, 0x2p-16384), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x8p-972, atan2(@as(f80, 0x8p-972), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-128, atan2(@as(f80, -0x4p-128), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-1024, atan2(@as(f80, -0x4p-1024), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-16384, atan2(@as(f80, -0x4p-16384), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x2p-16384, atan2(@as(f80, -0x2p-16384), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-972, atan2(@as(f80, -0x8p-972), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x8p-152, atan2(@as(f80, 0x8p-152), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x4p-1076, atan2(@as(f80, 0x4p-1076), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x8p-16448, atan2(@as(f80, 0x8p-16448), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-152, atan2(@as(f80, -0x8p-152), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-1076, atan2(@as(f80, -0x4p-1076), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-16448, atan2(@as(f80, -0x8p-16448), @as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-128), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1024), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-16384), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x2p-16384), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-972), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-128), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1024), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-16384), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x2p-16384), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-972), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-152), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x4p-1076), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, atan2(@as(f80, 0x8p-16448), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-152), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x4p-1076), @as(f80, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d4p+0, atan2(@as(f80, -0x8p-16448), @as(f80, -0x1p+0)));

    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x1p+0), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-128), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-1024), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16384), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x2p-16384), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-972), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-152), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-1076), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-16448), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16448), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16496), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0xf.fffffp+124), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x1p+0), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x0p+0), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-128), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-1024), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16384), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x2p-16384), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-972), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-152), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-1076), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-16448), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16448), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16496), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0xf.fffffp+124), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), std.math.inf(f128)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), std.math.inf(f128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(std.math.inf(f128), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(-std.math.inf(f128), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x1p+0), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x0p+0), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-128), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-152), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.fffffp+124), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), -std.math.inf(f128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x1p+0), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x0p+0), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-128), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-152), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.fffffp+124), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), -std.math.inf(f128)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(std.math.inf(f128), std.math.inf(f128)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(-std.math.inf(f128), std.math.inf(f128)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(std.math.inf(f128), -std.math.inf(f128)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(-std.math.inf(f128), -std.math.inf(f128)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x0p+0), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x0p+0), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x0p+0), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x0p+0), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x0p+0), @as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x0p+0), @as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x1p+0), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x1p+0), @as(f128, -0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x1p+0), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x1p+0), @as(f128, -0x0p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffff00000007fffff80000004p-900, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.fffff0000000000ffffffp-16260, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0xf.fffff00000000000000000000008p-16260, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0xf.fffff00000003fffffc0000005p-900, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.ffffffffffff801p-15364, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0xf.ffffffffffff8000000000000008p-15364, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0xc.90fdaa22168c034c4c6628b80fp-4, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0xf.fffffffffffffffp+16380)));
    // try std.testing.expectEqual(0xc.90fdaa22168c23444c6628b80dc8p-4, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xc.90fdaa22168c23544c6628b80dcp-4, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xc.90fdaa22168c434c4c6628b80c8p-4, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.ffffffffffffc00ffffffffffcp-15364, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0xf.ffffffffffffbffffffffffffc08p-15364, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a489e4e5327a2828p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469ece5327a28294p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469dce5327a28294p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a449e4e5327a282a8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xf.fffff00000007fffff80000004p-900, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xf.fffff0000000000ffffffp-16260, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0xf.fffff00000000000000000000008p-16260, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0xf.fffff00000003fffffc0000005p-900, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xf.ffffffffffff801p-15364, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0xf.ffffffffffff8000000000000008p-15364, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0xc.90fdaa22168c034c4c6628b80fp-4, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0xf.fffffffffffffffp+16380)));
    // try std.testing.expectEqual(-0xc.90fdaa22168c23444c6628b80dc8p-4, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xc.90fdaa22168c23544c6628b80dcp-4, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0xc.90fdaa22168c434c4c6628b80c8p-4, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0xf.ffffffffffffc00ffffffffffcp-15364, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0xf.ffffffffffffbffffffffffffc08p-15364, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a489e4e5327a2828p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469ece5327a28294p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469dce5327a28294p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a449e4e5327a282a8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffp+124), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffff8p+1020), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffp+16380), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffp+124), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffff8p+1020), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffp+16380), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0xa.4bc7d1934f7092419a87f2a457d8p-4, atan2(@as(f128, 0xcp-4), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0xa.4bc7d1934f7092419a87f2a457d8p-4, atan2(@as(f128, -0xcp-4), @as(f128, 0x1p+0)));
    // try std.testing.expectEqual(0x2.7f82ed6f50abffaef9710b03bdf2p+0, atan2(@as(f128, 0xcp-4), @as(f128, -0x1p+0)));
    // try std.testing.expectEqual(-0x2.7f82ed6f50abffaef9710b03bdf2p+0, atan2(@as(f128, -0xcp-4), @as(f128, -0x1p+0)));
    // try std.testing.expectEqual(0x1.91ef0ddcd8c8c7af32bf86fefb5cp+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164ap-12)));
    // try std.testing.expectEqual(0x1.91ef0ddd2ab449d869c37032d24p+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.301648p-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c7667bf2f32a78efp+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164840e171ap-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c7690b4f0474312cp+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164840e1719p-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c766906dd3b4c6b1p+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164840e1719f8p-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c76690bfbf36efe8p+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164840e1719f7ep-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c7669080482ae66cp+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164840e1719f7f8ca8198f1d3fp-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c7669080482ae66cp+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164840e1719f7f8ca8198f1d3ep-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c7669080482ae66cp+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164840e1719f7f8ca8198f1d8p-12)));
    try std.testing.expectEqual(0x1.91ef0ddd2052c7669080482ae66cp+0, atan2(@as(f128, 0x6.4p-4), @as(f128, 0x1.30164840e1719f7f8ca8198f1dp-12)));
    // try std.testing.expectEqual(0xf.b437a72087797cf2cd8d81111208p-4, atan2(@as(f128, 0x1.64p+0), @as(f128, 0xe.ep-4)));
    try std.testing.expectEqual(-0x1.cdaa9db2d1b4480f9a874bb84e72p+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9df15fde0913df719335d86fp+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1ep-12)));
    // try std.testing.expectEqual(-0x1.cdaa9de22c46dd6a08ed2c333239p+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9de22c46df5e7a3b2571f0a5p+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1de0e51248p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9de22c46de9520a4c629b208p+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46de955f32efe8d9dfp+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46de9526e102fad22cp+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46de9526e102fad22cp+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9de22c46de9526e102fad22bp+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    try std.testing.expectEqual(-0x1.cdaa9de22c46de9526e102fad23ap+0, atan2(@as(f128, -0x1.effe8p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    // try std.testing.expectEqual(-0x1.cdaa9d778c4525be29e0416fa252p+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db61a6ead0c075aa78c1d4dp+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1ep-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da6e6d78f68821ef0a0a971p+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da6e6d7915cf36b1c2c2c89p+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1de0e51248p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da6e6d7909399d576a70acbp+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da6e6d79093d863a02c7c3bp+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da6e6d79093a011b3726a46p+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da6e6d79093a011b3726a46p+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da6e6d79093a011b3726a44p+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da6e6d79093a011b3726a54p+0, atan2(@as(f128, -0x1.effe82p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    // try std.testing.expectEqual(-0x1.cdaa9d786fcfa2afd2e41f4d8586p+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1d8p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9db6fdf92adb3efd640a23eep+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1ep-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620d01e29ac2cea0eep+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620ef653e6f54698fcp+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1de0e51248p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620e2cfa514cf8542fp+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620e2d38df767ea32ep+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620e2d008d89c3c9bfp+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620e2d008d89c3c9cp+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620e2d008d89c3c9bep+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620e2d008d89c3c9cdp+0, atan2(@as(f128, -0x1.effe81f852716p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d5a76cce52e85bp+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db6fdf929011384455c4cc5p+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1ep-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620b27b722145353aep+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1de0e51244p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620d1c286e46cb4badp+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1de0e51248p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c52ced89e7d06e6p+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c530d66c80355e5p+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c52d514db487c76p+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c52d514db487c77p+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c52d514db487c75p+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c52d514db487c84p+0, atan2(@as(f128, -0x1.effe81f852717p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    // try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d61df7ac272703p+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1d8p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9db6fdf929018a0f23a3f83bp+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1ep-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620b282dacf27ef281p+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620d1c9ef924f6ea81p+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1de0e51248p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5345637ca8a5b9p+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5383f1a62ef4b8p+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c534b9fb9741b4ap+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c534b9fb9741b4ap+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c534b9fb9741b48p+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c534b9fb9741b58p+0, atan2(@as(f128, -0x1.effe81f852716ffcp-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    // try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d5e2b23d3d07afp+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1d8p-12)));
    try std.testing.expectEqual(-0x1.cdaa9db6fdf929014ec9b480228p+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620b27f26783692317p+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620d1c63b3b5e11b17p+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c530a1e0d92d64fp+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c5348ac3719254ep+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c53105a4a5e4bep+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c53105a4a5e4bep+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c53105a4a5e4bdep+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c53105a4a5e4beep+0, atan2(@as(f128, -0x1.effe81f852716ffep-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    // try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d61c33dac76549p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1d8p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9db6fdf92901884b52427e92p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620b282be9211de3cp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620d1c9d355395dbcp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1de0e51248p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c53439fab4796f8p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c53822dd4cde5f7p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c89p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c89p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c87p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c97p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac291p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    // try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d61c33dac76549p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1d8p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9db6fdf92901884b52427e92p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1ep-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620b282be9211de3cp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620d1c9d355395dbcp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1de0e51248p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c53439fab4796f8p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c53822dd4cde5f7p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c89p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c89p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c87p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c96p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac292p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    // try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d61c33dac7654bp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1d8p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9db6fdf92901884b52427e94p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1ep-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620b282be9211de3c2p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620d1c9d355395dbc2p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1de0e51248p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c53439fab4796fap+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c53822dd4cde5f9p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c8bp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c8bp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c89p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c99p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac28p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    // try std.testing.expectEqual(-0x1.cdaa9d786fcfa0d61c33dac7653cp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1d8p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9db6fdf92901884b52427e85p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1ep-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620b282be9211de3b3p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1de0e51244p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620d1c9d355395dbb3p+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1de0e51248p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c53439fab4796ebp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1de0e5124664p-12)));
    try std.testing.expectEqual(-0x1.cdaa9da7ca620c53822dd4cde5eap+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1de0e51246648p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c7cp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca483cp-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c7cp+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca484p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c7ap+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca48p-12)));
    // try std.testing.expectEqual(-0x1.cdaa9da7ca620c5349dbe8130c8ap+0, atan2(@as(f128, -0x1.effe81f852716ffc0f3eeb1ac3p-8), @as(f128, -0x7.57d1de0e51246640cc2340ca4ap-12)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x1.000002p+0), @as(f128, 0x1.000002p+0)));
    // try std.testing.expectEqual(0xc.90fdba22167c234c5710d362b87p-4, atan2(@as(f128, 0x1.000002p+0), @as(f128, 0x1p+0)));
    // try std.testing.expectEqual(0xc.90fdba22167ba34c5710d363bc7p-4, atan2(@as(f128, 0x1.000002p+0), @as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fdba22167c1b4c5710d362c87p-4, atan2(@as(f128, 0x1.000002p+0), @as(f128, 0x1.00000000000001p+0)));
    try std.testing.expectEqual(0xc.90fd9a22169c234c41bb7e0d6318p-4, atan2(@as(f128, 0x1p+0), @as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x1p+0), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168ba34c4c6628b811cp-4, atan2(@as(f128, 0x1p+0), @as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c1b4c4c6628b80dc8p-4, atan2(@as(f128, 0x1p+0), @as(f128, 0x1.00000000000001p+0)));
    try std.testing.expectEqual(0xc.90fd9a22169ca34c41bb7e0c5f18p-4, atan2(@as(f128, 0x1.0000000000001p+0), @as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168ca34c4c6628b809cp-4, atan2(@as(f128, 0x1.0000000000001p+0), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x1.0000000000001p+0), @as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c9b4c4c6628b809c8p-4, atan2(@as(f128, 0x1.0000000000001p+0), @as(f128, 0x1.00000000000001p+0)));
    try std.testing.expectEqual(0xc.90fd9a22169c2b4c41bb7e0d531p-4, atan2(@as(f128, 0x1.00000000000001p+0), @as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c2b4c4c6628b80dcp-4, atan2(@as(f128, 0x1.00000000000001p+0), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168bab4c4c6628b811cp-4, atan2(@as(f128, 0x1.00000000000001p+0), @as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x1.00000000000001p+0), @as(f128, 0x1.00000000000001p+0)));
    // try std.testing.expectEqual(0x1.9c22ce44a7229d114c2b882266fap-4, atan2(@as(f128, 0x4.c3841p-4), @as(f128, 0x2.f2f308p+0)));
    try std.testing.expectEqual(-0x1.1dd4c4e2645769d1f7ebdc32a451p-4, atan2(@as(f128, -0xe.cf143p-40), @as(f128, 0xd.3de7ap-36)));
    try std.testing.expectEqual(0x2.7c1782a75e16b743e48c247c62cap-4, atan2(@as(f128, 0x5.576cf8p-4), @as(f128, 0x2.21e65p+0)));
    try std.testing.expectEqual(-0x2.1dbac4fa4bfeb74f6140009955a6p-4, atan2(@as(f128, -0x4.29411p-4), @as(f128, 0x1.f4755cp+0)));
    try std.testing.expectEqual(-0x1.921fb6b9a118c89590d474178551p+0, atan2(@as(f128, -0xa.b4101p+20), @as(f128, -0xf.9c4c8p-4)));
    // try std.testing.expectEqual(0x9.23e97736442d915917b21858b148p-8, atan2(@as(f128, 0x4.251bb8p-4), @as(f128, 0x7.40ac68p+0)));
    // try std.testing.expectEqual(0x1.fd0a44d0aba440f30193e8545bc2p-4, atan2(@as(f128, 0x1.47239ep+68), @as(f128, 0xa.3ac3cp+68)));
    try std.testing.expectEqual(-0x1.de89352a0e839633c32d65e25422p-4, atan2(@as(f128, -0x6.b0794p-4), @as(f128, 0x3.8ff10cp+0)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x0p+0), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x8.000008000008000008000008p-280, atan2(@as(f128, -0x8p-152), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x7.15e7b61fff36dfff36dfff36ep-852, atan2(@as(f128, -0x7.15e7af0a1780cp-724), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.e5617957818bbb3ab867fdf781dp+0, atan2(@as(f128, 0x3.f16f1p+0), @as(f128, -0x1.546056p+0)));
    try std.testing.expectEqual(-0x1.921fb4fc926936de7a5c2816052ep+0, atan2(@as(f128, -0x1.9e657cp-24), @as(f128, 0x7.40bb4p-52)));
    // try std.testing.expectEqual(-0x3.f2b37e0ca7ac31a1c0615ac92a5cp-4, atan2(@as(f128, -0x3.f39e9p+48), @as(f128, 0xf.b02ccp+48)));
    try std.testing.expectEqual(0xc.92e26b8fc63dbaedf472f75dee78p-4, atan2(@as(f128, 0x6.f2aca8p-56), @as(f128, 0x6.f107d8p-56)));
    // try std.testing.expectEqual(0xc.92e274c80fb7ce899797f0d7d6b8p-4, atan2(@as(f128, 0x6.f2aca8p-56), @as(f128, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e270ff3b358490f84e1323b17p-4, atan2(@as(f128, 0x6.f2aca8p-56), @as(f128, 0x6.f107d348a52ep-56)));
    // try std.testing.expectEqual(0xc.92e26259ab2ddbf0371a8fb67b7p-4, atan2(@as(f128, 0x6.f2acap-56), @as(f128, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e26b91f4a7f20f1e8cafdcb858p-4, atan2(@as(f128, 0x6.f2acap-56), @as(f128, 0x6.f107dp-56)));
    try std.testing.expectEqual(0xc.92e267c92025a70e72e20b0139ap-4, atan2(@as(f128, 0x6.f2acap-56), @as(f128, 0x6.f107d348a52ep-56)));
    try std.testing.expectEqual(0xc.92e26ae105624fb58b68c52a4fe8p-4, atan2(@as(f128, 0x6.f2aca7683a51cp-56), @as(f128, 0x6.f107d8p-56)));
    try std.testing.expectEqual(0xc.92e274194edc6380dab95d83d958p-4, atan2(@as(f128, 0x6.f2aca7683a51cp-56), @as(f128, 0x6.f107dp-56)));
    // try std.testing.expectEqual(0xc.92e270507a5a1974a9dd2e786798p-4, atan2(@as(f128, 0x6.f2aca7683a51cp-56), @as(f128, 0x6.f107d348a52ep-56)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x8p-152), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0x3.423123d8009a8p-200, atan2(@as(f128, 0x1.a11891ec004d4p-348), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x1.a11891ec004d4p-348), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0xd.334f8c6d4a1b32001a10fa8adce8p-4, atan2(@as(f128, 0x1.a11891ec004d4p-348), @as(f128, 0x1.814830510be26p-348)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x8p-152), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0x3.706ddacf17c52p-440, atan2(@as(f128, 0x1.b836ed678be29p-588), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x1.b836ed678be29p-588), @as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0xc.932e874124780ca850c92ebb72cp-4, atan2(@as(f128, 0x1.b836ed678be29p-588), @as(f128, 0x1.b7be6f5a03a8cp-588)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x8p-152), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0x1.a83f842ef3f73p-484, atan2(@as(f128, 0xd.41fc21779fb98p-636), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0xd.41fc21779fb98p-636), @as(f128, 0x0p+0)));
    // try std.testing.expectEqual(0xc.941e0644c21b3c15d3dbbf4b3bap-4, atan2(@as(f128, 0xd.41fc21779fb98p-636), @as(f128, 0xd.3ccec5333bfp-636)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x8p-152), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x5.e53b26a270a29eb9f77ef8ef7af8p-13316)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x0p+0), @as(f128, 0x5.e53b26a270a29eb9f77ef8ef7af8p-13316)));
    try std.testing.expectEqual(0x8p-928, atan2(@as(f128, 0x4p-1076), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x4p-1076), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x5.e53b26a270a29eb9f77ef8ef7af8p-13316)));
    try std.testing.expectEqual(0x2.83beb54428c259p-13168, atan2(@as(f128, 0x1.41df5aa214612c8p-13316), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x1.41df5aa214612c8p-13316), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x5.077d6a885184b2p-12244, atan2(@as(f128, 0x1.41df5aa214612c8p-13316), @as(f128, 0x4p-1076)));
    // try std.testing.expectEqual(0x3.5ca813fb6ec1edbb4ddcc1a43934p-4, atan2(@as(f128, 0x1.41df5aa214612c8p-13316), @as(f128, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edbfbbaf43070b86p-4, atan2(@as(f128, 0x1.41df5aa214612c8p-13316), @as(f128, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edbea4efe3fb5d9ap-4, atan2(@as(f128, 0x1.41df5aa214612c8p-13316), @as(f128, 0x5.e53b26a270a29eb9f77ef8ef7af8p-13316)));
    try std.testing.expectEqual(0x2.83beb54428c258fcp-13168, atan2(@as(f128, 0x1.41df5aa214612c7ep-13316), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x1.41df5aa214612c7ep-13316), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x5.077d6a885184b1f8p-12244, atan2(@as(f128, 0x1.41df5aa214612c7ep-13316), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edb61cc3c20a3d86p-4, atan2(@as(f128, 0x1.41df5aa214612c7ep-13316), @as(f128, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edba8a96436d0fd8p-4, atan2(@as(f128, 0x1.41df5aa214612c7ep-13316), @as(f128, 0x5.e53b26a270a29eb8p-13316)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edb973d6e46161ecp-4, atan2(@as(f128, 0x1.41df5aa214612c7ep-13316), @as(f128, 0x5.e53b26a270a29eb9f77ef8ef7af8p-13316)));
    try std.testing.expectEqual(0x2.83beb54428c258fc033f4d5bd1p-13168, atan2(@as(f128, 0x1.41df5aa214612c7e019fa6ade88p-13316), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x1.41df5aa214612c7e019fa6ade88p-13316), @as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x5.077d6a885184b1f8067e9ab7a2p-12244, atan2(@as(f128, 0x1.41df5aa214612c7e019fa6ade88p-13316), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edb620fabe7bf832p-4, atan2(@as(f128, 0x1.41df5aa214612c7e019fa6ade88p-13316), @as(f128, 0x5.e53b26a270a29ecp-13316)));
    try std.testing.expectEqual(0x3.5ca813fb6ec1edba8ecd3fdeca84p-4, atan2(@as(f128, 0x1.41df5aa214612c7e019fa6ade88p-13316), @as(f128, 0x5.e53b26a270a29eb8p-13316)));
    // try std.testing.expectEqual(0x3.5ca813fb6ec1edb9780de0d31c98p-4, atan2(@as(f128, 0x1.41df5aa214612c7e019fa6ade88p-13316), @as(f128, 0x5.e53b26a270a29eb9f77ef8ef7af8p-13316)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x4p-128), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p-896, atan2(@as(f128, 0x4p-1024), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x4p-1024), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffffff54p-56, atan2(@as(f128, 0x4p-1024), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p-16256, atan2(@as(f128, 0x4p-16384), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p-15360, atan2(@as(f128, 0x4p-16384), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x4p-16384), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.1b6e192ebbe446c6d19aa220a39bp+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-15416, atan2(@as(f128, 0x4p-16384), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x8p-16260, atan2(@as(f128, 0x2p-16384), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x8p-15364, atan2(@as(f128, 0x2p-16384), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x7.6b19c1586ed3da2b7f222f65e1d4p-4, atan2(@as(f128, 0x2p-16384), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x2p-16384), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x4p-15416, atan2(@as(f128, 0x2p-16384), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x2p-844, atan2(@as(f128, 0x8p-972), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d17c69898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x8p-972), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a300d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x4p-16384)));
    // try std.testing.expectEqual(0x2.08d15159c9bec20c417ee80d5fd6p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x2.ad8dce72feb5cb305b276737a554p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18c69898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x4p-128), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-896, atan2(@as(f128, -0x4p-1024), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x4p-1024), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x7.ffffffffffffffffffffffffff54p-56, atan2(@as(f128, -0x4p-1024), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-16256, atan2(@as(f128, -0x4p-16384), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1p-15360, atan2(@as(f128, -0x4p-16384), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x4p-16384), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.1b6e192ebbe446c6d19aa220a39bp+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x8p-15416, atan2(@as(f128, -0x4p-16384), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x8p-16260, atan2(@as(f128, -0x2p-16384), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x8p-15364, atan2(@as(f128, -0x2p-16384), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x7.6b19c1586ed3da2b7f222f65e1d4p-4, atan2(@as(f128, -0x2p-16384), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x2p-16384), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x4p-15416, atan2(@as(f128, -0x2p-16384), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x2p-844, atan2(@as(f128, -0x8p-972), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d17c69898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x8p-972), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a300d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x4p-16384)));
    // try std.testing.expectEqual(-0x2.08d15159c9bec20c417ee80d5fd6p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x2.ad8dce72feb5cb305b276737a554p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18c69898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x8p-152), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x8p-928, atan2(@as(f128, 0x4p-1076), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x4p-1076), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1p-16296, atan2(@as(f128, 0x8p-16448), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2p-15372, atan2(@as(f128, 0x8p-16448), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x8p-16448), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.1b6e192ebbe446c6d19aa220a39bp+0, atan2(@as(f128, 0x8p-16448), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d10469898cc51701b8p+0, atan2(@as(f128, 0x8p-16448), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x8p-16300, atan2(@as(f128, 0x4p-16448), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p-15372, atan2(@as(f128, 0x4p-16448), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x7.6b19c1586ed3da2b7f222f65e1d4p-4, atan2(@as(f128, 0x4p-16448), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x4p-16448), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d08469898cc51701b8p+0, atan2(@as(f128, 0x4p-16448), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x8p-16348, atan2(@as(f128, 0x4p-16496), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p-15420, atan2(@as(f128, 0x4p-16496), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x7.ffffffffffffffffffffffff5554p-52, atan2(@as(f128, 0x4p-16496), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffaaaa8p-52, atan2(@as(f128, 0x4p-16496), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, 0x4p-16496), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x8p-16448)));
    // try std.testing.expectEqual(0x2.08d15159c9bec20c417ee80d5fd6p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d20469898cc51701b8p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x2.ad8dce72feb5cb305b276737a554p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d28469898cc51701b8p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x3.243f6a8885a288d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x3.243f6a8885a208d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x8p-152), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-928, atan2(@as(f128, -0x4p-1076), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x4p-1076), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x1p-16296, atan2(@as(f128, -0x8p-16448), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x2p-15372, atan2(@as(f128, -0x8p-16448), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x8p-16448), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.1b6e192ebbe446c6d19aa220a39bp+0, atan2(@as(f128, -0x8p-16448), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d10469898cc51701b8p+0, atan2(@as(f128, -0x8p-16448), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-16300, atan2(@as(f128, -0x4p-16448), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1p-15372, atan2(@as(f128, -0x4p-16448), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x7.6b19c1586ed3da2b7f222f65e1d4p-4, atan2(@as(f128, -0x4p-16448), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x4p-16448), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d08469898cc51701b8p+0, atan2(@as(f128, -0x4p-16448), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-16348, atan2(@as(f128, -0x4p-16496), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1p-15420, atan2(@as(f128, -0x4p-16496), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x7.ffffffffffffffffffffffff5554p-52, atan2(@as(f128, -0x4p-16496), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffaaaa8p-52, atan2(@as(f128, -0x4p-16496), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0xc.90fdaa22168c234c4c6628b80dcp-4, atan2(@as(f128, -0x4p-16496), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x8p-16448)));
    // try std.testing.expectEqual(-0x2.08d15159c9bec20c417ee80d5fd6p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d20469898cc51701b8p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x2.ad8dce72feb5cb305b276737a554p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d28469898cc51701b8p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x3.243f6a8885a288d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x3.243f6a8885a208d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x2.5b2f8fe6643a469e4e5327a28294p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb34442d184698c376fc1ac63p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x8p-876, atan2(@as(f128, 0x4p-1024), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d17469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x8p-16236, atan2(@as(f128, 0x4p-16384), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1p-15308, atan2(@as(f128, 0x4p-16384), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18467898cc51701b8p+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18468898cc51701b8p+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b7p+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x4p-16236, atan2(@as(f128, 0x2p-16384), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x8p-15312, atan2(@as(f128, 0x2p-16384), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18465898cc51701b8p+0, atan2(@as(f128, 0x2p-16384), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18467898cc51701b8p+0, atan2(@as(f128, 0x2p-16384), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b6p+0, atan2(@as(f128, 0x2p-16384), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1p-820, atan2(@as(f128, 0x8p-972), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc5170138p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0x4p-16496)));
    // try std.testing.expectEqual(0x1.921fb74442d1846986e21a6c570ep+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d19469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846b898cc51701b8p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846a898cc51701b8p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b9p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846d898cc51701b8p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846b898cc51701b8p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701bap+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc5170238p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x1.921fb34442d184698c376fc1ac63p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-876, atan2(@as(f128, -0x4p-1024), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d17469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x8p-16236, atan2(@as(f128, -0x4p-16384), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1p-15308, atan2(@as(f128, -0x4p-16384), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18467898cc51701b8p+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18468898cc51701b8p+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b7p+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x4p-16236, atan2(@as(f128, -0x2p-16384), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x8p-15312, atan2(@as(f128, -0x2p-16384), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18465898cc51701b8p+0, atan2(@as(f128, -0x2p-16384), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18467898cc51701b8p+0, atan2(@as(f128, -0x2p-16384), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b6p+0, atan2(@as(f128, -0x2p-16384), @as(f128, 0x4p-16496)));
    try std.testing.expectEqual(-0x1p-820, atan2(@as(f128, -0x8p-972), @as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc5170138p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0x4p-16496)));
    // try std.testing.expectEqual(-0x1.921fb74442d1846986e21a6c570ep+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d19469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846b898cc51701b8p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846a898cc51701b8p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b9p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d1846d898cc51701b8p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d1846b898cc51701b8p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701bap+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x8p-152)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc5170238p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x4p-1076)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x8p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x4p-16448)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x4p-16496)));
    try std.testing.expectEqual(0x1.fffffffffffd55555555555bbbbcp-24, atan2(@as(f128, 0x8p-152), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p-948, atan2(@as(f128, 0x4p-1076), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0xf.fffffffffffffffffffffffffaa8p-56, atan2(@as(f128, 0x4p-1076), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-108, atan2(@as(f128, 0x4p-1076), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x2p-16320, atan2(@as(f128, 0x8p-16448), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x2p-15424, atan2(@as(f128, 0x8p-16448), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x2p-64, atan2(@as(f128, 0x8p-16448), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x4p-64, atan2(@as(f128, 0x8p-16448), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1p-15476, atan2(@as(f128, 0x8p-16448), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p-16320, atan2(@as(f128, 0x4p-16448), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p-15424, atan2(@as(f128, 0x4p-16448), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1p-64, atan2(@as(f128, 0x4p-16448), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-64, atan2(@as(f128, 0x4p-16448), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-15480, atan2(@as(f128, 0x4p-16448), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1p-16368, atan2(@as(f128, 0x4p-16496), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1p-15472, atan2(@as(f128, 0x4p-16496), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1p-112, atan2(@as(f128, 0x4p-16496), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x2p-112, atan2(@as(f128, 0x4p-16496), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x8p-15528, atan2(@as(f128, 0x4p-16496), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x3.243f688885a308d315c434d8ae1cp+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a2f8d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e02fp+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x3.243f6a8885a308d113198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308cf13198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x3.243f6a8885a308d213198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d113198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e036ep+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x1.fffffffffffd55555555555bbbbcp-24, atan2(@as(f128, -0x8p-152), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-948, atan2(@as(f128, -0x4p-1076), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0xf.fffffffffffffffffffffffffaa8p-56, atan2(@as(f128, -0x4p-1076), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x8p-108, atan2(@as(f128, -0x4p-1076), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x2p-16320, atan2(@as(f128, -0x8p-16448), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x2p-15424, atan2(@as(f128, -0x8p-16448), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x2p-64, atan2(@as(f128, -0x8p-16448), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x4p-64, atan2(@as(f128, -0x8p-16448), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x1p-15476, atan2(@as(f128, -0x8p-16448), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-16320, atan2(@as(f128, -0x4p-16448), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1p-15424, atan2(@as(f128, -0x4p-16448), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1p-64, atan2(@as(f128, -0x4p-16448), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x2p-64, atan2(@as(f128, -0x4p-16448), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x8p-15480, atan2(@as(f128, -0x4p-16448), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x1p-16368, atan2(@as(f128, -0x4p-16496), @as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x1p-15472, atan2(@as(f128, -0x4p-16496), @as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x1p-112, atan2(@as(f128, -0x4p-16496), @as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x2p-112, atan2(@as(f128, -0x4p-16496), @as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x8p-15528, atan2(@as(f128, -0x4p-16496), @as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x3.243f688885a308d315c434d8ae1cp+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a2f8d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x1.921fb54442d18469898cc51701b8p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e02fp+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d113198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308cf13198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d213198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d113198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x4p-128)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x4p-1024)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x4p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e036ep+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x2p-16384)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x1p+0), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x1p+0), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x1p+0), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x1p+0), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x1p+0), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x1p+0), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x1p+0), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x1p+0), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x1p+0), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x1p+0), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x1.000001000001000001000001p-128, atan2(@as(f128, 0x1p+0), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.000000000000080000000000004p-1024, atan2(@as(f128, 0x1p+0), @as(f128, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(0x1.0000000000000001p-16384, atan2(@as(f128, 0x1p+0), @as(f128, 0xf.fffffffffffffffp+16380)));
    // try std.testing.expectEqual(0x1p-16384, atan2(@as(f128, 0x1p+0), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x1.000000000000040000000000005p-1024, atan2(@as(f128, 0x1p+0), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x1.000001000001000001000001p-128, atan2(@as(f128, -0x1p+0), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x1.000000000000080000000000004p-1024, atan2(@as(f128, -0x1p+0), @as(f128, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(-0x1.0000000000000001p-16384, atan2(@as(f128, -0x1p+0), @as(f128, 0xf.fffffffffffffffp+16380)));
    // try std.testing.expectEqual(-0x1p-16384, atan2(@as(f128, -0x1p+0), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x1.000000000000040000000000005p-1024, atan2(@as(f128, -0x1p+0), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4.000004000004000004000004p-256, atan2(@as(f128, 0x4p-128), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x4.00000000000020000000000001p-1152, atan2(@as(f128, 0x4p-128), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-128), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x4.000000000000100000000000014p-1152, atan2(@as(f128, 0x4p-128), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4.000004000004000004000004p-1152, atan2(@as(f128, 0x4p-1024), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x4.00000000000020000000000001p-2048, atan2(@as(f128, 0x4p-1024), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-1024), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-1024), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x4.000000000000100000000000014p-2048, atan2(@as(f128, 0x4p-1024), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16384), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x2p-16384), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x2p-16384), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x2p-16384), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x2p-16384), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x2p-16384), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x8.000008000008000008000008p-1100, atan2(@as(f128, 0x8p-972), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x8.00000000000040000000000002p-1996, atan2(@as(f128, 0x8p-972), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-972), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x8.000000000000200000000000028p-1996, atan2(@as(f128, 0x8p-972), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x4.000004000004000004000004p-256, atan2(@as(f128, -0x4p-128), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x4.00000000000020000000000001p-1152, atan2(@as(f128, -0x4p-128), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-128), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x4.000000000000100000000000014p-1152, atan2(@as(f128, -0x4p-128), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x4.000004000004000004000004p-1152, atan2(@as(f128, -0x4p-1024), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x4.00000000000020000000000001p-2048, atan2(@as(f128, -0x4p-1024), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-1024), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-1024), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x4.000000000000100000000000014p-2048, atan2(@as(f128, -0x4p-1024), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16384), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x2p-16384), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x2p-16384), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x2p-16384), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x2p-16384), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x2p-16384), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x8.000008000008000008000008p-1100, atan2(@as(f128, -0x8p-972), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x8.00000000000040000000000002p-1996, atan2(@as(f128, -0x8p-972), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-972), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x8.000000000000200000000000028p-1996, atan2(@as(f128, -0x8p-972), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x8.000008000008000008000008p-280, atan2(@as(f128, 0x8p-152), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x8.00000000000040000000000002p-1176, atan2(@as(f128, 0x8p-152), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-152), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x8.000000000000200000000000028p-1176, atan2(@as(f128, 0x8p-152), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4.000004000004000004000004p-1204, atan2(@as(f128, 0x4p-1076), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x4.00000000000020000000000001p-2100, atan2(@as(f128, 0x4p-1076), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-1076), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x4.000000000000100000000000014p-2100, atan2(@as(f128, 0x4p-1076), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-16448), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-16448), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-16448), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-16448), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x8p-16448), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16448), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16448), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16448), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16448), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16448), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16496), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16496), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16496), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16496), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x0p+0, atan2(@as(f128, 0x4p-16496), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x8.000008000008000008000008p-280, atan2(@as(f128, -0x8p-152), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x8.00000000000040000000000002p-1176, atan2(@as(f128, -0x8p-152), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-152), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x8.000000000000200000000000028p-1176, atan2(@as(f128, -0x8p-152), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x4.000004000004000004000004p-1204, atan2(@as(f128, -0x4p-1076), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x4.00000000000020000000000001p-2100, atan2(@as(f128, -0x4p-1076), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-1076), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-1076), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x4.000000000000100000000000014p-2100, atan2(@as(f128, -0x4p-1076), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-16448), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-16448), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-16448), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-16448), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x8p-16448), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16448), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16448), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16448), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16448), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16448), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16496), @as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16496), @as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16496), @as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16496), @as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(-0x0p+0, atan2(@as(f128, -0x4p-16496), @as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x4p-128, atan2(@as(f128, 0x4p-128), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x4p-1024, atan2(@as(f128, 0x4p-1024), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x4p-16384, atan2(@as(f128, 0x4p-16384), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x2p-16384, atan2(@as(f128, 0x2p-16384), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x8p-972, atan2(@as(f128, 0x8p-972), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-128, atan2(@as(f128, -0x4p-128), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-1024, atan2(@as(f128, -0x4p-1024), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-16384, atan2(@as(f128, -0x4p-16384), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x2p-16384, atan2(@as(f128, -0x2p-16384), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-972, atan2(@as(f128, -0x8p-972), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x8p-152, atan2(@as(f128, 0x8p-152), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x4p-1076, atan2(@as(f128, 0x4p-1076), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x8p-16448, atan2(@as(f128, 0x8p-16448), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x4p-16448, atan2(@as(f128, 0x4p-16448), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x4p-16496, atan2(@as(f128, 0x4p-16496), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-152, atan2(@as(f128, -0x8p-152), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-1076, atan2(@as(f128, -0x4p-1076), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x8p-16448, atan2(@as(f128, -0x8p-16448), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-16448, atan2(@as(f128, -0x4p-16448), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x4p-16496, atan2(@as(f128, -0x4p-16496), @as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-128), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1024), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16384), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x2p-16384), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-972), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-128), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1024), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16384), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x2p-16384), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-972), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-152), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-1076), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x8p-16448), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16448), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, 0x4p-16496), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-152), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-1076), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x8p-16448), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16448), @as(f128, -0x1p+0)));
    try std.testing.expectEqual(-0x3.243f6a8885a308d313198a2e037p+0, atan2(@as(f128, -0x4p-16496), @as(f128, -0x1p+0)));
}
