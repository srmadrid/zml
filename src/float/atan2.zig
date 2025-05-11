const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const dla = @import("dla.zig");
const atnat2 = @import("atnat2.zig");
const ldbl128 = @import("ldbl128.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn atan2(y: anytype, x: anytype) EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))) {
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
                        return atan2_32(cast(f32, y, .{}), cast(f32, x, .{}));
                    },
                    f64 => {
                        // glibc/sysdeps/ieee754/dbl-64/e_atan2.c
                        return atan2_64(cast(f64, y, .{}), cast(f64, x, .{}));
                    },
                    f80 => return cast(f80, atan2_128(cast(f128, y, .{}), cast(f128, x, .{})), .{}),
                    f128 => {
                        // glibc/sysdeps/ieee754/ldbl-128/e_atan2l.c
                        return atan2_128(cast(f128, y, .{}), cast(f128, x, .{}));
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
            const ret: f64 = float.copysign(z, y);
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
                return float.copysign(z, y);
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
            return float.copysign(z, y);
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
            return float.copysign(z, y);
        }

        var i: i32 = cast(i32, (TWO52 + 256 * u) - TWO52, .{});
        i -= 16;
        v = (u - @as(f64, @bitCast(atnat2.cij[@intCast(i)][0]))) + du;

        const zz: f64 = @as(f64, @bitCast(atnat2.hpi1)) - v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][2])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][3])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][4])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][5])) + v * @as(f64, @bitCast(atnat2.cij[@intCast(i)][6]))))));
        const t1: f64 = @as(f64, @bitCast(atnat2.hpi)) - @as(f64, @bitCast(atnat2.cij[@intCast(i)][1]));
        const z: f64 = t1 + zz;
        // Max ULP is 0.503.
        return float.copysign(z, y);
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
            return float.copysign(z, y);
        }

        var i: i32 = cast(i32, (TWO52 + 256 * u) - TWO52, .{});
        i -= 16;
        v = (u - @as(f64, @bitCast(atnat2.cij[@intCast(i)][0]))) + du;
        const zz: f64 = @as(f64, @bitCast(atnat2.hpi1)) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][2])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][3])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][4])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][5])) + v * @as(f64, @bitCast(atnat2.cij[@intCast(i)][6]))))));
        const t1: f64 = @as(f64, @bitCast(atnat2.hpi)) + @as(f64, @bitCast(atnat2.cij[@intCast(i)][1]));
        const z: f64 = t1 + zz;
        // Max ULP is 0.503.
        return float.copysign(z, y);
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
        return float.copysign(z, y);
    }

    var i: i32 = cast(i32, (TWO52 + 256 * u) - TWO52, .{});
    i -= 16;
    v = (u - @as(f64, @bitCast(atnat2.cij[@intCast(i)][0]))) + du;
    const zz: f64 = @as(f64, @bitCast(atnat2.opi1)) - v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][2])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][3])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][4])) + v * (@as(f64, @bitCast(atnat2.cij[@intCast(i)][5])) + v * @as(f64, @bitCast(atnat2.cij[@intCast(i)][6]))))));
    const t1: f64 = @as(f64, @bitCast(atnat2.opi)) - @as(f64, @bitCast(atnat2.cij[@intCast(i)][1]));
    const z: f64 = t1 + zz;
    // Max ULP is 0.502.
    return float.copysign(z, y);
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

    if (((hx - 0x3fff000000000000) | lx) == 0) return float.atan(y); // x = 1.0

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
        z = float.atan(@abs(y / x));
    }

    switch (m) {
        0 => return z, // atan(+,+)
        1 => return -z, // atan(-,+)

        2 => return pi - (z - pi_lo), // atan(+,-)
        else => return (z - pi_lo) - pi, // case 3, atan(-,-)
    }
}
