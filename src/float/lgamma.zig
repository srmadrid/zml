const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const mul_split = @import("mul_split.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const erf_data = @import("erf_data.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub fn lgamma(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return lgamma(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            var local_signgam: i32 = undefined;
            const y: @TypeOf(x) = switch (@TypeOf(x)) {
                f16 => cast(f16, lgamma_r32(cast(f32, x, .{}), &local_signgam), .{}),
                f32 => lgamma_r32(x, &local_signgam),
                f64 => lgamma_r64(x, &local_signgam),
                f80 => cast(f80, lgamma_r128(cast(f128, x, .{}), &local_signgam), .{}),
                f128 => lgamma_r128(x, &local_signgam),
                else => unreachable,
            };
            return y;
        },
        else => unreachable,
    }
}

fn as_r7(x: f64, c: []const f64) f64 {
    return (((x - c[0]) * (x - c[1])) * ((x - c[2]) * (x - c[3]))) * (((x - c[4]) * (x - c[5])) * ((x - c[6])));
}

fn as_r8(x: f64, c: []const f64) f64 {
    return (((x - c[0]) * (x - c[1])) * ((x - c[2]) * (x - c[3]))) * (((x - c[4]) * (x - c[5])) * ((x - c[6]) * (x - c[7])));
}

fn as_sinpi(x: f64) f64 {
    const c: [8]f64 = .{
        0x1p+2,                -0x1.de9e64df22ea4p+1,  0x1.472be122401f8p+0,
        -0x1.d4fcd82df91bp-3,  0x1.9f05c97e0aab2p-6,   -0x1.f3091c427b611p-10,
        0x1.b22c9bfdca547p-14, -0x1.15484325ef569p-18,
    };
    const xx: f64 = x - 0.5;
    const x2: f64 = xx * xx;
    const x4: f64 = x2 * x2;
    const x8: f64 = x4 * x4;
    return (0.25 - x2) * ((c[0] + x2 * c[1]) + x4 * (c[2] + x2 * c[3]) + x8 * ((c[4] + x2 * c[5]) + x4 * (c[6] + x2 * c[7])));
}

fn as_ln(x: f64) f64 {
    var t: u64 = @bitCast(x);
    const e: i32 = @bitCast(@as(u32, @truncate((t >> 52) -% 0x3ff)));
    const c: [8]f64 = .{
        0x1.fffffffffff24p-1,  -0x1.ffffffffd1d67p-2, 0x1.55555537802dep-2,
        -0x1.ffffeca81b866p-3, 0x1.999611761d772p-3,  -0x1.54f3e581b61bfp-3,
        0x1.1e642b4cb5143p-3,  -0x1.9115a5af1e1edp-4,
    };
    const il: [16]f64 = .{
        0x1.59caeec280116p-57, 0x1.f0a30c01162aap-5, 0x1.e27076e2af2ebp-4,
        0x1.5ff3070a793d6p-3,  0x1.c8ff7c79a9a2p-3,  0x1.1675cababa60fp-2,
        0x1.4618bc21c5ec2p-2,  0x1.739d7f6bbd007p-2, 0x1.9f323ecbf984dp-2,
        0x1.c8ff7c79a9a21p-2,  0x1.f128f5faf06ecp-2, 0x1.0be72e4252a83p-1,
        0x1.1e85f5e7040d1p-1,  0x1.307d7334f10bep-1, 0x1.41d8fe84672afp-1,
        0x1.52a2d265bc5abp-1,
    };
    const ix: [16]f64 = .{
        0x1p+0,               0x1.e1e1e1e1e1e1ep-1, 0x1.c71c71c71c71cp-1,
        0x1.af286bca1af28p-1, 0x1.999999999999ap-1, 0x1.8618618618618p-1,
        0x1.745d1745d1746p-1, 0x1.642c8590b2164p-1, 0x1.5555555555555p-1,
        0x1.47ae147ae147bp-1, 0x1.3b13b13b13b14p-1, 0x1.2f684bda12f68p-1,
        0x1.2492492492492p-1, 0x1.1a7b9611a7b96p-1, 0x1.1111111111111p-1,
        0x1.0842108421084p-1,
    };
    const i: i32 = @bitCast(@as(u32, @truncate((t >> 48) & 0xf)));
    t = (t & (~@as(u64, 0) >> 12)) | (0x3ff << 52);
    const z: f64 = ix[@intCast(i)] * @as(f64, @bitCast(t)) - 1;
    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    return cast(f64, e, .{}) * 0x1.62e42fefa39efp-1 + il[@intCast(i)] + z * ((c[0] + z * c[1]) + z2 * (c[2] + z * c[3]) + z4 * ((c[4] + z * c[5]) + z2 * (c[6] + z * c[7])));
}

fn lgamma_r32(x: f32, signgamp: *i32) f32 {
    const tb: [35]struct { x: f32, f: f32, df: f32 } = .{
        .{ .x = -0x1.efc2a2p+14, .f = -0x1.222dbcp+18, .df = -0x1p-7 },
        .{ .x = -0x1.627346p+7, .f = -0x1.73235ep+9, .df = -0x1p-16 },
        .{ .x = -0x1.08b14p+4, .f = -0x1.f0cbe6p+4, .df = -0x1p-21 },
        .{ .x = -0x1.69d628p+3, .f = -0x1.0eac2ap+4, .df = -0x1p-21 },
        .{ .x = -0x1.904902p+2, .f = -0x1.65532cp+2, .df = 0x1p-23 },
        .{ .x = -0x1.9272d2p+1, .f = -0x1.170b98p-8, .df = 0x1p-33 },
        .{ .x = -0x1.625edap+1, .f = 0x1.6a6c4ap-5, .df = -0x1p-30 },
        .{ .x = -0x1.5fc2aep+1, .f = 0x1.c0a484p-11, .df = -0x1p-36 },
        .{ .x = -0x1.5fb43ep+1, .f = 0x1.5b697p-17, .df = 0x1p-42 },
        .{ .x = -0x1.5fa20cp+1, .f = -0x1.132f7ap-10, .df = 0x1p-35 },
        .{ .x = -0x1.580c1ep+1, .f = -0x1.5787c6p-4, .df = 0x1p-29 },
        .{ .x = -0x1.3a7fcap+1, .f = -0x1.e4cf24p-24, .df = -0x1p-49 },
        .{ .x = -0x1.c2f04p-30, .f = 0x1.43a6f6p+4, .df = 0x1p-21 },
        .{ .x = -0x1.ade594p-30, .f = 0x1.446ab2p+4, .df = -0x1p-21 },
        .{ .x = -0x1.437e74p-40, .f = 0x1.b7dec2p+4, .df = -0x1p-21 },
        .{ .x = -0x1.d85bfep-43, .f = 0x1.d31592p+4, .df = -0x1p-21 },
        .{ .x = -0x1.f51c8ep-49, .f = 0x1.0a572ap+5, .df = -0x1p-20 },
        .{ .x = -0x1.108a5ap-66, .f = 0x1.6d7b18p+5, .df = -0x1p-20 },
        .{ .x = -0x1.ecf3fep-73, .f = 0x1.8f8e5ap+5, .df = -0x1p-20 },
        .{ .x = -0x1.25cb66p-123, .f = 0x1.547a44p+6, .df = -0x1p-19 },
        .{ .x = 0x1.ecf3fep-73, .f = 0x1.8f8e5ap+5, .df = -0x1p-20 },
        .{ .x = 0x1.108a5ap-66, .f = 0x1.6d7b18p+5, .df = -0x1p-20 },
        .{ .x = 0x1.a68bbcp-42, .f = 0x1.c9c6e8p+4, .df = 0x1p-21 },
        .{ .x = 0x1.ddfd06p-12, .f = 0x1.ec5ba8p+2, .df = -0x1p-23 },
        .{ .x = 0x1.f8a754p-9, .f = 0x1.63acc2p+2, .df = 0x1p-23 },
        .{ .x = 0x1.8d16b2p+5, .f = 0x1.1e4b4ep+7, .df = 0x1p-18 },
        .{ .x = 0x1.359e0ep+10, .f = 0x1.d9ad02p+12, .df = -0x1p-13 },
        .{ .x = 0x1.a82a2cp+13, .f = 0x1.c38036p+16, .df = 0x1p-9 },
        .{ .x = 0x1.62c646p+14, .f = 0x1.9075bep+17, .df = -0x1p-8 },
        .{ .x = 0x1.7f298p+31, .f = 0x1.f44946p+35, .df = -0x1p+10 },
        .{ .x = 0x1.a45ea4p+33, .f = 0x1.25dcbcp+38, .df = -0x1p+13 },
        .{ .x = 0x1.f9413ep+76, .f = 0x1.9d5ab4p+82, .df = -0x1p+57 },
        .{ .x = 0x1.dcbbaap+99, .f = 0x1.fc5772p+105, .df = 0x1p+80 },
        .{ .x = 0x1.58ace8p+112, .f = 0x1.9e4f66p+118, .df = -0x1p+93 },
        .{ .x = 0x1.87bdfp+115, .f = 0x1.e465aep+121, .df = 0x1p+96 },
    };

    const fx: f32 = float.floor(x);
    const ax: f32 = float.abs(x);
    var t: u32 = @bitCast(ax);
    if (t >= (0xff << 23)) {
        @branchHint(.unlikely);
        signgamp.* = 1;
        if (t == (0xff << 23))
            return std.math.inf(f32);

        return x + x; // nan
    }

    if (fx == x) {
        @branchHint(.unlikely);
        if (x <= 0) {
            signgamp.* = if (@as(u32, @bitCast(x)) >> 31 != 0) -1 else 1;
            return 1 / @as(f32, 0);
        }
        if (x == 1 or x == 2) {
            signgamp.* = 1;
            return 0;
        }
    }

    // Check the value of fx to avoid a spurious invalid exception.
    // Note that for a binary32 |x| >= 2^23, x is necessarily an integer,
    // and we already dealed with negative integers, thus now:
    // -2^23 < x < +Inf and x is not a negative integer nor 0, 1, 2.
    if (fx >= 0) {
        @branchHint(.likely);
        signgamp.* = 1;
    } else {
        // gamma(x) is negative in (-2n-1,-2n), thus when fx is odd.
        signgamp.* = 1 - (((cast(i32, fx, .{})) & 1) << 1);
    }

    const z: f64 = cast(f64, ax, .{});
    var f: f64 = undefined;
    if (ax < 0x1.52p-1) {
        @branchHint(.unlikely);
        const rn: [8]f64 = .{
            -0x1.505bdf4b65acp+4,  -0x1.51c80eb47e068p+2,
            0x1.0000000007cb8p+0,  -0x1.4ac529250a1fcp+1,
            -0x1.a8c99dbe1621ap+0, -0x1.4abdcc74115eap+0,
            -0x1.1b87fe5a5b923p+0, -0x1.05b8a4d47ff64p+0,
        };
        const c0: f64 = 0x1.0fc0fad268c4dp+2;
        const rd: [8]f64 = .{
            -0x1.4db2cfe9a5265p+5, -0x1.062e99d1c4f27p+3,
            -0x1.c81bc2ecf25f6p+1, -0x1.108e55c10091bp+1,
            -0x1.7dd25af0b83d4p+0, -0x1.36bf1880125fcp+0,
            -0x1.1379fc8023d9cp+0, -0x1.03712e41525d2p+0,
        };
        const s: f64 = cast(f64, x, .{});
        f = (c0 * s) * as_r8(s, &rn) / as_r8(s, &rd) - as_ln(z);
    } else {
        if (ax > 0x1.afc1ap+1) {
            if (x > 0x1.895f1cp+121) {
                @branchHint(.unlikely);
                return 0x1p127 * 0x1p127;
            }

            // |x|>=2**23, must be -integer
            if (x < 0 and ax > 0x1p+23) {
                @branchHint(.unlikely);
                return ax / @as(f32, 0);
            }

            const lz: f64 = as_ln(z);
            f = (z - 0.5) * (lz - 1) + 0x1.acfe390c97d69p-2;
            if (ax < 0x1.0p+20) {
                const iz: f64 = 1 / z;
                const iz2: f64 = iz * iz;
                if (ax > 1198) {
                    f += iz * 1 / 12;
                } else if (ax > 0x1.279a7p+6) {
                    const c: [2]f64 = .{ 0x1.555555547fbadp-4, -0x1.6c0fd270c465p-9 };
                    f += iz * (c[0] + iz2 * c[1]);
                } else if (ax > 0x1.555556p+3) {
                    const c: [4]f64 = .{
                        0x1.555555554de0bp-4,  -0x1.6c16bdc45944fp-9,
                        0x1.a0077f300ecb3p-11, -0x1.2e9cfff3b29c2p-11,
                    };
                    const iz4: f64 = iz2 * iz2;
                    f += iz * ((c[0] + iz2 * c[1]) + iz4 * (c[2] + iz2 * c[3]));
                } else {
                    const c: [8]f64 = .{
                        0x1.5555555551286p-4,  -0x1.6c16c0e7c4cf4p-9,
                        0x1.a0193267fe6f2p-11, -0x1.37e87ec19cb45p-11,
                        0x1.b40011dfff081p-11, -0x1.c16c8946b19b6p-10,
                        0x1.e9f47ace150d8p-9,  -0x1.4f5843a71a338p-8,
                    };
                    const iz4: f64 = iz2 * iz2;
                    const iz8: f64 = iz4 * iz4;
                    const p: f64 = ((c[0] + iz2 * c[1]) + iz4 * (c[2] + iz2 * c[3])) + iz8 * ((c[4] + iz2 * c[5]) + iz4 * (c[6] + iz2 * c[7]));
                    f += iz * p;
                }
            }
            if (x < 0) {
                f = 0x1.250d048e7a1bdp+0 - f - lz;
                const lp: f64 = as_ln(as_sinpi(x - fx));
                f -= lp;
            }
        } else {
            const rn: [7]f64 = .{
                -0x1.667923ff14df7p+5, -0x1.2d35f25ad8f64p+3,
                -0x1.b8c9eab9d5bd3p+1, -0x1.7a4a97f494127p+0,
                -0x1.3a6c8295b4445p-1, -0x1.da44e8b810024p-3,
                -0x1.9061e81c77e4ap-5,
            };

            if (x < 0) {
                const ni: i32 = cast(i32, float.floor(-2 * x), .{});
                if ((ni & 1) == 0 and cast(f32, ni, .{}) == -2 * x)
                    return 1 / @as(f32, 0);
            }

            const c0: f64 = 0x1.3cc0e6a0106b3p+2;
            const rd: [8]f64 = .{
                -0x1.491a899e84c52p+6, -0x1.d202961b9e098p+3,
                -0x1.4ced68c631ed6p+2, -0x1.2589eedf40738p+1,
                -0x1.1302e3337271p+0,  -0x1.c36b802f26dffp-2,
                -0x1.3ded448acc39dp-3, -0x1.bffc491078eafp-6,
            };
            f = (z - 1) * (z - 2) * c0 * as_r7(z, &rn) / as_r8(z, &rd);

            if (x < 0) {
                if (t < 0x40301b93 and t > 0x402f95c2) {
                    @branchHint(.unlikely);
                    const h: f64 = (cast(f64, x, .{}) + 0x1.5fb410a1bd901p+1) - 0x1.a19a96d2e6f85p-54;
                    const h2: f64 = h * h;
                    const h4: f64 = h2 * h2;
                    const c: [8]f64 = .{
                        -0x1.ea12da904b18cp+0,  0x1.3267f3c265a54p+3,
                        -0x1.4185ac30cadb3p+4,  0x1.f504accc3f2e4p+5,
                        -0x1.8588444c679b4p+7,  0x1.43740491dc22p+9,
                        -0x1.12400ea23f9e6p+11, 0x1.dac829f365795p+12,
                    };
                    f = h * ((c[0] + h * c[1]) + h2 * (c[2] + h * c[3]) + h4 * ((c[4] + h * c[5]) + h2 * (c[6] + h * c[7])));
                } else if (t > 0x401ceccb and t < 0x401d95ca) {
                    @branchHint(.unlikely);
                    const h: f64 = (cast(f64, x, .{}) + 0x1.3a7fc9600f86cp+1) + 0x1.55f64f98af8dp-55;
                    const h2: f64 = h * h;
                    const h4: f64 = h2 * h2;
                    const c: [7]f64 = .{
                        0x1.83fe966af535fp+0, 0x1.36eebb002f61ap+2,
                        0x1.694a60589a0b3p+0, 0x1.1718d7aedb0b5p+3,
                        0x1.733a045eca0d3p+2, 0x1.8d4297421205bp+4,
                        0x1.7feea5fb29965p+4,
                    };
                    f = h * ((c[0] + h * c[1]) + h2 * (c[2] + h * c[3]) + h4 * ((c[4] + h * c[5]) + h2 * (c[6])));
                } else if (t > 0x40492009 and t < 0x404940ef) {
                    @branchHint(.unlikely);
                    const h: f64 = (cast(f64, x, .{}) + 0x1.9260dbc9e59afp+1) + 0x1.f717cd335a7b3p-53;
                    const h2: f64 = h * h;
                    const h4: f64 = h2 * h2;
                    const c: [7]f64 = .{
                        0x1.f20a65f2fac55p+2,  0x1.9d4d297715105p+4,
                        0x1.c1137124d5b21p+6,  0x1.267203d24de38p+9,
                        0x1.99a63399a0b44p+11, 0x1.2941214faaf0cp+14,
                        0x1.bb912c0c9cdd1p+16,
                    };
                    f = h * ((c[0] + h * c[1]) + h2 * (c[2] + h * c[3]) + h4 * ((c[4] + h * c[5]) + h2 * (c[6])));
                } else {
                    f = 0x1.250d048e7a1bdp+0 - f;
                    const lp: f64 = as_ln(as_sinpi(cast(f64, x - fx, .{})) * z);
                    f -= lp;
                }
            }
        }
    }

    const tl: u64 = (@as(u64, @bitCast(f)) + 5) & 0xfffffff;
    const r: f32 = cast(f32, f, .{});
    if (tl <= 31) {
        @branchHint(.unlikely);
        t = @bitCast(x);
        var i: u32 = 0;
        while (i < 32) {
            if (t == @as(u32, @bitCast(tb[i].x)))
                return tb[i].f + tb[i].df;

            i += 1;
        }
    }
    return r;
}

// Compute the product of 1 + (T / (X + X_EPS)), 1 + (T / (X + X_EPS +
// 1)), ..., 1 + (T / (X + X_EPS + N - 1)), minus 1.  X is such that
// all the values X + 1, ..., X + N - 1 are exactly representable, and
// X_EPS / X is small enough that factors quadratic in it can be
// neglected.
fn lgamma_product64(t: f64, x: f64, x_eps: f64, n: i32) f64 {
    var ret: f64 = 0;
    var ret_eps: f64 = 0;
    var i: i32 = 0;
    while (i < n) {
        const xi: f64 = x + cast(f64, i, .{});
        const quot: f64 = t / xi;
        var mhi: f64 = undefined;
        var mlo: f64 = undefined;
        mul_split.mul_split64(&mhi, &mlo, quot, xi);
        const quot_lo: f64 = (t - mhi - mlo) / xi - t * x_eps / (xi * xi);
        // We want (1 + RET + RET_EPS) * (1 + QUOT + QUOT_LO) - 1.
        var rhi: f64 = undefined;
        var rlo: f64 = undefined;
        mul_split.mul_split64(&rhi, &rlo, ret, quot);
        const rpq: f64 = ret + quot;
        const rpq_eps: f64 = (ret - rpq) + quot;
        const nret: f64 = rpq + rhi;
        const nret_eps: f64 = (rpq - nret) + rhi;
        ret_eps += (rpq_eps + nret_eps + rlo + ret_eps * quot + quot_lo + quot_lo * (ret + ret_eps));
        ret = nret;

        i += 1;
    }

    return ret + ret_eps;
}

// Compute sin (pi * X) for -0.25 <= X <= 0.5.
fn lg_sinpi64(x: f64) f64 {
    if (x <= 0.25) {
        return float.sin(std.math.pi * x);
    } else {
        return float.cos(std.math.pi * (0.5 - x));
    }
}

// Compute cos (pi * X) for -0.25 <= X <= 0.5.
fn lg_cospi64(x: f64) f64 {
    if (x <= 0.25) {
        return float.cos(std.math.pi * x);
    } else {
        return float.sin(std.math.pi * (0.5 - x));
    }
}

// Compute cot (pi * X) for -0.25 <= X <= 0.5.
fn lg_cotpi64(x: f64) f64 {
    return lg_cospi64(x) / lg_sinpi64(x);
}

// Compute lgamma of a negative argument -28 < X < -2, setting
// *SIGNGAMP accordingly.
fn lgamma_neg64(x: f64, signgamp: *i32) f64 {
    const lgamma_zeros: [52][2]f64 =
        .{
            .{ -0x2.74ff92c01f0d8p+0, -0x2.abec9f315f1ap-56 },
            .{ -0x2.bf6821437b202p+0, 0x6.866a5b4b9be14p-56 },
            .{ -0x3.24c1b793cb35ep+0, -0xf.b8be699ad3d98p-56 },
            .{ -0x3.f48e2a8f85fcap+0, -0x1.70d4561291237p-56 },
            .{ -0x4.0a139e1665604p+0, 0xf.3c60f4f21e7fp-56 },
            .{ -0x4.fdd5de9bbabf4p+0, 0xa.ef2f55bf89678p-56 },
            .{ -0x5.021a95fc2db64p+0, -0x3.2a4c56e595394p-56 },
            .{ -0x5.ffa4bd647d034p+0, -0x1.7dd4ed62cbd32p-52 },
            .{ -0x6.005ac9625f234p+0, 0x4.9f83d2692e9c8p-56 },
            .{ -0x6.fff2fddae1bcp+0, 0xc.29d949a3dc03p-60 },
            .{ -0x7.000cff7b7f87cp+0, 0x1.20bb7d2324678p-52 },
            .{ -0x7.fffe5fe05673cp+0, -0x3.ca9e82b522b0cp-56 },
            .{ -0x8.0001a01459fc8p+0, -0x1.f60cb3cec1cedp-52 },
            .{ -0x8.ffffd1c425e8p+0, -0xf.fc864e9574928p-56 },
            .{ -0x9.00002e3bb47d8p+0, -0x6.d6d843fedc35p-56 },
            .{ -0x9.fffffb606bep+0, 0x2.32f9d51885afap-52 },
            .{ -0xa.0000049f93bb8p+0, -0x1.927b45d95e154p-52 },
            .{ -0xa.ffffff9466eap+0, 0xe.4c92532d5243p-56 },
            .{ -0xb.0000006b9915p+0, -0x3.15d965a6ffea4p-52 },
            .{ -0xb.fffffff708938p+0, -0x7.387de41acc3d4p-56 },
            .{ -0xc.00000008f76c8p+0, 0x8.cea983f0fdafp-56 },
            .{ -0xc.ffffffff4f6ep+0, 0x3.09e80685a0038p-52 },
            .{ -0xd.00000000b092p+0, -0x3.09c06683dd1bap-52 },
            .{ -0xd.fffffffff3638p+0, 0x3.a5461e7b5c1f6p-52 },
            .{ -0xe.000000000c9c8p+0, -0x3.a545e94e75ec6p-52 },
            .{ -0xe.ffffffffff29p+0, 0x3.f9f399fb10cfcp-52 },
            .{ -0xf.0000000000d7p+0, -0x3.f9f399bd0e42p-52 },
            .{ -0xf.fffffffffff28p+0, -0xc.060c6621f513p-56 },
            .{ -0x1.000000000000dp+4, -0x7.3f9f399da1424p-52 },
            .{ -0x1.0ffffffffffffp+4, -0x3.569c47e7a93e2p-52 },
            .{ -0x1.1000000000001p+4, 0x3.569c47e7a9778p-52 },
            .{ -0x1.2p+4, 0xb.413c31dcbecdp-56 },
            .{ -0x1.2p+4, -0xb.413c31dcbeca8p-56 },
            .{ -0x1.3p+4, 0x9.7a4da340a0ab8p-60 },
            .{ -0x1.3p+4, -0x9.7a4da340a0ab8p-60 },
            .{ -0x1.4p+4, 0x7.950ae90080894p-64 },
            .{ -0x1.4p+4, -0x7.950ae90080894p-64 },
            .{ -0x1.5p+4, 0x5.c6e3bdb73d5c8p-68 },
            .{ -0x1.5p+4, -0x5.c6e3bdb73d5c8p-68 },
            .{ -0x1.6p+4, 0x4.338e5b6dfe14cp-72 },
            .{ -0x1.6p+4, -0x4.338e5b6dfe14cp-72 },
            .{ -0x1.7p+4, 0x2.ec368262c7034p-76 },
            .{ -0x1.7p+4, -0x2.ec368262c7034p-76 },
            .{ -0x1.8p+4, 0x1.f2cf01972f578p-80 },
            .{ -0x1.8p+4, -0x1.f2cf01972f578p-80 },
            .{ -0x1.9p+4, 0x1.3f3ccdd165fa9p-84 },
            .{ -0x1.9p+4, -0x1.3f3ccdd165fa9p-84 },
            .{ -0x1.ap+4, 0xc.4742fe35272dp-92 },
            .{ -0x1.ap+4, -0xc.4742fe35272dp-92 },
            .{ -0x1.bp+4, 0x7.46ac70b733a8cp-96 },
            .{ -0x1.bp+4, -0x7.46ac70b733a8cp-96 },
            .{ -0x1.cp+4, 0x4.2862898d42174p-100 },
        };

    const e_hi: f64 = 0x2.b7e151628aed2p+0;
    const e_lo: f64 = 0xa.6abf7158809dp-56;

    // Coefficients B_2k / 2k(2k-1) of x^-(2k-1) in Stirling's
    // approximation to lgamma function.
    const lgamma_coeff: [12]f64 = .{
        0x1.5555555555555p-4,
        -0xb.60b60b60b60b8p-12,
        0x3.4034034034034p-12,
        -0x2.7027027027028p-12,
        0x3.72a3c5631fe46p-12,
        -0x7.daac36664f1f4p-12,
        0x1.a41a41a41a41ap-8,
        -0x7.90a1b2c3d4e6p-8,
        0x2.dfd2c703c0dp-4,
        -0x1.6476701181f3ap+0,
        0xd.672219167003p+0,
        -0x9.cd9292e6660d8p+4,
    };

    // Polynomial approximations to (|gamma(x)|-1)(x-n)/(x-x0), where n is
    // the integer end-point of the half-integer interval containing x and
    // x0 is the zero of lgamma in that half-integer interval.  Each
    // polynomial is expressed in terms of x-xm, where xm is the midpoint
    // of the interval for which the polynomial applies.
    const poly_coeff: [101]f64 = .{
        // Interval [-2.125, -2] (polynomial degree 10).
        -0x1.0b71c5c54d42fp+0,
        -0xc.73a1dc05f3758p-4,
        -0x1.ec84140851911p-4,
        -0xe.37c9da23847e8p-4,
        -0x1.03cd87cdc0ac6p-4,
        -0xe.ae9aedce12eep-4,
        0x9.b11a1780cfd48p-8,
        -0xe.f25fc460bdebp-4,
        0x2.6e984c61ca912p-4,
        -0xf.83fea1c6d35p-4,
        0x4.760c8c8909758p-4,
        // Interval [-2.25, -2.125] (polynomial degree 11).
        -0xf.2930890d7d678p-4,
        -0xc.a5cfde054eaa8p-4,
        0x3.9c9e0fdebd99cp-4,
        -0x1.02a5ad35619d9p+0,
        0x9.6e9b1167c164p-4,
        -0x1.4d8332eba090ap+0,
        0x1.1c0c94b1b2b6p+0,
        -0x1.c9a70d138c74ep+0,
        0x1.d7d9cf1d4c196p+0,
        -0x2.91fbf4cd6abacp+0,
        0x2.f6751f74b8ff8p+0,
        -0x3.e1bb7b09e3e76p+0,
        // Interval [-2.375, -2.25] (polynomial degree 12).
        -0xd.7d28d505d618p-4,
        -0xe.69649a3040958p-4,
        0xb.0d74a2827cd6p-4,
        -0x1.924b09228a86ep+0,
        0x1.d49b12bcf6175p+0,
        -0x3.0898bb530d314p+0,
        0x4.207a6be8fda4cp+0,
        -0x6.39eef56d4e9p+0,
        0x8.e2e42acbccec8p+0,
        -0xd.0d91c1e596a68p+0,
        0x1.2e20d7099c585p+4,
        -0x1.c4eb6691b4ca9p+4,
        0x2.96a1a11fd85fep+4,
        // Interval [-2.5, -2.375] (polynomial degree 13).
        -0xb.74ea1bcfff948p-4,
        -0x1.2a82bd590c376p+0,
        0x1.88020f828b81p+0,
        -0x3.32279f040d7aep+0,
        0x5.57ac8252ce868p+0,
        -0x9.c2aedd093125p+0,
        0x1.12c132716e94cp+4,
        -0x1.ea94dfa5c0a6dp+4,
        0x3.66b61abfe858cp+4,
        -0x6.0cfceb62a26e4p+4,
        0xa.beeba09403bd8p+4,
        -0x1.3188d9b1b288cp+8,
        0x2.37f774dd14c44p+8,
        -0x3.fdf0a64cd7136p+8,
        // Interval [-2.625, -2.5] (polynomial degree 13).
        -0x3.d10108c27ebbp-4,
        0x1.cd557caff7d2fp+0,
        0x3.819b4856d36cep+0,
        0x6.8505cbacfc42p+0,
        0xb.c1b2e6567a4dp+0,
        0x1.50a53a3ce6c73p+4,
        0x2.57adffbb1ec0cp+4,
        0x4.2b15549cf400cp+4,
        0x7.698cfd82b3e18p+4,
        0xd.2decde217755p+4,
        0x1.7699a624d07b9p+8,
        0x2.98ecf617abbfcp+8,
        0x4.d5244d44d60b4p+8,
        0x8.e962bf7395988p+8,
        // Interval [-2.75, -2.625] (polynomial degree 12).
        -0x6.b5d252a56e8a8p-4,
        0x1.28d60383da3a6p+0,
        0x1.db6513ada89bep+0,
        0x2.e217118fa8c02p+0,
        0x4.450112c651348p+0,
        0x6.4af990f589b8cp+0,
        0x9.2db5963d7a238p+0,
        0xd.62c03647da19p+0,
        0x1.379f81f6416afp+4,
        0x1.c5618b4fdb96p+4,
        0x2.9342d0af2ac4ep+4,
        0x3.d9cdf56d2b186p+4,
        0x5.ab9f91d5a27a4p+4,
        // Interval [-2.875, -2.75] (polynomial degree 11).
        -0x8.a41b1e4f36ff8p-4,
        0xc.da87d3b69dbe8p-4,
        0x1.1474ad5c36709p+0,
        0x1.761ecb90c8c5cp+0,
        0x1.d279bff588826p+0,
        0x2.4e5d003fb36a8p+0,
        0x2.d575575566842p+0,
        0x3.85152b0d17756p+0,
        0x4.5213d921ca13p+0,
        0x5.55da7dfcf69c4p+0,
        0x6.acef729b9404p+0,
        0x8.483cc21dd0668p+0,
        // Interval [-3, -2.875] (polynomial degree 11).
        -0xa.046d667e468f8p-4,
        0x9.70b88dcc006cp-4,
        0xa.a8a39421c94dp-4,
        0xd.2f4d1363f98ep-4,
        0xd.ca9aa19975b7p-4,
        0xf.cf09c2f54404p-4,
        0x1.04b1365a9adfcp+0,
        0x1.22b54ef213798p+0,
        0x1.2c52c25206bf5p+0,
        0x1.4aa3d798aace4p+0,
        0x1.5c3f278b504e3p+0,
        0x1.7e08292cc347bp+0,
    };

    const poly_deg: [8]u32 = .{
        10,
        11,
        12,
        13,
        13,
        12,
        11,
        11,
    };

    const poly_end: [8]u32 = .{
        10,
        22,
        35,
        49,
        63,
        76,
        88,
        100,
    };

    // Determine the half-integer region X lies in, handle exact
    // integers and determine the sign of the result.
    var i: i32 = cast(i32, float.floor(-2 * x), .{});
    if ((i & 1) == 0 and cast(f64, i, .{}) == -2 * x)
        return 1 / @as(f64, 0);

    const xn: f64 = cast(f64, if ((i & 1) == 0) @divFloor(-i, 2) else @divFloor(-i - 1, 2), .{});
    i -= 4;
    signgamp.* = if ((i & 2) == 0) -1 else 1;

    // Expand around the zero X0 = X0_HI + X0_LO.
    const x0_hi: f64 = lgamma_zeros[@intCast(i)][0];
    const x0_lo: f64 = lgamma_zeros[@intCast(i)][1];
    const xdiff: f64 = x - x0_hi - x0_lo;

    // For arguments in the range -3 to -2, use polynomial
    // approximations to an adjusted version of the gamma function.
    if (i < 2) {
        var j: i32 = cast(i32, float.floor(-8 * x) - 16, .{});
        const xm: f64 = (-33 - 2 * cast(f64, j, .{})) * 0.0625;
        const x_adj: f64 = x - xm;
        const deg: u32 = poly_deg[@intCast(j)];
        const end: u32 = poly_end[@intCast(j)];
        var g: f64 = poly_coeff[end];
        j = 1;
        while (j <= deg) {
            g = g * x_adj + poly_coeff[end - cast(u32, j, .{})];
            j += 1;
        }

        return float.log1p(g * xdiff / (x - xn));
    }

    // The result we want is log (sinpi (X0) / sinpi (X))
    //  + log (gamma (1 - X0) / gamma (1 - X)).
    const x_idiff: f64 = float.abs(xn - x);
    const x0_idiff: f64 = float.abs(xn - x0_hi - x0_lo);
    var log_sinpi_ratio: f64 = undefined;
    if (x0_idiff < x_idiff * 0.5) {
        // Use log not log1p to avoid inaccuracy from log1p of arguments
        // close to -1.
        log_sinpi_ratio = float.log(lg_sinpi64(x0_idiff) / lg_sinpi64(x_idiff));
    } else {
        // Use log1p not log to avoid inaccuracy from log of arguments
        // close to 1.  X0DIFF2 has positive sign if X0 is further from
        // XN than X is from XN, negative sign otherwise.
        const x0diff2: f64 = (if ((i & 1) == 0) xdiff else -xdiff) * 0.5;
        const sx0d2: f64 = lg_sinpi64(x0diff2);
        const cx0d2: f64 = lg_cospi64(x0diff2);
        log_sinpi_ratio = float.log1p(2 * sx0d2 * (-sx0d2 + cx0d2 * lg_cotpi64(x_idiff)));
    }

    var y0: f64 = 1 - x0_hi;
    var y0_eps: f64 = -x0_hi + (1 - y0) - x0_lo;
    var y: f64 = 1 - x;
    var y_eps: f64 = -x + (1 - y);
    // We now wish to compute LOG_GAMMA_RATIO
    // = log (gamma (Y0 + Y0_EPS) / gamma (Y + Y_EPS)).  XDIFF
    // accurately approximates the difference Y0 + Y0_EPS - Y -
    // Y_EPS.  Use Stirling's approximation.  First, we may need to
    // adjust into the range where Stirling's approximation is
    // sufficiently accurate.
    var log_gamma_adj: f64 = 0;
    if (i < 6) {
        const n_up: i32 = @divFloor(7 - i, 2);
        const ny0: f64 = y0 + cast(f64, n_up, .{});
        const ny0_eps: f64 = y0 - (ny0 - cast(f64, n_up, .{})) + y0_eps;
        y0 = ny0;
        y0_eps = ny0_eps;
        const ny: f64 = y + cast(f64, n_up, .{});
        const ny_eps: f64 = y - (ny - cast(f64, n_up, .{})) + y_eps;
        y = ny;
        y_eps = ny_eps;
        const prodm1: f64 = lgamma_product64(xdiff, y - cast(f64, n_up, .{}), y_eps, n_up);
        log_gamma_adj = -float.log1p(prodm1);
    }

    const log_gamma_high: f64 = (xdiff * float.log1p((y0 - e_hi - e_lo + y0_eps) / e_hi) + (y - 0.5 + y_eps) * float.log1p(xdiff / y) + log_gamma_adj);
    // Compute the sum of (B_2k / 2k(2k-1))(Y0^-(2k-1) - Y^-(2k-1)).
    const y0r: f64 = 1 / y0;
    const yr: f64 = 1 / y;
    const y0r2: f64 = y0r * y0r;
    const yr2: f64 = yr * yr;
    const rdiff: f64 = -xdiff / (y * y0);
    var bterm: [12]f64 = undefined;
    var dlast: f64 = rdiff;
    var elast: f64 = rdiff * yr * (yr + y0r);
    bterm[0] = dlast * lgamma_coeff[0];
    var j: u32 = 1;
    while (j < 12) {
        const dnext: f64 = dlast * y0r2 + elast;
        const enext: f64 = elast * yr2;
        bterm[j] = dnext * lgamma_coeff[j];
        dlast = dnext;
        elast = enext;

        j += 1;
    }
    var log_gamma_low: f64 = 0;
    j = 0;
    while (j < 12) {
        log_gamma_low += bterm[11 - j];

        j += 1;
    }

    const log_gamma_ratio: f64 = log_gamma_high + log_gamma_low;

    return log_sinpi_ratio + log_gamma_ratio;
}

const two52: f64 = 4.50359962737049600000e+15; // 0x43300000, 0x00000000
const pi: f64 = 3.14159265358979311600e+00; // 0x400921fb, 0x54442d18
const a0: f64 = 7.72156649015328655494e-02; // 0x3fb3c467, 0xe37db0c8
const a1: f64 = 3.22467033424113591611e-01; // 0x3fd4a34c, 0xc4a60fad
const a2: f64 = 6.73523010531292681824e-02; // 0x3fb13e00, 0x1a5562a7
const a3: f64 = 2.05808084325167332806e-02; // 0x3f951322, 0xac92547b
const a4: f64 = 7.38555086081402883957e-03; // 0x3f7e404f, 0xb68fefe8
const a5: f64 = 2.89051383673415629091e-03; // 0x3f67add8, 0xccb7926b
const a6: f64 = 1.19270763183362067845e-03; // 0x3f538a94, 0x116f3f5d
const a7: f64 = 5.10069792153511336608e-04; // 0x3f40b6c6, 0x89b99c00
const a8: f64 = 2.20862790713908385557e-04; // 0x3f2cf2ec, 0xed10e54d
const a9: f64 = 1.08011567247583939954e-04; // 0x3f1c5088, 0x987dfb07
const a10: f64 = 2.52144565451257326939e-05; // 0x3efa7074, 0x428cfa52
const a11: f64 = 4.48640949618915160150e-05; // 0x3f07858e, 0x90a45837
const tc: f64 = 1.46163214496836224576e+00; // 0x3ff762d8, 0x6356be3f
const tf: f64 = -1.21486290535849611461e-01; // 0xbfbf19b9, 0xbcc38a42
// tt = -(tail of tf)
const tt: f64 = -3.63867699703950536541e-18; // 0xbc50c7ca, 0xa48a971f
const t0: f64 = 4.83836122723810047042e-01; // 0x3fdef72b, 0xc8ee38a2
const t1: f64 = -1.47587722994593911752e-01; // 0xbfc2e427, 0x8dc6c509
const t2: f64 = 6.46249402391333854778e-02; // 0x3fb08b42, 0x94d5419b
const t3: f64 = -3.27885410759859649565e-02; // 0xbfa0c9a8, 0xdf35b713
const t4: f64 = 1.79706750811820387126e-02; // 0x3f9266e7, 0x970af9ec
const t5: f64 = -1.03142241298341437450e-02; // 0xbf851f9f, 0xba91ec6a
const t6: f64 = 6.10053870246291332635e-03; // 0x3f78fce0, 0xe370e344
const t7: f64 = -3.68452016781138256760e-03; // 0xbf6e2eff, 0xb3e914d7
const t8: f64 = 2.25964780900612472250e-03; // 0x3f6282d3, 0x2e15c915
const t9: f64 = -1.40346469989232843813e-03; // 0xbf56fe8e, 0xbf2d1af1
const t10: f64 = 8.81081882437654011382e-04; // 0x3f4cdf0c, 0xef61a8e9
const t11: f64 = -5.38595305356740546715e-04; // 0xbf41a610, 0x9c73e0ec
const t12: f64 = 3.15632070903625950361e-04; // 0x3f34af6d, 0x6c0ebbf7
const t13: f64 = -3.12754168375120860518e-04; // 0xbf347f24, 0xecc38c38
const t14: f64 = 3.35529192635519073543e-04; // 0x3f35fd3e, 0xe8c2d3f4
const U0: f64 = -7.72156649015328655494e-02; // 0xbfb3c467, 0xe37db0c8
const U1: f64 = 6.32827064025093366517e-01; // 0x3fe4401e, 0x8b005dff
const U2: f64 = 1.45492250137234768737e+00; // 0x3ff7475c, 0xd119bd6f
const U3: f64 = 9.77717527963372745603e-01; // 0x3fef4976, 0x44ea8450
const U4: f64 = 2.28963728064692451092e-01; // 0x3fcd4eae, 0xf6010924
const U5: f64 = 1.33810918536787660377e-02; // 0x3f8b678b, 0xbf2bab09
const v1: f64 = 2.45597793713041134822e+00; // 0x4003a5d7, 0xc2bd619c
const v2: f64 = 2.12848976379893395361e+00; // 0x40010725, 0xa42b18f5
const v3: f64 = 7.69285150456672783825e-01; // 0x3fe89dfb, 0xe45050af
const v4: f64 = 1.04222645593369134254e-01; // 0x3fbaae55, 0xd6537c88
const v5: f64 = 3.21709242282423911810e-03; // 0x3f6a5abb, 0x57d0cf61
const s0: f64 = -7.72156649015328655494e-02; // 0xbfb3c467, 0xe37db0c8
const s1: f64 = 2.14982415960608852501e-01; // 0x3fcb848b, 0x36e20878
const s2: f64 = 3.25778796408930981787e-01; // 0x3fd4d98f, 0x4f139f59
const s3: f64 = 1.46350472652464452805e-01; // 0x3fc2bb9c, 0xbee5f2f7
const s4: f64 = 2.66422703033638609560e-02; // 0x3f9b481c, 0x7e939961
const s5: f64 = 1.84028451407337715652e-03; // 0x3f5e26b6, 0x7368f239
const s6: f64 = 3.19475326584100867617e-05; // 0x3f00bfec, 0xdd17e945
const r1: f64 = 1.39200533467621045958e+00; // 0x3ff645a7, 0x62c4ab74
const r2: f64 = 7.21935547567138069525e-01; // 0x3fe71a18, 0x93d3dcdc
const r3: f64 = 1.71933865632803078993e-01; // 0x3fc601ed, 0xccfbdf27
const r4: f64 = 1.86459191715652901344e-02; // 0x3f9317ea, 0x742ed475
const r5: f64 = 7.77942496381893596434e-04; // 0x3f497dda, 0xca41a95b
const r6: f64 = 7.32668430744625636189e-06; // 0x3edebaf7, 0xa5b38140
const w0: f64 = 4.18938533204672725052e-01; // 0x3fdacfe3, 0x90c97d69
const w1: f64 = 8.33333333333329678849e-02; // 0x3fb55555, 0x5555553b
const w2: f64 = -2.77777777728775536470e-03; // 0xbf66c16c, 0x16b02e5c
const w3: f64 = 7.93650558643019558500e-04; // 0x3f4a019f, 0x98cf38b6
const w4: f64 = -5.95187557450339963135e-04; // 0xbf4380cb, 0x8c0fe741
const w5: f64 = 8.36339918996282139126e-04; // 0x3f4b67ba, 0x4cdad5d1
const w6: f64 = -1.63092934096575273989e-03; // 0xbf5ab89d, 0x0b9e43e4

fn sin_pi64(x: f64) f64 {
    var ix: i32 = undefined;
    dbl64.getHighWord(&ix, x);
    ix &= 0x7fffffff;
    if (ix < 0x3fd00000) return float.sin(pi * x);
    var y: f64 = -x; // x is assume negative

    // argument reduction, make sure inexact flag not raised if input
    // is an integer
    var z: f64 = float.floor(y);
    var n: i32 = undefined;
    if (z != y) { // inexact anyway
        y *= 0.5;
        y = 2.0 * (y - float.floor(y)); // y = |x| mod 2.0
        n = cast(i32, y * 4.0, .{});
    } else {
        if (ix >= 0x43400000) {
            y = 0;
            n = 0; // y must be even
        } else {
            if (ix < 0x43300000) z = y + two52; // exact
            dbl64.getLowWord(&n, z);
            n &= 1;
            y = cast(f64, n, .{});
            n <<= 2;
        }
    }

    switch (n) {
        0 => y = float.sin(pi * y),
        1, 2 => y = float.cos(pi * (0.5 - y)),
        3, 4 => y = float.sin(pi * (1 - y)),
        5, 6 => y = -float.cos(pi * (y - 1.5)),
        else => y = float.sin(pi * (y - 2.0)),
    }

    return -y;
}

pub fn lgamma_r64(x: f64, signgamp: *i32) f64 {
    var hx: i32 = undefined;
    var lx: i32 = undefined;
    dbl64.extractWords(&hx, &lx, x);

    // purge off +-inf, NaN, +-0, and negative arguments
    signgamp.* = 1;
    const ix: i32 = hx & 0x7fffffff;
    if (ix >= 0x7ff00000) {
        @branchHint(.unlikely);
        return x * x;
    }

    if ((ix | lx) == 0) {
        @branchHint(.unlikely);
        if (hx < 0)
            signgamp.* = -1;

        return 1 / float.abs(x);
    }

    if (ix < 0x3b900000) {
        @branchHint(.unlikely);
        // |x|<2**-70, return -log(|x|)
        if (hx < 0) {
            signgamp.* = -1;

            return -float.log(-x);
        } else {
            return -float.log(x);
        }
    }

    var nadj: f64 = undefined;
    var xx: f64 = x;
    if (hx < 0) {
        if (ix >= 0x43300000) {
            @branchHint(.unlikely);
            // |x|>=2**52, must be -integer
            return float.abs(x) / @as(f64, 0);
        }

        if (x < -2.0 and x > -28.0)
            return lgamma_neg64(x, signgamp);

        const t: f64 = sin_pi64(x);

        if (t == 0) return 1 / float.abs(t); // -integer

        nadj = float.log(pi / float.abs(t * x));
        if (t < 0) signgamp.* = -1;
        xx = -x;
    }

    // purge off 1 and 2
    var r: f64 = undefined;
    if ((((ix - 0x3ff00000) | lx) == 0) or (((ix - 0x40000000) | lx) == 0)) {
        r = 0;
        // for x < 2.0
    } else if (ix < 0x40000000) {
        var y: f64 = undefined;
        var i: i32 = undefined;
        if (ix <= 0x3feccccc) { // lgamma(x) = lgamma(x+1)-log(x)
            r = -float.log(xx);
            if (ix >= 0x3fe76944) {
                y = 1 - xx;
                i = 0;
            } else if (ix >= 0x3fcda661) {
                y = xx - (tc - 1);
                i = 1;
            } else {
                y = xx;
                i = 2;
            }
        } else {
            r = 0;
            if (ix >= 0x3ffbb4c3) { // [1.7316,2]
                y = 2.0 - xx;
                i = 0;
            } else if (ix >= 0x3ff3b4c4) { // [1.23,1.73]
                y = xx - tc;
                i = 1;
            } else {
                y = xx - 1;
                i = 2;
            }
        }

        switch (i) {
            0 => {
                const z: f64 = y * y;
                const p1: f64 = a0 + z * (a2 + z * (a4 + z * (a6 + z * (a8 + z * a10))));
                const p2: f64 = z * (a1 + z * (a3 + z * (a5 + z * (a7 + z * (a9 + z * a11)))));
                const p: f64 = y * p1 + p2;
                r += (p - 0.5 * y);
            },
            1 => {
                const z: f64 = y * y;
                const w: f64 = z * y;
                const p1: f64 = t0 + w * (t3 + w * (t6 + w * (t9 + w * t12))); // parallel comp
                const p2: f64 = t1 + w * (t4 + w * (t7 + w * (t10 + w * t13)));
                const p3: f64 = t2 + w * (t5 + w * (t8 + w * (t11 + w * t14)));
                const p: f64 = z * p1 - (tt - w * (p2 + y * p3));
                r += (tf + p);
            },
            2 => {
                const p1: f64 = y * (U0 + y * (U1 + y * (U2 + y * (U3 + y * (U4 + y * U5)))));
                const p2: f64 = 1 + y * (v1 + y * (v2 + y * (v3 + y * (v4 + y * v5))));
                r += (-0.5 * y + p1 / p2);
            },
            else => unreachable,
        }
    } else if (ix < 0x40200000) { // x < 8.0
        const i: i32 = cast(i32, xx, .{});
        const y: f64 = xx - cast(f64, i, .{});
        const p: f64 = y * (s0 + y * (s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6))))));
        const q: f64 = 1 + y * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6)))));
        r = 0.5 * y + p / q;
        var z: f64 = 1; // lgamma(1+s) = log(s) + lgamma(s)
        sw: switch (i) {
            7 => {
                z *= (y + 6.0);
                continue :sw 6;
            },
            6 => {
                z *= (y + 5.0);
                continue :sw 5;
            },
            5 => {
                z *= (y + 4.0);
                continue :sw 4;
            },
            4 => {
                z *= (y + 3.0);
                continue :sw 3;
            },
            3 => {
                z *= (y + 2.0);
                continue :sw 2;
            },
            else => r += float.log(z),
        }
        // 8.0 <= x < 2**58
    } else if (ix < 0x43900000) {
        const t: f64 = float.log(xx);
        const z: f64 = 1 / xx;
        const y: f64 = z * z;
        const w: f64 = w0 + z * (w1 + y * (w2 + y * (w3 + y * (w4 + y * (w5 + y * w6)))));
        r = (xx - 0.5) * (t - 1) + w;
    } else {
        // 2**58 <= x <= inf
        r = xx * (float.log(xx) - 1);
    }

    if (hx < 0) r = nadj - r;

    return r;
}

// Compute the product of 1 + (T / (X + X_EPS)), 1 + (T / (X + X_EPS +
// 1)), ..., 1 + (T / (X + X_EPS + N - 1)), minus 1.  X is such that
// all the values X + 1, ..., X + N - 1 are exactly representable, and
// X_EPS / X is small enough that factors quadratic in it can be
// neglected.
fn lgamma_product128(t: f128, x: f128, x_eps: f128, n: i32) f128 {
    var ret: f128 = 0;
    var ret_eps: f128 = 0;
    var i: i32 = 0;
    while (i < n) {
        const xi: f128 = x + cast(f128, i, .{});
        const quot: f128 = t / xi;
        var mhi: f128 = undefined;
        var mlo: f128 = undefined;
        mul_split.mul_split128(&mhi, &mlo, quot, xi);
        const quot_lo: f128 = (t - mhi - mlo) / xi - t * x_eps / (xi * xi);
        // We want (1 + RET + RET_EPS) * (1 + QUOT + QUOT_LO) - 1.
        var rhi: f128 = undefined;
        var rlo: f128 = undefined;
        mul_split.mul_split128(&rhi, &rlo, ret, quot);
        const rpq: f128 = ret + quot;
        const rpq_eps: f128 = (ret - rpq) + quot;
        const nret: f128 = rpq + rhi;
        const nret_eps: f128 = (rpq - nret) + rhi;
        ret_eps += (rpq_eps + nret_eps + rlo + ret_eps * quot + quot_lo + quot_lo * (ret + ret_eps));
        ret = nret;

        i += 1;
    }

    return ret + ret_eps;
}

// Compute sin (pi * X) for -0.25 <= X <= 0.5.
fn lg_sinpi128(x: f128) f128 {
    if (x <= 0.25) {
        return float.sin(std.math.pi * x);
    } else {
        return float.cos(std.math.pi * (0.5 - x));
    }
}

// Compute cos (pi * X) for -0.25 <= X <= 0.5.
fn lg_cospi128(x: f128) f128 {
    if (x <= 0.25) {
        return float.cos(std.math.pi * x);
    } else {
        return float.sin(std.math.pi * (0.5 - x));
    }
}

// Compute cot (pi * X) for -0.25 <= X <= 0.5.
fn lg_cotpi128(x: f128) f128 {
    return lg_cospi128(x) / lg_sinpi128(x);
}

// Compute lgamma of a negative argument -50 < X < -2, setting
// *SIGNGAMP accordingly.
fn lgamma_neg128(x: f128, signgamp: *i32) f128 {
    const lgamma_zeros: [96][2]f128 = .{
        .{ -0x2.74ff92c01f0d82abec9f315f1a08p+0, 0xe.d3ccb7fb2658634a2b9f6b2ba81p-116 },
        .{ -0x2.bf6821437b20197995a4b4641eaep+0, -0xb.f4b00b4829f961e428533e6ad048p-116 },
        .{ -0x3.24c1b793cb35efb8be699ad3d9bap+0, -0x6.5454cb7fac60e3f16d9d7840c2ep-116 },
        .{ -0x3.f48e2a8f85fca170d4561291236cp+0, -0xc.320a4887d1cb4c711828a75d5758p-116 },
        .{ -0x4.0a139e16656030c39f0b0de18114p+0, 0x1.53e84029416e1242006b2b3d1cfp-112 },
        .{ -0x4.fdd5de9bbabf3510d0aa40769884p+0, -0x1.01d7d78125286f78d1e501f14966p-112 },
        .{ -0x5.021a95fc2db6432a4c56e595394cp+0, -0x1.ecc6af0430d4fe5746fa7233356fp-112 },
        .{ -0x5.ffa4bd647d0357dd4ed62cbd31ecp+0, -0x1.f8e3f8e5deba2d67dbd70dd96ce1p-112 },
        .{ -0x6.005ac9625f233b607c2d96d16384p+0, -0x1.cb86ac569340cf1e5f24df7aab7bp-112 },
        .{ -0x6.fff2fddae1bbff3d626b65c23fd4p+0, 0x1.e0bfcff5c457ebcf4d3ad9674167p-112 },
        .{ -0x7.000cff7b7f87adf4482dcdb98784p+0, 0x1.54d99e35a74d6407b80292df199fp-112 },
        .{ -0x7.fffe5fe05673c3ca9e82b522b0ccp+0, 0x1.62d177c832e0eb42c2faffd1b145p-112 },
        .{ -0x8.0001a01459fc9f60cb3cec1cec88p+0, 0x2.8998835ac7277f7bcef67c47f188p-112 },
        .{ -0x8.ffffd1c425e80ffc864e95749258p+0, -0x1.e7e20210e7f81cf781b44e9d2b02p-112 },
        .{ -0x9.00002e3bb47d86d6d843fedc352p+0, 0x2.14852f613a16291751d2ab751f7ep-112 },
        .{ -0x9.fffffb606bdfdcd062ae77a50548p+0, 0x3.962d1490cc2e8f031c7007eaa1ap-116 },
        .{ -0xa.0000049f93bb9927b45d95e1544p+0, -0x1.e03086db9146a9287bd4f2172d5ap-112 },
        .{ -0xa.ffffff9466e9f1b36dacd2adbd18p+0, -0xd.05a4e458062f3f95345a4d9c9b6p-116 },
        .{ -0xb.0000006b9915315d965a6ffea41p+0, 0x1.b415c6fff233e7b7fdc3a094246fp-112 },
        .{ -0xb.fffffff7089387387de41acc3d4p+0, 0x3.687427c6373bd74a10306e10a28ep-112 },
        .{ -0xc.00000008f76c7731567c0f0250fp+0, -0x3.87920df5675833859190eb128ef6p-112 },
        .{ -0xc.ffffffff4f6dcf617f97a5ffc758p+0, 0x2.ab72d76f32eaee2d1a42ed515d3ap-116 },
        .{ -0xd.00000000b092309c06683dd1b9p+0, -0x3.e3700857a15c19ac5a611de9688ap-112 },
        .{ -0xd.fffffffff36345ab9e184a3e09dp+0, -0x1.176dc48e47f62d917973dd44e553p-112 },
        .{ -0xe.000000000c9cba545e94e75ec57p+0, -0x1.8f753e2501e757a17cf2ecbeeb89p-112 },
        .{ -0xe.ffffffffff28c060c6604ef3037p+0, -0x1.f89d37357c9e3dc17c6c6e63becap-112 },
        .{ -0xf.0000000000d73f9f399bd0e420f8p+0, -0x5.e9ee31b0b890744fc0e3fbc01048p-116 },
        .{ -0xf.fffffffffff28c060c6621f512e8p+0, 0xd.1b2eec9d960bd9adc5be5f5fa5p-116 },
        .{ -0x1.000000000000d73f9f399da1424cp+4, 0x6.c46e0e88305d2800f0e414c506a8p-116 },
        .{ -0x1.0ffffffffffff3569c47e7a93e1cp+4, -0x4.6a08a2e008a998ebabb8087efa2cp-112 },
        .{ -0x1.1000000000000ca963b818568887p+4, -0x6.ca5a3a64ec15db0a95caf2c9ffb4p-112 },
        .{ -0x1.1fffffffffffff4bec3ce234132dp+4, -0x8.b2b726187c841cb92cd5221e444p-116 },
        .{ -0x1.20000000000000b413c31dcbeca5p+4, 0x3.c4d005344b6cd0e7231120294abcp-112 },
        .{ -0x1.2ffffffffffffff685b25cbf5f54p+4, -0x5.ced932e38485f7dd296b8fa41448p-112 },
        .{ -0x1.30000000000000097a4da340a0acp+4, 0x7.e484e0e0ffe38d406ebebe112f88p-112 },
        .{ -0x1.3fffffffffffffff86af516ff7f7p+4, -0x6.bd67e720d57854502b7db75e1718p-112 },
        .{ -0x1.40000000000000007950ae900809p+4, 0x6.bec33375cac025d9c073168c5d9p-112 },
        .{ -0x1.4ffffffffffffffffa391c4248c3p+4, 0x5.c63022b62b5484ba346524db607p-112 },
        .{ -0x1.500000000000000005c6e3bdb73dp+4, -0x5.c62f55ed5322b2685c5e9a51e6a8p-112 },
        .{ -0x1.5fffffffffffffffffbcc71a492p+4, -0x1.eb5aeb96c74d7ad25e060528fb5p-112 },
        .{ -0x1.6000000000000000004338e5b6ep+4, 0x1.eb5aec04b2f2eb663e4e3d8a018cp-112 },
        .{ -0x1.6ffffffffffffffffffd13c97d9dp+4, -0x3.8fcc4d08d6fe5aa56ab04307ce7ep-112 },
        .{ -0x1.70000000000000000002ec368263p+4, 0x3.8fcc4d090cee2f5d0b69a99c353cp-112 },
        .{ -0x1.7fffffffffffffffffffe0d30fe7p+4, 0x7.2f577cca4b4c8cb1dc14001ac5ecp-112 },
        .{ -0x1.800000000000000000001f2cf019p+4, -0x7.2f577cca4b3442e35f0040b3b9e8p-112 },
        .{ -0x1.8ffffffffffffffffffffec0c332p+4, -0x2.e9a0572b1bb5b95f346a92d67a6p-112 },
        .{ -0x1.90000000000000000000013f3ccep+4, 0x2.e9a0572b1bb5c371ddb3561705ap-112 },
        .{ -0x1.9ffffffffffffffffffffff3b8bdp+4, -0x1.cad8d32e386fd783e97296d63dcbp-116 },
        .{ -0x1.a0000000000000000000000c4743p+4, 0x1.cad8d32e386fd7c1ab8c1fe34c0ep-116 },
        .{ -0x1.afffffffffffffffffffffff8b95p+4, -0x3.8f48cc5737d5979c39db806c5406p-112 },
        .{ -0x1.b00000000000000000000000746bp+4, 0x3.8f48cc5737d5979c3b3a6bda06f6p-112 },
        .{ -0x1.bffffffffffffffffffffffffbd8p+4, 0x6.2898d42174dcf171470d8c8c6028p-112 },
        .{ -0x1.c000000000000000000000000428p+4, -0x6.2898d42174dcf171470d18ba412cp-112 },
        .{ -0x1.cfffffffffffffffffffffffffdbp+4, -0x4.c0ce9794ea50a839e311320bde94p-112 },
        .{ -0x1.d000000000000000000000000025p+4, 0x4.c0ce9794ea50a839e311322f7cf8p-112 },
        .{ -0x1.dfffffffffffffffffffffffffffp+4, 0x3.932c5047d60e60caded4c298a174p-112 },
        .{ -0x1.e000000000000000000000000001p+4, -0x3.932c5047d60e60caded4c298973ap-112 },
        .{ -0x1.fp+4, 0xa.1a6973c1fade2170f7237d35fe3p-116 },
        .{ -0x1.fp+4, -0xa.1a6973c1fade2170f7237d35fe08p-116 },
        .{ -0x2p+4, 0x5.0d34b9e0fd6f10b87b91be9aff1p-120 },
        .{ -0x2p+4, -0x5.0d34b9e0fd6f10b87b91be9aff0cp-120 },
        .{ -0x2.1p+4, 0x2.73024a9ba1aa36a7059bff52e844p-124 },
        .{ -0x2.1p+4, -0x2.73024a9ba1aa36a7059bff52e844p-124 },
        .{ -0x2.2p+4, 0x1.2710231c0fd7a13f8a2b4af9d6b7p-128 },
        .{ -0x2.2p+4, -0x1.2710231c0fd7a13f8a2b4af9d6b7p-128 },
        .{ -0x2.3p+4, 0x8.6e2ce38b6c8f9419e3fad3f0312p-136 },
        .{ -0x2.3p+4, -0x8.6e2ce38b6c8f9419e3fad3f0312p-136 },
        .{ -0x2.4p+4, 0x3.bf30652185952560d71a254e4eb8p-140 },
        .{ -0x2.4p+4, -0x3.bf30652185952560d71a254e4eb8p-140 },
        .{ -0x2.5p+4, 0x1.9ec8d1c94e85af4c78b15c3d89d3p-144 },
        .{ -0x2.5p+4, -0x1.9ec8d1c94e85af4c78b15c3d89d3p-144 },
        .{ -0x2.6p+4, 0xa.ea565ce061d57489e9b85276274p-152 },
        .{ -0x2.6p+4, -0xa.ea565ce061d57489e9b85276274p-152 },
        .{ -0x2.7p+4, 0x4.7a6512692eb37804111dabad30ecp-156 },
        .{ -0x2.7p+4, -0x4.7a6512692eb37804111dabad30ecp-156 },
        .{ -0x2.8p+4, 0x1.ca8ed42a12ae3001a07244abad2bp-160 },
        .{ -0x2.8p+4, -0x1.ca8ed42a12ae3001a07244abad2bp-160 },
        .{ -0x2.9p+4, 0xb.2f30e1ce812063f12e7e8d8d96e8p-168 },
        .{ -0x2.9p+4, -0xb.2f30e1ce812063f12e7e8d8d96e8p-168 },
        .{ -0x2.ap+4, 0x4.42bd49d4c37a0db136489772e428p-172 },
        .{ -0x2.ap+4, -0x4.42bd49d4c37a0db136489772e428p-172 },
        .{ -0x2.bp+4, 0x1.95db45257e5122dcbae56def372p-176 },
        .{ -0x2.bp+4, -0x1.95db45257e5122dcbae56def372p-176 },
        .{ -0x2.cp+4, 0x9.3958d81ff63527ecf993f3fb6f48p-184 },
        .{ -0x2.cp+4, -0x9.3958d81ff63527ecf993f3fb6f48p-184 },
        .{ -0x2.dp+4, 0x3.47970e4440c8f1c058bd238c9958p-188 },
        .{ -0x2.dp+4, -0x3.47970e4440c8f1c058bd238c9958p-188 },
        .{ -0x2.ep+4, 0x1.240804f65951062ca46e4f25c608p-192 },
        .{ -0x2.ep+4, -0x1.240804f65951062ca46e4f25c608p-192 },
        .{ -0x2.fp+4, 0x6.36a382849fae6de2d15362d8a394p-200 },
        .{ -0x2.fp+4, -0x6.36a382849fae6de2d15362d8a394p-200 },
        .{ -0x3p+4, 0x2.123680d6dfe4cf4b9b1bcb9d8bdcp-204 },
        .{ -0x3p+4, -0x2.123680d6dfe4cf4b9b1bcb9d8bdcp-204 },
        .{ -0x3.1p+4, 0xa.d21786ff5842eca51fea0870919p-212 },
        .{ -0x3.1p+4, -0xa.d21786ff5842eca51fea0870919p-212 },
        .{ -0x3.2p+4, 0x3.766dedc259af040be140a68a6c04p-216 },
    };

    const e_hi: f128 = 0x2.b7e151628aed2a6abf7158809cf4p+0;
    const e_lo: f128 = 0xf.3c762e7160f38b4da56a784d9048p-116;

    // Coefficients B_2k / 2k(2k-1) of x^-(2k-1) in Stirling's
    // approximation to lgamma function.
    const lgamma_coeff: [27]f128 = .{
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
        0xa.8d1044d3708d1c219ee4fdc446ap+16,
        -0xe.8844d8a169abbc406169abbc406p+20,
        0x1.6d29a0f6433b79890cede62433b8p+28,
        -0x2.88a233b3c8cddaba9809357125d8p+32,
        0x5.0dde6f27500939a85c40939a85c4p+36,
        -0xb.4005bde03d4642a243581714af68p+40,
        0x1.bc8cd6f8f1f755c78753cdb5d5c9p+48,
        -0x4.bbebb143bb94de5a0284fa7ec424p+52,
        0xe.2e1337f5af0bed90b6b0a352d4fp+56,
        -0x2.e78250162b62405ad3e4bfe61b38p+64,
        0xa.5f7eef9e71ac7c80326ab4cc8bfp+68,
        -0x2.83be0395e550213369924971b21ap+76,
        0xa.8ebfe48da17dd999790760b0cep+80,
    };

    // Polynomial approximations to (|gamma(x)|-1)(x-n)/(x-x0), where n is
    // the integer end-point of the half-integer interval containing x and
    // x0 is the zero of lgamma in that half-integer interval.  Each
    // polynomial is expressed in terms of x-xm, where xm is the midpoint
    // of the interval for which the polynomial applies.
    const poly_coeff: [208]f128 = .{
        // Interval [-2.125, -2] (polynomial degree 23).
        -0x1.0b71c5c54d42eb6c17f30b7aa8f5p+0,
        -0xc.73a1dc05f34951602554c6d7506p-4,
        -0x1.ec841408528b51473e6c425ee5ffp-4,
        -0xe.37c9da26fc3c9a3c1844c8c7f1cp-4,
        -0x1.03cd87c519305703b021fa33f827p-4,
        -0xe.ae9ada65e09aa7f1c75216128f58p-4,
        0x9.b11855a4864b5731cf85736015a8p-8,
        -0xe.f28c133e697a95c28607c9701dep-4,
        0x2.6ec14a1c586a72a7cc33ee569d6ap-4,
        -0xf.57cab973e14464a262fc24723c38p-4,
        0x4.5b0fc25f16e52997b2886bbae808p-4,
        -0xf.f50e59f1a9b56e76e988dac9ccf8p-4,
        0x6.5f5eae15e9a93369e1d85146c6fcp-4,
        -0x1.0d2422daac459e33e0994325ed23p+0,
        0x8.82000a0e7401fb1117a0e6606928p-4,
        -0x1.1f492f178a3f1b19f58a2ca68e55p+0,
        0xa.cb545f949899a04c160b19389abp-4,
        -0x1.36165a1b155ba3db3d1b77caf498p+0,
        0xd.44c5d5576f74302e5cf79e183eep-4,
        -0x1.51f22e0cdd33d3d481e326c02f3ep+0,
        0xf.f73a349c08244ac389c007779bfp-4,
        -0x1.73317bf626156ba716747c4ca866p+0,
        0x1.379c3c97b9bc71e1c1c4802dd657p+0,
        -0x1.a72a351c54f902d483052000f5dfp+0,
        // Interval [-2.25, -2.125] (polynomial degree 24).
        -0xf.2930890d7d675a80c36afb0fd5e8p-4,
        -0xc.a5cfde054eab5c6770daeca577f8p-4,
        0x3.9c9e0fdebb07cdf89c61d41c9238p-4,
        -0x1.02a5ad35605fcf4af65a6dbacb84p+0,
        0x9.6e9b1185bb48be9de1918e00a2e8p-4,
        -0x1.4d8332f3cfbfa116fd611e9ce90dp+0,
        0x1.1c0c8cb4d9f4b1d490e1a41fae4dp+0,
        -0x1.c9a6f5ae9130cd0299e293a42714p+0,
        0x1.d7e9307fd58a2ea997f29573a112p+0,
        -0x2.921cb3473d96178ca2a11d2a8d46p+0,
        0x2.e8d59113b6f3409ff8db226e9988p+0,
        -0x3.cbab931625a1ae2b26756817f264p+0,
        0x4.7d9f0f05d5296d18663ca003912p+0,
        -0x5.ade9cba12a14ea485667b7135bbp+0,
        0x6.dc983a5da74fb48e767b7fec0a3p+0,
        -0x8.8d9ed454ae31d9e138dd8ee0d1a8p+0,
        0xa.6fa099d4e7c202e0c0fd6ed8492p+0,
        -0xc.ebc552a8090a0f0115e92d4ebbc8p+0,
        0xf.d695e4772c0d829b53fba9ca5568p+0,
        -0x1.38c32ae38e5e9eb79b2a4c5570a9p+4,
        0x1.8035145646cfab49306d0999a51bp+4,
        -0x1.d930adbb03dd342a4c2a8c4e1af6p+4,
        0x2.45c2edb1b4943ddb3686cd9c6524p+4,
        -0x2.e818ebbfafe2f916fa21abf7756p+4,
        0x3.9804ce51d0fb9a430a711fd7307p+4,
        // Interval [-2.375, -2.25] (polynomial degree 25).
        -0xd.7d28d505d6181218a25f31d5e45p-4,
        -0xe.69649a3040985140cdf946829fap-4,
        0xb.0d74a2827d053a8d44595012484p-4,
        -0x1.924b0922853617cac181afbc08ddp+0,
        0x1.d49b12bccf0a568582e2d3c410f3p+0,
        -0x3.0898bb7d8c4093e636279c791244p+0,
        0x4.207a6cac711cb53868e8a5057eep+0,
        -0x6.39ee63ea4fb1dcab0c9144bf3ddcp+0,
        0x8.e2e2556a797b649bf3f53bd26718p+0,
        -0xd.0e83ac82552ef12af508589e7a8p+0,
        0x1.2e4525e0ce6670563c6484a82b05p+4,
        -0x1.b8e350d6a8f2b222fa390a57c23dp+4,
        0x2.805cd69b919087d8a80295892c2cp+4,
        -0x3.a42585424a1b7e64c71743ab014p+4,
        0x5.4b4f409f98de49f7bfb03c05f984p+4,
        -0x7.b3c5827fbe934bc820d6832fb9fcp+4,
        0xb.33b7b90cc96c425526e0d0866e7p+4,
        -0x1.04b77047ac4f59ee3775ca10df0dp+8,
        0x1.7b366f5e94a34f41386eac086313p+8,
        -0x2.2797338429385c9849ca6355bfc2p+8,
        0x3.225273cf92a27c9aac1b35511256p+8,
        -0x4.8f078aa48afe6cb3a4e89690f898p+8,
        0x6.9f311d7b6654fc1d0b5195141d04p+8,
        -0x9.a0c297b6b4621619ca9bacc48ed8p+8,
        0xe.ce1f06b6f90d92138232a76e4cap+8,
        -0x1.5b0e6806fa064daf011613e43b17p+12,
        // Interval [-2.5, -2.375] (polynomial degree 27).
        -0xb.74ea1bcfff94b2c01afba9daa7d8p-4,
        -0x1.2a82bd590c37538cab143308de4dp+0,
        0x1.88020f828b966fec66b8649fd6fcp+0,
        -0x3.32279f040eb694970e9db24863dcp+0,
        0x5.57ac82517767e68a721005853864p+0,
        -0x9.c2aedcfe22833de43834a0a6cc4p+0,
        0x1.12c132f1f5577f99e1a0ed3538e1p+4,
        -0x1.ea94e26628a3de3597f7bb55a948p+4,
        0x3.66b4ac4fa582f58b59f96b2f7c7p+4,
        -0x6.0cf746a9cf4cba8c39afcc73fc84p+4,
        0xa.c102ef2c20d75a342197df7fedf8p+4,
        -0x1.31ebff06e8f14626782df58db3b6p+8,
        0x2.1fd6f0c0e710994e059b9dbdb1fep+8,
        -0x3.c6d76040407f447f8b5074f07706p+8,
        0x6.b6d18e0d8feb4c2ef5af6a40ed18p+8,
        -0xb.efaf542c529f91e34217f24ae6a8p+8,
        0x1.53852d873210e7070f5d9eb2296p+12,
        -0x2.5b977c0ddc6d540717173ac29fc8p+12,
        0x4.310d452ae05100eff1e02343a724p+12,
        -0x7.73a5d8f20c4f986a7dd1912b2968p+12,
        0xd.3f5ea2484f3fca15eab1f4d1a218p+12,
        -0x1.78d18aac156d1d93a2ffe7e08d3fp+16,
        0x2.9df49ca75e5b567f5ea3e47106cp+16,
        -0x4.a7149af8961a08aa7c3233b5bb94p+16,
        0x8.3db10ffa742c707c25197d989798p+16,
        -0xe.a26d6dd023cadd02041a049ec368p+16,
        0x1.c825d90514e7c57c7fa5316f947cp+20,
        -0x3.34bb81e5a0952df8ca1abdc6684cp+20,
        // Interval [-2.625, -2.5] (polynomial degree 28).
        -0x3.d10108c27ebafad533c20eac32bp-4,
        0x1.cd557caff7d2b2085f41dbec5106p+0,
        0x3.819b4856d399520dad9776ea2cacp+0,
        0x6.8505cbad03dc34c5e42e8b12eb78p+0,
        0xb.c1b2e653a9e38f82b399c94e7f08p+0,
        0x1.50a53a38f148138105124df65419p+4,
        0x2.57ae00cbe5232cbeeed34d89727ap+4,
        0x4.2b156301b8604db85a601544bfp+4,
        0x7.6989ed23ca3ca7579b3462592b5cp+4,
        0xd.2dd2976557939517f831f5552cc8p+4,
        0x1.76e1c3430eb860969bce40cd494p+8,
        0x2.9a77bf5488742466db3a2c7c1ec6p+8,
        0x4.a0d62ed7266e8eb36f725a8ebcep+8,
        0x8.3a6184dd3021067df2f8b91e99c8p+8,
        0xe.a0ade1538245bf55d39d7e436b1p+8,
        0x1.a01359fae8617b5826dd74428e9p+12,
        0x2.e3b0a32caae77251169acaca1ad4p+12,
        0x5.2301257c81589f62b38fb5993ee8p+12,
        0x9.21c9275db253d4e719b73b18cb9p+12,
        0x1.03c104bc96141cda3f3fa4b112bcp+16,
        0x1.cdc8ed65119196a08b0c78f1445p+16,
        0x3.34f31d2eaacf34382cdb0073572ap+16,
        0x5.b37628cadf12bf0000907d0ef294p+16,
        0xa.22d8b332c0b1e6a616f425dfe5ap+16,
        0x1.205b01444804c3ff922cd78b4c42p+20,
        0x1.fe8f0cea9d1e0ff25be2470b4318p+20,
        0x3.8872aebeb368399aee02b39340aep+20,
        0x6.ebd560d351e84e26a4381f5b293cp+20,
        0xc.c3644d094b0dae2fbcbf682cd428p+20,
        // Interval [-2.75, -2.625] (polynomial degree 26).
        -0x6.b5d252a56e8a75458a27ed1c2dd4p-4,
        0x1.28d60383da3ac721aed3c5794da9p+0,
        0x1.db6513ada8a66ea77d87d9a8827bp+0,
        0x2.e217118f9d348a27f7506a707e6ep+0,
        0x4.450112c5cbf725a0fb9802396c9p+0,
        0x6.4af99151eae7810a75df2a0303c4p+0,
        0x9.2db598b4a97a7f69aeef32aec758p+0,
        0xd.62bef9c22471f5ee47ea1b9c0b5p+0,
        0x1.379f294e412bd62328326d4222f9p+4,
        0x1.c5827349d8865f1e8825c37c31c6p+4,
        0x2.93a7e7a75b7568cc8cbe8c016c12p+4,
        0x3.bf9bb882afe57edb383d41879d3ap+4,
        0x5.73c737828cee095c43a5566731c8p+4,
        0x7.ee4653493a7f81e0442062b3823cp+4,
        0xb.891c6b83fc8b55bd973b5d962d6p+4,
        0x1.0c775d7de3bf9b246c0208e0207ep+8,
        0x1.867ee43ec4bd4f4fd56abc05110ap+8,
        0x2.37fe9ba6695821e9822d8c8af0a6p+8,
        0x3.3a2c667e37c942f182cd3223a936p+8,
        0x4.b1b500eb59f3f782c7ccec88754p+8,
        0x6.d3efd3b65b3d0d8488d30b79fa4cp+8,
        0x9.ee8224e65bed5ced8b75eaec609p+8,
        0xe.72416e510cca77d53fc615c1f3dp+8,
        0x1.4fb538b0a2dfe567a8904b7e0445p+12,
        0x1.e7f56a9266cf525a5b8cf4cb76cep+12,
        0x2.f0365c983f68c597ee49d099cce8p+12,
        0x4.53aa229e1b9f5b5e59625265951p+12,
        // Interval [-2.875, -2.75] (polynomial degree 24).
        -0x8.a41b1e4f36ff88dc820815607d68p-4,
        0xc.da87d3b69dc0f2f9c6f368b8ca1p-4,
        0x1.1474ad5c36158a7bea04fd2f98c6p+0,
        0x1.761ecb90c555df6555b7dba955b6p+0,
        0x1.d279bff9ae291caf6c4b4bcb3202p+0,
        0x2.4e5d00559a6e2b9b5d7fe1f6689cp+0,
        0x2.d57545a75cee8743ae2b17bc8d24p+0,
        0x3.8514eee3aac88b89bec2307021bap+0,
        0x4.5235e3b6e1891ffeb87fed9f8a24p+0,
        0x5.562acdb10eef3c9a773b3e27a864p+0,
        0x6.8ec8965c76efe03c26bff60b1194p+0,
        0x8.15251aca144877af32658399f9b8p+0,
        0x9.f08d56aba174d844138af782c0f8p+0,
        0xc.3dbbeda2679e8a1346ccc3f6da88p+0,
        0xf.0f5bfd5eacc26db308ffa0556fa8p+0,
        0x1.28a6ccd84476fbc713d6bab49ac9p+4,
        0x1.6d0a3ae2a3b1c8ff400641a3a21fp+4,
        0x1.c15701b28637f87acfb6a91d33b5p+4,
        0x2.28fbe0eccf472089b017651ca55ep+4,
        0x2.a8a453004f6e8ffaacd1603bc3dp+4,
        0x3.45ae4d9e1e7cd1a5dba0e4ec7f6cp+4,
        0x4.065fbfacb7fad3e473cb577a61e8p+4,
        0x4.f3d1473020927acac1944734a39p+4,
        0x6.54bb091245815a36fb74e314dd18p+4,
        0x7.d7f445129f7fb6c055e582d3f6ep+4,
        // Interval [-3, -2.875] (polynomial degree 23).
        -0xa.046d667e468f3e44dcae1afcc648p-4,
        0x9.70b88dcc006c214d8d996fdf5ccp-4,
        0xa.a8a39421c86d3ff24931a0929fp-4,
        0xd.2f4d1363f324da2b357c8b6ec94p-4,
        0xd.ca9aa1a3a5c00de11bf60499a97p-4,
        0xf.cf09c31eeb52a45dfa7ebe3778dp-4,
        0x1.04b133a39ed8a09691205660468bp+0,
        0x1.22b547a06edda944fcb12fd9b5ecp+0,
        0x1.2c57fce7db86a91df09602d344b3p+0,
        0x1.4aade4894708f84795212fe257eep+0,
        0x1.579c8b7b67ec4afed5b28c8bf787p+0,
        0x1.776820e7fc80ae5284239733078ap+0,
        0x1.883ab28c7301fde4ca6b8ec26ec8p+0,
        0x1.aa2ef6e1ae52eb42c9ee83b206e3p+0,
        0x1.bf4ad50f0a9a9311300cf0c51ee7p+0,
        0x1.e40206e0e96b1da463814dde0d09p+0,
        0x1.fdcbcffef3a21b29719c2bd9feb1p+0,
        0x2.25e2e8948939c4d42cf108fae4bep+0,
        0x2.44ce14d2b59c1c0e6bf2cfa81018p+0,
        0x2.70ee80bbd0387162be4861c43622p+0,
        0x2.954b64d2c2ebf3489b949c74476p+0,
        0x2.c616e133a811c1c9446105208656p+0,
        0x3.05a69dfe1a9ba1079f90fcf26bd4p+0,
        0x3.410d2ad16a0506de29736e6aafdap+0,
    };

    const poly_deg: [8]u32 = .{
        23,
        24,
        25,
        27,
        28,
        26,
        24,
        23,
    };

    const poly_end: [8]u32 = .{
        23,
        48,
        74,
        102,
        131,
        158,
        183,
        207,
    };

    // Determine the half-integer region X lies in, handle exact
    // integers and determine the sign of the result.
    var i: i32 = cast(i32, float.floor(-2 * x), .{});
    if ((i & 1) == 0 and cast(f128, i, .{}) == -2 * x)
        return 1 / @as(f128, 0);

    const xn: f128 = cast(f128, if ((i & 1) == 0) @divFloor(-i, 2) else @divFloor(-i - 1, 2), .{});
    i -= 4;
    signgamp.* = if ((i & 2) == 0) -1 else 1;

    // Expand around the zero X0 = X0_HI + X0_LO.
    const x0_hi: f128 = lgamma_zeros[@intCast(i)][0];
    const x0_lo: f128 = lgamma_zeros[@intCast(i)][1];
    const xdiff: f128 = x - x0_hi - x0_lo;

    // For arguments in the range -3 to -2, use polynomial
    // approximations to an adjusted version of the gamma function.
    if (i < 2) {
        var j: i32 = cast(i32, float.floor(-8 * x) - 16, .{});
        const xm: f128 = (-33 - 2 * cast(f128, j, .{})) * 0.0625;
        const x_adj: f128 = x - xm;
        const deg: u32 = poly_deg[@intCast(j)];
        const end: u32 = poly_end[@intCast(j)];
        var g: f128 = poly_coeff[end];
        j = 1;
        while (j <= deg) {
            g = g * x_adj + poly_coeff[@intCast(cast(i32, end, .{}) - j)];

            j += 1;
        }

        return float.log1p(g * xdiff / (x - xn));
    }

    // The result we want is log (sinpi (X0) / sinpi (X))
    //  + log (gamma (1 - X0) / gamma (1 - X)).
    const x_idiff: f128 = float.abs(xn - x);
    const x0_idiff: f128 = float.abs(xn - x0_hi - x0_lo);
    var log_sinpi_ratio: f128 = undefined;
    if (x0_idiff < x_idiff * 0.5) {
        // Use log not log1p to avoid inaccuracy from log1p of arguments
        // close to -1.
        log_sinpi_ratio = float.log(lg_sinpi128(x0_idiff) / lg_sinpi128(x_idiff));
    } else {
        // Use log1p not log to avoid inaccuracy from log of arguments
        // close to 1.  X0DIFF2 has positive sign if X0 is further from
        // XN than X is from XN, negative sign otherwise.
        const x0diff2 = (if ((i & 1) == 0) xdiff else -xdiff) * 0.5;
        const sx0d2 = lg_sinpi128(x0diff2);
        const cx0d2 = lg_cospi128(x0diff2);
        log_sinpi_ratio = float.log1p(2 * sx0d2 * (-sx0d2 + cx0d2 * lg_cotpi128(x_idiff)));
    }

    var y0: f128 = 1 - x0_hi;
    var y0_eps: f128 = -x0_hi + (1 - y0) - x0_lo;
    var y: f128 = 1 - x;
    var y_eps: f128 = -x + (1 - y);
    // We now wish to compute LOG_GAMMA_RATIO
    //  = log (gamma (Y0 + Y0_EPS) / gamma (Y + Y_EPS)).  XDIFF
    // accurately approximates the difference Y0 + Y0_EPS - Y -
    // Y_EPS.  Use Stirling's approximation.  First, we may need to
    // adjust into the range where Stirling's approximation is
    // sufficiently accurate.
    var log_gamma_adj: f128 = 0;
    if (i < 20) {
        const n_up: i32 = @divFloor(21 - i, 2);
        const ny0: f128 = y0 + cast(f128, n_up, .{});
        const ny0_eps: f128 = y0 - (ny0 - cast(f128, n_up, .{})) + y0_eps;
        y0 = ny0;
        y0_eps = ny0_eps;
        const ny: f128 = y + cast(f128, n_up, .{});
        const ny_eps: f128 = y - (ny - cast(f128, n_up, .{})) + y_eps;
        y = ny;
        y_eps = ny_eps;
        const prodm1: f128 = lgamma_product128(xdiff, y - cast(f128, n_up, .{}), y_eps, n_up);
        log_gamma_adj = -float.log1p(prodm1);
    }

    const log_gamma_high: f128 = (xdiff * float.log1p((y0 - e_hi - e_lo + y0_eps) / e_hi) + (y - 0.5 + y_eps) * float.log1p(xdiff / y) + log_gamma_adj);
    // Compute the sum of (B_2k / 2k(2k-1))(Y0^-(2k-1) - Y^-(2k-1)).
    const y0r: f128 = 1 / y0;
    const yr: f128 = 1 / y;
    const y0r2: f128 = y0r * y0r;
    const yr2: f128 = yr * yr;
    const rdiff: f128 = -xdiff / (y * y0);
    var bterm: [27]f128 = undefined;
    var dlast: f128 = rdiff;
    var elast: f128 = rdiff * yr * (yr + y0r);
    bterm[0] = dlast * lgamma_coeff[0];
    var j: u32 = 1;
    while (j < 27) {
        const dnext: f128 = dlast * y0r2 + elast;
        const enext: f128 = elast * yr2;
        bterm[j] = dnext * lgamma_coeff[j];
        dlast = dnext;
        elast = enext;

        j += 1;
    }

    var log_gamma_low: f128 = 0;
    j = 0;
    while (j < 27) {
        log_gamma_low += bterm[26 - j];

        j += 1;
    }

    const log_gamma_ratio: f128 = log_gamma_high + log_gamma_low;
    return log_sinpi_ratio + log_gamma_ratio;
}

pub fn lgamma_r128(x: f128, signgamp: *i32) f128 {
    const MAXLGM: f128 = 1.0485738685148938358098967157129705071571e4928;
    const huge: f128 = std.math.floatMax(f128);
    // log gamma(x) = ( x - 0.5 ) * log(x) - x + LS2PI + 1/x P(1/x^2)
    // 1/x <= 0.0741 (x >= 13.495...)
    // Peak relative error 1.5e-36
    const ls2pi: f128 = 9.1893853320467274178032973640561763986140e-1;
    const RASY: [13]f128 = .{
        8.333333333333333333333333333310437112111e-2,
        -2.777777777777777777777774789556228296902e-3,
        7.936507936507936507795933938448586499183e-4,
        -5.952380952380952041799269756378148574045e-4,
        8.417508417507928904209891117498524452523e-4,
        -1.917526917481263997778542329739806086290e-3,
        6.410256381217852504446848671499409919280e-3,
        -2.955064066900961649768101034477363301626e-2,
        1.796402955865634243663453415388336954675e-1,
        -1.391522089007758553455753477688592767741e0,
        1.326130089598399157988112385013829305510e1,
        -1.420412699593782497803472576479997819149e2,
        1.218058922427762808938869872528846787020e3,
    };
    // log gamma(x+13) = log gamma(13) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 12.5 <= x+13 <= 13.5
    // Peak relative error 1.1e-36
    const lgam13a: f128 = 1.9987213134765625e1;
    const lgam13b: f128 = 1.3608962611495173623870550785125024484248e-6;
    const RN13: [8]f128 = .{
        8.591478354823578150238226576156275285700e11,
        2.347931159756482741018258864137297157668e11,
        2.555408396679352028680662433943000804616e10,
        1.408581709264464345480765758902967123937e9,
        4.126759849752613822953004114044451046321e7,
        6.133298899622688505854211579222889943778e5,
        3.929248056293651597987893340755876578072e3,
        6.850783280018706668924952057996075215223e0,
    };
    const RD13: [7]f128 = .{
        3.401225382297342302296607039352935541669e11,
        8.756765276918037910363513243563234551784e10,
        8.873913342866613213078554180987647243903e9,
        4.483797255342763263361893016049310017973e8,
        1.178186288833066430952276702931512870676e7,
        1.519928623743264797939103740132278337476e5,
        7.989298844938119228411117593338850892311e2,
    };
    // log gamma(x+12) = log gamma(12) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 11.5 <= x+12 <= 12.5
    // Peak relative error 4.1e-36
    const lgam12a: f128 = 1.75023040771484375e1;
    const lgam12b: f128 = 3.7687254483392876529072161996717039575982e-6;
    const RN12: [8]f128 = .{
        4.709859662695606986110997348630997559137e11,
        1.398713878079497115037857470168777995230e11,
        1.654654931821564315970930093932954900867e10,
        9.916279414876676861193649489207282144036e8,
        3.159604070526036074112008954113411389879e7,
        5.109099197547205212294747623977502492861e5,
        3.563054878276102790183396740969279826988e3,
        6.769610657004672719224614163196946862747e0,
    };
    const RD12: [7]f128 = .{
        1.928167007860968063912467318985802726613e11,
        5.383198282277806237247492369072266389233e10,
        5.915693215338294477444809323037871058363e9,
        3.241438287570196713148310560147925781342e8,
        9.236680081763754597872713592701048455890e6,
        1.292246897881650919242713651166596478850e5,
        7.366532445427159272584194816076600211171e2,
    };
    // log gamma(x+11) = log gamma(11) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 10.5 <= x+11 <= 11.5
    // Peak relative error 1.8e-35
    const lgam11a: f128 = 1.5104400634765625e1;
    const lgam11b: f128 = 1.1938309890295225709329251070371882250744e-5;
    const RN11: [8]f128 = .{
        2.446960438029415837384622675816736622795e11,
        7.955444974446413315803799763901729640350e10,
        1.030555327949159293591618473447420338444e10,
        6.765022131195302709153994345470493334946e8,
        2.361892792609204855279723576041468347494e7,
        4.186623629779479136428005806072176490125e5,
        3.202506022088912768601325534149383594049e3,
        6.681356101133728289358838690666225691363e0,
    };
    const RD11: [7]f128 = .{
        1.040483786179428590683912396379079477432e11,
        3.172251138489229497223696648369823779729e10,
        3.806961885984850433709295832245848084614e9,
        2.278070344022934913730015420611609620171e8,
        7.089478198662651683977290023829391596481e6,
        1.083246385105903533237139380509590158658e5,
        6.744420991491385145885727942219463243597e2,
    };
    // log gamma(x+10) = log gamma(10) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 9.5 <= x+10 <= 10.5
    // Peak relative error 5.4e-37
    const lgam10a: f128 = 1.280181884765625e1;
    const lgam10b: f128 = 8.6324252196112077178745667061642811492557e-6;
    const RN10: [8]f128 = .{
        -1.239059737177249934158597996648808363783e14,
        -4.725899566371458992365624673357356908719e13,
        -7.283906268647083312042059082837754850808e12,
        -5.802855515464011422171165179767478794637e11,
        -2.532349691157548788382820303182745897298e10,
        -5.884260178023777312587193693477072061820e8,
        -6.437774864512125749845840472131829114906e6,
        -2.350975266781548931856017239843273049384e4,
    };
    const RD10: [8]f128 = .{
        -5.502645997581822567468347817182347679552e13,
        -1.970266640239849804162284805400136473801e13,
        -2.819677689615038489384974042561531409392e12,
        -2.056105863694742752589691183194061265094e11,
        -8.053670086493258693186307810815819662078e9,
        -1.632090155573373286153427982504851867131e8,
        -1.483575879240631280658077826889223634921e6,
        -4.002806669713232271615885826373550502510e3,
    };
    // log gamma(x+9) = log gamma(9) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 8.5 <= x+9 <= 9.5
    // Peak relative error 3.6e-36
    const lgam9a: f128 = 1.06045989990234375e1;
    const lgam9b: f128 = 3.9037218127284172274007216547549861681400e-6;
    const RN9: [8]f128 = .{
        -4.936332264202687973364500998984608306189e13,
        -2.101372682623700967335206138517766274855e13,
        -3.615893404644823888655732817505129444195e12,
        -3.217104993800878891194322691860075472926e11,
        -1.568465330337375725685439173603032921399e10,
        -4.073317518162025744377629219101510217761e8,
        -4.983232096406156139324846656819246974500e6,
        -2.036280038903695980912289722995505277253e4,
    };
    const RD9: [8]f128 = .{
        -2.306006080437656357167128541231915480393e13,
        -9.183606842453274924895648863832233799950e12,
        -1.461857965935942962087907301194381010380e12,
        -1.185728254682789754150068652663124298303e11,
        -5.166285094703468567389566085480783070037e9,
        -1.164573656694603024184768200787835094317e8,
        -1.177343939483908678474886454113163527909e6,
        -3.529391059783109732159524500029157638736e3,
    };
    // log gamma(x+8) = log gamma(8) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 7.5 <= x+8 <= 8.5
    // Peak relative error 2.4e-37
    const lgam8a: f128 = 8.525146484375e0;
    const lgam8b: f128 = 1.4876690414300165531036347125050759667737e-5;
    const RN8: [9]f128 = .{
        6.600775438203423546565361176829139703289e11,
        3.406361267593790705240802723914281025800e11,
        7.222460928505293914746983300555538432830e10,
        8.102984106025088123058747466840656458342e9,
        5.157620015986282905232150979772409345927e8,
        1.851445288272645829028129389609068641517e7,
        3.489261702223124354745894067468953756656e5,
        2.892095396706665774434217489775617756014e3,
        6.596977510622195827183948478627058738034e0,
    };
    const RD8: [8]f128 = .{
        3.274776546520735414638114828622673016920e11,
        1.581811207929065544043963828487733970107e11,
        3.108725655667825188135393076860104546416e10,
        3.193055010502912617128480163681842165730e9,
        1.830871482669835106357529710116211541839e8,
        5.790862854275238129848491555068073485086e6,
        9.305213264307921522842678835618803553589e4,
        6.216974105861848386918949336819572333622e2,
    };
    // log gamma(x+7) = log gamma(7) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 6.5 <= x+7 <= 7.5
    // Peak relative error 3.2e-36
    const lgam7a: f128 = 6.5792388916015625e0;
    const lgam7b: f128 = 1.2320408538495060178292903945321122583007e-5;
    const RN7: [9]f128 = .{
        2.065019306969459407636744543358209942213e11,
        1.226919919023736909889724951708796532847e11,
        2.996157990374348596472241776917953749106e10,
        3.873001919306801037344727168434909521030e9,
        2.841575255593761593270885753992732145094e8,
        1.176342515359431913664715324652399565551e7,
        2.558097039684188723597519300356028511547e5,
        2.448525238332609439023786244782810774702e3,
        6.460280377802030953041566617300902020435e0,
    };
    const RD7: [8]f128 = .{
        1.102646614598516998880874785339049304483e11,
        6.099297512712715445879759589407189290040e10,
        1.372898136289611312713283201112060238351e10,
        1.615306270420293159907951633566635172343e9,
        1.061114435798489135996614242842561967459e8,
        3.845638971184305248268608902030718674691e6,
        7.081730675423444975703917836972720495507e4,
        5.423122582741398226693137276201344096370e2,
    };
    // log gamma(x+6) = log gamma(6) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 5.5 <= x+6 <= 6.5
    // Peak relative error 6.2e-37
    const lgam6a: f128 = 4.7874908447265625e0;
    const lgam6b: f128 = 8.9805548349424770093452324304839959231517e-7;
    const RN6: [9]f128 = .{
        -3.538412754670746879119162116819571823643e13,
        -2.613432593406849155765698121483394257148e13,
        -8.020670732770461579558867891923784753062e12,
        -1.322227822931250045347591780332435433420e12,
        -1.262809382777272476572558806855377129513e11,
        -7.015006277027660872284922325741197022467e9,
        -2.149320689089020841076532186783055727299e8,
        -3.167210585700002703820077565539658995316e6,
        -1.576834867378554185210279285358586385266e4,
    };
    const RD6: [9]f128 = .{
        -2.073955870771283609792355579558899389085e13,
        -1.421592856111673959642750863283919318175e13,
        -4.012134994918353924219048850264207074949e12,
        -6.013361045800992316498238470888523722431e11,
        -5.145382510136622274784240527039643430628e10,
        -2.510575820013409711678540476918249524123e9,
        -6.564058379709759600836745035871373240904e7,
        -7.861511116647120540275354855221373571536e5,
        -2.821943442729620524365661338459579270561e3,
    };
    // log gamma(x+5) = log gamma(5) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 4.5 <= x+5 <= 5.5
    // Peak relative error 3.4e-37
    const lgam5a: f128 = 3.17803955078125e0;
    const lgam5b: f128 = 1.4279566695619646941601297055408873990961e-5;
    const RN5: [10]f128 = .{
        2.010952885441805899580403215533972172098e11,
        1.916132681242540921354921906708215338584e11,
        7.679102403710581712903937970163206882492e10,
        1.680514903671382470108010973615268125169e10,
        2.181011222911537259440775283277711588410e9,
        1.705361119398837808244780667539728356096e8,
        7.792391565652481864976147945997033946360e6,
        1.910741381027985291688667214472560023819e5,
        2.088138241893612679762260077783794329559e3,
        6.330318119566998299106803922739066556550e0,
    };
    const RD5: [9]f128 = .{
        1.335189758138651840605141370223112376176e11,
        1.174130445739492885895466097516530211283e11,
        4.308006619274572338118732154886328519910e10,
        8.547402888692578655814445003283720677468e9,
        9.934628078575618309542580800421370730906e8,
        6.847107420092173812998096295422311820672e7,
        2.698552646016599923609773122139463150403e6,
        5.526516251532464176412113632726150253215e4,
        4.772343321713697385780533022595450486932e2,
    };
    // log gamma(x+4) = log gamma(4) +  x P(x)/Q(x)
    // -0.5 <= x <= 0.5
    // 3.5 <= x+4 <= 4.5
    // Peak relative error 6.7e-37
    const lgam4a: f128 = 1.791748046875e0;
    const lgam4b: f128 = 1.1422353055000812477358380702272722990692e-5;
    const RN4: [10]f128 = .{
        -1.026583408246155508572442242188887829208e13,
        -1.306476685384622809290193031208776258809e13,
        -7.051088602207062164232806511992978915508e12,
        -2.100849457735620004967624442027793656108e12,
        -3.767473790774546963588549871673843260569e11,
        -4.156387497364909963498394522336575984206e10,
        -2.764021460668011732047778992419118757746e9,
        -1.036617204107109779944986471142938641399e8,
        -1.895730886640349026257780896972598305443e6,
        -1.180509051468390914200720003907727988201e4,
    };
    const RD4: [10]f128 = .{
        -8.172669122056002077809119378047536240889e12,
        -9.477592426087986751343695251801814226960e12,
        -4.629448850139318158743900253637212801682e12,
        -1.237965465892012573255370078308035272942e12,
        -1.971624313506929845158062177061297598956e11,
        -1.905434843346570533229942397763361493610e10,
        -1.089409357680461419743730978512856675984e9,
        -3.416703082301143192939774401370222822430e7,
        -4.981791914177103793218433195857635265295e5,
        -2.192507743896742751483055798411231453733e3,
    };
    // log gamma(x+3) = log gamma(3) +  x P(x)/Q(x)
    // -0.25 <= x <= 0.5
    // 2.75 <= x+3 <= 3.5
    // Peak relative error 6.0e-37
    const lgam3a: f128 = 6.93145751953125e-1;
    const lgam3b: f128 = 1.4286068203094172321214581765680755001344e-6;
    const RN3: [10]f128 = .{
        -4.813901815114776281494823863935820876670e11,
        -8.425592975288250400493910291066881992620e11,
        -6.228685507402467503655405482985516909157e11,
        -2.531972054436786351403749276956707260499e11,
        -6.170200796658926701311867484296426831687e10,
        -9.211477458528156048231908798456365081135e9,
        -8.251806236175037114064561038908691305583e8,
        -4.147886355917831049939930101151160447495e7,
        -1.010851868928346082547075956946476932162e6,
        -8.333374463411801009783402800801201603736e3,
    };
    const RD3: [10]f128 = .{
        -5.216713843111675050627304523368029262450e11,
        -8.014292925418308759369583419234079164391e11,
        -5.180106858220030014546267824392678611990e11,
        -1.830406975497439003897734969120997840011e11,
        -3.845274631904879621945745960119924118925e10,
        -4.891033385370523863288908070309417710903e9,
        -3.670172254411328640353855768698287474282e8,
        -1.505316381525727713026364396635522516989e7,
        -2.856327162923716881454613540575964890347e5,
        -1.622140448015769906847567212766206894547e3,
    };
    // log gamma(x+2.5) = log gamma(2.5) +  x P(x)/Q(x)
    // -0.125 <= x <= 0.25
    // 2.375 <= x+2.5 <= 2.75
    const lgam2r5a: f128 = 2.8466796875e-1;
    const lgam2r5b: f128 = 1.4901722919159632494669682701924320137696e-5;
    const RN2r5: [9]f128 = .{
        -4.676454313888335499356699817678862233205e9,
        -9.361888347911187924389905984624216340639e9,
        -7.695353600835685037920815799526540237703e9,
        -3.364370100981509060441853085968900734521e9,
        -8.449902011848163568670361316804900559863e8,
        -1.225249050950801905108001246436783022179e8,
        -9.732972931077110161639900388121650470926e6,
        -3.695711763932153505623248207576425983573e5,
        -4.717341584067827676530426007495274711306e3,
    };
    const RD2r5: [9]f128 = .{
        -6.650657966618993679456019224416926875619e9,
        -1.099511409330635807899718829033488771623e10,
        -7.482546968307837168164311101447116903148e9,
        -2.702967190056506495988922973755870557217e9,
        -5.570008176482922704972943389590409280950e8,
        -6.536934032192792470926310043166993233231e7,
        -4.101991193844953082400035444146067511725e6,
        -1.174082735875715802334430481065526664020e5,
        -9.932840389994157592102947657277692978511e2,
    };
    // log gamma(x+2) = x P(x)/Q(x)
    // -0.125 <= x <= +0.375
    // 1.875 <= x+2 <= 2.375
    // Peak relative error 4.6e-36
    const RN2: [10]f128 = .{
        -3.716661929737318153526921358113793421524e9,
        -1.138816715030710406922819131397532331321e10,
        -1.421017419363526524544402598734013569950e10,
        -9.510432842542519665483662502132010331451e9,
        -3.747528562099410197957514973274474767329e9,
        -8.923565763363912474488712255317033616626e8,
        -1.261396653700237624185350402781338231697e8,
        -9.918402520255661797735331317081425749014e6,
        -3.753996255897143855113273724233104768831e5,
        -4.778761333044147141559311805999540765612e3,
    };
    const RD2: [10]f128 = .{
        -8.790916836764308497770359421351673950111e9,
        -2.023108608053212516399197678553737477486e10,
        -1.958067901852022239294231785363504458367e10,
        -1.035515043621003101254252481625188704529e10,
        -3.253884432621336737640841276619272224476e9,
        -6.186383531162456814954947669274235815544e8,
        -6.932557847749518463038934953605969951466e7,
        -4.240731768287359608773351626528479703758e6,
        -1.197343995089189188078944689846348116630e5,
        -1.004622911670588064824904487064114090920e3,
    };
    // log gamma(x+1.75) = log gamma(1.75) +  x P(x)/Q(x)
    // -0.125 <= x <= +0.125
    // 1.625 <= x+1.75 <= 1.875
    // Peak relative error 9.2e-37
    const lgam1r75a: f128 = -8.441162109375e-2;
    const lgam1r75b: f128 = 1.0500073264444042213965868602268256157604e-5;
    const RN1r75: [9]f128 = .{
        -5.221061693929833937710891646275798251513e7,
        -2.052466337474314812817883030472496436993e8,
        -2.952718275974940270675670705084125640069e8,
        -2.132294039648116684922965964126389017840e8,
        -8.554103077186505960591321962207519908489e7,
        -1.940250901348870867323943119132071960050e7,
        -2.379394147112756860769336400290402208435e6,
        -1.384060879999526222029386539622255797389e5,
        -2.698453601378319296159355612094598695530e3,
    };
    const RD1r75: [9]f128 = .{
        -2.109754689501705828789976311354395393605e8,
        -5.036651829232895725959911504899241062286e8,
        -4.954234699418689764943486770327295098084e8,
        -2.589558042412676610775157783898195339410e8,
        -7.731476117252958268044969614034776883031e7,
        -1.316721702252481296030801191240867486965e7,
        -1.201296501404876774861190604303728810836e6,
        -5.007966406976106636109459072523610273928e4,
        -6.155817990560743422008969155276229018209e2,
    };
    // log gamma(x+x0) = y0 +  x^2 P(x)/Q(x)
    // -0.0867 <= x <= +0.1634
    // 1.374932... <= x+x0 <= 1.625032...
    // Peak relative error 4.0e-36
    const x0a: f128 = 1.4616241455078125;
    const x0b: f128 = 7.9994605498412626595423257213002588621246e-6;
    const y0a: f128 = -1.21490478515625e-1;
    const y0b: f128 = 4.1879797753919044854428223084178486438269e-6;
    const RN1r5: [9]f128 = .{
        6.827103657233705798067415468881313128066e5,
        1.910041815932269464714909706705242148108e6,
        2.194344176925978377083808566251427771951e6,
        1.332921400100891472195055269688876427962e6,
        4.589080973377307211815655093824787123508e5,
        8.900334161263456942727083580232613796141e4,
        9.053840838306019753209127312097612455236e3,
        4.053367147553353374151852319743594873771e2,
        5.040631576303952022968949605613514584950e0,
    };
    const RD1r5: [9]f128 = .{
        1.411036368843183477558773688484699813355e6,
        4.378121767236251950226362443134306184849e6,
        5.682322855631723455425929877581697918168e6,
        3.999065731556977782435009349967042222375e6,
        1.653651390456781293163585493620758410333e6,
        4.067774359067489605179546964969435858311e5,
        5.741463295366557346748361781768833633256e4,
        4.226404539738182992856094681115746692030e3,
        1.316980975410327975566999780608618774469e2,
    };
    // log gamma(x+1.25) = log gamma(1.25) +  x P(x)/Q(x)
    // -.125 <= x <= +.125
    // 1.125 <= x+1.25 <= 1.375
    // Peak relative error = 4.9e-36
    const lgam1r25a: f128 = -9.82818603515625e-2;
    const lgam1r25b: f128 = 1.0023929749338536146197303364159774377296e-5;
    const RN1r25: [10]f128 = .{
        -9.054787275312026472896002240379580536760e4,
        -8.685076892989927640126560802094680794471e4,
        2.797898965448019916967849727279076547109e5,
        6.175520827134342734546868356396008898299e5,
        5.179626599589134831538516906517372619641e5,
        2.253076616239043944538380039205558242161e5,
        5.312653119599957228630544772499197307195e4,
        6.434329437514083776052669599834938898255e3,
        3.385414416983114598582554037612347549220e2,
        4.907821957946273805080625052510832015792e0,
    };
    const RD1r25: [9]f128 = .{
        3.980939377333448005389084785896660309000e5,
        1.429634893085231519692365775184490465542e6,
        2.145438946455476062850151428438668234336e6,
        1.743786661358280837020848127465970357893e6,
        8.316364251289743923178092656080441655273e5,
        2.355732939106812496699621491135458324294e5,
        3.822267399625696880571810137601310855419e4,
        3.228463206479133236028576845538387620856e3,
        1.152133170470059555646301189220117965514e2,
    };
    // log gamma(x + 1) = x P(x)/Q(x)
    // 0.0 <= x <= +0.125
    // 1.0 <= x+1 <= 1.125
    // Peak relative error 1.1e-35
    const RN1: [9]f128 = .{
        -9.987560186094800756471055681088744738818e3,
        -2.506039379419574361949680225279376329742e4,
        -1.386770737662176516403363873617457652991e4,
        1.439445846078103202928677244188837130744e4,
        2.159612048879650471489449668295139990693e4,
        1.047439813638144485276023138173676047079e4,
        2.250316398054332592560412486630769139961e3,
        1.958510425467720733041971651126443864041e2,
        4.516830313569454663374271993200291219855e0,
    };
    const RD1: [8]f128 = .{
        1.730299573175751778863269333703788214547e4,
        6.807080914851328611903744668028014678148e4,
        1.090071629101496938655806063184092302439e5,
        9.124354356415154289343303999616003884080e4,
        4.262071638655772404431164427024003253954e4,
        1.096981664067373953673982635805821283581e4,
        1.431229503796575892151252708527595787588e3,
        7.734110684303689320830401788262295992921e1,
    };
    // log gamma(x + 1) = x P(x)/Q(x)
    // -0.125 <= x <= 0
    // 0.875 <= x+1 <= 1.0
    // Peak relative error 7.0e-37
    const RNr9: [9]f128 = .{
        4.441379198241760069548832023257571176884e5,
        1.273072988367176540909122090089580368732e6,
        9.732422305818501557502584486510048387724e5,
        -5.040539994443998275271644292272870348684e5,
        -1.208719055525609446357448132109723786736e6,
        -7.434275365370936547146540554419058907156e5,
        -2.075642969983377738209203358199008185741e5,
        -2.565534860781128618589288075109372218042e4,
        -1.032901669542994124131223797515913955938e3,
    };
    const RDr9: [9]f128 = .{
        -7.694488331323118759486182246005193998007e5,
        -3.301918855321234414232308938454112213751e6,
        -5.856830900232338906742924836032279404702e6,
        -5.540672519616151584486240871424021377540e6,
        -3.006530901041386626148342989181721176919e6,
        -9.350378280513062139466966374330795935163e5,
        -1.566179100031063346901755685375732739511e5,
        -1.205016539620260779274902967231510804992e4,
        -2.724583156305709733221564484006088794284e2,
    };

    signgamp.* = 1;

    if (!std.math.isFinite(x))
        return x * x;

    if (x == 0) {
        if (std.math.signbit(x))
            signgamp.* = -1;
    }

    if (x < 0) {
        if (x < -2 and x > -50)
            return lgamma_neg128(x, signgamp);

        const q: f128 = -x;
        var p: f128 = float.floor(q);
        if (p == q)
            return (1 / float.abs(p - p));

        const halfp: f128 = p * 0.5;
        if (halfp == float.floor(halfp)) {
            signgamp.* = -1;
        } else {
            signgamp.* = 1;
        }

        if (q < 0x1p-120)
            return -float.log(q);

        var z: f128 = q - p;
        if (z > 0.5) {
            p += 1;
            z = p - q;
        }

        z = q * float.sin(std.math.pi * z);
        var i: i32 = undefined;
        const w: f128 = lgamma_r128(q, &i);
        z = float.log(std.math.pi / z) - w;
        return z;
    }

    if (x < 13.5) {
        var p: f128 = 0;
        const nx: f128 = float.floor(x + 0.5);
        const nn: i32 = cast(i32, nx, .{});
        switch (nn) {
            0 => {
                // log gamma (x + 1) = log(x) + log gamma(x)
                if (x < 0x1p-120) {
                    return -float.log(x);
                } else if (x <= 0.125) {
                    p = x * erf_data.neval(x, &RN1, 8) / erf_data.deval(x, &RD1, 7);
                } else if (x <= 0.375) {
                    const z: f128 = x - 0.25;
                    p = z * erf_data.neval(z, &RN1r25, 9) / erf_data.deval(z, &RD1r25, 8);
                    p += lgam1r25b;
                    p += lgam1r25a;
                } else if (x <= 0.625) {
                    var z: f128 = x + (1 - x0a);
                    z = z - x0b;
                    p = erf_data.neval(z, &RN1r5, 8) / erf_data.deval(z, &RD1r5, 8);
                    p = p * z * z;
                    p = p + y0b;
                    p = p + y0a;
                } else if (x <= 0.875) {
                    const z: f128 = x - 0.75;
                    p = z * erf_data.neval(z, &RN1r75, 8) / erf_data.deval(z, &RD1r75, 8);
                    p += lgam1r75b;
                    p += lgam1r75a;
                } else {
                    const z: f128 = x - 1;
                    p = z * erf_data.neval(z, &RN2, 9) / erf_data.deval(z, &RD2, 9);
                }
                p = p - float.log(x);
            },
            1 => {
                if (x < 0.875) {
                    if (x <= 0.625) {
                        var z: f128 = x + (1 - x0a);
                        z = z - x0b;
                        p = erf_data.neval(z, &RN1r5, 8) / erf_data.deval(z, &RD1r5, 8);
                        p = p * z * z;
                        p = p + y0b;
                        p = p + y0a;
                    } else if (x <= 0.875) {
                        const z: f128 = x - 0.75;
                        p = z * erf_data.neval(z, &RN1r75, 8) / erf_data.deval(z, &RD1r75, 8);
                        p += lgam1r75b;
                        p += lgam1r75a;
                    } else {
                        const z: f128 = x - 1;
                        p = z * erf_data.neval(z, &RN2, 9) / erf_data.deval(z, &RD2, 9);
                    }
                    p = p - float.log(x);
                } else if (x < 1) {
                    const z: f128 = x - 1;
                    p = z * erf_data.neval(z, &RNr9, 8) / erf_data.deval(z, &RDr9, 8);
                } else if (x == 1) {
                    p = 0;
                } else if (x <= 1.125) {
                    const z: f128 = x - 1;
                    p = z * erf_data.neval(z, &RN1, 8) / erf_data.deval(z, &RD1, 7);
                } else if (x <= 1.375) {
                    const z: f128 = x - 1.25;
                    p = z * erf_data.neval(z, &RN1r25, 9) / erf_data.deval(z, &RD1r25, 8);
                    p += lgam1r25b;
                    p += lgam1r25a;
                } else {
                    // 1.375 <= x+x0 <= 1.625
                    var z: f128 = x - x0a;
                    z = z - x0b;
                    p = erf_data.neval(z, &RN1r5, 8) / erf_data.deval(z, &RD1r5, 8);
                    p = p * z * z;
                    p = p + y0b;
                    p = p + y0a;
                }
            },
            2 => {
                if (x < 1.625) {
                    var z: f128 = x - x0a;
                    z = z - x0b;
                    p = erf_data.neval(z, &RN1r5, 8) / erf_data.deval(z, &RD1r5, 8);
                    p = p * z * z;
                    p = p + y0b;
                    p = p + y0a;
                } else if (x < 1.875) {
                    const z: f128 = x - 1.75;
                    p = z * erf_data.neval(z, &RN1r75, 8) / erf_data.deval(z, &RD1r75, 8);
                    p += lgam1r75b;
                    p += lgam1r75a;
                } else if (x == 2) {
                    p = 0;
                } else if (x < 2.375) {
                    const z: f128 = x - 2;
                    p = z * erf_data.neval(z, &RN2, 9) / erf_data.deval(z, &RD2, 9);
                } else {
                    const z: f128 = x - 2.5;
                    p = z * erf_data.neval(z, &RN2r5, 8) / erf_data.deval(z, &RD2r5, 8);
                    p += lgam2r5b;
                    p += lgam2r5a;
                }
            },
            3 => {
                if (x < 2.75) {
                    const z: f128 = x - 2.5;
                    p = z * erf_data.neval(z, &RN2r5, 8) / erf_data.deval(z, &RD2r5, 8);
                    p += lgam2r5b;
                    p += lgam2r5a;
                } else {
                    const z: f128 = x - 3;
                    p = z * erf_data.neval(z, &RN3, 9) / erf_data.deval(z, &RD3, 9);
                    p += lgam3b;
                    p += lgam3a;
                }
            },
            4 => {
                const z: f128 = x - 4;
                p = z * erf_data.neval(z, &RN4, 9) / erf_data.deval(z, &RD4, 9);
                p += lgam4b;
                p += lgam4a;
            },
            5 => {
                const z: f128 = x - 5;
                p = z * erf_data.neval(z, &RN5, 9) / erf_data.deval(z, &RD5, 8);
                p += lgam5b;
                p += lgam5a;
            },
            6 => {
                const z: f128 = x - 6;
                p = z * erf_data.neval(z, &RN6, 8) / erf_data.deval(z, &RD6, 8);
                p += lgam6b;
                p += lgam6a;
            },
            7 => {
                const z: f128 = x - 7;
                p = z * erf_data.neval(z, &RN7, 8) / erf_data.deval(z, &RD7, 7);
                p += lgam7b;
                p += lgam7a;
            },
            8 => {
                const z: f128 = x - 8;
                p = z * erf_data.neval(z, &RN8, 8) / erf_data.deval(z, &RD8, 7);
                p += lgam8b;
                p += lgam8a;
            },
            9 => {
                const z: f128 = x - 9;
                p = z * erf_data.neval(z, &RN9, 7) / erf_data.deval(z, &RD9, 7);
                p += lgam9b;
                p += lgam9a;
            },
            10 => {
                const z: f128 = x - 10;
                p = z * erf_data.neval(z, &RN10, 7) / erf_data.deval(z, &RD10, 7);
                p += lgam10b;
                p += lgam10a;
            },
            11 => {
                const z: f128 = x - 11;
                p = z * erf_data.neval(z, &RN11, 7) / erf_data.deval(z, &RD11, 6);
                p += lgam11b;
                p += lgam11a;
            },
            12 => {
                const z: f128 = x - 12;
                p = z * erf_data.neval(z, &RN12, 7) / erf_data.deval(z, &RD12, 6);
                p += lgam12b;
                p += lgam12a;
            },
            13 => {
                const z: f128 = x - 13;
                p = z * erf_data.neval(z, &RN13, 7) / erf_data.deval(z, &RD13, 6);
                p += lgam13b;
                p += lgam13a;
            },
            else => unreachable,
        }
        return p;
    }

    if (x > MAXLGM)
        return (cast(f128, signgamp.*, .{}) * huge * huge);

    if (x > 0x1p120)
        return x * (float.log(x) - 1);

    var q: f128 = ls2pi - x;
    q = (x - 0.5) * float.log(x) + q;
    if (x > 1.0e18)
        return q;

    const p: f128 = 1 / (x * x);
    q += erf_data.neval(p, &RASY, 12) / x;
    return q;
}
