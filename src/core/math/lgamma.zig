const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
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

    const fx: f32 = math.floor(x);
    const ax: f32 = math.abs(x);
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
                const ni: i32 = cast(i32, math.floor(-2 * x), .{});
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
        return math.sin(std.math.pi * x);
    } else {
        return math.cos(std.math.pi * (0.5 - x));
    }
}

// Compute cos (pi * X) for -0.25 <= X <= 0.5.
fn lg_cospi64(x: f64) f64 {
    if (x <= 0.25) {
        return math.cos(std.math.pi * x);
    } else {
        return math.sin(std.math.pi * (0.5 - x));
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
    var i: i32 = cast(i32, math.floor(-2 * x), .{});
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
        var j: i32 = cast(i32, math.floor(-8 * x) - 16, .{});
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

        return math.log1p(g * xdiff / (x - xn));
    }

    // The result we want is log (sinpi (X0) / sinpi (X))
    //  + log (gamma (1 - X0) / gamma (1 - X)).
    const x_idiff: f64 = math.abs(xn - x);
    const x0_idiff: f64 = math.abs(xn - x0_hi - x0_lo);
    var log_sinpi_ratio: f64 = undefined;
    if (x0_idiff < x_idiff * 0.5) {
        // Use log not log1p to avoid inaccuracy from log1p of arguments
        // close to -1.
        log_sinpi_ratio = math.log(lg_sinpi64(x0_idiff) / lg_sinpi64(x_idiff));
    } else {
        // Use log1p not log to avoid inaccuracy from log of arguments
        // close to 1.  X0DIFF2 has positive sign if X0 is further from
        // XN than X is from XN, negative sign otherwise.
        const x0diff2: f64 = (if ((i & 1) == 0) xdiff else -xdiff) * 0.5;
        const sx0d2: f64 = lg_sinpi64(x0diff2);
        const cx0d2: f64 = lg_cospi64(x0diff2);
        log_sinpi_ratio = math.log1p(2 * sx0d2 * (-sx0d2 + cx0d2 * lg_cotpi64(x_idiff)));
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
        log_gamma_adj = -math.log1p(prodm1);
    }

    const log_gamma_high: f64 = (xdiff * math.log1p((y0 - e_hi - e_lo + y0_eps) / e_hi) + (y - 0.5 + y_eps) * math.log1p(xdiff / y) + log_gamma_adj);
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
    if (ix < 0x3fd00000) return math.sin(pi * x);
    var y: f64 = -x; // x is assume negative

    // argument reduction, make sure inexact flag not raised if input
    // is an integer
    var z: f64 = math.floor(y);
    var n: i32 = undefined;
    if (z != y) { // inexact anyway
        y *= 0.5;
        y = 2.0 * (y - math.floor(y)); // y = |x| mod 2.0
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
        0 => y = math.sin(pi * y),
        1, 2 => y = math.cos(pi * (0.5 - y)),
        3, 4 => y = math.sin(pi * (1 - y)),
        5, 6 => y = -math.cos(pi * (y - 1.5)),
        else => y = math.sin(pi * (y - 2.0)),
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

        return 1 / math.abs(x);
    }

    if (ix < 0x3b900000) {
        @branchHint(.unlikely);
        // |x|<2**-70, return -log(|x|)
        if (hx < 0) {
            signgamp.* = -1;

            return -math.log(-x);
        } else {
            return -math.log(x);
        }
    }

    var nadj: f64 = undefined;
    var xx: f64 = x;
    if (hx < 0) {
        if (ix >= 0x43300000) {
            @branchHint(.unlikely);
            // |x|>=2**52, must be -integer
            return math.abs(x) / @as(f64, 0);
        }

        if (x < -2.0 and x > -28.0)
            return lgamma_neg64(x, signgamp);

        const t: f64 = sin_pi64(x);

        if (t == 0) return 1 / math.abs(t); // -integer

        nadj = math.log(pi / math.abs(t * x));
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
            r = -math.log(xx);
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
            else => r += math.log(z),
        }
        // 8.0 <= x < 2**58
    } else if (ix < 0x43900000) {
        const t: f64 = math.log(xx);
        const z: f64 = 1 / xx;
        const y: f64 = z * z;
        const w: f64 = w0 + z * (w1 + y * (w2 + y * (w3 + y * (w4 + y * (w5 + y * w6)))));
        r = (xx - 0.5) * (t - 1) + w;
    } else {
        // 2**58 <= x <= inf
        r = xx * (math.log(xx) - 1);
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
        return math.sin(std.math.pi * x);
    } else {
        return math.cos(std.math.pi * (0.5 - x));
    }
}

// Compute cos (pi * X) for -0.25 <= X <= 0.5.
fn lg_cospi128(x: f128) f128 {
    if (x <= 0.25) {
        return math.cos(std.math.pi * x);
    } else {
        return math.sin(std.math.pi * (0.5 - x));
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
    var i: i32 = cast(i32, math.floor(-2 * x), .{});
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
        var j: i32 = cast(i32, math.floor(-8 * x) - 16, .{});
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

        return math.log1p(g * xdiff / (x - xn));
    }

    // The result we want is log (sinpi (X0) / sinpi (X))
    //  + log (gamma (1 - X0) / gamma (1 - X)).
    const x_idiff: f128 = math.abs(xn - x);
    const x0_idiff: f128 = math.abs(xn - x0_hi - x0_lo);
    var log_sinpi_ratio: f128 = undefined;
    if (x0_idiff < x_idiff * 0.5) {
        // Use log not log1p to avoid inaccuracy from log1p of arguments
        // close to -1.
        log_sinpi_ratio = math.log(lg_sinpi128(x0_idiff) / lg_sinpi128(x_idiff));
    } else {
        // Use log1p not log to avoid inaccuracy from log of arguments
        // close to 1.  X0DIFF2 has positive sign if X0 is further from
        // XN than X is from XN, negative sign otherwise.
        const x0diff2 = (if ((i & 1) == 0) xdiff else -xdiff) * 0.5;
        const sx0d2 = lg_sinpi128(x0diff2);
        const cx0d2 = lg_cospi128(x0diff2);
        log_sinpi_ratio = math.log1p(2 * sx0d2 * (-sx0d2 + cx0d2 * lg_cotpi128(x_idiff)));
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
        log_gamma_adj = -math.log1p(prodm1);
    }

    const log_gamma_high: f128 = (xdiff * math.log1p((y0 - e_hi - e_lo + y0_eps) / e_hi) + (y - 0.5 + y_eps) * math.log1p(xdiff / y) + log_gamma_adj);
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
        var p: f128 = math.floor(q);
        if (p == q)
            return (1 / math.abs(p - p));

        const halfp: f128 = p * 0.5;
        if (halfp == math.floor(halfp)) {
            signgamp.* = -1;
        } else {
            signgamp.* = 1;
        }

        if (q < 0x1p-120)
            return -math.log(q);

        var z: f128 = q - p;
        if (z > 0.5) {
            p += 1;
            z = p - q;
        }

        z = q * math.sin(std.math.pi * z);
        var i: i32 = undefined;
        const w: f128 = lgamma_r128(q, &i);
        z = math.log(std.math.pi / z) - w;
        return z;
    }

    if (x < 13.5) {
        var p: f128 = 0;
        const nx: f128 = math.floor(x + 0.5);
        const nn: i32 = cast(i32, nx, .{});
        switch (nn) {
            0 => {
                // log gamma (x + 1) = log(x) + log gamma(x)
                if (x < 0x1p-120) {
                    return -math.log(x);
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
                p = p - math.log(x);
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
                    p = p - math.log(x);
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
        return x * (math.log(x) - 1);

    var q: f128 = ls2pi - x;
    q = (x - 0.5) * math.log(x) + q;
    if (x > 1.0e18)
        return q;

    const p: f128 = 1 / (x * x);
    q += erf_data.neval(p, &RASY, 12) / x;
    return q;
}

test lgamma {
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x0p+0, lgamma(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0xb.17218p-4, lgamma(@as(f32, 0x3p+0)));
    try std.testing.expectEqual(0x9.28682p-4, lgamma(@as(f32, 0x8p-4)));
    try std.testing.expectEqual(0x4.2c831p-4, lgamma(@as(f32, 0xb.33334p-4)));
    try std.testing.expectEqual(0x4.2c8328p-4, lgamma(@as(f32, 0xb.33333p-4)));
    try std.testing.expectEqual(-0x1.5db13cp-4, lgamma(@as(f32, 0x1.333334p+0)));
    try std.testing.expectEqual(-0x1.5db134p-4, lgamma(@as(f32, 0x1.333332p+0)));
    try std.testing.expectEqual(0x8.8bdd4p+60, lgamma(@as(f32, 0x3.8p+56)));
    try std.testing.expectEqual(0x3.72d03p+0, lgamma(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x3.7c0e1p+0, lgamma(@as(f32, -0x8p-8)));
    try std.testing.expectEqual(0x6.ee5008p+0, lgamma(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x6.ee99fp+0, lgamma(@as(f32, -0x4p-12)));
    try std.testing.expectEqual(0xa.65ae4p+0, lgamma(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0xa.65b09p+0, lgamma(@as(f32, -0x2p-16)));
    try std.testing.expectEqual(0xd.dce9dp+0, lgamma(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0xd.dce9fp+0, lgamma(@as(f32, -0x1p-20)));
    try std.testing.expectEqual(0x1.154246p+4, lgamma(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x1.154246p+4, lgamma(@as(f32, -0x8p-28)));
    try std.testing.expectEqual(0x1.4cb5ecp+4, lgamma(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x1.4cb5ecp+4, lgamma(@as(f32, -0x4p-32)));
    try std.testing.expectEqual(0x1.bb9d3cp+4, lgamma(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x1.bb9d3cp+4, lgamma(@as(f32, -0x1p-40)));
    try std.testing.expectEqual(0x2.2a848cp+4, lgamma(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x2.2a848cp+4, lgamma(@as(f32, -0x4p-52)));
    try std.testing.expectEqual(0x2.996bd8p+4, lgamma(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x2.996bd8p+4, lgamma(@as(f32, -0x1p-60)));
    try std.testing.expectEqual(0x2.c5c86p+4, lgamma(@as(f32, 0x1p-64)));
    try std.testing.expectEqual(0x2.c5c86p+4, lgamma(@as(f32, -0x1p-64)));
    try std.testing.expectEqual(0x3.085328p+4, lgamma(@as(f32, 0x4p-72)));
    try std.testing.expectEqual(0x3.085328p+4, lgamma(@as(f32, -0x4p-72)));
    try std.testing.expectEqual(0x4.550918p+4, lgamma(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x4.550918p+4, lgamma(@as(f32, -0x1p-100)));
    try std.testing.expectEqual(0x5.75628p+4, lgamma(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x5.75628p+4, lgamma(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x6.74768p+4, lgamma(@as(f32, -0x8p-152)));
    try std.testing.expectEqual(-0x7.d809fp-4, lgamma(@as(f32, -0x3.ec4298p+0)));
    try std.testing.expectEqual(0xf.ffff1p+124, lgamma(@as(f32, 0x3.12be0cp+120)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0x3.12be6p+120)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f32), lgamma(@as(f32, 0xf.fffffp+124)));
    try std.testing.expectEqual(-0x3.511bccp-20, lgamma(@as(f32, -0x3.f48e28p+0)));
    try std.testing.expectEqual(0x1.dd4b54p-20, lgamma(@as(f32, -0x3.f48e2cp+0)));
    try std.testing.expectEqual(-0x3.4a0c54p-24, lgamma(@as(f32, -0x3.24c1b8p+0)));
    try std.testing.expectEqual(-0x7.78a01p+24, lgamma(@as(f32, -0x7.fffff8p+20)));
    try std.testing.expectEqual(0x1.96ee68p+0, lgamma(@as(f32, -0x4p-4)));
    try std.testing.expectEqual(0x1.43f89ap+0, lgamma(@as(f32, -0x8p-4)));
    try std.testing.expectEqual(0x1.93616p+0, lgamma(@as(f32, -0xcp-4)));
    try std.testing.expectEqual(0x1.5dce78p+0, lgamma(@as(f32, -0x1.4p+0)));
    try std.testing.expectEqual(0xd.c2c0bp-4, lgamma(@as(f32, -0x1.8p+0)));
    try std.testing.expectEqual(0x1.041e66p+0, lgamma(@as(f32, -0x1.cp+0)));
    try std.testing.expectEqual(0x2.bec33cp+0, lgamma(@as(f32, -0x2.08p+0)));
    try std.testing.expectEqual(0x2.07061p+0, lgamma(@as(f32, -0x2.1p+0)));
    try std.testing.expectEqual(0x1.99a9fep+0, lgamma(@as(f32, -0x2.18p+0)));
    try std.testing.expectEqual(0x1.4b32e6p+0, lgamma(@as(f32, -0x2.2p+0)));
    try std.testing.expectEqual(0x1.0e0298p+0, lgamma(@as(f32, -0x2.28p+0)));
    try std.testing.expectEqual(0xd.c0af4p-4, lgamma(@as(f32, -0x2.3p+0)));
    try std.testing.expectEqual(0xb.21412p-4, lgamma(@as(f32, -0x2.38p+0)));
    try std.testing.expectEqual(0x8.e3559p-4, lgamma(@as(f32, -0x2.4p+0)));
    try std.testing.expectEqual(0x6.f371cp-4, lgamma(@as(f32, -0x2.48p+0)));
    try std.testing.expectEqual(0x5.448598p-4, lgamma(@as(f32, -0x2.5p+0)));
    try std.testing.expectEqual(0x3.cd82f8p-4, lgamma(@as(f32, -0x2.58p+0)));
    try std.testing.expectEqual(0x2.8804acp-4, lgamma(@as(f32, -0x2.6p+0)));
    try std.testing.expectEqual(0x1.6f830ep-4, lgamma(@as(f32, -0x2.68p+0)));
    try std.testing.expectEqual(0x8.0d79bp-8, lgamma(@as(f32, -0x2.7p+0)));
    try std.testing.expectEqual(-0x4.60fecp-8, lgamma(@as(f32, -0x2.78p+0)));
    try std.testing.expectEqual(-0xe.65fdp-8, lgamma(@as(f32, -0x2.8p+0)));
    try std.testing.expectEqual(-0x1.60773ep-4, lgamma(@as(f32, -0x2.88p+0)));
    try std.testing.expectEqual(-0x1.b3f01cp-4, lgamma(@as(f32, -0x2.9p+0)));
    try std.testing.expectEqual(-0x1.df9732p-4, lgamma(@as(f32, -0x2.98p+0)));
    try std.testing.expectEqual(-0x1.e15352p-4, lgamma(@as(f32, -0x2.ap+0)));
    try std.testing.expectEqual(-0x1.b5f706p-4, lgamma(@as(f32, -0x2.a8p+0)));
    try std.testing.expectEqual(-0x1.58f3aap-4, lgamma(@as(f32, -0x2.bp+0)));
    try std.testing.expectEqual(-0xc.3dd14p-8, lgamma(@as(f32, -0x2.b8p+0)));
    try std.testing.expectEqual(0x1.261e6ep-8, lgamma(@as(f32, -0x2.cp+0)));
    try std.testing.expectEqual(0x1.36e062p-4, lgamma(@as(f32, -0x2.c8p+0)));
    try std.testing.expectEqual(0x2.bd204p-4, lgamma(@as(f32, -0x2.dp+0)));
    try std.testing.expectEqual(0x4.c3b23p-4, lgamma(@as(f32, -0x2.d8p+0)));
    try std.testing.expectEqual(0x7.7e1cp-4, lgamma(@as(f32, -0x2.ep+0)));
    try std.testing.expectEqual(0xb.4d46bp-4, lgamma(@as(f32, -0x2.e8p+0)));
    try std.testing.expectEqual(0x1.10b1c8p+0, lgamma(@as(f32, -0x2.fp+0)));
    try std.testing.expectEqual(0x1.b6f672p+0, lgamma(@as(f32, -0x2.f8p+0)));
    try std.testing.expectEqual(0x1.a2dd72p+0, lgamma(@as(f32, -0x3.08p+0)));
    try std.testing.expectEqual(0xe.88019p-4, lgamma(@as(f32, -0x3.1p+0)));
    try std.testing.expectEqual(0x7.88aafp-4, lgamma(@as(f32, -0x3.18p+0)));
    try std.testing.expectEqual(0x2.780efp-4, lgamma(@as(f32, -0x3.2p+0)));
    try std.testing.expectEqual(-0x1.83b7aep-4, lgamma(@as(f32, -0x3.28p+0)));
    try std.testing.expectEqual(-0x4.cb8ccp-4, lgamma(@as(f32, -0x3.3p+0)));
    try std.testing.expectEqual(-0x7.92f0fp-4, lgamma(@as(f32, -0x3.38p+0)));
    try std.testing.expectEqual(-0x9.f86fcp-4, lgamma(@as(f32, -0x3.4p+0)));
    try std.testing.expectEqual(-0xc.0f85ep-4, lgamma(@as(f32, -0x3.48p+0)));
    try std.testing.expectEqual(-0xd.e5453p-4, lgamma(@as(f32, -0x3.5p+0)));
    try std.testing.expectEqual(-0xf.82bdbp-4, lgamma(@as(f32, -0x3.58p+0)));
    try std.testing.expectEqual(-0x1.0ee564p+0, lgamma(@as(f32, -0x3.6p+0)));
    try std.testing.expectEqual(-0x1.22c984p+0, lgamma(@as(f32, -0x3.68p+0)));
    try std.testing.expectEqual(-0x1.340abcp+0, lgamma(@as(f32, -0x3.7p+0)));
    try std.testing.expectEqual(-0x1.42ca4cp+0, lgamma(@as(f32, -0x3.78p+0)));
    try std.testing.expectEqual(-0x1.4f1b1p+0, lgamma(@as(f32, -0x3.8p+0)));
    try std.testing.expectEqual(-0x1.590312p+0, lgamma(@as(f32, -0x3.88p+0)));
    try std.testing.expectEqual(-0x1.607c0ap+0, lgamma(@as(f32, -0x3.9p+0)));
    try std.testing.expectEqual(-0x1.6572dap+0, lgamma(@as(f32, -0x3.98p+0)));
    try std.testing.expectEqual(-0x1.67c606p+0, lgamma(@as(f32, -0x3.ap+0)));
    try std.testing.expectEqual(-0x1.6742cep+0, lgamma(@as(f32, -0x3.a8p+0)));
    try std.testing.expectEqual(-0x1.63a05ap+0, lgamma(@as(f32, -0x3.bp+0)));
    try std.testing.expectEqual(-0x1.5c77fcp+0, lgamma(@as(f32, -0x3.b8p+0)));
    try std.testing.expectEqual(-0x1.513878p+0, lgamma(@as(f32, -0x3.cp+0)));
    try std.testing.expectEqual(-0x1.41107p+0, lgamma(@as(f32, -0x3.c8p+0)));
    try std.testing.expectEqual(-0x1.2ac7d6p+0, lgamma(@as(f32, -0x3.dp+0)));
    try std.testing.expectEqual(-0x1.0c75b6p+0, lgamma(@as(f32, -0x3.d8p+0)));
    try std.testing.expectEqual(-0xe.2e1c1p-4, lgamma(@as(f32, -0x3.ep+0)));
    try std.testing.expectEqual(-0xa.7fd7cp-4, lgamma(@as(f32, -0x3.e8p+0)));
    try std.testing.expectEqual(-0x4.e2a518p-4, lgamma(@as(f32, -0x3.fp+0)));
    try std.testing.expectEqual(0x5.614458p-4, lgamma(@as(f32, -0x3.f8p+0)));
    try std.testing.expectEqual(-0x2.11f044p+0, lgamma(@as(f32, -0x4.4p+0)));
    try std.testing.expectEqual(-0x2.d02648p+0, lgamma(@as(f32, -0x4.8p+0)));
    try std.testing.expectEqual(-0x2.e01b08p+0, lgamma(@as(f32, -0x4.cp+0)));
    try std.testing.expectEqual(-0x3.ba71e8p+0, lgamma(@as(f32, -0x5.4p+0)));
    try std.testing.expectEqual(-0x4.8490a8p+0, lgamma(@as(f32, -0x5.8p+0)));
    try std.testing.expectEqual(-0x4.9fe698p+0, lgamma(@as(f32, -0x5.cp+0)));
    try std.testing.expectEqual(-0x5.8f95f8p+0, lgamma(@as(f32, -0x6.4p+0)));
    try std.testing.expectEqual(-0x6.63bf1p+0, lgamma(@as(f32, -0x6.8p+0)));
    try std.testing.expectEqual(-0x6.88be6p+0, lgamma(@as(f32, -0x6.cp+0)));
    try std.testing.expectEqual(-0x7.8ab8ep+0, lgamma(@as(f32, -0x7.4p+0)));
    try std.testing.expectEqual(-0x8.678fcp+0, lgamma(@as(f32, -0x7.8p+0)));
    try std.testing.expectEqual(-0x8.94f4p+0, lgamma(@as(f32, -0x7.cp+0)));
    try std.testing.expectEqual(-0x9.a6efdp+0, lgamma(@as(f32, -0x8.4p+0)));
    try std.testing.expectEqual(-0xa.8b6b2p+0, lgamma(@as(f32, -0x8.8p+0)));
    try std.testing.expectEqual(-0xa.c03b1p+0, lgamma(@as(f32, -0x8.cp+0)));
    try std.testing.expectEqual(-0xb.e070cp+0, lgamma(@as(f32, -0x9.4p+0)));
    try std.testing.expectEqual(-0xc.cbbfdp+0, lgamma(@as(f32, -0x9.8p+0)));
    try std.testing.expectEqual(-0xd.07361p+0, lgamma(@as(f32, -0x9.cp+0)));
    try std.testing.expectEqual(-0xe.34393p+0, lgamma(@as(f32, -0xa.4p+0)));
    try std.testing.expectEqual(-0xf.25b38p+0, lgamma(@as(f32, -0xa.8p+0)));
    try std.testing.expectEqual(-0xf.672fep+0, lgamma(@as(f32, -0xa.cp+0)));
    try std.testing.expectEqual(-0x1.09fd68p+4, lgamma(@as(f32, -0xb.4p+0)));
    try std.testing.expectEqual(-0x1.196f12p+4, lgamma(@as(f32, -0xb.8p+0)));
    try std.testing.expectEqual(-0x1.1ddefp+4, lgamma(@as(f32, -0xb.cp+0)));
    try std.testing.expectEqual(-0x1.32140ap+4, lgamma(@as(f32, -0xc.4p+0)));
    try std.testing.expectEqual(-0x1.41d876p+4, lgamma(@as(f32, -0xc.8p+0)));
    try std.testing.expectEqual(-0x1.46996ep+4, lgamma(@as(f32, -0xc.cp+0)));
    try std.testing.expectEqual(-0x1.5b6c18p+4, lgamma(@as(f32, -0xd.4p+0)));
    try std.testing.expectEqual(-0x1.6b7d14p+4, lgamma(@as(f32, -0xd.8p+0)));
    try std.testing.expectEqual(-0x1.708936p+4, lgamma(@as(f32, -0xd.cp+0)));
    try std.testing.expectEqual(-0x1.85ee2ap+4, lgamma(@as(f32, -0xe.4p+0)));
    try std.testing.expectEqual(-0x1.964664p+4, lgamma(@as(f32, -0xe.8p+0)));
    try std.testing.expectEqual(-0x1.9b988ap+4, lgamma(@as(f32, -0xe.cp+0)));
    try std.testing.expectEqual(-0x1.b1860cp+4, lgamma(@as(f32, -0xf.4p+0)));
    try std.testing.expectEqual(-0x1.c220dep+4, lgamma(@as(f32, -0xf.8p+0)));
    try std.testing.expectEqual(-0x1.c7b48ep+4, lgamma(@as(f32, -0xf.cp+0)));
    try std.testing.expectEqual(-0x1.de2212p+4, lgamma(@as(f32, -0x1.04p+4)));
    try std.testing.expectEqual(-0x1.eefb6ep+4, lgamma(@as(f32, -0x1.08p+4)));
    try std.testing.expectEqual(-0x1.f4ccb8p+4, lgamma(@as(f32, -0x1.0cp+4)));
    try std.testing.expectEqual(-0x2.0bb2b8p+4, lgamma(@as(f32, -0x1.14p+4)));
    try std.testing.expectEqual(-0x2.1cc7p+4, lgamma(@as(f32, -0x1.18p+4)));
    try std.testing.expectEqual(-0x2.22d264p+4, lgamma(@as(f32, -0x1.1cp+4)));
    try std.testing.expectEqual(-0x2.3a2a2cp+4, lgamma(@as(f32, -0x1.24p+4)));
    try std.testing.expectEqual(-0x2.4b7634p+4, lgamma(@as(f32, -0x1.28p+4)));
    try std.testing.expectEqual(-0x2.51b89p+4, lgamma(@as(f32, -0x1.2cp+4)));
    try std.testing.expectEqual(-0x2.697c24p+4, lgamma(@as(f32, -0x1.34p+4)));
    try std.testing.expectEqual(-0x2.7afd04p+4, lgamma(@as(f32, -0x1.38p+4)));
    try std.testing.expectEqual(-0x2.81739p+4, lgamma(@as(f32, -0x1.3cp+4)));
    try std.testing.expectEqual(-0x2.999d8cp+4, lgamma(@as(f32, -0x1.44p+4)));
    try std.testing.expectEqual(-0x2.ab50acp+4, lgamma(@as(f32, -0x1.48p+4)));
    try std.testing.expectEqual(-0x2.b1f8dcp+4, lgamma(@as(f32, -0x1.4cp+4)));
    try std.testing.expectEqual(-0x2.ca846p+4, lgamma(@as(f32, -0x1.54p+4)));
    try std.testing.expectEqual(-0x2.dc676cp+4, lgamma(@as(f32, -0x1.58p+4)));
    try std.testing.expectEqual(-0x2.e33ef8p+4, lgamma(@as(f32, -0x1.5cp+4)));
    try std.testing.expectEqual(-0x2.fc2794p+4, lgamma(@as(f32, -0x1.64p+4)));
    try std.testing.expectEqual(-0x3.0e386p+4, lgamma(@as(f32, -0x1.68p+4)));
    try std.testing.expectEqual(-0x3.153d3p+4, lgamma(@as(f32, -0x1.6cp+4)));
    try std.testing.expectEqual(-0x3.2e7ed8p+4, lgamma(@as(f32, -0x1.74p+4)));
    try std.testing.expectEqual(-0x3.40bb74p+4, lgamma(@as(f32, -0x1.78p+4)));
    try std.testing.expectEqual(-0x3.47eb9cp+4, lgamma(@as(f32, -0x1.7cp+4)));
    try std.testing.expectEqual(-0x3.618298p+4, lgamma(@as(f32, -0x1.84p+4)));
    try std.testing.expectEqual(-0x3.73e938p+4, lgamma(@as(f32, -0x1.88p+4)));
    try std.testing.expectEqual(-0x3.7b42f4p+4, lgamma(@as(f32, -0x1.8cp+4)));
    try std.testing.expectEqual(-0x3.952bdcp+4, lgamma(@as(f32, -0x1.94p+4)));
    try std.testing.expectEqual(-0x3.a7bad8p+4, lgamma(@as(f32, -0x1.98p+4)));
    try std.testing.expectEqual(-0x3.af3c8cp+4, lgamma(@as(f32, -0x1.9cp+4)));
    try std.testing.expectEqual(-0x3.c97438p+4, lgamma(@as(f32, -0x1.a4p+4)));
    try std.testing.expectEqual(-0x3.dc2a08p+4, lgamma(@as(f32, -0x1.a8p+4)));
    try std.testing.expectEqual(-0x3.e3d23p+4, lgamma(@as(f32, -0x1.acp+4)));
    try std.testing.expectEqual(-0x3.fe55b8p+4, lgamma(@as(f32, -0x1.b4p+4)));
    try std.testing.expectEqual(-0x4.1130fp+4, lgamma(@as(f32, -0x1.b8p+4)));
    try std.testing.expectEqual(-0x4.18fe28p+4, lgamma(@as(f32, -0x1.bcp+4)));
    try std.testing.expectEqual(-0x4.33cad8p+4, lgamma(@as(f32, -0x1.c4p+4)));
    try std.testing.expectEqual(-0x4.46ca28p+4, lgamma(@as(f32, -0x1.c8p+4)));
    try std.testing.expectEqual(-0x4.4ebb2p+4, lgamma(@as(f32, -0x1.ccp+4)));
    try std.testing.expectEqual(-0x4.69ce7p+4, lgamma(@as(f32, -0x1.d4p+4)));
    try std.testing.expectEqual(-0x4.7cf098p+4, lgamma(@as(f32, -0x1.d8p+4)));
    try std.testing.expectEqual(-0x4.850428p+4, lgamma(@as(f32, -0x1.dcp+4)));
    try std.testing.expectEqual(-0x4.a05bcp+4, lgamma(@as(f32, -0x1.e4p+4)));
    try std.testing.expectEqual(-0x4.b39fap+4, lgamma(@as(f32, -0x1.e8p+4)));
    try std.testing.expectEqual(-0x4.bbd4ap+4, lgamma(@as(f32, -0x1.ecp+4)));
    try std.testing.expectEqual(-0x4.d76e4p+4, lgamma(@as(f32, -0x1.f4p+4)));
    try std.testing.expectEqual(-0x4.ead2cp+4, lgamma(@as(f32, -0x1.f8p+4)));
    try std.testing.expectEqual(-0x4.f32828p+4, lgamma(@as(f32, -0x1.fcp+4)));
    try std.testing.expectEqual(-0x5.0f01c8p+4, lgamma(@as(f32, -0x2.04p+4)));
    try std.testing.expectEqual(-0x5.2285e8p+4, lgamma(@as(f32, -0x2.08p+4)));
    try std.testing.expectEqual(-0x5.2afabp+4, lgamma(@as(f32, -0x2.0cp+4)));
    try std.testing.expectEqual(-0x5.47126p+4, lgamma(@as(f32, -0x2.14p+4)));
    try std.testing.expectEqual(-0x5.5ab538p+4, lgamma(@as(f32, -0x2.18p+4)));
    try std.testing.expectEqual(-0x5.63487p+4, lgamma(@as(f32, -0x2.1cp+4)));
    try std.testing.expectEqual(-0x5.7f9c6p+4, lgamma(@as(f32, -0x2.24p+4)));
    try std.testing.expectEqual(-0x5.935cf8p+4, lgamma(@as(f32, -0x2.28p+4)));
    try std.testing.expectEqual(-0x5.9c0dc8p+4, lgamma(@as(f32, -0x2.2cp+4)));
    try std.testing.expectEqual(-0x5.b89c38p+4, lgamma(@as(f32, -0x2.34p+4)));
    try std.testing.expectEqual(-0x5.cc79c8p+4, lgamma(@as(f32, -0x2.38p+4)));
    try std.testing.expectEqual(-0x5.d5475p+4, lgamma(@as(f32, -0x2.3cp+4)));
    try std.testing.expectEqual(-0x5.f20ea8p+4, lgamma(@as(f32, -0x2.44p+4)));
    try std.testing.expectEqual(-0x6.06086p+4, lgamma(@as(f32, -0x2.48p+4)));
    try std.testing.expectEqual(-0x6.0ef1ep+4, lgamma(@as(f32, -0x2.4cp+4)));
    try std.testing.expectEqual(-0x6.2bf09p+4, lgamma(@as(f32, -0x2.54p+4)));
    try std.testing.expectEqual(-0x6.4005bp+4, lgamma(@as(f32, -0x2.58p+4)));
    try std.testing.expectEqual(-0x6.490a68p+4, lgamma(@as(f32, -0x2.5cp+4)));
    try std.testing.expectEqual(-0x6.663ef8p+4, lgamma(@as(f32, -0x2.64p+4)));
    try std.testing.expectEqual(-0x6.7a6ec8p+4, lgamma(@as(f32, -0x2.68p+4)));
    try std.testing.expectEqual(-0x6.838ep+4, lgamma(@as(f32, -0x2.6cp+4)));
    try std.testing.expectEqual(-0x6.a0f718p+4, lgamma(@as(f32, -0x2.74p+4)));
    try std.testing.expectEqual(-0x6.b540e8p+4, lgamma(@as(f32, -0x2.78p+4)));
    try std.testing.expectEqual(-0x6.be79f8p+4, lgamma(@as(f32, -0x2.7cp+4)));
    try std.testing.expectEqual(-0x6.dc1648p+4, lgamma(@as(f32, -0x2.84p+4)));
    try std.testing.expectEqual(-0x6.f0797p+4, lgamma(@as(f32, -0x2.88p+4)));
    try std.testing.expectEqual(-0x6.f9cbb8p+4, lgamma(@as(f32, -0x2.8cp+4)));
    try std.testing.expectEqual(-0x7.1799f8p+4, lgamma(@as(f32, -0x2.94p+4)));
    try std.testing.expectEqual(-0x7.2c15ep+4, lgamma(@as(f32, -0x2.98p+4)));
    try std.testing.expectEqual(-0x7.3580cp+4, lgamma(@as(f32, -0x2.9cp+4)));
    try std.testing.expectEqual(-0x7.537fc8p+4, lgamma(@as(f32, -0x2.a4p+4)));
    try std.testing.expectEqual(-0x7.6813d8p+4, lgamma(@as(f32, -0x2.a8p+4)));
    try std.testing.expectEqual(-0x7.7196cp+4, lgamma(@as(f32, -0x2.acp+4)));
    try std.testing.expectEqual(-0x7.8fc56p+4, lgamma(@as(f32, -0x2.b4p+4)));
    try std.testing.expectEqual(-0x7.a4711p+4, lgamma(@as(f32, -0x2.b8p+4)));
    try std.testing.expectEqual(-0x7.ae0b7p+4, lgamma(@as(f32, -0x2.bcp+4)));
    try std.testing.expectEqual(-0x7.cc68ap+4, lgamma(@as(f32, -0x2.c4p+4)));
    try std.testing.expectEqual(-0x7.e12b68p+4, lgamma(@as(f32, -0x2.c8p+4)));
    try std.testing.expectEqual(-0x7.eadcb8p+4, lgamma(@as(f32, -0x2.ccp+4)));
    try std.testing.expectEqual(-0x8.09677p+4, lgamma(@as(f32, -0x2.d4p+4)));
    try std.testing.expectEqual(-0x8.1e40cp+4, lgamma(@as(f32, -0x2.d8p+4)));
    try std.testing.expectEqual(-0x8.28088p+4, lgamma(@as(f32, -0x2.dcp+4)));
    try std.testing.expectEqual(-0x8.46bfcp+4, lgamma(@as(f32, -0x2.e4p+4)));
    try std.testing.expectEqual(-0x8.5baf2p+4, lgamma(@as(f32, -0x2.e8p+4)));
    try std.testing.expectEqual(-0x8.658cep+4, lgamma(@as(f32, -0x2.ecp+4)));
    try std.testing.expectEqual(-0x8.846fbp+4, lgamma(@as(f32, -0x2.f4p+4)));
    try std.testing.expectEqual(-0x8.9974bp+4, lgamma(@as(f32, -0x2.f8p+4)));
    try std.testing.expectEqual(-0x8.a367fp+4, lgamma(@as(f32, -0x2.fcp+4)));
    try std.testing.expectEqual(-0x8.c2756p+4, lgamma(@as(f32, -0x3.04p+4)));
    try std.testing.expectEqual(-0x8.d78f9p+4, lgamma(@as(f32, -0x3.08p+4)));
    try std.testing.expectEqual(-0x8.e197ep+4, lgamma(@as(f32, -0x3.0cp+4)));
    try std.testing.expectEqual(-0x9.00cf2p+4, lgamma(@as(f32, -0x3.14p+4)));
    try std.testing.expectEqual(-0x9.15fe1p+4, lgamma(@as(f32, -0x3.18p+4)));
    try std.testing.expectEqual(-0x9.201bp+4, lgamma(@as(f32, -0x3.1cp+4)));
    try std.testing.expectEqual(-0x9.3f7b3p+4, lgamma(@as(f32, -0x3.24p+4)));
    try std.testing.expectEqual(-0x9.54be7p+4, lgamma(@as(f32, -0x3.28p+4)));
    try std.testing.expectEqual(-0x9.5eefap+4, lgamma(@as(f32, -0x3.2cp+4)));
    try std.testing.expectEqual(-0x9.7e78p+4, lgamma(@as(f32, -0x3.34p+4)));
    try std.testing.expectEqual(-0x9.93cf3p+4, lgamma(@as(f32, -0x3.38p+4)));
    try std.testing.expectEqual(-0x9.9e143p+4, lgamma(@as(f32, -0x3.3cp+4)));
    try std.testing.expectEqual(-0x9.bdc3fp+4, lgamma(@as(f32, -0x3.44p+4)));
    try std.testing.expectEqual(-0x9.d32ebp+4, lgamma(@as(f32, -0x3.48p+4)));
    try std.testing.expectEqual(-0x9.dd872p+4, lgamma(@as(f32, -0x3.4cp+4)));
    try std.testing.expectEqual(-0x9.fd5d8p+4, lgamma(@as(f32, -0x3.54p+4)));
    try std.testing.expectEqual(-0xa.12db7p+4, lgamma(@as(f32, -0x3.58p+4)));
    try std.testing.expectEqual(-0xa.1d47p+4, lgamma(@as(f32, -0x3.5cp+4)));
    try std.testing.expectEqual(-0xa.3d435p+4, lgamma(@as(f32, -0x3.64p+4)));
    try std.testing.expectEqual(-0xa.52d41p+4, lgamma(@as(f32, -0x3.68p+4)));
    try std.testing.expectEqual(-0xa.5d526p+4, lgamma(@as(f32, -0x3.6cp+4)));
    try std.testing.expectEqual(-0xa.7d73fp+4, lgamma(@as(f32, -0x3.74p+4)));
    try std.testing.expectEqual(-0xa.93173p+4, lgamma(@as(f32, -0x3.78p+4)));
    try std.testing.expectEqual(-0xa.9da7ep+4, lgamma(@as(f32, -0x3.7cp+4)));
    try std.testing.expectEqual(-0xa.bdeep+4, lgamma(@as(f32, -0x3.84p+4)));
    try std.testing.expectEqual(-0xa.d3a37p+4, lgamma(@as(f32, -0x3.88p+4)));
    try std.testing.expectEqual(-0xa.de463p+4, lgamma(@as(f32, -0x3.8cp+4)));
    try std.testing.expectEqual(-0xa.feb04p+4, lgamma(@as(f32, -0x3.94p+4)));
    try std.testing.expectEqual(-0xb.14779p+4, lgamma(@as(f32, -0x3.98p+4)));
    try std.testing.expectEqual(-0xb.1f2c1p+4, lgamma(@as(f32, -0x3.9cp+4)));
    try std.testing.expectEqual(-0xb.3fb98p+4, lgamma(@as(f32, -0x3.a4p+4)));
    try std.testing.expectEqual(-0xb.55924p+4, lgamma(@as(f32, -0x3.a8p+4)));
    try std.testing.expectEqual(-0xb.60585p+4, lgamma(@as(f32, -0x3.acp+4)));
    try std.testing.expectEqual(-0xb.81086p+4, lgamma(@as(f32, -0x3.b4p+4)));
    try std.testing.expectEqual(-0xb.96f27p+4, lgamma(@as(f32, -0x3.b8p+4)));
    try std.testing.expectEqual(-0xb.a1c9ap+4, lgamma(@as(f32, -0x3.bcp+4)));
    try std.testing.expectEqual(-0xb.c29bep+4, lgamma(@as(f32, -0x3.c4p+4)));
    try std.testing.expectEqual(-0xb.d896ep+4, lgamma(@as(f32, -0x3.c8p+4)));
    try std.testing.expectEqual(-0xb.e37efp+4, lgamma(@as(f32, -0x3.ccp+4)));
    try std.testing.expectEqual(0x1.0a2b24p+4, lgamma(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0xf.f1402p+0, lgamma(@as(f32, -0x1.000002p+0)));
    try std.testing.expectEqual(0xf.3fce1p+0, lgamma(@as(f32, -0x1.fffffep+0)));
    try std.testing.expectEqual(0xe.8e5bfp+0, lgamma(@as(f32, -0x2.000004p+0)));
    try std.testing.expectEqual(0xd.751d5p+0, lgamma(@as(f32, -0x2.fffffcp+0)));
    try std.testing.expectEqual(0xd.751d5p+0, lgamma(@as(f32, -0x3.000004p+0)));
    try std.testing.expectEqual(0xc.12392p+0, lgamma(@as(f32, -0x3.fffffcp+0)));
    try std.testing.expectEqual(0xb.60c7p+0, lgamma(@as(f32, -0x4.000008p+0)));
    try std.testing.expectEqual(0x9.c4c2fp+0, lgamma(@as(f32, -0x4.fffff8p+0)));
    try std.testing.expectEqual(0x9.c4c2ep+0, lgamma(@as(f32, -0x5.000008p+0)));
    try std.testing.expectEqual(0x7.fa1238p+0, lgamma(@as(f32, -0x5.fffff8p+0)));
    try std.testing.expectEqual(0x7.fa1218p+0, lgamma(@as(f32, -0x6.000008p+0)));
    try std.testing.expectEqual(0x6.07eb1p+0, lgamma(@as(f32, -0x6.fffff8p+0)));
    try std.testing.expectEqual(0x6.07eafp+0, lgamma(@as(f32, -0x7.000008p+0)));
    try std.testing.expectEqual(0x3.f394c8p+0, lgamma(@as(f32, -0x7.fffff8p+0)));
    try std.testing.expectEqual(0x3.42227cp+0, lgamma(@as(f32, -0x8.00001p+0)));
    try std.testing.expectEqual(0x1.0fa572p+0, lgamma(@as(f32, -0x8.fffffp+0)));
    try std.testing.expectEqual(0x1.0fa52ap+0, lgamma(@as(f32, -0x9.00001p+0)));
    try std.testing.expectEqual(-0x1.3dd0c4p+0, lgamma(@as(f32, -0x9.fffffp+0)));
    try std.testing.expectEqual(-0x1.3dd10ep+0, lgamma(@as(f32, -0xa.00001p+0)));
    try std.testing.expectEqual(-0x3.a3ad38p+0, lgamma(@as(f32, -0xa.fffffp+0)));
    try std.testing.expectEqual(-0x3.a3ad88p+0, lgamma(@as(f32, -0xb.00001p+0)));
    try std.testing.expectEqual(-0x6.1fd01p+0, lgamma(@as(f32, -0xb.fffffp+0)));
    try std.testing.expectEqual(-0x6.1fd06p+0, lgamma(@as(f32, -0xc.00001p+0)));
    try std.testing.expectEqual(-0x8.b0709p+0, lgamma(@as(f32, -0xc.fffffp+0)));
    try std.testing.expectEqual(-0x8.b070ep+0, lgamma(@as(f32, -0xd.00001p+0)));
    try std.testing.expectEqual(-0xb.5409dp+0, lgamma(@as(f32, -0xd.fffffp+0)));
    try std.testing.expectEqual(-0xb.540a3p+0, lgamma(@as(f32, -0xe.00001p+0)));
    try std.testing.expectEqual(-0xe.094cap+0, lgamma(@as(f32, -0xe.fffffp+0)));
    try std.testing.expectEqual(-0xe.094cfp+0, lgamma(@as(f32, -0xf.00001p+0)));
    try std.testing.expectEqual(-0x1.0cf15p+4, lgamma(@as(f32, -0xf.fffffp+0)));
    try std.testing.expectEqual(-0x1.18087ap+4, lgamma(@as(f32, -0x1.000002p+4)));
    try std.testing.expectEqual(-0x1.455d46p+4, lgamma(@as(f32, -0x1.0ffffep+4)));
    try std.testing.expectEqual(-0x1.455d52p+4, lgamma(@as(f32, -0x1.100002p+4)));
    try std.testing.expectEqual(-0x1.739c3cp+4, lgamma(@as(f32, -0x1.1ffffep+4)));
    try std.testing.expectEqual(-0x1.739c48p+4, lgamma(@as(f32, -0x1.200002p+4)));
    try std.testing.expectEqual(-0x1.a2b8a8p+4, lgamma(@as(f32, -0x1.2ffffep+4)));
    try std.testing.expectEqual(-0x1.a2b8b4p+4, lgamma(@as(f32, -0x1.300002p+4)));
    try std.testing.expectEqual(-0x1.d2a72cp+4, lgamma(@as(f32, -0x1.3ffffep+4)));
    try std.testing.expectEqual(-0x1.d2a738p+4, lgamma(@as(f32, -0x1.400002p+4)));
    try std.testing.expectEqual(-0x2.035d88p+4, lgamma(@as(f32, -0x1.4ffffep+4)));
    try std.testing.expectEqual(-0x2.035d98p+4, lgamma(@as(f32, -0x1.500002p+4)));
    try std.testing.expectEqual(-0x2.34d274p+4, lgamma(@as(f32, -0x1.5ffffep+4)));
    try std.testing.expectEqual(-0x2.34d28p+4, lgamma(@as(f32, -0x1.600002p+4)));
    try std.testing.expectEqual(-0x2.66fd7p+4, lgamma(@as(f32, -0x1.6ffffep+4)));
    try std.testing.expectEqual(-0x2.66fd7cp+4, lgamma(@as(f32, -0x1.700002p+4)));
    try std.testing.expectEqual(-0x2.99d6bcp+4, lgamma(@as(f32, -0x1.7ffffep+4)));
    try std.testing.expectEqual(-0x2.99d6ccp+4, lgamma(@as(f32, -0x1.800002p+4)));
    try std.testing.expectEqual(-0x2.cd574p+4, lgamma(@as(f32, -0x1.8ffffep+4)));
    try std.testing.expectEqual(-0x2.cd575p+4, lgamma(@as(f32, -0x1.900002p+4)));
    try std.testing.expectEqual(-0x3.01786cp+4, lgamma(@as(f32, -0x1.9ffffep+4)));
    try std.testing.expectEqual(-0x3.017878p+4, lgamma(@as(f32, -0x1.a00002p+4)));
    try std.testing.expectEqual(-0x3.36342cp+4, lgamma(@as(f32, -0x1.affffep+4)));
    try std.testing.expectEqual(-0x3.363438p+4, lgamma(@as(f32, -0x1.b00002p+4)));
    try std.testing.expectEqual(-0x3.6b84ep+4, lgamma(@as(f32, -0x1.bffffep+4)));
    try std.testing.expectEqual(-0x3.6b84ecp+4, lgamma(@as(f32, -0x1.c00002p+4)));
    try std.testing.expectEqual(-0x3.a1655p+4, lgamma(@as(f32, -0x1.cffffep+4)));
    try std.testing.expectEqual(-0x3.a1656p+4, lgamma(@as(f32, -0x1.d00002p+4)));
    try std.testing.expectEqual(-0x3.d7d0ap+4, lgamma(@as(f32, -0x1.dffffep+4)));
    try std.testing.expectEqual(-0x3.d7d0acp+4, lgamma(@as(f32, -0x1.e00002p+4)));
    try std.testing.expectEqual(-0x4.0ec24p+4, lgamma(@as(f32, -0x1.effffep+4)));
    try std.testing.expectEqual(-0x4.0ec248p+4, lgamma(@as(f32, -0x1.f00002p+4)));
    try std.testing.expectEqual(-0x4.4635ep+4, lgamma(@as(f32, -0x1.fffffep+4)));
    try std.testing.expectEqual(-0x4.514d18p+4, lgamma(@as(f32, -0x2.000004p+4)));
    try std.testing.expectEqual(-0x4.893ebp+4, lgamma(@as(f32, -0x2.0ffffcp+4)));
    try std.testing.expectEqual(-0x4.893ec8p+4, lgamma(@as(f32, -0x2.100004p+4)));
    try std.testing.expectEqual(-0x4.c1aaa8p+4, lgamma(@as(f32, -0x2.1ffffcp+4)));
    try std.testing.expectEqual(-0x4.c1aac8p+4, lgamma(@as(f32, -0x2.200004p+4)));
    try std.testing.expectEqual(-0x4.fa8d6p+4, lgamma(@as(f32, -0x2.2ffffcp+4)));
    try std.testing.expectEqual(-0x4.fa8d78p+4, lgamma(@as(f32, -0x2.300004p+4)));
    try std.testing.expectEqual(-0x5.33e378p+4, lgamma(@as(f32, -0x2.3ffffcp+4)));
    try std.testing.expectEqual(-0x5.33e39p+4, lgamma(@as(f32, -0x2.400004p+4)));
    try std.testing.expectEqual(-0x5.6da9c8p+4, lgamma(@as(f32, -0x2.4ffffcp+4)));
    try std.testing.expectEqual(-0x5.6da9ep+4, lgamma(@as(f32, -0x2.500004p+4)));
    try std.testing.expectEqual(-0x5.a7dd58p+4, lgamma(@as(f32, -0x2.5ffffcp+4)));
    try std.testing.expectEqual(-0x5.a7dd7p+4, lgamma(@as(f32, -0x2.600004p+4)));
    try std.testing.expectEqual(-0x5.e27b48p+4, lgamma(@as(f32, -0x2.6ffffcp+4)));
    try std.testing.expectEqual(-0x5.e27b68p+4, lgamma(@as(f32, -0x2.700004p+4)));
    try std.testing.expectEqual(-0x6.1d80fp+4, lgamma(@as(f32, -0x2.7ffffcp+4)));
    try std.testing.expectEqual(-0x6.1d8108p+4, lgamma(@as(f32, -0x2.800004p+4)));
    try std.testing.expectEqual(-0x6.58ebb8p+4, lgamma(@as(f32, -0x2.8ffffcp+4)));
    try std.testing.expectEqual(-0x6.58ebd8p+4, lgamma(@as(f32, -0x2.900004p+4)));
    try std.testing.expectEqual(-0x6.94b938p+4, lgamma(@as(f32, -0x2.9ffffcp+4)));
    try std.testing.expectEqual(-0x6.94b958p+4, lgamma(@as(f32, -0x2.a00004p+4)));
    try std.testing.expectEqual(-0x6.d0e718p+4, lgamma(@as(f32, -0x2.affffcp+4)));
    try std.testing.expectEqual(-0x6.d0e738p+4, lgamma(@as(f32, -0x2.b00004p+4)));
    try std.testing.expectEqual(-0x7.0d732p+4, lgamma(@as(f32, -0x2.bffffcp+4)));
    try std.testing.expectEqual(-0x7.0d734p+4, lgamma(@as(f32, -0x2.c00004p+4)));
    try std.testing.expectEqual(-0x7.4a5b38p+4, lgamma(@as(f32, -0x2.cffffcp+4)));
    try std.testing.expectEqual(-0x7.4a5b58p+4, lgamma(@as(f32, -0x2.d00004p+4)));
    try std.testing.expectEqual(-0x7.879d58p+4, lgamma(@as(f32, -0x2.dffffcp+4)));
    try std.testing.expectEqual(-0x7.879d7p+4, lgamma(@as(f32, -0x2.e00004p+4)));
    try std.testing.expectEqual(-0x7.c53788p+4, lgamma(@as(f32, -0x2.effffcp+4)));
    try std.testing.expectEqual(-0x7.c537a8p+4, lgamma(@as(f32, -0x2.f00004p+4)));
    try std.testing.expectEqual(-0x8.0328p+4, lgamma(@as(f32, -0x2.fffffcp+4)));
    try std.testing.expectEqual(-0x8.03282p+4, lgamma(@as(f32, -0x3.000004p+4)));
    try std.testing.expectEqual(-0x8.416cep+4, lgamma(@as(f32, -0x3.0ffffcp+4)));
    try std.testing.expectEqual(-0x8.416dp+4, lgamma(@as(f32, -0x3.100004p+4)));
    try std.testing.expectEqual(-0x8.80048p+4, lgamma(@as(f32, -0x3.1ffffcp+4)));
    try std.testing.expectEqual(-0x8.8004ap+4, lgamma(@as(f32, -0x3.200004p+4)));
    try std.testing.expectEqual(-0x8.beed4p+4, lgamma(@as(f32, -0x3.2ffffcp+4)));
    try std.testing.expectEqual(-0x8.beed6p+4, lgamma(@as(f32, -0x3.300004p+4)));
    try std.testing.expectEqual(-0x8.fe259p+4, lgamma(@as(f32, -0x3.3ffffcp+4)));
    try std.testing.expectEqual(-0x8.fe25bp+4, lgamma(@as(f32, -0x3.400004p+4)));
    try std.testing.expectEqual(-0x9.3dabep+4, lgamma(@as(f32, -0x3.4ffffcp+4)));
    try std.testing.expectEqual(-0x9.3dacp+4, lgamma(@as(f32, -0x3.500004p+4)));
    try std.testing.expectEqual(-0x9.7d7ecp+4, lgamma(@as(f32, -0x3.5ffffcp+4)));
    try std.testing.expectEqual(-0x9.7d7eep+4, lgamma(@as(f32, -0x3.600004p+4)));
    try std.testing.expectEqual(-0x9.bd9cdp+4, lgamma(@as(f32, -0x3.6ffffcp+4)));
    try std.testing.expectEqual(-0x9.bd9cfp+4, lgamma(@as(f32, -0x3.700004p+4)));
    try std.testing.expectEqual(-0x9.fe04ap+4, lgamma(@as(f32, -0x3.7ffffcp+4)));
    try std.testing.expectEqual(-0x9.fe04cp+4, lgamma(@as(f32, -0x3.800004p+4)));
    try std.testing.expectEqual(-0xa.3eb5p+4, lgamma(@as(f32, -0x3.8ffffcp+4)));
    try std.testing.expectEqual(-0xa.3eb52p+4, lgamma(@as(f32, -0x3.900004p+4)));
    try std.testing.expectEqual(-0xa.7fac9p+4, lgamma(@as(f32, -0x3.9ffffcp+4)));
    try std.testing.expectEqual(-0xa.7facbp+4, lgamma(@as(f32, -0x3.a00004p+4)));
    try std.testing.expectEqual(-0xa.c0ea2p+4, lgamma(@as(f32, -0x3.affffcp+4)));
    try std.testing.expectEqual(-0xa.c0ea4p+4, lgamma(@as(f32, -0x3.b00004p+4)));
    try std.testing.expectEqual(-0xb.026c9p+4, lgamma(@as(f32, -0x3.bffffcp+4)));
    try std.testing.expectEqual(-0xb.026cbp+4, lgamma(@as(f32, -0x3.c00004p+4)));
    try std.testing.expectEqual(0x4.2b2b5p-24, lgamma(@as(f32, -0x2.74ff9p+0)));
    try std.testing.expectEqual(-0x1.e4cf24p-24, lgamma(@as(f32, -0x2.74ff94p+0)));
    try std.testing.expectEqual(-0x2.6b417p-24, lgamma(@as(f32, -0x2.bf682p+0)));
    try std.testing.expectEqual(0x5.3d0a3p-24, lgamma(@as(f32, -0x2.bf6824p+0)));
    try std.testing.expectEqual(0x1.bd69b6p-20, lgamma(@as(f32, -0x3.24c1b4p+0)));
    try std.testing.expectEqual(-0x3.4a0c54p-24, lgamma(@as(f32, -0x3.24c1b8p+0)));
    try std.testing.expectEqual(-0x3.511bccp-20, lgamma(@as(f32, -0x3.f48e28p+0)));
    try std.testing.expectEqual(0x1.dd4b54p-20, lgamma(@as(f32, -0x3.f48e2cp+0)));
    try std.testing.expectEqual(0xa.3165dp-20, lgamma(@as(f32, -0x4.0a1398p+0)));
    try std.testing.expectEqual(-0x3.33cb58p-20, lgamma(@as(f32, -0x4.0a13ap+0)));
    try std.testing.expectEqual(-0x3.02165cp-16, lgamma(@as(f32, -0x4.fdd5d8p+0)));
    try std.testing.expectEqual(0xa.22e78p-20, lgamma(@as(f32, -0x4.fdd5ep+0)));
    try std.testing.expectEqual(0x2.e258fp-16, lgamma(@as(f32, -0x5.021a9p+0)));
    try std.testing.expectEqual(-0xf.89067p-20, lgamma(@as(f32, -0x5.021a98p+0)));
    try std.testing.expectEqual(-0xf.15ee1p-16, lgamma(@as(f32, -0x5.ffa4b8p+0)));
    try std.testing.expectEqual(0x7.4bb0fp-16, lgamma(@as(f32, -0x5.ffa4cp+0)));
    try std.testing.expectEqual(0x3.e9df58p-16, lgamma(@as(f32, -0x6.005ac8p+0)));
    try std.testing.expectEqual(-0x1.2b35eep-12, lgamma(@as(f32, -0x6.005adp+0)));
    try std.testing.expectEqual(-0x7.313b9p-12, lgamma(@as(f32, -0x6.fff2f8p+0)));
    try std.testing.expectEqual(0x2.a3598cp-12, lgamma(@as(f32, -0x6.fff3p+0)));
    try std.testing.expectEqual(0x9.39801p-12, lgamma(@as(f32, -0x7.000cf8p+0)));
    try std.testing.expectEqual(-0xa.32834p-16, lgamma(@as(f32, -0x7.000dp+0)));
    try std.testing.expectEqual(-0x4.cccb88p-8, lgamma(@as(f32, -0x7.fffe58p+0)));
    try std.testing.expectEqual(0x1.37b06p-12, lgamma(@as(f32, -0x7.fffe6p+0)));
    try std.testing.expectEqual(0xc.86027p-16, lgamma(@as(f32, -0x8.0001ap+0)));
    try std.testing.expectEqual(-0x9.9cf5ep-8, lgamma(@as(f32, -0x8.0001bp+0)));
    try std.testing.expectEqual(-0x9.98ed1p-8, lgamma(@as(f32, -0x8.ffffdp+0)));
    try std.testing.expectEqual(0x5.e337e8p-4, lgamma(@as(f32, -0x8.ffffep+0)));
    try std.testing.expectEqual(0x5.e32ee8p-4, lgamma(@as(f32, -0x9.00002p+0)));
    try std.testing.expectEqual(-0x9.99c53p-8, lgamma(@as(f32, -0x9.00003p+0)));
    try std.testing.expectEqual(-0x1.3dd0c4p+0, lgamma(@as(f32, -0x9.fffffp+0)));
    try std.testing.expectEqual(-0x1.3dd10ep+0, lgamma(@as(f32, -0xa.00001p+0)));
    try std.testing.expectEqual(-0x3.a3ad38p+0, lgamma(@as(f32, -0xa.fffffp+0)));
    try std.testing.expectEqual(-0x3.a3ad88p+0, lgamma(@as(f32, -0xb.00001p+0)));
    try std.testing.expectEqual(-0x6.1fd01p+0, lgamma(@as(f32, -0xb.fffffp+0)));
    try std.testing.expectEqual(-0x6.1fd06p+0, lgamma(@as(f32, -0xc.00001p+0)));
    try std.testing.expectEqual(-0x8.b0709p+0, lgamma(@as(f32, -0xc.fffffp+0)));
    try std.testing.expectEqual(-0x8.b070ep+0, lgamma(@as(f32, -0xd.00001p+0)));
    try std.testing.expectEqual(-0xb.5409dp+0, lgamma(@as(f32, -0xd.fffffp+0)));
    try std.testing.expectEqual(-0xb.540a3p+0, lgamma(@as(f32, -0xe.00001p+0)));
    try std.testing.expectEqual(-0xe.094cap+0, lgamma(@as(f32, -0xe.fffffp+0)));
    try std.testing.expectEqual(-0xe.094cfp+0, lgamma(@as(f32, -0xf.00001p+0)));
    try std.testing.expectEqual(-0x1.0cf15p+4, lgamma(@as(f32, -0xf.fffffp+0)));
    try std.testing.expectEqual(-0x1.18087ap+4, lgamma(@as(f32, -0x1.000002p+4)));
    try std.testing.expectEqual(-0x1.455d46p+4, lgamma(@as(f32, -0x1.0ffffep+4)));
    try std.testing.expectEqual(-0x1.455d52p+4, lgamma(@as(f32, -0x1.100002p+4)));
    try std.testing.expectEqual(-0x1.739c3cp+4, lgamma(@as(f32, -0x1.1ffffep+4)));
    try std.testing.expectEqual(-0x1.739c48p+4, lgamma(@as(f32, -0x1.200002p+4)));
    try std.testing.expectEqual(-0x1.a2b8a8p+4, lgamma(@as(f32, -0x1.2ffffep+4)));
    try std.testing.expectEqual(-0x1.a2b8b4p+4, lgamma(@as(f32, -0x1.300002p+4)));
    try std.testing.expectEqual(-0x1.d2a72cp+4, lgamma(@as(f32, -0x1.3ffffep+4)));
    try std.testing.expectEqual(-0x1.d2a738p+4, lgamma(@as(f32, -0x1.400002p+4)));
    try std.testing.expectEqual(-0x2.035d88p+4, lgamma(@as(f32, -0x1.4ffffep+4)));
    try std.testing.expectEqual(-0x2.035d98p+4, lgamma(@as(f32, -0x1.500002p+4)));
    try std.testing.expectEqual(-0x2.34d274p+4, lgamma(@as(f32, -0x1.5ffffep+4)));
    try std.testing.expectEqual(-0x2.34d28p+4, lgamma(@as(f32, -0x1.600002p+4)));
    try std.testing.expectEqual(-0x2.66fd7p+4, lgamma(@as(f32, -0x1.6ffffep+4)));
    try std.testing.expectEqual(-0x2.66fd7cp+4, lgamma(@as(f32, -0x1.700002p+4)));
    try std.testing.expectEqual(-0x2.99d6bcp+4, lgamma(@as(f32, -0x1.7ffffep+4)));
    try std.testing.expectEqual(-0x2.99d6ccp+4, lgamma(@as(f32, -0x1.800002p+4)));
    try std.testing.expectEqual(-0x2.cd574p+4, lgamma(@as(f32, -0x1.8ffffep+4)));
    try std.testing.expectEqual(-0x2.cd575p+4, lgamma(@as(f32, -0x1.900002p+4)));
    try std.testing.expectEqual(-0x3.01786cp+4, lgamma(@as(f32, -0x1.9ffffep+4)));
    try std.testing.expectEqual(-0x3.017878p+4, lgamma(@as(f32, -0x1.a00002p+4)));
    try std.testing.expectEqual(-0x3.36342cp+4, lgamma(@as(f32, -0x1.affffep+4)));
    try std.testing.expectEqual(-0x3.363438p+4, lgamma(@as(f32, -0x1.b00002p+4)));
    try std.testing.expectEqual(-0x3.6b84ep+4, lgamma(@as(f32, -0x1.bffffep+4)));
    try std.testing.expectEqual(-0x3.6b84ecp+4, lgamma(@as(f32, -0x1.c00002p+4)));
    try std.testing.expectEqual(-0x3.a1655p+4, lgamma(@as(f32, -0x1.cffffep+4)));
    try std.testing.expectEqual(-0x3.a1656p+4, lgamma(@as(f32, -0x1.d00002p+4)));
    try std.testing.expectEqual(-0x3.d7d0ap+4, lgamma(@as(f32, -0x1.dffffep+4)));
    try std.testing.expectEqual(-0x3.d7d0acp+4, lgamma(@as(f32, -0x1.e00002p+4)));
    try std.testing.expectEqual(0x9.a8106p+0, lgamma(@as(f32, 0x8.8d2d5p+0)));
    try std.testing.expectEqual(0x3.2125f4p+56, lgamma(@as(f32, 0x1.6a324ap+52)));
    try std.testing.expectEqual(0xb.70d43p+0, lgamma(@as(f32, 0x9.62f59p+0)));
    try std.testing.expectEqual(0xe.b6cd6p+0, lgamma(@as(f32, 0xa.d55d7p+0)));
    try std.testing.expectEqual(0xe.b6cd4p+0, lgamma(@as(f32, 0xa.d55d6p+0)));
    try std.testing.expectEqual(0xa.41bp+0, lgamma(@as(f32, 0x8.d6315p+0)));
    try std.testing.expectEqual(0xf.88427p+0, lgamma(@as(f32, 0xb.2e679p+0)));
    try std.testing.expectEqual(0xf.1d4fdp+0, lgamma(@as(f32, 0xb.01191p+0)));
    try std.testing.expectEqual(0xf.76b51p+0, lgamma(@as(f32, 0xb.26fdap+0)));
    try std.testing.expectEqual(0xf.cbb4fp+0, lgamma(@as(f32, 0xb.4ad0ap+0)));
    try std.testing.expectEqual(0xe.0ed27p+24, lgamma(@as(f32, 0xe.7a678p+20)));
    try std.testing.expectEqual(0x1.d9db4cp+0, lgamma(@as(f32, -0x2.dea4ccp-4)));
    try std.testing.expectEqual(0x1.da47d6p+0, lgamma(@as(f32, -0x2.dd306p-4)));
    try std.testing.expectEqual(0xf.f273ep-4, lgamma(@as(f32, -0x1.bdc8bp+0)));
    try std.testing.expectEqual(0x1.950848p+0, lgamma(@as(f32, -0x4.0a82e8p-4)));
    try std.testing.expectEqual(0xf.ccp-4, lgamma(@as(f32, -0x1.bca67ap+0)));
    try std.testing.expectEqual(-0xb.a18b3p-4, lgamma(@as(f32, -0x3.464468p+0)));
    try std.testing.expectEqual(-0xb.a18c3p-4, lgamma(@as(f32, -0x3.46446cp+0)));
    try std.testing.expectEqual(-0xe.aa753p-8, lgamma(@as(f32, -0x3.f3d2c4p+0)));
    try std.testing.expectEqual(-0xe.aa27bp-8, lgamma(@as(f32, -0x3.f3d2c8p+0)));

    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), lgamma(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x0p+0, lgamma(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0xb.17217f7d1cf78p-4, lgamma(@as(f64, 0x3p+0)));
    try std.testing.expectEqual(0x9.28682473d0de8p-4, lgamma(@as(f64, 0x8p-4)));
    try std.testing.expectEqual(0x4.2c8312a971bcp-4, lgamma(@as(f64, 0xb.33334p-4)));
    // try std.testing.expectEqual(0x4.2c83262ea9194p-4, lgamma(@as(f64, 0xb.33333p-4)));
    try std.testing.expectEqual(0x4.2c832247379c4p-4, lgamma(@as(f64, 0xb.3333333333338p-4)));
    // try std.testing.expectEqual(0x4.2c832247379ccp-4, lgamma(@as(f64, 0xb.333333333333p-4)));
    try std.testing.expectEqual(-0x1.5db13c7af7432p-4, lgamma(@as(f64, 0x1.333334p+0)));
    try std.testing.expectEqual(-0x1.5db1333b26a22p-4, lgamma(@as(f64, 0x1.333332p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70cbp-4, lgamma(@as(f64, 0x1.3333333333334p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c6p-4, lgamma(@as(f64, 0x1.3333333333333p+0)));
    // try std.testing.expectEqual(0x8.8bdd41bf44848p+60, lgamma(@as(f64, 0x3.8p+56)));
    try std.testing.expectEqual(0x3.72d02ef880f8cp+0, lgamma(@as(f64, 0x8p-8)));
    // try std.testing.expectEqual(0x3.7c0e0ff92f04ap+0, lgamma(@as(f64, -0x8p-8)));
    try std.testing.expectEqual(0x6.ee500bbb72644p+0, lgamma(@as(f64, 0x4p-12)));
    // try std.testing.expectEqual(0x6.ee99edf298bep+0, lgamma(@as(f64, -0x4p-12)));
    try std.testing.expectEqual(0xa.65ae3fffc5928p+0, lgamma(@as(f64, 0x2p-16)));
    // try std.testing.expectEqual(0xa.65b08f116527p+0, lgamma(@as(f64, -0x2p-16)));
    try std.testing.expectEqual(0xd.dce9d6201e8ap+0, lgamma(@as(f64, 0x1p-20)));
    // try std.testing.expectEqual(0xd.dce9e898ab868p+0, lgamma(@as(f64, -0x1p-20)));
    try std.testing.expectEqual(0x1.1542456e99b0fp+4, lgamma(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x1.15424577d5f77p+4, lgamma(@as(f64, -0x8p-28)));
    // try std.testing.expectEqual(0x1.4cb5ecf08473fp+4, lgamma(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x1.4cb5ecf0ce562p+4, lgamma(@as(f64, -0x4p-32)));
    try std.testing.expectEqual(0x1.bb9d3beb8c7d7p+4, lgamma(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x1.bb9d3beb8c8ffp+4, lgamma(@as(f64, -0x1p-40)));
    try std.testing.expectEqual(0x2.2a848ae66fa86p+4, lgamma(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x2.2a848ae66fa86p+4, lgamma(@as(f64, -0x4p-52)));
    try std.testing.expectEqual(0x2.996bd9e152cap+4, lgamma(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x2.996bd9e152cap+4, lgamma(@as(f64, -0x1p-60)));
    try std.testing.expectEqual(0x2.c5c85fdf473dep+4, lgamma(@as(f64, 0x1p-64)));
    try std.testing.expectEqual(0x2.c5c85fdf473dep+4, lgamma(@as(f64, -0x1p-64)));
    try std.testing.expectEqual(0x3.085328dc35ebcp+4, lgamma(@as(f64, 0x4p-72)));
    try std.testing.expectEqual(0x3.085328dc35ebcp+4, lgamma(@as(f64, -0x4p-72)));
    try std.testing.expectEqual(0x4.550915ccdf50cp+4, lgamma(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x4.550915ccdf50cp+4, lgamma(@as(f64, -0x1p-100)));
    try std.testing.expectEqual(0x5.75627cbf9441cp+4, lgamma(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x5.75627cbf9441cp+4, lgamma(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x8.aa122b99bea18p+4, lgamma(@as(f64, 0x1p-200)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x8.aa122b99bea18p+4, lgamma(@as(f64, -0x1p-200)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.5a92d6d005c94p+8, lgamma(@as(f64, 0x1p-500)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.5a92d6d005c94p+8, lgamma(@as(f64, -0x1p-500)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x2.b525ada00b928p+8, lgamma(@as(f64, 0x1p-1000)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.b525ada00b928p+8, lgamma(@as(f64, -0x1p-1000)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x2.c4657baf579a4p+8, lgamma(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.c4657baf579a4p+8, lgamma(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, -0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dcp+4, lgamma(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386p+8, lgamma(@as(f64, -0x4p-1076)));
    // try std.testing.expectEqual(-0x7.d809ecd340fcp-4, lgamma(@as(f64, -0x3.ec4298p+0)));
    try std.testing.expectEqual(0xf.ffff142236928p+124, lgamma(@as(f64, 0x3.12be0cp+120)));
    try std.testing.expectEqual(0x1.00000ceb5ee8ap+128, lgamma(@as(f64, 0x3.12be6p+120)));
    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff88p+1020, lgamma(@as(f64, 0x5.d53649e2d4674p+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), lgamma(@as(f64, 0x5.d53649e2d46c8p+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), lgamma(@as(f64, 0x5.d53649e2d46ap+1012)));
    try std.testing.expectEqual(0xf.ffffffffffff8p+1020, lgamma(@as(f64, 0x5.d53649e2d469cp+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), lgamma(@as(f64, 0x5.d53649e2d46ap+1012)));
    try std.testing.expectEqual(0xf.ffffffffffff8p+1020, lgamma(@as(f64, 0x5.d53649e2d469cp+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), lgamma(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), lgamma(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), lgamma(@as(f64, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x5.7b90ba32fdbcp+132, lgamma(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(std.math.inf(f64), lgamma(@as(f64, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(-0x3.511bca412890ap-20, lgamma(@as(f64, -0x3.f48e28p+0)));
    // try std.testing.expectEqual(0x1.dd4b54ca863c2p-20, lgamma(@as(f64, -0x3.f48e2cp+0)));
    try std.testing.expectEqual(-0x1.ddc0336980b58p-52, lgamma(@as(f64, -0x3.f48e2a8f85fcap+0)));
    // try std.testing.expectEqual(-0x3.4a0c544eeb21ap-24, lgamma(@as(f64, -0x3.24c1b8p+0)));
    // try std.testing.expectEqual(-0x7.78a013681f5b8p+24, lgamma(@as(f64, -0x7.fffff8p+20)));
    try std.testing.expectEqual(-0x2.30b2cde569e24p+56, lgamma(@as(f64, -0xf.ffffffffffff8p+48)));
    try std.testing.expectEqual(-0x1.52e42ff102e65p+36, lgamma(@as(f64, -0x1.000000008p+32)));
    try std.testing.expectEqual(-0x1.52e42ff265ca8p+36, lgamma(@as(f64, -0x1.000000018p+32)));
    try std.testing.expectEqual(0x1.96ee685defb2dp+0, lgamma(@as(f64, -0x4p-4)));
    try std.testing.expectEqual(0x1.43f89a3f0edd6p+0, lgamma(@as(f64, -0x8p-4)));
    try std.testing.expectEqual(0x1.93616060ea5ep+0, lgamma(@as(f64, -0xcp-4)));
    // try std.testing.expectEqual(0x1.5dce78ceba7e9p+0, lgamma(@as(f64, -0x1.4p+0)));
    // try std.testing.expectEqual(0xd.c2c0a8c107c3p-4, lgamma(@as(f64, -0x1.8p+0)));
    try std.testing.expectEqual(0x1.041e656d68578p+0, lgamma(@as(f64, -0x1.cp+0)));
    try std.testing.expectEqual(0x2.bec33c279fa7ep+0, lgamma(@as(f64, -0x2.08p+0)));
    try std.testing.expectEqual(0x2.07060e6e8471ap+0, lgamma(@as(f64, -0x2.1p+0)));
    // try std.testing.expectEqual(0x1.99a9fdaac9a14p+0, lgamma(@as(f64, -0x2.18p+0)));
    try std.testing.expectEqual(0x1.4b32e6350c0ccp+0, lgamma(@as(f64, -0x2.2p+0)));
    try std.testing.expectEqual(0x1.0e029711cf8ddp+0, lgamma(@as(f64, -0x2.28p+0)));
    // try std.testing.expectEqual(0xd.c0af3f35d3ca8p-4, lgamma(@as(f64, -0x2.3p+0)));
    try std.testing.expectEqual(0xb.214127b24186p-4, lgamma(@as(f64, -0x2.38p+0)));
    // try std.testing.expectEqual(0x8.e355968bdbc3p-4, lgamma(@as(f64, -0x2.4p+0)));
    // try std.testing.expectEqual(0x6.f371c281277c8p-4, lgamma(@as(f64, -0x2.48p+0)));
    try std.testing.expectEqual(0x5.44859a67747f4p-4, lgamma(@as(f64, -0x2.5p+0)));
    try std.testing.expectEqual(0x3.cd82f61be0058p-4, lgamma(@as(f64, -0x2.58p+0)));
    // try std.testing.expectEqual(0x2.8804abda16ecap-4, lgamma(@as(f64, -0x2.6p+0)));
    try std.testing.expectEqual(0x1.6f830ebd2f0cbp-4, lgamma(@as(f64, -0x2.68p+0)));
    try std.testing.expectEqual(0x8.0d79aed68897p-8, lgamma(@as(f64, -0x2.7p+0)));
    try std.testing.expectEqual(-0x4.60febffedb54p-8, lgamma(@as(f64, -0x2.78p+0)));
    // try std.testing.expectEqual(-0xe.65fcfaf6878bp-8, lgamma(@as(f64, -0x2.8p+0)));
    // try std.testing.expectEqual(-0x1.60773dc36dfb4p-4, lgamma(@as(f64, -0x2.88p+0)));
    // try std.testing.expectEqual(-0x1.b3f01b8343f32p-4, lgamma(@as(f64, -0x2.9p+0)));
    // try std.testing.expectEqual(-0x1.df97311d4f4d8p-4, lgamma(@as(f64, -0x2.98p+0)));
    // try std.testing.expectEqual(-0x1.e15351cbe648ep-4, lgamma(@as(f64, -0x2.ap+0)));
    try std.testing.expectEqual(-0x1.b5f70616016fbp-4, lgamma(@as(f64, -0x2.a8p+0)));
    try std.testing.expectEqual(-0x1.58f3a915176d1p-4, lgamma(@as(f64, -0x2.bp+0)));
    try std.testing.expectEqual(-0xc.3dd1386983f58p-8, lgamma(@as(f64, -0x2.b8p+0)));
    // try std.testing.expectEqual(0x1.261e6d250cf63p-8, lgamma(@as(f64, -0x2.cp+0)));
    try std.testing.expectEqual(0x1.36e062f87a4ddp-4, lgamma(@as(f64, -0x2.c8p+0)));
    try std.testing.expectEqual(0x2.bd203eea3bb2ap-4, lgamma(@as(f64, -0x2.dp+0)));
    // try std.testing.expectEqual(0x4.c3b22d7ab0718p-4, lgamma(@as(f64, -0x2.d8p+0)));
    // try std.testing.expectEqual(0x7.7e1bfe9fdd9f4p-4, lgamma(@as(f64, -0x2.ep+0)));
    // try std.testing.expectEqual(0xb.4d46adb8bb958p-4, lgamma(@as(f64, -0x2.e8p+0)));
    try std.testing.expectEqual(0x1.10b1c8eb41e02p+0, lgamma(@as(f64, -0x2.fp+0)));
    try std.testing.expectEqual(0x1.b6f672f371762p+0, lgamma(@as(f64, -0x2.f8p+0)));
    try std.testing.expectEqual(0x1.a2dd71c565b74p+0, lgamma(@as(f64, -0x3.08p+0)));
    // try std.testing.expectEqual(0xe.88018878064ap-4, lgamma(@as(f64, -0x3.1p+0)));
    // try std.testing.expectEqual(0x7.88aaf3c5b63dp-4, lgamma(@as(f64, -0x3.18p+0)));
    try std.testing.expectEqual(0x2.780ef1ecfd4bcp-4, lgamma(@as(f64, -0x3.2p+0)));
    // try std.testing.expectEqual(-0x1.83b7ade05f104p-4, lgamma(@as(f64, -0x3.28p+0)));
    // try std.testing.expectEqual(-0x4.cb8cc177ba558p-4, lgamma(@as(f64, -0x3.3p+0)));
    try std.testing.expectEqual(-0x7.92f0f0407d53cp-4, lgamma(@as(f64, -0x3.38p+0)));
    // try std.testing.expectEqual(-0x9.f86fc0dd02fp-4, lgamma(@as(f64, -0x3.4p+0)));
    // try std.testing.expectEqual(-0xc.0f85e0da3243p-4, lgamma(@as(f64, -0x3.48p+0)));
    // try std.testing.expectEqual(-0xd.e54537e890f78p-4, lgamma(@as(f64, -0x3.5p+0)));
    try std.testing.expectEqual(-0xf.82bdb76fac928p-4, lgamma(@as(f64, -0x3.58p+0)));
    // try std.testing.expectEqual(-0x1.0ee5645b59b4cp+0, lgamma(@as(f64, -0x3.6p+0)));
    // try std.testing.expectEqual(-0x1.22c983fd69436p+0, lgamma(@as(f64, -0x3.68p+0)));
    // try std.testing.expectEqual(-0x1.340abce0a1f63p+0, lgamma(@as(f64, -0x3.7p+0)));
    // try std.testing.expectEqual(-0x1.42ca4c5b0ef64p+0, lgamma(@as(f64, -0x3.78p+0)));
    try std.testing.expectEqual(-0x1.4f1b0fe64a5d8p+0, lgamma(@as(f64, -0x3.8p+0)));
    try std.testing.expectEqual(-0x1.59031291fea94p+0, lgamma(@as(f64, -0x3.88p+0)));
    // try std.testing.expectEqual(-0x1.607c0a4453978p+0, lgamma(@as(f64, -0x3.9p+0)));
    // try std.testing.expectEqual(-0x1.6572da73cb38bp+0, lgamma(@as(f64, -0x3.98p+0)));
    try std.testing.expectEqual(-0x1.67c606af08bap+0, lgamma(@as(f64, -0x3.ap+0)));
    // try std.testing.expectEqual(-0x1.6742cd4618f51p+0, lgamma(@as(f64, -0x3.a8p+0)));
    // try std.testing.expectEqual(-0x1.63a05923d4971p+0, lgamma(@as(f64, -0x3.bp+0)));
    try std.testing.expectEqual(-0x1.5c77fc83c60b4p+0, lgamma(@as(f64, -0x3.b8p+0)));
    // try std.testing.expectEqual(-0x1.513878cce057fp+0, lgamma(@as(f64, -0x3.cp+0)));
    try std.testing.expectEqual(-0x1.41106fd92d20bp+0, lgamma(@as(f64, -0x3.c8p+0)));
    // try std.testing.expectEqual(-0x1.2ac7d6f6b00a3p+0, lgamma(@as(f64, -0x3.dp+0)));
    try std.testing.expectEqual(-0x1.0c75b5ade1a5ep+0, lgamma(@as(f64, -0x3.d8p+0)));
    // try std.testing.expectEqual(-0xe.2e1c140b222ep-4, lgamma(@as(f64, -0x3.ep+0)));
    // try std.testing.expectEqual(-0xa.7fd7bc9e5b2e8p-4, lgamma(@as(f64, -0x3.e8p+0)));
    // try std.testing.expectEqual(-0x4.e2a516e3ce8c4p-4, lgamma(@as(f64, -0x3.fp+0)));
    try std.testing.expectEqual(0x5.61445b27ef2f8p-4, lgamma(@as(f64, -0x3.f8p+0)));
    try std.testing.expectEqual(-0x2.11f0445d7c7f4p+0, lgamma(@as(f64, -0x4.4p+0)));
    // try std.testing.expectEqual(-0x2.d026474418ef6p+0, lgamma(@as(f64, -0x4.8p+0)));
    try std.testing.expectEqual(-0x2.e01b099dd31eap+0, lgamma(@as(f64, -0x4.cp+0)));
    // try std.testing.expectEqual(-0x3.ba71e6fbceb68p+0, lgamma(@as(f64, -0x5.4p+0)));
    // try std.testing.expectEqual(-0x4.8490a63c2e094p+0, lgamma(@as(f64, -0x5.8p+0)));
    try std.testing.expectEqual(-0x4.9fe6996865fd8p+0, lgamma(@as(f64, -0x5.cp+0)));
    try std.testing.expectEqual(-0x5.8f95f609dcbep+0, lgamma(@as(f64, -0x6.4p+0)));
    try std.testing.expectEqual(-0x6.63bf13aa8dc4p+0, lgamma(@as(f64, -0x6.8p+0)));
    try std.testing.expectEqual(-0x6.88be607932f0cp+0, lgamma(@as(f64, -0x6.cp+0)));
    // try std.testing.expectEqual(-0x7.8ab8df93f8e2cp+0, lgamma(@as(f64, -0x7.4p+0)));
    try std.testing.expectEqual(-0x8.678fc2dc64f88p+0, lgamma(@as(f64, -0x7.8p+0)));
    try std.testing.expectEqual(-0x8.94f3f99bb4bdp+0, lgamma(@as(f64, -0x7.cp+0)));
    try std.testing.expectEqual(-0x9.a6efce3f0c5ep+0, lgamma(@as(f64, -0x8.4p+0)));
    try std.testing.expectEqual(-0xa.8b6b2323e318p+0, lgamma(@as(f64, -0x8.8p+0)));
    try std.testing.expectEqual(-0xa.c03b140e0f968p+0, lgamma(@as(f64, -0x8.cp+0)));
    try std.testing.expectEqual(-0xb.e070bc16c1b68p+0, lgamma(@as(f64, -0x9.4p+0)));
    // try std.testing.expectEqual(-0xc.cbbfcbeca7aep+0, lgamma(@as(f64, -0x9.8p+0)));
    try std.testing.expectEqual(-0xd.0736112f6db28p+0, lgamma(@as(f64, -0x9.cp+0)));
    try std.testing.expectEqual(-0xe.343934d8f3a18p+0, lgamma(@as(f64, -0xa.4p+0)));
    // try std.testing.expectEqual(-0xf.25b38682cbb5p+0, lgamma(@as(f64, -0xa.8p+0)));
    // try std.testing.expectEqual(-0xf.672fe40267958p+0, lgamma(@as(f64, -0xa.cp+0)));
    try std.testing.expectEqual(-0x1.09fd673bdc937p+4, lgamma(@as(f64, -0xb.4p+0)));
    try std.testing.expectEqual(-0x1.196f12e453063p+4, lgamma(@as(f64, -0xb.8p+0)));
    // try std.testing.expectEqual(-0x1.1ddeefa04e20ep+4, lgamma(@as(f64, -0xb.cp+0)));
    try std.testing.expectEqual(-0x1.32140999470e3p+4, lgamma(@as(f64, -0xc.4p+0)));
    try std.testing.expectEqual(-0x1.41d87554b103ap+4, lgamma(@as(f64, -0xc.8p+0)));
    try std.testing.expectEqual(-0x1.46996e9ff5e8ep+4, lgamma(@as(f64, -0xc.cp+0)));
    // try std.testing.expectEqual(-0x1.5b6c176a914d9p+4, lgamma(@as(f64, -0xd.4p+0)));
    try std.testing.expectEqual(-0x1.6b7d13453aefdp+4, lgamma(@as(f64, -0xd.8p+0)));
    try std.testing.expectEqual(-0x1.70893507e7aacp+4, lgamma(@as(f64, -0xd.cp+0)));
    try std.testing.expectEqual(-0x1.85ee2af24d7d1p+4, lgamma(@as(f64, -0xe.4p+0)));
    // try std.testing.expectEqual(-0x1.9646635d59cf1p+4, lgamma(@as(f64, -0xe.8p+0)));
    // try std.testing.expectEqual(-0x1.9b9889f00a16bp+4, lgamma(@as(f64, -0xe.cp+0)));
    try std.testing.expectEqual(-0x1.b1860b9f9cf35p+4, lgamma(@as(f64, -0xf.4p+0)));
    try std.testing.expectEqual(-0x1.c220de6eff08dp+4, lgamma(@as(f64, -0xf.8p+0)));
    try std.testing.expectEqual(-0x1.c7b48e949c3d3p+4, lgamma(@as(f64, -0xf.cp+0)));
    try std.testing.expectEqual(-0x1.de2212eef35f3p+4, lgamma(@as(f64, -0x1.04p+4)));
    try std.testing.expectEqual(-0x1.eefb6ed92d5d8p+4, lgamma(@as(f64, -0x1.08p+4)));
    try std.testing.expectEqual(-0x1.f4ccb75a4248p+4, lgamma(@as(f64, -0x1.0cp+4)));
    // try std.testing.expectEqual(-0x2.0bb2b66649904p+4, lgamma(@as(f64, -0x1.14p+4)));
    try std.testing.expectEqual(-0x2.1cc701ffd028p+4, lgamma(@as(f64, -0x1.18p+4)));
    try std.testing.expectEqual(-0x2.22d2642bdb692p+4, lgamma(@as(f64, -0x1.1cp+4)));
    // try std.testing.expectEqual(-0x2.3a2a2c33d815ep+4, lgamma(@as(f64, -0x1.24p+4)));
    // try std.testing.expectEqual(-0x2.4b76325cc89aap+4, lgamma(@as(f64, -0x1.28p+4)));
    // try std.testing.expectEqual(-0x2.51b88f97694ccp+4, lgamma(@as(f64, -0x1.2cp+4)));
    try std.testing.expectEqual(-0x2.697c23520ea4ep+4, lgamma(@as(f64, -0x1.34p+4)));
    try std.testing.expectEqual(-0x2.7afd03ae5b994p+4, lgamma(@as(f64, -0x1.38p+4)));
    try std.testing.expectEqual(-0x2.81738ebf2dd88p+4, lgamma(@as(f64, -0x1.3cp+4)));
    // try std.testing.expectEqual(-0x2.999d8a3dc8772p+4, lgamma(@as(f64, -0x1.44p+4)));
    try std.testing.expectEqual(-0x2.ab50acb9fbd4ep+4, lgamma(@as(f64, -0x1.48p+4)));
    try std.testing.expectEqual(-0x2.b1f8ddf5bf30ap+4, lgamma(@as(f64, -0x1.4cp+4)));
    try std.testing.expectEqual(-0x2.ca8460bab0c94p+4, lgamma(@as(f64, -0x1.54p+4)));
    try std.testing.expectEqual(-0x2.dc676b66a8902p+4, lgamma(@as(f64, -0x1.58p+4)));
    try std.testing.expectEqual(-0x2.e33ef7090df6p+4, lgamma(@as(f64, -0x1.5cp+4)));
    try std.testing.expectEqual(-0x2.fc27921a70bb4p+4, lgamma(@as(f64, -0x1.64p+4)));
    // try std.testing.expectEqual(-0x3.0e3860d473066p+4, lgamma(@as(f64, -0x1.68p+4)));
    try std.testing.expectEqual(-0x3.153d2f0ea92fp+4, lgamma(@as(f64, -0x1.6cp+4)));
    try std.testing.expectEqual(-0x3.2e7ed62745dbp+4, lgamma(@as(f64, -0x1.74p+4)));
    try std.testing.expectEqual(-0x3.40bb73b417caep+4, lgamma(@as(f64, -0x1.78p+4)));
    try std.testing.expectEqual(-0x3.47eb9a13a5e8ap+4, lgamma(@as(f64, -0x1.7cp+4)));
    // try std.testing.expectEqual(-0x3.6182974be0d1p+4, lgamma(@as(f64, -0x1.84p+4)));
    // try std.testing.expectEqual(-0x3.73e93790ff62ap+4, lgamma(@as(f64, -0x1.88p+4)));
    // try std.testing.expectEqual(-0x3.7b42f37904236p+4, lgamma(@as(f64, -0x1.8cp+4)));
    // try std.testing.expectEqual(-0x3.952bdce9557fcp+4, lgamma(@as(f64, -0x1.94p+4)));
    try std.testing.expectEqual(-0x3.a7bad8102447ap+4, lgamma(@as(f64, -0x1.98p+4)));
    try std.testing.expectEqual(-0x3.af3c8a0f6e392p+4, lgamma(@as(f64, -0x1.9cp+4)));
    // try std.testing.expectEqual(-0x3.c974390b28308p+4, lgamma(@as(f64, -0x1.a4p+4)));
    try std.testing.expectEqual(-0x3.dc2a0760eba4p+4, lgamma(@as(f64, -0x1.a8p+4)));
    try std.testing.expectEqual(-0x3.e3d22f3b711cap+4, lgamma(@as(f64, -0x1.acp+4)));
    try std.testing.expectEqual(-0x3.fe55b8d8334aap+4, lgamma(@as(f64, -0x1.b4p+4)));
    try std.testing.expectEqual(-0x4.1130ef485a82cp+4, lgamma(@as(f64, -0x1.b8p+4)));
    try std.testing.expectEqual(-0x4.18fe289399754p+4, lgamma(@as(f64, -0x1.bcp+4)));
    try std.testing.expectEqual(-0x4.33cad742071ep+4, lgamma(@as(f64, -0x1.c4p+4)));
    try std.testing.expectEqual(-0x4.46ca244f93cf4p+4, lgamma(@as(f64, -0x1.c8p+4)));
    // try std.testing.expectEqual(-0x4.4ebb238830304p+4, lgamma(@as(f64, -0x1.ccp+4)));
    // try std.testing.expectEqual(-0x4.69ce718eca03p+4, lgamma(@as(f64, -0x1.d4p+4)));
    // try std.testing.expectEqual(-0x4.7cf09ab73358p+4, lgamma(@as(f64, -0x1.d8p+4)));
    try std.testing.expectEqual(-0x4.85042abb5d4fcp+4, lgamma(@as(f64, -0x1.dcp+4)));
    try std.testing.expectEqual(-0x4.a05bbd6dcca64p+4, lgamma(@as(f64, -0x1.e4p+4)));
    try std.testing.expectEqual(-0x4.b39f9ce3ffeb4p+4, lgamma(@as(f64, -0x1.e8p+4)));
    // try std.testing.expectEqual(-0x4.bbd49cc22d718p+4, lgamma(@as(f64, -0x1.ecp+4)));
    // try std.testing.expectEqual(-0x4.d76e40569b13cp+4, lgamma(@as(f64, -0x1.f4p+4)));
    try std.testing.expectEqual(-0x4.ead2c3080f2ecp+4, lgamma(@as(f64, -0x1.f8p+4)));
    try std.testing.expectEqual(-0x4.f3282414b3f08p+4, lgamma(@as(f64, -0x1.fcp+4)));
    // try std.testing.expectEqual(-0x5.0f01c7fe77b5p+4, lgamma(@as(f64, -0x2.04p+4)));
    try std.testing.expectEqual(-0x5.2285ebd6e2b7cp+4, lgamma(@as(f64, -0x2.08p+4)));
    try std.testing.expectEqual(-0x5.2afaaffe44d84p+4, lgamma(@as(f64, -0x2.0cp+4)));
    try std.testing.expectEqual(-0x5.471263b9b93bcp+4, lgamma(@as(f64, -0x2.14p+4)));
    try std.testing.expectEqual(-0x5.5ab5361c05df8p+4, lgamma(@as(f64, -0x2.18p+4)));
    try std.testing.expectEqual(-0x5.63486e673f348p+4, lgamma(@as(f64, -0x2.1cp+4)));
    try std.testing.expectEqual(-0x5.7f9c5ea615044p+4, lgamma(@as(f64, -0x2.24p+4)));
    // try std.testing.expectEqual(-0x5.935cfb12d92d4p+4, lgamma(@as(f64, -0x2.28p+4)));
    // try std.testing.expectEqual(-0x5.9c0dc658a126cp+4, lgamma(@as(f64, -0x2.2cp+4)));
    try std.testing.expectEqual(-0x5.b89c3a80e9aecp+4, lgamma(@as(f64, -0x2.34p+4)));
    try std.testing.expectEqual(-0x5.cc79c963ef6b8p+4, lgamma(@as(f64, -0x2.38p+4)));
    try std.testing.expectEqual(-0x5.d547531f08744p+4, lgamma(@as(f64, -0x2.3cp+4)));
    // try std.testing.expectEqual(-0x5.f20eab1178fe4p+4, lgamma(@as(f64, -0x2.44p+4)));
    try std.testing.expectEqual(-0x6.060860b0fb0e4p+4, lgamma(@as(f64, -0x2.48p+4)));
    // try std.testing.expectEqual(-0x6.0ef1dff71ff2p+4, lgamma(@as(f64, -0x2.4cp+4)));
    try std.testing.expectEqual(-0x6.2bf09212ee61p+4, lgamma(@as(f64, -0x2.54p+4)));
    try std.testing.expectEqual(-0x6.4005ad9c060ecp+4, lgamma(@as(f64, -0x2.58p+4)));
    // try std.testing.expectEqual(-0x6.490a643105b94p+4, lgamma(@as(f64, -0x2.5cp+4)));
    try std.testing.expectEqual(-0x6.663efb8d432c4p+4, lgamma(@as(f64, -0x2.64p+4)));
    // try std.testing.expectEqual(-0x6.7a6ec639b9ba8p+4, lgamma(@as(f64, -0x2.68p+4)));
    try std.testing.expectEqual(-0x6.838dffbb1b634p+4, lgamma(@as(f64, -0x2.6cp+4)));
    try std.testing.expectEqual(-0x6.a0f71a8eb113cp+4, lgamma(@as(f64, -0x2.74p+4)));
    // try std.testing.expectEqual(-0x6.b540e6e0fb638p+4, lgamma(@as(f64, -0x2.78p+4)));
    try std.testing.expectEqual(-0x6.be79f80712a5cp+4, lgamma(@as(f64, -0x2.7cp+4)));
    try std.testing.expectEqual(-0x6.dc1646398c9cp+4, lgamma(@as(f64, -0x2.84p+4)));
    // try std.testing.expectEqual(-0x6.f0796f4c3252cp+4, lgamma(@as(f64, -0x2.88p+4)));
    // try std.testing.expectEqual(-0x6.f9cbb53dffc1cp+4, lgamma(@as(f64, -0x2.8cp+4)));
    // try std.testing.expectEqual(-0x7.1799f71c2b61p+4, lgamma(@as(f64, -0x2.94p+4)));
    // try std.testing.expectEqual(-0x7.2c15e00240c7cp+4, lgamma(@as(f64, -0x2.98p+4)));
    try std.testing.expectEqual(-0x7.3580bfb9dce54p+4, lgamma(@as(f64, -0x2.9cp+4)));
    // try std.testing.expectEqual(-0x7.537fc4c9f7584p+4, lgamma(@as(f64, -0x2.a4p+4)));
    // try std.testing.expectEqual(-0x7.6813d7fea637p+4, lgamma(@as(f64, -0x2.a8p+4)));
    // try std.testing.expectEqual(-0x7.7196bdbc4617cp+4, lgamma(@as(f64, -0x2.acp+4)));
    try std.testing.expectEqual(-0x7.8fc563ae0f088p+4, lgamma(@as(f64, -0x2.b4p+4)));
    // try std.testing.expectEqual(-0x7.a471129172194p+4, lgamma(@as(f64, -0x2.b8p+4)));
    // try std.testing.expectEqual(-0x7.ae0b715b59528p+4, lgamma(@as(f64, -0x2.bcp+4)));
    try std.testing.expectEqual(-0x7.cc68a310de778p+4, lgamma(@as(f64, -0x2.c4p+4)));
    // try std.testing.expectEqual(-0x7.e12b6570af28p+4, lgamma(@as(f64, -0x2.c8p+4)));
    // try std.testing.expectEqual(-0x7.eadcb69e9c3ap+4, lgamma(@as(f64, -0x2.ccp+4)));
    // try std.testing.expectEqual(-0x8.09676b4afe7ap+4, lgamma(@as(f64, -0x2.d4p+4)));
    try std.testing.expectEqual(-0x8.1e40bef5c77ep+4, lgamma(@as(f64, -0x2.d8p+4)));
    // try std.testing.expectEqual(-0x8.280881c698b38p+4, lgamma(@as(f64, -0x2.dcp+4)));
    try std.testing.expectEqual(-0x8.46bfbc20675dp+4, lgamma(@as(f64, -0x2.e4p+4)));
    // try std.testing.expectEqual(-0x8.5baf248219bbp+4, lgamma(@as(f64, -0x2.e8p+4)));
    try std.testing.expectEqual(-0x8.658cddba91e7p+4, lgamma(@as(f64, -0x2.ecp+4)));
    try std.testing.expectEqual(-0x8.846fab3fa6868p+4, lgamma(@as(f64, -0x2.f4p+4)));
    try std.testing.expectEqual(-0x8.9974b10693918p+4, lgamma(@as(f64, -0x2.f8p+4)));
    try std.testing.expectEqual(-0x8.a367ea98497p+4, lgamma(@as(f64, -0x2.fcp+4)));
    try std.testing.expectEqual(-0x8.c27562e15186p+4, lgamma(@as(f64, -0x3.04p+4)));
    // try std.testing.expectEqual(-0x8.d78f93aaaba48p+4, lgamma(@as(f64, -0x3.08p+4)));
    // try std.testing.expectEqual(-0x8.e197dc624cdfp+4, lgamma(@as(f64, -0x3.0cp+4)));
    try std.testing.expectEqual(-0x9.00cf208467dbp+4, lgamma(@as(f64, -0x3.14p+4)));
    // try std.testing.expectEqual(-0x9.15fe0e8f86fcp+4, lgamma(@as(f64, -0x3.18p+4)));
    try std.testing.expectEqual(-0x9.201af9c9b1dbp+4, lgamma(@as(f64, -0x3.1cp+4)));
    // try std.testing.expectEqual(-0x9.3f7b33c4bae9p+4, lgamma(@as(f64, -0x3.24p+4)));
    try std.testing.expectEqual(-0x9.54be75ac78c8p+4, lgamma(@as(f64, -0x3.28p+4)));
    try std.testing.expectEqual(-0x9.5eef9b1085f78p+4, lgamma(@as(f64, -0x3.2cp+4)));
    // try std.testing.expectEqual(-0x9.7e77fd48cb95p+4, lgamma(@as(f64, -0x3.34p+4)));
    try std.testing.expectEqual(-0x9.93cf2dc25ffa8p+4, lgamma(@as(f64, -0x3.38p+4)));
    try std.testing.expectEqual(-0x9.9e142902892b8p+4, lgamma(@as(f64, -0x3.3cp+4)));
    try std.testing.expectEqual(-0x9.bdc3edc4d93p+4, lgamma(@as(f64, -0x3.44p+4)));
    // try std.testing.expectEqual(-0x9.d32eab63afc8p+4, lgamma(@as(f64, -0x3.48p+4)));
    try std.testing.expectEqual(-0x9.dd871c0210b98p+4, lgamma(@as(f64, -0x3.4cp+4)));
    try std.testing.expectEqual(-0x9.fd5d85111f54p+4, lgamma(@as(f64, -0x3.54p+4)));
    try std.testing.expectEqual(-0xa.12db720f2fc88p+4, lgamma(@as(f64, -0x3.58p+4)));
    try std.testing.expectEqual(-0xa.1d46fb272de5p+4, lgamma(@as(f64, -0x3.5cp+4)));
    try std.testing.expectEqual(-0xa.3d43515179cbp+4, lgamma(@as(f64, -0x3.64p+4)));
    try std.testing.expectEqual(-0xa.52d4135bb8p+4, lgamma(@as(f64, -0x3.68p+4)));
    // try std.testing.expectEqual(-0xa.5d525b6f696ep+4, lgamma(@as(f64, -0x3.6cp+4)));
    try std.testing.expectEqual(-0xa.7d73ee2cd7a9p+4, lgamma(@as(f64, -0x3.74p+4)));
    // try std.testing.expectEqual(-0xa.93172e335d758p+4, lgamma(@as(f64, -0x3.78p+4)));
    try std.testing.expectEqual(-0xa.9da7defc939c8p+4, lgamma(@as(f64, -0x3.7cp+4)));
    try std.testing.expectEqual(-0xa.bdee0413128f8p+4, lgamma(@as(f64, -0x3.84p+4)));
    // try std.testing.expectEqual(-0xa.d3a36e1cae66p+4, lgamma(@as(f64, -0x3.88p+4)));
    try std.testing.expectEqual(-0xa.de46346151a98p+4, lgamma(@as(f64, -0x3.8cp+4)));
    try std.testing.expectEqual(-0xa.feb0478fe5788p+4, lgamma(@as(f64, -0x3.94p+4)));
    try std.testing.expectEqual(-0xb.14778a90c23ep+4, lgamma(@as(f64, -0x3.98p+4)));
    // try std.testing.expectEqual(-0xb.1f2c15fa353b8p+4, lgamma(@as(f64, -0x3.9cp+4)));
    try std.testing.expectEqual(-0xb.3fb978a9e018p+4, lgamma(@as(f64, -0x3.a4p+4)));
    try std.testing.expectEqual(-0xb.5592465d023f8p+4, lgamma(@as(f64, -0x3.a8p+4)));
    try std.testing.expectEqual(-0xb.605849524a7p+4, lgamma(@as(f64, -0x3.acp+4)));
    // try std.testing.expectEqual(-0xb.8108624c51a7p+4, lgamma(@as(f64, -0x3.b4p+4)));
    try std.testing.expectEqual(-0xb.96f26f0fac7cp+4, lgamma(@as(f64, -0x3.b8p+4)));
    // try std.testing.expectEqual(-0xb.a1c99e9224b88p+4, lgamma(@as(f64, -0x3.bcp+4)));
    try std.testing.expectEqual(-0xb.c29bd9bb401fp+4, lgamma(@as(f64, -0x3.c4p+4)));
    // try std.testing.expectEqual(-0xb.d896dc6e2c3cp+4, lgamma(@as(f64, -0x3.c8p+4)));
    // try std.testing.expectEqual(-0xb.e37eeff88b8ep+4, lgamma(@as(f64, -0x3.ccp+4)));
    try std.testing.expectEqual(0x1.0a2b23fa7e70dp+4, lgamma(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff4p+4, lgamma(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0xf.f140266b6279p+0, lgamma(@as(f64, -0x1.000002p+0)));
    try std.testing.expectEqual(0x2.40b2cde569e24p+4, lgamma(@as(f64, -0x1.0000000000001p+0)));
    try std.testing.expectEqual(0xf.3fce11247f0a8p+0, lgamma(@as(f64, -0x1.fffffep+0)));
    try std.testing.expectEqual(0x2.359bac65ecc56p+4, lgamma(@as(f64, -0x1.fffffffffffffp+0)));
    try std.testing.expectEqual(0xe.8e5bf3a347bb8p+0, lgamma(@as(f64, -0x2.000004p+0)));
    try std.testing.expectEqual(0x2.2a848ae66fa86p+4, lgamma(@as(f64, -0x2.0000000000002p+0)));
    try std.testing.expectEqual(0xd.751d54afa9a2p+0, lgamma(@as(f64, -0x2.fffffcp+0)));
    try std.testing.expectEqual(0x2.18f0a06bc2a56p+4, lgamma(@as(f64, -0x2.ffffffffffffep+0)));
    try std.testing.expectEqual(0xd.751d4aa322368p+0, lgamma(@as(f64, -0x3.000004p+0)));
    // try std.testing.expectEqual(0x2.18f0a06bc2a54p+4, lgamma(@as(f64, -0x3.0000000000002p+0)));
    try std.testing.expectEqual(0xc.123925c006038p+0, lgamma(@as(f64, -0x3.fffffcp+0)));
    try std.testing.expectEqual(0x2.02c25d6cc86b6p+4, lgamma(@as(f64, -0x3.ffffffffffffep+0)));
    // try std.testing.expectEqual(0xb.60c6fbb5695c8p+0, lgamma(@as(f64, -0x4.000008p+0)));
    try std.testing.expectEqual(0x1.f7ab3bed4b4e6p+4, lgamma(@as(f64, -0x4.0000000000004p+0)));
    // try std.testing.expectEqual(0x9.c4c2f5e938fb8p+0, lgamma(@as(f64, -0x4.fffff8p+0)));
    try std.testing.expectEqual(0x1.ddeaf9f55dc14p+4, lgamma(@as(f64, -0x4.ffffffffffffcp+0)));
    // try std.testing.expectEqual(0x9.c4c2da9cf6f1p+0, lgamma(@as(f64, -0x5.000008p+0)));
    try std.testing.expectEqual(0x1.ddeaf9f55dc13p+4, lgamma(@as(f64, -0x5.0000000000004p+0)));
    try std.testing.expectEqual(0x7.fa12379bec518p+0, lgamma(@as(f64, -0x5.fffff8p+0)));
    // try std.testing.expectEqual(0x1.c13fedfb33a14p+4, lgamma(@as(f64, -0x5.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x7.fa1219a4ff9c8p+0, lgamma(@as(f64, -0x6.000008p+0)));
    // try std.testing.expectEqual(0x1.c13fedfb33a13p+4, lgamma(@as(f64, -0x6.0000000000004p+0)));
    try std.testing.expectEqual(0x6.07eb0ddd58f5cp+0, lgamma(@as(f64, -0x6.fffff8p+0)));
    try std.testing.expectEqual(0x1.a21d7b4d0146ep+4, lgamma(@as(f64, -0x6.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x6.07eaed9d47ae8p+0, lgamma(@as(f64, -0x7.000008p+0)));
    // try std.testing.expectEqual(0x1.a21d7b4d0146dp+4, lgamma(@as(f64, -0x7.0000000000004p+0)));
    try std.testing.expectEqual(0x3.f394c6f5e387cp+0, lgamma(@as(f64, -0x7.fffff8p+0)));
    try std.testing.expectEqual(0x1.80d816ce89fp+4, lgamma(@as(f64, -0x7.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x3.42227b9df8fep+0, lgamma(@as(f64, -0x8.00001p+0)));
    try std.testing.expectEqual(0x1.75c0f54f0cd2fp+4, lgamma(@as(f64, -0x8.0000000000008p+0)));
    try std.testing.expectEqual(0x1.0fa5728f979e9p+0, lgamma(@as(f64, -0x8.fffffp+0)));
    try std.testing.expectEqual(0x1.52992059b2cdp+4, lgamma(@as(f64, -0x8.ffffffffffff8p+0)));
    // try std.testing.expectEqual(0x1.0fa52a813c2c7p+0, lgamma(@as(f64, -0x9.00001p+0)));
    // try std.testing.expectEqual(0x1.52992059b2ccep+4, lgamma(@as(f64, -0x9.0000000000008p+0)));
    try std.testing.expectEqual(-0x1.3dd0c34d79694p+0, lgamma(@as(f64, -0x9.fffffp+0)));
    try std.testing.expectEqual(0x1.2dc1bce24822dp+4, lgamma(@as(f64, -0x9.ffffffffffff8p+0)));
    try std.testing.expectEqual(-0x1.3dd10e8f080e9p+0, lgamma(@as(f64, -0xa.00001p+0)));
    try std.testing.expectEqual(0x1.2dc1bce24822bp+4, lgamma(@as(f64, -0xa.0000000000008p+0)));
    try std.testing.expectEqual(-0x3.a3ad38c9033a6p+0, lgamma(@as(f64, -0xa.fffffp+0)));
    try std.testing.expectEqual(0x1.0763f57349b44p+4, lgamma(@as(f64, -0xa.ffffffffffff8p+0)));
    try std.testing.expectEqual(-0x3.a3ad86f34c0e4p+0, lgamma(@as(f64, -0xb.00001p+0)));
    try std.testing.expectEqual(0x1.0763f57349b41p+4, lgamma(@as(f64, -0xb.0000000000008p+0)));
    try std.testing.expectEqual(-0x6.1fd00f0e21b3cp+0, lgamma(@as(f64, -0xb.fffffp+0)));
    try std.testing.expectEqual(0xd.fa1c7f9a2774p+0, lgamma(@as(f64, -0xb.ffffffffffff8p+0)));
    // try std.testing.expectEqual(-0x6.1fd05fe315324p+0, lgamma(@as(f64, -0xc.00001p+0)));
    // try std.testing.expectEqual(0xd.fa1c7f9a27718p+0, lgamma(@as(f64, -0xc.0000000000008p+0)));
    try std.testing.expectEqual(-0x8.b07093393f8cp+0, lgamma(@as(f64, -0xc.fffffp+0)));
    try std.testing.expectEqual(0xb.697bfa33f5eap+0, lgamma(@as(f64, -0xc.ffffffffffff8p+0)));
    try std.testing.expectEqual(-0x8.b070e6845a6dp+0, lgamma(@as(f64, -0xd.00001p+0)));
    try std.testing.expectEqual(0xb.697bfa33f5e78p+0, lgamma(@as(f64, -0xd.0000000000008p+0)));
    try std.testing.expectEqual(-0xb.5409d4efa4b7p+0, lgamma(@as(f64, -0xd.fffffp+0)));
    try std.testing.expectEqual(0x8.c5e2b758fe75p+0, lgamma(@as(f64, -0xd.ffffffffffff8p+0)));
    // try std.testing.expectEqual(-0xb.540a2a83e42a8p+0, lgamma(@as(f64, -0xe.00001p+0)));
    try std.testing.expectEqual(0x8.c5e2b758fe728p+0, lgamma(@as(f64, -0xe.0000000000008p+0)));
    // try std.testing.expectEqual(-0xe.094c9b083ca98p+0, lgamma(@as(f64, -0xe.fffffp+0)));
    try std.testing.expectEqual(0x6.109ff02f55714p+0, lgamma(@as(f64, -0xe.ffffffffffff8p+0)));
    try std.testing.expectEqual(-0xe.094cf2be9e3e8p+0, lgamma(@as(f64, -0xf.00001p+0)));
    // try std.testing.expectEqual(0x6.109ff02f556e8p+0, lgamma(@as(f64, -0xf.0000000000008p+0)));
    try std.testing.expectEqual(-0x1.0cf14f9e783e7p+4, lgamma(@as(f64, -0xf.fffffp+0)));
    // try std.testing.expectEqual(0x3.4ad790500e338p+0, lgamma(@as(f64, -0xf.ffffffffffff8p+0)));
    // try std.testing.expectEqual(-0x1.180879870e33ep+4, lgamma(@as(f64, -0x1.000002p+4)));
    try std.testing.expectEqual(0x2.996578583c5fcp+0, lgamma(@as(f64, -0x1.0000000000001p+4)));
    try std.testing.expectEqual(-0x1.455d45b618e1fp+4, lgamma(@as(f64, -0x1.0ffffep+4)));
    // try std.testing.expectEqual(-0x3.be7ffe71389ccp-4, lgamma(@as(f64, -0x1.0ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x1.455d51292150ep+4, lgamma(@as(f64, -0x1.100002p+4)));
    try std.testing.expectEqual(-0x3.be7ffe7138f86p-4, lgamma(@as(f64, -0x1.1000000000001p+4)));
    // try std.testing.expectEqual(-0x1.739c3c0e7e3dcp+4, lgamma(@as(f64, -0x1.1ffffep+4)));
    // try std.testing.expectEqual(-0x3.1fd7673485ba8p+0, lgamma(@as(f64, -0x1.1ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.739c47ba6a3afp+4, lgamma(@as(f64, -0x1.200002p+4)));
    try std.testing.expectEqual(-0x3.1fd7673485c06p+0, lgamma(@as(f64, -0x1.2000000000001p+4)));
    // try std.testing.expectEqual(-0x1.a2b8a7ff951d5p+4, lgamma(@as(f64, -0x1.2ffffep+4)));
    try std.testing.expectEqual(-0x6.119e27f51c2p+0, lgamma(@as(f64, -0x1.2ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.a2b8b3e16627ep+4, lgamma(@as(f64, -0x1.300002p+4)));
    try std.testing.expectEqual(-0x6.119e27f51c26p+0, lgamma(@as(f64, -0x1.3000000000001p+4)));
    // try std.testing.expectEqual(-0x1.d2a72cdce34acp+4, lgamma(@as(f64, -0x1.3ffffep+4)));
    try std.testing.expectEqual(-0x9.108677639892p+0, lgamma(@as(f64, -0x1.3ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.d2a738f1e7889p+4, lgamma(@as(f64, -0x1.400002p+4)));
    try std.testing.expectEqual(-0x9.108677639898p+0, lgamma(@as(f64, -0x1.4000000000001p+4)));
    try std.testing.expectEqual(-0x2.035d89ed6122p+4, lgamma(@as(f64, -0x1.4ffffep+4)));
    // try std.testing.expectEqual(-0xc.1bec49f18e68p+0, lgamma(@as(f64, -0x1.4ffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.035d9633286cp+4, lgamma(@as(f64, -0x1.500002p+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e6e8p+0, lgamma(@as(f64, -0x1.5000000000001p+4)));
    try std.testing.expectEqual(-0x2.34d272c496dcp+4, lgamma(@as(f64, -0x1.5ffffep+4)));
    try std.testing.expectEqual(-0xf.333ad8d94721p+0, lgamma(@as(f64, -0x1.5ffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.34d27f38e9c8ep+4, lgamma(@as(f64, -0x1.600002p+4)));
    // try std.testing.expectEqual(-0xf.333ad8d947278p+0, lgamma(@as(f64, -0x1.6000000000001p+4)));
    // try std.testing.expectEqual(-0x2.66fd6ea9f77b8p+4, lgamma(@as(f64, -0x1.6ffffep+4)));
    try std.testing.expectEqual(-0x1.255ea98937d9fp+4, lgamma(@as(f64, -0x1.6ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x2.66fd7b4acff92p+4, lgamma(@as(f64, -0x1.700002p+4)));
    try std.testing.expectEqual(-0x1.255ea98937da5p+4, lgamma(@as(f64, -0x1.7000000000001p+4)));
    try std.testing.expectEqual(-0x2.99d6bd8dc68p+4, lgamma(@as(f64, -0x1.7ffffep+4)));
    try std.testing.expectEqual(-0x1.5837f8825c33ep+4, lgamma(@as(f64, -0x1.7ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x2.99d6ca5949a84p+4, lgamma(@as(f64, -0x1.800002p+4)));
    // try std.testing.expectEqual(-0x1.5837f8825c345p+4, lgamma(@as(f64, -0x1.8000000000001p+4)));
    try std.testing.expectEqual(-0x2.cd57416926b92p+4, lgamma(@as(f64, -0x1.8ffffep+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374e5p+4, lgamma(@as(f64, -0x1.8ffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.cd574e5d9fa3ep+4, lgamma(@as(f64, -0x1.900002p+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374ebp+4, lgamma(@as(f64, -0x1.9000000000001p+4)));
    try std.testing.expectEqual(-0x3.01786b2b55b3ap+4, lgamma(@as(f64, -0x1.9ffffep+4)));
    try std.testing.expectEqual(-0x1.bfd9a6481783ep+4, lgamma(@as(f64, -0x1.9ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x3.0178784731148p+4, lgamma(@as(f64, -0x1.a00002p+4)));
    // try std.testing.expectEqual(-0x1.bfd9a64817845p+4, lgamma(@as(f64, -0x1.a000000000001p+4)));
    try std.testing.expectEqual(-0x3.36342a886637ep+4, lgamma(@as(f64, -0x1.affffep+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8dp+4, lgamma(@as(f64, -0x1.affffffffffffp+4)));
    try std.testing.expectEqual(-0x3.363437ca2ea26p+4, lgamma(@as(f64, -0x1.b00002p+4)));
    // try std.testing.expectEqual(-0x1.f49565b81e8d7p+4, lgamma(@as(f64, -0x1.b000000000001p+4)));
    try std.testing.expectEqual(-0x3.6b84e02349a7ap+4, lgamma(@as(f64, -0x1.bffffep+4)));
    try std.testing.expectEqual(-0x2.29e61b654b214p+4, lgamma(@as(f64, -0x1.bffffffffffffp+4)));
    try std.testing.expectEqual(-0x3.6b84ed89a45b2p+4, lgamma(@as(f64, -0x1.c00002p+4)));
    // try std.testing.expectEqual(-0x2.29e61b654b21cp+4, lgamma(@as(f64, -0x1.c000000000001p+4)));
    try std.testing.expectEqual(-0x3.a16551a93dea6p+4, lgamma(@as(f64, -0x1.cffffep+4)));
    // try std.testing.expectEqual(-0x2.5fc68cfce71d8p+4, lgamma(@as(f64, -0x1.cffffffffffffp+4)));
    try std.testing.expectEqual(-0x3.a1655f32e810cp+4, lgamma(@as(f64, -0x1.d00002p+4)));
    // try std.testing.expectEqual(-0x2.5fc68cfce71dep+4, lgamma(@as(f64, -0x1.d000000000001p+4)));
    try std.testing.expectEqual(-0x3.d7d09f8a44868p+4, lgamma(@as(f64, -0x1.dffffep+4)));
    try std.testing.expectEqual(-0x2.9631daeefecacp+4, lgamma(@as(f64, -0x1.dffffffffffffp+4)));
    // try std.testing.expectEqual(-0x3.d7d0ad3610cfp+4, lgamma(@as(f64, -0x1.e00002p+4)));
    // try std.testing.expectEqual(-0x2.9631daeefecb2p+4, lgamma(@as(f64, -0x1.e000000000001p+4)));
    try std.testing.expectEqual(-0x4.0ec23c0ae2bc4p+4, lgamma(@as(f64, -0x1.effffep+4)));
    // try std.testing.expectEqual(-0x2.cd23778021216p+4, lgamma(@as(f64, -0x1.effffffffffffp+4)));
    try std.testing.expectEqual(-0x4.0ec249d7b746cp+4, lgamma(@as(f64, -0x1.f00002p+4)));
    // try std.testing.expectEqual(-0x2.cd2377802121ep+4, lgamma(@as(f64, -0x1.f000000000001p+4)));
    // try std.testing.expectEqual(-0x4.4635e378544dp+4, lgamma(@as(f64, -0x1.fffffep+4)));
    // try std.testing.expectEqual(-0x3.04971efd92b24p+4, lgamma(@as(f64, -0x1.fffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.514d19db0f01p+4, lgamma(@as(f64, -0x2.000004p+4)));
    // try std.testing.expectEqual(-0x3.0fae407d0fcfep+4, lgamma(@as(f64, -0x2.0000000000002p+4)));
    // try std.testing.expectEqual(-0x4.893eafcc099b4p+4, lgamma(@as(f64, -0x2.0ffffcp+4)));
    // try std.testing.expectEqual(-0x3.479ff266bb40ap+4, lgamma(@as(f64, -0x2.0fffffffffffep+4)));
    // try std.testing.expectEqual(-0x4.893ecbe3c2344p+4, lgamma(@as(f64, -0x2.100004p+4)));
    try std.testing.expectEqual(-0x3.479ff266bb418p+4, lgamma(@as(f64, -0x2.1000000000002p+4)));
    try std.testing.expectEqual(-0x4.c1aaa8b15d99p+4, lgamma(@as(f64, -0x2.1ffffcp+4)));
    try std.testing.expectEqual(-0x3.800beb6a2d5c8p+4, lgamma(@as(f64, -0x2.1fffffffffffep+4)));
    try std.testing.expectEqual(-0x4.c1aac505526e4p+4, lgamma(@as(f64, -0x2.200004p+4)));
    // try std.testing.expectEqual(-0x3.800beb6a2d5d6p+4, lgamma(@as(f64, -0x2.2000000000002p+4)));
    // try std.testing.expectEqual(-0x4.fa8d5d3a3bac4p+4, lgamma(@as(f64, -0x2.2ffffcp+4)));
    // try std.testing.expectEqual(-0x3.b8eea0104d442p+4, lgamma(@as(f64, -0x2.2fffffffffffep+4)));
    // try std.testing.expectEqual(-0x4.fa8d79c8b429cp+4, lgamma(@as(f64, -0x2.300004p+4)));
    try std.testing.expectEqual(-0x3.b8eea0104d45p+4, lgamma(@as(f64, -0x2.3000000000002p+4)));
    // try std.testing.expectEqual(-0x5.33e375121e254p+4, lgamma(@as(f64, -0x2.3ffffcp+4)));
    // try std.testing.expectEqual(-0x3.f244b804a1842p+4, lgamma(@as(f64, -0x2.3fffffffffffep+4)));
    try std.testing.expectEqual(-0x5.33e391d97a30cp+4, lgamma(@as(f64, -0x2.400004p+4)));
    // try std.testing.expectEqual(-0x3.f244b804a185p+4, lgamma(@as(f64, -0x2.4000000000002p+4)));
    try std.testing.expectEqual(-0x5.6da9c6d2e6bb8p+4, lgamma(@as(f64, -0x2.4ffffcp+4)));
    try std.testing.expectEqual(-0x4.2c0b09e117138p+4, lgamma(@as(f64, -0x2.4fffffffffffep+4)));
    // try std.testing.expectEqual(-0x5.6da9e3d19cb94p+4, lgamma(@as(f64, -0x2.500004p+4)));
    try std.testing.expectEqual(-0x4.2c0b09e117148p+4, lgamma(@as(f64, -0x2.5000000000002p+4)));
    try std.testing.expectEqual(-0x5.a7dd54437ab8p+4, lgamma(@as(f64, -0x2.5ffffcp+4)));
    try std.testing.expectEqual(-0x4.663e976c9d97p+4, lgamma(@as(f64, -0x2.5fffffffffffep+4)));
    try std.testing.expectEqual(-0x5.a7dd717815c34p+4, lgamma(@as(f64, -0x2.600004p+4)));
    // try std.testing.expectEqual(-0x4.663e976c9d97cp+4, lgamma(@as(f64, -0x2.6000000000002p+4)));
    try std.testing.expectEqual(-0x5.e27b46fa492f8p+4, lgamma(@as(f64, -0x2.6ffffcp+4)));
    try std.testing.expectEqual(-0x4.a0dc8a3dadb28p+4, lgamma(@as(f64, -0x2.6fffffffffffep+4)));
    try std.testing.expectEqual(-0x5.e27b64636783p+4, lgamma(@as(f64, -0x2.700004p+4)));
    try std.testing.expectEqual(-0x4.a0dc8a3dadb38p+4, lgamma(@as(f64, -0x2.7000000000002p+4)));
    // try std.testing.expectEqual(-0x6.1d80ed571479cp+4, lgamma(@as(f64, -0x2.7ffffcp+4)));
    // try std.testing.expectEqual(-0x4.dbe230b41296cp+4, lgamma(@as(f64, -0x2.7fffffffffffep+4)));
    // try std.testing.expectEqual(-0x6.1d810af366008p+4, lgamma(@as(f64, -0x2.800004p+4)));
    // try std.testing.expectEqual(-0x4.dbe230b412978p+4, lgamma(@as(f64, -0x2.8000000000002p+4)));
    // try std.testing.expectEqual(-0x6.58ebb7c93810cp+4, lgamma(@as(f64, -0x2.8ffffcp+4)));
    // try std.testing.expectEqual(-0x5.174cfb3f2fef4p+4, lgamma(@as(f64, -0x2.8fffffffffffep+4)));
    // try std.testing.expectEqual(-0x6.58ebd5977d1acp+4, lgamma(@as(f64, -0x2.900004p+4)));
    // try std.testing.expectEqual(-0x5.174cfb3f2ff04p+4, lgamma(@as(f64, -0x2.9000000000002p+4)));
    try std.testing.expectEqual(-0x6.94b936593305p+4, lgamma(@as(f64, -0x2.9ffffcp+4)));
    try std.testing.expectEqual(-0x5.531a79e78c698p+4, lgamma(@as(f64, -0x2.9fffffffffffep+4)));
    try std.testing.expectEqual(-0x6.94b954583b1bp+4, lgamma(@as(f64, -0x2.a00004p+4)));
    try std.testing.expectEqual(-0x5.531a79e78c6a8p+4, lgamma(@as(f64, -0x2.a000000000002p+4)));
    // try std.testing.expectEqual(-0x6.d0e7166d8c7dcp+4, lgamma(@as(f64, -0x2.affffcp+4)));
    try std.testing.expectEqual(-0x5.8f485a13b641cp+4, lgamma(@as(f64, -0x2.afffffffffffep+4)));
    try std.testing.expectEqual(-0x6.d0e7349c35524p+4, lgamma(@as(f64, -0x2.b00004p+4)));
    // try std.testing.expectEqual(-0x5.8f485a13b642cp+4, lgamma(@as(f64, -0x2.b000000000002p+4)));
    // try std.testing.expectEqual(-0x7.0d7320c43f54cp+4, lgamma(@as(f64, -0x2.bffffcp+4)));
    try std.testing.expectEqual(-0x5.cbd46481aeea4p+4, lgamma(@as(f64, -0x2.bfffffffffffep+4)));
    try std.testing.expectEqual(-0x7.0d733f2173cc4p+4, lgamma(@as(f64, -0x2.c00004p+4)));
    try std.testing.expectEqual(-0x5.cbd46481aeeb4p+4, lgamma(@as(f64, -0x2.c000000000002p+4)));
    try std.testing.expectEqual(-0x7.4a5b379ac57cp+4, lgamma(@as(f64, -0x2.cffffcp+4)));
    try std.testing.expectEqual(-0x6.08bc7b6ef67d8p+4, lgamma(@as(f64, -0x2.cfffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.4a5b56257ccb8p+4, lgamma(@as(f64, -0x2.d00004p+4)));
    try std.testing.expectEqual(-0x6.08bc7b6ef67e8p+4, lgamma(@as(f64, -0x2.d000000000002p+4)));
    try std.testing.expectEqual(-0x7.879d54ffa3388p+4, lgamma(@as(f64, -0x2.dffffcp+4)));
    // try std.testing.expectEqual(-0x6.45fe98ea17028p+4, lgamma(@as(f64, -0x2.dfffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.879d73b6e018cp+4, lgamma(@as(f64, -0x2.e00004p+4)));
    // try std.testing.expectEqual(-0x6.45fe98ea17034p+4, lgamma(@as(f64, -0x2.e000000000002p+4)));
    // try std.testing.expectEqual(-0x7.c5378948fb918p+4, lgamma(@as(f64, -0x2.effffcp+4)));
    // try std.testing.expectEqual(-0x6.8398cd4938e3cp+4, lgamma(@as(f64, -0x2.efffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.c537a82bcb824p+4, lgamma(@as(f64, -0x2.f00004p+4)));
    // try std.testing.expectEqual(-0x6.8398cd4938e4cp+4, lgamma(@as(f64, -0x2.f000000000002p+4)));
    try std.testing.expectEqual(-0x8.0327f9ac47b3p+4, lgamma(@as(f64, -0x2.fffffcp+4)));
    // try std.testing.expectEqual(-0x6.c1893dc1da5acp+4, lgamma(@as(f64, -0x2.ffffffffffffep+4)));
    try std.testing.expectEqual(-0x8.032818b9c24e8p+4, lgamma(@as(f64, -0x3.000004p+4)));
    // try std.testing.expectEqual(-0x6.c1893dc1da5bcp+4, lgamma(@as(f64, -0x3.0000000000002p+4)));
    try std.testing.expectEqual(-0x8.416cdef3c687p+4, lgamma(@as(f64, -0x3.0ffffcp+4)));
    // try std.testing.expectEqual(-0x6.ffce231e3f0f8p+4, lgamma(@as(f64, -0x3.0fffffffffffep+4)));
    // try std.testing.expectEqual(-0x8.416cfe2b0ce38p+4, lgamma(@as(f64, -0x3.100004p+4)));
    try std.testing.expectEqual(-0x6.ffce231e3f104p+4, lgamma(@as(f64, -0x3.1000000000002p+4)));
    try std.testing.expectEqual(-0x8.8004844ea3ddp+4, lgamma(@as(f64, -0x3.1ffffcp+4)));
    // try std.testing.expectEqual(-0x7.3e65c88d9746cp+4, lgamma(@as(f64, -0x3.1fffffffffffep+4)));
    // try std.testing.expectEqual(-0x8.8004a3aedffc8p+4, lgamma(@as(f64, -0x3.200004p+4)));
    // try std.testing.expectEqual(-0x7.3e65c88d9747cp+4, lgamma(@as(f64, -0x3.2000000000002p+4)));
    try std.testing.expectEqual(-0x8.beed463931cbp+4, lgamma(@as(f64, -0x3.2ffffcp+4)));
    // try std.testing.expectEqual(-0x7.7d4e8a8c3948cp+4, lgamma(@as(f64, -0x3.2fffffffffffep+4)));
    try std.testing.expectEqual(-0x8.beed65c196128p+4, lgamma(@as(f64, -0x3.300004p+4)));
    // try std.testing.expectEqual(-0x7.7d4e8a8c3949cp+4, lgamma(@as(f64, -0x3.3000000000002p+4)));
    try std.testing.expectEqual(-0x8.fe25917adde28p+4, lgamma(@as(f64, -0x3.3ffffcp+4)));
    // try std.testing.expectEqual(-0x7.bc86d5e1969b4p+4, lgamma(@as(f64, -0x3.3fffffffffffep+4)));
    try std.testing.expectEqual(-0x8.fe25b12aa4ap+4, lgamma(@as(f64, -0x3.400004p+4)));
    // try std.testing.expectEqual(-0x7.bc86d5e1969c4p+4, lgamma(@as(f64, -0x3.4000000000002p+4)));
    // try std.testing.expectEqual(-0x9.3dabe237d03e8p+4, lgamma(@as(f64, -0x3.4ffffcp+4)));
    // try std.testing.expectEqual(-0x7.fc0d26b1db14cp+4, lgamma(@as(f64, -0x3.4fffffffffffep+4)));
    try std.testing.expectEqual(-0x9.3dac020e3b36p+4, lgamma(@as(f64, -0x3.500004p+4)));
    // try std.testing.expectEqual(-0x7.fc0d26b1db15cp+4, lgamma(@as(f64, -0x3.5000000000002p+4)));
    try std.testing.expectEqual(-0x9.7d7ec3145dep+4, lgamma(@as(f64, -0x3.5ffffcp+4)));
    // try std.testing.expectEqual(-0x8.3be007a15f3a8p+4, lgamma(@as(f64, -0x3.5fffffffffffep+4)));
    // try std.testing.expectEqual(-0x9.7d7ee310b5e1p+4, lgamma(@as(f64, -0x3.600004p+4)));
    try std.testing.expectEqual(-0x8.3be007a15f3b8p+4, lgamma(@as(f64, -0x3.6000000000002p+4)));
    // try std.testing.expectEqual(-0x9.bd9ccc68ab9a8p+4, lgamma(@as(f64, -0x3.6ffffcp+4)));
    // try std.testing.expectEqual(-0x8.7bfe11084b368p+4, lgamma(@as(f64, -0x3.6fffffffffffep+4)));
    try std.testing.expectEqual(-0x9.bd9cec8a401ep+4, lgamma(@as(f64, -0x3.700004p+4)));
    // try std.testing.expectEqual(-0x8.7bfe11084b378p+4, lgamma(@as(f64, -0x3.7000000000002p+4)));
    // try std.testing.expectEqual(-0x9.fe04a3830c278p+4, lgamma(@as(f64, -0x3.7ffffcp+4)));
    // try std.testing.expectEqual(-0x8.bc65e834f4e8p+4, lgamma(@as(f64, -0x3.7fffffffffffep+4)));
    // try std.testing.expectEqual(-0x9.fe04c3c932f38p+4, lgamma(@as(f64, -0x3.800004p+4)));
    // try std.testing.expectEqual(-0x8.bc65e834f4e9p+4, lgamma(@as(f64, -0x3.8000000000002p+4)));
    // try std.testing.expectEqual(-0xa.3eb4f9f7cb8cp+4, lgamma(@as(f64, -0x3.8ffffcp+4)));
    // try std.testing.expectEqual(-0x8.fd163ebbab51p+4, lgamma(@as(f64, -0x3.8fffffffffffep+4)));
    // try std.testing.expectEqual(-0xa.3eb51a61e0618p+4, lgamma(@as(f64, -0x3.900004p+4)));
    // try std.testing.expectEqual(-0x8.fd163ebbab52p+4, lgamma(@as(f64, -0x3.9000000000002p+4)));
    // try std.testing.expectEqual(-0xa.7fac8cfd3cecp+4, lgamma(@as(f64, -0x3.9ffffcp+4)));
    // try std.testing.expectEqual(-0x9.3e0dd1d2c46a8p+4, lgamma(@as(f64, -0x3.9fffffffffffep+4)));
    try std.testing.expectEqual(-0xa.7facad8aa134p+4, lgamma(@as(f64, -0x3.a00004p+4)));
    // try std.testing.expectEqual(-0x9.3e0dd1d2c46b8p+4, lgamma(@as(f64, -0x3.a000000000002p+4)));
    // try std.testing.expectEqual(-0xa.c0ea24d2fe738p+4, lgamma(@as(f64, -0x3.affffcp+4)));
    // try std.testing.expectEqual(-0x9.7f4b69b9e11p+4, lgamma(@as(f64, -0x3.afffffffffffep+4)));
    // try std.testing.expectEqual(-0xa.c0ea458318f88p+4, lgamma(@as(f64, -0x3.b00004p+4)));
    try std.testing.expectEqual(-0x9.7f4b69b9e1118p+4, lgamma(@as(f64, -0x3.b000000000002p+4)));
    try std.testing.expectEqual(-0xb.026c9433822c8p+4, lgamma(@as(f64, -0x3.bffffcp+4)));
    // try std.testing.expectEqual(-0x9.c0cdd92b75da8p+4, lgamma(@as(f64, -0x3.bfffffffffffep+4)));
    // try std.testing.expectEqual(-0xb.026cb505bed38p+4, lgamma(@as(f64, -0x3.c00004p+4)));
    // try std.testing.expectEqual(-0x9.c0cdd92b75db8p+4, lgamma(@as(f64, -0x3.c000000000002p+4)));
    // try std.testing.expectEqual(0x4.2b2b52b5464fp-24, lgamma(@as(f64, -0x2.74ff9p+0)));
    // try std.testing.expectEqual(-0x1.e4cf2421a71b2p-24, lgamma(@as(f64, -0x2.74ff94p+0)));
    try std.testing.expectEqual(0x4.0c8edb47fa1b4p-56, lgamma(@as(f64, -0x2.74ff92c01f0d8p+0)));
    try std.testing.expectEqual(-0x2.c7343f216ac92p-52, lgamma(@as(f64, -0x2.74ff92c01f0dap+0)));
    // try std.testing.expectEqual(-0x2.6b416efc56fe4p-24, lgamma(@as(f64, -0x2.bf682p+0)));
    // try std.testing.expectEqual(0x5.3d0a33adaf4f4p-24, lgamma(@as(f64, -0x2.bf6824p+0)));
    try std.testing.expectEqual(-0x3.0c498b9ac27bep-52, lgamma(@as(f64, -0x2.bf6821437b2p+0)));
    // try std.testing.expectEqual(0xc.7dc2985d3b448p-56, lgamma(@as(f64, -0x2.bf6821437b202p+0)));
    // try std.testing.expectEqual(0x1.bd69b50d51b15p-20, lgamma(@as(f64, -0x3.24c1b4p+0)));
    // try std.testing.expectEqual(-0x3.4a0c544eeb21ap-24, lgamma(@as(f64, -0x3.24c1b8p+0)));
    // try std.testing.expectEqual(0x7.a58178eb9e988p-52, lgamma(@as(f64, -0x3.24c1b793cb35ep+0)));
    try std.testing.expectEqual(-0x7.ead1b6ac3791cp-52, lgamma(@as(f64, -0x3.24c1b793cb36p+0)));
    // try std.testing.expectEqual(-0x3.511bca412890ap-20, lgamma(@as(f64, -0x3.f48e28p+0)));
    // try std.testing.expectEqual(0x1.dd4b54ca863c2p-20, lgamma(@as(f64, -0x3.f48e2cp+0)));
    try std.testing.expectEqual(-0x1.ddc0336980b58p-52, lgamma(@as(f64, -0x3.f48e2a8f85fcap+0)));
    // try std.testing.expectEqual(0x2.7957af96f2c1p-48, lgamma(@as(f64, -0x3.f48e2a8f85fccp+0)));
    try std.testing.expectEqual(0xa.3165c90424948p-20, lgamma(@as(f64, -0x4.0a1398p+0)));
    try std.testing.expectEqual(-0x3.33cb5626dc332p-20, lgamma(@as(f64, -0x4.0a13ap+0)));
    try std.testing.expectEqual(0x5.1a6a378191448p-48, lgamma(@as(f64, -0x4.0a139e16656p+0)));
    try std.testing.expectEqual(-0x1.982d05a2f456bp-48, lgamma(@as(f64, -0x4.0a139e1665604p+0)));
    // try std.testing.expectEqual(-0x3.02165b2aa6efp-16, lgamma(@as(f64, -0x4.fdd5d8p+0)));
    try std.testing.expectEqual(0xa.22e7861540cap-20, lgamma(@as(f64, -0x4.fdd5ep+0)));
    try std.testing.expectEqual(-0x1.8280d0ba86861p-44, lgamma(@as(f64, -0x4.fdd5de9bbabfp+0)));
    try std.testing.expectEqual(0x4.fa3d33517a9ecp-48, lgamma(@as(f64, -0x4.fdd5de9bbabf4p+0)));
    try std.testing.expectEqual(0x2.e258f12a679eep-16, lgamma(@as(f64, -0x5.021a9p+0)));
    try std.testing.expectEqual(-0xf.89066929e3b18p-20, lgamma(@as(f64, -0x5.021a98p+0)));
    try std.testing.expectEqual(0x1.867827fdc0e93p-48, lgamma(@as(f64, -0x5.021a95fc2db64p+0)));
    // try std.testing.expectEqual(-0x1.d50b5e02beb78p-44, lgamma(@as(f64, -0x5.021a95fc2db68p+0)));
    try std.testing.expectEqual(-0xf.15ee1077e22dp-16, lgamma(@as(f64, -0x5.ffa4b8p+0)));
    // try std.testing.expectEqual(0x7.4bb0ef1ad813cp-16, lgamma(@as(f64, -0x5.ffa4cp+0)));
    try std.testing.expectEqual(-0x4.2c4d3e7ff052p-44, lgamma(@as(f64, -0x5.ffa4bd647d034p+0)));
    // try std.testing.expectEqual(0x7.04ae139d3fb74p-44, lgamma(@as(f64, -0x5.ffa4bd647d038p+0)));
    // try std.testing.expectEqual(0x3.e9df593e904f8p-16, lgamma(@as(f64, -0x6.005ac8p+0)));
    try std.testing.expectEqual(-0x1.2b35eea26dc94p-12, lgamma(@as(f64, -0x6.005adp+0)));
    try std.testing.expectEqual(0xa.7dd3bd697d2c8p-44, lgamma(@as(f64, -0x6.005ac9625f23p+0)));
    try std.testing.expectEqual(-0xd.11e91b3ff8f48p-48, lgamma(@as(f64, -0x6.005ac9625f234p+0)));
    // try std.testing.expectEqual(-0x7.313b929690048p-12, lgamma(@as(f64, -0x6.fff2f8p+0)));
    try std.testing.expectEqual(0x2.a3598cd9f522ap-12, lgamma(@as(f64, -0x6.fff3p+0)));
    // try std.testing.expectEqual(-0x4.dc097be5d1cc4p-40, lgamma(@as(f64, -0x6.fff2fddae1bbcp+0)));
    // try std.testing.expectEqual(0xe.f46d8dcca9e8p-48, lgamma(@as(f64, -0x6.fff2fddae1bcp+0)));
    try std.testing.expectEqual(0x9.39801333caa38p-12, lgamma(@as(f64, -0x7.000cf8p+0)));
    // try std.testing.expectEqual(-0xa.32834623023ep-16, lgamma(@as(f64, -0x7.000dp+0)));
    try std.testing.expectEqual(0x3.89727e62d4844p-40, lgamma(@as(f64, -0x7.000cff7b7f878p+0)));
    try std.testing.expectEqual(-0x1.638f6c2b4fb95p-40, lgamma(@as(f64, -0x7.000cff7b7f87cp+0)));
    try std.testing.expectEqual(-0x4.cccb8849515a8p-8, lgamma(@as(f64, -0x7.fffe58p+0)));
    // try std.testing.expectEqual(0x1.37b05f6d428dap-12, lgamma(@as(f64, -0x7.fffe6p+0)));
    try std.testing.expectEqual(-0x2.551849c02b7ep-40, lgamma(@as(f64, -0x7.fffe5fe05673cp+0)));
    try std.testing.expectEqual(0x2.509d5b2dadf1ep-36, lgamma(@as(f64, -0x7.fffe5fe05674p+0)));
    try std.testing.expectEqual(0xc.8602745a4491p-16, lgamma(@as(f64, -0x8.0001ap+0)));
    // try std.testing.expectEqual(-0x9.9cf5dfb6141fp-8, lgamma(@as(f64, -0x8.0001bp+0)));
    try std.testing.expectEqual(0x1.34e935f3e5a5dp-36, lgamma(@as(f64, -0x8.0001a01459fc8p+0)));
    try std.testing.expectEqual(-0x3.b73909c155552p-36, lgamma(@as(f64, -0x8.0001a01459fdp+0)));
    try std.testing.expectEqual(-0x9.98ed0cd062e4p-8, lgamma(@as(f64, -0x8.ffffdp+0)));
    // try std.testing.expectEqual(0x5.e337e9ef84f0cp-4, lgamma(@as(f64, -0x8.ffffep+0)));
    // try std.testing.expectEqual(-0x5.88479ad476d48p-36, lgamma(@as(f64, -0x8.ffffd1c425e8p+0)));
    // try std.testing.expectEqual(0x2.6c3945e214p-32, lgamma(@as(f64, -0x8.ffffd1c425e88p+0)));
    try std.testing.expectEqual(0x5.e32ee82416adcp-4, lgamma(@as(f64, -0x9.00002p+0)));
    // try std.testing.expectEqual(-0x9.99c537e2b9298p-8, lgamma(@as(f64, -0x9.00003p+0)));
    try std.testing.expectEqual(0x2.5debd4969bb28p-36, lgamma(@as(f64, -0x9.00002e3bb47d8p+0)));
    // try std.testing.expectEqual(-0x2.9ee383255c86ep-32, lgamma(@as(f64, -0x9.00002e3bb47ep+0)));
    try std.testing.expectEqual(-0x1.3dd0c34d79694p+0, lgamma(@as(f64, -0x9.fffffp+0)));
    // try std.testing.expectEqual(-0x1.41334d2c3ccaap-28, lgamma(@as(f64, -0x9.fffffb606bdf8p+0)));
    // try std.testing.expectEqual(0x7.9c48d283217d8p-32, lgamma(@as(f64, -0x9.fffffb606bep+0)));
    try std.testing.expectEqual(-0x1.3dd10e8f080e9p+0, lgamma(@as(f64, -0xa.00001p+0)));
    // try std.testing.expectEqual(0x5.70ddf269e6d68p-32, lgamma(@as(f64, -0xa.0000049f93bb8p+0)));
    // try std.testing.expectEqual(-0x1.63ea466b9e05cp-28, lgamma(@as(f64, -0xa.0000049f93bcp+0)));
    try std.testing.expectEqual(-0x3.a3ad38c9033a6p+0, lgamma(@as(f64, -0xa.fffffp+0)));
    try std.testing.expectEqual(-0x1.0e8528e5ba92dp-24, lgamma(@as(f64, -0xa.ffffff9466e98p+0)));
    // try std.testing.expectEqual(0x2.205541c47450ep-28, lgamma(@as(f64, -0xa.ffffff9466eap+0)));
    try std.testing.expectEqual(-0x3.a3ad86f34c0e4p+0, lgamma(@as(f64, -0xb.00001p+0)));
    // try std.testing.expectEqual(0x7.573b06696043p-28, lgamma(@as(f64, -0xb.0000006b9915p+0)));
    // try std.testing.expectEqual(-0xb.b16d1e1508e78p-28, lgamma(@as(f64, -0xb.0000006b99158p+0)));
    try std.testing.expectEqual(-0x6.1fd00f0e21b3cp+0, lgamma(@as(f64, -0xb.fffffp+0)));
    try std.testing.expectEqual(-0xc.e27c4f01cf53p-28, lgamma(@as(f64, -0xb.fffffff708938p+0)));
    // try std.testing.expectEqual(0xd.785692eee5fd8p-24, lgamma(@as(f64, -0xb.fffffff70894p+0)));
    // try std.testing.expectEqual(-0x6.1fd05fe315324p+0, lgamma(@as(f64, -0xc.00001p+0)));
    // try std.testing.expectEqual(0xd.4b0a2023492bp-24, lgamma(@as(f64, -0xc.00000008f76cp+0)));
    try std.testing.expectEqual(-0xf.b743a42616368p-28, lgamma(@as(f64, -0xc.00000008f76c8p+0)));
    try std.testing.expectEqual(-0x8.b07093393f8cp+0, lgamma(@as(f64, -0xc.fffffp+0)));
    try std.testing.expectEqual(-0x7.316d886018814p-20, lgamma(@as(f64, -0xc.ffffffff4f6d8p+0)));
    // try std.testing.expectEqual(0x4.67d7d4d0a161p-20, lgamma(@as(f64, -0xc.ffffffff4f6ep+0)));
    try std.testing.expectEqual(-0x8.b070e6845a6dp+0, lgamma(@as(f64, -0xd.00001p+0)));
    // try std.testing.expectEqual(0x4.679e61ad5163p-20, lgamma(@as(f64, -0xd.00000000b092p+0)));
    // try std.testing.expectEqual(-0x7.31a6fbad0e0ccp-20, lgamma(@as(f64, -0xd.00000000b0928p+0)));
    try std.testing.expectEqual(-0xb.5409d4efa4b7p+0, lgamma(@as(f64, -0xd.fffffp+0)));
    // try std.testing.expectEqual(-0x5.861824905c09p-16, lgamma(@as(f64, -0xd.fffffffff363p+0)));
    // try std.testing.expectEqual(0x4.a000dfad124b4p-16, lgamma(@as(f64, -0xd.fffffffff3638p+0)));
    // try std.testing.expectEqual(-0xb.540a2a83e42a8p+0, lgamma(@as(f64, -0xe.00001p+0)));
    try std.testing.expectEqual(0x4.a0009c38d0a84p-16, lgamma(@as(f64, -0xe.000000000c9c8p+0)));
    // try std.testing.expectEqual(-0x5.861868074a4e4p-16, lgamma(@as(f64, -0xe.000000000c9dp+0)));
    // try std.testing.expectEqual(-0xe.094c9b083ca98p+0, lgamma(@as(f64, -0xe.fffffp+0)));
    try std.testing.expectEqual(-0x4.c8585a763b9d4p-12, lgamma(@as(f64, -0xe.ffffffffff288p+0)));
    try std.testing.expectEqual(0x4.bb5f60f986f8cp-12, lgamma(@as(f64, -0xe.ffffffffff29p+0)));
    try std.testing.expectEqual(-0xe.094cf2be9e3e8p+0, lgamma(@as(f64, -0xf.00001p+0)));
    try std.testing.expectEqual(0x4.bb5f60afdccccp-12, lgamma(@as(f64, -0xf.0000000000d7p+0)));
    try std.testing.expectEqual(-0x4.c8585ac011a48p-12, lgamma(@as(f64, -0xf.0000000000d78p+0)));
    try std.testing.expectEqual(-0x1.0cf14f9e783e7p+4, lgamma(@as(f64, -0xf.fffffp+0)));
    // try std.testing.expectEqual(-0xe.466b0623a18dp-12, lgamma(@as(f64, -0xf.fffffffffff28p+0)));
    // try std.testing.expectEqual(0x8.c4f2f20afce3p-8, lgamma(@as(f64, -0xf.fffffffffff3p+0)));
    // try std.testing.expectEqual(-0x1.180879870e33ep+4, lgamma(@as(f64, -0x1.000002p+4)));
    // try std.testing.expectEqual(0x8.c4f2f20ab3ffp-8, lgamma(@as(f64, -0x1.000000000000dp+4)));
    try std.testing.expectEqual(-0xa.33ca82bb399ep-8, lgamma(@as(f64, -0x1.000000000000ep+4)));
    try std.testing.expectEqual(-0x1.455d45b618e1fp+4, lgamma(@as(f64, -0x1.0ffffep+4)));
    // try std.testing.expectEqual(-0x3.be7ffe71389ccp-4, lgamma(@as(f64, -0x1.0ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x1.455d51292150ep+4, lgamma(@as(f64, -0x1.100002p+4)));
    try std.testing.expectEqual(-0x3.be7ffe7138f86p-4, lgamma(@as(f64, -0x1.1000000000001p+4)));
    // try std.testing.expectEqual(-0x1.739c3c0e7e3dcp+4, lgamma(@as(f64, -0x1.1ffffep+4)));
    // try std.testing.expectEqual(-0x3.1fd7673485ba8p+0, lgamma(@as(f64, -0x1.1ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.739c47ba6a3afp+4, lgamma(@as(f64, -0x1.200002p+4)));
    try std.testing.expectEqual(-0x3.1fd7673485c06p+0, lgamma(@as(f64, -0x1.2000000000001p+4)));
    // try std.testing.expectEqual(-0x1.a2b8a7ff951d5p+4, lgamma(@as(f64, -0x1.2ffffep+4)));
    try std.testing.expectEqual(-0x6.119e27f51c2p+0, lgamma(@as(f64, -0x1.2ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.a2b8b3e16627ep+4, lgamma(@as(f64, -0x1.300002p+4)));
    try std.testing.expectEqual(-0x6.119e27f51c26p+0, lgamma(@as(f64, -0x1.3000000000001p+4)));
    // try std.testing.expectEqual(-0x1.d2a72cdce34acp+4, lgamma(@as(f64, -0x1.3ffffep+4)));
    try std.testing.expectEqual(-0x9.108677639892p+0, lgamma(@as(f64, -0x1.3ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.d2a738f1e7889p+4, lgamma(@as(f64, -0x1.400002p+4)));
    try std.testing.expectEqual(-0x9.108677639898p+0, lgamma(@as(f64, -0x1.4000000000001p+4)));
    try std.testing.expectEqual(-0x2.035d89ed6122p+4, lgamma(@as(f64, -0x1.4ffffep+4)));
    // try std.testing.expectEqual(-0xc.1bec49f18e68p+0, lgamma(@as(f64, -0x1.4ffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.035d9633286cp+4, lgamma(@as(f64, -0x1.500002p+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e6e8p+0, lgamma(@as(f64, -0x1.5000000000001p+4)));
    try std.testing.expectEqual(-0x2.34d272c496dcp+4, lgamma(@as(f64, -0x1.5ffffep+4)));
    try std.testing.expectEqual(-0xf.333ad8d94721p+0, lgamma(@as(f64, -0x1.5ffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.34d27f38e9c8ep+4, lgamma(@as(f64, -0x1.600002p+4)));
    // try std.testing.expectEqual(-0xf.333ad8d947278p+0, lgamma(@as(f64, -0x1.6000000000001p+4)));
    // try std.testing.expectEqual(-0x2.66fd6ea9f77b8p+4, lgamma(@as(f64, -0x1.6ffffep+4)));
    try std.testing.expectEqual(-0x1.255ea98937d9fp+4, lgamma(@as(f64, -0x1.6ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x2.66fd7b4acff92p+4, lgamma(@as(f64, -0x1.700002p+4)));
    try std.testing.expectEqual(-0x1.255ea98937da5p+4, lgamma(@as(f64, -0x1.7000000000001p+4)));
    try std.testing.expectEqual(-0x2.99d6bd8dc68p+4, lgamma(@as(f64, -0x1.7ffffep+4)));
    try std.testing.expectEqual(-0x1.5837f8825c33ep+4, lgamma(@as(f64, -0x1.7ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x2.99d6ca5949a84p+4, lgamma(@as(f64, -0x1.800002p+4)));
    // try std.testing.expectEqual(-0x1.5837f8825c345p+4, lgamma(@as(f64, -0x1.8000000000001p+4)));
    try std.testing.expectEqual(-0x2.cd57416926b92p+4, lgamma(@as(f64, -0x1.8ffffep+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374e5p+4, lgamma(@as(f64, -0x1.8ffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.cd574e5d9fa3ep+4, lgamma(@as(f64, -0x1.900002p+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374ebp+4, lgamma(@as(f64, -0x1.9000000000001p+4)));
    try std.testing.expectEqual(-0x3.01786b2b55b3ap+4, lgamma(@as(f64, -0x1.9ffffep+4)));
    try std.testing.expectEqual(-0x1.bfd9a6481783ep+4, lgamma(@as(f64, -0x1.9ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x3.0178784731148p+4, lgamma(@as(f64, -0x1.a00002p+4)));
    // try std.testing.expectEqual(-0x1.bfd9a64817845p+4, lgamma(@as(f64, -0x1.a000000000001p+4)));
    try std.testing.expectEqual(-0x3.36342a886637ep+4, lgamma(@as(f64, -0x1.affffep+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8dp+4, lgamma(@as(f64, -0x1.affffffffffffp+4)));
    try std.testing.expectEqual(-0x3.363437ca2ea26p+4, lgamma(@as(f64, -0x1.b00002p+4)));
    // try std.testing.expectEqual(-0x1.f49565b81e8d7p+4, lgamma(@as(f64, -0x1.b000000000001p+4)));
    try std.testing.expectEqual(-0x3.6b84e02349a7ap+4, lgamma(@as(f64, -0x1.bffffep+4)));
    try std.testing.expectEqual(-0x2.29e61b654b214p+4, lgamma(@as(f64, -0x1.bffffffffffffp+4)));
    try std.testing.expectEqual(-0x3.6b84ed89a45b2p+4, lgamma(@as(f64, -0x1.c00002p+4)));
    // try std.testing.expectEqual(-0x2.29e61b654b21cp+4, lgamma(@as(f64, -0x1.c000000000001p+4)));
    try std.testing.expectEqual(-0x3.a16551a93dea6p+4, lgamma(@as(f64, -0x1.cffffep+4)));
    // try std.testing.expectEqual(-0x2.5fc68cfce71d8p+4, lgamma(@as(f64, -0x1.cffffffffffffp+4)));
    try std.testing.expectEqual(-0x3.a1655f32e810cp+4, lgamma(@as(f64, -0x1.d00002p+4)));
    // try std.testing.expectEqual(-0x2.5fc68cfce71dep+4, lgamma(@as(f64, -0x1.d000000000001p+4)));
    try std.testing.expectEqual(-0x3.d7d09f8a44868p+4, lgamma(@as(f64, -0x1.dffffep+4)));
    try std.testing.expectEqual(-0x2.9631daeefecacp+4, lgamma(@as(f64, -0x1.dffffffffffffp+4)));
    // try std.testing.expectEqual(-0x3.d7d0ad3610cfp+4, lgamma(@as(f64, -0x1.e00002p+4)));
    // try std.testing.expectEqual(-0x2.9631daeefecb2p+4, lgamma(@as(f64, -0x1.e000000000001p+4)));
    try std.testing.expectEqual(0x9.a81063e7978p+0, lgamma(@as(f64, 0x8.8d2d5p+0)));
    // try std.testing.expectEqual(0x3.2125f40f9a1bep+56, lgamma(@as(f64, 0x1.6a324ap+52)));
    // try std.testing.expectEqual(0xb.70d4369f5b4c8p+0, lgamma(@as(f64, 0x9.62f59p+0)));
    try std.testing.expectEqual(0xe.b6cd62d45ad4p+0, lgamma(@as(f64, 0xa.d55d7p+0)));
    try std.testing.expectEqual(0xe.b6cd3d7503be8p+0, lgamma(@as(f64, 0xa.d55d6p+0)));
    // try std.testing.expectEqual(0xe.b6cd57db84cap+0, lgamma(@as(f64, 0xa.d55d6b4d78e28p+0)));
    // try std.testing.expectEqual(0xa.41afffa8a98e8p+0, lgamma(@as(f64, 0x8.d6315p+0)));
    // try std.testing.expectEqual(0xf.8842748a38e78p+0, lgamma(@as(f64, 0xb.2e679p+0)));
    try std.testing.expectEqual(0xf.1d4fd446695d8p+0, lgamma(@as(f64, 0xb.01191p+0)));
    // try std.testing.expectEqual(0xf.76b5167078378p+0, lgamma(@as(f64, 0xb.26fdap+0)));
    // try std.testing.expectEqual(0xf.cbb4eb9c9f4ep+0, lgamma(@as(f64, 0xb.4ad0ap+0)));
    try std.testing.expectEqual(0xe.0ed26f91598ep+24, lgamma(@as(f64, 0xe.7a678p+20)));
    // try std.testing.expectEqual(0x1.d9db4ca962b42p+0, lgamma(@as(f64, -0x2.dea4ccp-4)));
    try std.testing.expectEqual(0x1.da47d6051ae6cp+0, lgamma(@as(f64, -0x2.dd306p-4)));
    // try std.testing.expectEqual(0xf.f273df313426p-4, lgamma(@as(f64, -0x1.bdc8bp+0)));
    // try std.testing.expectEqual(0x1.950848252d48cp+0, lgamma(@as(f64, -0x4.0a82e8p-4)));
    try std.testing.expectEqual(0xf.cc00043a75098p-4, lgamma(@as(f64, -0x1.bca67ap+0)));
    try std.testing.expectEqual(-0xb.a18b329b453fp-4, lgamma(@as(f64, -0x3.464468p+0)));
    // try std.testing.expectEqual(-0xb.a18c341739da8p-4, lgamma(@as(f64, -0x3.46446cp+0)));
    // try std.testing.expectEqual(-0xb.a18c21a49017p-4, lgamma(@as(f64, -0x3.46446bb6a23aap+0)));
    try std.testing.expectEqual(-0xe.aa75345fa6408p-8, lgamma(@as(f64, -0x3.f3d2c4p+0)));
    // try std.testing.expectEqual(-0xe.aa27b7e3f86d8p-8, lgamma(@as(f64, -0x3.f3d2c8p+0)));
    // try std.testing.expectEqual(-0xe.aa7484b49666p-8, lgamma(@as(f64, -0x3.f3d2c40911814p+0)));

    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd08p+1032, lgamma(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), lgamma(@as(f80, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x0p+0, lgamma(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0xb.17217f7d1cf79acp-4, lgamma(@as(f80, 0x3p+0)));
    try std.testing.expectEqual(0x9.28682473d0de85fp-4, lgamma(@as(f80, 0x8p-4)));
    try std.testing.expectEqual(0x4.2c8312a971bbf728p-4, lgamma(@as(f80, 0xb.33334p-4)));
    try std.testing.expectEqual(0x4.2c83262ea9195468p-4, lgamma(@as(f80, 0xb.33333p-4)));
    try std.testing.expectEqual(0x4.2c832247379c436p-4, lgamma(@as(f80, 0xb.3333333333338p-4)));
    try std.testing.expectEqual(0x4.2c832247379cdf9p-4, lgamma(@as(f80, 0xb.333333333333p-4)));
    try std.testing.expectEqual(0x4.2c832247379ca108p-4, lgamma(@as(f80, 0xb.333333333333334p-4)));
    try std.testing.expectEqual(0x4.2c832247379ca118p-4, lgamma(@as(f80, 0xb.333333333333333p-4)));
    try std.testing.expectEqual(-0x1.5db13c7af7431d54p-4, lgamma(@as(f80, 0x1.333334p+0)));
    try std.testing.expectEqual(-0x1.5db1333b26a21d94p-4, lgamma(@as(f80, 0x1.333332p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70cadfep-4, lgamma(@as(f80, 0x1.3333333333334p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c63fep-4, lgamma(@as(f80, 0x1.3333333333333p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c72cep-4, lgamma(@as(f80, 0x1.3333333333333334p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c72c6p-4, lgamma(@as(f80, 0x1.3333333333333332p+0)));
    try std.testing.expectEqual(0x8.8bdd41bf4484605p+60, lgamma(@as(f80, 0x3.8p+56)));
    try std.testing.expectEqual(0x3.72d02ef880f8c918p+0, lgamma(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x3.7c0e0ff92f049588p+0, lgamma(@as(f80, -0x8p-8)));
    try std.testing.expectEqual(0x6.ee500bbb72645fdp+0, lgamma(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x6.ee99edf298bdfe38p+0, lgamma(@as(f80, -0x4p-12)));
    try std.testing.expectEqual(0xa.65ae3fffc592bd6p+0, lgamma(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0xa.65b08f1165271d6p+0, lgamma(@as(f80, -0x2p-16)));
    try std.testing.expectEqual(0xd.dce9d6201e89d6cp+0, lgamma(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0xd.dce9e898ab86468p+0, lgamma(@as(f80, -0x1p-20)));
    try std.testing.expectEqual(0x1.1542456e99b0f24ap+4, lgamma(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x1.15424577d5f77082p+4, lgamma(@as(f80, -0x8p-28)));
    try std.testing.expectEqual(0x1.4cb5ecf08473ea2ap+4, lgamma(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x1.4cb5ecf0ce561e1cp+4, lgamma(@as(f80, -0x4p-32)));
    try std.testing.expectEqual(0x1.bb9d3beb8c7d73e6p+4, lgamma(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x1.bb9d3beb8c8fec74p+4, lgamma(@as(f80, -0x1p-40)));
    try std.testing.expectEqual(0x2.2a848ae66fa859e8p+4, lgamma(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x2.2a848ae66fa85e88p+4, lgamma(@as(f80, -0x4p-52)));
    try std.testing.expectEqual(0x2.996bd9e152ca0844p+4, lgamma(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x2.996bd9e152ca0844p+4, lgamma(@as(f80, -0x1p-60)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6bp+4, lgamma(@as(f80, 0x1p-64)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6bp+4, lgamma(@as(f80, -0x1p-64)));
    try std.testing.expectEqual(0x3.085328dc35ebb45p+4, lgamma(@as(f80, 0x4p-72)));
    try std.testing.expectEqual(0x3.085328dc35ebb45p+4, lgamma(@as(f80, -0x4p-72)));
    try std.testing.expectEqual(0x4.550915ccdf50b87p+4, lgamma(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x4.550915ccdf50b87p+4, lgamma(@as(f80, -0x1p-100)));
    try std.testing.expectEqual(0x5.75627cbf9441de28p+4, lgamma(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x5.75627cbf9441de28p+4, lgamma(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x8.aa122b99bea170ep+4, lgamma(@as(f80, 0x1p-200)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x8.aa122b99bea170ep+4, lgamma(@as(f80, -0x1p-200)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.5a92d6d005c939a4p+8, lgamma(@as(f80, 0x1p-500)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.5a92d6d005c939a4p+8, lgamma(@as(f80, -0x1p-500)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2.b525ada00b927348p+8, lgamma(@as(f80, 0x1p-1000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.b525ada00b927348p+8, lgamma(@as(f80, -0x1p-1000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2.c4657baf579a47bcp+8, lgamma(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.c4657baf579a47bcp+8, lgamma(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0xd.89bc642039dc406p+8, lgamma(@as(f80, 0x1p-5000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0xd.89bc642039dc406p+8, lgamma(@as(f80, -0x1p-5000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.b1378c84073b880cp+12, lgamma(@as(f80, 0x1p-10000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.b1378c84073b880cp+12, lgamma(@as(f80, -0x1p-10000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x2.c5b2319c4843accp+12, lgamma(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x2.c5b2319c4843accp+12, lgamma(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x2.c86ce2daa80dcdbp+12, lgamma(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x2.c86ce2daa80dcdbp+12, lgamma(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x2.c86ce2daa80dcdbp+12, lgamma(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d1p+4, lgamma(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c74p+8, lgamma(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x2.c86ce2daa80dcdbp+12, lgamma(@as(f80, -0x8p-16448)));
    try std.testing.expectEqual(-0x7.d809ecd340fc16d8p-4, lgamma(@as(f80, -0x3.ec4298p+0)));
    try std.testing.expectEqual(0xf.ffff14223692bc4p+124, lgamma(@as(f80, 0x3.12be0cp+120)));
    try std.testing.expectEqual(0x1.00000ceb5ee8a07p+128, lgamma(@as(f80, 0x3.12be6p+120)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0xf.fffffffffff895bp+1020, lgamma(@as(f80, 0x5.d53649e2d4674p+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.000000000000701ap+1024, lgamma(@as(f80, 0x5.d53649e2d46c8p+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.0000000000000238p+1024, lgamma(@as(f80, 0x5.d53649e2d46ap+1012)));
    try std.testing.expectEqual(0xf.ffffffffffff73cp+1020, lgamma(@as(f80, 0x5.d53649e2d469cp+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffc01p+1020, lgamma(@as(f80, 0x5.d53649e2d469dbc8p+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffcp+1020, lgamma(@as(f80, 0x5.d53649e2d469dbcp+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.0000000000000238p+1024, lgamma(@as(f80, 0x5.d53649e2d46ap+1012)));
    try std.testing.expectEqual(0xf.ffffffffffff73cp+1020, lgamma(@as(f80, 0x5.d53649e2d469cp+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffc01p+1020, lgamma(@as(f80, 0x5.d53649e2d469dbc8p+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffcp+1020, lgamma(@as(f80, 0x5.d53649e2d469dbcp+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd08p+1032, lgamma(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0xf.ffffffffffffff1p+16380, lgamma(@as(f80, 0x5.c6aa645fffef5f5p+16368)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd08p+1032, lgamma(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), lgamma(@as(f80, 0x5.c6aa645fffef5ff8p+16368)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd08p+1032, lgamma(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), lgamma(@as(f80, 0x5.c6aa645fffef5fbp+16368)));
    try std.testing.expectEqual(std.math.inf(f80), lgamma(@as(f80, 0x5.c6aa645fffef5fa8p+16368)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16dp+132, lgamma(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd08p+1032, lgamma(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f80), lgamma(@as(f80, 0x5.c6aa645fffef5fbp+16368)));
    try std.testing.expectEqual(std.math.inf(f80), lgamma(@as(f80, 0x5.c6aa645fffef5fa8p+16368)));
    try std.testing.expectEqual(-0x3.511bca412890969p-20, lgamma(@as(f80, -0x3.f48e28p+0)));
    try std.testing.expectEqual(0x1.dd4b54ca863c1a48p-20, lgamma(@as(f80, -0x3.f48e2cp+0)));
    try std.testing.expectEqual(-0x1.ddc0336980b584d2p-52, lgamma(@as(f80, -0x3.f48e2a8f85fcap+0)));
    try std.testing.expectEqual(-0x3.4a0c544eeb21a028p-24, lgamma(@as(f80, -0x3.24c1b8p+0)));
    try std.testing.expectEqual(-0x7.78a013681f5b969p+24, lgamma(@as(f80, -0x7.fffff8p+20)));
    try std.testing.expectEqual(-0x2.30b2cde569e24b34p+56, lgamma(@as(f80, -0xf.ffffffffffff8p+48)));
    try std.testing.expectEqual(-0x1.55589f2fe510778ap+68, lgamma(@as(f80, -0x7.fffffffffffffff8p+60)));
    try std.testing.expectEqual(-0x1.52e42ff102e64be2p+36, lgamma(@as(f80, -0x1.000000008p+32)));
    try std.testing.expectEqual(-0x1.52e42ff265ca7bd2p+36, lgamma(@as(f80, -0x1.000000018p+32)));
    try std.testing.expectEqual(0x1.96ee685defb2cf08p+0, lgamma(@as(f80, -0x4p-4)));
    try std.testing.expectEqual(0x1.43f89a3f0edd620ap+0, lgamma(@as(f80, -0x8p-4)));
    try std.testing.expectEqual(0x1.93616060ea5dfbc4p+0, lgamma(@as(f80, -0xcp-4)));
    try std.testing.expectEqual(0x1.5dce78ceba7e8bbp+0, lgamma(@as(f80, -0x1.4p+0)));
    try std.testing.expectEqual(0xd.c2c0a8c107c324p-4, lgamma(@as(f80, -0x1.8p+0)));
    try std.testing.expectEqual(0x1.041e656d685779d4p+0, lgamma(@as(f80, -0x1.cp+0)));
    try std.testing.expectEqual(0x2.bec33c279fa7df5p+0, lgamma(@as(f80, -0x2.08p+0)));
    try std.testing.expectEqual(0x2.07060e6e8471a488p+0, lgamma(@as(f80, -0x2.1p+0)));
    try std.testing.expectEqual(0x1.99a9fdaac9a13816p+0, lgamma(@as(f80, -0x2.18p+0)));
    try std.testing.expectEqual(0x1.4b32e6350c0cbcfcp+0, lgamma(@as(f80, -0x2.2p+0)));
    try std.testing.expectEqual(0x1.0e029711cf8dcadcp+0, lgamma(@as(f80, -0x2.28p+0)));
    try std.testing.expectEqual(0xd.c0af3f35d3ca5ffp-4, lgamma(@as(f80, -0x2.3p+0)));
    try std.testing.expectEqual(0xb.214127b24185c3ap-4, lgamma(@as(f80, -0x2.38p+0)));
    try std.testing.expectEqual(0x8.e355968bdbc2c1ap-4, lgamma(@as(f80, -0x2.4p+0)));
    try std.testing.expectEqual(0x6.f371c281277c8f58p-4, lgamma(@as(f80, -0x2.48p+0)));
    try std.testing.expectEqual(0x5.44859a67747f55dp-4, lgamma(@as(f80, -0x2.5p+0)));
    try std.testing.expectEqual(0x3.cd82f61be0057224p-4, lgamma(@as(f80, -0x2.58p+0)));
    try std.testing.expectEqual(0x2.8804abda16ec96fcp-4, lgamma(@as(f80, -0x2.6p+0)));
    try std.testing.expectEqual(0x1.6f830ebd2f0cb62cp-4, lgamma(@as(f80, -0x2.68p+0)));
    try std.testing.expectEqual(0x8.0d79aed6889706dp-8, lgamma(@as(f80, -0x2.7p+0)));
    try std.testing.expectEqual(-0x4.60febffedb540e98p-8, lgamma(@as(f80, -0x2.78p+0)));
    try std.testing.expectEqual(-0xe.65fcfaf6878ac47p-8, lgamma(@as(f80, -0x2.8p+0)));
    try std.testing.expectEqual(-0x1.60773dc36dfb3a28p-4, lgamma(@as(f80, -0x2.88p+0)));
    try std.testing.expectEqual(-0x1.b3f01b8343f3228ep-4, lgamma(@as(f80, -0x2.9p+0)));
    try std.testing.expectEqual(-0x1.df97311d4f4d7d7ap-4, lgamma(@as(f80, -0x2.98p+0)));
    try std.testing.expectEqual(-0x1.e15351cbe648e7a6p-4, lgamma(@as(f80, -0x2.ap+0)));
    try std.testing.expectEqual(-0x1.b5f70616016fabf4p-4, lgamma(@as(f80, -0x2.a8p+0)));
    try std.testing.expectEqual(-0x1.58f3a915176d0a5ep-4, lgamma(@as(f80, -0x2.bp+0)));
    try std.testing.expectEqual(-0xc.3dd1386983f5bc3p-8, lgamma(@as(f80, -0x2.b8p+0)));
    try std.testing.expectEqual(0x1.261e6d250cf634acp-8, lgamma(@as(f80, -0x2.cp+0)));
    try std.testing.expectEqual(0x1.36e062f87a4dd0cap-4, lgamma(@as(f80, -0x2.c8p+0)));
    try std.testing.expectEqual(0x2.bd203eea3bb29668p-4, lgamma(@as(f80, -0x2.dp+0)));
    try std.testing.expectEqual(0x4.c3b22d7ab0718b98p-4, lgamma(@as(f80, -0x2.d8p+0)));
    try std.testing.expectEqual(0x7.7e1bfe9fdd9f4e88p-4, lgamma(@as(f80, -0x2.ep+0)));
    try std.testing.expectEqual(0xb.4d46adb8bb95b0ep-4, lgamma(@as(f80, -0x2.e8p+0)));
    try std.testing.expectEqual(0x1.10b1c8eb41e01f2cp+0, lgamma(@as(f80, -0x2.fp+0)));
    try std.testing.expectEqual(0x1.b6f672f371761ee2p+0, lgamma(@as(f80, -0x2.f8p+0)));
    try std.testing.expectEqual(0x1.a2dd71c565b73f6ep+0, lgamma(@as(f80, -0x3.08p+0)));
    try std.testing.expectEqual(0xe.88018878064a0a8p-4, lgamma(@as(f80, -0x3.1p+0)));
    try std.testing.expectEqual(0x7.88aaf3c5b63ce8bp-4, lgamma(@as(f80, -0x3.18p+0)));
    try std.testing.expectEqual(0x2.780ef1ecfd4bca08p-4, lgamma(@as(f80, -0x3.2p+0)));
    try std.testing.expectEqual(-0x1.83b7ade05f104574p-4, lgamma(@as(f80, -0x3.28p+0)));
    try std.testing.expectEqual(-0x4.cb8cc177ba556a8p-4, lgamma(@as(f80, -0x3.3p+0)));
    try std.testing.expectEqual(-0x7.92f0f0407d53cff8p-4, lgamma(@as(f80, -0x3.38p+0)));
    try std.testing.expectEqual(-0x9.f86fc0dd02f005fp-4, lgamma(@as(f80, -0x3.4p+0)));
    try std.testing.expectEqual(-0xc.0f85e0da3242c1dp-4, lgamma(@as(f80, -0x3.48p+0)));
    try std.testing.expectEqual(-0xd.e54537e890f7a84p-4, lgamma(@as(f80, -0x3.5p+0)));
    try std.testing.expectEqual(-0xf.82bdb76fac924fcp-4, lgamma(@as(f80, -0x3.58p+0)));
    try std.testing.expectEqual(-0x1.0ee5645b59b4c5fp+0, lgamma(@as(f80, -0x3.6p+0)));
    try std.testing.expectEqual(-0x1.22c983fd69436638p+0, lgamma(@as(f80, -0x3.68p+0)));
    try std.testing.expectEqual(-0x1.340abce0a1f62ff2p+0, lgamma(@as(f80, -0x3.7p+0)));
    try std.testing.expectEqual(-0x1.42ca4c5b0ef6441ep+0, lgamma(@as(f80, -0x3.78p+0)));
    try std.testing.expectEqual(-0x1.4f1b0fe64a5d866p+0, lgamma(@as(f80, -0x3.8p+0)));
    try std.testing.expectEqual(-0x1.59031291fea94166p+0, lgamma(@as(f80, -0x3.88p+0)));
    try std.testing.expectEqual(-0x1.607c0a445397833p+0, lgamma(@as(f80, -0x3.9p+0)));
    try std.testing.expectEqual(-0x1.6572da73cb38af5p+0, lgamma(@as(f80, -0x3.98p+0)));
    try std.testing.expectEqual(-0x1.67c606af08b9f924p+0, lgamma(@as(f80, -0x3.ap+0)));
    try std.testing.expectEqual(-0x1.6742cd4618f50d22p+0, lgamma(@as(f80, -0x3.a8p+0)));
    try std.testing.expectEqual(-0x1.63a05923d49717ap+0, lgamma(@as(f80, -0x3.bp+0)));
    try std.testing.expectEqual(-0x1.5c77fc83c60b4488p+0, lgamma(@as(f80, -0x3.b8p+0)));
    try std.testing.expectEqual(-0x1.513878cce057f69ap+0, lgamma(@as(f80, -0x3.cp+0)));
    try std.testing.expectEqual(-0x1.41106fd92d20b088p+0, lgamma(@as(f80, -0x3.c8p+0)));
    try std.testing.expectEqual(-0x1.2ac7d6f6b00a2856p+0, lgamma(@as(f80, -0x3.dp+0)));
    try std.testing.expectEqual(-0x1.0c75b5ade1a5e5d8p+0, lgamma(@as(f80, -0x3.d8p+0)));
    try std.testing.expectEqual(-0xe.2e1c140b222dc37p-4, lgamma(@as(f80, -0x3.ep+0)));
    try std.testing.expectEqual(-0xa.7fd7bc9e5b2e8a7p-4, lgamma(@as(f80, -0x3.e8p+0)));
    try std.testing.expectEqual(-0x4.e2a516e3ce8c2598p-4, lgamma(@as(f80, -0x3.fp+0)));
    try std.testing.expectEqual(0x5.61445b27ef2f9af8p-4, lgamma(@as(f80, -0x3.f8p+0)));
    try std.testing.expectEqual(-0x2.11f0445d7c7f46e8p+0, lgamma(@as(f80, -0x4.4p+0)));
    try std.testing.expectEqual(-0x2.d026474418ef5fap+0, lgamma(@as(f80, -0x4.8p+0)));
    try std.testing.expectEqual(-0x2.e01b099dd31e9184p+0, lgamma(@as(f80, -0x4.cp+0)));
    try std.testing.expectEqual(-0x3.ba71e6fbceb6724cp+0, lgamma(@as(f80, -0x5.4p+0)));
    try std.testing.expectEqual(-0x4.8490a63c2e095cfp+0, lgamma(@as(f80, -0x5.8p+0)));
    try std.testing.expectEqual(-0x4.9fe6996865fd9f5p+0, lgamma(@as(f80, -0x5.cp+0)));
    try std.testing.expectEqual(-0x5.8f95f609dcbdec58p+0, lgamma(@as(f80, -0x6.4p+0)));
    try std.testing.expectEqual(-0x6.63bf13aa8dc4031p+0, lgamma(@as(f80, -0x6.8p+0)));
    try std.testing.expectEqual(-0x6.88be607932f0a858p+0, lgamma(@as(f80, -0x6.cp+0)));
    try std.testing.expectEqual(-0x7.8ab8df93f8e2d0a8p+0, lgamma(@as(f80, -0x7.4p+0)));
    try std.testing.expectEqual(-0x8.678fc2dc64f8699p+0, lgamma(@as(f80, -0x7.8p+0)));
    try std.testing.expectEqual(-0x8.94f3f99bb4bcf32p+0, lgamma(@as(f80, -0x7.cp+0)));
    try std.testing.expectEqual(-0x9.a6efce3f0c5dfdcp+0, lgamma(@as(f80, -0x8.4p+0)));
    try std.testing.expectEqual(-0xa.8b6b2323e31829cp+0, lgamma(@as(f80, -0x8.8p+0)));
    try std.testing.expectEqual(-0xa.c03b140e0f96abcp+0, lgamma(@as(f80, -0x8.cp+0)));
    try std.testing.expectEqual(-0xb.e070bc16c1b6b16p+0, lgamma(@as(f80, -0x9.4p+0)));
    try std.testing.expectEqual(-0xc.cbbfcbeca7ae3e5p+0, lgamma(@as(f80, -0x9.8p+0)));
    try std.testing.expectEqual(-0xd.0736112f6db281bp+0, lgamma(@as(f80, -0x9.cp+0)));
    try std.testing.expectEqual(-0xe.343934d8f3a1738p+0, lgamma(@as(f80, -0xa.4p+0)));
    try std.testing.expectEqual(-0xf.25b38682cbb4e36p+0, lgamma(@as(f80, -0xa.8p+0)));
    try std.testing.expectEqual(-0xf.672fe4026795b13p+0, lgamma(@as(f80, -0xa.cp+0)));
    try std.testing.expectEqual(-0x1.09fd673bdc93709cp+4, lgamma(@as(f80, -0xb.4p+0)));
    try std.testing.expectEqual(-0x1.196f12e4530636aep+4, lgamma(@as(f80, -0xb.8p+0)));
    try std.testing.expectEqual(-0x1.1ddeefa04e20d89p+4, lgamma(@as(f80, -0xb.cp+0)));
    try std.testing.expectEqual(-0x1.32140999470e301p+4, lgamma(@as(f80, -0xc.4p+0)));
    try std.testing.expectEqual(-0x1.41d87554b103a5eap+4, lgamma(@as(f80, -0xc.8p+0)));
    try std.testing.expectEqual(-0x1.46996e9ff5e8e79p+4, lgamma(@as(f80, -0xc.cp+0)));
    try std.testing.expectEqual(-0x1.5b6c176a914d9642p+4, lgamma(@as(f80, -0xd.4p+0)));
    try std.testing.expectEqual(-0x1.6b7d13453aefce14p+4, lgamma(@as(f80, -0xd.8p+0)));
    try std.testing.expectEqual(-0x1.70893507e7aac336p+4, lgamma(@as(f80, -0xd.cp+0)));
    try std.testing.expectEqual(-0x1.85ee2af24d7d0a88p+4, lgamma(@as(f80, -0xe.4p+0)));
    try std.testing.expectEqual(-0x1.9646635d59cf13f4p+4, lgamma(@as(f80, -0xe.8p+0)));
    try std.testing.expectEqual(-0x1.9b9889f00a16b6dap+4, lgamma(@as(f80, -0xe.cp+0)));
    try std.testing.expectEqual(-0x1.b1860b9f9cf34edap+4, lgamma(@as(f80, -0xf.4p+0)));
    try std.testing.expectEqual(-0x1.c220de6eff08d03cp+4, lgamma(@as(f80, -0xf.8p+0)));
    try std.testing.expectEqual(-0x1.c7b48e949c3d3428p+4, lgamma(@as(f80, -0xf.cp+0)));
    try std.testing.expectEqual(-0x1.de2212eef35f350cp+4, lgamma(@as(f80, -0x1.04p+4)));
    try std.testing.expectEqual(-0x1.eefb6ed92d5d7aa8p+4, lgamma(@as(f80, -0x1.08p+4)));
    try std.testing.expectEqual(-0x1.f4ccb75a4247fa76p+4, lgamma(@as(f80, -0x1.0cp+4)));
    try std.testing.expectEqual(-0x2.0bb2b6664990308p+4, lgamma(@as(f80, -0x1.14p+4)));
    try std.testing.expectEqual(-0x2.1cc701ffd0280dccp+4, lgamma(@as(f80, -0x1.18p+4)));
    try std.testing.expectEqual(-0x2.22d2642bdb692f18p+4, lgamma(@as(f80, -0x1.1cp+4)));
    try std.testing.expectEqual(-0x2.3a2a2c33d815da18p+4, lgamma(@as(f80, -0x1.24p+4)));
    try std.testing.expectEqual(-0x2.4b76325cc89a90ap+4, lgamma(@as(f80, -0x1.28p+4)));
    try std.testing.expectEqual(-0x2.51b88f97694cb15p+4, lgamma(@as(f80, -0x1.2cp+4)));
    try std.testing.expectEqual(-0x2.697c23520ea4d9a8p+4, lgamma(@as(f80, -0x1.34p+4)));
    try std.testing.expectEqual(-0x2.7afd03ae5b99459cp+4, lgamma(@as(f80, -0x1.38p+4)));
    try std.testing.expectEqual(-0x2.81738ebf2dd88e14p+4, lgamma(@as(f80, -0x1.3cp+4)));
    try std.testing.expectEqual(-0x2.999d8a3dc87714dp+4, lgamma(@as(f80, -0x1.44p+4)));
    try std.testing.expectEqual(-0x2.ab50acb9fbd4e958p+4, lgamma(@as(f80, -0x1.48p+4)));
    try std.testing.expectEqual(-0x2.b1f8ddf5bf30a558p+4, lgamma(@as(f80, -0x1.4cp+4)));
    try std.testing.expectEqual(-0x2.ca8460bab0c94ca4p+4, lgamma(@as(f80, -0x1.54p+4)));
    try std.testing.expectEqual(-0x2.dc676b66a89013e8p+4, lgamma(@as(f80, -0x1.58p+4)));
    try std.testing.expectEqual(-0x2.e33ef7090df5fe34p+4, lgamma(@as(f80, -0x1.5cp+4)));
    try std.testing.expectEqual(-0x2.fc27921a70bb365p+4, lgamma(@as(f80, -0x1.64p+4)));
    try std.testing.expectEqual(-0x3.0e3860d4730664e8p+4, lgamma(@as(f80, -0x1.68p+4)));
    try std.testing.expectEqual(-0x3.153d2f0ea92f085p+4, lgamma(@as(f80, -0x1.6cp+4)));
    try std.testing.expectEqual(-0x3.2e7ed62745db0594p+4, lgamma(@as(f80, -0x1.74p+4)));
    try std.testing.expectEqual(-0x3.40bb73b417cadap+4, lgamma(@as(f80, -0x1.78p+4)));
    try std.testing.expectEqual(-0x3.47eb9a13a5e8a568p+4, lgamma(@as(f80, -0x1.7cp+4)));
    try std.testing.expectEqual(-0x3.6182974be0d0f664p+4, lgamma(@as(f80, -0x1.84p+4)));
    try std.testing.expectEqual(-0x3.73e93790ff62911p+4, lgamma(@as(f80, -0x1.88p+4)));
    try std.testing.expectEqual(-0x3.7b42f379042362d4p+4, lgamma(@as(f80, -0x1.8cp+4)));
    try std.testing.expectEqual(-0x3.952bdce9557fca24p+4, lgamma(@as(f80, -0x1.94p+4)));
    try std.testing.expectEqual(-0x3.a7bad810244797a8p+4, lgamma(@as(f80, -0x1.98p+4)));
    try std.testing.expectEqual(-0x3.af3c8a0f6e392338p+4, lgamma(@as(f80, -0x1.9cp+4)));
    try std.testing.expectEqual(-0x3.c974390b28307048p+4, lgamma(@as(f80, -0x1.a4p+4)));
    try std.testing.expectEqual(-0x3.dc2a0760eba3f578p+4, lgamma(@as(f80, -0x1.a8p+4)));
    try std.testing.expectEqual(-0x3.e3d22f3b711ca1c8p+4, lgamma(@as(f80, -0x1.acp+4)));
    try std.testing.expectEqual(-0x3.fe55b8d8334a9a3p+4, lgamma(@as(f80, -0x1.b4p+4)));
    try std.testing.expectEqual(-0x4.1130ef485a82c8b8p+4, lgamma(@as(f80, -0x1.b8p+4)));
    try std.testing.expectEqual(-0x4.18fe289399753798p+4, lgamma(@as(f80, -0x1.bcp+4)));
    try std.testing.expectEqual(-0x4.33cad742071e19ep+4, lgamma(@as(f80, -0x1.c4p+4)));
    try std.testing.expectEqual(-0x4.46ca244f93cf3498p+4, lgamma(@as(f80, -0x1.c8p+4)));
    try std.testing.expectEqual(-0x4.4ebb238830305bep+4, lgamma(@as(f80, -0x1.ccp+4)));
    try std.testing.expectEqual(-0x4.69ce718eca02e1d8p+4, lgamma(@as(f80, -0x1.d4p+4)));
    try std.testing.expectEqual(-0x4.7cf09ab733581fd8p+4, lgamma(@as(f80, -0x1.d8p+4)));
    try std.testing.expectEqual(-0x4.85042abb5d4fb7ap+4, lgamma(@as(f80, -0x1.dcp+4)));
    try std.testing.expectEqual(-0x4.a05bbd6dcca6218p+4, lgamma(@as(f80, -0x1.e4p+4)));
    try std.testing.expectEqual(-0x4.b39f9ce3ffeb5bc8p+4, lgamma(@as(f80, -0x1.e8p+4)));
    try std.testing.expectEqual(-0x4.bbd49cc22d716e58p+4, lgamma(@as(f80, -0x1.ecp+4)));
    try std.testing.expectEqual(-0x4.d76e40569b13cc88p+4, lgamma(@as(f80, -0x1.f4p+4)));
    try std.testing.expectEqual(-0x4.ead2c3080f2ed0bp+4, lgamma(@as(f80, -0x1.f8p+4)));
    try std.testing.expectEqual(-0x4.f3282414b3f08778p+4, lgamma(@as(f80, -0x1.fcp+4)));
    try std.testing.expectEqual(-0x5.0f01c7fe77b50a18p+4, lgamma(@as(f80, -0x2.04p+4)));
    try std.testing.expectEqual(-0x5.2285ebd6e2b7ae78p+4, lgamma(@as(f80, -0x2.08p+4)));
    try std.testing.expectEqual(-0x5.2afaaffe44d821e8p+4, lgamma(@as(f80, -0x2.0cp+4)));
    try std.testing.expectEqual(-0x5.471263b9b93bcb18p+4, lgamma(@as(f80, -0x2.14p+4)));
    try std.testing.expectEqual(-0x5.5ab5361c05df6c6p+4, lgamma(@as(f80, -0x2.18p+4)));
    try std.testing.expectEqual(-0x5.63486e673f3485e8p+4, lgamma(@as(f80, -0x2.1cp+4)));
    try std.testing.expectEqual(-0x5.7f9c5ea615044f3p+4, lgamma(@as(f80, -0x2.24p+4)));
    try std.testing.expectEqual(-0x5.935cfb12d92d5f7p+4, lgamma(@as(f80, -0x2.28p+4)));
    try std.testing.expectEqual(-0x5.9c0dc658a126dfep+4, lgamma(@as(f80, -0x2.2cp+4)));
    try std.testing.expectEqual(-0x5.b89c3a80e9aed74p+4, lgamma(@as(f80, -0x2.34p+4)));
    try std.testing.expectEqual(-0x5.cc79c963ef6b8bbp+4, lgamma(@as(f80, -0x2.38p+4)));
    try std.testing.expectEqual(-0x5.d547531f08742a18p+4, lgamma(@as(f80, -0x2.3cp+4)));
    try std.testing.expectEqual(-0x5.f20eab1178fe58f8p+4, lgamma(@as(f80, -0x2.44p+4)));
    try std.testing.expectEqual(-0x6.060860b0fb0e2cep+4, lgamma(@as(f80, -0x2.48p+4)));
    try std.testing.expectEqual(-0x6.0ef1dff71ff1f428p+4, lgamma(@as(f80, -0x2.4cp+4)));
    try std.testing.expectEqual(-0x6.2bf09212ee611d8p+4, lgamma(@as(f80, -0x2.54p+4)));
    try std.testing.expectEqual(-0x6.4005ad9c060ea6bp+4, lgamma(@as(f80, -0x2.58p+4)));
    try std.testing.expectEqual(-0x6.490a643105b9512p+4, lgamma(@as(f80, -0x2.5cp+4)));
    try std.testing.expectEqual(-0x6.663efb8d432c3718p+4, lgamma(@as(f80, -0x2.64p+4)));
    try std.testing.expectEqual(-0x6.7a6ec639b9ba9dd8p+4, lgamma(@as(f80, -0x2.68p+4)));
    try std.testing.expectEqual(-0x6.838dffbb1b634938p+4, lgamma(@as(f80, -0x2.6cp+4)));
    try std.testing.expectEqual(-0x6.a0f71a8eb113b98p+4, lgamma(@as(f80, -0x2.74p+4)));
    try std.testing.expectEqual(-0x6.b540e6e0fb63724p+4, lgamma(@as(f80, -0x2.78p+4)));
    try std.testing.expectEqual(-0x6.be79f80712a5bap+4, lgamma(@as(f80, -0x2.7cp+4)));
    try std.testing.expectEqual(-0x6.dc1646398c9c01bp+4, lgamma(@as(f80, -0x2.84p+4)));
    try std.testing.expectEqual(-0x6.f0796f4c3252a5p+4, lgamma(@as(f80, -0x2.88p+4)));
    try std.testing.expectEqual(-0x6.f9cbb53dffc1b6dp+4, lgamma(@as(f80, -0x2.8cp+4)));
    try std.testing.expectEqual(-0x7.1799f71c2b60e7fp+4, lgamma(@as(f80, -0x2.94p+4)));
    try std.testing.expectEqual(-0x7.2c15e00240c7b3ep+4, lgamma(@as(f80, -0x2.98p+4)));
    try std.testing.expectEqual(-0x7.3580bfb9dce5422p+4, lgamma(@as(f80, -0x2.9cp+4)));
    try std.testing.expectEqual(-0x7.537fc4c9f7583cbp+4, lgamma(@as(f80, -0x2.a4p+4)));
    try std.testing.expectEqual(-0x7.6813d7fea636e348p+4, lgamma(@as(f80, -0x2.a8p+4)));
    try std.testing.expectEqual(-0x7.7196bdbc4617c0f8p+4, lgamma(@as(f80, -0x2.acp+4)));
    try std.testing.expectEqual(-0x7.8fc563ae0f0867fp+4, lgamma(@as(f80, -0x2.b4p+4)));
    try std.testing.expectEqual(-0x7.a4711291721933cp+4, lgamma(@as(f80, -0x2.b8p+4)));
    try std.testing.expectEqual(-0x7.ae0b715b59528ffp+4, lgamma(@as(f80, -0x2.bcp+4)));
    try std.testing.expectEqual(-0x7.cc68a310de776628p+4, lgamma(@as(f80, -0x2.c4p+4)));
    try std.testing.expectEqual(-0x7.e12b6570af281508p+4, lgamma(@as(f80, -0x2.c8p+4)));
    try std.testing.expectEqual(-0x7.eadcb69e9c3a0108p+4, lgamma(@as(f80, -0x2.ccp+4)));
    try std.testing.expectEqual(-0x8.09676b4afe7a218p+4, lgamma(@as(f80, -0x2.d4p+4)));
    try std.testing.expectEqual(-0x8.1e40bef5c77e16cp+4, lgamma(@as(f80, -0x2.d8p+4)));
    try std.testing.expectEqual(-0x8.280881c698b34ffp+4, lgamma(@as(f80, -0x2.dcp+4)));
    try std.testing.expectEqual(-0x8.46bfbc20675ce02p+4, lgamma(@as(f80, -0x2.e4p+4)));
    try std.testing.expectEqual(-0x8.5baf248219baddap+4, lgamma(@as(f80, -0x2.e8p+4)));
    try std.testing.expectEqual(-0x8.658cddba91e6ebdp+4, lgamma(@as(f80, -0x2.ecp+4)));
    try std.testing.expectEqual(-0x8.846fab3fa68668p+4, lgamma(@as(f80, -0x2.f4p+4)));
    try std.testing.expectEqual(-0x8.9974b1069391725p+4, lgamma(@as(f80, -0x2.f8p+4)));
    try std.testing.expectEqual(-0x8.a367ea98496fe48p+4, lgamma(@as(f80, -0x2.fcp+4)));
    try std.testing.expectEqual(-0x8.c27562e15186009p+4, lgamma(@as(f80, -0x3.04p+4)));
    try std.testing.expectEqual(-0x8.d78f93aaaba45acp+4, lgamma(@as(f80, -0x3.08p+4)));
    try std.testing.expectEqual(-0x8.e197dc624cded55p+4, lgamma(@as(f80, -0x3.0cp+4)));
    try std.testing.expectEqual(-0x9.00cf208467db157p+4, lgamma(@as(f80, -0x3.14p+4)));
    try std.testing.expectEqual(-0x9.15fe0e8f86fc0fcp+4, lgamma(@as(f80, -0x3.18p+4)));
    try std.testing.expectEqual(-0x9.201af9c9b1dafebp+4, lgamma(@as(f80, -0x3.1cp+4)));
    try std.testing.expectEqual(-0x9.3f7b33c4bae8e66p+4, lgamma(@as(f80, -0x3.24p+4)));
    try std.testing.expectEqual(-0x9.54be75ac78c7db2p+4, lgamma(@as(f80, -0x3.28p+4)));
    try std.testing.expectEqual(-0x9.5eef9b1085f7a45p+4, lgamma(@as(f80, -0x3.2cp+4)));
    try std.testing.expectEqual(-0x9.7e77fd48cb94c5ep+4, lgamma(@as(f80, -0x3.34p+4)));
    try std.testing.expectEqual(-0x9.93cf2dc25ffa932p+4, lgamma(@as(f80, -0x3.38p+4)));
    try std.testing.expectEqual(-0x9.9e142902892baa5p+4, lgamma(@as(f80, -0x3.3cp+4)));
    try std.testing.expectEqual(-0x9.bdc3edc4d92fc7p+4, lgamma(@as(f80, -0x3.44p+4)));
    try std.testing.expectEqual(-0x9.d32eab63afc830ep+4, lgamma(@as(f80, -0x3.48p+4)));
    try std.testing.expectEqual(-0x9.dd871c0210b9a67p+4, lgamma(@as(f80, -0x3.4cp+4)));
    try std.testing.expectEqual(-0x9.fd5d85111f54064p+4, lgamma(@as(f80, -0x3.54p+4)));
    try std.testing.expectEqual(-0xa.12db720f2fc8a7p+4, lgamma(@as(f80, -0x3.58p+4)));
    try std.testing.expectEqual(-0xa.1d46fb272de50cdp+4, lgamma(@as(f80, -0x3.5cp+4)));
    try std.testing.expectEqual(-0xa.3d43515179cb224p+4, lgamma(@as(f80, -0x3.64p+4)));
    try std.testing.expectEqual(-0xa.52d4135bb7ffc89p+4, lgamma(@as(f80, -0x3.68p+4)));
    try std.testing.expectEqual(-0xa.5d525b6f696dc1p+4, lgamma(@as(f80, -0x3.6cp+4)));
    try std.testing.expectEqual(-0xa.7d73ee2cd7a8c8ap+4, lgamma(@as(f80, -0x3.74p+4)));
    try std.testing.expectEqual(-0xa.93172e335d7555fp+4, lgamma(@as(f80, -0x3.78p+4)));
    try std.testing.expectEqual(-0xa.9da7defc939ca96p+4, lgamma(@as(f80, -0x3.7cp+4)));
    try std.testing.expectEqual(-0xa.bdee0413128f557p+4, lgamma(@as(f80, -0x3.84p+4)));
    try std.testing.expectEqual(-0xa.d3a36e1cae65cd4p+4, lgamma(@as(f80, -0x3.88p+4)));
    try std.testing.expectEqual(-0xa.de46346151a9651p+4, lgamma(@as(f80, -0x3.8cp+4)));
    try std.testing.expectEqual(-0xa.feb0478fe57870ep+4, lgamma(@as(f80, -0x3.94p+4)));
    try std.testing.expectEqual(-0xb.14778a90c23de92p+4, lgamma(@as(f80, -0x3.98p+4)));
    try std.testing.expectEqual(-0xb.1f2c15fa353b6f3p+4, lgamma(@as(f80, -0x3.9cp+4)));
    try std.testing.expectEqual(-0xb.3fb978a9e017ff8p+4, lgamma(@as(f80, -0x3.a4p+4)));
    try std.testing.expectEqual(-0xb.5592465d023fa8bp+4, lgamma(@as(f80, -0x3.a8p+4)));
    try std.testing.expectEqual(-0xb.605849524a70202p+4, lgamma(@as(f80, -0x3.acp+4)));
    try std.testing.expectEqual(-0xb.8108624c51a6e6dp+4, lgamma(@as(f80, -0x3.b4p+4)));
    try std.testing.expectEqual(-0xb.96f26f0fac7bfc1p+4, lgamma(@as(f80, -0x3.b8p+4)));
    try std.testing.expectEqual(-0xb.a1c99e9224b8975p+4, lgamma(@as(f80, -0x3.bcp+4)));
    try std.testing.expectEqual(-0xb.c29bd9bb401ef0ap+4, lgamma(@as(f80, -0x3.c4p+4)));
    try std.testing.expectEqual(-0xb.d896dc6e2c3c335p+4, lgamma(@as(f80, -0x3.c8p+4)));
    try std.testing.expectEqual(-0xb.e37eeff88b8ddd1p+4, lgamma(@as(f80, -0x3.ccp+4)));
    try std.testing.expectEqual(0x1.0a2b23fa7e70cd72p+4, lgamma(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff434p+4, lgamma(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6bp+4, lgamma(@as(f80, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0xf.f140266b6278ffap+0, lgamma(@as(f80, -0x1.000002p+0)));
    try std.testing.expectEqual(0x2.40b2cde569e24b04p+4, lgamma(@as(f80, -0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.bab13e5fca20ef14p+4, lgamma(@as(f80, -0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0xf.3fce11247f0a78p+0, lgamma(@as(f80, -0x1.fffffep+0)));
    try std.testing.expectEqual(0x2.359bac65ecc554cp+4, lgamma(@as(f80, -0x1.fffffffffffffp+0)));
    try std.testing.expectEqual(0x2.af9a1ce04d03f778p+4, lgamma(@as(f80, -0x1.fffffffffffffffep+0)));
    try std.testing.expectEqual(0xe.8e5bf3a347bbb1ep+0, lgamma(@as(f80, -0x2.000004p+0)));
    try std.testing.expectEqual(0x2.2a848ae66fa85a6p+4, lgamma(@as(f80, -0x2.0000000000002p+0)));
    try std.testing.expectEqual(0x2.a482fb60cfe6ffep+4, lgamma(@as(f80, -0x2.0000000000000004p+0)));
    try std.testing.expectEqual(0xd.751d54afa9a2256p+0, lgamma(@as(f80, -0x2.fffffcp+0)));
    try std.testing.expectEqual(0x2.18f0a06bc2a55424p+4, lgamma(@as(f80, -0x2.ffffffffffffep+0)));
    try std.testing.expectEqual(0x2.92ef10e622e3f548p+4, lgamma(@as(f80, -0x2.fffffffffffffffcp+0)));
    try std.testing.expectEqual(0xd.751d4aa3223696ap+0, lgamma(@as(f80, -0x3.000004p+0)));
    try std.testing.expectEqual(0x2.18f0a06bc2a54f2p+4, lgamma(@as(f80, -0x3.0000000000002p+0)));
    try std.testing.expectEqual(0x2.92ef10e622e3f548p+4, lgamma(@as(f80, -0x3.0000000000000004p+0)));
    try std.testing.expectEqual(0xc.123925c00603b21p+0, lgamma(@as(f80, -0x3.fffffcp+0)));
    try std.testing.expectEqual(0x2.02c25d6cc86b657p+4, lgamma(@as(f80, -0x3.ffffffffffffep+0)));
    try std.testing.expectEqual(0x2.7cc0cde728aa0614p+4, lgamma(@as(f80, -0x3.fffffffffffffffcp+0)));
    try std.testing.expectEqual(0xb.60c6fbb5695c876p+0, lgamma(@as(f80, -0x4.000008p+0)));
    try std.testing.expectEqual(0x1.f7ab3bed4b4e64cap+4, lgamma(@as(f80, -0x4.0000000000004p+0)));
    try std.testing.expectEqual(0x2.71a9ac67ab8d0e78p+4, lgamma(@as(f80, -0x4.0000000000000008p+0)));
    try std.testing.expectEqual(0x9.c4c2f5e938fb4f8p+0, lgamma(@as(f80, -0x4.fffff8p+0)));
    try std.testing.expectEqual(0x1.ddeaf9f55dc13e3ap+4, lgamma(@as(f80, -0x4.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x2.57e96a6fbdffdb0cp+4, lgamma(@as(f80, -0x4.fffffffffffffff8p+0)));
    try std.testing.expectEqual(0x9.c4c2da9cf6f0fedp+0, lgamma(@as(f80, -0x5.000008p+0)));
    try std.testing.expectEqual(0x1.ddeaf9f55dc13094p+4, lgamma(@as(f80, -0x5.0000000000004p+0)));
    try std.testing.expectEqual(0x2.57e96a6fbdffdb0cp+4, lgamma(@as(f80, -0x5.0000000000000008p+0)));
    try std.testing.expectEqual(0x7.fa12379bec516538p+0, lgamma(@as(f80, -0x5.fffff8p+0)));
    try std.testing.expectEqual(0x1.c13fedfb33a13cb2p+4, lgamma(@as(f80, -0x5.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x2.3b3e5e7593dfd8dcp+4, lgamma(@as(f80, -0x5.fffffffffffffff8p+0)));
    try std.testing.expectEqual(0x7.fa1219a4ff9c69ep+0, lgamma(@as(f80, -0x6.000008p+0)));
    try std.testing.expectEqual(0x1.c13fedfb33a12db6p+4, lgamma(@as(f80, -0x6.0000000000004p+0)));
    try std.testing.expectEqual(0x2.3b3e5e7593dfd8d8p+4, lgamma(@as(f80, -0x6.0000000000000008p+0)));
    try std.testing.expectEqual(0x6.07eb0ddd58f5bbbp+0, lgamma(@as(f80, -0x6.fffff8p+0)));
    try std.testing.expectEqual(0x1.a21d7b4d0146e5fp+4, lgamma(@as(f80, -0x6.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x2.1c1bebc761858188p+4, lgamma(@as(f80, -0x6.fffffffffffffff8p+0)));
    try std.testing.expectEqual(0x6.07eaed9d47ae7738p+0, lgamma(@as(f80, -0x7.000008p+0)));
    try std.testing.expectEqual(0x1.a21d7b4d0146d5dp+4, lgamma(@as(f80, -0x7.0000000000004p+0)));
    try std.testing.expectEqual(0x2.1c1bebc761858184p+4, lgamma(@as(f80, -0x7.0000000000000008p+0)));
    try std.testing.expectEqual(0x3.f394c6f5e387cebp+0, lgamma(@as(f80, -0x7.fffff8p+0)));
    try std.testing.expectEqual(0x1.80d816ce89efffap+4, lgamma(@as(f80, -0x7.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x1.fad68748ea2e9ab6p+4, lgamma(@as(f80, -0x7.fffffffffffffff8p+0)));
    try std.testing.expectEqual(0x3.42227b9df8fdfa1cp+0, lgamma(@as(f80, -0x8.00001p+0)));
    try std.testing.expectEqual(0x1.75c0f54f0cd2ee54p+4, lgamma(@as(f80, -0x8.0000000000008p+0)));
    try std.testing.expectEqual(0x1.efbf65c96d11a318p+4, lgamma(@as(f80, -0x8.000000000000001p+0)));
    try std.testing.expectEqual(0x1.0fa5728f979e8bdp+0, lgamma(@as(f80, -0x8.fffffp+0)));
    try std.testing.expectEqual(0x1.52992059b2ccfc4ap+4, lgamma(@as(f80, -0x8.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.cc9790d4130b8deep+4, lgamma(@as(f80, -0x8.fffffffffffffffp+0)));
    try std.testing.expectEqual(0x1.0fa52a813c2c749ep+0, lgamma(@as(f80, -0x9.00001p+0)));
    try std.testing.expectEqual(0x1.52992059b2ccd842p+4, lgamma(@as(f80, -0x9.0000000000008p+0)));
    try std.testing.expectEqual(0x1.cc9790d4130b8deap+4, lgamma(@as(f80, -0x9.000000000000001p+0)));
    try std.testing.expectEqual(-0x1.3dd0c34d79694344p+0, lgamma(@as(f80, -0x9.fffffp+0)));
    try std.testing.expectEqual(0x1.2dc1bce24822d21p+4, lgamma(@as(f80, -0x9.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.a7c02d5ca86162e8p+4, lgamma(@as(f80, -0x9.fffffffffffffffp+0)));
    try std.testing.expectEqual(-0x1.3dd10e8f080e8daap+0, lgamma(@as(f80, -0xa.00001p+0)));
    try std.testing.expectEqual(0x1.2dc1bce24822ac7p+4, lgamma(@as(f80, -0xa.0000000000008p+0)));
    try std.testing.expectEqual(0x1.a7c02d5ca86162e4p+4, lgamma(@as(f80, -0xa.000000000000001p+0)));
    try std.testing.expectEqual(-0x3.a3ad38c9033a659cp+0, lgamma(@as(f80, -0xa.fffffp+0)));
    try std.testing.expectEqual(0x1.0763f57349b43b5cp+4, lgamma(@as(f80, -0xa.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.816265eda9f2cb7ap+4, lgamma(@as(f80, -0xa.fffffffffffffffp+0)));
    try std.testing.expectEqual(-0x3.a3ad86f34c0e3ba4p+0, lgamma(@as(f80, -0xb.00001p+0)));
    try std.testing.expectEqual(0x1.0763f57349b41446p+4, lgamma(@as(f80, -0xb.0000000000008p+0)));
    try std.testing.expectEqual(0x1.816265eda9f2cb74p+4, lgamma(@as(f80, -0xb.000000000000001p+0)));
    try std.testing.expectEqual(-0x6.1fd00f0e21b3c988p+0, lgamma(@as(f80, -0xb.fffffp+0)));
    try std.testing.expectEqual(0xd.fa1c7f9a2774239p+0, lgamma(@as(f80, -0xb.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.59a0387402b5d1acp+4, lgamma(@as(f80, -0xb.fffffffffffffffp+0)));
    try std.testing.expectEqual(-0x6.1fd05fe315324a38p+0, lgamma(@as(f80, -0xc.00001p+0)));
    try std.testing.expectEqual(0xd.fa1c7f9a27719cfp+0, lgamma(@as(f80, -0xc.0000000000008p+0)));
    try std.testing.expectEqual(0x1.59a0387402b5d1a8p+4, lgamma(@as(f80, -0xc.000000000000001p+0)));
    try std.testing.expectEqual(-0x8.b07093393f8bec6p+0, lgamma(@as(f80, -0xc.fffffp+0)));
    try std.testing.expectEqual(0xb.697bfa33f5ea0d9p+0, lgamma(@as(f80, -0xc.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.3096301d9f9d2fbp+4, lgamma(@as(f80, -0xc.fffffffffffffffp+0)));
    try std.testing.expectEqual(-0x8.b070e6845a6ce34p+0, lgamma(@as(f80, -0xd.00001p+0)));
    try std.testing.expectEqual(0xb.697bfa33f5e7734p+0, lgamma(@as(f80, -0xd.0000000000008p+0)));
    try std.testing.expectEqual(0x1.3096301d9f9d2faap+4, lgamma(@as(f80, -0xd.000000000000001p+0)));
    try std.testing.expectEqual(-0xb.5409d4efa4b70f9p+0, lgamma(@as(f80, -0xd.fffffp+0)));
    try std.testing.expectEqual(0x8.c5e2b758fe7527dp+0, lgamma(@as(f80, -0xd.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.065c9beff025e0cp+4, lgamma(@as(f80, -0xd.fffffffffffffffp+0)));
    try std.testing.expectEqual(-0xb.540a2a83e42a4f9p+0, lgamma(@as(f80, -0xe.00001p+0)));
    try std.testing.expectEqual(0x8.c5e2b758fe727b2p+0, lgamma(@as(f80, -0xe.0000000000008p+0)));
    try std.testing.expectEqual(0x1.065c9beff025e0bap+4, lgamma(@as(f80, -0xe.000000000000001p+0)));
    try std.testing.expectEqual(-0xe.094c9b083ca94dp+0, lgamma(@as(f80, -0xe.fffffp+0)));
    try std.testing.expectEqual(0x6.109ff02f55715028p+0, lgamma(@as(f80, -0xe.ffffffffffff8p+0)));
    try std.testing.expectEqual(0xd.b086f7d5595a2bep+0, lgamma(@as(f80, -0xe.fffffffffffffffp+0)));
    try std.testing.expectEqual(-0xe.094cf2be9e3eaf2p+0, lgamma(@as(f80, -0xf.00001p+0)));
    try std.testing.expectEqual(0x6.109ff02f556e9278p+0, lgamma(@as(f80, -0xf.0000000000008p+0)));
    try std.testing.expectEqual(0xd.b086f7d5595a2b9p+0, lgamma(@as(f80, -0xf.000000000000001p+0)));
    try std.testing.expectEqual(-0x1.0cf14f9e783e6b3cp+4, lgamma(@as(f80, -0xf.fffffp+0)));
    try std.testing.expectEqual(0x3.4ad790500e33717cp+0, lgamma(@as(f80, -0xf.ffffffffffff8p+0)));
    try std.testing.expectEqual(0xa.eabe97f6121c453p+0, lgamma(@as(f80, -0xf.fffffffffffffffp+0)));
    try std.testing.expectEqual(-0x1.180879870e33e356p+4, lgamma(@as(f80, -0x1.000002p+4)));
    try std.testing.expectEqual(0x2.996578583c5fc344p+0, lgamma(@as(f80, -0x1.0000000000001p+4)));
    try std.testing.expectEqual(0xa.394c7ffe404ccbp+0, lgamma(@as(f80, -0x1.0000000000000002p+4)));
    try std.testing.expectEqual(-0x1.455d45b618e1f038p+4, lgamma(@as(f80, -0x1.0ffffep+4)));
    try std.testing.expectEqual(-0x3.be7ffe71389cc268p-4, lgamma(@as(f80, -0x1.0ffffffffffffp+4)));
    try std.testing.expectEqual(0x7.63ff07bef05d91d8p+0, lgamma(@as(f80, -0x1.0ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.455d51292150d8bap+4, lgamma(@as(f80, -0x1.100002p+4)));
    try std.testing.expectEqual(-0x3.be7ffe7138f85aacp-4, lgamma(@as(f80, -0x1.1000000000001p+4)));
    try std.testing.expectEqual(0x7.63ff07bef05d912p+0, lgamma(@as(f80, -0x1.1000000000000002p+4)));
    try std.testing.expectEqual(-0x1.739c3c0e7e3dc748p+4, lgamma(@as(f80, -0x1.1ffffep+4)));
    try std.testing.expectEqual(-0x3.1fd7673485ba8a88p+0, lgamma(@as(f80, -0x1.1ffffffffffffp+4)));
    try std.testing.expectEqual(0x4.800fa0717e2cc54p+0, lgamma(@as(f80, -0x1.1ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.739c47ba6a3ae8acp+4, lgamma(@as(f80, -0x1.200002p+4)));
    try std.testing.expectEqual(-0x3.1fd7673485c0607cp+0, lgamma(@as(f80, -0x1.2000000000001p+4)));
    try std.testing.expectEqual(0x4.800fa0717e2cc48p+0, lgamma(@as(f80, -0x1.2000000000000002p+4)));
    try std.testing.expectEqual(-0x1.a2b8a7ff951d4cd8p+4, lgamma(@as(f80, -0x1.2ffffep+4)));
    try std.testing.expectEqual(-0x6.119e27f51c200b5p+0, lgamma(@as(f80, -0x1.2ffffffffffffp+4)));
    try std.testing.expectEqual(0x1.8e48dfb0e7c736fep+0, lgamma(@as(f80, -0x1.2ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.a2b8b3e16627e78p+4, lgamma(@as(f80, -0x1.300002p+4)));
    try std.testing.expectEqual(-0x6.119e27f51c25fc38p+0, lgamma(@as(f80, -0x1.3000000000001p+4)));
    try std.testing.expectEqual(0x1.8e48dfb0e7c7364p+0, lgamma(@as(f80, -0x1.3000000000000002p+4)));
    try std.testing.expectEqual(-0x1.d2a72cdce34ac164p+4, lgamma(@as(f80, -0x1.3ffffep+4)));
    try std.testing.expectEqual(-0x9.108677639892289p+0, lgamma(@as(f80, -0x1.3ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.709f6fbd94aaf306p+0, lgamma(@as(f80, -0x1.3ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.d2a738f1e7888f4p+4, lgamma(@as(f80, -0x1.400002p+4)));
    try std.testing.expectEqual(-0x9.108677639898331p+0, lgamma(@as(f80, -0x1.4000000000001p+4)));
    try std.testing.expectEqual(-0x1.709f6fbd94aaf3c8p+0, lgamma(@as(f80, -0x1.4000000000000002p+4)));
    try std.testing.expectEqual(-0x2.035d89ed6121f85cp+4, lgamma(@as(f80, -0x1.4ffffep+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e683b1p+0, lgamma(@as(f80, -0x1.4ffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.7c05424b8a8111cp+0, lgamma(@as(f80, -0x1.4ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.035d9633286bf6f8p+4, lgamma(@as(f80, -0x1.500002p+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e6e5ep+0, lgamma(@as(f80, -0x1.5000000000001p+4)));
    try std.testing.expectEqual(-0x4.7c05424b8a811288p+0, lgamma(@as(f80, -0x1.5000000000000002p+4)));
    try std.testing.expectEqual(-0x2.34d272c496dc021cp+4, lgamma(@as(f80, -0x1.5ffffep+4)));
    try std.testing.expectEqual(-0xf.333ad8d94721201p+0, lgamma(@as(f80, -0x1.5ffffffffffffp+4)));
    try std.testing.expectEqual(-0x7.9353d133433a0268p+0, lgamma(@as(f80, -0x1.5ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.34d27f38e9c8e974p+4, lgamma(@as(f80, -0x1.600002p+4)));
    try std.testing.expectEqual(-0xf.333ad8d947275a4p+0, lgamma(@as(f80, -0x1.6000000000001p+4)));
    try std.testing.expectEqual(-0x7.9353d133433a033p+0, lgamma(@as(f80, -0x1.6000000000000002p+4)));
    try std.testing.expectEqual(-0x2.66fd6ea9f77b79a8p+4, lgamma(@as(f80, -0x1.6ffffep+4)));
    try std.testing.expectEqual(-0x1.255ea98937d9f162p+4, lgamma(@as(f80, -0x1.6ffffffffffffp+4)));
    try std.testing.expectEqual(-0xa.b60390ed79b8038p+0, lgamma(@as(f80, -0x1.6ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.66fd7b4acff91314p+4, lgamma(@as(f80, -0x1.700002p+4)));
    try std.testing.expectEqual(-0x1.255ea98937da5668p+4, lgamma(@as(f80, -0x1.7000000000001p+4)));
    try std.testing.expectEqual(-0xa.b60390ed79b8045p+0, lgamma(@as(f80, -0x1.7000000000000002p+4)));
    try std.testing.expectEqual(-0x2.99d6bd8dc680078p+4, lgamma(@as(f80, -0x1.7ffffep+4)));
    try std.testing.expectEqual(-0x1.5837f8825c33e21ep+4, lgamma(@as(f80, -0x1.7ffffffffffffp+4)));
    try std.testing.expectEqual(-0xd.e398807fbf571ap+0, lgamma(@as(f80, -0x1.7ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.99d6ca5949a84b98p+4, lgamma(@as(f80, -0x1.800002p+4)));
    try std.testing.expectEqual(-0x1.5837f8825c34487ap+4, lgamma(@as(f80, -0x1.8000000000001p+4)));
    try std.testing.expectEqual(-0xd.e398807fbf571adp+0, lgamma(@as(f80, -0x1.8000000000000002p+4)));
    try std.testing.expectEqual(-0x2.cd57416926b9198cp+4, lgamma(@as(f80, -0x1.8ffffep+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374e485p+4, lgamma(@as(f80, -0x1.8ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.11ba0bf7d70fd876p+4, lgamma(@as(f80, -0x1.8ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.cd574e5d9fa3edp+4, lgamma(@as(f80, -0x1.900002p+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374eaff4p+4, lgamma(@as(f80, -0x1.9000000000001p+4)));
    try std.testing.expectEqual(-0x1.11ba0bf7d70fd882p+4, lgamma(@as(f80, -0x1.9000000000000002p+4)));
    try std.testing.expectEqual(-0x3.01786b2b55b39354p+4, lgamma(@as(f80, -0x1.9ffffep+4)));
    try std.testing.expectEqual(-0x1.bfd9a6481783e14ap+4, lgamma(@as(f80, -0x1.9ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.45db35cdb745720ep+4, lgamma(@as(f80, -0x1.9ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.0178784731148e2cp+4, lgamma(@as(f80, -0x1.a00002p+4)));
    try std.testing.expectEqual(-0x1.bfd9a64817844a2ap+4, lgamma(@as(f80, -0x1.a000000000001p+4)));
    try std.testing.expectEqual(-0x1.45db35cdb745721ap+4, lgamma(@as(f80, -0x1.a000000000000002p+4)));
    try std.testing.expectEqual(-0x3.36342a886637ea3cp+4, lgamma(@as(f80, -0x1.affffep+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8d007ap+4, lgamma(@as(f80, -0x1.affffffffffffp+4)));
    try std.testing.expectEqual(-0x1.7a96f53dbe4e91d4p+4, lgamma(@as(f80, -0x1.affffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.363437ca2ea26058p+4, lgamma(@as(f80, -0x1.b00002p+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8d6a88p+4, lgamma(@as(f80, -0x1.b000000000001p+4)));
    try std.testing.expectEqual(-0x1.7a96f53dbe4e91ep+4, lgamma(@as(f80, -0x1.b000000000000002p+4)));
    try std.testing.expectEqual(-0x3.6b84e02349a7940cp+4, lgamma(@as(f80, -0x1.bffffep+4)));
    try std.testing.expectEqual(-0x2.29e61b654b21467p+4, lgamma(@as(f80, -0x1.bffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.afe7aaeaeae2d85ep+4, lgamma(@as(f80, -0x1.bffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.6b84ed89a45b2eb8p+4, lgamma(@as(f80, -0x1.c00002p+4)));
    try std.testing.expectEqual(-0x2.29e61b654b21b1a4p+4, lgamma(@as(f80, -0x1.c000000000001p+4)));
    try std.testing.expectEqual(-0x1.afe7aaeaeae2d86ap+4, lgamma(@as(f80, -0x1.c000000000000002p+4)));
    try std.testing.expectEqual(-0x3.a16551a93dea66acp+4, lgamma(@as(f80, -0x1.cffffep+4)));
    try std.testing.expectEqual(-0x2.5fc68cfce71d836p+4, lgamma(@as(f80, -0x1.cffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.e5c81c8286df15d8p+4, lgamma(@as(f80, -0x1.cffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.a1655f32e810c39p+4, lgamma(@as(f80, -0x1.d00002p+4)));
    try std.testing.expectEqual(-0x2.5fc68cfce71defacp+4, lgamma(@as(f80, -0x1.d000000000001p+4)));
    try std.testing.expectEqual(-0x1.e5c81c8286df15e6p+4, lgamma(@as(f80, -0x1.d000000000000002p+4)));
    try std.testing.expectEqual(-0x3.d7d09f8a4486822p+4, lgamma(@as(f80, -0x1.dffffep+4)));
    try std.testing.expectEqual(-0x2.9631daeefecab874p+4, lgamma(@as(f80, -0x1.dffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.1c336a749e8c4b74p+4, lgamma(@as(f80, -0x1.dffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.d7d0ad3610cf0124p+4, lgamma(@as(f80, -0x1.e00002p+4)));
    try std.testing.expectEqual(-0x2.9631daeefecb25dp+4, lgamma(@as(f80, -0x1.e000000000001p+4)));
    try std.testing.expectEqual(-0x2.1c336a749e8c4b84p+4, lgamma(@as(f80, -0x1.e000000000000002p+4)));
    try std.testing.expectEqual(-0x4.0ec23c0ae2bc2538p+4, lgamma(@as(f80, -0x1.effffep+4)));
    try std.testing.expectEqual(-0x2.cd23778021216bdp+4, lgamma(@as(f80, -0x1.effffffffffffp+4)));
    try std.testing.expectEqual(-0x2.53250705c0e2ff58p+4, lgamma(@as(f80, -0x1.effffffffffffffep+4)));
    try std.testing.expectEqual(-0x4.0ec249d7b746b4cp+4, lgamma(@as(f80, -0x1.f00002p+4)));
    try std.testing.expectEqual(-0x2.cd2377802121da38p+4, lgamma(@as(f80, -0x1.f000000000001p+4)));
    try std.testing.expectEqual(-0x2.53250705c0e2ff64p+4, lgamma(@as(f80, -0x1.f000000000000002p+4)));
    try std.testing.expectEqual(-0x4.4635e378544cf34p+4, lgamma(@as(f80, -0x1.fffffep+4)));
    try std.testing.expectEqual(-0x3.04971efd92b24158p+4, lgamma(@as(f80, -0x1.fffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.8a98ae833273d55cp+4, lgamma(@as(f80, -0x1.fffffffffffffffep+4)));
    try std.testing.expectEqual(-0x4.514d19db0f00e278p+4, lgamma(@as(f80, -0x2.000004p+4)));
    try std.testing.expectEqual(-0x3.0fae407d0fcfe00cp+4, lgamma(@as(f80, -0x2.0000000000002p+4)));
    try std.testing.expectEqual(-0x2.95afd002af90cd0cp+4, lgamma(@as(f80, -0x2.0000000000000004p+4)));
    try std.testing.expectEqual(-0x4.893eafcc099b56ep+4, lgamma(@as(f80, -0x2.0ffffcp+4)));
    try std.testing.expectEqual(-0x3.479ff266bb40a24cp+4, lgamma(@as(f80, -0x2.0fffffffffffep+4)));
    try std.testing.expectEqual(-0x2.cda181ec5b026ef8p+4, lgamma(@as(f80, -0x2.0ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x4.893ecbe3c23456ep+4, lgamma(@as(f80, -0x2.100004p+4)));
    try std.testing.expectEqual(-0x3.479ff266bb41830cp+4, lgamma(@as(f80, -0x2.1000000000002p+4)));
    try std.testing.expectEqual(-0x2.cda181ec5b026f14p+4, lgamma(@as(f80, -0x2.1000000000000004p+4)));
    try std.testing.expectEqual(-0x4.c1aaa8b15d9907ap+4, lgamma(@as(f80, -0x2.1ffffcp+4)));
    try std.testing.expectEqual(-0x3.800beb6a2d5c8c94p+4, lgamma(@as(f80, -0x2.1fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.060d7aefcd1e5a3p+4, lgamma(@as(f80, -0x2.1ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x4.c1aac505526e43ep+4, lgamma(@as(f80, -0x2.200004p+4)));
    try std.testing.expectEqual(-0x3.800beb6a2d5d6f34p+4, lgamma(@as(f80, -0x2.2000000000002p+4)));
    try std.testing.expectEqual(-0x3.060d7aefcd1e5a4cp+4, lgamma(@as(f80, -0x2.2000000000000004p+4)));
    try std.testing.expectEqual(-0x4.fa8d5d3a3bac5a6p+4, lgamma(@as(f80, -0x2.2ffffcp+4)));
    try std.testing.expectEqual(-0x3.b8eea0104d44166cp+4, lgamma(@as(f80, -0x2.2fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.3ef02f95ed05e4fp+4, lgamma(@as(f80, -0x2.2ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x4.fa8d79c8b429d12p+4, lgamma(@as(f80, -0x2.300004p+4)));
    try std.testing.expectEqual(-0x3.b8eea0104d44fadcp+4, lgamma(@as(f80, -0x2.3000000000002p+4)));
    try std.testing.expectEqual(-0x3.3ef02f95ed05e50cp+4, lgamma(@as(f80, -0x2.3000000000000004p+4)));
    try std.testing.expectEqual(-0x5.33e375121e252908p+4, lgamma(@as(f80, -0x2.3ffffcp+4)));
    try std.testing.expectEqual(-0x3.f244b804a18419ecp+4, lgamma(@as(f80, -0x2.3fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.7846478a4145e954p+4, lgamma(@as(f80, -0x2.3ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x5.33e391d97a30d8bp+4, lgamma(@as(f80, -0x2.400004p+4)));
    try std.testing.expectEqual(-0x3.f244b804a1850024p+4, lgamma(@as(f80, -0x2.4000000000002p+4)));
    try std.testing.expectEqual(-0x3.7846478a4145e97p+4, lgamma(@as(f80, -0x2.4000000000000004p+4)));
    try std.testing.expectEqual(-0x5.6da9c6d2e6bb76c8p+4, lgamma(@as(f80, -0x2.4ffffcp+4)));
    try std.testing.expectEqual(-0x4.2c0b09e11713938p+4, lgamma(@as(f80, -0x2.4fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.b20c9966b6d563c4p+4, lgamma(@as(f80, -0x2.4ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x5.6da9e3d19cb94ffp+4, lgamma(@as(f80, -0x2.500004p+4)));
    try std.testing.expectEqual(-0x4.2c0b09e117147b7p+4, lgamma(@as(f80, -0x2.5000000000002p+4)));
    try std.testing.expectEqual(-0x3.b20c9966b6d563ep+4, lgamma(@as(f80, -0x2.5000000000000004p+4)));
    try std.testing.expectEqual(-0x5.a7dd54437ab7f3fp+4, lgamma(@as(f80, -0x2.5ffffcp+4)));
    try std.testing.expectEqual(-0x4.663e976c9d96e32p+4, lgamma(@as(f80, -0x2.5fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.ec4026f23d58b44p+4, lgamma(@as(f80, -0x2.5ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x5.a7dd717815c3466p+4, lgamma(@as(f80, -0x2.600004p+4)));
    try std.testing.expectEqual(-0x4.663e976c9d97ccc8p+4, lgamma(@as(f80, -0x2.6000000000002p+4)));
    try std.testing.expectEqual(-0x3.ec4026f23d58b46p+4, lgamma(@as(f80, -0x2.6000000000000004p+4)));
    try std.testing.expectEqual(-0x5.e27b46fa492f70b8p+4, lgamma(@as(f80, -0x2.6ffffcp+4)));
    try std.testing.expectEqual(-0x4.a0dc8a3dadb28ee8p+4, lgamma(@as(f80, -0x2.6fffffffffffep+4)));
    try std.testing.expectEqual(-0x4.26de19c34d7460d8p+4, lgamma(@as(f80, -0x2.6ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x5.e27b64636782f7bp+4, lgamma(@as(f80, -0x2.700004p+4)));
    try std.testing.expectEqual(-0x4.a0dc8a3dadb37a3p+4, lgamma(@as(f80, -0x2.7000000000002p+4)));
    try std.testing.expectEqual(-0x4.26de19c34d7460fp+4, lgamma(@as(f80, -0x2.7000000000000004p+4)));
    try std.testing.expectEqual(-0x6.1d80ed571479dcep+4, lgamma(@as(f80, -0x2.7ffffcp+4)));
    try std.testing.expectEqual(-0x4.dbe230b41296a858p+4, lgamma(@as(f80, -0x2.7fffffffffffep+4)));
    try std.testing.expectEqual(-0x4.61e3c039b2587b1p+4, lgamma(@as(f80, -0x2.7ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x6.1d810af366009708p+4, lgamma(@as(f80, -0x2.800004p+4)));
    try std.testing.expectEqual(-0x4.dbe230b412979538p+4, lgamma(@as(f80, -0x2.8000000000002p+4)));
    try std.testing.expectEqual(-0x4.61e3c039b2587b3p+4, lgamma(@as(f80, -0x2.8000000000000004p+4)));
    try std.testing.expectEqual(-0x6.58ebb7c93810d52p+4, lgamma(@as(f80, -0x2.8ffffcp+4)));
    try std.testing.expectEqual(-0x5.174cfb3f2fef42e8p+4, lgamma(@as(f80, -0x2.8fffffffffffep+4)));
    try std.testing.expectEqual(-0x4.9d4e8ac4cfb11668p+4, lgamma(@as(f80, -0x2.8ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x6.58ebd5977d1aae78p+4, lgamma(@as(f80, -0x2.900004p+4)));
    try std.testing.expectEqual(-0x5.174cfb3f2ff03158p+4, lgamma(@as(f80, -0x2.9000000000002p+4)));
    try std.testing.expectEqual(-0x4.9d4e8ac4cfb11688p+4, lgamma(@as(f80, -0x2.9000000000000004p+4)));
    try std.testing.expectEqual(-0x6.94b93659330503bp+4, lgamma(@as(f80, -0x2.9ffffcp+4)));
    try std.testing.expectEqual(-0x5.531a79e78c699ba8p+4, lgamma(@as(f80, -0x2.9fffffffffffep+4)));
    try std.testing.expectEqual(-0x4.d91c096d2c2b6ffp+4, lgamma(@as(f80, -0x2.9ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x6.94b954583b1b0ddp+4, lgamma(@as(f80, -0x2.a00004p+4)));
    try std.testing.expectEqual(-0x5.531a79e78c6a8bap+4, lgamma(@as(f80, -0x2.a000000000002p+4)));
    try std.testing.expectEqual(-0x4.d91c096d2c2b701p+4, lgamma(@as(f80, -0x2.a000000000000004p+4)));
    try std.testing.expectEqual(-0x6.d0e7166d8c7dd2a8p+4, lgamma(@as(f80, -0x2.affffcp+4)));
    try std.testing.expectEqual(-0x5.8f485a13b641bd18p+4, lgamma(@as(f80, -0x2.afffffffffffep+4)));
    try std.testing.expectEqual(-0x5.1549e99956039218p+4, lgamma(@as(f80, -0x2.affffffffffffffcp+4)));
    try std.testing.expectEqual(-0x6.d0e7349c35525fcp+4, lgamma(@as(f80, -0x2.b00004p+4)));
    try std.testing.expectEqual(-0x5.8f485a13b642ae88p+4, lgamma(@as(f80, -0x2.b000000000002p+4)));
    try std.testing.expectEqual(-0x5.1549e99956039238p+4, lgamma(@as(f80, -0x2.b000000000000004p+4)));
    try std.testing.expectEqual(-0x7.0d7320c43f54d4p+4, lgamma(@as(f80, -0x2.bffffcp+4)));
    try std.testing.expectEqual(-0x5.cbd46481aeea43p+4, lgamma(@as(f80, -0x2.bfffffffffffep+4)));
    try std.testing.expectEqual(-0x5.51d5f4074eac18cp+4, lgamma(@as(f80, -0x2.bffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x7.0d733f2173cc49d8p+4, lgamma(@as(f80, -0x2.c00004p+4)));
    try std.testing.expectEqual(-0x5.cbd46481aeeb35e8p+4, lgamma(@as(f80, -0x2.c000000000002p+4)));
    try std.testing.expectEqual(-0x5.51d5f4074eac18ep+4, lgamma(@as(f80, -0x2.c000000000000004p+4)));
    try std.testing.expectEqual(-0x7.4a5b379ac57bf5a8p+4, lgamma(@as(f80, -0x2.cffffcp+4)));
    try std.testing.expectEqual(-0x6.08bc7b6ef67d8ae8p+4, lgamma(@as(f80, -0x2.cfffffffffffep+4)));
    try std.testing.expectEqual(-0x5.8ebe0af4963f6158p+4, lgamma(@as(f80, -0x2.cffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x7.4a5b56257ccb99p+4, lgamma(@as(f80, -0x2.d00004p+4)));
    try std.testing.expectEqual(-0x6.08bc7b6ef67e7f38p+4, lgamma(@as(f80, -0x2.d000000000002p+4)));
    try std.testing.expectEqual(-0x5.8ebe0af4963f6178p+4, lgamma(@as(f80, -0x2.d000000000000004p+4)));
    try std.testing.expectEqual(-0x7.879d54ffa33864dp+4, lgamma(@as(f80, -0x2.dffffcp+4)));
    try std.testing.expectEqual(-0x6.45fe98ea170261ep+4, lgamma(@as(f80, -0x2.dfffffffffffep+4)));
    try std.testing.expectEqual(-0x5.cc00286fb6c43908p+4, lgamma(@as(f80, -0x2.dffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x7.879d73b6e018ba4p+4, lgamma(@as(f80, -0x2.e00004p+4)));
    try std.testing.expectEqual(-0x6.45fe98ea17035798p+4, lgamma(@as(f80, -0x2.e000000000002p+4)));
    try std.testing.expectEqual(-0x5.cc00286fb6c43928p+4, lgamma(@as(f80, -0x2.e000000000000004p+4)));
    try std.testing.expectEqual(-0x7.c5378948fb919718p+4, lgamma(@as(f80, -0x2.effffcp+4)));
    try std.testing.expectEqual(-0x6.8398cd4938e3cde8p+4, lgamma(@as(f80, -0x2.efffffffffffep+4)));
    try std.testing.expectEqual(-0x6.099a5cced8a5a5b8p+4, lgamma(@as(f80, -0x2.effffffffffffffcp+4)));
    try std.testing.expectEqual(-0x7.c537a82bcb8243bp+4, lgamma(@as(f80, -0x2.f00004p+4)));
    try std.testing.expectEqual(-0x6.8398cd4938e4c4f8p+4, lgamma(@as(f80, -0x2.f000000000002p+4)));
    try std.testing.expectEqual(-0x6.099a5cced8a5a5d8p+4, lgamma(@as(f80, -0x2.f000000000000004p+4)));
    try std.testing.expectEqual(-0x8.0327f9ac47b31c9p+4, lgamma(@as(f80, -0x2.fffffcp+4)));
    try std.testing.expectEqual(-0x6.c1893dc1da5ab638p+4, lgamma(@as(f80, -0x2.ffffffffffffep+4)));
    try std.testing.expectEqual(-0x6.478acd477a1c8ecp+4, lgamma(@as(f80, -0x2.fffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x8.032818b9c24e73dp+4, lgamma(@as(f80, -0x3.000004p+4)));
    try std.testing.expectEqual(-0x6.c1893dc1da5baea8p+4, lgamma(@as(f80, -0x3.0000000000002p+4)));
    try std.testing.expectEqual(-0x6.478acd477a1c8ed8p+4, lgamma(@as(f80, -0x3.0000000000000004p+4)));
    try std.testing.expectEqual(-0x8.416cdef3c687166p+4, lgamma(@as(f80, -0x3.0ffffcp+4)));
    try std.testing.expectEqual(-0x6.ffce231e3f0f644p+4, lgamma(@as(f80, -0x3.0fffffffffffep+4)));
    try std.testing.expectEqual(-0x6.85cfb2a3ded13d68p+4, lgamma(@as(f80, -0x3.0ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x8.416cfe2b0ce3bcp+4, lgamma(@as(f80, -0x3.100004p+4)));
    try std.testing.expectEqual(-0x6.ffce231e3f105df8p+4, lgamma(@as(f80, -0x3.1000000000002p+4)));
    try std.testing.expectEqual(-0x6.85cfb2a3ded13d88p+4, lgamma(@as(f80, -0x3.1000000000000004p+4)));
    try std.testing.expectEqual(-0x8.8004844ea3dd201p+4, lgamma(@as(f80, -0x3.1ffffcp+4)));
    try std.testing.expectEqual(-0x7.3e65c88d9746c208p+4, lgamma(@as(f80, -0x3.1fffffffffffep+4)));
    try std.testing.expectEqual(-0x6.c467581337089bd8p+4, lgamma(@as(f80, -0x3.1ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x8.8004a3aedffc55p+4, lgamma(@as(f80, -0x3.200004p+4)));
    try std.testing.expectEqual(-0x7.3e65c88d9747bd1p+4, lgamma(@as(f80, -0x3.2000000000002p+4)));
    try std.testing.expectEqual(-0x6.c467581337089bf8p+4, lgamma(@as(f80, -0x3.2000000000000004p+4)));
    try std.testing.expectEqual(-0x8.beed463931cafd9p+4, lgamma(@as(f80, -0x3.2ffffcp+4)));
    try std.testing.expectEqual(-0x7.7d4e8a8c3948bfap+4, lgamma(@as(f80, -0x3.2fffffffffffep+4)));
    try std.testing.expectEqual(-0x7.03501a11d90a9a08p+4, lgamma(@as(f80, -0x3.2ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x8.beed65c196125abp+4, lgamma(@as(f80, -0x3.300004p+4)));
    try std.testing.expectEqual(-0x7.7d4e8a8c3949bbep+4, lgamma(@as(f80, -0x3.3000000000002p+4)));
    try std.testing.expectEqual(-0x7.03501a11d90a9a28p+4, lgamma(@as(f80, -0x3.3000000000000004p+4)));
    try std.testing.expectEqual(-0x8.fe25917adde26efp+4, lgamma(@as(f80, -0x3.3ffffcp+4)));
    try std.testing.expectEqual(-0x7.bc86d5e1969b5038p+4, lgamma(@as(f80, -0x3.3fffffffffffep+4)));
    try std.testing.expectEqual(-0x7.42886567365d2b4p+4, lgamma(@as(f80, -0x3.3ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x8.fe25b12aa49ff38p+4, lgamma(@as(f80, -0x3.400004p+4)));
    try std.testing.expectEqual(-0x7.bc86d5e1969c4dbp+4, lgamma(@as(f80, -0x3.4000000000002p+4)));
    try std.testing.expectEqual(-0x7.42886567365d2b6p+4, lgamma(@as(f80, -0x3.4000000000000004p+4)));
    try std.testing.expectEqual(-0x9.3dabe237d03ebd8p+4, lgamma(@as(f80, -0x3.4ffffcp+4)));
    try std.testing.expectEqual(-0x7.fc0d26b1db14a5p+4, lgamma(@as(f80, -0x3.4fffffffffffep+4)));
    try std.testing.expectEqual(-0x7.820eb6377ad680a8p+4, lgamma(@as(f80, -0x3.4ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x9.3dac020e3b36386p+4, lgamma(@as(f80, -0x3.500004p+4)));
    try std.testing.expectEqual(-0x7.fc0d26b1db15a3b8p+4, lgamma(@as(f80, -0x3.5000000000002p+4)));
    try std.testing.expectEqual(-0x7.820eb6377ad680c8p+4, lgamma(@as(f80, -0x3.5000000000000004p+4)));
    try std.testing.expectEqual(-0x9.7d7ec3145de00c1p+4, lgamma(@as(f80, -0x3.5ffffcp+4)));
    try std.testing.expectEqual(-0x8.3be007a15f3abbdp+4, lgamma(@as(f80, -0x3.5fffffffffffep+4)));
    try std.testing.expectEqual(-0x7.c1e19726fefc9808p+4, lgamma(@as(f80, -0x3.5ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x9.7d7ee310b5e1023p+4, lgamma(@as(f80, -0x3.600004p+4)));
    try std.testing.expectEqual(-0x8.3be007a15f3bbbbp+4, lgamma(@as(f80, -0x3.6000000000002p+4)));
    try std.testing.expectEqual(-0x7.c1e19726fefc9828p+4, lgamma(@as(f80, -0x3.6000000000000004p+4)));
    try std.testing.expectEqual(-0x9.bd9ccc68ab9aa23p+4, lgamma(@as(f80, -0x3.6ffffcp+4)));
    try std.testing.expectEqual(-0x8.7bfe11084b36861p+4, lgamma(@as(f80, -0x3.6fffffffffffep+4)));
    try std.testing.expectEqual(-0x8.01ffa08deaf862ep+4, lgamma(@as(f80, -0x3.6ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x9.bd9cec8a401dec1p+4, lgamma(@as(f80, -0x3.700004p+4)));
    try std.testing.expectEqual(-0x8.7bfe11084b37872p+4, lgamma(@as(f80, -0x3.7000000000002p+4)));
    try std.testing.expectEqual(-0x8.01ffa08deaf863p+4, lgamma(@as(f80, -0x3.7000000000000004p+4)));
    try std.testing.expectEqual(-0x9.fe04a3830c27439p+4, lgamma(@as(f80, -0x3.7ffffcp+4)));
    try std.testing.expectEqual(-0x8.bc65e834f4e7c3ap+4, lgamma(@as(f80, -0x3.7fffffffffffep+4)));
    try std.testing.expectEqual(-0x8.426777ba94a9a1p+4, lgamma(@as(f80, -0x3.7ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x9.fe04c3c932f3b21p+4, lgamma(@as(f80, -0x3.800004p+4)));
    try std.testing.expectEqual(-0x8.bc65e834f4e8c5dp+4, lgamma(@as(f80, -0x3.8000000000002p+4)));
    try std.testing.expectEqual(-0x8.426777ba94a9a12p+4, lgamma(@as(f80, -0x3.8000000000000004p+4)));
    try std.testing.expectEqual(-0xa.3eb4f9f7cb8c1f4p+4, lgamma(@as(f80, -0x3.8ffffcp+4)));
    try std.testing.expectEqual(-0x8.fd163ebbab51269p+4, lgamma(@as(f80, -0x3.8fffffffffffep+4)));
    try std.testing.expectEqual(-0x8.8317ce414b13048p+4, lgamma(@as(f80, -0x3.8ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0xa.3eb51a61e061893p+4, lgamma(@as(f80, -0x3.900004p+4)));
    try std.testing.expectEqual(-0x8.fd163ebbab5229ep+4, lgamma(@as(f80, -0x3.9000000000002p+4)));
    try std.testing.expectEqual(-0x8.8317ce414b1304ap+4, lgamma(@as(f80, -0x3.9000000000000004p+4)));
    try std.testing.expectEqual(-0xa.7fac8cfd3cebe97p+4, lgamma(@as(f80, -0x3.9ffffcp+4)));
    try std.testing.expectEqual(-0x9.3e0dd1d2c46a5b1p+4, lgamma(@as(f80, -0x3.9fffffffffffep+4)));
    try std.testing.expectEqual(-0x8.c40f6158642c399p+4, lgamma(@as(f80, -0x3.9ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0xa.7facad8aa13415ap+4, lgamma(@as(f80, -0x3.a00004p+4)));
    try std.testing.expectEqual(-0x9.3e0dd1d2c46b5f8p+4, lgamma(@as(f80, -0x3.a000000000002p+4)));
    try std.testing.expectEqual(-0x8.c40f6158642c39bp+4, lgamma(@as(f80, -0x3.a000000000000004p+4)));
    try std.testing.expectEqual(-0xa.c0ea24d2fe73637p+4, lgamma(@as(f80, -0x3.affffcp+4)));
    try std.testing.expectEqual(-0x9.7f4b69b9e1103d6p+4, lgamma(@as(f80, -0x3.afffffffffffep+4)));
    try std.testing.expectEqual(-0x9.054cf93f80d21c7p+4, lgamma(@as(f80, -0x3.affffffffffffffcp+4)));
    try std.testing.expectEqual(-0xa.c0ea458318f84e9p+4, lgamma(@as(f80, -0x3.b00004p+4)));
    try std.testing.expectEqual(-0x9.7f4b69b9e11142fp+4, lgamma(@as(f80, -0x3.b000000000002p+4)));
    try std.testing.expectEqual(-0x9.054cf93f80d21c9p+4, lgamma(@as(f80, -0x3.b000000000000004p+4)));
    try std.testing.expectEqual(-0xb.026c9433822c768p+4, lgamma(@as(f80, -0x3.bffffcp+4)));
    try std.testing.expectEqual(-0x9.c0cdd92b75da6a1p+4, lgamma(@as(f80, -0x3.bfffffffffffep+4)));
    try std.testing.expectEqual(-0x9.46cf68b1159c49bp+4, lgamma(@as(f80, -0x3.bffffffffffffffcp+4)));
    try std.testing.expectEqual(-0xb.026cb505bed383cp+4, lgamma(@as(f80, -0x3.c00004p+4)));
    try std.testing.expectEqual(-0x9.c0cdd92b75db70bp+4, lgamma(@as(f80, -0x3.c000000000002p+4)));
    try std.testing.expectEqual(-0x9.46cf68b1159c49dp+4, lgamma(@as(f80, -0x3.c000000000000004p+4)));
    try std.testing.expectEqual(0x4.2b2b52b5464eed28p-24, lgamma(@as(f80, -0x2.74ff9p+0)));
    try std.testing.expectEqual(-0x1.e4cf2421a71b195p-24, lgamma(@as(f80, -0x2.74ff94p+0)));
    try std.testing.expectEqual(0x4.0c8edb47fa1b3508p-56, lgamma(@as(f80, -0x2.74ff92c01f0d8p+0)));
    try std.testing.expectEqual(-0x2.c7343f216ac923dp-52, lgamma(@as(f80, -0x2.74ff92c01f0dap+0)));
    try std.testing.expectEqual(0x5.f29bbbdec3d4a8d8p-64, lgamma(@as(f80, -0x2.74ff92c01f0d82a8p+0)));
    try std.testing.expectEqual(-0x1.d5e9dcd11030bbap-68, lgamma(@as(f80, -0x2.74ff92c01f0d82acp+0)));
    try std.testing.expectEqual(-0x2.6b416efc56fe3eb8p-24, lgamma(@as(f80, -0x2.bf682p+0)));
    try std.testing.expectEqual(0x5.3d0a33adaf4f5878p-24, lgamma(@as(f80, -0x2.bf6824p+0)));
    try std.testing.expectEqual(-0x3.0c498b9ac27bd8d4p-52, lgamma(@as(f80, -0x2.bf6821437b2p+0)));
    try std.testing.expectEqual(0xc.7dc2985d3b44557p-56, lgamma(@as(f80, -0x2.bf6821437b202p+0)));
    try std.testing.expectEqual(-0x3.088b212f3705dc7cp-64, lgamma(@as(f80, -0x2.bf6821437b201978p+0)));
    try std.testing.expectEqual(0x4.9fc04911f55d35cp-64, lgamma(@as(f80, -0x2.bf6821437b20197cp+0)));
    try std.testing.expectEqual(0x1.bd69b50d51b1488p-20, lgamma(@as(f80, -0x3.24c1b4p+0)));
    try std.testing.expectEqual(-0x3.4a0c544eeb21a028p-24, lgamma(@as(f80, -0x3.24c1b8p+0)));
    try std.testing.expectEqual(0x7.a58178eb9e987768p-52, lgamma(@as(f80, -0x3.24c1b793cb35ep+0)));
    try std.testing.expectEqual(-0x7.ead1b6ac3791da08p-52, lgamma(@as(f80, -0x3.24c1b793cb36p+0)));
    try std.testing.expectEqual(0x5.c9c4ac92390bb718p-64, lgamma(@as(f80, -0x3.24c1b793cb35efb8p+0)));
    try std.testing.expectEqual(-0x1.956e1b29d7349244p-60, lgamma(@as(f80, -0x3.24c1b793cb35efbcp+0)));
    try std.testing.expectEqual(-0x3.511bca412890969p-20, lgamma(@as(f80, -0x3.f48e28p+0)));
    try std.testing.expectEqual(0x1.dd4b54ca863c1a48p-20, lgamma(@as(f80, -0x3.f48e2cp+0)));
    try std.testing.expectEqual(-0x1.ddc0336980b584d2p-52, lgamma(@as(f80, -0x3.f48e2a8f85fcap+0)));
    try std.testing.expectEqual(0x2.7957af96f2c10fp-48, lgamma(@as(f80, -0x3.f48e2a8f85fccp+0)));
    try std.testing.expectEqual(-0x1.130ae5c4f54dbe92p-60, lgamma(@as(f80, -0x3.f48e2a8f85fca17p+0)));
    try std.testing.expectEqual(0x4.1b5c7fd62043e838p-60, lgamma(@as(f80, -0x3.f48e2a8f85fca174p+0)));
    try std.testing.expectEqual(0xa.3165c90424948cfp-20, lgamma(@as(f80, -0x4.0a1398p+0)));
    try std.testing.expectEqual(-0x3.33cb5626dc331eccp-20, lgamma(@as(f80, -0x4.0a13ap+0)));
    try std.testing.expectEqual(0x5.1a6a37819144766p-48, lgamma(@as(f80, -0x4.0a139e16656p+0)));
    try std.testing.expectEqual(-0x1.982d05a2f456b4bp-48, lgamma(@as(f80, -0x4.0a139e1665604p+0)));
    try std.testing.expectEqual(0x6.103eebf7b96ec358p-60, lgamma(@as(f80, -0x4.0a139e16656030cp+0)));
    try std.testing.expectEqual(-0x7.54ef8e5151b2567p-60, lgamma(@as(f80, -0x4.0a139e16656030c8p+0)));
    try std.testing.expectEqual(-0x3.02165b2aa6eef1f4p-16, lgamma(@as(f80, -0x4.fdd5d8p+0)));
    try std.testing.expectEqual(0xa.22e7861540c9fcep-20, lgamma(@as(f80, -0x4.fdd5ep+0)));
    try std.testing.expectEqual(-0x1.8280d0ba86860c98p-44, lgamma(@as(f80, -0x4.fdd5de9bbabfp+0)));
    try std.testing.expectEqual(0x4.fa3d33517a9ecdc8p-48, lgamma(@as(f80, -0x4.fdd5de9bbabf4p+0)));
    try std.testing.expectEqual(-0x5.efcf1ba2a53065f8p-60, lgamma(@as(f80, -0x4.fdd5de9bbabf351p+0)));
    try std.testing.expectEqual(0x3.454c56251230ebfp-56, lgamma(@as(f80, -0x4.fdd5de9bbabf3518p+0)));
    try std.testing.expectEqual(0x2.e258f12a679ed408p-16, lgamma(@as(f80, -0x5.021a9p+0)));
    try std.testing.expectEqual(-0xf.89066929e3b181p-20, lgamma(@as(f80, -0x5.021a98p+0)));
    try std.testing.expectEqual(0x1.867827fdc0e929bcp-48, lgamma(@as(f80, -0x5.021a95fc2db64p+0)));
    try std.testing.expectEqual(-0x1.d50b5e02beb77b12p-44, lgamma(@as(f80, -0x5.021a95fc2db68p+0)));
    try std.testing.expectEqual(0x1.1b82d6b2b33045c6p-56, lgamma(@as(f80, -0x5.021a95fc2db64328p+0)));
    try std.testing.expectEqual(-0x2.bf62ea52828ff32cp-56, lgamma(@as(f80, -0x5.021a95fc2db6433p+0)));
    try std.testing.expectEqual(-0xf.15ee1077e22d21cp-16, lgamma(@as(f80, -0x5.ffa4b8p+0)));
    try std.testing.expectEqual(0x7.4bb0ef1ad813da38p-16, lgamma(@as(f80, -0x5.ffa4cp+0)));
    try std.testing.expectEqual(-0x4.2c4d3e7ff051f43p-44, lgamma(@as(f80, -0x5.ffa4bd647d034p+0)));
    try std.testing.expectEqual(0x7.04ae139d3fb74038p-44, lgamma(@as(f80, -0x5.ffa4bd647d038p+0)));
    try std.testing.expectEqual(-0xe.d9cc85177f957fap-56, lgamma(@as(f80, -0x5.ffa4bd647d0357d8p+0)));
    try std.testing.expectEqual(0x7.882a1f22de7c7118p-56, lgamma(@as(f80, -0x5.ffa4bd647d0357ep+0)));
    try std.testing.expectEqual(0x3.e9df593e904f847cp-16, lgamma(@as(f80, -0x6.005ac8p+0)));
    try std.testing.expectEqual(-0x1.2b35eea26dc93cd2p-12, lgamma(@as(f80, -0x6.005adp+0)));
    try std.testing.expectEqual(0xa.7dd3bd697d2c7b9p-44, lgamma(@as(f80, -0x6.005ac9625f23p+0)));
    try std.testing.expectEqual(-0xd.11e91b3ff8f4b94p-48, lgamma(@as(f80, -0x6.005ac9625f234p+0)));
    try std.testing.expectEqual(0x1.5f103a1b00a48702p-56, lgamma(@as(f80, -0x6.005ac9625f233b6p+0)));
    try std.testing.expectEqual(-0x1.53ed4641ff204b2cp-52, lgamma(@as(f80, -0x6.005ac9625f233b68p+0)));
    try std.testing.expectEqual(-0x7.313b92969004729p-12, lgamma(@as(f80, -0x6.fff2f8p+0)));
    try std.testing.expectEqual(0x2.a3598cd9f522a41p-12, lgamma(@as(f80, -0x6.fff3p+0)));
    try std.testing.expectEqual(-0x4.dc097be5d1cc3e08p-40, lgamma(@as(f80, -0x6.fff2fddae1bbcp+0)));
    try std.testing.expectEqual(0xe.f46d8dcca9e8197p-48, lgamma(@as(f80, -0x6.fff2fddae1bcp+0)));
    try std.testing.expectEqual(-0x6.9ebebbccaa51db7p-52, lgamma(@as(f80, -0x6.fff2fddae1bbff38p+0)));
    try std.testing.expectEqual(0x3.373d171aaa3ac8cp-52, lgamma(@as(f80, -0x6.fff2fddae1bbff4p+0)));
    try std.testing.expectEqual(0x9.39801333caa3622p-12, lgamma(@as(f80, -0x7.000cf8p+0)));
    try std.testing.expectEqual(-0xa.32834623023dc45p-16, lgamma(@as(f80, -0x7.000dp+0)));
    try std.testing.expectEqual(0x3.89727e62d4843a64p-40, lgamma(@as(f80, -0x7.000cff7b7f878p+0)));
    try std.testing.expectEqual(-0x1.638f6c2b4fb95152p-40, lgamma(@as(f80, -0x7.000cff7b7f87cp+0)));
    try std.testing.expectEqual(0x5.45e474b8c68eb4e8p-52, lgamma(@as(f80, -0x7.000cff7b7f87adfp+0)));
    try std.testing.expectEqual(-0x4.941f6063775a1f68p-52, lgamma(@as(f80, -0x7.000cff7b7f87adf8p+0)));
    try std.testing.expectEqual(-0x4.cccb8849515a9e48p-8, lgamma(@as(f80, -0x7.fffe58p+0)));
    try std.testing.expectEqual(0x1.37b05f6d428d9a9ap-12, lgamma(@as(f80, -0x7.fffe6p+0)));
    try std.testing.expectEqual(-0x2.551849c02b7e0c5cp-40, lgamma(@as(f80, -0x7.fffe5fe05673cp+0)));
    try std.testing.expectEqual(0x2.509d5b2dadf1ea4p-36, lgamma(@as(f80, -0x7.fffe5fe05674p+0)));
    try std.testing.expectEqual(-0x1.9c7a33ad9478c56cp-48, lgamma(@as(f80, -0x7.fffe5fe05673c3c8p+0)));
    try std.testing.expectEqual(0x3.4f638be5777634c4p-48, lgamma(@as(f80, -0x7.fffe5fe05673c3dp+0)));
    try std.testing.expectEqual(0xc.8602745a44910cdp-16, lgamma(@as(f80, -0x8.0001ap+0)));
    try std.testing.expectEqual(-0x9.9cf5dfb6141ef53p-8, lgamma(@as(f80, -0x8.0001bp+0)));
    try std.testing.expectEqual(0x1.34e935f3e5a5cd1p-36, lgamma(@as(f80, -0x8.0001a01459fc8p+0)));
    try std.testing.expectEqual(-0x3.b73909c1555516bp-36, lgamma(@as(f80, -0x8.0001a01459fdp+0)));
    try std.testing.expectEqual(0x7.d0d61593cad19968p-52, lgamma(@as(f80, -0x8.0001a01459fc9f6p+0)));
    try std.testing.expectEqual(-0x9.5b371e11feb316fp-48, lgamma(@as(f80, -0x8.0001a01459fc9f7p+0)));
    try std.testing.expectEqual(-0x9.98ed0cd062e3fd4p-8, lgamma(@as(f80, -0x8.ffffdp+0)));
    try std.testing.expectEqual(0x5.e337e9ef84f0aaap-4, lgamma(@as(f80, -0x8.ffffep+0)));
    try std.testing.expectEqual(-0x5.88479ad476d496a8p-36, lgamma(@as(f80, -0x8.ffffd1c425e8p+0)));
    try std.testing.expectEqual(0x2.6c3945e213fff55cp-32, lgamma(@as(f80, -0x8.ffffd1c425e88p+0)));
    try std.testing.expectEqual(-0x4.55973b8ddaa3ac88p-44, lgamma(@as(f80, -0x8.ffffd1c425e80ffp+0)));
    try std.testing.expectEqual(0x1.33e4438b1b9d4bdcp-44, lgamma(@as(f80, -0x8.ffffd1c425e81p+0)));
    try std.testing.expectEqual(0x5.e32ee82416adc46p-4, lgamma(@as(f80, -0x9.00002p+0)));
    try std.testing.expectEqual(-0x9.99c537e2b92992bp-8, lgamma(@as(f80, -0x9.00003p+0)));
    try std.testing.expectEqual(0x2.5debd4969bb286fp-36, lgamma(@as(f80, -0x9.00002e3bb47d8p+0)));
    try std.testing.expectEqual(-0x2.9ee383255c86df28p-32, lgamma(@as(f80, -0x9.00002e3bb47ep+0)));
    try std.testing.expectEqual(0x2.5e69b52fd9e19d14p-44, lgamma(@as(f80, -0x9.00002e3bb47d86dp+0)));
    try std.testing.expectEqual(-0x3.2b1acbb48b0afdbp-44, lgamma(@as(f80, -0x9.00002e3bb47d86ep+0)));
    try std.testing.expectEqual(-0x1.3dd0c34d79694344p+0, lgamma(@as(f80, -0x9.fffffp+0)));
    try std.testing.expectEqual(-0x1.41334d2c3ccaa62ap-28, lgamma(@as(f80, -0x9.fffffb606bdf8p+0)));
    try std.testing.expectEqual(0x7.9c48d283217d793p-32, lgamma(@as(f80, -0x9.fffffb606bep+0)));
    try std.testing.expectEqual(-0x1.55818a2b42ba2174p-44, lgamma(@as(f80, -0x9.fffffb606bdfdcdp+0)));
    try std.testing.expectEqual(0x3.60979c1bc0b3232cp-40, lgamma(@as(f80, -0x9.fffffb606bdfdcep+0)));
    try std.testing.expectEqual(-0x1.3dd10e8f080e8daap+0, lgamma(@as(f80, -0xa.00001p+0)));
    try std.testing.expectEqual(0x5.70ddf269e6d667ap-32, lgamma(@as(f80, -0xa.0000049f93bb8p+0)));
    try std.testing.expectEqual(-0x1.63ea466b9e05b9e2p-28, lgamma(@as(f80, -0xa.0000049f93bcp+0)));
    try std.testing.expectEqual(0x1.aa9c2e2b1029c57ep-40, lgamma(@as(f80, -0xa.0000049f93bb992p+0)));
    try std.testing.expectEqual(-0x1.cb541d167c13dafcp-40, lgamma(@as(f80, -0xa.0000049f93bb993p+0)));
    try std.testing.expectEqual(-0x3.a3ad38c9033a659cp+0, lgamma(@as(f80, -0xa.fffffp+0)));
    try std.testing.expectEqual(-0x1.0e8528e5ba92d3d6p-24, lgamma(@as(f80, -0xa.ffffff9466e98p+0)));
    try std.testing.expectEqual(0x2.205541c47450d1d4p-28, lgamma(@as(f80, -0xa.ffffff9466eap+0)));
    try std.testing.expectEqual(-0x8.28300f9cbbc503bp-40, lgamma(@as(f80, -0xa.ffffff9466e9f1bp+0)));
    try std.testing.expectEqual(0x1.de91fa23a9940dd2p-36, lgamma(@as(f80, -0xa.ffffff9466e9f1cp+0)));
    try std.testing.expectEqual(-0x3.a3ad86f34c0e3ba4p+0, lgamma(@as(f80, -0xb.00001p+0)));
    try std.testing.expectEqual(0x7.573b0669604304ap-28, lgamma(@as(f80, -0xb.0000006b9915p+0)));
    try std.testing.expectEqual(-0xb.b16d1e1508e7a9cp-28, lgamma(@as(f80, -0xb.0000006b99158p+0)));
    try std.testing.expectEqual(0x2.053cabc3adfebe3p-36, lgamma(@as(f80, -0xb.0000006b9915315p+0)));
    try std.testing.expectEqual(-0x5.bd8591f162c0dab8p-40, lgamma(@as(f80, -0xb.0000006b9915316p+0)));
    try std.testing.expectEqual(-0x6.1fd00f0e21b3c988p+0, lgamma(@as(f80, -0xb.fffffp+0)));
    try std.testing.expectEqual(-0xc.e27c4f01cf53284p-28, lgamma(@as(f80, -0xb.fffffff708938p+0)));
    try std.testing.expectEqual(0xd.785692eee5fd5cp-24, lgamma(@as(f80, -0xb.fffffff70894p+0)));
    try std.testing.expectEqual(-0xf.272276e2f7d5552p-36, lgamma(@as(f80, -0xb.fffffff70893873p+0)));
    try std.testing.expectEqual(0xd.65d9840e2817355p-36, lgamma(@as(f80, -0xb.fffffff70893874p+0)));
    try std.testing.expectEqual(-0x6.1fd05fe315324a38p+0, lgamma(@as(f80, -0xc.00001p+0)));
    try std.testing.expectEqual(0xd.4b0a2023492b1c3p-24, lgamma(@as(f80, -0xc.00000008f76cp+0)));
    try std.testing.expectEqual(-0xf.b743a426163665bp-28, lgamma(@as(f80, -0xc.00000008f76c8p+0)));
    try std.testing.expectEqual(0x2.6322ea559f93a0b8p-36, lgamma(@as(f80, -0xc.00000008f76c773p+0)));
    try std.testing.expectEqual(-0x1.a29d91aa27903fb6p-32, lgamma(@as(f80, -0xc.00000008f76c774p+0)));
    try std.testing.expectEqual(-0x8.b07093393f8bec6p+0, lgamma(@as(f80, -0xc.fffffp+0)));
    try std.testing.expectEqual(-0x7.316d886018815098p-20, lgamma(@as(f80, -0xc.ffffffff4f6d8p+0)));
    try std.testing.expectEqual(0x4.67d7d4d0a160ff8p-20, lgamma(@as(f80, -0xc.ffffffff4f6ep+0)));
    try std.testing.expectEqual(-0x2.2c25e6e64d1da5ecp-32, lgamma(@as(f80, -0xc.ffffffff4f6dcf6p+0)));
    try std.testing.expectEqual(0x1.50666d9a11231798p-28, lgamma(@as(f80, -0xc.ffffffff4f6dcf7p+0)));
    try std.testing.expectEqual(-0x8.b070e6845a6ce34p+0, lgamma(@as(f80, -0xd.00001p+0)));
    try std.testing.expectEqual(0x4.679e61ad5162fc8p-20, lgamma(@as(f80, -0xd.00000000b092p+0)));
    try std.testing.expectEqual(-0x7.31a6fbad0e0cc41p-20, lgamma(@as(f80, -0xd.00000000b0928p+0)));
    try std.testing.expectEqual(0x1.16f33a7d23d6cb18p-28, lgamma(@as(f80, -0xd.00000000b092309p+0)));
    try std.testing.expectEqual(-0x5.c35919086cfd4ec8p-32, lgamma(@as(f80, -0xd.00000000b09230ap+0)));
    try std.testing.expectEqual(-0xb.5409d4efa4b70f9p+0, lgamma(@as(f80, -0xd.fffffp+0)));
    try std.testing.expectEqual(-0x5.861824905c091e7p-16, lgamma(@as(f80, -0xd.fffffffff363p+0)));
    try std.testing.expectEqual(0x4.a000dfad124b37a8p-16, lgamma(@as(f80, -0xd.fffffffff3638p+0)));
    try std.testing.expectEqual(-0xe.bcf83d656a15decp-28, lgamma(@as(f80, -0xd.fffffffff36345ap+0)));
    try std.testing.expectEqual(0x5.8f42e4c2cdc7cbcp-28, lgamma(@as(f80, -0xd.fffffffff36345bp+0)));
    try std.testing.expectEqual(-0xb.540a2a83e42a4f9p+0, lgamma(@as(f80, -0xe.00001p+0)));
    try std.testing.expectEqual(0x4.a0009c38d0a82858p-16, lgamma(@as(f80, -0xe.000000000c9c8p+0)));
    try std.testing.expectEqual(-0x5.861868074a4e2958p-16, lgamma(@as(f80, -0xe.000000000c9dp+0)));
    try std.testing.expectEqual(0x5.8b0b8d2a481f47p-28, lgamma(@as(f80, -0xe.000000000c9cba5p+0)));
    try std.testing.expectEqual(-0xe.c12f950349025abp-28, lgamma(@as(f80, -0xe.000000000c9cba6p+0)));
    try std.testing.expectEqual(-0xe.094c9b083ca94dp+0, lgamma(@as(f80, -0xe.fffffp+0)));
    try std.testing.expectEqual(-0x4.c8585a763b9d58p-12, lgamma(@as(f80, -0xe.ffffffffff288p+0)));
    try std.testing.expectEqual(0x4.bb5f60f986f8a8ep-12, lgamma(@as(f80, -0xe.ffffffffff29p+0)));
    try std.testing.expectEqual(-0xe.beef09380560f81p-28, lgamma(@as(f80, -0xe.ffffffffff28c06p+0)));
    try std.testing.expectEqual(0x1.21b8928708bc37b6p-20, lgamma(@as(f80, -0xe.ffffffffff28c07p+0)));
    try std.testing.expectEqual(-0xe.094cf2be9e3eaf2p+0, lgamma(@as(f80, -0xf.00001p+0)));
    try std.testing.expectEqual(0x4.bb5f60afdcccb468p-12, lgamma(@as(f80, -0xf.0000000000d7p+0)));
    try std.testing.expectEqual(-0x4.c8585ac011a47d4p-12, lgamma(@as(f80, -0xf.0000000000d78p+0)));
    try std.testing.expectEqual(0x1.21b848c7158f27a4p-20, lgamma(@as(f80, -0xf.0000000000d73f9p+0)));
    try std.testing.expectEqual(-0xe.bf38c930add7227p-28, lgamma(@as(f80, -0xf.0000000000d73fap+0)));
    try std.testing.expectEqual(-0x1.0cf14f9e783e6b3cp+4, lgamma(@as(f80, -0xf.fffffp+0)));
    try std.testing.expectEqual(-0xe.466b0623a18cfb1p-12, lgamma(@as(f80, -0xf.fffffffffff28p+0)));
    try std.testing.expectEqual(0x8.c4f2f20afce33e2p-8, lgamma(@as(f80, -0xf.fffffffffff3p+0)));
    try std.testing.expectEqual(-0x7.318a3fab1e86e0bp-20, lgamma(@as(f80, -0xf.fffffffffff28cp+0)));
    try std.testing.expectEqual(0xb.d5eff885a06ba07p-20, lgamma(@as(f80, -0xf.fffffffffff28c1p+0)));
    try std.testing.expectEqual(-0x1.180879870e33e356p+4, lgamma(@as(f80, -0x1.000002p+4)));
    try std.testing.expectEqual(0x8.c4f2f20ab3ff0edp-8, lgamma(@as(f80, -0x1.000000000000dp+4)));
    try std.testing.expectEqual(-0xa.33ca82bb399dc63p-8, lgamma(@as(f80, -0x1.000000000000ep+4)));
    try std.testing.expectEqual(0x1.edd80cde02fd7df4p-16, lgamma(@as(f80, -0x1.000000000000d73ep+4)));
    try std.testing.expectEqual(-0x7.318a4462081fae6p-20, lgamma(@as(f80, -0x1.000000000000d74p+4)));
    try std.testing.expectEqual(-0x1.455d45b618e1f038p+4, lgamma(@as(f80, -0x1.0ffffep+4)));
    try std.testing.expectEqual(-0x3.be7ffe71389cc268p-4, lgamma(@as(f80, -0x1.0ffffffffffffp+4)));
    try std.testing.expectEqual(-0xc.57773dac63c8289p-16, lgamma(@as(f80, -0x1.0ffffffffffff356p+4)));
    try std.testing.expectEqual(0x1.c19a5332b053695p-12, lgamma(@as(f80, -0x1.0ffffffffffff358p+4)));
    try std.testing.expectEqual(-0x1.455d51292150d8bap+4, lgamma(@as(f80, -0x1.100002p+4)));
    try std.testing.expectEqual(-0x3.be7ffe7138f85aacp-4, lgamma(@as(f80, -0x1.1000000000001p+4)));
    try std.testing.expectEqual(0x1.c19a533267df77f2p-12, lgamma(@as(f80, -0x1.1000000000000ca8p+4)));
    try std.testing.expectEqual(-0xc.57773db0ebbe6efp-16, lgamma(@as(f80, -0x1.1000000000000caap+4)));
    try std.testing.expectEqual(-0x1.739c3c0e7e3dc748p+4, lgamma(@as(f80, -0x1.1ffffep+4)));
    try std.testing.expectEqual(-0x3.1fd7673485ba8a88p+0, lgamma(@as(f80, -0x1.1ffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.b80fd7d902af06e4p-8, lgamma(@as(f80, -0x1.1fffffffffffff4ap+4)));
    try std.testing.expectEqual(0x1.c19a53328e26a91cp-12, lgamma(@as(f80, -0x1.1fffffffffffff4cp+4)));
    try std.testing.expectEqual(-0x1.739c47ba6a3ae8acp+4, lgamma(@as(f80, -0x1.200002p+4)));
    try std.testing.expectEqual(-0x3.1fd7673485c0607cp+0, lgamma(@as(f80, -0x1.2000000000001p+4)));
    try std.testing.expectEqual(0x1.c19a53328a0c3826p-12, lgamma(@as(f80, -0x1.20000000000000b4p+4)));
    try std.testing.expectEqual(-0x2.b80fd7d902f168bp-8, lgamma(@as(f80, -0x1.20000000000000b6p+4)));
    try std.testing.expectEqual(-0x1.a2b8a7ff951d4cd8p+4, lgamma(@as(f80, -0x1.2ffffep+4)));
    try std.testing.expectEqual(-0x6.119e27f51c200b5p+0, lgamma(@as(f80, -0x1.2ffffffffffffp+4)));
    try std.testing.expectEqual(-0xd.bb3fcdf10bfe34bp-8, lgamma(@as(f80, -0x1.2ffffffffffffff6p+4)));
    try std.testing.expectEqual(0x2.b64afc1442844c48p-4, lgamma(@as(f80, -0x1.2ffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x1.a2b8b3e16627e78p+4, lgamma(@as(f80, -0x1.300002p+4)));
    try std.testing.expectEqual(-0x6.119e27f51c25fc38p+0, lgamma(@as(f80, -0x1.3000000000001p+4)));
    try std.testing.expectEqual(0x2.b64afc1442841ccp-4, lgamma(@as(f80, -0x1.3000000000000008p+4)));
    try std.testing.expectEqual(-0xd.bb3fcdf10c01eb4p-8, lgamma(@as(f80, -0x1.300000000000000ap+4)));
    try std.testing.expectEqual(-0x1.d2a72cdce34ac164p+4, lgamma(@as(f80, -0x1.3ffffep+4)));
    try std.testing.expectEqual(-0x9.108677639892289p+0, lgamma(@as(f80, -0x1.3ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.709f6fbd94aaf306p+0, lgamma(@as(f80, -0x1.3ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.d2a738f1e7888f4p+4, lgamma(@as(f80, -0x1.400002p+4)));
    try std.testing.expectEqual(-0x9.108677639898331p+0, lgamma(@as(f80, -0x1.4000000000001p+4)));
    try std.testing.expectEqual(-0x1.709f6fbd94aaf3c8p+0, lgamma(@as(f80, -0x1.4000000000000002p+4)));
    try std.testing.expectEqual(-0x2.035d89ed6121f85cp+4, lgamma(@as(f80, -0x1.4ffffep+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e683b1p+0, lgamma(@as(f80, -0x1.4ffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.7c05424b8a8111cp+0, lgamma(@as(f80, -0x1.4ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.035d9633286bf6f8p+4, lgamma(@as(f80, -0x1.500002p+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e6e5ep+0, lgamma(@as(f80, -0x1.5000000000001p+4)));
    try std.testing.expectEqual(-0x4.7c05424b8a811288p+0, lgamma(@as(f80, -0x1.5000000000000002p+4)));
    try std.testing.expectEqual(-0x2.34d272c496dc021cp+4, lgamma(@as(f80, -0x1.5ffffep+4)));
    try std.testing.expectEqual(-0xf.333ad8d94721201p+0, lgamma(@as(f80, -0x1.5ffffffffffffp+4)));
    try std.testing.expectEqual(-0x7.9353d133433a0268p+0, lgamma(@as(f80, -0x1.5ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.34d27f38e9c8e974p+4, lgamma(@as(f80, -0x1.600002p+4)));
    try std.testing.expectEqual(-0xf.333ad8d947275a4p+0, lgamma(@as(f80, -0x1.6000000000001p+4)));
    try std.testing.expectEqual(-0x7.9353d133433a033p+0, lgamma(@as(f80, -0x1.6000000000000002p+4)));
    try std.testing.expectEqual(-0x2.66fd6ea9f77b79a8p+4, lgamma(@as(f80, -0x1.6ffffep+4)));
    try std.testing.expectEqual(-0x1.255ea98937d9f162p+4, lgamma(@as(f80, -0x1.6ffffffffffffp+4)));
    try std.testing.expectEqual(-0xa.b60390ed79b8038p+0, lgamma(@as(f80, -0x1.6ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.66fd7b4acff91314p+4, lgamma(@as(f80, -0x1.700002p+4)));
    try std.testing.expectEqual(-0x1.255ea98937da5668p+4, lgamma(@as(f80, -0x1.7000000000001p+4)));
    try std.testing.expectEqual(-0xa.b60390ed79b8045p+0, lgamma(@as(f80, -0x1.7000000000000002p+4)));
    try std.testing.expectEqual(-0x2.99d6bd8dc680078p+4, lgamma(@as(f80, -0x1.7ffffep+4)));
    try std.testing.expectEqual(-0x1.5837f8825c33e21ep+4, lgamma(@as(f80, -0x1.7ffffffffffffp+4)));
    try std.testing.expectEqual(-0xd.e398807fbf571ap+0, lgamma(@as(f80, -0x1.7ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.99d6ca5949a84b98p+4, lgamma(@as(f80, -0x1.800002p+4)));
    try std.testing.expectEqual(-0x1.5837f8825c34487ap+4, lgamma(@as(f80, -0x1.8000000000001p+4)));
    try std.testing.expectEqual(-0xd.e398807fbf571adp+0, lgamma(@as(f80, -0x1.8000000000000002p+4)));
    try std.testing.expectEqual(-0x2.cd57416926b9198cp+4, lgamma(@as(f80, -0x1.8ffffep+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374e485p+4, lgamma(@as(f80, -0x1.8ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.11ba0bf7d70fd876p+4, lgamma(@as(f80, -0x1.8ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.cd574e5d9fa3edp+4, lgamma(@as(f80, -0x1.900002p+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374eaff4p+4, lgamma(@as(f80, -0x1.9000000000001p+4)));
    try std.testing.expectEqual(-0x1.11ba0bf7d70fd882p+4, lgamma(@as(f80, -0x1.9000000000000002p+4)));
    try std.testing.expectEqual(-0x3.01786b2b55b39354p+4, lgamma(@as(f80, -0x1.9ffffep+4)));
    try std.testing.expectEqual(-0x1.bfd9a6481783e14ap+4, lgamma(@as(f80, -0x1.9ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.45db35cdb745720ep+4, lgamma(@as(f80, -0x1.9ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.0178784731148e2cp+4, lgamma(@as(f80, -0x1.a00002p+4)));
    try std.testing.expectEqual(-0x1.bfd9a64817844a2ap+4, lgamma(@as(f80, -0x1.a000000000001p+4)));
    try std.testing.expectEqual(-0x1.45db35cdb745721ap+4, lgamma(@as(f80, -0x1.a000000000000002p+4)));
    try std.testing.expectEqual(-0x3.36342a886637ea3cp+4, lgamma(@as(f80, -0x1.affffep+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8d007ap+4, lgamma(@as(f80, -0x1.affffffffffffp+4)));
    try std.testing.expectEqual(-0x1.7a96f53dbe4e91d4p+4, lgamma(@as(f80, -0x1.affffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.363437ca2ea26058p+4, lgamma(@as(f80, -0x1.b00002p+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8d6a88p+4, lgamma(@as(f80, -0x1.b000000000001p+4)));
    try std.testing.expectEqual(-0x1.7a96f53dbe4e91ep+4, lgamma(@as(f80, -0x1.b000000000000002p+4)));
    try std.testing.expectEqual(-0x3.6b84e02349a7940cp+4, lgamma(@as(f80, -0x1.bffffep+4)));
    try std.testing.expectEqual(-0x2.29e61b654b21467p+4, lgamma(@as(f80, -0x1.bffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.afe7aaeaeae2d85ep+4, lgamma(@as(f80, -0x1.bffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.6b84ed89a45b2eb8p+4, lgamma(@as(f80, -0x1.c00002p+4)));
    try std.testing.expectEqual(-0x2.29e61b654b21b1a4p+4, lgamma(@as(f80, -0x1.c000000000001p+4)));
    try std.testing.expectEqual(-0x1.afe7aaeaeae2d86ap+4, lgamma(@as(f80, -0x1.c000000000000002p+4)));
    try std.testing.expectEqual(-0x3.a16551a93dea66acp+4, lgamma(@as(f80, -0x1.cffffep+4)));
    try std.testing.expectEqual(-0x2.5fc68cfce71d836p+4, lgamma(@as(f80, -0x1.cffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.e5c81c8286df15d8p+4, lgamma(@as(f80, -0x1.cffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.a1655f32e810c39p+4, lgamma(@as(f80, -0x1.d00002p+4)));
    try std.testing.expectEqual(-0x2.5fc68cfce71defacp+4, lgamma(@as(f80, -0x1.d000000000001p+4)));
    try std.testing.expectEqual(-0x1.e5c81c8286df15e6p+4, lgamma(@as(f80, -0x1.d000000000000002p+4)));
    try std.testing.expectEqual(-0x3.d7d09f8a4486822p+4, lgamma(@as(f80, -0x1.dffffep+4)));
    try std.testing.expectEqual(-0x2.9631daeefecab874p+4, lgamma(@as(f80, -0x1.dffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.1c336a749e8c4b74p+4, lgamma(@as(f80, -0x1.dffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.d7d0ad3610cf0124p+4, lgamma(@as(f80, -0x1.e00002p+4)));
    try std.testing.expectEqual(-0x2.9631daeefecb25dp+4, lgamma(@as(f80, -0x1.e000000000001p+4)));
    try std.testing.expectEqual(-0x2.1c336a749e8c4b84p+4, lgamma(@as(f80, -0x1.e000000000000002p+4)));
    try std.testing.expectEqual(0x9.a81063e79780375p+0, lgamma(@as(f80, 0x8.8d2d5p+0)));
    try std.testing.expectEqual(0x3.2125f40f9a1beba4p+56, lgamma(@as(f80, 0x1.6a324ap+52)));
    try std.testing.expectEqual(0xb.70d4369f5b4c557p+0, lgamma(@as(f80, 0x9.62f59p+0)));
    try std.testing.expectEqual(0xe.b6cd62d45ad40ddp+0, lgamma(@as(f80, 0xa.d55d7p+0)));
    try std.testing.expectEqual(0xe.b6cd3d7503be73bp+0, lgamma(@as(f80, 0xa.d55d6p+0)));
    try std.testing.expectEqual(0xe.b6cd57db84c9ef4p+0, lgamma(@as(f80, 0xa.d55d6b4d78e28p+0)));
    try std.testing.expectEqual(0xa.41afffa8a98e845p+0, lgamma(@as(f80, 0x8.d6315p+0)));
    try std.testing.expectEqual(0xf.8842748a38e7a7p+0, lgamma(@as(f80, 0xb.2e679p+0)));
    try std.testing.expectEqual(0xf.1d4fd446695d45fp+0, lgamma(@as(f80, 0xb.01191p+0)));
    try std.testing.expectEqual(0xf.76b5167078375cp+0, lgamma(@as(f80, 0xb.26fdap+0)));
    try std.testing.expectEqual(0xf.cbb4eb9c9f4ddefp+0, lgamma(@as(f80, 0xb.4ad0ap+0)));
    try std.testing.expectEqual(0xe.0ed26f91598df35p+24, lgamma(@as(f80, 0xe.7a678p+20)));
    try std.testing.expectEqual(0x1.d9db4ca962b419ep+0, lgamma(@as(f80, -0x2.dea4ccp-4)));
    try std.testing.expectEqual(0x1.da47d6051ae6bf5ep+0, lgamma(@as(f80, -0x2.dd306p-4)));
    try std.testing.expectEqual(0xf.f273df313425f4ep-4, lgamma(@as(f80, -0x1.bdc8bp+0)));
    try std.testing.expectEqual(0x1.950848252d48c05ap+0, lgamma(@as(f80, -0x4.0a82e8p-4)));
    try std.testing.expectEqual(0xf.cc00043a75099f4p-4, lgamma(@as(f80, -0x1.bca67ap+0)));
    try std.testing.expectEqual(-0xb.a18b329b453f2e7p-4, lgamma(@as(f80, -0x3.464468p+0)));
    try std.testing.expectEqual(-0xb.a18c341739da8b3p-4, lgamma(@as(f80, -0x3.46446cp+0)));
    try std.testing.expectEqual(-0xb.a18c21a49016c03p-4, lgamma(@as(f80, -0x3.46446bb6a23aap+0)));
    try std.testing.expectEqual(-0xe.aa75345fa640643p-8, lgamma(@as(f80, -0x3.f3d2c4p+0)));
    try std.testing.expectEqual(-0xe.aa27b7e3f86d49ap-8, lgamma(@as(f80, -0x3.f3d2c8p+0)));
    try std.testing.expectEqual(-0xe.aa7484b49666213p-8, lgamma(@as(f80, -0x3.f3d2c40911814p+0)));

    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd080e48fd4262096p+1032, lgamma(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), lgamma(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(std.math.inf(f128), lgamma(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x2.c4c85fdf473ddb98060f5143178p+1032, lgamma(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x0p+0, lgamma(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0xb.17217f7d1cf79abc9e3b39803f3p-4, lgamma(@as(f128, 0x3p+0)));
    try std.testing.expectEqual(0x9.28682473d0de85eafcab635421f8p-4, lgamma(@as(f128, 0x8p-4)));
    try std.testing.expectEqual(0x4.2c8312a971bbf7287f1c24c96238p-4, lgamma(@as(f128, 0xb.33334p-4)));
    try std.testing.expectEqual(0x4.2c83262ea91954655f5ec6068384p-4, lgamma(@as(f128, 0xb.33333p-4)));
    // try std.testing.expectEqual(0x4.2c832247379c4363b0be5aa54848p-4, lgamma(@as(f128, 0xb.3333333333338p-4)));
    try std.testing.expectEqual(0x4.2c832247379cdf8d6c1618623c58p-4, lgamma(@as(f128, 0xb.333333333333p-4)));
    try std.testing.expectEqual(0x4.2c832247379ca106b69376ea1344p-4, lgamma(@as(f128, 0xb.333333333333334p-4)));
    // try std.testing.expectEqual(0x4.2c832247379ca11a3bcae1e1cae4p-4, lgamma(@as(f128, 0xb.333333333333333p-4)));
    // try std.testing.expectEqual(0x4.2c832247379ca11654596616a624p-4, lgamma(@as(f128, 0xb.3333333333333333333333333338p-4)));
    try std.testing.expectEqual(0x4.2c832247379ca11654596616a62cp-4, lgamma(@as(f128, 0xb.333333333333333333333333333p-4)));
    try std.testing.expectEqual(0x4.2c832247379ca11654596616a53p-4, lgamma(@as(f128, 0xb.33333333333333333333333334p-4)));
    try std.testing.expectEqual(0x4.2c832247379ca11654596616aa1p-4, lgamma(@as(f128, 0xb.3333333333333333333333333p-4)));
    try std.testing.expectEqual(-0x1.5db13c7af7431d54a91acd0484e1p-4, lgamma(@as(f128, 0x1.333334p+0)));
    try std.testing.expectEqual(-0x1.5db1333b26a21d93053dff519c4p-4, lgamma(@as(f128, 0x1.333332p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70cadfd0f2a4555835dp-4, lgamma(@as(f128, 0x1.3333333333334p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c63fe8a632b0ceafap-4, lgamma(@as(f128, 0x1.3333333333333p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c72cebe78070ff51bp-4, lgamma(@as(f128, 0x1.3333333333333334p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c72c57ea76e2cac08p-4, lgamma(@as(f128, 0x1.3333333333333332p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c72cb0b57c9e83e4ap-4, lgamma(@as(f128, 0x1.3333333333333333333333333334p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c72cb0b57c9e83e46p-4, lgamma(@as(f128, 0x1.3333333333333333333333333333p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c72cb0b57c9e83faap-4, lgamma(@as(f128, 0x1.333333333333333333333333338p+0)));
    try std.testing.expectEqual(-0x1.5db138c7d70c72cb0b57c9e83d5ap-4, lgamma(@as(f128, 0x1.33333333333333333333333333p+0)));
    // try std.testing.expectEqual(0x8.8bdd41bf44846050819264e2d57p+60, lgamma(@as(f128, 0x3.8p+56)));
    // try std.testing.expectEqual(0x3.72d02ef880f8c917fc232be05e9p+0, lgamma(@as(f128, 0x8p-8)));
    // try std.testing.expectEqual(0x3.7c0e0ff92f04958709ad5a1ae648p+0, lgamma(@as(f128, -0x8p-8)));
    // try std.testing.expectEqual(0x6.ee500bbb72645fcecb166c9d8c64p+0, lgamma(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x6.ee99edf298bdfe3b9118d8828c58p+0, lgamma(@as(f128, -0x4p-12)));
    try std.testing.expectEqual(0xa.65ae3fffc592bd634ed0d8487708p+0, lgamma(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0xa.65b08f1165271d5bc46c11c53e68p+0, lgamma(@as(f128, -0x2p-16)));
    try std.testing.expectEqual(0xd.dce9d6201e89d6bd62b2e7a79c18p+0, lgamma(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0xd.dce9e898ab8646804e122fd4c418p+0, lgamma(@as(f128, -0x1p-20)));
    try std.testing.expectEqual(0x1.1542456e99b0f24ab2b908b14804p+4, lgamma(@as(f128, 0x8p-28)));
    // try std.testing.expectEqual(0x1.15424577d5f770828dc71d4bb9a7p+4, lgamma(@as(f128, -0x8p-28)));
    // try std.testing.expectEqual(0x1.4cb5ecf08473ea2a0dabf1e4d0e9p+4, lgamma(@as(f128, 0x4p-32)));
    // try std.testing.expectEqual(0x1.4cb5ecf0ce561e1bcc8455ba6e63p+4, lgamma(@as(f128, -0x4p-32)));
    try std.testing.expectEqual(0x1.bb9d3beb8c7d73e6fa81731862c3p+4, lgamma(@as(f128, 0x1p-40)));
    // try std.testing.expectEqual(0x1.bb9d3beb8c8fec73f6f12931575ep+4, lgamma(@as(f128, -0x1p-40)));
    try std.testing.expectEqual(0x2.2a848ae66fa859e9c54803444a0cp+4, lgamma(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x2.2a848ae66fa85e87e8871f31d048p+4, lgamma(@as(f128, -0x4p-52)));
    // try std.testing.expectEqual(0x2.996bd9e152ca0843a1517996911ep+4, lgamma(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0x2.996bd9e152ca0844c8da495d8c8p+4, lgamma(@as(f128, -0x1p-60)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6af1e5287e1d7fp+4, lgamma(@as(f128, 0x1p-64)));
    // try std.testing.expectEqual(0x2.c5c85fdf473de6af30cb14de47a6p+4, lgamma(@as(f128, -0x1p-64)));
    // try std.testing.expectEqual(0x3.085328dc35ebb44f931f409f1868p+4, lgamma(@as(f128, 0x4p-72)));
    try std.testing.expectEqual(0x3.085328dc35ebb44f936922d30a26p+4, lgamma(@as(f128, -0x4p-72)));
    try std.testing.expectEqual(0x4.550915ccdf50b871adcf2276181cp+4, lgamma(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x4.550915ccdf50b871adcf22761944p+4, lgamma(@as(f128, -0x1p-100)));
    try std.testing.expectEqual(0x5.75627cbf9441de28d5e1264d1f18p+4, lgamma(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x5.75627cbf9441de28d5e1264d1f18p+4, lgamma(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x8.aa122b99bea170e35b9e44ec316p+4, lgamma(@as(f128, 0x1p-200)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x8.aa122b99bea170e35b9e44ec316p+4, lgamma(@as(f128, -0x1p-200)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.5a92d6d005c939a38650bac4e7b7p+8, lgamma(@as(f128, 0x1p-500)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.5a92d6d005c939a38650bac4e7b7p+8, lgamma(@as(f128, -0x1p-500)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2.b525ada00b9273470ca17589cf6ep+8, lgamma(@as(f128, 0x1p-1000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.b525ada00b9273470ca17589cf6ep+8, lgamma(@as(f128, -0x1p-1000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2.c4657baf579a47bbcffb06f8dfc4p+8, lgamma(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.c4657baf579a47bbcffb06f8dfc4p+8, lgamma(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0xd.89bc642039dc40633f274bb10d2p+8, lgamma(@as(f128, 0x1p-5000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0xd.89bc642039dc40633f274bb10d2p+8, lgamma(@as(f128, -0x1p-5000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.b1378c84073b880c67e4e97621a4p+12, lgamma(@as(f128, 0x1p-10000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.b1378c84073b880c67e4e97621a4p+12, lgamma(@as(f128, -0x1p-10000)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x2.c5b2319c4843acbff21591e99cccp+12, lgamma(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x2.c5b2319c4843acbff21591e99cccp+12, lgamma(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x2.c86ce2daa80dcdaf0680827cc35ap+12, lgamma(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x2.c86ce2daa80dcdaf0680827cc35ap+12, lgamma(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x2.c86ce2daa80dcdaf0680827cc35ap+12, lgamma(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x2.c877f9fc278aeaa6a13d20b7fcdcp+12, lgamma(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x2.ca8c50440f005913a49acbd2c4e8p+12, lgamma(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x6.74767f33d1dc1d0fc8187877a4c8p+4, lgamma(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x2.e870a88dae386c72b4fd4773c092p+8, lgamma(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x2.c86ce2daa80dcdaf0680827cc35ap+12, lgamma(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x2.c877f9fc278aeaa6a13d20b7fcdcp+12, lgamma(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x2.ca8c50440f005913a49acbd2c4e8p+12, lgamma(@as(f128, -0x4p-16496)));
    // try std.testing.expectEqual(-0x7.d809ecd340fc16da6722ad116694p-4, lgamma(@as(f128, -0x3.ec4298p+0)));
    try std.testing.expectEqual(0xf.ffff14223692bc3c374a35f59b5p+124, lgamma(@as(f128, 0x3.12be0cp+120)));
    try std.testing.expectEqual(0x1.00000ceb5ee8a070db2fe7db5d9p+128, lgamma(@as(f128, 0x3.12be6p+120)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0xf.fffffffffff895ade04ea9c1c858p+1020, lgamma(@as(f128, 0x5.d53649e2d4674p+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x1.000000000000701a0eb2451958d2p+1024, lgamma(@as(f128, 0x5.d53649e2d46c8p+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0x1.0000000000000238eb5387b923bp+1024, lgamma(@as(f128, 0x5.d53649e2d46ap+1012)));
    try std.testing.expectEqual(0xf.ffffffffffff73c0163a7fc51948p+1020, lgamma(@as(f128, 0x5.d53649e2d469cp+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffc010a6fe7fb2849p+1020, lgamma(@as(f128, 0x5.d53649e2d469dbc8p+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffbffaad2a9ff30ae8p+1020, lgamma(@as(f128, 0x5.d53649e2d469dbcp+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffbfffffffffffc39p+1020, lgamma(@as(f128, 0x5.d53649e2d469dbc1f01e99fd52p+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    // try std.testing.expectEqual(0x1.0000000000000238eb5387b923bp+1024, lgamma(@as(f128, 0x5.d53649e2d46ap+1012)));
    try std.testing.expectEqual(0xf.ffffffffffff73c0163a7fc51948p+1020, lgamma(@as(f128, 0x5.d53649e2d469cp+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffc010a6fe7fb2849p+1020, lgamma(@as(f128, 0x5.d53649e2d469dbc8p+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffbffaad2a9ff30ae8p+1020, lgamma(@as(f128, 0x5.d53649e2d469dbcp+1012)));
    try std.testing.expectEqual(0xf.ffffffffffffc0000000000036fp+1020, lgamma(@as(f128, 0x5.d53649e2d469dbc1f01e99fd7cp+1012)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd080e48fd4262096p+1032, lgamma(@as(f128, 0xf.ffffffffffff8p+1020)));
    // try std.testing.expectEqual(0xf.ffffffffffffff093d65feafaa5p+16380, lgamma(@as(f128, 0x5.c6aa645fffef5f5p+16368)));
    try std.testing.expectEqual(0x2.c4c85fdf473ddb98060f5143178p+1032, lgamma(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd080e48fd4262096p+1032, lgamma(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), lgamma(@as(f128, 0x5.c6aa645fffef5ff8p+16368)));
    try std.testing.expectEqual(0x2.c4c85fdf473ddb98060f5143178p+1032, lgamma(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd080e48fd4262096p+1032, lgamma(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), lgamma(@as(f128, 0x5.c6aa645fffef5fbp+16368)));
    try std.testing.expectEqual(0xf.fffffffffffffffd06ecf74e1a6p+16380, lgamma(@as(f128, 0x5.c6aa645fffef5fa8p+16368)));
    try std.testing.expectEqual(0xf.ffffffffffffffffffffffffff8p+16380, lgamma(@as(f128, 0x5.c6aa645fffef5fa912b9b480f7acp+16368)));
    try std.testing.expectEqual(0x2.c4c85fdf473ddb98060f5143178p+1032, lgamma(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    try std.testing.expectEqual(0x5.7b90ba32fdbc16cfd352e91badfcp+132, lgamma(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x2.c4c85fdf473dd080e48fd4262096p+1032, lgamma(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(std.math.inf(f128), lgamma(@as(f128, 0x5.c6aa645fffef5fbp+16368)));
    try std.testing.expectEqual(0xf.fffffffffffffffd06ecf74e1a6p+16380, lgamma(@as(f128, 0x5.c6aa645fffef5fa8p+16368)));
    try std.testing.expectEqual(std.math.inf(f128), lgamma(@as(f128, 0x5.c6aa645fffef5fa912b9b480f8p+16368)));
    try std.testing.expectEqual(0x2.c4c85fdf473ddb98060f5143178p+1032, lgamma(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
    // try std.testing.expectEqual(-0x3.511bca412890968ef5acdaae7dbep-20, lgamma(@as(f128, -0x3.f48e28p+0)));
    // try std.testing.expectEqual(0x1.dd4b54ca863c1a476cbd9fd337c3p-20, lgamma(@as(f128, -0x3.f48e2cp+0)));
    // try std.testing.expectEqual(-0x1.ddc0336980b584d18e3a66026b1p-52, lgamma(@as(f128, -0x3.f48e2a8f85fcap+0)));
    // try std.testing.expectEqual(-0x3.4a0c544eeb21a026dc79de4e099ap-24, lgamma(@as(f128, -0x3.24c1b8p+0)));
    try std.testing.expectEqual(-0x7.78a013681f5b968e1639b3340434p+24, lgamma(@as(f128, -0x7.fffff8p+20)));
    // try std.testing.expectEqual(-0x2.30b2cde569e24b3482adbc59e6aap+56, lgamma(@as(f128, -0xf.ffffffffffff8p+48)));
    // try std.testing.expectEqual(-0x1.55589f2fe510778a31db722e9284p+68, lgamma(@as(f128, -0x7.fffffffffffffff8p+60)));
    try std.testing.expectEqual(-0x8.f8f97a94a1c31ceeb9cc952b33dp+108, lgamma(@as(f128, -0x1.ffffffffffffffffffffffffff8p+104)));
    try std.testing.expectEqual(-0x4.ca1ea7c6bcac53b28539e9281ba4p+116, lgamma(@as(f128, -0xf.fffffffffffffffffffffffffff8p+108)));
    // try std.testing.expectEqual(-0x1.52e42ff102e64be289794d246255p+36, lgamma(@as(f128, -0x1.000000008p+32)));
    // try std.testing.expectEqual(-0x1.52e42ff265ca7bd24518407be41dp+36, lgamma(@as(f128, -0x1.000000018p+32)));
    // try std.testing.expectEqual(0x1.96ee685defb2cf07c13b52ad8c5fp+0, lgamma(@as(f128, -0x4p-4)));
    // try std.testing.expectEqual(0x1.43f89a3f0edd620a79ae69cd4613p+0, lgamma(@as(f128, -0x8p-4)));
    // try std.testing.expectEqual(0x1.93616060ea5dfbc406c13494046bp+0, lgamma(@as(f128, -0xcp-4)));
    // try std.testing.expectEqual(0x1.5dce78ceba7e8baf758b14c78cebp+0, lgamma(@as(f128, -0x1.4p+0)));
    // try std.testing.expectEqual(0xd.c2c0a8c107c323f9f78901044cap-4, lgamma(@as(f128, -0x1.8p+0)));
    try std.testing.expectEqual(0x1.041e656d685779d4a3f404f4e635p+0, lgamma(@as(f128, -0x1.cp+0)));
    try std.testing.expectEqual(0x2.bec33c279fa7df4e0daf52f683dcp+0, lgamma(@as(f128, -0x2.08p+0)));
    try std.testing.expectEqual(0x2.07060e6e8471a4872889bc43cbacp+0, lgamma(@as(f128, -0x2.1p+0)));
    // try std.testing.expectEqual(0x1.99a9fdaac9a1381519be768d0a24p+0, lgamma(@as(f128, -0x2.18p+0)));
    try std.testing.expectEqual(0x1.4b32e6350c0cbcfce3355e8d3eb7p+0, lgamma(@as(f128, -0x2.2p+0)));
    try std.testing.expectEqual(0x1.0e029711cf8dcadbfb31b31203bcp+0, lgamma(@as(f128, -0x2.28p+0)));
    // try std.testing.expectEqual(0xd.c0af3f35d3ca5ff45faa2778d698p-4, lgamma(@as(f128, -0x2.3p+0)));
    try std.testing.expectEqual(0xb.214127b24185c3a55f714cce42dp-4, lgamma(@as(f128, -0x2.38p+0)));
    try std.testing.expectEqual(0x8.e355968bdbc2c19c11f614d8a5ap-4, lgamma(@as(f128, -0x2.4p+0)));
    // try std.testing.expectEqual(0x6.f371c281277c8f59db2107586bf8p-4, lgamma(@as(f128, -0x2.48p+0)));
    try std.testing.expectEqual(0x5.44859a67747f55d25257b423b27p-4, lgamma(@as(f128, -0x2.5p+0)));
    try std.testing.expectEqual(0x3.cd82f61be0057224635e100a5774p-4, lgamma(@as(f128, -0x2.58p+0)));
    try std.testing.expectEqual(0x2.8804abda16ec96fcd236c3350162p-4, lgamma(@as(f128, -0x2.6p+0)));
    // try std.testing.expectEqual(0x1.6f830ebd2f0cb62bd9edb09f5f4bp-4, lgamma(@as(f128, -0x2.68p+0)));
    // try std.testing.expectEqual(0x8.0d79aed6889706c84e242cc5979p-8, lgamma(@as(f128, -0x2.7p+0)));
    try std.testing.expectEqual(-0x4.60febffedb540e956d2cd1b5a1ap-8, lgamma(@as(f128, -0x2.78p+0)));
    // try std.testing.expectEqual(-0xe.65fcfaf6878ac4761b616dbe9c28p-8, lgamma(@as(f128, -0x2.8p+0)));
    // try std.testing.expectEqual(-0x1.60773dc36dfb3a2737aebb096e94p-4, lgamma(@as(f128, -0x2.88p+0)));
    // try std.testing.expectEqual(-0x1.b3f01b8343f3228d295d6a35e13ep-4, lgamma(@as(f128, -0x2.9p+0)));
    // try std.testing.expectEqual(-0x1.df97311d4f4d7d7a72d1c691228cp-4, lgamma(@as(f128, -0x2.98p+0)));
    // try std.testing.expectEqual(-0x1.e15351cbe648e7a5981795472499p-4, lgamma(@as(f128, -0x2.ap+0)));
    // try std.testing.expectEqual(-0x1.b5f70616016fabf3429fe8803652p-4, lgamma(@as(f128, -0x2.a8p+0)));
    try std.testing.expectEqual(-0x1.58f3a915176d0a5efef66be2fcfbp-4, lgamma(@as(f128, -0x2.bp+0)));
    // try std.testing.expectEqual(-0xc.3dd1386983f5bc2ded655fb6d13p-8, lgamma(@as(f128, -0x2.b8p+0)));
    // try std.testing.expectEqual(0x1.261e6d250cf634ac23728ff074a3p-8, lgamma(@as(f128, -0x2.cp+0)));
    // try std.testing.expectEqual(0x1.36e062f87a4dd0c9524322e6ec4ap-4, lgamma(@as(f128, -0x2.c8p+0)));
    // try std.testing.expectEqual(0x2.bd203eea3bb29666abf95fe7845ep-4, lgamma(@as(f128, -0x2.dp+0)));
    try std.testing.expectEqual(0x4.c3b22d7ab0718b98997c98b7068cp-4, lgamma(@as(f128, -0x2.d8p+0)));
    try std.testing.expectEqual(0x7.7e1bfe9fdd9f4e8993dbb3f56a64p-4, lgamma(@as(f128, -0x2.ep+0)));
    // try std.testing.expectEqual(0xb.4d46adb8bb95b0de76fefcf48cfp-4, lgamma(@as(f128, -0x2.e8p+0)));
    // try std.testing.expectEqual(0x1.10b1c8eb41e01f2bbe3cd6a4e165p+0, lgamma(@as(f128, -0x2.fp+0)));
    try std.testing.expectEqual(0x1.b6f672f371761ee1bd1431bd6852p+0, lgamma(@as(f128, -0x2.f8p+0)));
    try std.testing.expectEqual(0x1.a2dd71c565b73f6d228bcaa1eadbp+0, lgamma(@as(f128, -0x3.08p+0)));
    // try std.testing.expectEqual(0xe.88018878064a0a862ef5d058f41p-4, lgamma(@as(f128, -0x3.1p+0)));
    // try std.testing.expectEqual(0x7.88aaf3c5b63ce8b3765e44615c58p-4, lgamma(@as(f128, -0x3.18p+0)));
    try std.testing.expectEqual(0x2.780ef1ecfd4bca081f12f293bdd4p-4, lgamma(@as(f128, -0x3.2p+0)));
    // try std.testing.expectEqual(-0x1.83b7ade05f1045749d53035a3a6ap-4, lgamma(@as(f128, -0x3.28p+0)));
    // try std.testing.expectEqual(-0x4.cb8cc177ba556a81c83a394ed6b4p-4, lgamma(@as(f128, -0x3.3p+0)));
    try std.testing.expectEqual(-0x7.92f0f0407d53cff6a7a5bbe0ce18p-4, lgamma(@as(f128, -0x3.38p+0)));
    // try std.testing.expectEqual(-0x9.f86fc0dd02f005f7ad31696ed7p-4, lgamma(@as(f128, -0x3.4p+0)));
    try std.testing.expectEqual(-0xc.0f85e0da3242c1ceb0136cc9832p-4, lgamma(@as(f128, -0x3.48p+0)));
    try std.testing.expectEqual(-0xd.e54537e890f7a838809244d08ca8p-4, lgamma(@as(f128, -0x3.5p+0)));
    // try std.testing.expectEqual(-0xf.82bdb76fac924fc405f972015d18p-4, lgamma(@as(f128, -0x3.58p+0)));
    try std.testing.expectEqual(-0x1.0ee5645b59b4c5f0c17e2103b3c3p+0, lgamma(@as(f128, -0x3.6p+0)));
    try std.testing.expectEqual(-0x1.22c983fd694366382ba9b0e9d6c4p+0, lgamma(@as(f128, -0x3.68p+0)));
    // try std.testing.expectEqual(-0x1.340abce0a1f62ff17b32ac1e2f9cp+0, lgamma(@as(f128, -0x3.7p+0)));
    // try std.testing.expectEqual(-0x1.42ca4c5b0ef6441d453f9c2aa8d2p+0, lgamma(@as(f128, -0x3.78p+0)));
    // try std.testing.expectEqual(-0x1.4f1b0fe64a5d865fa2cc44a4e0c5p+0, lgamma(@as(f128, -0x3.8p+0)));
    // try std.testing.expectEqual(-0x1.59031291fea941652f12d85eed1dp+0, lgamma(@as(f128, -0x3.88p+0)));
    // try std.testing.expectEqual(-0x1.607c0a44539783305e737806a597p+0, lgamma(@as(f128, -0x3.9p+0)));
    // try std.testing.expectEqual(-0x1.6572da73cb38af4f45f25dd2654dp+0, lgamma(@as(f128, -0x3.98p+0)));
    try std.testing.expectEqual(-0x1.67c606af08b9f923cfea3c97ad0ap+0, lgamma(@as(f128, -0x3.ap+0)));
    // try std.testing.expectEqual(-0x1.6742cd4618f50d225aa4764e79cep+0, lgamma(@as(f128, -0x3.a8p+0)));
    try std.testing.expectEqual(-0x1.63a05923d497179fdbe61c0d2c86p+0, lgamma(@as(f128, -0x3.bp+0)));
    // try std.testing.expectEqual(-0x1.5c77fc83c60b44881903014e4e94p+0, lgamma(@as(f128, -0x3.b8p+0)));
    // try std.testing.expectEqual(-0x1.513878cce057f69a43a658ab143ap+0, lgamma(@as(f128, -0x3.cp+0)));
    try std.testing.expectEqual(-0x1.41106fd92d20b08790993f0378c8p+0, lgamma(@as(f128, -0x3.c8p+0)));
    try std.testing.expectEqual(-0x1.2ac7d6f6b00a28569d1e5ad0ed63p+0, lgamma(@as(f128, -0x3.dp+0)));
    // try std.testing.expectEqual(-0x1.0c75b5ade1a5e5d7c50e9c38be78p+0, lgamma(@as(f128, -0x3.d8p+0)));
    // try std.testing.expectEqual(-0xe.2e1c140b222dc36eeeb710644fb8p-4, lgamma(@as(f128, -0x3.ep+0)));
    try std.testing.expectEqual(-0xa.7fd7bc9e5b2e8a6c5847d7ab97ep-4, lgamma(@as(f128, -0x3.e8p+0)));
    try std.testing.expectEqual(-0x4.e2a516e3ce8c25958fc0c743f628p-4, lgamma(@as(f128, -0x3.fp+0)));
    try std.testing.expectEqual(0x5.61445b27ef2f9af42e57f9354204p-4, lgamma(@as(f128, -0x3.f8p+0)));
    try std.testing.expectEqual(-0x2.11f0445d7c7f46e7ccff367e6afp+0, lgamma(@as(f128, -0x4.4p+0)));
    try std.testing.expectEqual(-0x2.d026474418ef5fa1211babb6e74ap+0, lgamma(@as(f128, -0x4.8p+0)));
    // try std.testing.expectEqual(-0x2.e01b099dd31e9182bf31d3e95e68p+0, lgamma(@as(f128, -0x4.cp+0)));
    // try std.testing.expectEqual(-0x3.ba71e6fbceb6724dd3e5f3728e6p+0, lgamma(@as(f128, -0x5.4p+0)));
    try std.testing.expectEqual(-0x4.8490a63c2e095cece2cff1b3e0fcp+0, lgamma(@as(f128, -0x5.8p+0)));
    // try std.testing.expectEqual(-0x4.9fe6996865fd9f4ddc2a8b04b60cp+0, lgamma(@as(f128, -0x5.cp+0)));
    // try std.testing.expectEqual(-0x5.8f95f609dcbdec55ff0dd66e952cp+0, lgamma(@as(f128, -0x6.4p+0)));
    // try std.testing.expectEqual(-0x6.63bf13aa8dc40311e8a61d305cb8p+0, lgamma(@as(f128, -0x6.8p+0)));
    try std.testing.expectEqual(-0x6.88be607932f0a85a34afcbd3bdd8p+0, lgamma(@as(f128, -0x6.cp+0)));
    // try std.testing.expectEqual(-0x7.8ab8df93f8e2d0ab3f5a4d49d3ep+0, lgamma(@as(f128, -0x7.4p+0)));
    // try std.testing.expectEqual(-0x8.678fc2dc64f8698ca2539c036558p+0, lgamma(@as(f128, -0x7.8p+0)));
    try std.testing.expectEqual(-0x8.94f3f99bb4bcf32586bcabb15d7p+0, lgamma(@as(f128, -0x7.cp+0)));
    // try std.testing.expectEqual(-0x9.a6efce3f0c5dfdc1db446d03ceep+0, lgamma(@as(f128, -0x8.4p+0)));
    // try std.testing.expectEqual(-0xa.8b6b2323e31829c0be636f82e6dp+0, lgamma(@as(f128, -0x8.8p+0)));
    try std.testing.expectEqual(-0xa.c03b140e0f96abc4c901806682f8p+0, lgamma(@as(f128, -0x8.cp+0)));
    try std.testing.expectEqual(-0xb.e070bc16c1b6b15d44a869cfc8a8p+0, lgamma(@as(f128, -0x9.4p+0)));
    try std.testing.expectEqual(-0xc.cbbfcbeca7ae3e5503d29e5934fp+0, lgamma(@as(f128, -0x9.8p+0)));
    try std.testing.expectEqual(-0xd.0736112f6db281b4a90d85ap+0, lgamma(@as(f128, -0x9.cp+0)));
    try std.testing.expectEqual(-0xe.343934d8f3a1737b4ce05d06fe28p+0, lgamma(@as(f128, -0xa.4p+0)));
    // try std.testing.expectEqual(-0xf.25b38682cbb4e366d49d0ee55c5p+0, lgamma(@as(f128, -0xa.8p+0)));
    // try std.testing.expectEqual(-0xf.672fe4026795b128b9dd0a7c2c48p+0, lgamma(@as(f128, -0xa.cp+0)));
    try std.testing.expectEqual(-0x1.09fd673bdc93709c0e0c3b597081p+4, lgamma(@as(f128, -0xb.4p+0)));
    try std.testing.expectEqual(-0x1.196f12e4530636addbb797998b7fp+4, lgamma(@as(f128, -0xb.8p+0)));
    try std.testing.expectEqual(-0x1.1ddeefa04e20d8902b3ea2985287p+4, lgamma(@as(f128, -0xb.cp+0)));
    try std.testing.expectEqual(-0x1.32140999470e300f73a257c054c6p+4, lgamma(@as(f128, -0xc.4p+0)));
    try std.testing.expectEqual(-0x1.41d87554b103a5e91b085102cc2bp+4, lgamma(@as(f128, -0xc.8p+0)));
    try std.testing.expectEqual(-0x1.46996e9ff5e8e7901aa2fd6c1ab3p+4, lgamma(@as(f128, -0xc.cp+0)));
    // try std.testing.expectEqual(-0x1.5b6c176a914d9642f7b1b82c4984p+4, lgamma(@as(f128, -0xd.4p+0)));
    try std.testing.expectEqual(-0x1.6b7d13453aefce149d2ee0493ce6p+4, lgamma(@as(f128, -0xd.8p+0)));
    try std.testing.expectEqual(-0x1.70893507e7aac335181780e3ca85p+4, lgamma(@as(f128, -0xd.cp+0)));
    // try std.testing.expectEqual(-0x1.85ee2af24d7d0a88e9ac08b57e7ap+4, lgamma(@as(f128, -0xe.4p+0)));
    try std.testing.expectEqual(-0x1.9646635d59cf13f4add1e2f07111p+4, lgamma(@as(f128, -0xe.8p+0)));
    try std.testing.expectEqual(-0x1.9b9889f00a16b6da301362abbacep+4, lgamma(@as(f128, -0xe.cp+0)));
    try std.testing.expectEqual(-0x1.b1860b9f9cf34eda33665e357553p+4, lgamma(@as(f128, -0xf.4p+0)));
    try std.testing.expectEqual(-0x1.c220de6eff08d03c1f90ec27cb49p+4, lgamma(@as(f128, -0xf.8p+0)));
    // try std.testing.expectEqual(-0x1.c7b48e949c3d3427fac367504d59p+4, lgamma(@as(f128, -0xf.cp+0)));
    // try std.testing.expectEqual(-0x1.de2212eef35f350cc51d00051d45p+4, lgamma(@as(f128, -0x1.04p+4)));
    try std.testing.expectEqual(-0x1.eefb6ed92d5d7aa845edc95ceb38p+4, lgamma(@as(f128, -0x1.08p+4)));
    // try std.testing.expectEqual(-0x1.f4ccb75a4247fa751ee3e945c0a9p+4, lgamma(@as(f128, -0x1.0cp+4)));
    try std.testing.expectEqual(-0x2.0bb2b66649903080e12e244c2314p+4, lgamma(@as(f128, -0x1.14p+4)));
    try std.testing.expectEqual(-0x2.1cc701ffd0280dccf6b051e1bddp+4, lgamma(@as(f128, -0x1.18p+4)));
    try std.testing.expectEqual(-0x2.22d2642bdb692f166219f8d5855p+4, lgamma(@as(f128, -0x1.1cp+4)));
    try std.testing.expectEqual(-0x2.3a2a2c33d815da18a695407a3c66p+4, lgamma(@as(f128, -0x1.24p+4)));
    try std.testing.expectEqual(-0x2.4b76325cc89a90a169e4cce7fdacp+4, lgamma(@as(f128, -0x1.28p+4)));
    try std.testing.expectEqual(-0x2.51b88f97694cb14e4f0e0fda961p+4, lgamma(@as(f128, -0x1.2cp+4)));
    try std.testing.expectEqual(-0x2.697c23520ea4d9a7157b930d7e24p+4, lgamma(@as(f128, -0x1.34p+4)));
    try std.testing.expectEqual(-0x2.7afd03ae5b99459b2483c87515bcp+4, lgamma(@as(f128, -0x1.38p+4)));
    try std.testing.expectEqual(-0x2.81738ebf2dd88e145b52ffe790acp+4, lgamma(@as(f128, -0x1.3cp+4)));
    try std.testing.expectEqual(-0x2.999d8a3dc87714cf45457fefbef4p+4, lgamma(@as(f128, -0x1.44p+4)));
    // try std.testing.expectEqual(-0x2.ab50acb9fbd4e957c1a582e20954p+4, lgamma(@as(f128, -0x1.48p+4)));
    try std.testing.expectEqual(-0x2.b1f8ddf5bf30a5572ac9d4d3adacp+4, lgamma(@as(f128, -0x1.4cp+4)));
    try std.testing.expectEqual(-0x2.ca8460bab0c94ca2c85fbc3f9742p+4, lgamma(@as(f128, -0x1.54p+4)));
    try std.testing.expectEqual(-0x2.dc676b66a89013e9bf50b6694c56p+4, lgamma(@as(f128, -0x1.58p+4)));
    try std.testing.expectEqual(-0x2.e33ef7090df5fe33e9103516b1eap+4, lgamma(@as(f128, -0x1.5cp+4)));
    try std.testing.expectEqual(-0x2.fc27921a70bb36502d28015982fap+4, lgamma(@as(f128, -0x1.64p+4)));
    try std.testing.expectEqual(-0x3.0e3860d4730664e8d52d272bcd34p+4, lgamma(@as(f128, -0x1.68p+4)));
    try std.testing.expectEqual(-0x3.153d2f0ea92f084fec38a601ebc8p+4, lgamma(@as(f128, -0x1.6cp+4)));
    // try std.testing.expectEqual(-0x3.2e7ed62745db05944c8a682cad48p+4, lgamma(@as(f128, -0x1.74p+4)));
    try std.testing.expectEqual(-0x3.40bb73b417cada01316c3455dd36p+4, lgamma(@as(f128, -0x1.78p+4)));
    try std.testing.expectEqual(-0x3.47eb9a13a5e8a56971e8d80730e2p+4, lgamma(@as(f128, -0x1.7cp+4)));
    try std.testing.expectEqual(-0x3.6182974be0d0f6629117525e7796p+4, lgamma(@as(f128, -0x1.84p+4)));
    // try std.testing.expectEqual(-0x3.73e93790ff62910f53a08bf641bap+4, lgamma(@as(f128, -0x1.88p+4)));
    // try std.testing.expectEqual(-0x3.7b42f379042362d245e912d820e4p+4, lgamma(@as(f128, -0x1.8cp+4)));
    try std.testing.expectEqual(-0x3.952bdce9557fca25fb42723f72cap+4, lgamma(@as(f128, -0x1.94p+4)));
    try std.testing.expectEqual(-0x3.a7bad810244797a9ffa322038a26p+4, lgamma(@as(f128, -0x1.98p+4)));
    // try std.testing.expectEqual(-0x3.af3c8a0f6e39233617cadb8d50b4p+4, lgamma(@as(f128, -0x1.9cp+4)));
    try std.testing.expectEqual(-0x3.c974390b2830704759a858601576p+4, lgamma(@as(f128, -0x1.a4p+4)));
    try std.testing.expectEqual(-0x3.dc2a0760eba3f5784050bda8ff22p+4, lgamma(@as(f128, -0x1.a8p+4)));
    try std.testing.expectEqual(-0x3.e3d22f3b711ca1c8360ee4bbdd26p+4, lgamma(@as(f128, -0x1.acp+4)));
    // try std.testing.expectEqual(-0x3.fe55b8d8334a9a313e1792cbd2fep+4, lgamma(@as(f128, -0x1.b4p+4)));
    // try std.testing.expectEqual(-0x4.1130ef485a82c8b7fa637c5a2f34p+4, lgamma(@as(f128, -0x1.b8p+4)));
    // try std.testing.expectEqual(-0x4.18fe28939975379956e6bd5ded14p+4, lgamma(@as(f128, -0x1.bcp+4)));
    try std.testing.expectEqual(-0x4.33cad742071e19dd9fabc343392p+4, lgamma(@as(f128, -0x1.c4p+4)));
    try std.testing.expectEqual(-0x4.46ca244f93cf3498a8fc081ce468p+4, lgamma(@as(f128, -0x1.c8p+4)));
    try std.testing.expectEqual(-0x4.4ebb238830305be106adc3210304p+4, lgamma(@as(f128, -0x1.ccp+4)));
    try std.testing.expectEqual(-0x4.69ce718eca02e1d407ee1c6c2144p+4, lgamma(@as(f128, -0x1.d4p+4)));
    try std.testing.expectEqual(-0x4.7cf09ab733581fd87d96251e54fp+4, lgamma(@as(f128, -0x1.d8p+4)));
    // try std.testing.expectEqual(-0x4.85042abb5d4fb79dfb39ce8c6d4p+4, lgamma(@as(f128, -0x1.dcp+4)));
    try std.testing.expectEqual(-0x4.a05bbd6dcca6217d8024a52bc07cp+4, lgamma(@as(f128, -0x1.e4p+4)));
    try std.testing.expectEqual(-0x4.b39f9ce3ffeb5bc483eeb5d7cc08p+4, lgamma(@as(f128, -0x1.e8p+4)));
    // try std.testing.expectEqual(-0x4.bbd49cc22d716e5745fee69530ecp+4, lgamma(@as(f128, -0x1.ecp+4)));
    try std.testing.expectEqual(-0x4.d76e40569b13cc8900ce9dace15cp+4, lgamma(@as(f128, -0x1.f4p+4)));
    try std.testing.expectEqual(-0x4.ead2c3080f2ed0ad0b3cf5b5ded4p+4, lgamma(@as(f128, -0x1.f8p+4)));
    try std.testing.expectEqual(-0x4.f3282414b3f0877aa37c836452ecp+4, lgamma(@as(f128, -0x1.fcp+4)));
    try std.testing.expectEqual(-0x5.0f01c7fe77b50a17ac1d2ecff474p+4, lgamma(@as(f128, -0x2.04p+4)));
    // try std.testing.expectEqual(-0x5.2285ebd6e2b7ae7a5991d2bf0704p+4, lgamma(@as(f128, -0x2.08p+4)));
    try std.testing.expectEqual(-0x5.2afaaffe44d821eaf73380b599d8p+4, lgamma(@as(f128, -0x2.0cp+4)));
    try std.testing.expectEqual(-0x5.471263b9b93bcb1aa33f2ff0cb78p+4, lgamma(@as(f128, -0x2.14p+4)));
    try std.testing.expectEqual(-0x5.5ab5361c05df6c623a508fedfa94p+4, lgamma(@as(f128, -0x2.18p+4)));
    try std.testing.expectEqual(-0x5.63486e673f3485e6bab34f13eac8p+4, lgamma(@as(f128, -0x2.1cp+4)));
    try std.testing.expectEqual(-0x5.7f9c5ea615044f2ed87c198429e8p+4, lgamma(@as(f128, -0x2.24p+4)));
    try std.testing.expectEqual(-0x5.935cfb12d92d5f7112ffef6e80ap+4, lgamma(@as(f128, -0x2.28p+4)));
    try std.testing.expectEqual(-0x5.9c0dc658a126dfe2fc435dba0efcp+4, lgamma(@as(f128, -0x2.2cp+4)));
    // try std.testing.expectEqual(-0x5.b89c3a80e9aed743e25e844a0ap+4, lgamma(@as(f128, -0x2.34p+4)));
    // try std.testing.expectEqual(-0x5.cc79c963ef6b8bad12d43a37c588p+4, lgamma(@as(f128, -0x2.38p+4)));
    // try std.testing.expectEqual(-0x5.d547531f08742a1a08bc04d1a654p+4, lgamma(@as(f128, -0x2.3cp+4)));
    try std.testing.expectEqual(-0x5.f20eab1178fe58f4345ac6091e6p+4, lgamma(@as(f128, -0x2.44p+4)));
    try std.testing.expectEqual(-0x6.060860b0fb0e2cdf94d9919f5f18p+4, lgamma(@as(f128, -0x2.48p+4)));
    try std.testing.expectEqual(-0x6.0ef1dff71ff1f424d893ba0ddaecp+4, lgamma(@as(f128, -0x2.4cp+4)));
    // try std.testing.expectEqual(-0x6.2bf09212ee611d815d24c6d465c4p+4, lgamma(@as(f128, -0x2.54p+4)));
    try std.testing.expectEqual(-0x6.4005ad9c060ea6b23e6be3ddf018p+4, lgamma(@as(f128, -0x2.58p+4)));
    try std.testing.expectEqual(-0x6.490a643105b9511ee42b1d09d58cp+4, lgamma(@as(f128, -0x2.5cp+4)));
    // try std.testing.expectEqual(-0x6.663efb8d432c3718b6caba7d7e44p+4, lgamma(@as(f128, -0x2.64p+4)));
    try std.testing.expectEqual(-0x6.7a6ec639b9ba9ddb69f071aab214p+4, lgamma(@as(f128, -0x2.68p+4)));
    try std.testing.expectEqual(-0x6.838dffbb1b634936974365590ffcp+4, lgamma(@as(f128, -0x2.6cp+4)));
    try std.testing.expectEqual(-0x6.a0f71a8eb113b9818cc418ac52bp+4, lgamma(@as(f128, -0x2.74p+4)));
    try std.testing.expectEqual(-0x6.b540e6e0fb63723c32d39cf12cfp+4, lgamma(@as(f128, -0x2.78p+4)));
    try std.testing.expectEqual(-0x6.be79f80712a5ba0185945e9a550cp+4, lgamma(@as(f128, -0x2.7cp+4)));
    try std.testing.expectEqual(-0x6.dc1646398c9c01b2adfced8afa8cp+4, lgamma(@as(f128, -0x2.84p+4)));
    try std.testing.expectEqual(-0x6.f0796f4c3252a4ff1f3bc50ceep+4, lgamma(@as(f128, -0x2.88p+4)));
    // try std.testing.expectEqual(-0x6.f9cbb53dffc1b6d177c75140b75p+4, lgamma(@as(f128, -0x2.8cp+4)));
    try std.testing.expectEqual(-0x7.1799f71c2b60e7ef15b309d7fabp+4, lgamma(@as(f128, -0x2.94p+4)));
    try std.testing.expectEqual(-0x7.2c15e00240c7b3dcab50d5328b4p+4, lgamma(@as(f128, -0x2.98p+4)));
    try std.testing.expectEqual(-0x7.3580bfb9dce5421d4ec6a4336554p+4, lgamma(@as(f128, -0x2.9cp+4)));
    // try std.testing.expectEqual(-0x7.537fc4c9f7583cb3b66dcf478a28p+4, lgamma(@as(f128, -0x2.a4p+4)));
    try std.testing.expectEqual(-0x7.6813d7fea636e34aeb094cbbe3ccp+4, lgamma(@as(f128, -0x2.a8p+4)));
    try std.testing.expectEqual(-0x7.7196bdbc4617c0faab028d91ea9cp+4, lgamma(@as(f128, -0x2.acp+4)));
    try std.testing.expectEqual(-0x7.8fc563ae0f0867eef6decb627e68p+4, lgamma(@as(f128, -0x2.b4p+4)));
    try std.testing.expectEqual(-0x7.a4711291721933c265ede8386848p+4, lgamma(@as(f128, -0x2.b8p+4)));
    try std.testing.expectEqual(-0x7.ae0b715b59528fef9d1e552e9d6cp+4, lgamma(@as(f128, -0x2.bcp+4)));
    try std.testing.expectEqual(-0x7.cc68a310de77662b791c45ffbf08p+4, lgamma(@as(f128, -0x2.c4p+4)));
    try std.testing.expectEqual(-0x7.e12b6570af28150a8754688bd444p+4, lgamma(@as(f128, -0x2.c8p+4)));
    try std.testing.expectEqual(-0x7.eadcb69e9c3a010a72e32fd184b4p+4, lgamma(@as(f128, -0x2.ccp+4)));
    try std.testing.expectEqual(-0x8.09676b4afe7a217d040b36c8c0b8p+4, lgamma(@as(f128, -0x2.d4p+4)));
    try std.testing.expectEqual(-0x8.1e40bef5c77e16c1471b14b08e6p+4, lgamma(@as(f128, -0x2.d8p+4)));
    try std.testing.expectEqual(-0x8.280881c698b34ff326df1e26cbep+4, lgamma(@as(f128, -0x2.dcp+4)));
    try std.testing.expectEqual(-0x8.46bfbc20675ce021b898f0e6e0c8p+4, lgamma(@as(f128, -0x2.e4p+4)));
    // try std.testing.expectEqual(-0x8.5baf248219badda0231bb6bd38fp+4, lgamma(@as(f128, -0x2.e8p+4)));
    // try std.testing.expectEqual(-0x8.658cddba91e6ebcb24bb5fbe939p+4, lgamma(@as(f128, -0x2.ecp+4)));
    try std.testing.expectEqual(-0x8.846fab3fa6866806ed8a8e60c3a8p+4, lgamma(@as(f128, -0x2.f4p+4)));
    // try std.testing.expectEqual(-0x8.9974b10693917254656a23fbfe4p+4, lgamma(@as(f128, -0x2.f8p+4)));
    try std.testing.expectEqual(-0x8.a367ea98496fe483be1e42eff598p+4, lgamma(@as(f128, -0x2.fcp+4)));
    // try std.testing.expectEqual(-0x8.c27562e1518600979bcb6443f218p+4, lgamma(@as(f128, -0x3.04p+4)));
    try std.testing.expectEqual(-0x8.d78f93aaaba45abd6695496748dp+4, lgamma(@as(f128, -0x3.08p+4)));
    try std.testing.expectEqual(-0x8.e197dc624cded54dba167d94edd8p+4, lgamma(@as(f128, -0x3.0cp+4)));
    try std.testing.expectEqual(-0x9.00cf208467db1573aecee045af18p+4, lgamma(@as(f128, -0x3.14p+4)));
    try std.testing.expectEqual(-0x9.15fe0e8f86fc0fc0f733bf71b91p+4, lgamma(@as(f128, -0x3.18p+4)));
    try std.testing.expectEqual(-0x9.201af9c9b1dafeb0561df89315bp+4, lgamma(@as(f128, -0x3.1cp+4)));
    try std.testing.expectEqual(-0x9.3f7b33c4bae8e6583d30fb1072cp+4, lgamma(@as(f128, -0x3.24p+4)));
    // try std.testing.expectEqual(-0x9.54be75ac78c7db1f1dfd1a8c3488p+4, lgamma(@as(f128, -0x3.28p+4)));
    try std.testing.expectEqual(-0x9.5eef9b1085f7a44a198c096dbcp+4, lgamma(@as(f128, -0x3.2cp+4)));
    try std.testing.expectEqual(-0x9.7e77fd48cb94c5e51babf495469p+4, lgamma(@as(f128, -0x3.34p+4)));
    try std.testing.expectEqual(-0x9.93cf2dc25ffa931dac7d1e7ae49p+4, lgamma(@as(f128, -0x3.38p+4)));
    // try std.testing.expectEqual(-0x9.9e142902892baa559fdec68a122p+4, lgamma(@as(f128, -0x3.3cp+4)));
    try std.testing.expectEqual(-0x9.bdc3edc4d92fc7031c1e2be27b08p+4, lgamma(@as(f128, -0x3.44p+4)));
    try std.testing.expectEqual(-0x9.d32eab63afc830d9c7813fd5078p+4, lgamma(@as(f128, -0x3.48p+4)));
    // try std.testing.expectEqual(-0x9.dd871c0210b9a670b4e9fc3cb508p+4, lgamma(@as(f128, -0x3.4cp+4)));
    // try std.testing.expectEqual(-0x9.fd5d85111f54063bc995d4479008p+4, lgamma(@as(f128, -0x3.54p+4)));
    try std.testing.expectEqual(-0xa.12db720f2fc8a706a263843d143p+4, lgamma(@as(f128, -0x3.58p+4)));
    try std.testing.expectEqual(-0xa.1d46fb272de50cd2f3ee6edbd848p+4, lgamma(@as(f128, -0x3.5cp+4)));
    // try std.testing.expectEqual(-0xa.3d43515179cb223cee1febb25c6p+4, lgamma(@as(f128, -0x3.64p+4)));
    try std.testing.expectEqual(-0xa.52d4135bb7ffc88b4370f9e251f8p+4, lgamma(@as(f128, -0x3.68p+4)));
    // try std.testing.expectEqual(-0xa.5d525b6f696dc102239723df41e8p+4, lgamma(@as(f128, -0x3.6cp+4)));
    try std.testing.expectEqual(-0xa.7d73ee2cd7a8c8a2803e4ba21c3p+4, lgamma(@as(f128, -0x3.74p+4)));
    // try std.testing.expectEqual(-0xa.93172e335d7555f720e70dbde228p+4, lgamma(@as(f128, -0x3.78p+4)));
    try std.testing.expectEqual(-0xa.9da7defc939ca9661d2ac26ea158p+4, lgamma(@as(f128, -0x3.7cp+4)));
    try std.testing.expectEqual(-0xa.bdee0413128f5571d773fb7c7d48p+4, lgamma(@as(f128, -0x3.84p+4)));
    // try std.testing.expectEqual(-0xa.d3a36e1cae65cd3e3f19796ec888p+4, lgamma(@as(f128, -0x3.88p+4)));
    try std.testing.expectEqual(-0xa.de46346151a96509ce5f57f6a478p+4, lgamma(@as(f128, -0x3.8cp+4)));
    // try std.testing.expectEqual(-0xa.feb0478fe57870dbb4f6896ee5p+4, lgamma(@as(f128, -0x3.94p+4)));
    try std.testing.expectEqual(-0xb.14778a90c23de920ab7eba6b5eb8p+4, lgamma(@as(f128, -0x3.98p+4)));
    // try std.testing.expectEqual(-0xb.1f2c15fa353b6f2fa787435f3688p+4, lgamma(@as(f128, -0x3.9cp+4)));
    try std.testing.expectEqual(-0xb.3fb978a9e017ff7a66b90ccbf7a8p+4, lgamma(@as(f128, -0x3.a4p+4)));
    // try std.testing.expectEqual(-0xb.5592465d023fa8b1d05f4ecdc718p+4, lgamma(@as(f128, -0x3.a8p+4)));
    try std.testing.expectEqual(-0xb.605849524a702018451f8fa126cp+4, lgamma(@as(f128, -0x3.acp+4)));
    try std.testing.expectEqual(-0xb.8108624c51a6e6d7dd3f95ae4298p+4, lgamma(@as(f128, -0x3.b4p+4)));
    // try std.testing.expectEqual(-0xb.96f26f0fac7bfc0981899572b19p+4, lgamma(@as(f128, -0x3.b8p+4)));
    // try std.testing.expectEqual(-0xb.a1c99e9224b8974b9b855feefdep+4, lgamma(@as(f128, -0x3.bcp+4)));
    try std.testing.expectEqual(-0xb.c29bd9bb401ef09ac5948b7e3aep+4, lgamma(@as(f128, -0x3.c4p+4)));
    // try std.testing.expectEqual(-0xb.d896dc6e2c3c334db65e596bd108p+4, lgamma(@as(f128, -0x3.c8p+4)));
    try std.testing.expectEqual(-0xb.e37eeff88b8ddd0b3590e5a68fp+4, lgamma(@as(f128, -0x3.ccp+4)));
    try std.testing.expectEqual(0x1.0a2b23fa7e70cd72a6f928bada49p+4, lgamma(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x2.4bc9ef64e6ff433f2a8e5128b4b6p+4, lgamma(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x2.c5c85fdf473de6af2e5287e1d7fp+4, lgamma(@as(f128, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x4.e535c94639c94b4d41d824619be4p+4, lgamma(@as(f128, -0xf.fffffffffffffffffffffffffff8p-4)));
    // try std.testing.expectEqual(0x4.9793dec9cdfe8612198485cf1a2cp+4, lgamma(@as(f128, -0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(0xf.f140266b6278ff9f51d8bd4f62bp+0, lgamma(@as(f128, -0x1.000002p+0)));
    try std.testing.expectEqual(0x2.40b2cde569e24b02148beb2bbdb2p+4, lgamma(@as(f128, -0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.bab13e5fca20ef145d692022ff42p+4, lgamma(@as(f128, -0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x4.da1ea7c6bcac53b28539e9281ba4p+4, lgamma(@as(f128, -0x1.0000000000000000000000000001p+0)));
    try std.testing.expectEqual(0x4.8c7cbd4a50e18e775ce64a9599e8p+4, lgamma(@as(f128, -0x1.000000000000000000000000008p+0)));
    try std.testing.expectEqual(0xf.3fce11247f0a77fcc417ebc455b8p+0, lgamma(@as(f128, -0x1.fffffep+0)));
    try std.testing.expectEqual(0x2.359bac65ecc554bfcf1de8f6dbe6p+4, lgamma(@as(f128, -0x1.fffffffffffffp+0)));
    try std.testing.expectEqual(0x2.af9a1ce04d03f779cbd9caf09f98p+4, lgamma(@as(f128, -0x1.fffffffffffffffep+0)));
    try std.testing.expectEqual(0x4.cf0786473f8f5c17c89badee9b64p+4, lgamma(@as(f128, -0x1.ffffffffffffffffffffffffffffp+0)));
    try std.testing.expectEqual(0x4.81659bcad3c496dca0480f5c19b4p+4, lgamma(@as(f128, -0x1.ffffffffffffffffffffffffff8p+0)));
    try std.testing.expectEqual(0xe.8e5bf3a347bbb1e1859aa88e0cfp+0, lgamma(@as(f128, -0x2.000004p+0)));
    try std.testing.expectEqual(0x2.2a848ae66fa85a605fb758366e44p+4, lgamma(@as(f128, -0x2.0000000000002p+0)));
    try std.testing.expectEqual(0x2.a482fb60cfe6ffdeb6a536ac6e7ap+4, lgamma(@as(f128, -0x2.0000000000000004p+0)));
    try std.testing.expectEqual(0x4.c3f064c7c272647d0bfd72b51b28p+4, lgamma(@as(f128, -0x2.0000000000000000000000000002p+0)));
    try std.testing.expectEqual(0x4.764e7a4b56a79f41e3a9d422995cp+4, lgamma(@as(f128, -0x2.00000000000000000000000001p+0)));
    try std.testing.expectEqual(0xd.751d54afa9a22560e6fd730a2a3p+0, lgamma(@as(f128, -0x2.fffffcp+0)));
    try std.testing.expectEqual(0x2.18f0a06bc2a554248e80dc15058p+4, lgamma(@as(f128, -0x2.ffffffffffffep+0)));
    try std.testing.expectEqual(0x2.92ef10e622e3f547d7d6bf3ab4a4p+4, lgamma(@as(f128, -0x2.fffffffffffffffcp+0)));
    try std.testing.expectEqual(0x4.b25c7a4d156f59e5a1bbd9dfcad4p+4, lgamma(@as(f128, -0x2.fffffffffffffffffffffffffffep+0)));
    try std.testing.expectEqual(0x4.64ba8fd0a9a494aa79683b4d492cp+4, lgamma(@as(f128, -0x2.ffffffffffffffffffffffffffp+0)));
    // try std.testing.expectEqual(0xd.751d4aa3223696a3c4450e957978p+0, lgamma(@as(f128, -0x3.000004p+0)));
    try std.testing.expectEqual(0x2.18f0a06bc2a54f1e4acb14b67348p+4, lgamma(@as(f128, -0x3.0000000000002p+0)));
    try std.testing.expectEqual(0x2.92ef10e622e3f547370e4881c8d2p+4, lgamma(@as(f128, -0x3.0000000000000004p+0)));
    try std.testing.expectEqual(0x4.b25c7a4d156f59e5a1bbd9dfcad4p+4, lgamma(@as(f128, -0x3.0000000000000000000000000002p+0)));
    try std.testing.expectEqual(0x4.64ba8fd0a9a494aa79683b4d4904p+4, lgamma(@as(f128, -0x3.00000000000000000000000001p+0)));
    try std.testing.expectEqual(0xc.123925c00603b209538b612fb7ap+0, lgamma(@as(f128, -0x3.fffffcp+0)));
    try std.testing.expectEqual(0x2.02c25d6cc86b656f154465a20502p+4, lgamma(@as(f128, -0x3.ffffffffffffep+0)));
    // try std.testing.expectEqual(0x2.7cc0cde728aa06126e9a48c7b426p+4, lgamma(@as(f128, -0x3.fffffffffffffffcp+0)));
    try std.testing.expectEqual(0x4.9c2e374e1b356ab0287f636cca54p+4, lgamma(@as(f128, -0x3.fffffffffffffffffffffffffffep+0)));
    // try std.testing.expectEqual(0x4.4e8c4cd1af6aa575002bc4da48b4p+4, lgamma(@as(f128, -0x3.ffffffffffffffffffffffffffp+0)));
    try std.testing.expectEqual(0xb.60c6fbb5695c876615d9b462c398p+0, lgamma(@as(f128, -0x4.000008p+0)));
    try std.testing.expectEqual(0x1.f7ab3bed4b4e64caf3157f5aaa98p+4, lgamma(@as(f128, -0x4.0000000000004p+0)));
    // try std.testing.expectEqual(0x2.71a9ac67ab8d0e7690cf5b78d22ap+4, lgamma(@as(f128, -0x4.0000000000000008p+0)));
    try std.testing.expectEqual(0x4.911715ce9e1873156be128334a14p+4, lgamma(@as(f128, -0x4.0000000000000000000000000004p+0)));
    // try std.testing.expectEqual(0x4.43752b52324dadda438d89a0c82cp+4, lgamma(@as(f128, -0x4.00000000000000000000000002p+0)));
    try std.testing.expectEqual(0x9.c4c2f5e938fb4f78265b70fa6d2p+0, lgamma(@as(f128, -0x4.fffff8p+0)));
    try std.testing.expectEqual(0x1.ddeaf9f55dc13e39495660933b62p+4, lgamma(@as(f128, -0x4.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x2.57e96a6fbdffdb0d2e026832e2f4p+4, lgamma(@as(f128, -0x4.fffffffffffffff8p+0)));
    try std.testing.expectEqual(0x4.7756d3d6b08b3faa6de9ade1e9ap+4, lgamma(@as(f128, -0x4.fffffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(0x4.29b4e95a44c07a6f45960f4f681cp+4, lgamma(@as(f128, -0x4.fffffffffffffffffffffffffep+0)));
    try std.testing.expectEqual(0x9.c4c2da9cf6f0fecaafe5d0803558p+0, lgamma(@as(f128, -0x5.000008p+0)));
    // try std.testing.expectEqual(0x1.ddeaf9f55dc130932851383c7d5bp+4, lgamma(@as(f128, -0x5.0000000000004p+0)));
    try std.testing.expectEqual(0x2.57e96a6fbdffdb0b793e478dd81cp+4, lgamma(@as(f128, -0x5.0000000000000008p+0)));
    try std.testing.expectEqual(0x4.7756d3d6b08b3faa6de9ade1e9ap+4, lgamma(@as(f128, -0x5.0000000000000000000000000004p+0)));
    try std.testing.expectEqual(0x4.29b4e95a44c07a6f45960f4f67bp+4, lgamma(@as(f128, -0x5.00000000000000000000000002p+0)));
    // try std.testing.expectEqual(0x7.fa12379bec516539476159244dd4p+0, lgamma(@as(f128, -0x5.fffff8p+0)));
    try std.testing.expectEqual(0x1.c13fedfb33a13cb1cd21372f157dp+4, lgamma(@as(f128, -0x5.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x2.3b3e5e7593dfd8db1c77e97967b6p+4, lgamma(@as(f128, -0x5.fffffffffffffff8p+0)));
    try std.testing.expectEqual(0x4.5aabc7dc866b3d784709d9d3190cp+4, lgamma(@as(f128, -0x5.fffffffffffffffffffffffffffcp+0)));
    // try std.testing.expectEqual(0x4.0d09dd601aa0783d1eb63b40979p+4, lgamma(@as(f128, -0x5.fffffffffffffffffffffffffep+0)));
    try std.testing.expectEqual(0x7.fa1219a4ff9c69e124ac82ef9dfp+0, lgamma(@as(f128, -0x6.000008p+0)));
    // try std.testing.expectEqual(0x1.c13fedfb33a12db656c6b9830221p+4, lgamma(@as(f128, -0x6.0000000000004p+0)));
    try std.testing.expectEqual(0x2.3b3e5e7593dfd8d93d091e29b234p+4, lgamma(@as(f128, -0x6.0000000000000008p+0)));
    try std.testing.expectEqual(0x4.5aabc7dc866b3d784709d9d3190cp+4, lgamma(@as(f128, -0x6.0000000000000000000000000004p+0)));
    try std.testing.expectEqual(0x4.0d09dd601aa0783d1eb63b409718p+4, lgamma(@as(f128, -0x6.00000000000000000000000002p+0)));
    try std.testing.expectEqual(0x6.07eb0ddd58f5bbb39faa2d8f8c7p+0, lgamma(@as(f128, -0x6.fffff8p+0)));
    try std.testing.expectEqual(0x1.a21d7b4d0146e5efa6dc800b47bp+4, lgamma(@as(f128, -0x6.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x2.1c1bebc761858186bf57c49ebe78p+4, lgamma(@as(f128, -0x6.fffffffffffffff8p+0)));
    try std.testing.expectEqual(0x4.3b89552e5410e623d7a0906626acp+4, lgamma(@as(f128, -0x6.fffffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(0x3.ede76ab1e84620e8af4cf1d3a532p+4, lgamma(@as(f128, -0x6.fffffffffffffffffffffffffep+0)));
    // try std.testing.expectEqual(0x6.07eaed9d47ae7736e9ad713a84fcp+0, lgamma(@as(f128, -0x7.000008p+0)));
    try std.testing.expectEqual(0x1.a21d7b4d0146d5cf9e38ddcceb2fp+4, lgamma(@as(f128, -0x7.0000000000004p+0)));
    // try std.testing.expectEqual(0x2.1c1bebc761858184bb56b02a76aep+4, lgamma(@as(f128, -0x7.0000000000000008p+0)));
    // try std.testing.expectEqual(0x4.3b89552e5410e623d7a0906626acp+4, lgamma(@as(f128, -0x7.0000000000000000000000000004p+0)));
    try std.testing.expectEqual(0x3.ede76ab1e84620e8af4cf1d3a4bp+4, lgamma(@as(f128, -0x7.00000000000000000000000002p+0)));
    try std.testing.expectEqual(0x3.f394c6f5e387ceb04254681d15ecp+0, lgamma(@as(f128, -0x7.fffff8p+0)));
    // try std.testing.expectEqual(0x1.80d816ce89efff9f7101ce5ec6f5p+4, lgamma(@as(f128, -0x7.ffffffffffffcp+0)));
    try std.testing.expectEqual(0x1.fad68748ea2e9ab6997d12f23dbbp+4, lgamma(@as(f128, -0x7.fffffffffffffff8p+0)));
    // try std.testing.expectEqual(0x4.1a43f0afdcb9ff53a1c5deb9a5ecp+4, lgamma(@as(f128, -0x7.fffffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(0x3.cca2063370ef3a18797240272478p+4, lgamma(@as(f128, -0x7.fffffffffffffffffffffffffep+0)));
    // try std.testing.expectEqual(0x3.42227b9df8fdfa1c5dea97787f46p+0, lgamma(@as(f128, -0x8.00001p+0)));
    try std.testing.expectEqual(0x1.75c0f54f0cd2ee54a76e1fc7c0b6p+4, lgamma(@as(f128, -0x8.0000000000008p+0)));
    // try std.testing.expectEqual(0x1.efbf65c96d11a318a6dd390a51cbp+4, lgamma(@as(f128, -0x8.000000000000001p+0)));
    try std.testing.expectEqual(0x4.0f2ccf305f9d07b8e527a38025acp+4, lgamma(@as(f128, -0x8.0000000000000000000000000008p+0)));
    // try std.testing.expectEqual(0x3.c18ae4b3f3d2427dbcd404eda36ap+4, lgamma(@as(f128, -0x8.00000000000000000000000004p+0)));
    // try std.testing.expectEqual(0x1.0fa5728f979e8bcff85a754cd032p+0, lgamma(@as(f128, -0x8.fffffp+0)));
    // try std.testing.expectEqual(0x1.52992059b2ccfc49726b162811fbp+4, lgamma(@as(f128, -0x8.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.cc9790d4130b8dee36cdf764b281p+4, lgamma(@as(f128, -0x8.fffffffffffffffp+0)));
    try std.testing.expectEqual(0x3.ec04fa3b0596f28a10a471d58508p+4, lgamma(@as(f128, -0x8.fffffffffffffffffffffffffff8p+0)));
    try std.testing.expectEqual(0x3.9e630fbe99cc2d4ee850d34303dcp+4, lgamma(@as(f128, -0x8.fffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(0x1.0fa52a813c2c749db1de5995236p+0, lgamma(@as(f128, -0x9.00001p+0)));
    // try std.testing.expectEqual(0x1.52992059b2ccd84244b20a8ee732p+4, lgamma(@as(f128, -0x9.0000000000008p+0)));
    try std.testing.expectEqual(0x1.cc9790d4130b8de9b5e840433f5cp+4, lgamma(@as(f128, -0x9.000000000000001p+0)));
    try std.testing.expectEqual(0x3.ec04fa3b0596f28a10a471d58506p+4, lgamma(@as(f128, -0x9.0000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x3.9e630fbe99cc2d4ee850d34302bcp+4, lgamma(@as(f128, -0x9.00000000000000000000000004p+0)));
    try std.testing.expectEqual(-0x1.3dd0c34d79694344018ee202113p+0, lgamma(@as(f128, -0x9.fffffp+0)));
    try std.testing.expectEqual(0x1.2dc1bce24822d21084a22d69fe18p+4, lgamma(@as(f128, -0x9.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.a7c02d5ca86162e895d1db736b66p+4, lgamma(@as(f128, -0x9.fffffffffffffffp+0)));
    // try std.testing.expectEqual(0x3.c72d96c39aecc784560ebc4aa454p+4, lgamma(@as(f128, -0x9.fffffffffffffffffffffffffff8p+0)));
    try std.testing.expectEqual(0x3.798bac472f2202492dbb1db8232ep+4, lgamma(@as(f128, -0x9.fffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(-0x1.3dd10e8f080e8da97df93de56ed2p+0, lgamma(@as(f128, -0xa.00001p+0)));
    // try std.testing.expectEqual(0x1.2dc1bce24822ac6fbd4f883739b5p+4, lgamma(@as(f128, -0xa.0000000000008p+0)));
    // try std.testing.expectEqual(0x1.a7c02d5ca86162e3e1b8f11ec50ep+4, lgamma(@as(f128, -0xa.000000000000001p+0)));
    // try std.testing.expectEqual(0x3.c72d96c39aecc784560ebc4aa45p+4, lgamma(@as(f128, -0xa.0000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x3.798bac472f2202492dbb1db82202p+4, lgamma(@as(f128, -0xa.00000000000000000000000004p+0)));
    try std.testing.expectEqual(-0x3.a3ad38c9033a659ac104c00477e4p+0, lgamma(@as(f128, -0xa.fffffp+0)));
    // try std.testing.expectEqual(0x1.0763f57349b43b5b3a7450b9687p+4, lgamma(@as(f128, -0xa.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.816265eda9f2cb79345e2d4e78a3p+4, lgamma(@as(f128, -0xa.fffffffffffffffp+0)));
    try std.testing.expectEqual(0x3.a0cfcf549c7e3014dd553cb15478p+4, lgamma(@as(f128, -0xa.fffffffffffffffffffffffffff8p+0)));
    try std.testing.expectEqual(0x3.532de4d830b36ad9b5019e1ed35ap+4, lgamma(@as(f128, -0xa.fffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(-0x3.a3ad86f34c0e3ba328367f78cabp+0, lgamma(@as(f128, -0xb.00001p+0)));
    try std.testing.expectEqual(0x1.0763f57349b41446160a65b52fbp+4, lgamma(@as(f128, -0xb.0000000000008p+0)));
    try std.testing.expectEqual(0x1.816265eda9f2cb7451b9a011181cp+4, lgamma(@as(f128, -0xb.000000000000001p+0)));
    try std.testing.expectEqual(0x3.a0cfcf549c7e3014dd553cb15476p+4, lgamma(@as(f128, -0xb.0000000000000000000000000008p+0)));
    // try std.testing.expectEqual(0x3.532de4d830b36ad9b5019e1ed22p+4, lgamma(@as(f128, -0xb.00000000000000000000000004p+0)));
    try std.testing.expectEqual(-0x6.1fd00f0e21b3c98569e28b729b24p+0, lgamma(@as(f128, -0xb.fffffp+0)));
    // try std.testing.expectEqual(0xd.fa1c7f9a277423901a0ec1bc24c8p+0, lgamma(@as(f128, -0xb.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.59a0387402b5d1ac6635735b7d26p+4, lgamma(@as(f128, -0xb.fffffffffffffffp+0)));
    // try std.testing.expectEqual(0x3.790da1daf5413647f9d72d6903a6p+4, lgamma(@as(f128, -0xb.fffffffffffffffffffffffffff8p+0)));
    // try std.testing.expectEqual(0x3.2b6bb75e8976710cd1838ed6828cp+4, lgamma(@as(f128, -0xb.fffffffffffffffffffffffffcp+0)));
    // try std.testing.expectEqual(-0x6.1fd05fe315324a387d5380a1660cp+0, lgamma(@as(f128, -0xc.00001p+0)));
    try std.testing.expectEqual(0xd.fa1c7f9a27719ce87e1abc234378p+0, lgamma(@as(f128, -0xc.0000000000008p+0)));
    // try std.testing.expectEqual(0x1.59a0387402b5d1a758e63b7371f5p+4, lgamma(@as(f128, -0xc.000000000000001p+0)));
    try std.testing.expectEqual(0x3.790da1daf5413647f9d72d6903a4p+4, lgamma(@as(f128, -0xc.0000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x3.2b6bb75e8976710cd1838ed6814ap+4, lgamma(@as(f128, -0xc.00000000000000000000000004p+0)));
    try std.testing.expectEqual(-0x8.b07093393f8bec5dcbeca94ad53p+0, lgamma(@as(f128, -0xc.fffffp+0)));
    try std.testing.expectEqual(0xb.697bfa33f5ea0d97e7debb452f2p+0, lgamma(@as(f128, -0xc.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.3096301d9f9d2faf6ceb107de666p+4, lgamma(@as(f128, -0xc.fffffffffffffffp+0)));
    try std.testing.expectEqual(0x3.500399849228944aecdb8f77bbacp+4, lgamma(@as(f128, -0xc.fffffffffffffffffffffffffff8p+0)));
    try std.testing.expectEqual(0x3.0261af08265dcf0fc487f0e53a96p+4, lgamma(@as(f128, -0xc.fffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(-0x8.b070e6845a6ce3384311f5033318p+0, lgamma(@as(f128, -0xd.00001p+0)));
    // try std.testing.expectEqual(0xb.697bfa33f5e7733f10d704713a18p+0, lgamma(@as(f128, -0xd.0000000000008p+0)));
    try std.testing.expectEqual(0x1.3096301d9f9d2faa3839626e78bep+4, lgamma(@as(f128, -0xd.000000000000001p+0)));
    try std.testing.expectEqual(0x3.500399849228944aecdb8f77bbaap+4, lgamma(@as(f128, -0xd.0000000000000000000000000008p+0)));
    // try std.testing.expectEqual(0x3.0261af08265dcf0fc487f0e5394ap+4, lgamma(@as(f128, -0xd.00000000000000000000000004p+0)));
    try std.testing.expectEqual(-0xb.5409d4efa4b70f8f3d8788779a88p+0, lgamma(@as(f128, -0xd.fffffp+0)));
    try std.testing.expectEqual(0x8.c5e2b758fe7527c9b9af95704e6p+0, lgamma(@as(f128, -0xd.ffffffffffff8p+0)));
    try std.testing.expectEqual(0x1.065c9beff025e0c0532cb069bceap+4, lgamma(@as(f128, -0xd.fffffffffffffffp+0)));
    try std.testing.expectEqual(0x3.25ca0556e2b1455bc0d40ad1490ap+4, lgamma(@as(f128, -0xd.fffffffffffffffffffffffffff8p+0)));
    try std.testing.expectEqual(0x2.d8281ada76e6802098806c3ec7fap+4, lgamma(@as(f128, -0xd.fffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(-0xb.540a2a83e42a4f8e47f4ba505p+0, lgamma(@as(f128, -0xe.00001p+0)));
    // try std.testing.expectEqual(0x8.c5e2b758fe727b27be159577c71p+0, lgamma(@as(f128, -0xe.0000000000008p+0)));
    try std.testing.expectEqual(0x1.065c9beff025e0baf9e8b935bcf9p+4, lgamma(@as(f128, -0xe.000000000000001p+0)));
    try std.testing.expectEqual(0x3.25ca0556e2b1455bc0d40ad14908p+4, lgamma(@as(f128, -0xe.0000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x2.d8281ada76e6802098806c3ec6a4p+4, lgamma(@as(f128, -0xe.00000000000000000000000004p+0)));
    // try std.testing.expectEqual(-0xe.094c9b083ca94d01fbdb43c57ae8p+0, lgamma(@as(f128, -0xe.fffffp+0)));
    try std.testing.expectEqual(0x6.109ff02f5571502bbea6eb8dca78p+0, lgamma(@as(f128, -0xe.ffffffffffff8p+0)));
    // try std.testing.expectEqual(0xd.b086f7d5595a2bdfc04ae541d32p+0, lgamma(@as(f128, -0xe.fffffffffffffffp+0)));
    try std.testing.expectEqual(0x2.fa75d8e448210759589af7aa9842p+4, lgamma(@as(f128, -0xe.fffffffffffffffffffffffffff8p+0)));
    // try std.testing.expectEqual(0x2.acd3ee67dc56421e304759181736p+4, lgamma(@as(f128, -0xe.fffffffffffffffffffffffffcp+0)));
    // try std.testing.expectEqual(-0xe.094cf2be9e3eaf232939b809f3p+0, lgamma(@as(f128, -0xf.00001p+0)));
    // try std.testing.expectEqual(0x6.109ff02f556e9278b1fbda843218p+0, lgamma(@as(f128, -0xf.0000000000008p+0)));
    // try std.testing.expectEqual(0xd.b086f7d5595a2b8809e94fdfb1e8p+0, lgamma(@as(f128, -0xf.000000000000001p+0)));
    // try std.testing.expectEqual(0x2.fa75d8e448210759589af7aa983ep+4, lgamma(@as(f128, -0xf.0000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x2.acd3ee67dc56421e3047591815d6p+4, lgamma(@as(f128, -0xf.00000000000000000000000004p+0)));
    // try std.testing.expectEqual(-0x1.0cf14f9e783e6b3b12314bccff56p+4, lgamma(@as(f128, -0xf.fffffp+0)));
    try std.testing.expectEqual(0x3.4ad790500e33717c97181d2dbaccp+0, lgamma(@as(f128, -0xf.ffffffffffff8p+0)));
    try std.testing.expectEqual(0xa.eabe97f6121c453198bc16e1c35p+0, lgamma(@as(f128, -0xf.fffffffffffffffp+0)));
    // try std.testing.expectEqual(0x2.ce1952e653ad28ee66220ac49744p+4, lgamma(@as(f128, -0xf.fffffffffffffffffffffffffff8p+0)));
    try std.testing.expectEqual(0x2.80776869e7e263b33dce6c32163cp+4, lgamma(@as(f128, -0xf.fffffffffffffffffffffffffcp+0)));
    try std.testing.expectEqual(-0x1.180879870e33e355b67293d3944bp+4, lgamma(@as(f128, -0x1.000002p+4)));
    // try std.testing.expectEqual(0x2.996578583c5fc3443a33d008884ep+0, lgamma(@as(f128, -0x1.0000000000001p+4)));
    try std.testing.expectEqual(0xa.394c7ffe404ccaff3d4603368d9p+0, lgamma(@as(f128, -0x1.0000000000000002p+4)));
    // try std.testing.expectEqual(0x2.c3023166d6903153a983cf8b1702p+4, lgamma(@as(f128, -0x1.0000000000000000000000000001p+4)));
    // try std.testing.expectEqual(0x2.756046ea6ac56c18813030f893e4p+4, lgamma(@as(f128, -0x1.000000000000000000000000008p+4)));
    // try std.testing.expectEqual(-0x1.455d45b618e1f038dddeea5dfff7p+4, lgamma(@as(f128, -0x1.0ffffep+4)));
    // try std.testing.expectEqual(-0x3.be7ffe71389cc26835a85ecbcda4p-4, lgamma(@as(f128, -0x1.0ffffffffffffp+4)));
    try std.testing.expectEqual(0x7.63ff07bef05d91d4a5f788c52c74p+0, lgamma(@as(f128, -0x1.0ffffffffffffffep+4)));
    try std.testing.expectEqual(0x2.95ad59e2e1913db5ab2497199ebp+4, lgamma(@as(f128, -0x1.0fffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0x2.480b6f6675c6787a82d0f8871e62p+4, lgamma(@as(f128, -0x1.0fffffffffffffffffffffffff8p+4)));
    // try std.testing.expectEqual(-0x1.455d51292150d8b93e426f65c468p+4, lgamma(@as(f128, -0x1.100002p+4)));
    // try std.testing.expectEqual(-0x3.be7ffe7138f85aabacec61e0bb4ep-4, lgamma(@as(f128, -0x1.1000000000001p+4)));
    // try std.testing.expectEqual(0x7.63ff07bef05d911d75709a3d2648p+0, lgamma(@as(f128, -0x1.1000000000000002p+4)));
    // try std.testing.expectEqual(0x2.95ad59e2e1913db5ab2497199eacp+4, lgamma(@as(f128, -0x1.1000000000000000000000000001p+4)));
    // try std.testing.expectEqual(0x2.480b6f6675c6787a82d0f8871b86p+4, lgamma(@as(f128, -0x1.100000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x1.739c3c0e7e3dc747f6c9173a7b13p+4, lgamma(@as(f128, -0x1.1ffffep+4)));
    try std.testing.expectEqual(-0x3.1fd7673485ba8a86b1e31b4b3ca6p+0, lgamma(@as(f128, -0x1.1ffffffffffffp+4)));
    // try std.testing.expectEqual(0x4.800fa0717e2cc53d5afd2c4a3a7cp+0, lgamma(@as(f128, -0x1.1ffffffffffffffep+4)));
    try std.testing.expectEqual(0x2.676e636e0a6e30ec1a032a357dcap+4, lgamma(@as(f128, -0x1.1fffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0x2.19cc78f19ea36bb0f1af8ba2fd82p+4, lgamma(@as(f128, -0x1.1fffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x1.739c47ba6a3ae8abe5a16e7d7a65p+4, lgamma(@as(f128, -0x1.200002p+4)));
    try std.testing.expectEqual(-0x3.1fd7673485c0607cb073cd43a7f2p+0, lgamma(@as(f128, -0x1.2000000000001p+4)));
    try std.testing.expectEqual(0x4.800fa0717e2cc4829c3d5a33fb6cp+0, lgamma(@as(f128, -0x1.2000000000000002p+4)));
    // try std.testing.expectEqual(0x2.676e636e0a6e30ec1a032a357dc4p+4, lgamma(@as(f128, -0x1.2000000000000000000000000001p+4)));
    // try std.testing.expectEqual(0x2.19cc78f19ea36bb0f1af8ba2fa98p+4, lgamma(@as(f128, -0x1.200000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x1.a2b8a7ff951d4cd8ff71bbc81688p+4, lgamma(@as(f128, -0x1.2ffffep+4)));
    // try std.testing.expectEqual(-0x6.119e27f51c200b4d7dd7ace1fa28p+0, lgamma(@as(f128, -0x1.2ffffffffffffp+4)));
    try std.testing.expectEqual(0x1.8e48dfb0e7c736fefad2b5a6035bp+0, lgamma(@as(f128, -0x1.2ffffffffffffffep+4)));
    try std.testing.expectEqual(0x2.3851f7620107d808190dfc0e98aap+4, lgamma(@as(f128, -0x1.2fffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0x1.eab00ce5953d12ccf0ba5d7c1868p+4, lgamma(@as(f128, -0x1.2fffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x1.a2b8b3e16627e7804ccde008eedcp+4, lgamma(@as(f128, -0x1.300002p+4)));
    try std.testing.expectEqual(-0x6.119e27f51c25fc36032500898dep+0, lgamma(@as(f128, -0x1.3000000000001p+4)));
    // try std.testing.expectEqual(0x1.8e48dfb0e7c73640ddc20bfb8e68p+0, lgamma(@as(f128, -0x1.3000000000000002p+4)));
    try std.testing.expectEqual(0x2.3851f7620107d808190dfc0e98a4p+4, lgamma(@as(f128, -0x1.3000000000000000000000000001p+4)));
    try std.testing.expectEqual(0x1.eab00ce5953d12ccf0ba5d7c157p+4, lgamma(@as(f128, -0x1.300000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x1.d2a72cdce34ac164fbae8c7684ddp+4, lgamma(@as(f128, -0x1.3ffffep+4)));
    try std.testing.expectEqual(-0x9.10867763989228882449ec5b3c48p+0, lgamma(@as(f128, -0x1.3ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.709f6fbd94aaf306ded2bd06724bp+0, lgamma(@as(f128, -0x1.3ffffffffffffffep+4)));
    // try std.testing.expectEqual(0x2.0863726b1940b567a1da0b4a37b6p+4, lgamma(@as(f128, -0x1.3fffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0x1.bac187eead75f02c79866cb7b77bp+4, lgamma(@as(f128, -0x1.3fffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x1.d2a738f1e7888f3f7c6994ba183ep+4, lgamma(@as(f128, -0x1.400002p+4)));
    // try std.testing.expectEqual(-0x9.108677639898330a4330d99c6998p+0, lgamma(@as(f128, -0x1.4000000000001p+4)));
    // try std.testing.expectEqual(-0x1.709f6fbd94aaf3c82f1699e41a71p+0, lgamma(@as(f128, -0x1.4000000000000002p+4)));
    try std.testing.expectEqual(0x2.0863726b1940b567a1da0b4a37bp+4, lgamma(@as(f128, -0x1.4000000000000000000000000001p+4)));
    try std.testing.expectEqual(0x1.bac187eead75f02c79866cb7b475p+4, lgamma(@as(f128, -0x1.400000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x2.035d89ed6121f85bdcd2763fe0bcp+4, lgamma(@as(f128, -0x1.4ffffep+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e683b14fbebdfbc5b28p+0, lgamma(@as(f128, -0x1.4ffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.7c05424b8a8111c2f3687fa4854p+0, lgamma(@as(f128, -0x1.4ffffffffffffffep+4)));
    // try std.testing.expectEqual(0x1.d7ad154239e3537bc82f2907f5p+4, lgamma(@as(f128, -0x1.4fffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0x1.8a0b2ac5ce188e409fdb8a7574cbp+4, lgamma(@as(f128, -0x1.4fffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x2.035d9633286bf6f969e3ff6bccfp+4, lgamma(@as(f128, -0x1.500002p+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e6e5df8a0eb2e83a0d8p+0, lgamma(@as(f128, -0x1.5000000000001p+4)));
    // try std.testing.expectEqual(-0x4.7c05424b8a8112874fdd1f8e5e28p+0, lgamma(@as(f128, -0x1.5000000000000002p+4)));
    try std.testing.expectEqual(0x1.d7ad154239e3537bc82f2907f4fap+4, lgamma(@as(f128, -0x1.5000000000000000000000000001p+4)));
    try std.testing.expectEqual(0x1.8a0b2ac5ce188e409fdb8a7571bap+4, lgamma(@as(f128, -0x1.500000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x2.34d272c496dc021c05680f598766p+4, lgamma(@as(f128, -0x1.5ffffep+4)));
    // try std.testing.expectEqual(-0xf.333ad8d94721201568ad5e5db99p+0, lgamma(@as(f128, -0x1.5ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x7.9353d133433a0264d487158bb564p+0, lgamma(@as(f128, -0x1.5ffffffffffffffep+4)));
    try std.testing.expectEqual(0x1.a6382c53be57c47192d76e3524e6p+4, lgamma(@as(f128, -0x1.5fffffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(0x1.589641d7528cff366a83cfa2a4b7p+4, lgamma(@as(f128, -0x1.5fffffffffffffffffffffffff8p+4)));
    // try std.testing.expectEqual(-0x2.34d27f38e9c8e973c1260ebe82ecp+4, lgamma(@as(f128, -0x1.600002p+4)));
    // try std.testing.expectEqual(-0xf.333ad8d947275a3edf210a3c4518p+0, lgamma(@as(f128, -0x1.6000000000001p+4)));
    // try std.testing.expectEqual(-0x7.9353d133433a032c19b5e4013138p+0, lgamma(@as(f128, -0x1.6000000000000002p+4)));
    try std.testing.expectEqual(0x1.a6382c53be57c47192d76e3524ep+4, lgamma(@as(f128, -0x1.6000000000000000000000000001p+4)));
    try std.testing.expectEqual(0x1.589641d7528cff366a83cfa2a19ap+4, lgamma(@as(f128, -0x1.600000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x2.66fd6ea9f77b79a6b027c2a9dfa2p+4, lgamma(@as(f128, -0x1.6ffffep+4)));
    // try std.testing.expectEqual(-0x1.255ea98937d9f1616b540f71866cp+4, lgamma(@as(f128, -0x1.6ffffffffffffp+4)));
    try std.testing.expectEqual(-0xa.b60390ed79b80386211aae4662d8p+0, lgamma(@as(f128, -0x1.6ffffffffffffffep+4)));
    try std.testing.expectEqual(0x1.740d30581aefe45f67cb6c506eeep+4, lgamma(@as(f128, -0x1.6fffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0x1.266b45dbaf251f243f77cdbdeec4p+4, lgamma(@as(f128, -0x1.6fffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x2.66fd7b4acff91314aecad54bcbe4p+4, lgamma(@as(f128, -0x1.700002p+4)));
    try std.testing.expectEqual(-0x1.255ea98937da56682f40dae18568p+4, lgamma(@as(f128, -0x1.7000000000001p+4)));
    try std.testing.expectEqual(-0xa.b60390ed79b804502ea287dd42d8p+0, lgamma(@as(f128, -0x1.7000000000000002p+4)));
    // try std.testing.expectEqual(0x1.740d30581aefe45f67cb6c506ee7p+4, lgamma(@as(f128, -0x1.7000000000000000000000000001p+4)));
    try std.testing.expectEqual(0x1.266b45dbaf251f243f77cdbdeb9cp+4, lgamma(@as(f128, -0x1.700000000000000000000000008p+4)));
    // try std.testing.expectEqual(-0x2.99d6bd8dc68007801753da9a4214p+4, lgamma(@as(f128, -0x1.7ffffep+4)));
    // try std.testing.expectEqual(-0x1.5837f8825c33e21e60c5af48acdp+4, lgamma(@as(f128, -0x1.7ffffffffffffp+4)));
    try std.testing.expectEqual(-0xd.e398807fbf5719fecd8a010e1e98p+0, lgamma(@as(f128, -0x1.7ffffffffffffffep+4)));
    // try std.testing.expectEqual(0x1.4133e15ef695f2f7c7af21ce9dddp+4, lgamma(@as(f128, -0x1.7fffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0xf.391f6e28acb2dbc9f5b833c1db88p+0, lgamma(@as(f128, -0x1.7fffffffffffffffffffffffff8p+4)));
    // try std.testing.expectEqual(-0x2.99d6ca5949a84b98c0bae097d5dap+4, lgamma(@as(f128, -0x1.800002p+4)));
    try std.testing.expectEqual(-0x1.5837f8825c34487a7a07d00e012p+4, lgamma(@as(f128, -0x1.8000000000001p+4)));
    try std.testing.expectEqual(-0xd.e398807fbf571acb85bc854fa94p+0, lgamma(@as(f128, -0x1.8000000000000002p+4)));
    try std.testing.expectEqual(0x1.4133e15ef695f2f7c7af21ce9dd6p+4, lgamma(@as(f128, -0x1.8000000000000000000000000001p+4)));
    try std.testing.expectEqual(0xf.391f6e28acb2dbc9f5b833c1a858p+0, lgamma(@as(f128, -0x1.800000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x2.cd57416926b9198c8d473083f362p+4, lgamma(@as(f128, -0x1.8ffffep+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374e485085aa667ac9ep+4, lgamma(@as(f128, -0x1.8ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.11ba0bf7d70fd875d44cb36bf4cp+4, lgamma(@as(f128, -0x1.8ffffffffffffffep+4)));
    try std.testing.expectEqual(0x1.0db35d6f1b7b8c21cbc02d2bdcf1p+4, lgamma(@as(f128, -0x1.8fffffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(0xc.01172f2afb0c6e6a36c8e995cd28p+0, lgamma(@as(f128, -0x1.8fffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x2.cd574e5d9fa3ed015fba57b06442p+4, lgamma(@as(f128, -0x1.900002p+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374eaff44d01022165dfp+4, lgamma(@as(f128, -0x1.9000000000001p+4)));
    try std.testing.expectEqual(-0x1.11ba0bf7d70fd882c8c59e3f6994p+4, lgamma(@as(f128, -0x1.9000000000000002p+4)));
    try std.testing.expectEqual(0x1.0db35d6f1b7b8c21cbc02d2bdcebp+4, lgamma(@as(f128, -0x1.9000000000000000000000000001p+4)));
    try std.testing.expectEqual(0xc.01172f2afb0c6e6a36c8e9959958p+0, lgamma(@as(f128, -0x1.900000000000000000000000008p+4)));
    // try std.testing.expectEqual(-0x3.01786b2b55b39354d0060d9af742p+4, lgamma(@as(f128, -0x1.9ffffep+4)));
    try std.testing.expectEqual(-0x1.bfd9a6481783e14ac56ba21bb97ap+4, lgamma(@as(f128, -0x1.9ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.45db35cdb745720d8a3551830bbfp+4, lgamma(@as(f128, -0x1.9ffffffffffffffep+4)));
    try std.testing.expectEqual(0xd.99233993b45f28a0226540114b7p+0, lgamma(@as(f128, -0x1.9fffffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(0x8.bf0491ccf7b2d4ed9d2b56e949dp+0, lgamma(@as(f128, -0x1.9fffffffffffffffffffffffff8p+4)));
    // try std.testing.expectEqual(-0x3.0178784731148e2c18b47a300152p+4, lgamma(@as(f128, -0x1.a00002p+4)));
    try std.testing.expectEqual(-0x1.bfd9a64817844a29a07378d606b4p+4, lgamma(@as(f128, -0x1.a000000000001p+4)));
    try std.testing.expectEqual(-0x1.45db35cdb745721aa610b27de309p+4, lgamma(@as(f128, -0x1.a000000000000002p+4)));
    // try std.testing.expectEqual(0xd.99233993b45f28a0226540114b08p+0, lgamma(@as(f128, -0x1.a000000000000000000000000001p+4)));
    try std.testing.expectEqual(0x8.bf0491ccf7b2d4ed9d2b56e9156p+0, lgamma(@as(f128, -0x1.a00000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x3.36342a886637ea3d1ee94bbf39f4p+4, lgamma(@as(f128, -0x1.affffep+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8d0079500a7f922f3p+4, lgamma(@as(f128, -0x1.affffffffffffp+4)));
    try std.testing.expectEqual(-0x1.7a96f53dbe4e91d3b60397455b8bp+4, lgamma(@as(f128, -0x1.affffffffffffffep+4)));
    try std.testing.expectEqual(0xa.4d67429343cd2c3c361898123bcp+0, lgamma(@as(f128, -0x1.afffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0x5.73489acc8720d889b0deaeea3a6cp+0, lgamma(@as(f128, -0x1.afffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x3.363437ca2ea26056c67a1202c958p+4, lgamma(@as(f128, -0x1.b00002p+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8d6a87935e305f72eep+4, lgamma(@as(f128, -0x1.b000000000001p+4)));
    // try std.testing.expectEqual(-0x1.7a96f53dbe4e91e0f7cc01bb7533p+4, lgamma(@as(f128, -0x1.b000000000000002p+4)));
    try std.testing.expectEqual(0xa.4d67429343cd2c3c361898123b58p+0, lgamma(@as(f128, -0x1.b000000000000000000000000001p+4)));
    try std.testing.expectEqual(0x5.73489acc8720d889b0deaeea0564p+0, lgamma(@as(f128, -0x1.b00000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x3.6b84e02349a7940af2a134eb8688p+4, lgamma(@as(f128, -0x1.bffffep+4)));
    try std.testing.expectEqual(-0x2.29e61b654b214670ef8bad28fd7cp+4, lgamma(@as(f128, -0x1.bffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.afe7aaeaeae2d85d8c6032930547p+4, lgamma(@as(f128, -0x1.bffffffffffffffep+4)));
    try std.testing.expectEqual(0x6.f85be7c07a88c39dabbc9a130dcp+0, lgamma(@as(f128, -0x1.bfffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(0x2.1e3d3ff9bddc6feb2682b0eb0cb4p+0, lgamma(@as(f128, -0x1.bfffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x3.6b84ed89a45b2eb6e36679911b66p+4, lgamma(@as(f128, -0x1.c00002p+4)));
    try std.testing.expectEqual(-0x2.29e61b654b21b1a3c52882888a5ep+4, lgamma(@as(f128, -0x1.c000000000001p+4)));
    try std.testing.expectEqual(-0x1.afe7aaeaeae2d86af2bae62db138p+4, lgamma(@as(f128, -0x1.c000000000000002p+4)));
    try std.testing.expectEqual(0x6.f85be7c07a88c39dabbc9a130d54p+0, lgamma(@as(f128, -0x1.c000000000000000000000000001p+4)));
    try std.testing.expectEqual(0x2.1e3d3ff9bddc6feb2682b0ead71ap+0, lgamma(@as(f128, -0x1.c00000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x3.a16551a93dea66ada032f329cee6p+4, lgamma(@as(f128, -0x1.cffffep+4)));
    try std.testing.expectEqual(-0x2.5fc68cfce71d835e7f01e235d532p+4, lgamma(@as(f128, -0x1.cffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.e5c81c8286df15d847f9b7129f35p+4, lgamma(@as(f128, -0x1.cffffffffffffffep+4)));
    try std.testing.expectEqual(0x3.9a54ce46bac4ebf0d7a8bc07c728p+0, lgamma(@as(f128, -0x1.cfffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.3fc9d98001e767c1ad912d2039ap+0, lgamma(@as(f128, -0x1.cfffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x3.a1655f32e810c38e8832afeceb82p+4, lgamma(@as(f128, -0x1.d00002p+4)));
    try std.testing.expectEqual(-0x2.5fc68cfce71defabd034c93d1b76p+4, lgamma(@as(f128, -0x1.d000000000001p+4)));
    try std.testing.expectEqual(-0x1.e5c81c8286df15e5d1a3dd6f801dp+4, lgamma(@as(f128, -0x1.d000000000000002p+4)));
    try std.testing.expectEqual(0x3.9a54ce46bac4ebf0d7a8bc07c6bcp+0, lgamma(@as(f128, -0x1.d000000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x1.3fc9d98001e767c1ad912d206fc7p+0, lgamma(@as(f128, -0x1.d00000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x3.d7d09f8a4486821f88b66a182d2cp+4, lgamma(@as(f128, -0x1.dffffep+4)));
    // try std.testing.expectEqual(-0x2.9631daeefecab8731b50a80d7dbp+4, lgamma(@as(f128, -0x1.dffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.1c336a749e8c4b755bbff461bf2cp+4, lgamma(@as(f128, -0x1.dffffffffffffffep+4)));
    try std.testing.expectEqual(0x3.39fef253ff1921e8a33d604b6a06p-4, lgamma(@as(f128, -0x1.dfffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.a67eb8a17cbac193fb06132349e4p+0, lgamma(@as(f128, -0x1.dfffffffffffffffffffffffff8p+4)));
    // try std.testing.expectEqual(-0x3.d7d0ad3610cf012292e53b0205f2p+4, lgamma(@as(f128, -0x1.e00002p+4)));
    // try std.testing.expectEqual(-0x2.9631daeefecb25d17d94a025d504p+4, lgamma(@as(f128, -0x1.e000000000001p+4)));
    // try std.testing.expectEqual(-0x2.1c336a749e8c4b83078c3ce0c236p+4, lgamma(@as(f128, -0x1.e000000000000002p+4)));
    // try std.testing.expectEqual(0x3.39fef253ff1921e8a33d604b633p-4, lgamma(@as(f128, -0x1.e000000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x4.a67eb8a17cbac193fb0613238094p+0, lgamma(@as(f128, -0x1.e00000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x4.0ec23c0ae2bc253963f0c770efd4p+4, lgamma(@as(f128, -0x1.effffep+4)));
    try std.testing.expectEqual(-0x2.cd23778021216bd128a5aa6dd404p+4, lgamma(@as(f128, -0x1.effffffffffffp+4)));
    // try std.testing.expectEqual(-0x2.53250705c0e2ff57799917ca5794p+4, lgamma(@as(f128, -0x1.effffffffffffffep+4)));
    // try std.testing.expectEqual(-0x3.3b79d9ece579ac045ba07108f0dap+0, lgamma(@as(f128, -0x1.efffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x8.159881b3a225ffb6e0da5a30f118p+0, lgamma(@as(f128, -0x1.efffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x4.0ec249d7b746b4c08f3395f5c198p+4, lgamma(@as(f128, -0x1.f00002p+4)));
    // try std.testing.expectEqual(-0x2.cd2377802121da37ccfa26a7339cp+4, lgamma(@as(f128, -0x1.f000000000001p+4)));
    try std.testing.expectEqual(-0x2.53250705c0e2ff65466da259decp+4, lgamma(@as(f128, -0x1.f000000000000002p+4)));
    try std.testing.expectEqual(-0x3.3b79d9ece579ac045ba07108f148p+0, lgamma(@as(f128, -0x1.f000000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x8.159881b3a225ffb6e0da5a31285p+0, lgamma(@as(f128, -0x1.f00000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x4.4635e378544cf33f13029a3b17b8p+4, lgamma(@as(f128, -0x1.fffffep+4)));
    try std.testing.expectEqual(-0x3.04971efd92b24156d7bcd28d553ep+4, lgamma(@as(f128, -0x1.fffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.8a98ae833273d55d18b03fe9d8dp+4, lgamma(@as(f128, -0x1.fffffffffffffffep+4)));
    try std.testing.expectEqual(-0x6.b2b451c3fe870c5f4d12f3010498p+0, lgamma(@as(f128, -0x1.ffffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0xb.8cd2f98abb336011d24cdc290498p+0, lgamma(@as(f128, -0x1.ffffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x4.514d19db0f00e277f11efebfce18p+4, lgamma(@as(f128, -0x2.000004p+4)));
    try std.testing.expectEqual(-0x3.0fae407d0fcfe00b8ad9c81c96a8p+4, lgamma(@as(f128, -0x2.0000000000002p+4)));
    try std.testing.expectEqual(-0x2.95afd002af90cd0cb88d4afaa3dp+4, lgamma(@as(f128, -0x2.0000000000000004p+4)));
    // try std.testing.expectEqual(-0x7.642669bbd056860b16f6a699093p+0, lgamma(@as(f128, -0x2.0000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0xc.3e4511828d02d9bd9c308fc15c18p+0, lgamma(@as(f128, -0x2.00000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x4.893eafcc099b56dd588d9a421c54p+4, lgamma(@as(f128, -0x2.0ffffcp+4)));
    try std.testing.expectEqual(-0x3.479ff266bb40a24ce71c64377f62p+4, lgamma(@as(f128, -0x2.0fffffffffffep+4)));
    try std.testing.expectEqual(-0x2.cda181ec5b026ef7a2d78c59b988p+4, lgamma(@as(f128, -0x2.0ffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0xa.e34188568770a67946a82d830b3p+0, lgamma(@as(f128, -0x2.0ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0xf.bd60301d441cfa2bcbe216aaef18p+0, lgamma(@as(f128, -0x2.0fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.893ecbe3c23456e3e9a55f8403fcp+4, lgamma(@as(f128, -0x2.100004p+4)));
    try std.testing.expectEqual(-0x3.479ff266bb41830aabe4646c2f0cp+4, lgamma(@as(f128, -0x2.1000000000002p+4)));
    // try std.testing.expectEqual(-0x2.cda181ec5b026f13ba902559c01ep+4, lgamma(@as(f128, -0x2.1000000000000004p+4)));
    try std.testing.expectEqual(-0xa.e34188568770a67946a82d830c18p+0, lgamma(@as(f128, -0x2.1000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0xf.bd60301d441cfa2bcbe216ab5f78p+0, lgamma(@as(f128, -0x2.10000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x4.c1aaa8b15d99079f60785038ed2cp+4, lgamma(@as(f128, -0x2.1ffffcp+4)));
    try std.testing.expectEqual(-0x3.800beb6a2d5c8c94b128e6f187p+4, lgamma(@as(f128, -0x2.1fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.060d7aefcd1e5a303fb6e1e693fep+4, lgamma(@as(f128, -0x2.1ffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0xe.6a01188da92f5a04f67f68329488p+0, lgamma(@as(f128, -0x2.1ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.3441fc05465dbadb77bb9515a78p+4, lgamma(@as(f128, -0x2.1fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.c1aac505526e43e22e13764fb458p+4, lgamma(@as(f128, -0x2.200004p+4)));
    try std.testing.expectEqual(-0x3.800beb6a2d5d6f3457d2c908188cp+4, lgamma(@as(f128, -0x2.2000000000002p+4)));
    // try std.testing.expectEqual(-0x3.060d7aefcd1e5a4c93abb722d6dp+4, lgamma(@as(f128, -0x2.2000000000000004p+4)));
    try std.testing.expectEqual(-0xe.6a01188da92f5a04f67f6832957p+0, lgamma(@as(f128, -0x2.2000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x1.3441fc05465dbadb77bb9515ae95p+4, lgamma(@as(f128, -0x2.20000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x4.fa8d5d3a3bac5a5d21bafa1a797cp+4, lgamma(@as(f128, -0x2.2ffffcp+4)));
    try std.testing.expectEqual(-0x3.b8eea0104d44166a0fe8c0a138e6p+4, lgamma(@as(f128, -0x2.2fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.3ef02f95ed05e4ef8fd5d187a502p+4, lgamma(@as(f128, -0x2.2ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x1.1f82c62efa7a805fbcc8ba417c2p+4, lgamma(@as(f128, -0x2.2ffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x1.6d24b0ab6645459ae51c58d3fa5p+4, lgamma(@as(f128, -0x2.2fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.fa8d79c8b429d12397d1db8b28d8p+4, lgamma(@as(f128, -0x2.300004p+4)));
    try std.testing.expectEqual(-0x3.b8eea0104d44faddd3d476d50c46p+4, lgamma(@as(f128, -0x2.3000000000002p+4)));
    // try std.testing.expectEqual(-0x3.3ef02f95ed05e50c1e4e4efe6b7cp+4, lgamma(@as(f128, -0x2.3000000000000004p+4)));
    try std.testing.expectEqual(-0x1.1f82c62efa7a805fbcc8ba417c2ep+4, lgamma(@as(f128, -0x2.3000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x1.6d24b0ab6645459ae51c58d40173p+4, lgamma(@as(f128, -0x2.30000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x5.33e375121e252906f743623074d8p+4, lgamma(@as(f128, -0x2.3ffffcp+4)));
    try std.testing.expectEqual(-0x3.f244b804a18419eacf6f8530a122p+4, lgamma(@as(f128, -0x2.3fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.7846478a4145e953c123b288d46p+4, lgamma(@as(f128, -0x2.3ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x1.58d8de234eba84c40a88625f1d46p+4, lgamma(@as(f128, -0x2.3ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.a67ac89fba8549ff32dc00f19b6ep+4, lgamma(@as(f128, -0x2.3fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x5.33e391d97a30d8b0fbcf15dc5f14p+4, lgamma(@as(f128, -0x2.400004p+4)));
    try std.testing.expectEqual(-0x3.f244b804a1850025afcd0280e64ap+4, lgamma(@as(f128, -0x2.4000000000002p+4)));
    // try std.testing.expectEqual(-0x3.7846478a4145e970887fbe387e6ap+4, lgamma(@as(f128, -0x2.4000000000000004p+4)));
    try std.testing.expectEqual(-0x1.58d8de234eba84c40a88625f1d54p+4, lgamma(@as(f128, -0x2.4000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x1.a67ac89fba8549ff32dc00f1a2ap+4, lgamma(@as(f128, -0x2.40000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x5.6da9c6d2e6bb76c4e0dab44bd54cp+4, lgamma(@as(f128, -0x2.4ffffcp+4)));
    try std.testing.expectEqual(-0x4.2c0b09e11713937c977995628abcp+4, lgamma(@as(f128, -0x2.4fffffffffffep+4)));
    // try std.testing.expectEqual(-0x3.b20c9966b6d563c2d5496fb3d2c2p+4, lgamma(@as(f128, -0x2.4ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x1.929f2fffc449ff333a5b189edd61p+4, lgamma(@as(f128, -0x2.4ffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x1.e0411a7c3014c46e62aeb7315b83p+4, lgamma(@as(f128, -0x2.4fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x5.6da9e3d19cb94ff25b3cbe8e617cp+4, lgamma(@as(f128, -0x2.500004p+4)));
    try std.testing.expectEqual(-0x4.2c0b09e117147b7247685ece7cdcp+4, lgamma(@as(f128, -0x2.5000000000002p+4)));
    try std.testing.expectEqual(-0x3.b20c9966b6d563dfd3ff6d8d004p+4, lgamma(@as(f128, -0x2.5000000000000004p+4)));
    try std.testing.expectEqual(-0x1.929f2fffc449ff333a5b189edd6fp+4, lgamma(@as(f128, -0x2.5000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x1.e0411a7c3014c46e62aeb73162c2p+4, lgamma(@as(f128, -0x2.50000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x5.a7dd54437ab7f3f0a6219412f1p+4, lgamma(@as(f128, -0x2.5ffffcp+4)));
    // try std.testing.expectEqual(-0x4.663e976c9d96e323c0d719b576d4p+4, lgamma(@as(f128, -0x2.5fffffffffffep+4)));
    try std.testing.expectEqual(-0x3.ec4026f23d58b44177ea52579672p+4, lgamma(@as(f128, -0x2.5ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x1.ccd2bd8b4acd4fb1f7ee81ff42c1p+4, lgamma(@as(f128, -0x2.5ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.1a74a807b69814ed20422091c0dcp+4, lgamma(@as(f128, -0x2.5fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x5.a7dd717815c346617f076b535634p+4, lgamma(@as(f128, -0x2.600004p+4)));
    try std.testing.expectEqual(-0x4.663e976c9d97ccc89931ad3c5b78p+4, lgamma(@as(f128, -0x2.6000000000002p+4)));
    try std.testing.expectEqual(-0x3.ec4026f23d58b45eac855daa074ep+4, lgamma(@as(f128, -0x2.6000000000000004p+4)));
    // try std.testing.expectEqual(-0x1.ccd2bd8b4acd4fb1f7ee81ff42dp+4, lgamma(@as(f128, -0x2.6000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x2.1a74a807b69814ed20422091c82ap+4, lgamma(@as(f128, -0x2.60000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x5.e27b46fa492f70b847dc589bbab8p+4, lgamma(@as(f128, -0x2.6ffffcp+4)));
    try std.testing.expectEqual(-0x4.a0dc8a3dadb28ee62af37e6eee48p+4, lgamma(@as(f128, -0x2.6fffffffffffep+4)));
    try std.testing.expectEqual(-0x4.26de19c34d7460d5d4e5e503ed1cp+4, lgamma(@as(f128, -0x2.6ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x2.0770b05c5ae8fc466f2bb8c5db1p+4, lgamma(@as(f128, -0x2.6ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.55129ad8c6b3c181977f57585924p+4, lgamma(@as(f128, -0x2.6fffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0x5.e27b64636782f7ac6925d6926bacp+4, lgamma(@as(f128, -0x2.700004p+4)));
    try std.testing.expectEqual(-0x4.a0dc8a3dadb37a2f1d8fb6101494p+4, lgamma(@as(f128, -0x2.7000000000002p+4)));
    // try std.testing.expectEqual(-0x4.26de19c34d7460f33e04388ae14p+4, lgamma(@as(f128, -0x2.7000000000000004p+4)));
    // try std.testing.expectEqual(-0x2.0770b05c5ae8fc466f2bb8c5db1ep+4, lgamma(@as(f128, -0x2.7000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x2.55129ad8c6b3c181977f5758607ep+4, lgamma(@as(f128, -0x2.70000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x6.1d80ed571479dcdf00b76483a94cp+4, lgamma(@as(f128, -0x2.7ffffcp+4)));
    try std.testing.expectEqual(-0x4.dbe230b41296a85491f8dda002acp+4, lgamma(@as(f128, -0x2.7fffffffffffep+4)));
    try std.testing.expectEqual(-0x4.61e3c039b2587b10ef1e776834b8p+4, lgamma(@as(f128, -0x2.7ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x2.427656d2bfcd1681a2fde4c3bc42p+4, lgamma(@as(f128, -0x2.7ffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x2.9018414f2b97dbbccb5183563a52p+4, lgamma(@as(f128, -0x2.7fffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0x6.1d810af366009706555fc67d155p+4, lgamma(@as(f128, -0x2.800004p+4)));
    try std.testing.expectEqual(-0x4.dbe230b4129795371e2eaedac29p+4, lgamma(@as(f128, -0x2.8000000000002p+4)));
    try std.testing.expectEqual(-0x4.61e3c039b2587b2e8b6ffe225c1p+4, lgamma(@as(f128, -0x2.8000000000000004p+4)));
    try std.testing.expectEqual(-0x2.427656d2bfcd1681a2fde4c3bc52p+4, lgamma(@as(f128, -0x2.8000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x2.9018414f2b97dbbccb51835641b8p+4, lgamma(@as(f128, -0x2.80000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x6.58ebb7c93810d51c01c1079e2b44p+4, lgamma(@as(f128, -0x2.8ffffcp+4)));
    // try std.testing.expectEqual(-0x5.174cfb3f2fef42e41dac5665aebp+4, lgamma(@as(f128, -0x2.8fffffffffffep+4)));
    try std.testing.expectEqual(-0x4.9d4e8ac4cfb116682fe4ab7f0c74p+4, lgamma(@as(f128, -0x2.8ffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x2.7de1215ddd25b1d8fcbdda6a301ap+4, lgamma(@as(f128, -0x2.8ffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x2.cb830bda48f07714251178fcae22p+4, lgamma(@as(f128, -0x2.8fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x6.58ebd5977d1aae7b88857eec360cp+4, lgamma(@as(f128, -0x2.900004p+4)));
    // try std.testing.expectEqual(-0x5.174cfb3f2ff0315645fb2161fe3p+4, lgamma(@as(f128, -0x2.9000000000002p+4)));
    try std.testing.expectEqual(-0x4.9d4e8ac4cfb11685fe29b5586bfcp+4, lgamma(@as(f128, -0x2.9000000000000004p+4)));
    try std.testing.expectEqual(-0x2.7de1215ddd25b1d8fcbdda6a3028p+4, lgamma(@as(f128, -0x2.9000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x2.cb830bda48f07714251178fcb594p+4, lgamma(@as(f128, -0x2.90000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x6.94b93659330503ad9f832ca10764p+4, lgamma(@as(f128, -0x2.9ffffcp+4)));
    // try std.testing.expectEqual(-0x5.531a79e78c699ba7a7c4b0d540dcp+4, lgamma(@as(f128, -0x2.9fffffffffffep+4)));
    // try std.testing.expectEqual(-0x4.d91c096d2c2b6feeadcc42e26de4p+4, lgamma(@as(f128, -0x2.9ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x2.b9aea00639a00b5f9306f7e5f30ep+4, lgamma(@as(f128, -0x2.9ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.07508a82a56ad09abb5a9678711p+4, lgamma(@as(f128, -0x2.9fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x6.94b954583b1b0dd0329e24d76afcp+4, lgamma(@as(f128, -0x2.a00004p+4)));
    // try std.testing.expectEqual(-0x5.531a79e78c6a8b9fe87501e9f1e4p+4, lgamma(@as(f128, -0x2.a000000000002p+4)));
    try std.testing.expectEqual(-0x4.d91c096d2c2b700cacd458ec9078p+4, lgamma(@as(f128, -0x2.a000000000000004p+4)));
    // try std.testing.expectEqual(-0x2.b9aea00639a00b5f9306f7e5f31cp+4, lgamma(@as(f128, -0x2.a000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x3.07508a82a56ad09abb5a9678789p+4, lgamma(@as(f128, -0x2.a0000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x6.d0e7166d8c7dd2a47e015b47ae7p+4, lgamma(@as(f128, -0x2.affffcp+4)));
    try std.testing.expectEqual(-0x5.8f485a13b641bd15df1413add478p+4, lgamma(@as(f128, -0x2.afffffffffffep+4)));
    try std.testing.expectEqual(-0x5.1549e9995603921b50455261b42p+4, lgamma(@as(f128, -0x2.affffffffffffffcp+4)));
    try std.testing.expectEqual(-0x2.f5dc803263782d8c4d5066a6b65p+4, lgamma(@as(f128, -0x2.affffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.437e6aaecf42f2c775a40539344cp+4, lgamma(@as(f128, -0x2.afffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x6.d0e7349c35525fc11d27ae73b1fp+4, lgamma(@as(f128, -0x2.b00004p+4)));
    try std.testing.expectEqual(-0x5.8f485a13b642ae8b25b87c92e4cp+4, lgamma(@as(f128, -0x2.b000000000002p+4)));
    try std.testing.expectEqual(-0x5.1549e999560392397eee26eed0cp+4, lgamma(@as(f128, -0x2.b000000000000004p+4)));
    try std.testing.expectEqual(-0x2.f5dc803263782d8c4d5066a6b66p+4, lgamma(@as(f128, -0x2.b000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x3.437e6aaecf42f2c775a405393bd8p+4, lgamma(@as(f128, -0x2.b0000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x7.0d7320c43f54d3ff63352f9ad55cp+4, lgamma(@as(f128, -0x2.bffffcp+4)));
    // try std.testing.expectEqual(-0x5.cbd46481aeea4300a27e66d16aap+4, lgamma(@as(f128, -0x2.bfffffffffffep+4)));
    try std.testing.expectEqual(-0x5.51d5f4074eac18c02af576f9a76p+4, lgamma(@as(f128, -0x2.bffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x3.32688aa05c20b4313f465cb306aap+4, lgamma(@as(f128, -0x2.bffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.800a751cc7eb796c6799fb4584ap+4, lgamma(@as(f128, -0x2.bfffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x7.0d733f2173cc49d63107f8ffe82cp+4, lgamma(@as(f128, -0x2.c00004p+4)));
    try std.testing.expectEqual(-0x5.cbd46481aeeb35ea463a1587ef44p+4, lgamma(@as(f128, -0x2.c000000000002p+4)));
    try std.testing.expectEqual(-0x5.51d5f4074eac18de8829ee6f7e3p+4, lgamma(@as(f128, -0x2.c000000000000004p+4)));
    try std.testing.expectEqual(-0x3.32688aa05c20b4313f465cb306b8p+4, lgamma(@as(f128, -0x2.c000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x3.800a751cc7eb796c6799fb458c36p+4, lgamma(@as(f128, -0x2.c0000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x7.4a5b379ac57bf5a943c938e5a458p+4, lgamma(@as(f128, -0x2.cffffcp+4)));
    // try std.testing.expectEqual(-0x6.08bc7b6ef67d8ae469985cc20bp+4, lgamma(@as(f128, -0x2.cfffffffffffep+4)));
    try std.testing.expectEqual(-0x5.8ebe0af4963f6159e6aeb6dee71p+4, lgamma(@as(f128, -0x2.cffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x3.6f50a18da3b3fccb11c108af07c6p+4, lgamma(@as(f128, -0x2.cffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x3.bcf28c0a0f7ec2063a14a74185b6p+4, lgamma(@as(f128, -0x2.cfffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0x7.4a5b56257ccb9902e9e8349618e4p+4, lgamma(@as(f128, -0x2.d00004p+4)));
    try std.testing.expectEqual(-0x6.08bc7b6ef67e7f3a2415778f5114p+4, lgamma(@as(f128, -0x2.d000000000002p+4)));
    // try std.testing.expectEqual(-0x5.8ebe0af4963f61787166068240bcp+4, lgamma(@as(f128, -0x2.d000000000000004p+4)));
    try std.testing.expectEqual(-0x3.6f50a18da3b3fccb11c108af07d6p+4, lgamma(@as(f128, -0x2.d000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x3.bcf28c0a0f7ec2063a14a7418d5ap+4, lgamma(@as(f128, -0x2.d0000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x7.879d54ffa33864ceab27276f7cd4p+4, lgamma(@as(f128, -0x2.dffffcp+4)));
    try std.testing.expectEqual(-0x6.45fe98ea170261df3affd1873614p+4, lgamma(@as(f128, -0x2.dfffffffffffep+4)));
    try std.testing.expectEqual(-0x5.cc00286fb6c43906b8162ba41228p+4, lgamma(@as(f128, -0x2.dffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x3.ac92bf08c438d477f96b45cd3dfep+4, lgamma(@as(f128, -0x2.dffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x3.fa34a985300399b321bee45fbbe8p+4, lgamma(@as(f128, -0x2.dfffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x7.879d73b6e018ba3e942b365ce218p+4, lgamma(@as(f128, -0x2.e00004p+4)));
    // try std.testing.expectEqual(-0x6.45fe98ea1703579922027d069268p+4, lgamma(@as(f128, -0x2.e000000000002p+4)));
    // try std.testing.expectEqual(-0x5.cc00286fb6c439256f530bf98214p+4, lgamma(@as(f128, -0x2.e000000000000004p+4)));
    try std.testing.expectEqual(-0x3.ac92bf08c438d477f96b45cd3e0ep+4, lgamma(@as(f128, -0x2.e000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x3.fa34a985300399b321bee45fc396p+4, lgamma(@as(f128, -0x2.e0000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x7.c5378948fb9197193badba79d9b4p+4, lgamma(@as(f128, -0x2.effffcp+4)));
    try std.testing.expectEqual(-0x6.8398cd4938e3cde4079bbd524398p+4, lgamma(@as(f128, -0x2.efffffffffffep+4)));
    try std.testing.expectEqual(-0x6.099a5cced8a5a5b9bb29ebdc0f58p+4, lgamma(@as(f128, -0x2.effffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x3.ea2cf367e61a412b12488e30ce4p+4, lgamma(@as(f128, -0x2.effffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x4.37cedde451e506663a9c2cc34c24p+4, lgamma(@as(f128, -0x2.efffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0x7.c537a82bcb8243af457b045a7028p+4, lgamma(@as(f128, -0x2.f00004p+4)));
    try std.testing.expectEqual(-0x6.8398cd4938e4c4fa87212202a56p+4, lgamma(@as(f128, -0x2.f000000000002p+4)));
    // try std.testing.expectEqual(-0x6.099a5cced8a5a5d89df9dc88a564p+4, lgamma(@as(f128, -0x2.f000000000000004p+4)));
    try std.testing.expectEqual(-0x3.ea2cf367e61a412b12488e30ce4ep+4, lgamma(@as(f128, -0x2.f000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x4.37cedde451e506663a9c2cc353dcp+4, lgamma(@as(f128, -0x2.f0000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x8.0327f9ac47b31c8d5f780da3bc68p+4, lgamma(@as(f128, -0x2.fffffcp+4)));
    try std.testing.expectEqual(-0x6.c1893dc1da5ab63bb9ab9862ea3cp+4, lgamma(@as(f128, -0x2.ffffffffffffep+4)));
    try std.testing.expectEqual(-0x6.478acd477a1c8ebc028f1c420b54p+4, lgamma(@as(f128, -0x2.fffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x4.281d63e087912a2d6f0313ec1f9p+4, lgamma(@as(f128, -0x2.fffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x4.75bf4e5cf35bef689756b27e9d7p+4, lgamma(@as(f128, -0x2.ffffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0x8.032818b9c24e73ce14094adffa6p+4, lgamma(@as(f128, -0x3.000004p+4)));
    try std.testing.expectEqual(-0x6.c1893dc1da5baea78e865268a158p+4, lgamma(@as(f128, -0x3.0000000000002p+4)));
    try std.testing.expectEqual(-0x6.478acd477a1c8edb1009b7994c0cp+4, lgamma(@as(f128, -0x3.0000000000000004p+4)));
    try std.testing.expectEqual(-0x4.281d63e087912a2d6f0313ec1fap+4, lgamma(@as(f128, -0x3.0000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x4.75bf4e5cf35bef689756b27ea534p+4, lgamma(@as(f128, -0x3.00000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x8.416cdef3c68716616cd6c95c6898p+4, lgamma(@as(f128, -0x3.0ffffcp+4)));
    try std.testing.expectEqual(-0x6.ffce231e3f0f643d6978f1c4a53p+4, lgamma(@as(f128, -0x3.0fffffffffffep+4)));
    try std.testing.expectEqual(-0x6.85cfb2a3ded13d64cc7bce74c114p+4, lgamma(@as(f128, -0x3.0ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x4.6662493cec45d8d64dd5a6c60454p+4, lgamma(@as(f128, -0x3.0ffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x4.b40433b958109e1176294558823p+4, lgamma(@as(f128, -0x3.0fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x8.416cfe2b0ce3bc002bf2bb5c03b8p+4, lgamma(@as(f128, -0x3.100004p+4)));
    try std.testing.expectEqual(-0x6.ffce231e3f105df79c5e1ebaafe4p+4, lgamma(@as(f128, -0x3.1000000000002p+4)));
    try std.testing.expectEqual(-0x6.85cfb2a3ded13d8403c22b1a5fd4p+4, lgamma(@as(f128, -0x3.1000000000000004p+4)));
    try std.testing.expectEqual(-0x4.6662493cec45d8d64dd5a6c60464p+4, lgamma(@as(f128, -0x3.1000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x4.b40433b958109e117629455889fcp+4, lgamma(@as(f128, -0x3.10000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x8.8004844ea3dd20089f685a7f9a28p+4, lgamma(@as(f128, -0x3.1ffffcp+4)));
    // try std.testing.expectEqual(-0x7.3e65c88d9746c20a4afbe430428p+4, lgamma(@as(f128, -0x3.1fffffffffffep+4)));
    try std.testing.expectEqual(-0x6.c467581337089bd5708e1d095428p+4, lgamma(@as(f128, -0x3.1ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x4.a4f9eeac447d37470662d6a2457cp+4, lgamma(@as(f128, -0x3.1ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x4.f29bd928b047fc822eb67534c354p+4, lgamma(@as(f128, -0x3.1fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x8.8004a3aedffc550387906dae126p+4, lgamma(@as(f128, -0x3.200004p+4)));
    // try std.testing.expectEqual(-0x7.3e65c88d9747bd0c2bf58c0794ep+4, lgamma(@as(f128, -0x3.2000000000002p+4)));
    // try std.testing.expectEqual(-0x6.c467581337089bf4d0ca3c3e4f14p+4, lgamma(@as(f128, -0x3.2000000000000004p+4)));
    // try std.testing.expectEqual(-0x4.a4f9eeac447d37470662d6a2458cp+4, lgamma(@as(f128, -0x3.2000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x4.f29bd928b047fc822eb67534cb2cp+4, lgamma(@as(f128, -0x3.20000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x8.beed463931cafd90ce384ebd5978p+4, lgamma(@as(f128, -0x3.2ffffcp+4)));
    // try std.testing.expectEqual(-0x7.7d4e8a8c3948bf9f12fc14d66a84p+4, lgamma(@as(f128, -0x3.2fffffffffffep+4)));
    try std.testing.expectEqual(-0x7.03501a11d90a9a0ac51ada3c08cp+4, lgamma(@as(f128, -0x3.2ffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x4.e3e2b0aae67f357c6f03a7e90e28p+4, lgamma(@as(f128, -0x3.2ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x5.31849b275249fab79757467b8bf8p+4, lgamma(@as(f128, -0x3.2fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x8.beed65c196125ab3de9d9e6720cp+4, lgamma(@as(f128, -0x3.300004p+4)));
    try std.testing.expectEqual(-0x7.7d4e8a8c3949bbe23536fdeefe28p+4, lgamma(@as(f128, -0x3.3000000000002p+4)));
    try std.testing.expectEqual(-0x7.03501a11d90a9a2a4d7f21992bdp+4, lgamma(@as(f128, -0x3.3000000000000004p+4)));
    try std.testing.expectEqual(-0x4.e3e2b0aae67f357c6f03a7e90e38p+4, lgamma(@as(f128, -0x3.3000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x5.31849b275249fab79757467b93dcp+4, lgamma(@as(f128, -0x3.30000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x8.fe25917adde26ef3cd95670ddd98p+4, lgamma(@as(f128, -0x3.3ffffcp+4)));
    // try std.testing.expectEqual(-0x7.bc86d5e1969b50340f5b8bb0da5cp+4, lgamma(@as(f128, -0x3.3fffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.42886567365d2b3d37a1b38c9ffcp+4, lgamma(@as(f128, -0x3.3ffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x5.231afc0043d1c6aef53bbc4d56ap+4, lgamma(@as(f128, -0x3.3ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x5.70bce67caf9c8bea1d8f5adfd46cp+4, lgamma(@as(f128, -0x3.3fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x8.fe25b12aa49ff3795435fc203e1p+4, lgamma(@as(f128, -0x3.400004p+4)));
    try std.testing.expectEqual(-0x7.bc86d5e1969c4db24547afdd1f3cp+4, lgamma(@as(f128, -0x3.4000000000002p+4)));
    // try std.testing.expectEqual(-0x7.42886567365d2b5ce76871112584p+4, lgamma(@as(f128, -0x3.4000000000000004p+4)));
    // try std.testing.expectEqual(-0x5.231afc0043d1c6aef53bbc4d56bp+4, lgamma(@as(f128, -0x3.4000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x5.70bce67caf9c8bea1d8f5adfdc58p+4, lgamma(@as(f128, -0x3.40000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x9.3dabe237d03ebd86fcf7bf96196p+4, lgamma(@as(f128, -0x3.4ffffcp+4)));
    // try std.testing.expectEqual(-0x7.fc0d26b1db14a5027bbf8934098cp+4, lgamma(@as(f128, -0x3.4fffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.820eb6377ad680a6219b6d7069cp+4, lgamma(@as(f128, -0x3.4ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x5.62a14cd0884b1c17f287932c4bdcp+4, lgamma(@as(f128, -0x3.4ffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x5.b043374cf415e1531adb31bec9a4p+4, lgamma(@as(f128, -0x3.4fffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0x9.3dac020e3b363863752d871afc18p+4, lgamma(@as(f128, -0x3.500004p+4)));
    // try std.testing.expectEqual(-0x7.fc0d26b1db15a3b5d37b6017da8p+4, lgamma(@as(f128, -0x3.5000000000002p+4)));
    // try std.testing.expectEqual(-0x7.820eb6377ad680c5f80664eb4638p+4, lgamma(@as(f128, -0x3.5000000000000004p+4)));
    try std.testing.expectEqual(-0x5.62a14cd0884b1c17f287932c4becp+4, lgamma(@as(f128, -0x3.5000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x5.b043374cf415e1531adb31bed198p+4, lgamma(@as(f128, -0x3.50000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x9.7d7ec3145de00c0a087938f3dc5p+4, lgamma(@as(f128, -0x3.5ffffcp+4)));
    try std.testing.expectEqual(-0x8.3be007a15f3abbcbc2fca1e3ff8p+4, lgamma(@as(f128, -0x3.5fffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.c1e19726fefc98070a07ee6c39ccp+4, lgamma(@as(f128, -0x3.5ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x5.a2742dc00c713378edea98e5bd18p+4, lgamma(@as(f128, -0x3.5ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x5.f016183c783bf8b4163e37783ad8p+4, lgamma(@as(f128, -0x3.5fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x9.7d7ee310b5e10228df915a274458p+4, lgamma(@as(f128, -0x3.600004p+4)));
    // try std.testing.expectEqual(-0x8.3be007a15f3bbbae830452dac6f8p+4, lgamma(@as(f128, -0x3.6000000000002p+4)));
    // try std.testing.expectEqual(-0x7.c1e19726fefc9827065fef6258a4p+4, lgamma(@as(f128, -0x3.6000000000000004p+4)));
    // try std.testing.expectEqual(-0x5.a2742dc00c713378edea98e5bd28p+4, lgamma(@as(f128, -0x3.6000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x5.f016183c783bf8b4163e377842d8p+4, lgamma(@as(f128, -0x3.60000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x9.bd9ccc68ab9aa22b549067817f38p+4, lgamma(@as(f128, -0x3.6ffffcp+4)));
    // try std.testing.expectEqual(-0x8.7bfe11084b36861147a44cae1aep+4, lgamma(@as(f128, -0x3.6fffffffffffep+4)));
    // try std.testing.expectEqual(-0x8.01ffa08deaf862e16e1aa72d0608p+4, lgamma(@as(f128, -0x3.6ffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x5.e2923726f86cfe53649b92d06d68p+4, lgamma(@as(f128, -0x3.6ffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x6.303421a36437c38e8cef3162eb24p+4, lgamma(@as(f128, -0x3.6fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x9.bd9cec8a401dec1250f5d987d4d8p+4, lgamma(@as(f128, -0x3.700004p+4)));
    try std.testing.expectEqual(-0x8.7bfe11084b37871debbe9be60c38p+4, lgamma(@as(f128, -0x3.7000000000002p+4)));
    // try std.testing.expectEqual(-0x8.01ffa08deaf863018faf2a76ed08p+4, lgamma(@as(f128, -0x3.7000000000000004p+4)));
    try std.testing.expectEqual(-0x5.e2923726f86cfe53649b92d06d78p+4, lgamma(@as(f128, -0x3.7000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x6.303421a36437c38e8cef3162f32cp+4, lgamma(@as(f128, -0x3.70000000000000000000000001p+4)));
    try std.testing.expectEqual(-0x9.fe04a3830c274393e4e68be74c1p+4, lgamma(@as(f128, -0x3.7ffffcp+4)));
    try std.testing.expectEqual(-0x8.bc65e834f4e7c3a3a3c3b57e6968p+4, lgamma(@as(f128, -0x3.7fffffffffffep+4)));
    try std.testing.expectEqual(-0x8.426777ba94a9a10601157db43p+4, lgamma(@as(f128, -0x3.7ffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x6.22fa0e53a21e3c7809df8de9e084p+4, lgamma(@as(f128, -0x3.7ffffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x6.709bf8d00de901b332332c7c5e4p+4, lgamma(@as(f128, -0x3.7fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x9.fe04c3c932f3b20d2a807c4fa728p+4, lgamma(@as(f128, -0x3.800004p+4)));
    try std.testing.expectEqual(-0x8.bc65e834f4e8c5d4da272948a3e8p+4, lgamma(@as(f128, -0x3.8000000000002p+4)));
    // try std.testing.expectEqual(-0x8.426777ba94a9a126473c4a22a948p+4, lgamma(@as(f128, -0x3.8000000000000004p+4)));
    try std.testing.expectEqual(-0x6.22fa0e53a21e3c7809df8de9e098p+4, lgamma(@as(f128, -0x3.8000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x6.709bf8d00de901b332332c7c665p+4, lgamma(@as(f128, -0x3.80000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0xa.3eb4f9f7cb8c1f3849128b37f388p+4, lgamma(@as(f128, -0x3.8ffffcp+4)));
    // try std.testing.expectEqual(-0x8.fd163ebbab51268f56d68e71a358p+4, lgamma(@as(f128, -0x3.8fffffffffffep+4)));
    // try std.testing.expectEqual(-0x8.8317ce414b1304815a554032a458p+4, lgamma(@as(f128, -0x3.8ffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x6.63aa64da58879ff3751654e615fcp+4, lgamma(@as(f128, -0x3.8ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x6.b14c4f56c452652e9d69f37893bp+4, lgamma(@as(f128, -0x3.8fffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0xa.3eb51a61e0618933cd9c24a727a8p+4, lgamma(@as(f128, -0x3.900004p+4)));
    // try std.testing.expectEqual(-0x8.fd163ebbab5229dffd81de4dd4d8p+4, lgamma(@as(f128, -0x3.9000000000002p+4)));
    // try std.testing.expectEqual(-0x8.8317ce414b1304a1c46a159c9fep+4, lgamma(@as(f128, -0x3.9000000000000004p+4)));
    try std.testing.expectEqual(-0x6.63aa64da58879ff3751654e6160cp+4, lgamma(@as(f128, -0x3.9000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x6.b14c4f56c452652e9d69f3789bccp+4, lgamma(@as(f128, -0x3.90000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0xa.7fac8cfd3cebe975b34284afbc28p+4, lgamma(@as(f128, -0x3.9ffffcp+4)));
    // try std.testing.expectEqual(-0x9.3e0dd1d2c46a5b17a2eafeb7fb5p+4, lgamma(@as(f128, -0x3.9fffffffffffep+4)));
    // try std.testing.expectEqual(-0x8.c40f6158642c3996d28cffebbe88p+4, lgamma(@as(f128, -0x3.9ffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x6.a4a1f7f171a0d508fef5ce004aa4p+4, lgamma(@as(f128, -0x3.9ffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x6.f243e26ddd6b9a4427496c92c854p+4, lgamma(@as(f128, -0x3.9fffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0xa.7facad8aa13415a62f06963c78p+4, lgamma(@as(f128, -0x3.a00004p+4)));
    // try std.testing.expectEqual(-0x9.3e0dd1d2c46b5f82c52c603be63p+4, lgamma(@as(f128, -0x3.a000000000002p+4)));
    // try std.testing.expectEqual(-0x8.c40f6158642c39b75ff14817efp+4, lgamma(@as(f128, -0x3.a000000000000004p+4)));
    // try std.testing.expectEqual(-0x6.a4a1f7f171a0d508fef5ce004ab4p+4, lgamma(@as(f128, -0x3.a000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x6.f243e26ddd6b9a4427496c92d078p+4, lgamma(@as(f128, -0x3.a0000000000000000000000001p+4)));
    try std.testing.expectEqual(-0xa.c0ea24d2fe73637149acc1a8d5a8p+4, lgamma(@as(f128, -0x3.affffcp+4)));
    // try std.testing.expectEqual(-0x9.7f4b69b9e1103d675b305b49b3a8p+4, lgamma(@as(f128, -0x3.afffffffffffep+4)));
    // try std.testing.expectEqual(-0x9.054cf93f80d21c71526a39c73a28p+4, lgamma(@as(f128, -0x3.affffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x6.e5df8fd88e46b7e3902e263b3b6cp+4, lgamma(@as(f128, -0x3.affffffffffffffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.33817a54fa117d1eb881c4cdb918p+4, lgamma(@as(f128, -0x3.afffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0xa.c0ea458318f84e8c139879496c5p+4, lgamma(@as(f128, -0x3.b00004p+4)));
    // try std.testing.expectEqual(-0x9.7f4b69b9e11142e82f57b4200f6p+4, lgamma(@as(f128, -0x3.b000000000002p+4)));
    // try std.testing.expectEqual(-0x9.054cf93f80d21c920284beb254fp+4, lgamma(@as(f128, -0x3.b000000000000004p+4)));
    // try std.testing.expectEqual(-0x6.e5df8fd88e46b7e3902e263b3b7cp+4, lgamma(@as(f128, -0x3.b000000000000000000000000002p+4)));
    try std.testing.expectEqual(-0x7.33817a54fa117d1eb881c4cdc144p+4, lgamma(@as(f128, -0x3.b0000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0xb.026c9433822c767deece73d0b428p+4, lgamma(@as(f128, -0x3.bffffcp+4)));
    try std.testing.expectEqual(-0x9.c0cdd92b75da6a16b41d5c5adc68p+4, lgamma(@as(f128, -0x3.bfffffffffffep+4)));
    // try std.testing.expectEqual(-0x9.46cf68b1159c49a922ceb24fda6p+4, lgamma(@as(f128, -0x3.bffffffffffffffcp+4)));
    // try std.testing.expectEqual(-0x7.2761ff4a2310e51b71a3afd4ecb4p+4, lgamma(@as(f128, -0x3.bffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(-0x7.7503e9c68edbaa5699f74e676a5cp+4, lgamma(@as(f128, -0x3.bfffffffffffffffffffffffffp+4)));
    // try std.testing.expectEqual(-0xb.026cb505bed383badae93f9807p+4, lgamma(@as(f128, -0x3.c00004p+4)));
    try std.testing.expectEqual(-0x9.c0cdd92b75db70a89955c642493p+4, lgamma(@as(f128, -0x3.c000000000002p+4)));
    try std.testing.expectEqual(-0x9.46cf68b1159c49c9f50b595d1748p+4, lgamma(@as(f128, -0x3.c000000000000004p+4)));
    // try std.testing.expectEqual(-0x7.2761ff4a2310e51b71a3afd4ecc4p+4, lgamma(@as(f128, -0x3.c000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x7.7503e9c68edbaa5699f74e67729p+4, lgamma(@as(f128, -0x3.c0000000000000000000000001p+4)));
    // try std.testing.expectEqual(0x4.2b2b52b5464eed2ad208f29fe2b8p-24, lgamma(@as(f128, -0x2.74ff9p+0)));
    // try std.testing.expectEqual(-0x1.e4cf2421a71b194f877fec38af2ep-24, lgamma(@as(f128, -0x2.74ff94p+0)));
    try std.testing.expectEqual(0x4.0c8edb47fa1b350b6db03f5fc5c8p-56, lgamma(@as(f128, -0x2.74ff92c01f0d8p+0)));
    try std.testing.expectEqual(-0x2.c7343f216ac923ce0f23c10783bcp-52, lgamma(@as(f128, -0x2.74ff92c01f0dap+0)));
    // try std.testing.expectEqual(0x5.f29bbbdec3d4a8d5f94b191f425p-64, lgamma(@as(f128, -0x2.74ff92c01f0d82a8p+0)));
    try std.testing.expectEqual(-0x1.d5e9dcd11030bba08ecffc3ee9e9p-68, lgamma(@as(f128, -0x2.74ff92c01f0d82acp+0)));
    // try std.testing.expectEqual(0x1.a06e577d39632215de23f7cd18afp-112, lgamma(@as(f128, -0x2.74ff92c01f0d82abec9f315f1a06p+0)));
    try std.testing.expectEqual(-0x1.678ed558b108b80cbdb57a893e71p-112, lgamma(@as(f128, -0x2.74ff92c01f0d82abec9f315f1a08p+0)));
    try std.testing.expectEqual(0xa.b865ddfef8a6b07db1b04ed01ef8p-112, lgamma(@as(f128, -0x2.74ff92c01f0d82abec9f315f1ap+0)));
    // try std.testing.expectEqual(-0x1.7946308cf63d4660d03b08dc56dap-104, lgamma(@as(f128, -0x2.74ff92c01f0d82abec9f315f1bp+0)));
    try std.testing.expectEqual(-0x2.6b416efc56fe3eb64a9ce8defc3ep-24, lgamma(@as(f128, -0x2.bf682p+0)));
    try std.testing.expectEqual(0x5.3d0a33adaf4f5876308316ab3e4p-24, lgamma(@as(f128, -0x2.bf6824p+0)));
    // try std.testing.expectEqual(-0x3.0c498b9ac27bd8d3bd09edf35378p-52, lgamma(@as(f128, -0x2.bf6821437b2p+0)));
    try std.testing.expectEqual(0xc.7dc2985d3b445569287e423f8a9p-56, lgamma(@as(f128, -0x2.bf6821437b202p+0)));
    // try std.testing.expectEqual(-0x3.088b212f3705dc7ce8ce2ecc3a8ep-64, lgamma(@as(f128, -0x2.bf6821437b201978p+0)));
    try std.testing.expectEqual(0x4.9fc04911f55d35c2fd558331628cp-64, lgamma(@as(f128, -0x2.bf6821437b20197cp+0)));
    try std.testing.expectEqual(-0x1.6e33a5f13f9e3c6bdc30fa8b829ap-112, lgamma(@as(f128, -0x2.bf6821437b20197995a4b4641eaep+0)));
    try std.testing.expectEqual(0x2.65f20f2f56934ca42d56d24de4c8p-112, lgamma(@as(f128, -0x2.bf6821437b20197995a4b4641ebp+0)));
    try std.testing.expectEqual(-0x1.4e870434044a73d0df195798692bp-104, lgamma(@as(f128, -0x2.bf6821437b20197995a4b4641ep+0)));
    // try std.testing.expectEqual(0x9.b8bd65c46ce50b725aa8ed447018p-108, lgamma(@as(f128, -0x2.bf6821437b20197995a4b4641fp+0)));
    // try std.testing.expectEqual(0x1.bd69b50d51b1488028c7a2b72e33p-20, lgamma(@as(f128, -0x3.24c1b4p+0)));
    // try std.testing.expectEqual(-0x3.4a0c544eeb21a026dc79de4e099ap-24, lgamma(@as(f128, -0x3.24c1b8p+0)));
    // try std.testing.expectEqual(0x7.a58178eb9e9877664321f5a1f468p-52, lgamma(@as(f128, -0x3.24c1b793cb35ep+0)));
    // try std.testing.expectEqual(-0x7.ead1b6ac3791da04d17fd2da7cacp-52, lgamma(@as(f128, -0x3.24c1b793cb36p+0)));
    try std.testing.expectEqual(0x5.c9c4ac92390bb71a2f0034baceacp-64, lgamma(@as(f128, -0x3.24c1b793cb35efb8p+0)));
    try std.testing.expectEqual(-0x1.956e1b29d7349243d86eea41d5dbp-60, lgamma(@as(f128, -0x3.24c1b793cb35efbcp+0)));
    // try std.testing.expectEqual(0x3.1413b11d6bffa548f8c9e81a7f92p-112, lgamma(@as(f128, -0x3.24c1b793cb35efb8be699ad3d9bap+0)));
    try std.testing.expectEqual(-0xc.7c3f7e7a6a2ac8e4ae516726fc2p-112, lgamma(@as(f128, -0x3.24c1b793cb35efb8be699ad3d9bcp+0)));
    try std.testing.expectEqual(0x5.aa824bfb463969abdeadb5b2f238p-104, lgamma(@as(f128, -0x3.24c1b793cb35efb8be699ad3d9p+0)));
    // try std.testing.expectEqual(-0x2.1da74bd0a4dbcd6af4dff1edd78ap-104, lgamma(@as(f128, -0x3.24c1b793cb35efb8be699ad3dap+0)));
    // try std.testing.expectEqual(-0x3.511bca412890968ef5acdaae7dbep-20, lgamma(@as(f128, -0x3.f48e28p+0)));
    // try std.testing.expectEqual(0x1.dd4b54ca863c1a476cbd9fd337c3p-20, lgamma(@as(f128, -0x3.f48e2cp+0)));
    // try std.testing.expectEqual(-0x1.ddc0336980b584d18e3a66026b1p-52, lgamma(@as(f128, -0x3.f48e2a8f85fcap+0)));
    try std.testing.expectEqual(0x2.7957af96f2c10efef97cd0a7f8b6p-48, lgamma(@as(f128, -0x3.f48e2a8f85fccp+0)));
    try std.testing.expectEqual(-0x1.130ae5c4f54dbe9194a6099cc513p-60, lgamma(@as(f128, -0x3.f48e2a8f85fca17p+0)));
    // try std.testing.expectEqual(0x4.1b5c7fd62043e8361bf7df13d68p-60, lgamma(@as(f128, -0x3.f48e2a8f85fca174p+0)));
    try std.testing.expectEqual(-0xf.cc078f3d044cc25d934ed1413dfp-112, lgamma(@as(f128, -0x3.f48e2a8f85fca170d4561291236cp+0)));
    try std.testing.expectEqual(0x1.9a7339d9ba8406f455a88d6bf278p-108, lgamma(@as(f128, -0x3.f48e2a8f85fca170d4561291236ep+0)));
    // try std.testing.expectEqual(-0x8.ce1a8304f16a153abbbecc4129p-104, lgamma(@as(f128, -0x3.f48e2a8f85fca170d456129123p+0)));
    // try std.testing.expectEqual(0xb.eb83136764dc8396bb2d07bf2f1p-104, lgamma(@as(f128, -0x3.f48e2a8f85fca170d456129124p+0)));
    // try std.testing.expectEqual(0xa.3165c90424948cf2db600b526e28p-20, lgamma(@as(f128, -0x4.0a1398p+0)));
    // try std.testing.expectEqual(-0x3.33cb5626dc331ecdd9f62339e71ap-20, lgamma(@as(f128, -0x4.0a13ap+0)));
    try std.testing.expectEqual(0x5.1a6a3781914476625cd69c566e5p-48, lgamma(@as(f128, -0x4.0a139e16656p+0)));
    try std.testing.expectEqual(-0x1.982d05a2f456b4af798342c40356p-48, lgamma(@as(f128, -0x4.0a139e1665604p+0)));
    // try std.testing.expectEqual(0x6.103eebf7b96ec358066dd1892148p-60, lgamma(@as(f128, -0x4.0a139e16656030cp+0)));
    // try std.testing.expectEqual(-0x7.54ef8e5151b25673cdaf3d854ca4p-60, lgamma(@as(f128, -0x4.0a139e16656030c8p+0)));
    try std.testing.expectEqual(0x4.7972c9f44c23cd27b57158763d9p-108, lgamma(@as(f128, -0x4.0a139e16656030c39f0b0de1811p+0)));
    // try std.testing.expectEqual(-0x2.39247330396cbffb9bff1408caf6p-108, lgamma(@as(f128, -0x4.0a139e16656030c39f0b0de18114p+0)));
    try std.testing.expectEqual(0x1.cbe99f07a7c6894a89574e2a4bb8p-100, lgamma(@as(f128, -0x4.0a139e16656030c39f0b0de18p+0)));
    // try std.testing.expectEqual(-0x1.8d61ff8a9b01bd471f60e8153e3ep-100, lgamma(@as(f128, -0x4.0a139e16656030c39f0b0de182p+0)));
    try std.testing.expectEqual(-0x3.02165b2aa6eef1f3030056865942p-16, lgamma(@as(f128, -0x4.fdd5d8p+0)));
    try std.testing.expectEqual(0xa.22e7861540c9fce321e0c06dbc9p-20, lgamma(@as(f128, -0x4.fdd5ep+0)));
    // try std.testing.expectEqual(-0x1.8280d0ba86860c975c015010a996p-44, lgamma(@as(f128, -0x4.fdd5de9bbabfp+0)));
    // try std.testing.expectEqual(0x4.fa3d33517a9ecdcbd38ddb02bfa8p-48, lgamma(@as(f128, -0x4.fdd5de9bbabf4p+0)));
    try std.testing.expectEqual(-0x5.efcf1ba2a53065f6ccd269593bcp-60, lgamma(@as(f128, -0x4.fdd5de9bbabf351p+0)));
    // try std.testing.expectEqual(0x3.454c56251230ebf0190e0a892456p-56, lgamma(@as(f128, -0x4.fdd5de9bbabf3518p+0)));
    // try std.testing.expectEqual(-0x7.55ff3704a7af9f64dfe331ebcf94p-108, lgamma(@as(f128, -0x4.fdd5de9bbabf3510d0aa40769884p+0)));
    try std.testing.expectEqual(0x1.5cc4b07f53c6fc793bcc51cff2ebp-104, lgamma(@as(f128, -0x4.fdd5de9bbabf3510d0aa40769888p+0)));
    // try std.testing.expectEqual(-0x3.c8c191553b0fbbe57111955dbed4p-100, lgamma(@as(f128, -0x4.fdd5de9bbabf3510d0aa407698p+0)));
    try std.testing.expectEqual(0xa.c8638e27b6fff796dd42921b01b8p-100, lgamma(@as(f128, -0x4.fdd5de9bbabf3510d0aa40769ap+0)));
    try std.testing.expectEqual(0x2.e258f12a679ed407ae94080315ccp-16, lgamma(@as(f128, -0x5.021a9p+0)));
    try std.testing.expectEqual(-0xf.89066929e3b180fd518e8b9d62f8p-20, lgamma(@as(f128, -0x5.021a98p+0)));
    // try std.testing.expectEqual(0x1.867827fdc0e929bcc699f977a672p-48, lgamma(@as(f128, -0x5.021a95fc2db64p+0)));
    try std.testing.expectEqual(-0x1.d50b5e02beb77b1150b98c01af96p-44, lgamma(@as(f128, -0x5.021a95fc2db68p+0)));
    // try std.testing.expectEqual(0x1.1b82d6b2b33045c61b8f03708c2ep-56, lgamma(@as(f128, -0x5.021a95fc2db64328p+0)));
    // try std.testing.expectEqual(-0x2.bf62ea52828ff32acab6b018a736p-56, lgamma(@as(f128, -0x5.021a95fc2db6433p+0)));
    // try std.testing.expectEqual(0xe.d75efeb9083d919f12877a9eb3fp-108, lgamma(@as(f128, -0x5.021a95fc2db6432a4c56e595394cp+0)));
    try std.testing.expectEqual(-0xf.ffcf0970a5c44e84d51ed6bcd858p-108, lgamma(@as(f128, -0x5.021a95fc2db6432a4c56e595395p+0)));
    // try std.testing.expectEqual(0xa.0e9b4ba43c72d93d432d73de60fp-100, lgamma(@as(f128, -0x5.021a95fc2db6432a4c56e59538p+0)));
    try std.testing.expectEqual(-0x5.5cfbb8709a8e16d4b0a5b4d1993cp-100, lgamma(@as(f128, -0x5.021a95fc2db6432a4c56e5953ap+0)));
    // try std.testing.expectEqual(-0xf.15ee1077e22d21b977289dc12a58p-16, lgamma(@as(f128, -0x5.ffa4b8p+0)));
    // try std.testing.expectEqual(0x7.4bb0ef1ad813da34a81bb0995568p-16, lgamma(@as(f128, -0x5.ffa4cp+0)));
    try std.testing.expectEqual(-0x4.2c4d3e7ff051f430d17064abe11p-44, lgamma(@as(f128, -0x5.ffa4bd647d034p+0)));
    // try std.testing.expectEqual(0x7.04ae139d3fb740351122ea1d2804p-44, lgamma(@as(f128, -0x5.ffa4bd647d038p+0)));
    try std.testing.expectEqual(-0xe.d9cc85177f957fa2719e31081b5p-56, lgamma(@as(f128, -0x5.ffa4bd647d0357d8p+0)));
    // try std.testing.expectEqual(0x7.882a1f22de7c7117c696484fb6ccp-56, lgamma(@as(f128, -0x5.ffa4bd647d0357ep+0)));
    try std.testing.expectEqual(-0x5.84998680d25d3b2fae7819c285d4p-104, lgamma(@as(f128, -0x5.ffa4bd647d0357dd4ed62cbd31ecp+0)));
    try std.testing.expectEqual(0x5.ac61cb9c5cabe658c0f85dbf47e4p-104, lgamma(@as(f128, -0x5.ffa4bd647d0357dd4ed62cbd31fp+0)));
    try std.testing.expectEqual(-0x5.660d59fa866bc057bd39817679ccp-96, lgamma(@as(f128, -0x5.ffa4bd647d0357dd4ed62cbd3p+0)));
    // try std.testing.expectEqual(0x3.2704f141118d06c7a7eba3bcb69cp-100, lgamma(@as(f128, -0x5.ffa4bd647d0357dd4ed62cbd32p+0)));
    // try std.testing.expectEqual(0x3.e9df593e904f847b411ee284216ap-16, lgamma(@as(f128, -0x6.005ac8p+0)));
    try std.testing.expectEqual(-0x1.2b35eea26dc93cd1810a63cdaf1dp-12, lgamma(@as(f128, -0x6.005adp+0)));
    try std.testing.expectEqual(0xa.7dd3bd697d2c7b894581895ac238p-44, lgamma(@as(f128, -0x6.005ac9625f23p+0)));
    // try std.testing.expectEqual(-0xd.11e91b3ff8f4b947413d3dd522d8p-48, lgamma(@as(f128, -0x6.005ac9625f234p+0)));
    // try std.testing.expectEqual(0x1.5f103a1b00a487010706431c9b3ap-56, lgamma(@as(f128, -0x6.005ac9625f233b6p+0)));
    // try std.testing.expectEqual(-0x1.53ed4641ff204b2bd9b67f3df5e9p-52, lgamma(@as(f128, -0x6.005ac9625f233b68p+0)));
    // try std.testing.expectEqual(0x5.131fd8e3456e59f544101180674cp-104, lgamma(@as(f128, -0x6.005ac9625f233b607c2d96d16384p+0)));
    // try std.testing.expectEqual(-0x6.3bd2763a33e6b2b51738389d2bccp-104, lgamma(@as(f128, -0x6.005ac9625f233b607c2d96d16388p+0)));
    // try std.testing.expectEqual(0x4.4dfcefd30e3ea82681da742fef2cp-96, lgamma(@as(f128, -0x6.005ac9625f233b607c2d96d162p+0)));
    try std.testing.expectEqual(-0x1.597c37bbae6bde2eabc9b0e72d53p-96, lgamma(@as(f128, -0x6.005ac9625f233b607c2d96d164p+0)));
    // try std.testing.expectEqual(-0x7.313b92969004728fd3b1bc642318p-12, lgamma(@as(f128, -0x6.fff2f8p+0)));
    try std.testing.expectEqual(0x2.a3598cd9f522a41184bfaa007bf2p-12, lgamma(@as(f128, -0x6.fff3p+0)));
    try std.testing.expectEqual(-0x4.dc097be5d1cc3e05e676e18044ap-40, lgamma(@as(f128, -0x6.fff2fddae1bbcp+0)));
    try std.testing.expectEqual(0xe.f46d8dcca9e8196c5247230ffccp-48, lgamma(@as(f128, -0x6.fff2fddae1bcp+0)));
    // try std.testing.expectEqual(-0x6.9ebebbccaa51db7332b540aeeea4p-52, lgamma(@as(f128, -0x6.fff2fddae1bbff38p+0)));
    // try std.testing.expectEqual(0x3.373d171aaa3ac8c09e470f086dbap-52, lgamma(@as(f128, -0x6.fff2fddae1bbff4p+0)));
    // try std.testing.expectEqual(-0x2.9beb0f7a288ab0d66ab8a8c438c4p-100, lgamma(@as(f128, -0x6.fff2fddae1bbff3d626b65c23fdp+0)));
    try std.testing.expectEqual(0x2.4f12d9f981bc274f10d21440cca4p-100, lgamma(@as(f128, -0x6.fff2fddae1bbff3d626b65c23fd4p+0)));
    try std.testing.expectEqual(-0x2.3d16f8d7e350a4a1d2659626f086p-92, lgamma(@as(f128, -0x6.fff2fddae1bbff3d626b65c23ep+0)));
    // try std.testing.expectEqual(0x3.867fbe1f1d2c770eb5fc833dbe12p-96, lgamma(@as(f128, -0x6.fff2fddae1bbff3d626b65c24p+0)));
    try std.testing.expectEqual(0x9.39801333caa3621e695e453c3f98p-12, lgamma(@as(f128, -0x7.000cf8p+0)));
    try std.testing.expectEqual(-0xa.32834623023dc44da1eb1b2327p-16, lgamma(@as(f128, -0x7.000dp+0)));
    // try std.testing.expectEqual(0x3.89727e62d4843a650676bf6a43a4p-40, lgamma(@as(f128, -0x7.000cff7b7f878p+0)));
    try std.testing.expectEqual(-0x1.638f6c2b4fb951525a3d1dd68618p-40, lgamma(@as(f128, -0x7.000cff7b7f87cp+0)));
    // try std.testing.expectEqual(0x5.45e474b8c68eb4e8057e7bcdb0f4p-52, lgamma(@as(f128, -0x7.000cff7b7f87adfp+0)));
    try std.testing.expectEqual(-0x4.941f6063775a1f6940ea2871a034p-52, lgamma(@as(f128, -0x7.000cff7b7f87adf8p+0)));
    // try std.testing.expectEqual(0x3.49444b311bfe7229d00b52c56e9ep-100, lgamma(@as(f128, -0x7.000cff7b7f87adf4482dcdb9878p+0)));
    // try std.testing.expectEqual(-0x1.a3bd9f5d02f5dca722ae4f13c9a3p-100, lgamma(@as(f128, -0x7.000cff7b7f87adf4482dcdb98784p+0)));
    try std.testing.expectEqual(0x1.dc29fc407cb79c0084d5a81fc999p-92, lgamma(@as(f128, -0x7.000cff7b7f87adf4482dcdb986p+0)));
    try std.testing.expectEqual(-0x9.a56f90692c28b67f48728e572078p-96, lgamma(@as(f128, -0x7.000cff7b7f87adf4482dcdb988p+0)));
    try std.testing.expectEqual(-0x4.cccb8849515a9e45ca27a76b35dcp-8, lgamma(@as(f128, -0x7.fffe58p+0)));
    try std.testing.expectEqual(0x1.37b05f6d428d9a997989792587b5p-12, lgamma(@as(f128, -0x7.fffe6p+0)));
    try std.testing.expectEqual(-0x2.551849c02b7e0c5decab28ca9bap-40, lgamma(@as(f128, -0x7.fffe5fe05673cp+0)));
    try std.testing.expectEqual(0x2.509d5b2dadf1ea40b619bc764a54p-36, lgamma(@as(f128, -0x7.fffe5fe05674p+0)));
    // try std.testing.expectEqual(-0x1.9c7a33ad9478c56c11c3c4797b4ep-48, lgamma(@as(f128, -0x7.fffe5fe05673c3c8p+0)));
    try std.testing.expectEqual(0x3.4f638be5777634c4771c01985932p-48, lgamma(@as(f128, -0x7.fffe5fe05673c3dp+0)));
    try std.testing.expectEqual(-0x1.9ba8efbb83ce91492a4f17354a03p-96, lgamma(@as(f128, -0x7.fffe5fe05673c3ca9e82b522b0c8p+0)));
    try std.testing.expectEqual(0xd.a45f00e0226d4a9fa6d172f7146p-100, lgamma(@as(f128, -0x7.fffe5fe05673c3ca9e82b522b0ccp+0)));
    // try std.testing.expectEqual(-0x7.ca450a517adbc7ac6571008b904p-92, lgamma(@as(f128, -0x7.fffe5fe05673c3ca9e82b522bp+0)));
    // try std.testing.expectEqual(0xb.e531f3fab4cf67ecc07075284a7p-92, lgamma(@as(f128, -0x7.fffe5fe05673c3ca9e82b522b2p+0)));
    // try std.testing.expectEqual(0xc.8602745a44910cd1b3729bc87778p-16, lgamma(@as(f128, -0x8.0001ap+0)));
    // try std.testing.expectEqual(-0x9.9cf5dfb6141ef53063d367d23668p-8, lgamma(@as(f128, -0x8.0001bp+0)));
    // try std.testing.expectEqual(0x1.34e935f3e5a5cd0f2efbeff26f35p-36, lgamma(@as(f128, -0x8.0001a01459fc8p+0)));
    // try std.testing.expectEqual(-0x3.b73909c1555516b085e6b72d9978p-36, lgamma(@as(f128, -0x8.0001a01459fdp+0)));
    // try std.testing.expectEqual(0x7.d0d61593cad19969bee6042f80f8p-52, lgamma(@as(f128, -0x8.0001a01459fc9f6p+0)));
    // try std.testing.expectEqual(-0x9.5b371e11feb316e9319036a2cc18p-48, lgamma(@as(f128, -0x8.0001a01459fc9f7p+0)));
    try std.testing.expectEqual(0x3.5c700de7e9cd30f840dc30cb631ap-96, lgamma(@as(f128, -0x8.0001a01459fc9f60cb3cec1cec8p+0)));
    try std.testing.expectEqual(-0x1.8fb231cdb3f8ba025de82c9ed996p-96, lgamma(@as(f128, -0x8.0001a01459fc9f60cb3cec1cec88p+0)));
    try std.testing.expectEqual(0x5.21e940941c62be0a22d22144e1cp-92, lgamma(@as(f128, -0x8.0001a01459fc9f60cb3cec1cecp+0)));
    // try std.testing.expectEqual(-0x2.23f28bd18d1cc99cad350a61697cp-88, lgamma(@as(f128, -0x8.0001a01459fc9f60cb3cec1cfp+0)));
    try std.testing.expectEqual(-0x9.98ed0cd062e3fd423095cc578b1p-8, lgamma(@as(f128, -0x8.ffffdp+0)));
    // try std.testing.expectEqual(0x5.e337e9ef84f0aaa1574e4105202cp-4, lgamma(@as(f128, -0x8.ffffep+0)));
    try std.testing.expectEqual(-0x5.88479ad476d496a4dd586c5c1bf8p-36, lgamma(@as(f128, -0x8.ffffd1c425e8p+0)));
    try std.testing.expectEqual(0x2.6c3945e213fff55a9528d2a68d04p-32, lgamma(@as(f128, -0x8.ffffd1c425e88p+0)));
    // try std.testing.expectEqual(-0x4.55973b8ddaa3ac8b9449bc0f0998p-44, lgamma(@as(f128, -0x8.ffffd1c425e80ffp+0)));
    // try std.testing.expectEqual(0x1.33e4438b1b9d4bdcab5c1cb2606ap-44, lgamma(@as(f128, -0x8.ffffd1c425e81p+0)));
    // try std.testing.expectEqual(-0xa.8d6d593f44a9c110c1a62526e19p-96, lgamma(@as(f128, -0x8.ffffd1c425e80ffc864e95749258p+0)));
    // try std.testing.expectEqual(0x2.1be6e9f8871b3c84f4d5e8289284p-92, lgamma(@as(f128, -0x8.ffffd1c425e80ffc864e9574926p+0)));
    // try std.testing.expectEqual(-0xd.04c6df3bc1b211003527adf10bb8p-88, lgamma(@as(f128, -0x8.ffffd1c425e80ffc864e95749p+0)));
    // try std.testing.expectEqual(0x9.21271d28197cb3afd25a79d46348p-88, lgamma(@as(f128, -0x8.ffffd1c425e80ffc864e957494p+0)));
    // try std.testing.expectEqual(0x5.e32ee82416adc45e301db9b6542cp-4, lgamma(@as(f128, -0x9.00002p+0)));
    // try std.testing.expectEqual(-0x9.99c537e2b92992ab902d09aa68ep-8, lgamma(@as(f128, -0x9.00003p+0)));
    try std.testing.expectEqual(0x2.5debd4969bb286f02ef187d99c92p-36, lgamma(@as(f128, -0x9.00002e3bb47d8p+0)));
    try std.testing.expectEqual(-0x2.9ee383255c86df27cc3863245bacp-32, lgamma(@as(f128, -0x9.00002e3bb47ep+0)));
    try std.testing.expectEqual(0x2.5e69b52fd9e19d1256f345b0cfbep-44, lgamma(@as(f128, -0x9.00002e3bb47d86dp+0)));
    // try std.testing.expectEqual(-0x3.2b1acbb48b0afdaeb19034b5fd9ep-44, lgamma(@as(f128, -0x9.00002e3bb47d86ep+0)));
    try std.testing.expectEqual(0x2.0c77b3317df50e2c02edb5cb24bap-92, lgamma(@as(f128, -0x9.00002e3bb47d86d6d843fedc3518p+0)));
    try std.testing.expectEqual(-0xb.84a8d40b492f4ceba172fed99c88p-96, lgamma(@as(f128, -0x9.00002e3bb47d86d6d843fedc352p+0)));
    // try std.testing.expectEqual(0x6.2f30682ce668d7673da9a4bdc4a8p-88, lgamma(@as(f128, -0x9.00002e3bb47d86d6d843fedc34p+0)));
    // try std.testing.expectEqual(-0xf.f6e19b64add7406eaa7d1bce4d1p-88, lgamma(@as(f128, -0x9.00002e3bb47d86d6d843fedc38p+0)));
    try std.testing.expectEqual(-0x1.3dd0c34d79694344018ee202113p+0, lgamma(@as(f128, -0x9.fffffp+0)));
    try std.testing.expectEqual(-0x1.41334d2c3ccaa629ce4c669f13ccp-28, lgamma(@as(f128, -0x9.fffffb606bdf8p+0)));
    try std.testing.expectEqual(0x7.9c48d283217d7932db943025f438p-32, lgamma(@as(f128, -0x9.fffffb606bep+0)));
    try std.testing.expectEqual(-0x1.55818a2b42ba217320b6a0d61e98p-44, lgamma(@as(f128, -0x9.fffffb606bdfdcdp+0)));
    // try std.testing.expectEqual(0x3.60979c1bc0b3232c9100d27de7ccp-40, lgamma(@as(f128, -0x9.fffffb606bdfdcep+0)));
    // try std.testing.expectEqual(-0x1.ae8e54cac61a239123ee45984ab3p-88, lgamma(@as(f128, -0x9.fffffb606bdfdcd062ae77a5054p+0)));
    try std.testing.expectEqual(0xc.698594717bb0c97a661c1e399808p-96, lgamma(@as(f128, -0x9.fffffb606bdfdcd062ae77a50548p+0)));
    // try std.testing.expectEqual(-0x4.6e54873ab758351b9e113f2f461cp-84, lgamma(@as(f128, -0x9.fffffb606bdfdcd062ae77a504p+0)));
    // try std.testing.expectEqual(0x9.696a4bbf05566db954940168b358p-84, lgamma(@as(f128, -0x9.fffffb606bdfdcd062ae77a508p+0)));
    try std.testing.expectEqual(-0x1.3dd10e8f080e8da97df93de56ed2p+0, lgamma(@as(f128, -0xa.00001p+0)));
    // try std.testing.expectEqual(0x5.70ddf269e6d667a1b2a416297d08p-32, lgamma(@as(f128, -0xa.0000049f93bb8p+0)));
    // try std.testing.expectEqual(-0x1.63ea466b9e05b9e2182cf19da66ap-28, lgamma(@as(f128, -0xa.0000049f93bcp+0)));
    try std.testing.expectEqual(0x1.aa9c2e2b1029c57ecdd5e15136d8p-40, lgamma(@as(f128, -0xa.0000049f93bb992p+0)));
    // try std.testing.expectEqual(-0x1.cb541d167c13dafc1bbdbfcaab97p-40, lgamma(@as(f128, -0xa.0000049f93bb993p+0)));
    try std.testing.expectEqual(0x6.7dca7cdca65e5e48553a99e1c884p-92, lgamma(@as(f128, -0xa.0000049f93bb9927b45d95e1544p+0)));
    // try std.testing.expectEqual(-0x1.531b7dd2fbd538fb486d457ffd2p-88, lgamma(@as(f128, -0xa.0000049f93bb9927b45d95e15448p+0)));
    // try std.testing.expectEqual(0xe.3f9dd4d3fc3edce2f35b8e589238p-88, lgamma(@as(f128, -0xa.0000049f93bb9927b45d95e154p+0)));
    // try std.testing.expectEqual(-0xc.f3c74fb8f21509303ecc8ea97cp-84, lgamma(@as(f128, -0xa.0000049f93bb9927b45d95e158p+0)));
    try std.testing.expectEqual(-0x3.a3ad38c9033a659ac104c00477e4p+0, lgamma(@as(f128, -0xa.fffffp+0)));
    // try std.testing.expectEqual(-0x1.0e8528e5ba92d3d57ef9a360b17ep-24, lgamma(@as(f128, -0xa.ffffff9466e98p+0)));
    // try std.testing.expectEqual(0x2.205541c47450d1d3de54560bdabep-28, lgamma(@as(f128, -0xa.ffffff9466eap+0)));
    // try std.testing.expectEqual(-0x8.28300f9cbbc503a9a39ca7d04828p-40, lgamma(@as(f128, -0xa.ffffff9466e9f1bp+0)));
    // try std.testing.expectEqual(0x1.de91fa23a9940dd2694ba76e4fe6p-36, lgamma(@as(f128, -0xa.ffffff9466e9f1cp+0)));
    try std.testing.expectEqual(-0x1.efb7e7a1e33d47532bf23ebbde56p-88, lgamma(@as(f128, -0xa.ffffff9466e9f1b36dacd2adbd18p+0)));
    // try std.testing.expectEqual(0x1.118eff148f83e9e172e04ca13dd5p-84, lgamma(@as(f128, -0xa.ffffff9466e9f1b36dacd2adbd2p+0)));
    try std.testing.expectEqual(-0x2.9c1eaa8fbde52da4dd753ce2f2a6p-80, lgamma(@as(f128, -0xa.ffffff9466e9f1b36dacd2adbcp+0)));
    // try std.testing.expectEqual(0x6.e83541e5afd8c5104f9a72184984p-80, lgamma(@as(f128, -0xa.ffffff9466e9f1b36dacd2adcp+0)));
    try std.testing.expectEqual(-0x3.a3ad86f34c0e3ba328367f78cabp+0, lgamma(@as(f128, -0xb.00001p+0)));
    // try std.testing.expectEqual(0x7.573b0669604304a200ed7fab9af4p-28, lgamma(@as(f128, -0xb.0000006b9915p+0)));
    try std.testing.expectEqual(-0xb.b16d1e1508e7a9bd7460e95f6afp-28, lgamma(@as(f128, -0xb.0000006b99158p+0)));
    try std.testing.expectEqual(0x2.053cabc3adfebe315489e4a8bd3ep-36, lgamma(@as(f128, -0xb.0000006b9915315p+0)));
    // try std.testing.expectEqual(-0x5.bd8591f162c0dabb0c2a6c755abp-40, lgamma(@as(f128, -0xb.0000006b9915316p+0)));
    try std.testing.expectEqual(0xe.fb1c8a928784bc973ad267a8c848p-88, lgamma(@as(f128, -0xb.0000006b9915315d965a6ffea408p+0)));
    try std.testing.expectEqual(-0x4.0d8b9c829ccafedbd8515060eb2p-88, lgamma(@as(f128, -0xb.0000006b9915315d965a6ffea41p+0)));
    try std.testing.expectEqual(0x2.203c4b1a7abd4780a4df789fbb98p-84, lgamma(@as(f128, -0xb.0000006b9915315d965a6ffea4p+0)));
    // try std.testing.expectEqual(-0x9.62504ed8ea7c09417f17ab0a119p-80, lgamma(@as(f128, -0xb.0000006b9915315d965a6ffea8p+0)));
    try std.testing.expectEqual(-0x6.1fd00f0e21b3c98569e28b729b24p+0, lgamma(@as(f128, -0xb.fffffp+0)));
    try std.testing.expectEqual(-0xc.e27c4f01cf5328473b68e6bd241p-28, lgamma(@as(f128, -0xb.fffffff708938p+0)));
    try std.testing.expectEqual(0xd.785692eee5fd5bfcead41cc15278p-24, lgamma(@as(f128, -0xb.fffffff70894p+0)));
    // try std.testing.expectEqual(-0xf.272276e2f7d5551ccfbcd35dbe1p-36, lgamma(@as(f128, -0xb.fffffff70893873p+0)));
    // try std.testing.expectEqual(0xd.65d9840e2817354de03bdbee8ae8p-36, lgamma(@as(f128, -0xb.fffffff70893874p+0)));
    // try std.testing.expectEqual(-0x8.31ab137078f08157c5a627294428p-84, lgamma(@as(f128, -0xb.fffffff7089387387de41acc3d38p+0)));
    // try std.testing.expectEqual(0x6.14d2ea08df7366befa4020c93208p-84, lgamma(@as(f128, -0xb.fffffff7089387387de41acc3d4p+0)));
    try std.testing.expectEqual(-0x2.34eedcb0ecf028dccedaec8d0426p-76, lgamma(@as(f128, -0xb.fffffff7089387387de41acc3cp+0)));
    // try std.testing.expectEqual(0x4.ee50220bbf41cb2e91b43d75bc74p-76, lgamma(@as(f128, -0xb.fffffff7089387387de41acc4p+0)));
    // try std.testing.expectEqual(-0x6.1fd05fe315324a387d5380a1660cp+0, lgamma(@as(f128, -0xc.00001p+0)));
    // try std.testing.expectEqual(0xd.4b0a2023492b1c2bf31822109db8p-24, lgamma(@as(f128, -0xc.00000008f76cp+0)));
    // try std.testing.expectEqual(-0xf.b743a426163665b0453dbafc72ep-28, lgamma(@as(f128, -0xc.00000008f76c8p+0)));
    // try std.testing.expectEqual(0x2.6322ea559f93a0b65b314a996df8p-36, lgamma(@as(f128, -0xc.00000008f76c773p+0)));
    // try std.testing.expectEqual(-0x1.a29d91aa27903fb6055269a5275bp-32, lgamma(@as(f128, -0xc.00000008f76c774p+0)));
    try std.testing.expectEqual(0x6.4c596ec141406827148aa9bdeb64p-84, lgamma(@as(f128, -0xc.00000008f76c7731567c0f0250fp+0)));
    try std.testing.expectEqual(-0x7.fa2493c5665b67fadd26ecb71368p-84, lgamma(@as(f128, -0xc.00000008f76c7731567c0f0250f8p+0)));
    // try std.testing.expectEqual(0x1.b28f1dba88e582cc217e7cd29c4p-76, lgamma(@as(f128, -0xc.00000008f76c7731567c0f025p+0)));
    // try std.testing.expectEqual(-0x5.70afe388cae86544d684f7c84998p-76, lgamma(@as(f128, -0xc.00000008f76c7731567c0f0254p+0)));
    try std.testing.expectEqual(-0x8.b07093393f8bec5dcbeca94ad53p+0, lgamma(@as(f128, -0xc.fffffp+0)));
    // try std.testing.expectEqual(-0x7.316d88601881509658502a6f3f3p-20, lgamma(@as(f128, -0xc.ffffffff4f6d8p+0)));
    // try std.testing.expectEqual(0x4.67d7d4d0a160ff7dc6f636e473bp-20, lgamma(@as(f128, -0xc.ffffffff4f6ep+0)));
    try std.testing.expectEqual(-0x2.2c25e6e64d1da5ede86337c40edap-32, lgamma(@as(f128, -0xc.ffffffff4f6dcf6p+0)));
    try std.testing.expectEqual(0x1.50666d9a112317971ea00308e623p-28, lgamma(@as(f128, -0xc.ffffffff4f6dcf7p+0)));
    try std.testing.expectEqual(-0xb.5b581a4ac393dd116537fa8b602p-80, lgamma(@as(f128, -0xc.ffffffff4f6dcf617f97a5ffc75p+0)));
    // try std.testing.expectEqual(0x3.dee458b96deb245de3658790cfe8p-84, lgamma(@as(f128, -0xc.ffffffff4f6dcf617f97a5ffc758p+0)));
    try std.testing.expectEqual(-0x4.d8d27bc90c37033727acbba9d3dcp-72, lgamma(@as(f128, -0xc.ffffffff4f6dcf617f97a5ffc4p+0)));
    // try std.testing.expectEqual(0xf.3d0b422210244746edfa7b5d241p-76, lgamma(@as(f128, -0xc.ffffffff4f6dcf617f97a5ffc8p+0)));
    try std.testing.expectEqual(-0x8.b070e6845a6ce3384311f5033318p+0, lgamma(@as(f128, -0xd.00001p+0)));
    // try std.testing.expectEqual(0x4.679e61ad5162fc7e1c654d564528p-20, lgamma(@as(f128, -0xd.00000000b092p+0)));
    try std.testing.expectEqual(-0x7.31a6fbad0e0cc4117020643e69bcp-20, lgamma(@as(f128, -0xd.00000000b0928p+0)));
    try std.testing.expectEqual(0x1.16f33a7d23d6cb18bb112232c1d8p-28, lgamma(@as(f128, -0xd.00000000b092309p+0)));
    // try std.testing.expectEqual(-0x5.c35919086cfd4ecafbcfe5a84b98p-32, lgamma(@as(f128, -0xd.00000000b09230ap+0)));
    // try std.testing.expectEqual(0x5.a339fee9d14554c80472b7f2bbdp-80, lgamma(@as(f128, -0xd.00000000b092309c06683dd1b9p+0)));
    // try std.testing.expectEqual(-0x5.f60c613fd4481b8619b8d1b38718p-80, lgamma(@as(f128, -0xd.00000000b092309c06683dd1b908p+0)));
    // try std.testing.expectEqual(0x1.78cc06041e82f35e8cdf5c0cba5p-72, lgamma(@as(f128, -0xd.00000000b092309c06683dd1b8p+0)));
    try std.testing.expectEqual(-0x4.53d72a10b443c4c879ef8da38468p-72, lgamma(@as(f128, -0xd.00000000b092309c06683dd1bcp+0)));
    try std.testing.expectEqual(-0xb.5409d4efa4b70f8f3d8788779a88p+0, lgamma(@as(f128, -0xd.fffffp+0)));
    // try std.testing.expectEqual(-0x5.861824905c091e728232d794138p-16, lgamma(@as(f128, -0xd.fffffffff363p+0)));
    try std.testing.expectEqual(0x4.a000dfad124b37a42c08da284184p-16, lgamma(@as(f128, -0xd.fffffffff3638p+0)));
    try std.testing.expectEqual(-0xe.bcf83d656a15debaeee43e4b325p-28, lgamma(@as(f128, -0xd.fffffffff36345ap+0)));
    // try std.testing.expectEqual(0x5.8f42e4c2cdc7cbbccabf0a7808f4p-28, lgamma(@as(f128, -0xd.fffffffff36345bp+0)));
    try std.testing.expectEqual(-0x1.627c8836779854634351d0f7c6ddp-76, lgamma(@as(f128, -0xd.fffffffff36345ab9e184a3e09dp+0)));
    // try std.testing.expectEqual(0x8.c3a10bc6dbc5b0028a081a21a7ap-76, lgamma(@as(f128, -0xd.fffffffff36345ab9e184a3e09d8p+0)));
    try std.testing.expectEqual(-0x2.4e05300f9b5ae55348c3229a5cep-68, lgamma(@as(f128, -0xd.fffffffff36345ab9e184a3e08p+0)));
    try std.testing.expectEqual(0x2.c50999ef0e541cdfaf9dd2520952p-68, lgamma(@as(f128, -0xd.fffffffff36345ab9e184a3e0cp+0)));
    try std.testing.expectEqual(-0xb.540a2a83e42a4f8e47f4ba505p+0, lgamma(@as(f128, -0xe.00001p+0)));
    try std.testing.expectEqual(0x4.a0009c38d0a8285ae87c2fd3240cp-16, lgamma(@as(f128, -0xe.000000000c9c8p+0)));
    // try std.testing.expectEqual(-0x5.861868074a4e2955c5b8093665ap-16, lgamma(@as(f128, -0xe.000000000c9dp+0)));
    // try std.testing.expectEqual(0x5.8b0b8d2a481f4700368f7fdea0bcp-28, lgamma(@as(f128, -0xe.000000000c9cba5p+0)));
    // try std.testing.expectEqual(-0xe.c12f950349025aab8304d77f03d8p-28, lgamma(@as(f128, -0xe.000000000c9cba6p+0)));
    // try std.testing.expectEqual(0x1.fac1bf7cf1f74c5fcd608a5ca6dcp-76, lgamma(@as(f128, -0xe.000000000c9cba545e94e75ec57p+0)));
    // try std.testing.expectEqual(-0x8.2b5bd485baaaaf39fb9cba2a67bp-76, lgamma(@as(f128, -0xe.000000000c9cba545e94e75ec578p+0)));
    // try std.testing.expectEqual(0x1.d4d41257f7f712821d19d1e48499p-68, lgamma(@as(f128, -0xe.000000000c9cba545e94e75ec4p+0)));
    try std.testing.expectEqual(-0x3.3e3ab7a95e59eb4a8f10d25cfe24p-68, lgamma(@as(f128, -0xe.000000000c9cba545e94e75ec8p+0)));
    // try std.testing.expectEqual(-0xe.094c9b083ca94d01fbdb43c57ae8p+0, lgamma(@as(f128, -0xe.fffffp+0)));
    // try std.testing.expectEqual(-0x4.c8585a763b9d58036e94236507a8p-12, lgamma(@as(f128, -0xe.ffffffffff288p+0)));
    // try std.testing.expectEqual(0x4.bb5f60f986f8a8e0b908fc5bb77p-12, lgamma(@as(f128, -0xe.ffffffffff29p+0)));
    try std.testing.expectEqual(-0xe.beef09380560f8096fc599fed51p-28, lgamma(@as(f128, -0xe.ffffffffff28c06p+0)));
    try std.testing.expectEqual(0x1.21b8928708bc37b5ecc9dcb97281p-20, lgamma(@as(f128, -0xe.ffffffffff28c07p+0)));
    // try std.testing.expectEqual(-0x2.58262de2adbf5f56b3ba66632876p-72, lgamma(@as(f128, -0xe.ffffffffff28c060c6604ef3037p+0)));
    // try std.testing.expectEqual(0x7.2b958cdd26656fdeb1f0835f3cdp-72, lgamma(@as(f128, -0xe.ffffffffff28c060c6604ef30378p+0)));
    // try std.testing.expectEqual(-0x4.18f2d06c4fd5905fd98608342c4p-64, lgamma(@as(f128, -0xe.ffffffffff28c060c6604ef3p+0)));
    try std.testing.expectEqual(0xa.8eb0cf39a3cd732a08cf71f5fbbp-68, lgamma(@as(f128, -0xe.ffffffffff28c060c6604ef304p+0)));
    // try std.testing.expectEqual(-0xe.094cf2be9e3eaf232939b809f3p+0, lgamma(@as(f128, -0xf.00001p+0)));
    try std.testing.expectEqual(0x4.bb5f60afdcccb46b4f271d7625a4p-12, lgamma(@as(f128, -0xf.0000000000d7p+0)));
    try std.testing.expectEqual(-0x4.c8585ac011a47d4389869bd07ddp-12, lgamma(@as(f128, -0xf.0000000000d78p+0)));
    // try std.testing.expectEqual(0x1.21b848c7158f27a4dd8cba8a9fafp-20, lgamma(@as(f128, -0xf.0000000000d73f9p+0)));
    try std.testing.expectEqual(-0xe.bf38c930add7226ecefaf98e3218p-28, lgamma(@as(f128, -0xf.0000000000d73fap+0)));
    // try std.testing.expectEqual(0x7.088d5a8137b6f804702dc202fff8p-76, lgamma(@as(f128, -0xf.0000000000d73f9f399bd0e420f8p+0)));
    try std.testing.expectEqual(-0x9.1332e518185fc14a40bf2730213p-72, lgamma(@as(f128, -0xf.0000000000d73f9f399bd0e421p+0)));
    // try std.testing.expectEqual(0x1.27644472ed630658b61d0b0eae39p-64, lgamma(@as(f128, -0xf.0000000000d73f9f399bd0e42p+0)));
    try std.testing.expectEqual(-0x3.9a7998ed288a9206ceb26aa19f5ep-64, lgamma(@as(f128, -0xf.0000000000d73f9f399bd0e424p+0)));
    // try std.testing.expectEqual(-0x1.0cf14f9e783e6b3b12314bccff56p+4, lgamma(@as(f128, -0xf.fffffp+0)));
    try std.testing.expectEqual(-0xe.466b0623a18cfb084ac2ebacb15p-12, lgamma(@as(f128, -0xf.fffffffffff28p+0)));
    // try std.testing.expectEqual(0x8.c4f2f20afce33e1bd4e089a038e8p-8, lgamma(@as(f128, -0xf.fffffffffff3p+0)));
    // try std.testing.expectEqual(-0x7.318a3fab1e86e0b05917d7632f18p-20, lgamma(@as(f128, -0xf.fffffffffff28cp+0)));
    // try std.testing.expectEqual(0xb.d5eff885a06ba0727b7eafe8fc7p-20, lgamma(@as(f128, -0xf.fffffffffff28c1p+0)));
    // try std.testing.expectEqual(-0x8.8a5563410902f2fc7c7a69503978p-68, lgamma(@as(f128, -0xf.fffffffffff28c060c6621f512ep+0)));
    // try std.testing.expectEqual(0xf.966577ef42f59f4988b610f9e87p-72, lgamma(@as(f128, -0xf.fffffffffff28c060c6621f512e8p+0)));
    // try std.testing.expectEqual(-0x3.73e1cc804007163b134e47452d6cp-60, lgamma(@as(f128, -0xf.fffffffffff28c060c6621f51p+0)));
    try std.testing.expectEqual(0x1.4dfc10dfbe920fecd6fb77f829cp-60, lgamma(@as(f128, -0xf.fffffffffff28c060c6621f514p+0)));
    try std.testing.expectEqual(-0x1.180879870e33e355b67293d3944bp+4, lgamma(@as(f128, -0x1.000002p+4)));
    // try std.testing.expectEqual(0x8.c4f2f20ab3ff0ed275259026a5dp-8, lgamma(@as(f128, -0x1.000000000000dp+4)));
    // try std.testing.expectEqual(-0xa.33ca82bb399dc62af456a083b24p-8, lgamma(@as(f128, -0x1.000000000000ep+4)));
    // try std.testing.expectEqual(0x1.edd80cde02fd7df4f903a50896b6p-16, lgamma(@as(f128, -0x1.000000000000d73ep+4)));
    try std.testing.expectEqual(-0x7.318a4462081fae5c7fba91eb6ecp-20, lgamma(@as(f128, -0x1.000000000000d74p+4)));
    // try std.testing.expectEqual(0x1.286b0c2ff32e03dda8ad8c34d8d9p-64, lgamma(@as(f128, -0x1.000000000000d73f9f399da1424bp+4)));
    // try std.testing.expectEqual(-0x8.0c6b280d2bb28498300670555558p-72, lgamma(@as(f128, -0x1.000000000000d73f9f399da1424cp+4)));
    try std.testing.expectEqual(0x5.a5b6b02f80d768265f4c0ec03bap-60, lgamma(@as(f128, -0x1.000000000000d73f9f399da142p+4)));
    // try std.testing.expectEqual(-0x3.de050a9081f64b6dc9abec14f97ep-60, lgamma(@as(f128, -0x1.000000000000d73f9f399da1428p+4)));
    // try std.testing.expectEqual(-0x1.455d45b618e1f038dddeea5dfff7p+4, lgamma(@as(f128, -0x1.0ffffep+4)));
    // try std.testing.expectEqual(-0x3.be7ffe71389cc26835a85ecbcda4p-4, lgamma(@as(f128, -0x1.0ffffffffffffp+4)));
    // try std.testing.expectEqual(-0xc.57773dac63c828891bc6e3f41378p-16, lgamma(@as(f128, -0x1.0ffffffffffff356p+4)));
    // try std.testing.expectEqual(0x1.c19a5332b053694fd8a4b888ace1p-12, lgamma(@as(f128, -0x1.0ffffffffffff358p+4)));
    // try std.testing.expectEqual(-0x5.93f933dffa74012ca203ed3155bcp-64, lgamma(@as(f128, -0x1.0ffffffffffff3569c47e7a93e1cp+4)));
    // try std.testing.expectEqual(0xe.a3f5b8f8053066eb84078217fbp-64, lgamma(@as(f128, -0x1.0ffffffffffff3569c47e7a93e1dp+4)));
    // try std.testing.expectEqual(-0x2.3bb21b1b7ff06cdb8124ca6f017p-56, lgamma(@as(f128, -0x1.0ffffffffffff3569c47e7a93ep+4)));
    // try std.testing.expectEqual(0x7.e0455b507fe1e3888b7147b9896cp-56, lgamma(@as(f128, -0x1.0ffffffffffff3569c47e7a93e8p+4)));
    // try std.testing.expectEqual(-0x1.455d51292150d8b93e426f65c468p+4, lgamma(@as(f128, -0x1.100002p+4)));
    // try std.testing.expectEqual(-0x3.be7ffe7138f85aabacec61e0bb4ep-4, lgamma(@as(f128, -0x1.1000000000001p+4)));
    try std.testing.expectEqual(0x1.c19a533267df77f20158487aad6fp-12, lgamma(@as(f128, -0x1.1000000000000ca8p+4)));
    // try std.testing.expectEqual(-0xc.57773db0ebbe6eed7f15eafde5fp-16, lgamma(@as(f128, -0x1.1000000000000caap+4)));
    // try std.testing.expectEqual(0x8.94adc5a9656a5944f73f433c03cp-64, lgamma(@as(f128, -0x1.1000000000000ca963b818568887p+4)));
    // try std.testing.expectEqual(-0xb.a341272e9af13edf98c4eab5553p-64, lgamma(@as(f128, -0x1.1000000000000ca963b818568888p+4)));
    try std.testing.expectEqual(0x9.61c363f9167ebaefaf2d7ab0fa2p-60, lgamma(@as(f128, -0x1.1000000000000ca963b81856888p+4)));
    // try std.testing.expectEqual(-0x9.85db402c6ec5b3470f30b672e988p-56, lgamma(@as(f128, -0x1.1000000000000ca963b8185689p+4)));
    try std.testing.expectEqual(-0x1.739c3c0e7e3dc747f6c9173a7b13p+4, lgamma(@as(f128, -0x1.1ffffep+4)));
    try std.testing.expectEqual(-0x3.1fd7673485ba8a86b1e31b4b3ca6p+0, lgamma(@as(f128, -0x1.1ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x2.b80fd7d902af06e3995bb05523dcp-8, lgamma(@as(f128, -0x1.1fffffffffffff4ap+4)));
    // try std.testing.expectEqual(0x1.c19a53328e26a91c6bdd01b59dafp-12, lgamma(@as(f128, -0x1.1fffffffffffff4cp+4)));
    // try std.testing.expectEqual(-0xc.5d86cd624ca79f6dceb2b1d20bb8p-64, lgamma(@as(f128, -0x1.1fffffffffffff4bec3ce234132dp+4)));
    try std.testing.expectEqual(0x1.5f9145d9cdb2fbf24dacd8acd86ep-56, lgamma(@as(f128, -0x1.1fffffffffffff4bec3ce234132ep+4)));
    try std.testing.expectEqual(-0x4.005578030d2343c85effa7bae274p-52, lgamma(@as(f128, -0x1.1fffffffffffff4bec3ce23413p+4)));
    // try std.testing.expectEqual(0x7.5f20ed3672db03e0763d1b82b2ecp-52, lgamma(@as(f128, -0x1.1fffffffffffff4bec3ce234138p+4)));
    try std.testing.expectEqual(-0x1.739c47ba6a3ae8abe5a16e7d7a65p+4, lgamma(@as(f128, -0x1.200002p+4)));
    try std.testing.expectEqual(-0x3.1fd7673485c0607cb073cd43a7f2p+0, lgamma(@as(f128, -0x1.2000000000001p+4)));
    // try std.testing.expectEqual(0x1.c19a53328a0c38256e1fdf0a2bfp-12, lgamma(@as(f128, -0x1.20000000000000b4p+4)));
    try std.testing.expectEqual(-0x2.b80fd7d902f168b1c90998bee862p-8, lgamma(@as(f128, -0x1.20000000000000b6p+4)));
    // try std.testing.expectEqual(0x1.16355d66125b301ee0b5e5281079p-56, lgamma(@as(f128, -0x1.20000000000000b413c31dcbeca4p+4)));
    try std.testing.expectEqual(-0x5.5b96f411da52dc950d7091977238p-60, lgamma(@as(f128, -0x1.20000000000000b413c31dcbeca5p+4)));
    // try std.testing.expectEqual(0x3.443ca24e8d26da6bf3ab93ba4fd4p-52, lgamma(@as(f128, -0x1.20000000000000b413c31dcbec8p+4)));
    // try std.testing.expectEqual(-0x8.1b39c2eaf2da5837e0d9887f7b3p-52, lgamma(@as(f128, -0x1.20000000000000b413c31dcbedp+4)));
    try std.testing.expectEqual(-0x1.a2b8a7ff951d4cd8ff71bbc81688p+4, lgamma(@as(f128, -0x1.2ffffep+4)));
    // try std.testing.expectEqual(-0x6.119e27f51c200b4d7dd7ace1fa28p+0, lgamma(@as(f128, -0x1.2ffffffffffffp+4)));
    // try std.testing.expectEqual(-0xd.bb3fcdf10bfe34aa839c1b1a191p-8, lgamma(@as(f128, -0x1.2ffffffffffffff6p+4)));
    // try std.testing.expectEqual(0x2.b64afc1442844c492a44cf5aae0cp-4, lgamma(@as(f128, -0x1.2ffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x9.ce0b9a4828db14dcde92cf4d3f9p-56, lgamma(@as(f128, -0x1.2ffffffffffffff685b25cbf5f54p+4)));
    // try std.testing.expectEqual(0x1.134ad96206724f00837cae8366b3p-52, lgamma(@as(f128, -0x1.2ffffffffffffff685b25cbf5f55p+4)));
    try std.testing.expectEqual(-0x8.e6b2cf7c97411e8ac8b53cb0d9fp-48, lgamma(@as(f128, -0x1.2ffffffffffffff685b25cbf5fp+4)));
    // try std.testing.expectEqual(0x4.9aa9c8b7b0a1adc41049a8430efp-48, lgamma(@as(f128, -0x1.2ffffffffffffff685b25cbf5f8p+4)));
    try std.testing.expectEqual(-0x1.a2b8b3e16627e7804ccde008eedcp+4, lgamma(@as(f128, -0x1.300002p+4)));
    try std.testing.expectEqual(-0x6.119e27f51c25fc36032500898dep+0, lgamma(@as(f128, -0x1.3000000000001p+4)));
    // try std.testing.expectEqual(0x2.b64afc1442841cc1e61a64bd716ep-4, lgamma(@as(f128, -0x1.3000000000000008p+4)));
    // try std.testing.expectEqual(-0xd.bb3fcdf10c01eb3bd6ec6f62d56p-8, lgamma(@as(f128, -0x1.300000000000000ap+4)));
    // try std.testing.expectEqual(0xd.afc10a2f38dbffb02355aa0679p-56, lgamma(@as(f128, -0x1.30000000000000097a4da340a0abp+4)));
    try std.testing.expectEqual(-0xd.52f82639572464437b5ce2b7abp-56, lgamma(@as(f128, -0x1.30000000000000097a4da340a0acp+4)));
    try std.testing.expectEqual(0x4.9724dc2bbf7374cb291b72d5b15cp-48, lgamma(@as(f128, -0x1.30000000000000097a4da340a08p+4)));
    try std.testing.expectEqual(-0x8.ea37bc08886f870af40ddcbb742p-48, lgamma(@as(f128, -0x1.30000000000000097a4da340a1p+4)));
    try std.testing.expectEqual(-0x1.d2a72cdce34ac164fbae8c7684ddp+4, lgamma(@as(f128, -0x1.3ffffep+4)));
    try std.testing.expectEqual(-0x9.10867763989228882449ec5b3c48p+0, lgamma(@as(f128, -0x1.3ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.709f6fbd94aaf306ded2bd06724bp+0, lgamma(@as(f128, -0x1.3ffffffffffffffep+4)));
    try std.testing.expectEqual(-0xe.38f646c46c4d27208e77961a34f8p-52, lgamma(@as(f128, -0x1.3fffffffffffffff86af516ff7f7p+4)));
    try std.testing.expectEqual(0x1.38a7135be47b86f556df0964a769p-48, lgamma(@as(f128, -0x1.3fffffffffffffff86af516ff7f8p+4)));
    try std.testing.expectEqual(-0xf.c00e11277e57c2d191e36fd60bfp-44, lgamma(@as(f128, -0x1.3fffffffffffffff86af516ff78p+4)));
    try std.testing.expectEqual(0x1.21a5ad19d3f1ea0c24626e1975bep-44, lgamma(@as(f128, -0x1.3fffffffffffffff86af516ff8p+4)));
    try std.testing.expectEqual(-0x1.d2a738f1e7888f3f7c6994ba183ep+4, lgamma(@as(f128, -0x1.400002p+4)));
    // try std.testing.expectEqual(-0x9.108677639898330a4330d99c6998p+0, lgamma(@as(f128, -0x1.4000000000001p+4)));
    // try std.testing.expectEqual(-0x1.709f6fbd94aaf3c82f1699e41a71p+0, lgamma(@as(f128, -0x1.4000000000000002p+4)));
    try std.testing.expectEqual(0x1.3879456d6785f6aa67429896d1f9p-48, lgamma(@as(f128, -0x1.40000000000000007950ae900808p+4)));
    // try std.testing.expectEqual(-0xe.3bd325ac3ba631da0a5d89e4cd2p-52, lgamma(@as(f128, -0x1.40000000000000007950ae900809p+4)));
    // try std.testing.expectEqual(0x1.21a2d03aec229137c9799e44027p-44, lgamma(@as(f128, -0x1.40000000000000007950ae9008p+4)));
    try std.testing.expectEqual(-0xf.c010ee0666271eab2ddbb3221fdp-44, lgamma(@as(f128, -0x1.40000000000000007950ae90088p+4)));
    try std.testing.expectEqual(-0x2.035d89ed6121f85bdcd2763fe0bcp+4, lgamma(@as(f128, -0x1.4ffffep+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e683b14fbebdfbc5b28p+0, lgamma(@as(f128, -0x1.4ffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.7c05424b8a8111c2f3687fa4854p+0, lgamma(@as(f128, -0x1.4ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.c526944f697ba10afa8327f58ac2p-44, lgamma(@as(f128, -0x1.4ffffffffffffffffa391c4248c2p+4)));
    // try std.testing.expectEqual(0xf.fe0e8e74f374bf119b7af74a1518p-48, lgamma(@as(f128, -0x1.4ffffffffffffffffa391c4248c3p+4)));
    try std.testing.expectEqual(-0xb.89114dc66e395aed6a7298d3ee8p-40, lgamma(@as(f128, -0x1.4ffffffffffffffffa391c42488p+4)));
    // try std.testing.expectEqual(0xa.9f2a9bef4dc759d87f75188b7328p-40, lgamma(@as(f128, -0x1.4ffffffffffffffffa391c4249p+4)));
    try std.testing.expectEqual(-0x2.035d9633286bf6f969e3ff6bccfp+4, lgamma(@as(f128, -0x1.500002p+4)));
    try std.testing.expectEqual(-0xc.1bec49f18e6e5df8a0eb2e83a0d8p+0, lgamma(@as(f128, -0x1.5000000000001p+4)));
    // try std.testing.expectEqual(-0x4.7c05424b8a8112874fdd1f8e5e28p+0, lgamma(@as(f128, -0x1.5000000000000002p+4)));
    try std.testing.expectEqual(0xf.fe0c5746b70b02a08b297e590878p-48, lgamma(@as(f128, -0x1.500000000000000005c6e3bdb73dp+4)));
    try std.testing.expectEqual(-0x1.c526b7c24d423cd82e6be483ea53p-44, lgamma(@as(f128, -0x1.500000000000000005c6e3bdb73ep+4)));
    // try std.testing.expectEqual(0xa.9f2a99b81f8af0337368ac1fde78p-40, lgamma(@as(f128, -0x1.500000000000000005c6e3bdb7p+4)));
    // try std.testing.expectEqual(-0xb.89114ffd9c75c4c38d9c2d39f97p-40, lgamma(@as(f128, -0x1.500000000000000005c6e3bdb78p+4)));
    try std.testing.expectEqual(-0x2.34d272c496dc021c05680f598766p+4, lgamma(@as(f128, -0x1.5ffffep+4)));
    // try std.testing.expectEqual(-0xf.333ad8d94721201568ad5e5db99p+0, lgamma(@as(f128, -0x1.5ffffffffffffp+4)));
    // try std.testing.expectEqual(-0x7.9353d133433a0264d487158bb564p+0, lgamma(@as(f128, -0x1.5ffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.4f359fa4fdbfce6f5bb1209ad1ecp-44, lgamma(@as(f128, -0x1.5fffffffffffffffffbcc71a492p+4)));
    // try std.testing.expectEqual(0x3.59f6f230f3b480081844e49e805p-40, lgamma(@as(f128, -0x1.5fffffffffffffffffbcc71a4921p+4)));
    // try std.testing.expectEqual(-0x7.a523cdf44d8687ba0b959d73c248p-36, lgamma(@as(f128, -0x1.5fffffffffffffffffbcc71a49p+4)));
    // try std.testing.expectEqual(0x1.6d22e937415b6230eb30e8e04412p-32, lgamma(@as(f128, -0x1.5fffffffffffffffffbcc71a498p+4)));
    // try std.testing.expectEqual(-0x2.34d27f38e9c8e973c1260ebe82ecp+4, lgamma(@as(f128, -0x1.600002p+4)));
    // try std.testing.expectEqual(-0xf.333ad8d947275a3edf210a3c4518p+0, lgamma(@as(f128, -0x1.6000000000001p+4)));
    // try std.testing.expectEqual(-0x7.9353d133433a032c19b5e4013138p+0, lgamma(@as(f128, -0x1.6000000000000002p+4)));
    try std.testing.expectEqual(0x3.59f6f216ca01e747521a9ee9b19cp-40, lgamma(@as(f128, -0x1.6000000000000000004338e5b6dfp+4)));
    try std.testing.expectEqual(-0x7.4f35a14798e95a81f87ef25b691p-44, lgamma(@as(f128, -0x1.6000000000000000004338e5b6ep+4)));
    // try std.testing.expectEqual(0x1.6d22e9372731af984f6414c9de1p-32, lgamma(@as(f128, -0x1.6000000000000000004338e5b68p+4)));
    // try std.testing.expectEqual(-0x7.a523cdf5f021b146e5779a13f85cp-36, lgamma(@as(f128, -0x1.6000000000000000004338e5b7p+4)));
    try std.testing.expectEqual(-0x2.66fd6ea9f77b79a6b027c2a9dfa2p+4, lgamma(@as(f128, -0x1.6ffffep+4)));
    // try std.testing.expectEqual(-0x1.255ea98937d9f1616b540f71866cp+4, lgamma(@as(f128, -0x1.6ffffffffffffp+4)));
    try std.testing.expectEqual(-0xa.b60390ed79b80386211aae4662d8p+0, lgamma(@as(f128, -0x1.6ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x1.37f86d6de495a16d879a63d238cbp-36, lgamma(@as(f128, -0x1.6ffffffffffffffffffd13c97d9dp+4)));
    try std.testing.expectEqual(0x4.41786010c9a0650a03802fb0197cp-36, lgamma(@as(f128, -0x1.6ffffffffffffffffffd13c97d9ep+4)));
    // try std.testing.expectEqual(-0x9.ff9bfb198d9cba2c3e54cf6d9fdp-32, lgamma(@as(f128, -0x1.6ffffffffffffffffffd13c97d8p+4)));
    // try std.testing.expectEqual(0x2.1cbea72e0da562cf14d9766a28d6p-28, lgamma(@as(f128, -0x1.6ffffffffffffffffffd13c97ep+4)));
    try std.testing.expectEqual(-0x2.66fd7b4acff91314aecad54bcbe4p+4, lgamma(@as(f128, -0x1.700002p+4)));
    try std.testing.expectEqual(-0x1.255ea98937da56682f40dae18568p+4, lgamma(@as(f128, -0x1.7000000000001p+4)));
    try std.testing.expectEqual(-0xa.b60390ed79b804502ea287dd42d8p+0, lgamma(@as(f128, -0x1.7000000000000002p+4)));
    try std.testing.expectEqual(0x4.41786010b72c10946e58f94676c4p-36, lgamma(@as(f128, -0x1.70000000000000000002ec368262p+4)));
    try std.testing.expectEqual(-0x1.37f86d6df709f5e32312067aa839p-36, lgamma(@as(f128, -0x1.70000000000000000002ec368263p+4)));
    try std.testing.expectEqual(0x2.1cbea72e0d92ee7aa1af18a3c99p-28, lgamma(@as(f128, -0x1.70000000000000000002ec3682p+4)));
    try std.testing.expectEqual(-0x9.ff9bfb198ec3ff73a37e0dc9f9dp-32, lgamma(@as(f128, -0x1.70000000000000000002ec36828p+4)));
    // try std.testing.expectEqual(-0x2.99d6bd8dc68007801753da9a4214p+4, lgamma(@as(f128, -0x1.7ffffep+4)));
    // try std.testing.expectEqual(-0x1.5837f8825c33e21e60c5af48acdp+4, lgamma(@as(f128, -0x1.7ffffffffffffp+4)));
    try std.testing.expectEqual(-0xd.e398807fbf5719fecd8a010e1e98p+0, lgamma(@as(f128, -0x1.7ffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x4.862b2dc3253947bc30ce36b569f4p-32, lgamma(@as(f128, -0x1.7fffffffffffffffffffe0d30fe6p+4)));
    try std.testing.expectEqual(0x3.affe0676a937375e99944a595968p-32, lgamma(@as(f128, -0x1.7fffffffffffffffffffe0d30fe7p+4)));
    try std.testing.expectEqual(-0x3.4a1a90952a99eda0ccd8a7fc6f4p-24, lgamma(@as(f128, -0x1.7fffffffffffffffffffe0d30f8p+4)));
    // try std.testing.expectEqual(0xd.0fa0475b683650543df5fbc1be68p-28, lgamma(@as(f128, -0x1.7fffffffffffffffffffe0d31p+4)));
    // try std.testing.expectEqual(-0x2.99d6ca5949a84b98c0bae097d5dap+4, lgamma(@as(f128, -0x1.800002p+4)));
    try std.testing.expectEqual(-0x1.5837f8825c34487a7a07d00e012p+4, lgamma(@as(f128, -0x1.8000000000001p+4)));
    try std.testing.expectEqual(-0xd.e398807fbf571acb85bc854fa94p+0, lgamma(@as(f128, -0x1.8000000000000002p+4)));
    // try std.testing.expectEqual(0x3.affe0676a92ac03fb5e4485c0afcp-32, lgamma(@as(f128, -0x1.800000000000000000001f2cf019p+4)));
    // try std.testing.expectEqual(-0x4.862b2dc32545bedb14e494cbfa8p-32, lgamma(@as(f128, -0x1.800000000000000000001f2cf01ap+4)));
    try std.testing.expectEqual(0xd.0fa0475b683588e2505aeb8960d8p-28, lgamma(@as(f128, -0x1.800000000000000000001f2cfp+4)));
    try std.testing.expectEqual(-0x3.4a1a90952a99fa17ebe5870c962ap-24, lgamma(@as(f128, -0x1.800000000000000000001f2cf08p+4)));
    try std.testing.expectEqual(-0x2.cd57416926b9198c8d473083f362p+4, lgamma(@as(f128, -0x1.8ffffep+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374e485085aa667ac9ep+4, lgamma(@as(f128, -0x1.8ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.11ba0bf7d70fd875d44cb36bf4cp+4, lgamma(@as(f128, -0x1.8ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.55ecffc0812aac5347cab8b401a8p-28, lgamma(@as(f128, -0x1.8ffffffffffffffffffffec0c332p+4)));
    // try std.testing.expectEqual(0xa.7eb36524b0d30bc5c23e9be4e6d8p-28, lgamma(@as(f128, -0x1.8ffffffffffffffffffffec0c333p+4)));
    // try std.testing.expectEqual(-0x2.83dd0d761876f08e6bfe35926ae6p-20, lgamma(@as(f128, -0x1.8ffffffffffffffffffffec0c3p+4)));
    try std.testing.expectEqual(0x3.e6736a6ff2727a302becd4cbe6a8p-20, lgamma(@as(f128, -0x1.8ffffffffffffffffffffec0c38p+4)));
    try std.testing.expectEqual(-0x2.cd574e5d9fa3ed015fba57b06442p+4, lgamma(@as(f128, -0x1.900002p+4)));
    try std.testing.expectEqual(-0x1.8bb87c72374eaff44d01022165dfp+4, lgamma(@as(f128, -0x1.9000000000001p+4)));
    try std.testing.expectEqual(-0x1.11ba0bf7d70fd882c8c59e3f6994p+4, lgamma(@as(f128, -0x1.9000000000000002p+4)));
    // try std.testing.expectEqual(0xa.7eb36524b0d303b1e7123f36cfb8p-28, lgamma(@as(f128, -0x1.90000000000000000000013f3ccdp+4)));
    // try std.testing.expectEqual(-0x2.55ecffc0812ab46722fd8f9e8e34p-28, lgamma(@as(f128, -0x1.90000000000000000000013f3ccep+4)));
    try std.testing.expectEqual(0x3.e6736a6ff2727a2818139b3367e2p-20, lgamma(@as(f128, -0x1.90000000000000000000013f3c8p+4)));
    try std.testing.expectEqual(-0x2.83dd0d761876f0967fdaac49246p-20, lgamma(@as(f128, -0x1.90000000000000000000013f3dp+4)));
    // try std.testing.expectEqual(-0x3.01786b2b55b39354d0060d9af742p+4, lgamma(@as(f128, -0x1.9ffffep+4)));
    try std.testing.expectEqual(-0x1.bfd9a6481783e14ac56ba21bb97ap+4, lgamma(@as(f128, -0x1.9ffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.45db35cdb745720d8a3551830bbfp+4, lgamma(@as(f128, -0x1.9ffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.55ecffc0812ab034f847f8a4c53ap-28, lgamma(@as(f128, -0x1.9ffffffffffffffffffffff3b8bdp+4)));
    // try std.testing.expectEqual(0x1.4b426a4f71edda0a8f5218c17e64p-20, lgamma(@as(f128, -0x1.9ffffffffffffffffffffff3b8bep+4)));
    try std.testing.expectEqual(-0x4.f7eda0c3cb5f0e25a9dac51cd54cp-16, lgamma(@as(f128, -0x1.9ffffffffffffffffffffff3b88p+4)));
    // try std.testing.expectEqual(0x5.74d73977b83c2eb193bf5b3f712p-16, lgamma(@as(f128, -0x1.9ffffffffffffffffffffff3b9p+4)));
    // try std.testing.expectEqual(-0x3.0178784731148e2c18b47a300152p+4, lgamma(@as(f128, -0x1.a00002p+4)));
    try std.testing.expectEqual(-0x1.bfd9a64817844a29a07378d606b4p+4, lgamma(@as(f128, -0x1.a000000000001p+4)));
    try std.testing.expectEqual(-0x1.45db35cdb745721aa610b27de309p+4, lgamma(@as(f128, -0x1.a000000000000002p+4)));
    try std.testing.expectEqual(0x1.4b426a4f71edda0a3ed7e6f8630fp-20, lgamma(@as(f128, -0x1.a0000000000000000000000c4742p+4)));
    try std.testing.expectEqual(-0x2.55ecffc0812ab08572804fadcaa2p-28, lgamma(@as(f128, -0x1.a0000000000000000000000c4743p+4)));
    // try std.testing.expectEqual(0x5.74d73977b83c2eb18eb7d32c53e4p-16, lgamma(@as(f128, -0x1.a0000000000000000000000c47p+4)));
    // try std.testing.expectEqual(-0x4.f7eda0c3cb5f0e25aee2819f600cp-16, lgamma(@as(f128, -0x1.a0000000000000000000000c478p+4)));
    try std.testing.expectEqual(-0x3.36342a886637ea3d1ee94bbf39f4p+4, lgamma(@as(f128, -0x1.affffep+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8d0079500a7f922f3p+4, lgamma(@as(f128, -0x1.affffffffffffp+4)));
    try std.testing.expectEqual(-0x1.7a96f53dbe4e91d3b60397455b8bp+4, lgamma(@as(f128, -0x1.affffffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.d3e5be4d445c73cfbd7fad6d1cacp-20, lgamma(@as(f128, -0x1.afffffffffffffffffffffff8b95p+4)));
    try std.testing.expectEqual(0x1.b5b3f8628ba7d3d7a8c361f02124p-16, lgamma(@as(f128, -0x1.afffffffffffffffffffffff8b96p+4)));
    try std.testing.expectEqual(-0x2.ea6c2c67704b8ff10edfef0ae1bcp-12, lgamma(@as(f128, -0x1.afffffffffffffffffffffff8b8p+4)));
    // try std.testing.expectEqual(0xe.b396b51ee93d97f1fb62ed4b0cc8p-12, lgamma(@as(f128, -0x1.afffffffffffffffffffffff8cp+4)));
    try std.testing.expectEqual(-0x3.363437ca2ea26056c67a1202c958p+4, lgamma(@as(f128, -0x1.b00002p+4)));
    try std.testing.expectEqual(-0x1.f49565b81e8d6a87935e305f72eep+4, lgamma(@as(f128, -0x1.b000000000001p+4)));
    // try std.testing.expectEqual(-0x1.7a96f53dbe4e91e0f7cc01bb7533p+4, lgamma(@as(f128, -0x1.b000000000000002p+4)));
    try std.testing.expectEqual(0x1.b5b3f8628ba7d3d7a893278fb757p-16, lgamma(@as(f128, -0x1.b00000000000000000000000746ap+4)));
    try std.testing.expectEqual(-0x7.d3e5be4d445c73cfc0835a149dacp-20, lgamma(@as(f128, -0x1.b00000000000000000000000746bp+4)));
    // try std.testing.expectEqual(0xe.b396b51ee93d97f1fb5fec63a4a8p-12, lgamma(@as(f128, -0x1.b0000000000000000000000074p+4)));
    try std.testing.expectEqual(-0x2.ea6c2c67704b8ff10ee2f342bbf6p-12, lgamma(@as(f128, -0x1.b00000000000000000000000748p+4)));
    try std.testing.expectEqual(-0x3.6b84e02349a7940af2a134eb8688p+4, lgamma(@as(f128, -0x1.bffffep+4)));
    try std.testing.expectEqual(-0x2.29e61b654b214670ef8bad28fd7cp+4, lgamma(@as(f128, -0x1.bffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.afe7aaeaeae2d85d8c6032930547p+4, lgamma(@as(f128, -0x1.bffffffffffffffep+4)));
    try std.testing.expectEqual(-0x2.5dc72642d59f49efcf8837264d78p-12, lgamma(@as(f128, -0x1.bffffffffffffffffffffffffbd7p+4)));
    // try std.testing.expectEqual(0x1.7b435490f5313482a060e6848b86p-12, lgamma(@as(f128, -0x1.bffffffffffffffffffffffffbd8p+4)));
    // try std.testing.expectEqual(-0x1.4400f5be284e30221b069be275e3p-4, lgamma(@as(f128, -0x1.bffffffffffffffffffffffffb8p+4)));
    try std.testing.expectEqual(0x9.e6f812486e02e7c8d79623bc733p-8, lgamma(@as(f128, -0x1.bffffffffffffffffffffffffcp+4)));
    try std.testing.expectEqual(-0x3.6b84ed89a45b2eb6e36679911b66p+4, lgamma(@as(f128, -0x1.c00002p+4)));
    try std.testing.expectEqual(-0x2.29e61b654b21b1a3c52882888a5ep+4, lgamma(@as(f128, -0x1.c000000000001p+4)));
    try std.testing.expectEqual(-0x1.afe7aaeaeae2d86af2bae62db138p+4, lgamma(@as(f128, -0x1.c000000000000002p+4)));
    try std.testing.expectEqual(0x1.7b435490f5313482a060caabd708p-12, lgamma(@as(f128, -0x1.c000000000000000000000000428p+4)));
    try std.testing.expectEqual(-0x2.5dc72642d59f49efcf885305b524p-12, lgamma(@as(f128, -0x1.c000000000000000000000000429p+4)));
    // try std.testing.expectEqual(0x9.e6f812486e02e7c8d796220fa7d8p-8, lgamma(@as(f128, -0x1.c0000000000000000000000004p+4)));
    try std.testing.expectEqual(-0x1.4400f5be284e30221b069c009c2fp-4, lgamma(@as(f128, -0x1.c00000000000000000000000048p+4)));
    try std.testing.expectEqual(-0x3.a16551a93dea66ada032f329cee6p+4, lgamma(@as(f128, -0x1.cffffep+4)));
    try std.testing.expectEqual(-0x2.5fc68cfce71d835e7f01e235d532p+4, lgamma(@as(f128, -0x1.cffffffffffffp+4)));
    try std.testing.expectEqual(-0x1.e5c81c8286df15d847f9b7129f35p+4, lgamma(@as(f128, -0x1.cffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x2.104f809e32bb022582a7f432eeacp-8, lgamma(@as(f128, -0x1.cfffffffffffffffffffffffffdbp+4)));
    // try std.testing.expectEqual(0x4.f34f0176c4a5abfbae3a2dbc31e8p-8, lgamma(@as(f128, -0x1.cfffffffffffffffffffffffffdcp+4)));
    try std.testing.expectEqual(-0x1.3fc9d98001e767c1ad912d2039ap+0, lgamma(@as(f128, -0x1.cfffffffffffffffffffffffff8p+4)));
    try std.testing.expectEqual(-0x3.a1655f32e810c38e8832afeceb82p+4, lgamma(@as(f128, -0x1.d00002p+4)));
    try std.testing.expectEqual(-0x2.5fc68cfce71defabd034c93d1b76p+4, lgamma(@as(f128, -0x1.d000000000001p+4)));
    try std.testing.expectEqual(-0x1.e5c81c8286df15e5d1a3dd6f801dp+4, lgamma(@as(f128, -0x1.d000000000000002p+4)));
    // try std.testing.expectEqual(0x4.f34f0176c4a5abfbae3a2dacf70cp-8, lgamma(@as(f128, -0x1.d000000000000000000000000024p+4)));
    try std.testing.expectEqual(-0x2.104f809e32bb022582a7f44295d8p-8, lgamma(@as(f128, -0x1.d000000000000000000000000025p+4)));
    // try std.testing.expectEqual(-0x1.3fc9d98001e767c1ad912d206fc7p+0, lgamma(@as(f128, -0x1.d00000000000000000000000008p+4)));
    try std.testing.expectEqual(-0x3.d7d09f8a4486821f88b66a182d2cp+4, lgamma(@as(f128, -0x1.dffffep+4)));
    // try std.testing.expectEqual(-0x2.9631daeefecab8731b50a80d7dbp+4, lgamma(@as(f128, -0x1.dffffffffffffp+4)));
    try std.testing.expectEqual(-0x2.1c336a749e8c4b755bbff461bf2cp+4, lgamma(@as(f128, -0x1.dffffffffffffffep+4)));
    // try std.testing.expectEqual(-0x7.dd228d291dde78d3fafdd934d1bcp-4, lgamma(@as(f128, -0x1.dffffffffffffffffffffffffffep+4)));
    try std.testing.expectEqual(0x3.39fef253ff1921e8a33d604b6a06p-4, lgamma(@as(f128, -0x1.dfffffffffffffffffffffffffffp+4)));
    try std.testing.expectEqual(-0x4.a67eb8a17cbac193fb06132349e4p+0, lgamma(@as(f128, -0x1.dfffffffffffffffffffffffff8p+4)));
    // try std.testing.expectEqual(-0x3.d7d0ad3610cf012292e53b0205f2p+4, lgamma(@as(f128, -0x1.e00002p+4)));
    // try std.testing.expectEqual(-0x2.9631daeefecb25d17d94a025d504p+4, lgamma(@as(f128, -0x1.e000000000001p+4)));
    // try std.testing.expectEqual(-0x2.1c336a749e8c4b83078c3ce0c236p+4, lgamma(@as(f128, -0x1.e000000000000002p+4)));
    // try std.testing.expectEqual(0x3.39fef253ff1921e8a33d604b633p-4, lgamma(@as(f128, -0x1.e000000000000000000000000001p+4)));
    // try std.testing.expectEqual(-0x7.dd228d291dde78d3fafdd934df68p-4, lgamma(@as(f128, -0x1.e000000000000000000000000002p+4)));
    // try std.testing.expectEqual(-0x4.a67eb8a17cbac193fb0613238094p+0, lgamma(@as(f128, -0x1.e00000000000000000000000008p+4)));
    try std.testing.expectEqual(0x9.a81063e797803748580495bd2f48p+0, lgamma(@as(f128, 0x8.8d2d5p+0)));
    try std.testing.expectEqual(0x3.2125f40f9a1beba2b9f1959dbd98p+56, lgamma(@as(f128, 0x1.6a324ap+52)));
    try std.testing.expectEqual(0xb.70d4369f5b4c5572c84c32a22198p+0, lgamma(@as(f128, 0x9.62f59p+0)));
    try std.testing.expectEqual(0xe.b6cd62d45ad40dd2814b1697eb7p+0, lgamma(@as(f128, 0xa.d55d7p+0)));
    try std.testing.expectEqual(0xe.b6cd3d7503be73b09b5064553898p+0, lgamma(@as(f128, 0xa.d55d6p+0)));
    try std.testing.expectEqual(0xe.b6cd57db84c9ef437a5fd131a98p+0, lgamma(@as(f128, 0xa.d55d6b4d78e28p+0)));
    try std.testing.expectEqual(0xa.41afffa8a98e8455472818ee0938p+0, lgamma(@as(f128, 0x8.d6315p+0)));
    try std.testing.expectEqual(0xf.8842748a38e7a706e0144479dfc8p+0, lgamma(@as(f128, 0xb.2e679p+0)));
    try std.testing.expectEqual(0xf.1d4fd446695d45f71085f9be1868p+0, lgamma(@as(f128, 0xb.01191p+0)));
    try std.testing.expectEqual(0xf.76b5167078375bfcf413bd552c88p+0, lgamma(@as(f128, 0xb.26fdap+0)));
    try std.testing.expectEqual(0xf.cbb4eb9c9f4ddef22be7eb70eddp+0, lgamma(@as(f128, 0xb.4ad0ap+0)));
    // try std.testing.expectEqual(0xe.0ed26f91598df34bb14f20fb465p+24, lgamma(@as(f128, 0xe.7a678p+20)));
    // try std.testing.expectEqual(0x1.d9db4ca962b419e05bba7e38076bp+0, lgamma(@as(f128, -0x2.dea4ccp-4)));
    try std.testing.expectEqual(0x1.da47d6051ae6bf5e4dbd9b3a4acp+0, lgamma(@as(f128, -0x2.dd306p-4)));
    // try std.testing.expectEqual(0xf.f273df313425f4e361e154f408f8p-4, lgamma(@as(f128, -0x1.bdc8bp+0)));
    // try std.testing.expectEqual(0x1.950848252d48c05ac1f462baa5b6p+0, lgamma(@as(f128, -0x4.0a82e8p-4)));
    try std.testing.expectEqual(0xf.cc00043a75099f3d7c46acf58bb8p-4, lgamma(@as(f128, -0x1.bca67ap+0)));
    // try std.testing.expectEqual(-0xb.a18b329b453f2e72c0bc36a2e3dp-4, lgamma(@as(f128, -0x3.464468p+0)));
    // try std.testing.expectEqual(-0xb.a18c341739da8b29bdd8519d4f9p-4, lgamma(@as(f128, -0x3.46446cp+0)));
    // try std.testing.expectEqual(-0xb.a18c21a49016c028c0b54f8a1f5p-4, lgamma(@as(f128, -0x3.46446bb6a23aap+0)));
    // try std.testing.expectEqual(-0xe.aa75345fa640642f79f7d11a867p-8, lgamma(@as(f128, -0x3.f3d2c4p+0)));
    // try std.testing.expectEqual(-0xe.aa27b7e3f86d49a3ece042406398p-8, lgamma(@as(f128, -0x3.f3d2c8p+0)));
    // try std.testing.expectEqual(-0xe.aa7484b49666212f34177cf52dd8p-8, lgamma(@as(f128, -0x3.f3d2c40911814p+0)));
}
