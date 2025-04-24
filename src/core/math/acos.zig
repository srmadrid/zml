const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const uasncs = @import("uasncs.zig");
const root_tbl = @import("root_tbl.zig");
const powtwo_tbl = @import("powtwo_tbl.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn acos(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return acos(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, acos32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_acosf.c
                    return acos32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_asin.c
                    return acos64(x);
                },
                f80 => return cast(f80, acos128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_acosl.c
                    return acos128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn as_special(x: f32) f32 {
    const pih: f32 = 0x1.921fb6p+1;
    const pil: f32 = -0x1p-24;
    const t: u32 = @bitCast(x);
    if (t == (0x7f << 23))
        return 0; // x=1

    if (t == (0x17f << 23))
        return pih + pil; // x=-1

    const ax: u32 = t << 1;
    if (ax > (0xff << 24))
        return x + x; // nan

    return (0 - 0) / (0 - 0);
}

inline fn poly12(z: f64, c: *const [12]f64) f64 {
    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    var c0: f64 = c[0] + z * c[1];
    const c2: f64 = c[2] + z * c[3];
    var c4: f64 = c[4] + z * c[5];
    const c6: f64 = c[6] + z * c[7];
    var c8: f64 = c[8] + z * c[9];
    const c10: f64 = c[10] + z * c[11];
    c0 += c2 * z2;
    c4 += c6 * z2;
    c8 += z2 * c10;
    c0 += z4 * (c4 + z4 * c8);
    return c0;
}

fn acos32(x: f32) f32 {
    const pi2: f64 = 0x1.921fb54442d18p+0;
    const o: [2]f64 = .{ 0, 0x1.921fb54442d18p+1 };
    const xs: f64 = cast(f64, x, .{});

    const t: u32 = @bitCast(x);
    const ax: u32 = t << 1;
    if (ax >= 0x7f << 24) {
        @branchHint(.unlikely);
        return as_special(x);
    }

    if (ax < 0x7ec2a1dc) { // |x| < 0x1.c2a1dcp-1
        @branchHint(.likely);
        const b: [16]f64 = .{
            0x1.fffffffd9ccb8p-1,  0x1.5555c94838007p-3,  0x1.32ded4b7c20fap-4,
            0x1.8566df703309ep-5,  -0x1.980c959bec9a3p-6, 0x1.56fbb04998344p-1,
            -0x1.403d8e4c49f52p+2, 0x1.b06c3e9f311eap+4,  -0x1.9ea97c4e2c21fp+6,
            0x1.200b8261cc61bp+8,  -0x1.2274c2799a5c7p+9, 0x1.a558a59cc19d3p+9,
            -0x1.aca4b6a529ffp+9,  0x1.228744703f813p+9,  -0x1.d7dbb0b322228p+7,
            0x1.5c2018c0c0105p+5,
        };
        if (ax <= 0x40000000) { // |x| < 2^-63
            @branchHint(.unlikely);
            return cast(f32, pi2, .{});
        }

        const z: f64 = xs;
        const z2: f64 = z * z;
        const z4: f64 = z2 * z2;
        const z8: f64 = z4 * z4;
        const z16: f64 = z8 * z8;
        const r: f64 = z * ((((b[0] + z2 * b[1]) + z4 * (b[2] + z2 * b[3])) + z8 * ((b[4] + z2 * b[5]) + z4 * (b[6] + z2 * b[7]))) + z16 * (((b[8] + z2 * b[9]) + z4 * (b[10] + z2 * b[11])) + z8 * ((b[12] + z2 * b[13]) + z4 * (b[14] + z2 * b[15]))));
        const ub: f32 = cast(f32, 0x1.921fb54574191p+0 - r, .{});
        const lb: f32 = cast(f32, 0x1.921fb543118ap+0 - r, .{});
        if (ub == lb)
            return ub;
    }
    if (ax < (0x7e << 24)) { // accurate path
        const c: [12]f64 = .{
            0x1.555555555529cp-3, 0x1.333333337e0ddp-4,  0x1.6db6db3b4465ep-5,
            0x1.f1c72e13ac306p-6, 0x1.6e89cebe06bc4p-6,  0x1.1c6dcf5289094p-6,
            0x1.c6dbbcc7c6315p-7, 0x1.8f8dc2615e996p-7,  0x1.a5833b7bf15e8p-8,
            0x1.43f44ace1665cp-6, -0x1.0fb17df881c73p-6, 0x1.07520c026b2d6p-5,
        };
        if (t == 0x328885a3)
            return 0x1.921fb6p+0 + 0x1p-25;

        if (t == 0x39826222)
            return 0x1.920f6ap+0 + 0x1p-25;

        const x2: f64 = xs * xs;
        const r: f64 = (pi2 - xs) - (xs * x2) * poly12(x2, &c);
        return cast(f32, r, .{});
    } else {
        const c: [12]f64 = .{ 0x1.6a09e667f3bcbp+0, 0x1.e2b7dddff2db9p-4, 0x1.b27247ab42dbcp-6, 0x1.02995cc4e0744p-7, 0x1.5ffb0276ec8eap-9, 0x1.033885a928decp-10, 0x1.911f2be23f8c7p-12, 0x1.4c3c55d2437fdp-13, 0x1.af477e1d7b461p-15, 0x1.abd6bdff67dcbp-15, -0x1.1717e86d0fa28p-16, 0x1.6ff526de46023p-16 };
        const bx: f64 = math.abs(xs);
        const z: f64 = 1 - bx;
        const s: f64 = math.copysign(math.sqrt(z), xs);
        const r: f64 = o[t >> 31] + s * poly12(z, &c);
        return cast(f32, r, .{});
    }
}

// acos with max ULP of ~0.523 based on random sampling.
fn acos64(x: f64) f64 {
    const u: [2]i32 = @bitCast(x);
    const m: i32 = u[uasncs.HIGH_HALF];
    var k: i32 = 0x7fffffff & m;
    if (k < 0x3c880000) { //-------------------  |x|<2.77556*10^-17 ----------------------
        return @as(f64, @bitCast(uasncs.hp0));
    } else if (k < 0x3fc00000) { //-----------------  2.77556*10^-17 <= |x| < 2^-3 --------------
        const x2: f64 = x * x;
        const t: f64 = (((((uasncs.f6 * x2 + uasncs.f5) * x2 + uasncs.f4) * x2 + uasncs.f3) * x2 + uasncs.f2) * x2 + uasncs.f1) * (x2 * x);
        const r: f64 = @as(f64, @bitCast(uasncs.hp0)) - x;
        const cor: f64 = (((@as(f64, @bitCast(uasncs.hp0)) - r) - x) + @as(f64, @bitCast(uasncs.hp1))) - t;
        const res: f64 = r + cor;
        // Max ULP is 0.502.
        return res;
    } else if (k < 0x3fe00000) { //----------------------  0.125 <= |x| < 0.5 --------------------
        var n: i32 = undefined;
        if (k < 0x3fd00000) {
            n = 11 * ((k & 0x000fffff) >> 15);
        } else {
            n = 11 * ((k & 0x000fffff) >> 14) + 352;
        }

        var xx: f64 = undefined;
        if (m > 0) {
            xx = x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        } else {
            xx = -x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        }

        var t: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 1)])) * xx;
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) +
            xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)]));
        t += p;
        const y: f64 = if (m > 0) (@as(f64, @bitCast(uasncs.hp0)) - @as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)]))) else (@as(f64, @bitCast(uasncs.hp0)) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)])));
        t = if (m > 0) (@as(f64, @bitCast(uasncs.hp1)) - t) else (@as(f64, @bitCast(uasncs.hp1)) + t);
        const res: f64 = y + t;
        // Max ULP is 0.51.
        return res;
    } else if (k < 0x3fe80000) { //--------------------------- 0.5 <= |x| < 0.75 ---------------------
        const n: i32 = 1056 + ((k & 0x000fe000) >> 11) * 3;

        var xx: f64 = undefined;
        if (m > 0) {
            xx = x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        } else {
            xx = -x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        }

        var t: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 1)])) * xx;
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) +
            xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])) +
                xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)]))))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)]));
        t += p;
        const y: f64 = if (m > 0) (@as(f64, @bitCast(uasncs.hp0)) - @as(f64, @bitCast(uasncs.asncs[@intCast(n + 9)]))) else (@as(f64, @bitCast(uasncs.hp0)) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 9)])));
        t = if (m > 0) (@as(f64, @bitCast(uasncs.hp1)) - t) else (@as(f64, @bitCast(uasncs.hp1)) + t);
        const res: f64 = y + t;
        // Max ULP is 0.523 based on random sampling.
        return res;
    } else if (k < 0x3fed8000) { //------------------------- 0.75 <= |x| < 0.921875 -------------
        const n: i32 = 992 + ((k & 0x000fe000) >> 13) * 13;

        var xx: f64 = undefined;
        if (m > 0) {
            xx = x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        } else {
            xx = -x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        }

        var t: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 1)])) * xx;
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) +
            xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)])) +
                xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)])))))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 9)]));
        t += p;
        const y: f64 = if (m > 0) (@as(f64, @bitCast(uasncs.hp0)) - @as(f64, @bitCast(uasncs.asncs[@intCast(n + 10)]))) else (@as(f64, @bitCast(uasncs.hp0)) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 10)])));
        t = if (m > 0) (@as(f64, @bitCast(uasncs.hp1)) - t) else (@as(f64, @bitCast(uasncs.hp1)) + t);
        const res: f64 = y + t;
        // Max ULP is 0.523 based on random sampling.
        return res;
    } else if (k < 0x3fee8000) { //-------------------0.921875 <= |x| < 0.953125 ------------------
        const n: i32 = 884 + ((k & 0x000fe000) >> 13) * 14;

        var xx: f64 = undefined;
        if (m > 0) {
            xx = x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        } else {
            xx = -x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        }

        var t: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 1)])) * xx;
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) +
            xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)])) +
                xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 9)]))))))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 10)]));
        t += p;
        const y: f64 = if (m > 0) (@as(f64, @bitCast(uasncs.hp0)) - @as(f64, @bitCast(uasncs.asncs[@intCast(n + 11)]))) else (@as(f64, @bitCast(uasncs.hp0)) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 11)])));
        t = if (m > 0) (@as(f64, @bitCast(uasncs.hp1)) - t) else (@as(f64, @bitCast(uasncs.hp1)) + t);
        const res: f64 = y + t;
        // Max ULP is 0.523 based on random sampling.
        return res;
    } else if (k < 0x3fef0000) { //--------------------0.953125 <= |x| < 0.96875 ----------------
        const n: i32 = 768 + ((k & 0x000fe000) >> 13) * 15;

        var xx: f64 = undefined;
        if (m > 0) {
            xx = x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        } else {
            xx = -x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        }

        var t: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 1)])) * xx;
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) +
            xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 9)])) +
                xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 10)])))))))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 11)]));
        t += p;
        const y: f64 = if (m > 0) (@as(f64, @bitCast(uasncs.hp0)) - @as(f64, @bitCast(uasncs.asncs[@intCast(n + 12)]))) else (@as(f64, @bitCast(uasncs.hp0)) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 12)])));
        t = if (m > 0) (@as(f64, @bitCast(uasncs.hp1)) - t) else (@as(f64, @bitCast(uasncs.hp1)) + t);
        const res: f64 = y + t;
        // Max ULP is 0.523 based on random sampling.
        return res;
    } else if (k < 0x3ff00000) { //-----------------0.96875 <= |x| < 1 ---------------------------
        const z: f64 = 0.5 * (if (m > 0) (1 - x) else (1 + x));
        const v: [2]i32 = @bitCast(z);
        k = v[uasncs.HIGH_HALF];
        var t: f64 = root_tbl.inroot[@intCast((k & 0x001fffff) >> 14)] * powtwo_tbl.powtwo[@intCast(511 - (k >> 21))];
        const r: f64 = 1 - t * t * z;
        t = t * (uasncs.rt0 + r * (uasncs.rt1 + r * (uasncs.rt2 + r * uasncs.rt3)));
        const c: f64 = t * z;
        t = c * (1.5 - 0.5 * t * c);
        const y: f64 = (uasncs.t27 * c + c) - uasncs.t27 * c;
        const cc: f64 = (z - y * y) / (t + y);
        const p: f64 = (((((uasncs.f6 * z + uasncs.f5) * z + uasncs.f4) * z + uasncs.f3) * z + uasncs.f2) * z + uasncs.f1) * z;
        if (m < 0) {
            const cor: f64 = (@as(f64, @bitCast(uasncs.hp1)) - cc) - (y + cc) * p;
            const res1: f64 = @as(f64, @bitCast(uasncs.hp0)) - y;
            const res: f64 = res1 + cor;
            // Max ULP is 0.501.
            return (res + res);
        } else {
            const cor: f64 = cc + p * (y + cc);
            const res: f64 = y + cor;
            // Max ULP is 0.515.
            return (res + res);
        }
    } else if (k == 0x3ff00000 and u[uasncs.LOW_HALF] == 0) { //---------------------------- |x|>=1 -----------------------
        return if (m > 0) 0 else 2 * @as(f64, @bitCast(uasncs.hp0));
    } else {
        return (x - x) / (x - x);
    }
}

fn acos128(x: f128) f128 {
    const pio2_hi: f128 = 1.5707963267948966192313216916397514420986;
    const pio2_lo: f128 = 4.3359050650618905123985220130216759843812e-35;
    // acos(0.5625 + x): f128 = acos(0.5625) + x rS(x) / sS(x)
    // -0.0625 <= x <= 0.0625
    // peak relative error 3.3e-35
    const rS0: f128 = 5.619049346208901520945464704848780243887e0;
    const rS1: f128 = -4.460504162777731472539175700169871920352e1;
    const rS2: f128 = 1.317669505315409261479577040530751477488e2;
    const rS3: f128 = -1.626532582423661989632442410808596009227e2;
    const rS4: f128 = 3.144806644195158614904369445440583873264e1;
    const rS5: f128 = 9.806674443470740708765165604769099559553e1;
    const rS6: f128 = -5.708468492052010816555762842394927806920e1;
    const rS7: f128 = -1.396540499232262112248553357962639431922e1;
    const rS8: f128 = 1.126243289311910363001762058295832610344e1;
    const rS9: f128 = 4.956179821329901954211277873774472383512e-1;
    const rS10: f128 = -3.313227657082367169241333738391762525780e-1;
    const sS0: f128 = -4.645814742084009935700221277307007679325e0;
    const sS1: f128 = 3.879074822457694323970438316317961918430e1;
    const sS2: f128 = -1.221986588013474694623973554726201001066e2;
    const sS3: f128 = 1.658821150347718105012079876756201905822e2;
    const sS4: f128 = -4.804379630977558197953176474426239748977e1;
    const sS5: f128 = -1.004296417397316948114344573811562952793e2;
    const sS6: f128 = 7.530281592861320234941101403870010111138e1;
    const sS7: f128 = 1.270735595411673647119592092304357226607e1;
    const sS8: f128 = -1.815144839646376500705105967064792930282e1;
    const sS9: f128 = -7.821597334910963922204235247786840828217e-2;
    const acosr5625: f128 = 9.7338991014954640492751132535550279812151e-1;
    const pimacosr5625: f128 = 2.1682027434402468335351320579240000860757e0;
    // acos(0.4375 + x): f128 = acos(0.4375) + x rS(x) / sS(x)
    // -0.0625 <= x <= 0.0625
    // peak relative error 2.1e-35
    const P0: f128 = 2.177690192235413635229046633751390484892e0;
    const P1: f128 = -2.848698225706605746657192566166142909573e1;
    const P2: f128 = 1.040076477655245590871244795403659880304e2;
    const P3: f128 = -1.400087608918906358323551402881238180553e2;
    const P4: f128 = 2.221047917671449176051896400503615543757e1;
    const P5: f128 = 9.643714856395587663736110523917499638702e1;
    const P6: f128 = -5.158406639829833829027457284942389079196e1;
    const P7: f128 = -1.578651828337585944715290382181219741813e1;
    const P8: f128 = 1.093632715903802870546857764647931045906e1;
    const P9: f128 = 5.448925479898460003048760932274085300103e-1;
    const P10: f128 = -3.315886001095605268470690485170092986337e-1;
    const Q0: f128 = -1.958219113487162405143608843774587557016e0;
    const Q1: f128 = 2.614577866876185080678907676023269360520e1;
    const Q2: f128 = -9.990858606464150981009763389881793660938e1;
    const Q3: f128 = 1.443958741356995763628660823395334281596e2;
    const Q4: f128 = -3.206441012484232867657763518369723873129e1;
    const Q5: f128 = -1.048560885341833443564920145642588991492e2;
    const Q6: f128 = 6.745883931909770880159915641984874746358e1;
    const Q7: f128 = 1.806809656342804436118449982647641392951e1;
    const Q8: f128 = -1.770150690652438294290020775359580915464e1;
    const Q9: f128 = -5.659156469628629327045433069052560211164e-1;
    const acosr4375: f128 = 1.1179797320499710475919903296900511518755e0;
    const pimacosr4375: f128 = 2.0236129215398221908706530535894517323217e0;
    // asin(x): f128 = x + x^3 pS(x^2) / qS(x^2)
    // 0 <= x <= 0.5
    // peak relative error 1.9e-35  */
    const pS0: f128 = -8.358099012470680544198472400254596543711e2;
    const pS1: f128 = 3.674973957689619490312782828051860366493e3;
    const pS2: f128 = -6.730729094812979665807581609853656623219e3;
    const pS3: f128 = 6.643843795209060298375552684423454077633e3;
    const pS4: f128 = -3.817341990928606692235481812252049415993e3;
    const pS5: f128 = 1.284635388402653715636722822195716476156e3;
    const pS6: f128 = -2.410736125231549204856567737329112037867e2;
    const pS7: f128 = 2.219191969382402856557594215833622156220e1;
    const pS8: f128 = -7.249056260830627156600112195061001036533e-1;
    const pS9: f128 = 1.055923570937755300061509030361395604448e-3;
    const qS0: f128 = -5.014859407482408326519083440151745519205e3;
    const qS1: f128 = 2.430653047950480068881028451580393430537e4;
    const qS2: f128 = -4.997904737193653607449250593976069726962e4;
    const qS3: f128 = 5.675712336110456923807959930107347511086e4;
    const qS4: f128 = -3.881523118339661268482937768522572588022e4;
    const qS5: f128 = 1.634202194895541569749717032234510811216e4;
    const qS6: f128 = -4.151452662440709301601820849901296953752e3;
    const qS7: f128 = 5.956050864057192019085175976175695342168e2;
    const qS8: f128 = -4.175375777334867025769346564600396877176e1;

    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const sign: u32 = u.w0;
    const ix: i32 = @as(i32, @bitCast(sign)) & 0x7fffffff;
    u.w0 = @bitCast(ix); // |x|
    if (ix >= 0x3fff0000) { // |x| >= 1
        if (ix == 0x3fff0000 and (u.w1 | u.w2 | u.w3) == 0) { // |x| == 1
            if ((sign & 0x80000000) == 0) {
                return 0; // acos(1) = 0
            } else {
                return (2 * pio2_hi) + (2 * pio2_lo); // acos(-1)= pi
            }
        }
        return (x - x) / (x - x); // acos(|x| > 1) is NaN
    } else if (ix < 0x3ffe0000) { // |x| < 0.5
        if (ix < 0x3f8e0000) { // |x| < 2**-113
            return pio2_hi + pio2_lo;
        }

        if (ix < 0x3ffde000) { // |x| < .4375
            // Arcsine of x.
            var z: f128 = x * x;
            const p: f128 = (((((((((pS9 * z + pS8) * z + pS7) * z + pS6) * z + pS5) * z + pS4) * z + pS3) * z + pS2) * z + pS1) * z + pS0) * z;
            const q: f128 = ((((((((z + qS8) * z + qS7) * z + qS6) * z + qS5) * z + qS4) * z + qS3) * z + qS2) * z + qS1) * z + qS0;
            const r: f128 = x + x * p / q;
            z = pio2_hi - (r - pio2_lo);
            return z;
        }

        // .4375 <= |x| < .5
        const t: f128 = @as(f128, @bitCast(u)) - 0.4375;
        const p: f128 = ((((((((((P10 * t + P9) * t + P8) * t + P7) * t + P6) * t + P5) * t + P4) * t + P3) * t + P2) * t + P1) * t + P0) * t;

        const q: f128 = (((((((((t + Q9) * t + Q8) * t + Q7) * t + Q6) * t + Q5) * t + Q4) * t + Q3) * t + Q2) * t + Q1) * t + Q0;
        var r: f128 = p / q;
        if ((sign & 0x80000000) != 0) {
            r = pimacosr4375 - r;
        } else {
            r = acosr4375 + r;
        }

        return r;
    } else if (ix < 0x3ffe4000) { // |x| < 0.625
        const t: f128 = @as(f128, @bitCast(u)) - 0.5625;
        const p: f128 = ((((((((((rS10 * t + rS9) * t + rS8) * t + rS7) * t + rS6) * t + rS5) * t + rS4) * t + rS3) * t + rS2) * t + rS1) * t + rS0) * t;

        const q: f128 = (((((((((t + sS9) * t + sS8) * t + sS7) * t + sS6) * t + sS5) * t + sS4) * t + sS3) * t + sS2) * t + sS1) * t + sS0;
        if ((sign & 0x80000000) != 0) {
            return pimacosr5625 - p / q;
        } else {
            return acosr5625 + p / q;
        }
    } else { // |x| >= .625
        const z: f128 = (1 - @as(f128, @bitCast(u))) * 0.5;
        const s: f128 = math.sqrt(z);
        // Compute an extended precision square root from
        // the Newton iteration  s -> 0.5 * (s + z / s).
        // The change w from s to the improved value is
        // w = 0.5 * (s + z / s) - s  = (s^2 + z)/2s - s = (z - s^2)/2s.
        // Express s = f1 + f2 where f1 * f1 is exactly representable.
        // w = (z - s^2)/2s = (z - f1^2 - 2 f1 f2 - f2^2)/2s .
        // s + w has extended precision.
        u = @bitCast(s);
        u.w2 = 0;
        u.w3 = 0;
        const f2: f128 = s - @as(f128, @bitCast(u));
        var w: f128 = z - @as(f128, @bitCast(u)) * @as(f128, @bitCast(u));
        w = w - 2 * @as(f128, @bitCast(u)) * f2;
        w = w - f2 * f2;
        w = w / (2 * s);
        // Arcsine of s.
        const p: f128 = (((((((((pS9 * z + pS8) * z + pS7) * z + pS6) * z + pS5) * z + pS4) * z + pS3) * z + pS2) * z + pS1) * z + pS0) * z;
        const q: f128 = ((((((((z + qS8) * z + qS7) * z + qS6) * z + qS5) * z + qS4) * z + qS3) * z + qS2) * z + qS1) * z + qS0;
        const r: f128 = s + (w + s * p / q);

        if ((sign & 0x80000000) != 0) {
            w = pio2_hi + (pio2_lo - r);
        } else {
            w = r;
        }

        return 2 * w;
    }
}

test acos {
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x0p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x3.243f6cp+0, acos(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x1.0c1524p+0, acos(@as(f32, 0x8p-4)));
    try std.testing.expectEqual(0x2.182a48p+0, acos(@as(f32, -0x8p-4)));
    try std.testing.expectEqual(0xb.9051dp-4, acos(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x1.70ef56p-56)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x1.70ef54p-56)));
    try std.testing.expectEqual(0x1.821d0ap+0, acos(@as(f32, 0x1p-4)));
    try std.testing.expectEqual(0x1.6a09e6p-12, acos(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c8p+0, acos(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6p-12, acos(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c8p+0, acos(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6cp+0, acos(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6p-12, acos(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c8p+0, acos(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6cp+0, acos(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6p-12, acos(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c8p+0, acos(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6cp+0, acos(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6p-12, acos(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c8p+0, acos(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6cp+0, acos(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6p-12, acos(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c8p+0, acos(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6cp+0, acos(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6p-12, acos(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c8p+0, acos(@as(f32, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6cp+0, acos(@as(f32, -0x1p+0)));
    try std.testing.expectEqual(0x1.8a1f6p+0, acos(@as(f32, 0x8p-8)));
    try std.testing.expectEqual(0x1.91dfb6p+0, acos(@as(f32, 0x4p-12)));
    try std.testing.expectEqual(0x1.921db6p+0, acos(@as(f32, 0x2p-16)));
    try std.testing.expectEqual(0x1.921fa6p+0, acos(@as(f32, 0x1p-20)));
    try std.testing.expectEqual(0x1.921fb4p+0, acos(@as(f32, 0x8p-28)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x4p-32)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x2p-36)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x1p-40)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x8p-48)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x4p-52)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x2p-56)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x1p-60)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x8p-68)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x4p-72)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x2p-76)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x1p-80)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x8p-88)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x4p-92)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x2p-96)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x1p-100)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x8p-108)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x4p-112)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x2p-116)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x1p-120)));
    try std.testing.expectEqual(0x1.9a200ap+0, acos(@as(f32, -0x8p-8)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x8p-28)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x8p-48)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x8p-68)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x8p-88)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x8p-108)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x8p-128)));
    try std.testing.expectEqual(0x1.b23ddep+0, acos(@as(f32, -0x2.0089a4p-4)));
    try std.testing.expectEqual(0x5.a2499p-4, acos(@as(f32, 0xf.04aeep-4)));
    try std.testing.expectEqual(0x1.321054p+0, acos(@as(f32, 0x5.dd2588p-4)));
    try std.testing.expectEqual(0x1.321054p+0, acos(@as(f32, 0x5.dd258p-4)));
    try std.testing.expectEqual(0x1.b59bcap+0, acos(@as(f32, -0x2.35f05p-4)));
    try std.testing.expectEqual(0x1.b59bcap+0, acos(@as(f32, -0x2.35f054p-4)));
    try std.testing.expectEqual(0x6.bc5e58p-4, acos(@as(f32, 0xe.9a5c1p-4)));
    try std.testing.expectEqual(0x6.bc5e8p-4, acos(@as(f32, 0xe.9a5cp-4)));
    try std.testing.expectEqual(0x7.e544bp-4, acos(@as(f32, 0xe.17514p-4)));
    try std.testing.expectEqual(0x7.e544dp-4, acos(@as(f32, 0xe.17513p-4)));
    try std.testing.expectEqual(0x1.532618p+0, acos(@as(f32, 0x3.e57824p-4)));
    try std.testing.expectEqual(0x1.532618p+0, acos(@as(f32, 0x3.e5782p-4)));
    try std.testing.expectEqual(0x1.7149c6p+0, acos(@as(f32, 0x2.0bee8p-4)));
    try std.testing.expectEqual(0x1.afd0cap+0, acos(@as(f32, -0x1.da00d8p-4)));
    try std.testing.expectEqual(0x3.76cf6p-12, acos(@as(f32, 0xf.ffffap-4)));
    try std.testing.expectEqual(0x3.bddd44p-12, acos(@as(f32, 0xf.ffff9p-4)));
    try std.testing.expectEqual(0x7.ffb55p-8, acos(@as(f32, 0xf.fe003p-4)));
    try std.testing.expectEqual(0x7.ffd55p-8, acos(@as(f32, 0xf.fe002p-4)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb6p+0, acos(@as(f32, -0x8p-152)));

    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x0p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, acos(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x1.0c152382d7366p+0, acos(@as(f64, 0x8p-4)));
    try std.testing.expectEqual(0x2.182a4705ae6ccp+0, acos(@as(f64, -0x8p-4)));
    try std.testing.expectEqual(0xb.9051c960ecaa8p-4, acos(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x1.70ef56p-56)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x1.70ef54p-56)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x1.70ef54646d497p-56)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x1.70ef54646d496p-56)));
    try std.testing.expectEqual(0x1.821d0965ad9b7p+0, acos(@as(f64, 0x1p-4)));
    try std.testing.expectEqual(0x1.6a09e6861f3abp-12, acos(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d412p+0, acos(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3abp-12, acos(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.000000000aaabp-16, acos(@as(f64, 0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d412p+0, acos(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, acos(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243e6a8885a3p+0, acos(@as(f64, -0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3abp-12, acos(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcep-24, acos(@as(f64, 0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d412p+0, acos(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, acos(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f691e7bbcap+0, acos(@as(f64, -0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3abp-12, acos(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4p-28, acos(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d412p+0, acos(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, acos(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a3p+0, acos(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3abp-12, acos(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4p-28, acos(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d412p+0, acos(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, acos(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a3p+0, acos(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3abp-12, acos(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4p-28, acos(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d412p+0, acos(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, acos(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a3p+0, acos(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3abp-12, acos(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4p-28, acos(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d412p+0, acos(@as(f64, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a3p+0, acos(@as(f64, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a3p+0, acos(@as(f64, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.8a1f5fe55274ap+0, acos(@as(f64, 0x8p-8)));
    try std.testing.expectEqual(0x1.91dfb5439826dp+0, acos(@as(f64, 0x4p-12)));
    try std.testing.expectEqual(0x1.921db54442d03p+0, acos(@as(f64, 0x2p-16)));
    try std.testing.expectEqual(0x1.921fa54442d18p+0, acos(@as(f64, 0x1p-20)));
    try std.testing.expectEqual(0x1.921fb4c442d18p+0, acos(@as(f64, 0x8p-28)));
    try std.testing.expectEqual(0x1.921fb54042d18p+0, acos(@as(f64, 0x4p-32)));
    try std.testing.expectEqual(0x1.921fb54422d18p+0, acos(@as(f64, 0x2p-36)));
    try std.testing.expectEqual(0x1.921fb54441d18p+0, acos(@as(f64, 0x1p-40)));
    try std.testing.expectEqual(0x1.921fb54442c98p+0, acos(@as(f64, 0x8p-48)));
    try std.testing.expectEqual(0x1.921fb54442d14p+0, acos(@as(f64, 0x4p-52)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x2p-56)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x1p-60)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x8p-68)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x4p-72)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x2p-76)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x1p-80)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x8p-88)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x4p-92)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x2p-96)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x1p-100)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x8p-108)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x4p-112)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x2p-116)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x1p-120)));
    try std.testing.expectEqual(0x1.9a200aa3332e7p+0, acos(@as(f64, -0x8p-8)));
    try std.testing.expectEqual(0x1.921fb5c442d18p+0, acos(@as(f64, -0x8p-28)));
    try std.testing.expectEqual(0x1.921fb54442d98p+0, acos(@as(f64, -0x8p-48)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x8p-68)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x8p-88)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x8p-108)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x8p-128)));
    try std.testing.expectEqual(0x1.b23ddd09f0cc1p+0, acos(@as(f64, -0x2.0089a4p-4)));
    try std.testing.expectEqual(0x5.a2498fffcffd4p-4, acos(@as(f64, 0xf.04aeep-4)));
    try std.testing.expectEqual(0x1.32105458cb00ep+0, acos(@as(f64, 0x5.dd2588p-4)));
    try std.testing.expectEqual(0x1.321054e25d71bp+0, acos(@as(f64, 0x5.dd258p-4)));
    try std.testing.expectEqual(0x1.321054e1f50c8p+0, acos(@as(f64, 0x5.dd258006121b8p-4)));
    try std.testing.expectEqual(0x1.b59bc9f3d809fp+0, acos(@as(f64, -0x2.35f05p-4)));
    try std.testing.expectEqual(0x1.b59bca3476b44p+0, acos(@as(f64, -0x2.35f054p-4)));
    try std.testing.expectEqual(0x1.b59bca12945d5p+0, acos(@as(f64, -0x2.35f051e70dbc4p-4)));
    try std.testing.expectEqual(0x6.bc5e5bb8473b8p-4, acos(@as(f64, 0xe.9a5c1p-4)));
    try std.testing.expectEqual(0x6.bc5e82df35ea8p-4, acos(@as(f64, 0xe.9a5cp-4)));
    try std.testing.expectEqual(0x6.bc5e61d72acccp-4, acos(@as(f64, 0xe.9a5c0d7fabbap-4)));
    try std.testing.expectEqual(0x6.bc5e61d72acep-4, acos(@as(f64, 0xe.9a5c0d7fabb98p-4)));
    try std.testing.expectEqual(0x7.e544b07f9332cp-4, acos(@as(f64, 0xe.17514p-4)));
    try std.testing.expectEqual(0x7.e544d2469d9fp-4, acos(@as(f64, 0xe.17513p-4)));
    try std.testing.expectEqual(0x7.e544c6955c77cp-4, acos(@as(f64, 0xe.17513589de7ap-4)));
    try std.testing.expectEqual(0x7.e544c6955c78cp-4, acos(@as(f64, 0xe.17513589de798p-4)));
    try std.testing.expectEqual(0x1.532617c14a05dp+0, acos(@as(f64, 0x3.e57824p-4)));
    try std.testing.expectEqual(0x1.532618034691ep+0, acos(@as(f64, 0x3.e5782p-4)));
    try std.testing.expectEqual(0x1.532617e527e23p+0, acos(@as(f64, 0x3.e57821d368ebap-4)));
    try std.testing.expectEqual(0x1.7149c5a449b95p+0, acos(@as(f64, 0x2.0bee8p-4)));
    try std.testing.expectEqual(0x1.afd0ca8858c9fp+0, acos(@as(f64, -0x1.da00d8p-4)));
    try std.testing.expectEqual(0x3.76cf5ec671462p-12, acos(@as(f64, 0xf.ffffap-4)));
    try std.testing.expectEqual(0x3.bddd445bc8fdep-12, acos(@as(f64, 0xf.ffff9p-4)));
    // try std.testing.expectEqual(0x3.8d25c52edd924p-12, acos(@as(f64, 0xf.ffff9b1a566bp-4)));
    try std.testing.expectEqual(0x7.ffb550aec7a58p-8, acos(@as(f64, 0xf.fe003p-4)));
    try std.testing.expectEqual(0x7.ffd552eedca5cp-8, acos(@as(f64, 0xf.fe002p-4)));
    try std.testing.expectEqual(0x7.ffc7175d8527p-8, acos(@as(f64, 0xf.fe00271d507fp-4)));
    try std.testing.expectEqual(0x7.ffc7175d8627p-8, acos(@as(f64, 0xf.fe00271d507e8p-4)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18p+0, acos(@as(f64, -0x4p-1076)));

    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x0p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, acos(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x1.0c152382d7365846p+0, acos(@as(f80, 0x8p-4)));
    try std.testing.expectEqual(0x2.182a4705ae6cb08cp+0, acos(@as(f80, -0x8p-4)));
    try std.testing.expectEqual(0xb.9051c960ecaa429p-4, acos(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x1.921fb54442d182f8p+0, acos(@as(f80, 0x1.70ef56p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f8p+0, acos(@as(f80, 0x1.70ef54p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f8p+0, acos(@as(f80, 0x1.70ef54646d497p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f8p+0, acos(@as(f80, 0x1.70ef54646d496p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f8p+0, acos(@as(f80, 0x1.70ef54646d496894p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f8p+0, acos(@as(f80, 0x1.70ef54646d496892p-56)));
    try std.testing.expectEqual(0x1.821d0965ad9b6b24p+0, acos(@as(f80, 0x1p-4)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd2p-12, acos(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d411528p+0, acos(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd2p-12, acos(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.000000000aaaaaaap-16, acos(@as(f80, 0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d411528p+0, acos(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, acos(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243e6a8885a2fe28p+0, acos(@as(f80, -0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd2p-12, acos(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bce734p-24, acos(@as(f80, 0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d411528p+0, acos(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, acos(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f691e7bbca0ep+0, acos(@as(f80, -0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd2p-12, acos(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4.00000000000002a8p-28, acos(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d411528p+0, acos(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, acos(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a308d4p+0, acos(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd2p-12, acos(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4.00000000000002a8p-28, acos(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908p-32, acos(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d411528p+0, acos(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, acos(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a308d4p+0, acos(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.243f6a871b99226cp+0, acos(@as(f80, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd2p-12, acos(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4.00000000000002a8p-28, acos(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908p-32, acos(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d411528p+0, acos(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, acos(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a308d4p+0, acos(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.243f6a871b99226cp+0, acos(@as(f80, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd2p-12, acos(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4.00000000000002a8p-28, acos(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908p-32, acos(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d411528p+0, acos(@as(f80, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d4p+0, acos(@as(f80, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a308d4p+0, acos(@as(f80, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.243f6a871b99226cp+0, acos(@as(f80, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x1.8a1f5fe55274a05ap+0, acos(@as(f80, 0x8p-8)));
    try std.testing.expectEqual(0x1.91dfb5439826d4f2p+0, acos(@as(f80, 0x4p-12)));
    try std.testing.expectEqual(0x1.921db54442d02f14p+0, acos(@as(f80, 0x2p-16)));
    try std.testing.expectEqual(0x1.921fa54442d18466p+0, acos(@as(f80, 0x1p-20)));
    try std.testing.expectEqual(0x1.921fb4c442d1846ap+0, acos(@as(f80, 0x8p-28)));
    try std.testing.expectEqual(0x1.921fb54042d1846ap+0, acos(@as(f80, 0x4p-32)));
    try std.testing.expectEqual(0x1.921fb54422d1846ap+0, acos(@as(f80, 0x2p-36)));
    try std.testing.expectEqual(0x1.921fb54441d1846ap+0, acos(@as(f80, 0x1p-40)));
    try std.testing.expectEqual(0x1.921fb54442c9846ap+0, acos(@as(f80, 0x8p-48)));
    try std.testing.expectEqual(0x1.921fb54442d1446ap+0, acos(@as(f80, 0x4p-52)));
    try std.testing.expectEqual(0x1.921fb54442d1826ap+0, acos(@as(f80, 0x2p-56)));
    try std.testing.expectEqual(0x1.921fb54442d1845ap+0, acos(@as(f80, 0x1p-60)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x8p-68)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x4p-72)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x2p-76)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x1p-80)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x8p-88)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x4p-92)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x2p-96)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x1p-100)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x8p-108)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x4p-112)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x2p-116)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x1p-120)));
    try std.testing.expectEqual(0x1.9a200aa3332e6878p+0, acos(@as(f80, -0x8p-8)));
    try std.testing.expectEqual(0x1.921fb5c442d1846ap+0, acos(@as(f80, -0x8p-28)));
    try std.testing.expectEqual(0x1.921fb54442d9846ap+0, acos(@as(f80, -0x8p-48)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x8p-68)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x8p-88)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x8p-108)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x8p-128)));
    try std.testing.expectEqual(0x1.b23ddd09f0cc16cap+0, acos(@as(f80, -0x2.0089a4p-4)));
    try std.testing.expectEqual(0x5.a2498fffcffd3be8p-4, acos(@as(f80, 0xf.04aeep-4)));
    try std.testing.expectEqual(0x1.32105458cb00d85ap+0, acos(@as(f80, 0x5.dd2588p-4)));
    try std.testing.expectEqual(0x1.321054e25d71b6fcp+0, acos(@as(f80, 0x5.dd258p-4)));
    try std.testing.expectEqual(0x1.321054e1f50c7ffap+0, acos(@as(f80, 0x5.dd258006121b8p-4)));
    try std.testing.expectEqual(0x1.b59bc9f3d809e8bep+0, acos(@as(f80, -0x2.35f05p-4)));
    try std.testing.expectEqual(0x1.b59bca3476b43f62p+0, acos(@as(f80, -0x2.35f054p-4)));
    try std.testing.expectEqual(0x1.b59bca12945d4ffep+0, acos(@as(f80, -0x2.35f051e70dbc4p-4)));
    try std.testing.expectEqual(0x6.bc5e5bb8473b8b2p-4, acos(@as(f80, 0xe.9a5c1p-4)));
    try std.testing.expectEqual(0x6.bc5e82df35ea6dap-4, acos(@as(f80, 0xe.9a5cp-4)));
    try std.testing.expectEqual(0x6.bc5e61d72accaa2p-4, acos(@as(f80, 0xe.9a5c0d7fabbap-4)));
    try std.testing.expectEqual(0x6.bc5e61d72acde358p-4, acos(@as(f80, 0xe.9a5c0d7fabb98p-4)));
    try std.testing.expectEqual(0x6.bc5e61d72acd7c68p-4, acos(@as(f80, 0xe.9a5c0d7fabb9aa1p-4)));
    try std.testing.expectEqual(0x7.e544b07f9332da58p-4, acos(@as(f80, 0xe.17514p-4)));
    try std.testing.expectEqual(0x7.e544d2469d9f1438p-4, acos(@as(f80, 0xe.17513p-4)));
    try std.testing.expectEqual(0x7.e544c6955c77c6p-4, acos(@as(f80, 0xe.17513589de7ap-4)));
    try std.testing.expectEqual(0x7.e544c6955c78d438p-4, acos(@as(f80, 0xe.17513589de798p-4)));
    try std.testing.expectEqual(0x7.e544c6955c785f7p-4, acos(@as(f80, 0xe.17513589de79b75p-4)));
    try std.testing.expectEqual(0x1.532617c14a05d60ep+0, acos(@as(f80, 0x3.e57824p-4)));
    try std.testing.expectEqual(0x1.532618034691e422p+0, acos(@as(f80, 0x3.e5782p-4)));
    try std.testing.expectEqual(0x1.532617e527e23p+0, acos(@as(f80, 0x3.e57821d368ebap-4)));
    try std.testing.expectEqual(0x1.7149c5a449b958p+0, acos(@as(f80, 0x2.0bee8p-4)));
    try std.testing.expectEqual(0x1.afd0ca8858c9ea46p+0, acos(@as(f80, -0x1.da00d8p-4)));
    try std.testing.expectEqual(0x3.76cf5ec671462a94p-12, acos(@as(f80, 0xf.ffffap-4)));
    try std.testing.expectEqual(0x3.bddd445bc8fdd904p-12, acos(@as(f80, 0xf.ffff9p-4)));
    try std.testing.expectEqual(0x3.8d25c52edd923p-12, acos(@as(f80, 0xf.ffff9b1a566bp-4)));
    try std.testing.expectEqual(0x7.ffb550aec7a56398p-8, acos(@as(f80, 0xf.fe003p-4)));
    try std.testing.expectEqual(0x7.ffd552eedca5adfp-8, acos(@as(f80, 0xf.fe002p-4)));
    try std.testing.expectEqual(0x7.ffc7175d8526e8dp-8, acos(@as(f80, 0xf.fe00271d507fp-4)));
    try std.testing.expectEqual(0x7.ffc7175d8626fa98p-8, acos(@as(f80, 0xf.fe00271d507e8p-4)));
    try std.testing.expectEqual(0x7.ffc7175d855b4c78p-8, acos(@as(f80, 0xf.fe00271d507ee5dp-4)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d1846ap+0, acos(@as(f80, -0x8p-16448)));

    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x0p+0)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x0p+0)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, acos(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x1.0c152382d73658465bb32e0f567bp+0, acos(@as(f128, 0x8p-4)));
    try std.testing.expectEqual(0x2.182a4705ae6cb08cb7665c1eacf6p+0, acos(@as(f128, -0x8p-4)));
    try std.testing.expectEqual(0xb.9051c960ecaa428dd6deb6696c7p-4, acos(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a36c51701b8p+0, acos(@as(f128, 0x1.70ef56p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a38c51701b8p+0, acos(@as(f128, 0x1.70ef54p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a3860a9b848p+0, acos(@as(f128, 0x1.70ef54646d497p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a3860a9b858p+0, acos(@as(f128, 0x1.70ef54646d496p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a3860a9b85p+0, acos(@as(f128, 0x1.70ef54646d496894p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a3860a9b85p+0, acos(@as(f128, 0x1.70ef54646d496892p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a3860a9b85p+0, acos(@as(f128, 0x1.70ef54646d496892137dfd73f5aap-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a3860a9b85p+0, acos(@as(f128, 0x1.70ef54646d496892137dfd73f5a9p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a3860a9b85p+0, acos(@as(f128, 0x1.70ef54646d496892137dfd73f6p-56)));
    try std.testing.expectEqual(0x1.921fb54442d182f89a3860a9b85p+0, acos(@as(f128, 0x1.70ef54646d496892137dfd73f58p-56)));
    try std.testing.expectEqual(0x1.821d0965ad9b6b237e01535f8604p+0, acos(@as(f128, 0x1p-4)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd17681ee6db02ap-12, acos(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d4115283602220f1c96p+0, acos(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd17681ee6db02ap-12, acos(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.000000000aaaaaaaabddddddde0cp-16, acos(@as(f128, 0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d4115283602220f1c96p+0, acos(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, acos(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243e6a8885a2fe28686ede502592p+0, acos(@as(f128, -0xf.fffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd17681ee6db02ap-12, acos(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bce73430d912615775p-24, acos(@as(f128, 0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d4115283602220f1c96p+0, acos(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, acos(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f691e7bbca0df563255fd2a5ep+0, acos(@as(f128, -0xf.fffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd17681ee6db02ap-12, acos(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4.00000000000002aaaaaaaaaaaabp-28, acos(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d4115283602220f1c96p+0, acos(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, acos(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a308d313195f8358c6p+0, acos(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd17681ee6db02ap-12, acos(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4.00000000000002aaaaaaaaaaaabp-28, acos(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908d1269144e99p-32, acos(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d4115283602220f1c96p+0, acos(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, acos(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a308d313195f8358c6p+0, acos(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.243f6a871b99226b1f5cc125324ap+0, acos(@as(f128, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd17681ee6db02ap-12, acos(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4.00000000000002aaaaaaaaaaaabp-28, acos(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908d1269144e99p-32, acos(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0xb.504f333f9de6484597d89b3754e8p-56, acos(@as(f128, 0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d4115283602220f1c96p+0, acos(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, acos(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a308d313195f8358c6p+0, acos(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.243f6a871b99226b1f5cc125324ap+0, acos(@as(f128, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a2fd82c3e64a901d28p+0, acos(@as(f128, -0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(0x0p+0, acos(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.6a09e6861f3aadd17681ee6db02ap-12, acos(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x4.00000000000002aaaaaaaaaaaabp-28, acos(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x1.6a09e667f3bcc908d1269144e99p-32, acos(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x1p-56, acos(@as(f128, 0xf.fffffffffffffffffffffffffff8p-4)));
    try std.testing.expectEqual(0xb.504f333f9de6484597d89b3754e8p-56, acos(@as(f128, 0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(0x3.2428c9ea1d4115283602220f1c96p+0, acos(@as(f128, -0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a308d313198a2e037p+0, acos(@as(f128, -0x1p+0)));
    try std.testing.expectEqual(0x3.243f6a4885a308d313195f8358c6p+0, acos(@as(f128, -0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.243f6a871b99226b1f5cc125324ap+0, acos(@as(f128, -0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x3.243f6a8885a307d313198a2e037p+0, acos(@as(f128, -0xf.fffffffffffffffffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.243f6a8885a2fd82c3e64a901d28p+0, acos(@as(f128, -0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(0x1.8a1f5fe55274a05ad2c29a890fc4p+0, acos(@as(f128, 0x8p-8)));
    try std.testing.expectEqual(0x1.91dfb5439826d4f211e796c22ac3p+0, acos(@as(f128, 0x4p-12)));
    try std.testing.expectEqual(0x1.921db54442d02f143435095b45f7p+0, acos(@as(f128, 0x2p-16)));
    try std.testing.expectEqual(0x1.921fa54442d18466dee21a6c55dap+0, acos(@as(f128, 0x1p-20)));
    try std.testing.expectEqual(0x1.921fb4c442d1846989876fc1ac63p+0, acos(@as(f128, 0x8p-28)));
    try std.testing.expectEqual(0x1.921fb54042d18469898cc50c570ep+0, acos(@as(f128, 0x4p-32)));
    try std.testing.expectEqual(0x1.921fb54422d18469898cc51701a3p+0, acos(@as(f128, 0x2p-36)));
    try std.testing.expectEqual(0x1.921fb54441d18469898cc51701b8p+0, acos(@as(f128, 0x1p-40)));
    try std.testing.expectEqual(0x1.921fb54442c98469898cc51701b8p+0, acos(@as(f128, 0x8p-48)));
    try std.testing.expectEqual(0x1.921fb54442d14469898cc51701b8p+0, acos(@as(f128, 0x4p-52)));
    try std.testing.expectEqual(0x1.921fb54442d18269898cc51701b8p+0, acos(@as(f128, 0x2p-56)));
    try std.testing.expectEqual(0x1.921fb54442d18459898cc51701b8p+0, acos(@as(f128, 0x1p-60)));
    try std.testing.expectEqual(0x1.921fb54442d18469098cc51701b8p+0, acos(@as(f128, 0x8p-68)));
    try std.testing.expectEqual(0x1.921fb54442d18469858cc51701b8p+0, acos(@as(f128, 0x4p-72)));
    try std.testing.expectEqual(0x1.921fb54442d18469896cc51701b8p+0, acos(@as(f128, 0x2p-76)));
    try std.testing.expectEqual(0x1.921fb54442d18469898bc51701b8p+0, acos(@as(f128, 0x1p-80)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cbd1701b8p+0, acos(@as(f128, 0x8p-88)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc4d701b8p+0, acos(@as(f128, 0x4p-92)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51501b8p+0, acos(@as(f128, 0x2p-96)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc516f1b8p+0, acos(@as(f128, 0x1p-100)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc5170138p+0, acos(@as(f128, 0x8p-108)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b4p+0, acos(@as(f128, 0x4p-112)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x2p-116)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x1p-120)));
    try std.testing.expectEqual(0x1.9a200aa3332e68784056efa4f3adp+0, acos(@as(f128, -0x8p-8)));
    try std.testing.expectEqual(0x1.921fb5c442d1846989921a6c570ep+0, acos(@as(f128, -0x8p-28)));
    try std.testing.expectEqual(0x1.921fb54442d98469898cc51701b8p+0, acos(@as(f128, -0x8p-48)));
    try std.testing.expectEqual(0x1.921fb54442d1846a098cc51701b8p+0, acos(@as(f128, -0x8p-68)));
    try std.testing.expectEqual(0x1.921fb54442d18469898ccd1701b8p+0, acos(@as(f128, -0x8p-88)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc5170238p+0, acos(@as(f128, -0x8p-108)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x8p-128)));
    try std.testing.expectEqual(0x1.b23ddd09f0cc16c92820303d9531p+0, acos(@as(f128, -0x2.0089a4p-4)));
    try std.testing.expectEqual(0x5.a2498fffcffd3be855770b48848p-4, acos(@as(f128, 0xf.04aeep-4)));
    try std.testing.expectEqual(0x1.32105458cb00d859a030e4b5837ap+0, acos(@as(f128, 0x5.dd2588p-4)));
    try std.testing.expectEqual(0x1.321054e25d71b6fb99d7a46a0848p+0, acos(@as(f128, 0x5.dd258p-4)));
    try std.testing.expectEqual(0x1.321054e1f50c7ffa4b9e7d998261p+0, acos(@as(f128, 0x5.dd258006121b8p-4)));
    try std.testing.expectEqual(0x1.b59bc9f3d809e8bdaeb280b0409fp+0, acos(@as(f128, -0x2.35f05p-4)));
    try std.testing.expectEqual(0x1.b59bca3476b43f6179a7770ed704p+0, acos(@as(f128, -0x2.35f054p-4)));
    try std.testing.expectEqual(0x1.b59bca12945d4ffe399019e670ecp+0, acos(@as(f128, -0x2.35f051e70dbc4p-4)));
    try std.testing.expectEqual(0x6.bc5e5bb8473b8b1d4f737c4445ccp-4, acos(@as(f128, 0xe.9a5c1p-4)));
    try std.testing.expectEqual(0x6.bc5e82df35ea6da2ad38b839a7f8p-4, acos(@as(f128, 0xe.9a5cp-4)));
    try std.testing.expectEqual(0x6.bc5e61d72accaa1effa645786978p-4, acos(@as(f128, 0xe.9a5c0d7fabbap-4)));
    try std.testing.expectEqual(0x6.bc5e61d72acde35677699ecc47dcp-4, acos(@as(f128, 0xe.9a5c0d7fabb98p-4)));
    try std.testing.expectEqual(0x6.bc5e61d72acd7c691d2e8d119e4cp-4, acos(@as(f128, 0xe.9a5c0d7fabb9aa1p-4)));
    try std.testing.expectEqual(0x7.e544b07f9332da597436a3d6aa38p-4, acos(@as(f128, 0xe.17514p-4)));
    try std.testing.expectEqual(0x7.e544d2469d9f143b4750eec2959p-4, acos(@as(f128, 0xe.17513p-4)));
    try std.testing.expectEqual(0x7.e544c6955c77c5fc56134aa509f8p-4, acos(@as(f128, 0xe.17513589de7ap-4)));
    try std.testing.expectEqual(0x7.e544c6955c78d434a8d1ef5ead5p-4, acos(@as(f128, 0xe.17513589de798p-4)));
    try std.testing.expectEqual(0x7.e544c6955c785f6f92104d6ff92cp-4, acos(@as(f128, 0xe.17513589de79b75p-4)));
    try std.testing.expectEqual(0x1.532617c14a05d60df1b559c6ae81p+0, acos(@as(f128, 0x3.e57824p-4)));
    try std.testing.expectEqual(0x1.532618034691e421e199dcbd888bp+0, acos(@as(f128, 0x3.e5782p-4)));
    try std.testing.expectEqual(0x1.532617e527e22fffe0ea49e76f54p+0, acos(@as(f128, 0x3.e57821d368ebap-4)));
    try std.testing.expectEqual(0x1.7149c5a449b957ffe712405f62fbp+0, acos(@as(f128, 0x2.0bee8p-4)));
    // try std.testing.expectEqual(0x1.afd0ca8858c9ea46ebc1be7c97dfp+0, acos(@as(f128, -0x1.da00d8p-4)));
    try std.testing.expectEqual(0x3.76cf5ec671462a93def0c25c152p-12, acos(@as(f128, 0xf.ffffap-4)));
    try std.testing.expectEqual(0x3.bddd445bc8fdd903a92245263ea4p-12, acos(@as(f128, 0xf.ffff9p-4)));
    try std.testing.expectEqual(0x3.8d25c52edd92300000026414cf62p-12, acos(@as(f128, 0xf.ffff9b1a566bp-4)));
    try std.testing.expectEqual(0x7.ffb550aec7a5639438e648655f7p-8, acos(@as(f128, 0xf.fe003p-4)));
    try std.testing.expectEqual(0x7.ffd552eedca5adef3be2b7553514p-8, acos(@as(f128, 0xf.fe002p-4)));
    try std.testing.expectEqual(0x7.ffc7175d8526e8d2a5d9ca757ae8p-8, acos(@as(f128, 0xf.fe00271d507fp-4)));
    try std.testing.expectEqual(0x7.ffc7175d8626fa9a9b2f299f3b2p-8, acos(@as(f128, 0xf.fe00271d507e8p-4)));
    try std.testing.expectEqual(0x7.ffc7175d855b4c75eeab24686aecp-8, acos(@as(f128, 0xf.fe00271d507ee5dp-4)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x4p-128)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x4p-1024)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x4p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x2p-16384)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x8p-972)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x8p-152)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x4p-1076)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x8p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x4p-16448)));
    try std.testing.expectEqual(0x1.921fb54442d18469898cc51701b8p+0, acos(@as(f128, -0x4p-16496)));
}
