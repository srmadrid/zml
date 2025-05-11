const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const uasncs = @import("uasncs.zig");
const root_tbl = @import("root_tbl.zig");
const powtwo_tbl = @import("powtwo_tbl.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn asin(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return asin(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, asin32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_asinf.c
                    return asin32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_asin.c
                    return asin64(x);
                },
                f80 => return cast(f80, asin128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_asinl.c
                    return asin128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn as_special(x: f32) f32 {
    const ax: u32 = @as(u32, @bitCast(x)) << 1;
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

fn asin32(x: f32) f32 {
    const pi2: f64 = 0x1.921fb54442d18p+0;

    const ax: u32 = @as(u32, @bitCast(x)) << 1;
    const xs: f64 = cast(f64, x, .{});
    if (ax > 0x7f << 24) {
        @branchHint(.unlikely);
        return as_special(x);
    }

    if (ax < 0x7ec29000) {
        @branchHint(.likely);
        if (ax < 115 << 24) {
            return @mulAdd(f32, x, 0x1p-25, x);
        }

        const b: [16]f64 = .{
            0x1.0000000000005p+0,  0x1.55557aeca105dp-3, 0x1.3314ec3db7d12p-4,
            0x1.775738a5a6f92p-5,  0x1.5d5f7ce1c8538p-8, 0x1.605c6d58740fp-2,
            -0x1.5728b732d73c6p+1, 0x1.f152170f151ebp+3, -0x1.f962ea3ca992ep+5,
            0x1.71971e17375ap+7,   -0x1.860512b4ba23p+8, 0x1.26a3b8d4bdb14p+9,
            -0x1.36f2ea5698b51p+9, 0x1.b3d722aebfa2ep+8, -0x1.6cf89703b1289p+7,
            0x1.1518af6a65e2dp+5,
        };
        const z: f64 = xs;
        const z2: f64 = z * z;
        const z4: f64 = z2 * z2;
        const z8: f64 = z4 * z4;
        const z16: f64 = z8 * z8;
        const r: f64 = z * ((((b[0] + z2 * b[1]) + z4 * (b[2] + z2 * b[3])) + z8 * ((b[4] + z2 * b[5]) + z4 * (b[6] + z2 * b[7]))) + z16 * (((b[8] + z2 * b[9]) + z4 * (b[10] + z2 * b[11])) + z8 * ((b[12] + z2 * b[13]) + z4 * (b[14] + z2 * b[15]))));
        const ub: f32 = cast(f32, r, .{});
        const lb: f32 = cast(f32, r - z * 0x1.efa8ebp-31, .{});
        if (ub == lb)
            return ub;

        return cast(f32, r, .{});
    }

    if (ax < (0x7e << 24)) {
        const c: [12]f64 = .{
            0x1.555555555529cp-3, 0x1.333333337e0ddp-4,  0x1.6db6db3b4465ep-5,
            0x1.f1c72e13ac306p-6, 0x1.6e89cebe06bc4p-6,  0x1.1c6dcf5289094p-6,
            0x1.c6dbbcc7c6315p-7, 0x1.8f8dc2615e996p-7,  0x1.a5833b7bf15e8p-8,
            0x1.43f44ace1665cp-6, -0x1.0fb17df881c73p-6, 0x1.07520c026b2d6p-5,
        };
        const z: f64 = xs;
        const z2: f64 = z * z;
        const c0: f64 = poly12(z2, &c);
        return cast(f32, z + (z * z2) * c0, .{});
    } else {
        if (ax == 0x7e55688a) {
            @branchHint(.unlikely);
            return float.copysign(@as(f32, 0x1.75b8a2p-1), x) + float.copysign(@as(f32, 0x1p-26), x);
        }

        if (ax == 0x7e107434) {
            @branchHint(.unlikely);
            return float.copysign(@as(f32, 0x1.1f4b64p-1), x) + float.copysign(@as(f32, 0x1p-26), x);
        }

        const bx: f64 = float.abs(xs);
        const z: f64 = 1 - bx;
        const s: f64 = float.sqrt(z);
        const c: [12]f64 = .{
            0x1.6a09e667f3bcbp+0,  0x1.e2b7dddff2db9p-4,   0x1.b27247ab42dbcp-6,
            0x1.02995cc4e0744p-7,  0x1.5ffb0276ec8eap-9,   0x1.033885a928decp-10,
            0x1.911f2be23f8c7p-12, 0x1.4c3c55d2437fdp-13,  0x1.af477e1d7b461p-15,
            0x1.abd6bdff67dcbp-15, -0x1.1717e86d0fa28p-16, 0x1.6ff526de46023p-16,
        };

        const r: f64 = pi2 - s * poly12(z, &c);
        return cast(f32, float.copysign(r, xs), .{});
    }
}

// asin with max ULP of ~0.516 based on random sampling.
fn asin64(x: f64) f64 {
    const u: [2]i32 = @bitCast(x);
    const m: i32 = u[uasncs.HIGH_HALF];
    var k: i32 = 0x7fffffff & m; // no sign
    if (k < 0x3e500000) {
        if (float.abs(x) < std.math.floatMin(f64)) {
            const vx: f64 = x * x;
            std.mem.doNotOptimizeAway(vx);
        }

        return x; // for x->0 => sin(x)=x
    } else if (k < 0x3fc00000) { //----------------------2^-26 <= |x| < 2^ -3    -----------------
        const x2: f64 = x * x;
        const t: f64 = (((((uasncs.f6 * x2 + uasncs.f5) * x2 + uasncs.f4) * x2 + uasncs.f3) * x2 + uasncs.f2) * x2 + uasncs.f1) * (x2 * x);
        const res: f64 = x + t; //  res=arcsin(x) according to Taylor series
        // Max ULP is 0.513.
        return res;
    } else if (k < 0x3fe00000) { //---------------------0.125 <= |x| < 0.5 -----------------------------
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
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)]));
        t += p;
        const res: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)])) + t;
        // Max ULP is 0.524.
        return if (m > 0) res else -res;
    } else if (k < 0x3fe80000) { //-------------------- 0.5 <= |x| < 0.75 -----------------------------
        const n: i32 = 1056 + ((k & 0x000fe000) >> 11) * 3;

        var xx: f64 = undefined;
        if (m > 0) {
            xx = x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        } else {
            xx = -x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        }

        var t: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 1)])) * xx;
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])) + xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)]))))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)]));
        t += p;
        const res: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 9)])) + t;
        // Max ULP is 0.505.
        return if (m > 0) res else -res;
    } else if (k < 0x3fed8000) { //--------------------- 0.75 <= |x|< 0.921875 ----------------------
        const n: i32 = 992 + ((k & 0x000fe000) >> 13) * 13;

        var xx: f64 = undefined;
        if (m > 0) {
            xx = x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        } else {
            xx = -x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        }

        var t: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 1)])) * xx;
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)])) + xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)])))))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 9)]));
        t += p;
        const res: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 10)])) + t;
        // Max ULP is 0.505.
        return if (m > 0) res else -res;
    } else if (k < 0x3fee8000) { //-------------------0.921875 <= |x| < 0.953125 ------------------------
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
        const res: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 11)])) + t;
        // Max ULP is 0.505.
        return if (m > 0) res else -res;
    } else if (k < 0x3fef0000) { //--------------------0.953125 <= |x| < 0.96875 ------------------------
        const n: i32 = 768 + ((k & 0x000fe000) >> 13) * 15;

        var xx: f64 = undefined;
        if (m > 0) {
            xx = x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        } else {
            xx = -x - @as(f64, @bitCast(uasncs.asncs[@intCast(n)]));
        }

        var t: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 1)])) * xx;
        const p: f64 = xx * xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 2)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 3)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 4)])) +
            xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 5)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 6)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 7)])) + xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 8)])) +
                xx * (@as(f64, @bitCast(uasncs.asncs[@intCast(n + 9)])) + xx * @as(f64, @bitCast(uasncs.asncs[@intCast(n + 10)])))))))))) + @as(f64, @bitCast(uasncs.asncs[@intCast(n + 11)]));
        t += p;
        const res: f64 = @as(f64, @bitCast(uasncs.asncs[@intCast(n + 12)])) + t;
        // Max ULP is 0.505.
        return if (m > 0) res else -res;
    } else if (k < 0x3ff00000) { //--------------------0.96875 <= |x| < 1 --------------------------------
        const z: f64 = 0.5 * (if (m > 0) 1 - x else 1 + x);
        const v: [2]i32 = @bitCast(z);
        k = v[uasncs.HIGH_HALF];
        var t: f64 = root_tbl.inroot[@intCast((k & 0x001fffff) >> 14)] * powtwo_tbl.powtwo[@intCast(511 - (k >> 21))];
        const r: f64 = 1 - t * t * z;
        t = t * (uasncs.rt0 + r * (uasncs.rt1 + r * (uasncs.rt2 + r * uasncs.rt3)));
        const c: f64 = t * z;
        t = c * (1.5 - 0.5 * t * c);
        const y: f64 = (c + uasncs.t24) - uasncs.t24;
        const cc: f64 = (z - y * y) / (t + y);
        const p: f64 = (((((uasncs.f6 * z + uasncs.f5) * z + uasncs.f4) * z + uasncs.f3) * z + uasncs.f2) * z + uasncs.f1) * z;
        const cor: f64 = (@as(f64, @bitCast(uasncs.hp1)) - 2 * cc) - 2 * (y + cc) * p;
        const res1: f64 = @as(f64, @bitCast(uasncs.hp0)) - 2 * y;
        const res: f64 = res1 + cor;
        // Max ULP is 0.5015.
        return if (m > 0) res else -res;
    } else if (k == 0x3ff00000 and u[uasncs.LOW_HALF] == 0) { //---------------------------- |x|>=1 -------------------------------
        return if (m > 0) @as(f64, @bitCast(uasncs.hp0)) else -@as(f64, @bitCast(uasncs.hp0));
    } else {
        return (x - x) / (x - x);
    }
}

fn asin128(x: f128) f128 {
    const huge: f128 = 1.0e+4932;
    const pio2_hi: f128 = 1.5707963267948966192313216916397514420986;
    const pio2_lo: f128 = 4.3359050650618905123985220130216759843812e-35;
    const pio4_hi: f128 = 7.8539816339744830961566084581987569936977e-1;
    // coefficient for R(x^2)
    // asin(x): f128 = x + x^3 pS(x^2) / qS(x^2)
    // 0 <= x <= 0.5
    // peak relative error 1.9e-35
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
    // asin(0.5625 + x): f128 = asin(0.5625) + x rS(x) / sS(x)
    // -0.0625 <= x <= 0.0625
    // peak relative error 3.3e-35
    const rS0: f128 = -5.619049346208901520945464704848780243887e0;
    const rS1: f128 = 4.460504162777731472539175700169871920352e1;
    const rS2: f128 = -1.317669505315409261479577040530751477488e2;
    const rS3: f128 = 1.626532582423661989632442410808596009227e2;
    const rS4: f128 = -3.144806644195158614904369445440583873264e1;
    const rS5: f128 = -9.806674443470740708765165604769099559553e1;
    const rS6: f128 = 5.708468492052010816555762842394927806920e1;
    const rS7: f128 = 1.396540499232262112248553357962639431922e1;
    const rS8: f128 = -1.126243289311910363001762058295832610344e1;
    const rS9: f128 = -4.956179821329901954211277873774472383512e-1;
    const rS10: f128 = 3.313227657082367169241333738391762525780e-1;
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
    const asinr5625: f128 = 5.9740641664535021430381036628424864397707e-1;

    var flag: i32 = 0;
    var u: ldbl128.ieee_f128_shape32 = @bitCast(x);
    const sign: u32 = u.w0;
    const ix: i32 = @as(i32, @bitCast(sign)) & 0x7fffffff;
    u.w0 = @bitCast(ix); // |x|
    var t: f128 = undefined;
    var w: f128 = undefined;
    if (ix >= 0x3fff0000) { // |x|>= 1
        if (ix == 0x3fff0000 and (u.w1 | u.w2 | u.w3) == 0) {
            // asin(1)=+-pi/2 with inexact
            return x * pio2_hi + x * pio2_lo;
        }
        return (x - x) / (x - x); // asin(|x|>1) is NaN
    } else if (ix < 0x3ffe0000) { // |x| < 0.5
        if (ix < 0x3fc60000) { // |x| < 2**-57
            if (float.abs(x) < std.math.floatMin(f128)) {
                const vx: f128 = x * x;
                std.mem.doNotOptimizeAway(vx);
            }

            const force_inexact: f128 = huge + x;
            std.mem.doNotOptimizeAway(force_inexact);
            return x; // return x with inexact if x!=0
        } else {
            t = x * x;
            // Mark to use pS, qS later on.
            flag = 1;
        }
    } else if (ix < 0x3ffe4000) { // 0.625
        t = @as(f128, @bitCast(u)) - 0.5625;
        const p: f128 = ((((((((((rS10 * t + rS9) * t + rS8) * t + rS7) * t + rS6) * t + rS5) * t + rS4) * t + rS3) * t + rS2) * t + rS1) * t + rS0) * t;

        const q: f128 = (((((((((t + sS9) * t + sS8) * t + sS7) * t + sS6) * t + sS5) * t + sS4) * t + sS3) * t + sS2) * t + sS1) * t + sS0;
        t = asinr5625 + p / q;

        if ((sign & 0x80000000) == 0) {
            return t;
        } else {
            return -t;
        }
    } else {
        // 1 > |x| >= 0.625
        w = 1 - @as(f128, @bitCast(u));
        t = w * 0.5;
    }

    var p: f128 = (((((((((pS9 * t + pS8) * t + pS7) * t + pS6) * t + pS5) * t + pS4) * t + pS3) * t + pS2) * t + pS1) * t + pS0) * t;

    var q: f128 = ((((((((t + qS8) * t + qS7) * t + qS6) * t + qS5) * t + qS4) * t + qS3) * t + qS2) * t + qS1) * t + qS0;

    if (flag != 0) { // 2^-57 < |x| < 0.5
        w = p / q;
        return x + x * w;
    }

    const s = float.sqrt(t);
    if (ix >= 0x3ffef333) { // |x| > 0.975
        w = p / q;
        t = pio2_hi - (2 * (s + s * w) - pio2_lo);
    } else {
        u = @bitCast(s);
        u.w3 = 0;
        u.w2 = 0;
        w = @bitCast(u);
        const c: f128 = (t - w * w) / (s + w);
        const r: f128 = p / q;
        p = 2 * s * r - (pio2_lo - 2 * c);
        q = pio4_hi - 2.0 * w;
        t = pio4_hi - (p - q);
    }

    if ((sign & 0x80000000) == 0) {
        return t;
    } else {
        return -t;
    }
}
