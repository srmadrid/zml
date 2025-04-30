const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const log2_data = @import("log2_data.zig");
const ldbl128 = @import("ldbl128.zig");
const erf_data = @import("erf_data.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn log2(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return log2(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, log2_32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_log2f.c
                    return log2_32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_log2.c
                    return log2_64(x);
                },
                f80 => return cast(f80, log2_128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_log2l.c
                    return log2_128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn log2_32(x: f32) f32 {
    var ix: u32 = @bitCast(x);
    // Fix sign of zero with downward rounding when x==1.
    if (ix == 0x3f800000) {
        @branchHint(.unlikely);
        return 0;
    }

    if (ix -% 0x00800000 >= 0x7f800000 -% 0x00800000) {
        @branchHint(.unlikely);
        // x < 0x1p-126 or inf or nan.
        if (ix * 2 == 0)
            return -1 / @as(f32, 0);

        if (ix == 0x7f800000) // log2(inf) == inf.
            return x;

        if ((ix & 0x80000000) != 0 or ix * 2 >= 0xff000000)
            return (x - x) / (x - x);

        // x is subnormal, normalize it.
        ix = @bitCast(x * 0x1p23);
        ix -%= 23 << 23;
    }

    // x = 2^k z; where z is in range [0x3f330000,2*0x3f330000] and exact.
    // The range is split into 16 subintervals.
    // The ith subinterval contains z and c is near its center.
    const tmp: u32 = ix -% 0x3f330000;
    const i: u32 = (tmp >> 19) % 16;
    const top: u32 = tmp & 0xff800000;
    const iz: u32 = ix -% top;
    var k: i32 = undefined;
    {
        @setRuntimeSafety(false);
        k = @as(i32, @intCast(tmp)) >> 23; // arithmetic shift
    }
    const invc: f64 = log2_data.tab_32[i].invc;
    const logc: f64 = log2_data.tab_32[i].logc;
    const z = cast(f64, @as(f32, @bitCast(iz)), .{});

    // log2(x) = log1p(z/c-1)/ln2 + log2(c) + k
    const r: f64 = z * invc - 1;
    const y0: f64 = logc + cast(f64, k, .{});

    // Pipelined polynomial evaluation to approximate log1p(r)/ln2.
    const r2: f64 = r * r;
    var y: f64 = log2_data.poly_32[1] * r + log2_data.poly_32[2];
    y = log2_data.poly_32[0] * r2 + y;
    const p: f64 = log2_data.poly_32[3] * r + y0;
    y = y * r2 + p;
    return cast(f32, y, .{});
}

// Top 16 bits of a double.
inline fn top16(x: f64) u32 {
    return @truncate(@as(u64, @bitCast(x)) >> 48);
}

fn log2_64(x: f64) f64 {
    var ix: u64 = @bitCast(x);
    const top: u32 = top16(x);
    const LO: u64 = @bitCast(@as(f64, 1.0 - 0x1.5b51p-5));
    const HI: u64 = @bitCast(@as(f64, 1.0 + 0x1.6ab2p-5));
    if (ix -% LO < HI -% LO) {
        @branchHint(.unlikely);
        // Handle close to 1.0 inputs separately.
        // Fix sign of zero with downward rounding when x==1.
        if (ix == @as(u64, @bitCast(@as(f64, 1.0)))) {
            @branchHint(.unlikely);
            return 0;
        }

        const r: f64 = x - 1.0;
        var hi: f64 = undefined;
        var lo: f64 = undefined;
        if (true) {
            hi = r * log2_data.invln2hi_64;
            lo = r * log2_data.invln2lo_64 + @mulAdd(f64, r, log2_data.invln2hi_64, -hi);
        } else {
            const rhi: f64 = @bitCast(@as(u64, @bitCast(r)) & -1 << 32);
            const rlo: f64 = r - rhi;
            hi = rhi * log2_data.invln2hi_64;
            lo = rlo * log2_data.invln2hi_64 + r * log2_data.invln2lo_64;
        }
        const r2: f64 = r * r; // rounding error: 0x1p-62.
        const r4: f64 = r2 * r2;
        // Worst-case error is less than 0.54 ULP (0.55 ULP without fma).
        const p: f64 = r2 * (log2_data.poly1_64[0] + r * log2_data.poly1_64[1]);
        var y: f64 = hi + p;
        lo += hi - y + p;
        lo += r4 * (log2_data.poly1_64[2] + r * log2_data.poly1_64[3] + r2 * (log2_data.poly1_64[4] + r * log2_data.poly1_64[5]) + r4 * (log2_data.poly1_64[6] + r * log2_data.poly1_64[7] + r2 * (log2_data.poly1_64[8] + r * log2_data.poly1_64[9])));
        y += lo;
        return y;
    }

    if (top -% 0x0010 >= 0x7ff0 -% 0x0010) {
        @branchHint(.unlikely);
        // x < 0x1p-1022 or inf or nan.
        if (ix * 2 == 0)
            return -1 / @as(f64, 0);

        if (ix == @as(u64, @bitCast(std.math.inf(f64)))) // log(inf) == inf.
            return x;

        if ((top & 0x8000) != 0 or (top & 0x7ff0) == 0x7ff0)
            return (x - x) / (x - x);

        // x is subnormal, normalize it.
        ix = @bitCast(x * 0x1p52);
        ix -%= 52 << 52;
    }

    // x = 2^k z; where z is in range [OFF,2*OFF) and exact.
    // The range is split into N subintervals.
    // The ith subinterval contains z and c is near its center.
    const tmp: u64 = ix -% 0x3fe6000000000000;
    const i: u64 = (tmp >> 46) % 64;
    var k: i64 = undefined;
    {
        @setRuntimeSafety(false);
        k = @as(i64, @intCast(tmp)) >> 52; // arithmetic shift
    }
    const iz: u64 = ix -% (tmp & 0xfff << 52);
    const invc: f64 = log2_data.tab_64[i].invc;
    const logc: f64 = log2_data.tab_64[i].logc;
    const z: f64 = @bitCast(iz);
    const kd: f64 = cast(f64, k, .{});

    // log2(x) = log2(z/c) + log2(c) + k.
    // r ~= z/c - 1, |r| < 1/(2*N).
    var r: f64 = undefined;
    var t1: f64 = undefined;
    var t2: f64 = undefined;
    if (true) {
        // rounding error: 0x1p-55/N.
        r = @mulAdd(f64, z, invc, -1.0);
        t1 = r * log2_data.invln2hi_64;
        t2 = r * log2_data.invln2lo_64 + @mulAdd(f64, r, log2_data.invln2hi_64, -t1);
    } else {
        // rounding error: 0x1p-55/N + 0x1p-65.
        r = (z - log2_data.tab2_64[i].chi - log2_data.tab2_64[i].clo) * invc;
        const rhi: f64 = @bitCast(@as(u64, @bitCast(r)) & -1 << 32);
        const rlo: f64 = r - rhi;
        t1 = rhi * log2_data.invln2hi_64;
        t2 = rlo * log2_data.invln2hi_64 + r * log2_data.invln2lo_64;
    }

    // hi + lo = r/ln2 + log2(c) + k.
    const t3: f64 = kd + logc;
    const hi: f64 = t3 + t1;
    const lo: f64 = t3 - hi + t1 + t2;

    // log2(r+1) = r/ln2 + r^2*poly(r).
    // Evaluation is optimized assuming superscalar pipelined execution.
    const r2: f64 = r * r; // rounding error: 0x1p-54/N^2.
    const r4: f64 = r2 * r2;
    // Worst-case error if |y| > 0x1p-4: 0.547 ULP (0.550 ULP without fma).
    // ~ 0.5 + 2/N/ln2 + abs-poly-error*0x1p56 ULP (+ 0.003 ULP without fma).
    const p: f64 = log2_data.poly_64[0] + r * log2_data.poly_64[1] + r2 * (log2_data.poly_64[2] + r * log2_data.poly_64[3]) + r4 * (log2_data.poly_64[4] + r * log2_data.poly_64[5]);
    const y: f64 = lo + r2 * p + hi;
    return y;
}

fn log2_128(x: f128) f128 {
    // Coefficients for ln(1+x) = x - x**2/2 + x**3 P(x)/Q(x)
    // 1/sqrt(2) <= x < sqrt(2)
    // Theoretical peak relative error = 5.3e-37,
    // relative peak error spread = 2.3e-14
    const P: [13]f128 = .{
        1.313572404063446165910279910527789794488e4,
        7.771154681358524243729929227226708890930e4,
        2.014652742082537582487669938141683759923e5,
        3.007007295140399532324943111654767187848e5,
        2.854829159639697837788887080758954924001e5,
        1.797628303815655343403735250238293741397e5,
        7.594356839258970405033155585486712125861e4,
        2.128857716871515081352991964243375186031e4,
        3.824952356185897735160588078446136783779e3,
        4.114517881637811823002128927449878962058e2,
        2.321125933898420063925789532045674660756e1,
        4.998469661968096229986658302195402690910e-1,
        1.538612243596254322971797716843006400388e-6,
    };
    const Q: [12]f128 = .{
        3.940717212190338497730839731583397586124e4,
        2.626900195321832660448791748036714883242e5,
        7.777690340007566932935753241556479363645e5,
        1.347518538384329112529391120390701166528e6,
        1.514882452993549494932585972882995548426e6,
        1.158019977462989115839826904108208787040e6,
        6.132189329546557743179177159925690841200e5,
        2.248234257620569139969141618556349415120e5,
        5.605842085972455027590989944010492125825e4,
        9.147150349299596453976674231612674085381e3,
        9.104928120962988414618126155557301584078e2,
        4.839208193348159620282142911143429644326e1,
    };

    // Coefficients for log(x) = z + z^3 P(z^2)/Q(z^2),
    // where z = 2(x-1)/(x+1)
    // 1/sqrt(2) <= x < sqrt(2)
    // Theoretical peak relative error = 1.1e-35,
    // relative peak error spread 1.1e-9
    const R: [6]f128 = .{
        1.418134209872192732479751274970992665513e5,
        -8.977257995689735303686582344659576526998e4,
        2.048819892795278657810231591630928516206e4,
        -2.024301798136027039250415126250455056397e3,
        8.057002716646055371965756206836056074715e1,
        -8.828896441624934385266096344596648080902e-1,
    };
    const S: [6]f128 = .{
        1.701761051846631278975701529965589676574e6,
        -1.332535117259762928288745111081235577029e6,
        4.001557694070773974936904547424676279307e5,
        -5.748542087379434595104154610899551484314e4,
        3.998526750980007367835804959888064681098e3,
        -1.186359407982897997337150403816839480438e2,
    };
    // log2(e) - 1
    const LOG2EA: f128 = 4.4269504088896340735992468100189213742664595e-1;
    // sqrt(2)/2
    const SQRTH: f128 = 7.071067811865475244008443621048490392848359e-1;

    // Test for domain
    var hx: i64 = undefined;
    var lx: i64 = undefined;
    ldbl128.getWords(&hx, &lx, x);
    if (((hx & 0x7fffffffffffffff) | lx) == 0)
        return (-1 / math.abs(x)); // log2l(+-0)=-inf

    if (hx < 0)
        return (x - x) / (x - x);

    if (hx >= 0x7fff000000000000)
        return (x + x);

    if (x == 1)
        return 0;

    // separate mantissa from exponent
    // Note, frexp is used so that denormal numbers
    // will be handled properly.
    var e: i32 = undefined;
    var xx: f128 = math.frexp(x, &e);

    // logarithm using log(x) = z + z**3 P(z)/Q(z),
    // where z = 2(x-1)/x+1)
    var done = false;
    var y: f128 = undefined;
    var z: f128 = undefined;
    if ((e > 2) or (e < -2)) {
        if (xx < SQRTH) { // 2( 2x-1 )/( 2x+1 )
            e -= 1;
            z = xx - 0.5;
            y = 0.5 * z + 0.5;
        } else { //  2 (x-1)/(x+1)
            z = xx - 0.5;
            z -= 0.5;
            y = 0.5 * xx + 0.5;
        }
        xx = z / y;
        z = xx * xx;
        y = xx * (z * erf_data.neval(z, &R, 5) / erf_data.deval(z, &S, 5));
        done = true;
    }

    // logarithm using log(1+x) = x - .5x**2 + x**3 P(x)/Q(x)
    if (!done) {
        if (xx < SQRTH) {
            e -= 1;
            xx = 2.0 * xx - 1; //  2x - 1
        } else {
            xx = xx - 1;
        }
        z = xx * xx;
        y = xx * (z * erf_data.neval(xx, &P, 12) / erf_data.deval(xx, &Q, 11));
        y = y - 0.5 * z;
    }

    // Multiply log of fraction by log2(e)
    // and base 2 exponent by 1
    z = y * LOG2EA;
    z += xx * LOG2EA;
    z += y;
    z += xx;
    z += cast(f128, e, .{});
    return z;
}

test log2 {
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x1.715478p+0, log2(@as(f32, 0x2.b7e154p+0)));
    try std.testing.expectEqual(0x1.715476p+0, log2(@as(f32, 0x2.b7e15p+0)));
    try std.testing.expectEqual(0x1p+0, log2(@as(f32, 0x2p+0)));
    try std.testing.expectEqual(0x4p+0, log2(@as(f32, 0x1p+4)));
    try std.testing.expectEqual(0x8p+0, log2(@as(f32, 0x1p+8)));
    try std.testing.expectEqual(-0x6.a3fe6p-4, log2(@as(f32, 0xcp-4)));
    try std.testing.expectEqual(0x2.e2a8e8p-24, log2(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x2.e2a8e8p-24, log2(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2.e2a8e8p-24, log2(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2.e2a8e8p-24, log2(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    try std.testing.expectEqual(0x2.e2a8e8p-24, log2(@as(f32, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    // try std.testing.expectEqual(-0x1.715478p-24, log2(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    // try std.testing.expectEqual(-0x1.715478p-24, log2(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    // try std.testing.expectEqual(-0x1.715478p-24, log2(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    // try std.testing.expectEqual(-0x1.715478p-24, log2(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f32, 0x1p+0)));
    // try std.testing.expectEqual(-0x1.715478p-24, log2(@as(f32, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x3.6a6848p-4, log2(@as(f32, 0x1.28d3b4p+0)));
    try std.testing.expectEqual(-0x1.b8bdeep-4, log2(@as(f32, 0xe.d99dap-4)));
    try std.testing.expectEqual(0x7.999f28p-4, log2(@as(f32, 0x1.63d204p+0)));
    try std.testing.expectEqual(0x7.999f08p-4, log2(@as(f32, 0x1.63d202p+0)));
    try std.testing.expectEqual(-0x3.75cc4p-8, log2(@as(f32, 0xf.d9ce1p-4)));
    try std.testing.expectEqual(-0x3.75cdb4p-8, log2(@as(f32, 0xf.d9cep-4)));
    try std.testing.expectEqual(0xa.59594p-8, log2(@as(f32, 0x1.07465cp+0)));
    try std.testing.expectEqual(0xa.59567p-8, log2(@as(f32, 0x1.07465ap+0)));
    try std.testing.expectEqual(-0x2.c10694p+4, log2(@as(f32, 0xf.4dfb4p-48)));
    try std.testing.expectEqual(0xe.a1dd4p-8, log2(@as(f32, 0x1.0a588ep+0)));
    try std.testing.expectEqual(-0x6.d35688p-4, log2(@as(f32, 0xb.e77c6p-4)));
    try std.testing.expectEqual(0x6.44f93p-4, log2(@as(f32, 0x1.4fe37ep+0)));
    try std.testing.expectEqual(0x9.d9a8cp+0, log2(@as(f32, 0x3.9b0754p+8)));
    try std.testing.expectEqual(-0x6.df8b3p-4, log2(@as(f32, 0xb.e132ap-4)));
    try std.testing.expectEqual(-0x7.e84208p-4, log2(@as(f32, 0xb.5bf83p-4)));
    try std.testing.expectEqual(-0x7.e84228p-4, log2(@as(f32, 0xb.5bf82p-4)));
    try std.testing.expectEqual(-0x7.b18b68p-4, log2(@as(f32, 0xb.7704ep-4)));
    try std.testing.expectEqual(-0x7.b18b88p-4, log2(@as(f32, 0xb.7704dp-4)));
    try std.testing.expectEqual(-0x7.f27148p-4, log2(@as(f32, 0xb.56f64p-4)));
    try std.testing.expectEqual(-0x7.f27168p-4, log2(@as(f32, 0xb.56f63p-4)));
    try std.testing.expectEqual(-0x7.f84a98p-4, log2(@as(f32, 0xb.54171p-4)));
    try std.testing.expectEqual(-0x7.f84ab8p-4, log2(@as(f32, 0xb.5417p-4)));
    try std.testing.expectEqual(-0x7.ep+4, log2(@as(f32, 0x4p-128)));
    try std.testing.expectEqual(-0x9.5p+4, log2(@as(f32, 0x8p-152)));
    try std.testing.expectEqual(0x8p+4, log2(@as(f32, 0xf.fffffp+124)));

    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.715477b6069b3p+0, log2(@as(f64, 0x2.b7e154p+0)));
    try std.testing.expectEqual(0x1.715475968cddcp+0, log2(@as(f64, 0x2.b7e15p+0)));
    try std.testing.expectEqual(0x1.71547652b82ffp+0, log2(@as(f64, 0x2.b7e151628aed4p+0)));
    try std.testing.expectEqual(0x1.71547652b82fep+0, log2(@as(f64, 0x2.b7e151628aed2p+0)));
    try std.testing.expectEqual(0x1p+0, log2(@as(f64, 0x2p+0)));
    try std.testing.expectEqual(0x4p+0, log2(@as(f64, 0x1p+4)));
    try std.testing.expectEqual(0x8p+0, log2(@as(f64, 0x1p+8)));
    try std.testing.expectEqual(-0x6.a3fe5c6042978p-4, log2(@as(f64, 0xcp-4)));
    try std.testing.expectEqual(0x2.e2a8e9c2c777p-24, log2(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c777p-24, log2(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fdp-52, log2(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c777p-24, log2(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fdp-52, log2(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c777p-24, log2(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fdp-52, log2(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c777p-24, log2(@as(f64, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fdp-52, log2(@as(f64, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b8p-24, log2(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b8p-24, log2(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17fp-56, log2(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b8p-24, log2(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17fp-56, log2(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b8p-24, log2(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17fp-56, log2(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f64, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b8p-24, log2(@as(f64, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17fp-56, log2(@as(f64, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x3.6a6846287159cp-4, log2(@as(f64, 0x1.28d3b4p+0)));
    try std.testing.expectEqual(-0x1.b8bdee5bd4282p-4, log2(@as(f64, 0xe.d99dap-4)));
    try std.testing.expectEqual(0x7.999f2a94857b4p-4, log2(@as(f64, 0x1.63d204p+0)));
    try std.testing.expectEqual(0x7.999f095d7e15p-4, log2(@as(f64, 0x1.63d202p+0)));
    try std.testing.expectEqual(0x7.999f16e03b55p-4, log2(@as(f64, 0x1.63d202d04392cp+0)));
    try std.testing.expectEqual(-0x3.75cc3f2233d66p-8, log2(@as(f64, 0xf.d9ce1p-4)));
    try std.testing.expectEqual(-0x3.75cdb3f0a0b66p-8, log2(@as(f64, 0xf.d9cep-4)));
    try std.testing.expectEqual(-0x3.75ccb13d89cc6p-8, log2(@as(f64, 0xf.d9ce0b1a50e08p-4)));
    try std.testing.expectEqual(0xa.5959448ade0cp-8, log2(@as(f64, 0x1.07465cp+0)));
    try std.testing.expectEqual(0xa.5956764b14aap-8, log2(@as(f64, 0x1.07465ap+0)));
    try std.testing.expectEqual(0xa.595912bb7fffp-8, log2(@as(f64, 0x1.07465bdc7e41cp+0)));
    try std.testing.expectEqual(0xa.595912bb7fe88p-8, log2(@as(f64, 0x1.07465bdc7e41bp+0)));
    try std.testing.expectEqual(-0x2.c106931f2bfdp+4, log2(@as(f64, 0xf.4dfb4p-48)));
    try std.testing.expectEqual(0xe.a1dd43a221dp-8, log2(@as(f64, 0x1.0a588ep+0)));
    try std.testing.expectEqual(-0x6.d35688edc44a4p-4, log2(@as(f64, 0xb.e77c6p-4)));
    try std.testing.expectEqual(0x6.44f92e0fda7dp-4, log2(@as(f64, 0x1.4fe37ep+0)));
    try std.testing.expectEqual(0x9.d9a8c6de34328p+0, log2(@as(f64, 0x3.9b0754p+8)));
    try std.testing.expectEqual(-0x6.df8b2c2c5ea4p-4, log2(@as(f64, 0xb.e132ap-4)));
    try std.testing.expectEqual(-0x7.e842050c531dp-4, log2(@as(f64, 0xb.5bf83p-4)));
    try std.testing.expectEqual(-0x7.e842258fcc5d8p-4, log2(@as(f64, 0xb.5bf82p-4)));
    try std.testing.expectEqual(-0x7.e8420994680dcp-4, log2(@as(f64, 0xb.5bf82dc51f028p-4)));
    try std.testing.expectEqual(-0x7.e8420994680ecp-4, log2(@as(f64, 0xb.5bf82dc51f02p-4)));
    try std.testing.expectEqual(-0x7.b18b6b68ffa24p-4, log2(@as(f64, 0xb.7704ep-4)));
    try std.testing.expectEqual(-0x7.b18b8b9fc309cp-4, log2(@as(f64, 0xb.7704dp-4)));
    try std.testing.expectEqual(-0x7.b18b723cc4c5p-4, log2(@as(f64, 0xb.7704dc9beb05p-4)));
    try std.testing.expectEqual(-0x7.f27149af9dc8cp-4, log2(@as(f64, 0xb.56f64p-4)));
    try std.testing.expectEqual(-0x7.f2716a4172a7p-4, log2(@as(f64, 0xb.56f63p-4)));
    try std.testing.expectEqual(-0x7.f27151a15d70cp-4, log2(@as(f64, 0xb.56f63c18e93fp-4)));
    try std.testing.expectEqual(-0x7.f27151a15d71cp-4, log2(@as(f64, 0xb.56f63c18e93e8p-4)));
    try std.testing.expectEqual(-0x7.f84a9424a2fbcp-4, log2(@as(f64, 0xb.54171p-4)));
    try std.testing.expectEqual(-0x7.f84ab4beb988cp-4, log2(@as(f64, 0xb.5417p-4)));
    try std.testing.expectEqual(-0x7.f84a998412a44p-4, log2(@as(f64, 0xb.54170d5cfa9p-4)));
    try std.testing.expectEqual(-0x7.f84a998412a54p-4, log2(@as(f64, 0xb.54170d5cfa8f8p-4)));
    try std.testing.expectEqual(-0x7.ep+4, log2(@as(f64, 0x4p-128)));
    try std.testing.expectEqual(-0x3.fep+8, log2(@as(f64, 0x4p-1024)));
    try std.testing.expectEqual(-0x3.c9p+8, log2(@as(f64, 0x8p-972)));
    try std.testing.expectEqual(-0x9.5p+4, log2(@as(f64, 0x8p-152)));
    try std.testing.expectEqual(-0x4.32p+8, log2(@as(f64, 0x4p-1076)));
    try std.testing.expectEqual(0x7.ffffffe8eab88p+4, log2(@as(f64, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x4p+8, log2(@as(f64, 0xf.ffffffffffff8p+1020)));

    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.715477b6069b2d4cp+0, log2(@as(f80, 0x2.b7e154p+0)));
    try std.testing.expectEqual(0x1.715475968cddc4e6p+0, log2(@as(f80, 0x2.b7e15p+0)));
    try std.testing.expectEqual(0x1.71547652b82fececp+0, log2(@as(f80, 0x2.b7e151628aed4p+0)));
    try std.testing.expectEqual(0x1.71547652b82fdbfp+0, log2(@as(f80, 0x2.b7e151628aed2p+0)));
    try std.testing.expectEqual(0x1.71547652b82fe178p+0, log2(@as(f80, 0x2.b7e151628aed2a6cp+0)));
    try std.testing.expectEqual(0x1.71547652b82fe176p+0, log2(@as(f80, 0x2.b7e151628aed2a68p+0)));
    try std.testing.expectEqual(0x1p+0, log2(@as(f80, 0x2p+0)));
    try std.testing.expectEqual(0x4p+0, log2(@as(f80, 0x1p+4)));
    try std.testing.expectEqual(0x8p+0, log2(@as(f80, 0x1p+8)));
    try std.testing.expectEqual(-0x6.a3fe5c6042978608p-4, log2(@as(f80, 0xcp-4)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f66p-24, log2(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f66p-24, log2(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fd5ecp-52, log2(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f66p-24, log2(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fd5ecp-52, log2(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.e2a8eca5705fc2ecp-64, log2(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f66p-24, log2(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fd5ecp-52, log2(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.e2a8eca5705fc2ecp-64, log2(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f66p-24, log2(@as(f80, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fd5ecp-52, log2(@as(f80, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.e2a8eca5705fc2ecp-64, log2(@as(f80, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85fp-24, log2(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85fp-24, log2(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f39ep-56, log2(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85fp-24, log2(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f39ep-56, log2(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.71547652b82fe178p-64, log2(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85fp-24, log2(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f39ep-56, log2(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.71547652b82fe178p-64, log2(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f80, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85fp-24, log2(@as(f80, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f39ep-56, log2(@as(f80, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.71547652b82fe178p-64, log2(@as(f80, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x3.6a6846287159cebp-4, log2(@as(f80, 0x1.28d3b4p+0)));
    try std.testing.expectEqual(-0x1.b8bdee5bd428254ep-4, log2(@as(f80, 0xe.d99dap-4)));
    try std.testing.expectEqual(0x7.999f2a94857b22ep-4, log2(@as(f80, 0x1.63d204p+0)));
    try std.testing.expectEqual(0x7.999f095d7e150198p-4, log2(@as(f80, 0x1.63d202p+0)));
    try std.testing.expectEqual(0x7.999f16e03b54e0f8p-4, log2(@as(f80, 0x1.63d202d04392cp+0)));
    try std.testing.expectEqual(-0x3.75cc3f2233d65c7p-8, log2(@as(f80, 0xf.d9ce1p-4)));
    try std.testing.expectEqual(-0x3.75cdb3f0a0b65598p-8, log2(@as(f80, 0xf.d9cep-4)));
    try std.testing.expectEqual(-0x3.75ccb13d89cc66fp-8, log2(@as(f80, 0xf.d9ce0b1a50e08p-4)));
    try std.testing.expectEqual(0xa.5959448ade0bdcbp-8, log2(@as(f80, 0x1.07465cp+0)));
    try std.testing.expectEqual(0xa.5956764b14a9d7ep-8, log2(@as(f80, 0x1.07465ap+0)));
    try std.testing.expectEqual(0xa.595912bb7fff0bap-8, log2(@as(f80, 0x1.07465bdc7e41cp+0)));
    try std.testing.expectEqual(0xa.595912bb7fe899ap-8, log2(@as(f80, 0x1.07465bdc7e41bp+0)));
    try std.testing.expectEqual(0xa.595912bb7fefdddp-8, log2(@as(f80, 0x1.07465bdc7e41b52ep+0)));
    try std.testing.expectEqual(-0x2.c106931f2bfd0af4p+4, log2(@as(f80, 0xf.4dfb4p-48)));
    try std.testing.expectEqual(0xe.a1dd43a221d02a3p-8, log2(@as(f80, 0x1.0a588ep+0)));
    try std.testing.expectEqual(-0x6.d35688edc44a496p-4, log2(@as(f80, 0xb.e77c6p-4)));
    try std.testing.expectEqual(0x6.44f92e0fda7d1b48p-4, log2(@as(f80, 0x1.4fe37ep+0)));
    try std.testing.expectEqual(0x9.d9a8c6de3432cp+0, log2(@as(f80, 0x3.9b0754p+8)));
    try std.testing.expectEqual(-0x6.df8b2c2c5ea4p-4, log2(@as(f80, 0xb.e132ap-4)));
    try std.testing.expectEqual(-0x7.e842050c531d023p-4, log2(@as(f80, 0xb.5bf83p-4)));
    try std.testing.expectEqual(-0x7.e842258fcc5d9f3p-4, log2(@as(f80, 0xb.5bf82p-4)));
    try std.testing.expectEqual(-0x7.e8420994680da578p-4, log2(@as(f80, 0xb.5bf82dc51f028p-4)));
    try std.testing.expectEqual(-0x7.e8420994680ea99p-4, log2(@as(f80, 0xb.5bf82dc51f02p-4)));
    try std.testing.expectEqual(-0x7.e8420994680ea2d8p-4, log2(@as(f80, 0xb.5bf82dc51f02035p-4)));
    try std.testing.expectEqual(-0x7.b18b6b68ffa23508p-4, log2(@as(f80, 0xb.7704ep-4)));
    try std.testing.expectEqual(-0x7.b18b8b9fc309de5p-4, log2(@as(f80, 0xb.7704dp-4)));
    try std.testing.expectEqual(-0x7.b18b723cc4c4fae8p-4, log2(@as(f80, 0xb.7704dc9beb05p-4)));
    try std.testing.expectEqual(-0x7.f27149af9dc8b0fp-4, log2(@as(f80, 0xb.56f64p-4)));
    try std.testing.expectEqual(-0x7.f2716a4172a70438p-4, log2(@as(f80, 0xb.56f63p-4)));
    try std.testing.expectEqual(-0x7.f27151a15d70d158p-4, log2(@as(f80, 0xb.56f63c18e93fp-4)));
    try std.testing.expectEqual(-0x7.f27151a15d71d5e8p-4, log2(@as(f80, 0xb.56f63c18e93e8p-4)));
    try std.testing.expectEqual(-0x7.f27151a15d70f868p-4, log2(@as(f80, 0xb.56f63c18e93eecdp-4)));
    try std.testing.expectEqual(-0x7.f84a9424a2fba58p-4, log2(@as(f80, 0xb.54171p-4)));
    try std.testing.expectEqual(-0x7.f84ab4beb988b358p-4, log2(@as(f80, 0xb.5417p-4)));
    try std.testing.expectEqual(-0x7.f84a998412a436a8p-4, log2(@as(f80, 0xb.54170d5cfa9p-4)));
    try std.testing.expectEqual(-0x7.f84a998412a53b78p-4, log2(@as(f80, 0xb.54170d5cfa8f8p-4)));
    try std.testing.expectEqual(-0x7.f84a998412a489dp-4, log2(@as(f80, 0xb.54170d5cfa8fd73p-4)));
    try std.testing.expectEqual(-0x7.f84a998412a489fp-4, log2(@as(f80, 0xb.54170d5cfa8fd72p-4)));
    try std.testing.expectEqual(-0x7.ep+4, log2(@as(f80, 0x4p-128)));
    try std.testing.expectEqual(-0x3.fep+8, log2(@as(f80, 0x4p-1024)));
    try std.testing.expectEqual(-0x3.ffep+12, log2(@as(f80, 0x4p-16384)));
    try std.testing.expectEqual(-0x3.fffp+12, log2(@as(f80, 0x2p-16384)));
    try std.testing.expectEqual(-0x3.c9p+8, log2(@as(f80, 0x8p-972)));
    try std.testing.expectEqual(-0x9.5p+4, log2(@as(f80, 0x8p-152)));
    try std.testing.expectEqual(-0x4.32p+8, log2(@as(f80, 0x4p-1076)));
    try std.testing.expectEqual(-0x4.03dp+12, log2(@as(f80, 0x8p-16448)));
    try std.testing.expectEqual(0x7.ffffffe8eab88f48p+4, log2(@as(f80, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.fffffffffffffff4p+8, log2(@as(f80, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x4p+12, log2(@as(f80, 0xf.fffffffffffffffp+16380)));

    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.715477b6069b2d4b70c9ce6a329p+0, log2(@as(f128, 0x2.b7e154p+0)));
    try std.testing.expectEqual(0x1.715475968cddc4e6c2b1cbd1a7c8p+0, log2(@as(f128, 0x2.b7e15p+0)));
    try std.testing.expectEqual(0x1.71547652b82fecebf2ef1ac7b6dep+0, log2(@as(f128, 0x2.b7e151628aed4p+0)));
    try std.testing.expectEqual(0x1.71547652b82fdbf024ffffda5e62p+0, log2(@as(f128, 0x2.b7e151628aed2p+0)));
    // try std.testing.expectEqual(0x1.71547652b82fe1782731bf3f6b29p+0, log2(@as(f128, 0x2.b7e151628aed2a6cp+0)));
    try std.testing.expectEqual(0x1.71547652b82fe17607b8015c0d7ep+0, log2(@as(f128, 0x2.b7e151628aed2a68p+0)));
    try std.testing.expectEqual(0x1.71547652b82fe1777d0ffda0d23bp+0, log2(@as(f128, 0x2.b7e151628aed2a6abf7158809cf6p+0)));
    try std.testing.expectEqual(0x1.71547652b82fe1777d0ffda0d23ap+0, log2(@as(f128, 0x2.b7e151628aed2a6abf7158809cf4p+0)));
    try std.testing.expectEqual(0x1.71547652b82fe1777d0ffda0d24p+0, log2(@as(f128, 0x2.b7e151628aed2a6abf7158809dp+0)));
    try std.testing.expectEqual(0x1.71547652b82fe1777d0ffda0d1b8p+0, log2(@as(f128, 0x2.b7e151628aed2a6abf7158809cp+0)));
    try std.testing.expectEqual(0x1p+0, log2(@as(f128, 0x2p+0)));
    try std.testing.expectEqual(0x4p+0, log2(@as(f128, 0x1p+4)));
    try std.testing.expectEqual(0x8p+0, log2(@as(f128, 0x1p+8)));
    // try std.testing.expectEqual(-0x6.a3fe5c6042978605ff4edf5f9744p-4, log2(@as(f128, 0xcp-4)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f65fd01efaf723e4p-24, log2(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f65fd01efaf723e4p-24, log2(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fd5ecd95d67df53aap-52, log2(@as(f128, 0x1.0000000000001p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f65fd01efaf723e4p-24, log2(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fd5ecd95d67df53aap-52, log2(@as(f128, 0x1.0000000000001p+0)));
    // try std.testing.expectEqual(0x2.e2a8eca5705fc2ec17770e9c3416p-64, log2(@as(f128, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f65fd01efaf723e4p-24, log2(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fd5ecd95d67df53aap-52, log2(@as(f128, 0x1.0000000000001p+0)));
    // try std.testing.expectEqual(0x2.e2a8eca5705fc2ec17770e9c3416p-64, log2(@as(f128, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0xb.8aa3b295c17f0bbbe87fed068efp-108, log2(@as(f128, 0x1.000000000000000000000000008p+0)));
    try std.testing.expectEqual(0x2.e2a8e9c2c776f65fd01efaf723e4p-24, log2(@as(f128, 0x1.000002p+0)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(0x1.71547652b82fd5ecd95d67df53aap-52, log2(@as(f128, 0x1.0000000000001p+0)));
    // try std.testing.expectEqual(0x2.e2a8eca5705fc2ec17770e9c3416p-64, log2(@as(f128, 0x1.0000000000000002p+0)));
    try std.testing.expectEqual(0x1.71547652b82fe1777d0ffda0d23ap-112, log2(@as(f128, 0x1.0000000000000000000000000001p+0)));
    try std.testing.expectEqual(0xb.8aa3b295c17f0bbbe87fed068efp-108, log2(@as(f128, 0x1.000000000000000000000000008p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85efbccdf68d2e98p-24, log2(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85efbccdf68d2e98p-24, log2(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f39e6774a440c8ef8p-56, log2(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85efbccdf68d2e98p-24, log2(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f39e6774a440c8ef8p-56, log2(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.71547652b82fe17835ba38ca2e52p-64, log2(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85efbccdf68d2e98p-24, log2(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f39e6774a440c8ef8p-56, log2(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.71547652b82fe17835ba38ca2e52p-64, log2(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0x5.c551d94ae0bf85ddf43ff68349a4p-108, log2(@as(f128, 0xf.fffffffffffffffffffffffffcp-4)));
    try std.testing.expectEqual(0x0p+0, log2(@as(f128, 0x1p+0)));
    try std.testing.expectEqual(-0x1.7154770b626b85efbccdf68d2e98p-24, log2(@as(f128, 0xf.fffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f39e6774a440c8ef8p-56, log2(@as(f128, 0xf.ffffffffffff8p-4)));
    try std.testing.expectEqual(-0x1.71547652b82fe17835ba38ca2e52p-64, log2(@as(f128, 0xf.fffffffffffffffp-4)));
    try std.testing.expectEqual(-0xb.8aa3b295c17f0bbbe87fed0691d8p-116, log2(@as(f128, 0xf.fffffffffffffffffffffffffff8p-4)));
    try std.testing.expectEqual(-0x5.c551d94ae0bf85ddf43ff68349a4p-108, log2(@as(f128, 0xf.fffffffffffffffffffffffffcp-4)));
    // try std.testing.expectEqual(0x3.6a6846287159ceb1f8d965deea72p-4, log2(@as(f128, 0x1.28d3b4p+0)));
    // try std.testing.expectEqual(-0x1.b8bdee5bd428254ebf7ead086211p-4, log2(@as(f128, 0xe.d99dap-4)));
    try std.testing.expectEqual(0x7.999f2a94857b22e23de191021e94p-4, log2(@as(f128, 0x1.63d204p+0)));
    // try std.testing.expectEqual(0x7.999f095d7e150199097308f7442p-4, log2(@as(f128, 0x1.63d202p+0)));
    try std.testing.expectEqual(0x7.999f16e03b54e0f55354326ba2c8p-4, log2(@as(f128, 0x1.63d202d04392cp+0)));
    try std.testing.expectEqual(-0x3.75cc3f2233d65c71f491713f54c8p-8, log2(@as(f128, 0xf.d9ce1p-4)));
    // try std.testing.expectEqual(-0x3.75cdb3f0a0b655972857a48ecdfap-8, log2(@as(f128, 0xf.d9cep-4)));
    try std.testing.expectEqual(-0x3.75ccb13d89cc66ee72a7e068c0eap-8, log2(@as(f128, 0xf.d9ce0b1a50e08p-4)));
    // try std.testing.expectEqual(0xa.5959448ade0bdcb61198f8dd2668p-8, log2(@as(f128, 0x1.07465cp+0)));
    try std.testing.expectEqual(0xa.5956764b14a9d7da921863b84bbp-8, log2(@as(f128, 0x1.07465ap+0)));
    // try std.testing.expectEqual(0xa.595912bb7fff0ba2d9d999b5e918p-8, log2(@as(f128, 0x1.07465bdc7e41cp+0)));
    // try std.testing.expectEqual(0xa.595912bb7fe899a4a195bab67568p-8, log2(@as(f128, 0x1.07465bdc7e41bp+0)));
    try std.testing.expectEqual(0xa.595912bb7fefddcbce0cb3878e3p-8, log2(@as(f128, 0x1.07465bdc7e41b52ep+0)));
    try std.testing.expectEqual(-0x2.c106931f2bfd0af427fc474396bp+4, log2(@as(f128, 0xf.4dfb4p-48)));
    try std.testing.expectEqual(0xe.a1dd43a221d02a32622e9cba02cp-8, log2(@as(f128, 0x1.0a588ep+0)));
    // try std.testing.expectEqual(-0x6.d35688edc44a495fd74b5e1b8d9cp-4, log2(@as(f128, 0xb.e77c6p-4)));
    try std.testing.expectEqual(0x6.44f92e0fda7d1b46e2bc2dcfa988p-4, log2(@as(f128, 0x1.4fe37ep+0)));
    try std.testing.expectEqual(0x9.d9a8c6de3432bff9b0fef9633f8p+0, log2(@as(f128, 0x3.9b0754p+8)));
    // try std.testing.expectEqual(-0x6.df8b2c2c5ea400001520bc941b08p-4, log2(@as(f128, 0xb.e132ap-4)));
    // try std.testing.expectEqual(-0x7.e842050c531d02307168b728f45cp-4, log2(@as(f128, 0xb.5bf83p-4)));
    // try std.testing.expectEqual(-0x7.e842258fcc5d9f2cd3d1ec3f0854p-4, log2(@as(f128, 0xb.5bf82p-4)));
    try std.testing.expectEqual(-0x7.e8420994680da57678fdcfd1bea8p-4, log2(@as(f128, 0xb.5bf82dc51f028p-4)));
    // try std.testing.expectEqual(-0x7.e8420994680ea992427e98066258p-4, log2(@as(f128, 0xb.5bf82dc51f02p-4)));
    // try std.testing.expectEqual(-0x7.e8420994680ea2d70a67a2d80578p-4, log2(@as(f128, 0xb.5bf82dc51f02035p-4)));
    // try std.testing.expectEqual(-0x7.b18b6b68ffa235098af8c5c4f5d4p-4, log2(@as(f128, 0xb.7704ep-4)));
    // try std.testing.expectEqual(-0x7.b18b8b9fc309de4f9564e0281104p-4, log2(@as(f128, 0xb.7704dp-4)));
    try std.testing.expectEqual(-0x7.b18b723cc4c4faeb8adda8e96a08p-4, log2(@as(f128, 0xb.7704dc9beb05p-4)));
    // try std.testing.expectEqual(-0x7.f27149af9dc8b0f1993d141d1a5p-4, log2(@as(f128, 0xb.56f64p-4)));
    // try std.testing.expectEqual(-0x7.f2716a4172a70437981d6d2faa0cp-4, log2(@as(f128, 0xb.56f63p-4)));
    // try std.testing.expectEqual(-0x7.f27151a15d70d15a62a4c18ae804p-4, log2(@as(f128, 0xb.56f63c18e93fp-4)));
    // try std.testing.expectEqual(-0x7.f27151a15d71d5e90939366e0698p-4, log2(@as(f128, 0xb.56f63c18e93e8p-4)));
    try std.testing.expectEqual(-0x7.f27151a15d70f86944dd429073p-4, log2(@as(f128, 0xb.56f63c18e93eecdp-4)));
    try std.testing.expectEqual(-0x7.f84a9424a2fba583bbaf3ba7a54p-4, log2(@as(f128, 0xb.54171p-4)));
    // try std.testing.expectEqual(-0x7.f84ab4beb988b35b6eb65af72e58p-4, log2(@as(f128, 0xb.5417p-4)));
    // try std.testing.expectEqual(-0x7.f84a998412a436a8365c4397d8f8p-4, log2(@as(f128, 0xb.54170d5cfa9p-4)));
    // try std.testing.expectEqual(-0x7.f84a998412a53b78ea4938468dfcp-4, log2(@as(f128, 0xb.54170d5cfa8f8p-4)));
    // try std.testing.expectEqual(-0x7.f84a998412a489d141bab11c0c3p-4, log2(@as(f128, 0xb.54170d5cfa8fd73p-4)));
    try std.testing.expectEqual(-0x7.f84a998412a489f1dbd12ebaa208p-4, log2(@as(f128, 0xb.54170d5cfa8fd72p-4)));
    // try std.testing.expectEqual(-0x7.f84a998412a489dce921cd17e614p-4, log2(@as(f128, 0xb.54170d5cfa8fd72a47d6bda19068p-4)));
    // try std.testing.expectEqual(-0x7.f84a998412a489dce921cd17decp-4, log2(@as(f128, 0xb.54170d5cfa8fd72a47d6bda194p-4)));
    try std.testing.expectEqual(-0x7.f84a998412a489dce921cd17e6e8p-4, log2(@as(f128, 0xb.54170d5cfa8fd72a47d6bda19p-4)));
    try std.testing.expectEqual(-0x7.ep+4, log2(@as(f128, 0x4p-128)));
    try std.testing.expectEqual(-0x3.fep+8, log2(@as(f128, 0x4p-1024)));
    try std.testing.expectEqual(-0x3.ffep+12, log2(@as(f128, 0x4p-16384)));
    try std.testing.expectEqual(-0x3.fffp+12, log2(@as(f128, 0x2p-16384)));
    try std.testing.expectEqual(-0x3.c9p+8, log2(@as(f128, 0x8p-972)));
    try std.testing.expectEqual(-0x9.5p+4, log2(@as(f128, 0x8p-152)));
    try std.testing.expectEqual(-0x4.32p+8, log2(@as(f128, 0x4p-1076)));
    try std.testing.expectEqual(-0x4.03dp+12, log2(@as(f128, 0x8p-16448)));
    try std.testing.expectEqual(-0x4.03ep+12, log2(@as(f128, 0x4p-16448)));
    try std.testing.expectEqual(-0x4.06ep+12, log2(@as(f128, 0x4p-16496)));
    try std.testing.expectEqual(0x7.ffffffe8eab88f49d947a104332p+4, log2(@as(f128, 0xf.fffffp+124)));
    try std.testing.expectEqual(0x3.fffffffffffffff4755c4d6a3e8p+8, log2(@as(f128, 0xf.ffffffffffff8p+1020)));
    try std.testing.expectEqual(0x3.ffffffffffffffffffe8eab89ad4p+12, log2(@as(f128, 0xf.fffffffffffffffp+16380)));
    try std.testing.expectEqual(0x4p+12, log2(@as(f128, 0xf.fffffffffffffffffffffffffff8p+16380)));
    try std.testing.expectEqual(0x3.fffffffffffffffa3aae26b51f4p+8, log2(@as(f128, 0xf.ffffffffffffbffffffffffffcp+1020)));
}
