const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
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
        return (-1 / float.abs(x)); // log2l(+-0)=-inf

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
    var xx: f128 = float.frexp(x, &e);

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
