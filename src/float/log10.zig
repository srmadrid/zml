const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");
const erf_data = @import("erf_data.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn log10(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return log10(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, log10_32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_log10f.c
                    return log10_32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_log10.c
                    return log10_64(x);
                },
                f80 => return cast(f80, log10_128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/e_log10l.c
                    return log10_128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn as_special(x: f32) f32 {
    const ux: u32 = @bitCast(x);
    if (ux == 0x7f800000)
        return x; // +inf

    const ax: u32 = ux << 1;
    if (ax == 0) {
        return -1 / @as(f32, 0); // -0
    }

    if (ax > 0xff000000)
        return x + x; // nan

    return (x - x) / (x - x);
}

fn log10_32(x: f32) f32 {
    const tr: [65]f64 = .{
        0x1p+0,         0x1.f81f82p-1,  0x1.f07c1fp-1,  0x1.e9131acp-1,
        0x1.e1e1e1ep-1, 0x1.dae6077p-1, 0x1.d41d41dp-1, 0x1.cd85689p-1,
        0x1.c71c71cp-1, 0x1.c0e0704p-1, 0x1.bacf915p-1, 0x1.b4e81b5p-1,
        0x1.af286bdp-1, 0x1.a98ef6p-1,  0x1.a41a41ap-1, 0x1.9ec8e95p-1,
        0x1.999999ap-1, 0x1.948b0fdp-1, 0x1.8f9c19p-1,  0x1.8acb90fp-1,
        0x1.8618618p-1, 0x1.8181818p-1, 0x1.7d05f41p-1, 0x1.78a4c81p-1,
        0x1.745d174p-1, 0x1.702e05cp-1, 0x1.6c16c17p-1, 0x1.6816817p-1,
        0x1.642c859p-1, 0x1.605816p-1,  0x1.5c9882cp-1, 0x1.58ed231p-1,
        0x1.5555555p-1, 0x1.51d07ebp-1, 0x1.4e5e0a7p-1, 0x1.4afd6ap-1,
        0x1.47ae148p-1, 0x1.446f865p-1, 0x1.4141414p-1, 0x1.3e22cbdp-1,
        0x1.3b13b14p-1, 0x1.3813814p-1, 0x1.3521cfbp-1, 0x1.323e34ap-1,
        0x1.2f684bep-1, 0x1.2c9fb4ep-1, 0x1.29e412ap-1, 0x1.27350b9p-1,
        0x1.2492492p-1, 0x1.21fb781p-1, 0x1.1f7047ep-1, 0x1.1cf06aep-1,
        0x1.1a7b961p-1, 0x1.1811812p-1, 0x1.15b1e5fp-1, 0x1.135c811p-1,
        0x1.1111111p-1, 0x1.0ecf56cp-1, 0x1.0c9715p-1,  0x1.0a6810ap-1,
        0x1.0842108p-1, 0x1.0624dd3p-1, 0x1.041041p-1,  0x1.0204081p-1,
        0.5,
    };
    const tl: [65]f64 = .{
        -0x1.d45fd6237ebe3p-47, 0x1.b947689311b6ep-8, 0x1.b5e909c96d7d5p-7,
        0x1.45f4f59ed2165p-6,   0x1.af5f92cbd8f1ep-6, 0x1.0ba01a606de8cp-5,
        0x1.3ed119b9a2b7bp-5,   0x1.714834298eec2p-5, 0x1.a30a9d98357fbp-5,
        0x1.d41d512670813p-5,   0x1.02428c0f65519p-4, 0x1.1a23444eecc3ep-4,
        0x1.31b30543f4cb4p-4,   0x1.48f3ed39bfd04p-4, 0x1.5fe8049a0e423p-4,
        0x1.769140a6aa008p-4,   0x1.8cf1836c98cb3p-4, 0x1.a30a9d55541a1p-4,
        0x1.b8de4d1ee823ep-4,   0x1.ce6e4202ca2e6p-4, 0x1.e3bc1accace07p-4,
        0x1.f8c9683b5abd4p-4,   0x1.06cbd68ca9a6ep-3, 0x1.11142f19df73p-3,
        0x1.1b3e71fa7a97fp-3,   0x1.254b4d37a46e3p-3, 0x1.2f3b6912cbf07p-3,
        0x1.390f683115886p-3,   0x1.42c7e7fffc5a8p-3, 0x1.4c65808c78d3cp-3,
        0x1.55e8c50751c55p-3,   0x1.5f52445dec3d8p-3, 0x1.68a288c3f12p-3,
        0x1.71da17bdf0d19p-3,   0x1.7af973608afd9p-3, 0x1.84011952a2579p-3,
        0x1.8cf1837a7ea6p-3,    0x1.95cb2891e43d6p-3, 0x1.9e8e7b0f869ep-3,
        0x1.a73beaa5db18dp-3,   0x1.afd3e394558d3p-3, 0x1.b856cf060d9f1p-3,
        0x1.c0c5134de1ffcp-3,   0x1.c91f1371bc99fp-3, 0x1.d1652ffcd3f53p-3,
        0x1.d997c6f635e75p-3,   0x1.e1b733ab90f3bp-3, 0x1.e9c3ceadac856p-3,
        0x1.f1bdeec43a305p-3,   0x1.f9a5e7a5fa3fep-3, 0x1.00be05ac02f2bp-2,
        0x1.04a054d81a2d4p-2,   0x1.087a0835957fbp-2, 0x1.0c4b457099517p-2,
        0x1.101431aa1fe51p-2,   0x1.13d4f08b98dd8p-2, 0x1.178da53edb892p-2,
        0x1.1b3e71e9f9d58p-2,   0x1.1ee777defdeedp-2, 0x1.2288d7b48e23bp-2,
        0x1.2622b0f52e49fp-2,   0x1.29b522a4c6314p-2, 0x1.2d404b0e30f8p-2,
        0x1.30c4478f3fbe5p-2,   0x1.34413509f7915p-2,
    };
    const st: [16]f32 = .{
        0x1p+0,        0x1.4p+3,      0x1.9p+6,       0x1.f4p+9,
        0x1.388p+13,   0x1.86ap+16,   0x1.e848p+19,   0x1.312dp+23,
        0x1.7d784p+26, 0x1.dcd65p+29, 0x1.2a05f2p+33, 0,
        0,             0,             0,              0,
    };
    const b: [3]f64 = .{ 0x1.bcb7b15c5a2f8p-2, -0x1.bcbb1dbb88ebap-3, 0x1.2871c39d521c6p-3 };
    const c: [7]f64 = .{
        0x1.bcb7b1526e50ep-2,  -0x1.bcb7b1526e53dp-3, 0x1.287a7636f3fa2p-3,
        -0x1.bcb7b146a14b3p-4, 0x1.63c627d5219cbp-4,  -0x1.2880736c8762dp-4,
        0x1.fc1ecf913961ap-5,
    };

    var ux: u32 = @bitCast(x);
    if (ux < (1 << 23) or ux >= 0x7f800000) {
        @branchHint(.unlikely);
        if (ux == 0 or ux >= 0x7f800000)
            return as_special(x);

        // subnormal
        const n: u32 = @clz(ux) - 8;
        ux <<= @as(u5, @intCast(n));
        ux -%= @as(u32, @intCast(n)) << 23;
    }

    const m: u32 = ux & ((1 << 23) - 1);
    const j: u32 = (m + (1 << (23 - 7))) >> (23 - 6);
    const ix: f64 = tr[j];
    const l: f64 = tl[j];
    var e: i32 = undefined;
    {
        @setRuntimeSafety(false);
        e = (@as(i32, @intCast(ux)) >> 23) - 127;
    }
    var je: u32 = undefined;
    {
        @setRuntimeSafety(false);
        je = @intCast(e + 1);
    }
    je = (je *% 0x4d104d4) >> 28;
    if (ux == @as(u32, @bitCast(st[je]))) {
        @branchHint(.unlikely);
        return cast(f32, je, .{});
    }

    var tz: f64 = @bitCast((cast(i64, m, .{}) | (1023 << 23)) << (52 - 23));
    const z: f64 = tz * ix - 1;
    const z2: f64 = z * z;
    var r: f64 = ((cast(f64, e, .{}) * 0x1.34413509f79ffp-2 + l) + z * b[0]) + z2 * (b[1] + z * b[2]);
    var ub: f32 = cast(f32, r, .{});
    const lb: f32 = cast(f32, r + 0x1.b008p-34, .{});
    if (ub != lb) {
        @branchHint(.unlikely);
        var f: f64 = z * ((c[0] + z * c[1]) + z2 * ((c[2] + z * c[3]) + z2 * (c[4] + z * c[5] + z2 * c[6])));
        f -= 0x1.0cee0ed4ca7e9p-54 * cast(f64, e, .{});
        f += l - tl[0];
        const el: f64 = cast(f64, e, .{}) * 0x1.34413509f7ap-2;
        r = el + f;
        ub = cast(f32, r, .{});
        tz = r;
        if (@as(u64, @bitCast(tz)) & ((1 << 28) - 1) == 0) {
            @branchHint(.unlikely);
            const dr: f64 = (el - r) + f;
            r += dr * 32;
            ub = cast(f32, r, .{});
        }
    }
    return ub;
}

fn log10_64(x: f64) f64 {
    const two54: f64 = 1.80143985094819840000e+16; // 0x4350000000000000
    const ivln10: f64 = 4.34294481903251816668e-01; // 0x3fdbcb7b1526e50e
    const log10_2hi: f64 = 3.01029995663611771306e-01; // 0x3fd34413509f6000
    const log10_2lo: f64 = 3.69423907715893078616e-13; // 0x3d59fef311f12b36

    var hx: i64 = undefined;
    dbl64.extractWords64(&hx, x);

    var k: i32 = 0;
    if (hx < 0x0010000000000000) { // x < 2**-1022
        if ((hx & 0x7fffffffffffffff) == 0) {
            @branchHint(.unlikely);
            return -two54 / float.abs(x); // log(+-0)=-inf
        }

        if (hx < 0) {
            @branchHint(.unlikely);
            return (x - x) / (x - x); // log(-#) = NaN
        }

        k -= 54;
        const tmp: f64 = x * two54; // subnormal number, scale up x
        dbl64.extractWords64(&hx, tmp);
    }

    // scale up resulted in a NaN number
    if (hx >= 0x7ff0000000000000) {
        @branchHint(.unlikely);
        return x + x;
    }

    k += cast(i32, (hx >> 52) - 1023, .{});
    var i: i64 = undefined;
    {
        @setRuntimeSafety(false);
        i = cast(i64, (@as(u64, @intCast(k)) & 0x8000000000000000) >> 63, .{});
    }
    hx = (hx & 0x000fffffffffffff) | ((0x3ff - i) << 52);
    var y: f64 = cast(f64, k + i, .{});
    if (y == 0)
        y = 0;

    var xx: f64 = x;
    dbl64.insertWords64(&xx, hx);
    const z: f64 = y * log10_2lo + ivln10 * float.log(xx);
    return z + y * log10_2hi;
}

fn log10_128(x: f128) f128 {
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

    // log10(2)
    const L102A: f128 = 0.3125;
    const L102B: f128 = -1.14700043360188047862611052755069732318101185e-2;
    // log10(e)
    const L10EA: f128 = 0.5;
    const L10EB: f128 = -6.570551809674817234887108108339491770560299e-2;
    // sqrt(2)/2
    const SQRTH: f128 = 7.071067811865475244008443621048490392848359e-1;

    // Test for domain
    var hx: i64 = undefined;
    var lx: i64 = undefined;
    ldbl128.getWords(&hx, &lx, x);
    if (((hx & 0x7fffffffffffffff) | lx) == 0)
        return (-1 / float.abs(x)); // log10l(+-0)=-inf

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

    // Multiply log of fraction by log10(e)
    // and base 2 exponent by log10(2).
    z = y * L10EB;
    z += xx * L10EB;
    z += cast(f128, e, .{}) * L102B;
    z += y * L10EA;
    z += xx * L10EA;
    z += cast(f128, e, .{}) * L102A;
    return z;
}
