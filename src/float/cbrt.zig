const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const classify = @import("classify.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn cbrt(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return cbrt(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, cbrt32(cast(f32, x))),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_cbrtf.c
                    return cbrt32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_cbrt.c
                    return cbrt64(x);
                },
                f80 => return cast(f80, cbrt128(cast(f128, x, .{})), .{}),
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_cbrtl.c
                    return cbrt128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn cbrt32(x: f32) f32 {
    const escale: [3]f64 = .{
        1.0,
        0x1.428a2f98d728bp+0, // 2^(1/3)
        0x1.965fea53d6e3dp+0, // 2^(2/3)
    };

    const u: u32 = @bitCast(x);
    var au: u32 = u << 1;
    const sgn: u32 = u >> 31;
    var e: u32 = au >> 24;
    if (au < 1 << 24 or au >= 0xff << 24) {
        @branchHint(.unlikely);
        if (au >= 0xff << 24)
            return x + x; // inf, nan

        if (au == 0)
            return x; // +-0

        const nz: i32 = @clz(au) - 7; // subnormal
        au <<= @as(u5, @intCast(nz));
        e -%= cast(u32, nz -% 1, .{});
    }

    const mant: u32 = au & 0xffffff;
    e +%= 899;
    const et: u32 = @divFloor(e, 3);
    const it: u32 = e % 3;
    var isc: u64 = @bitCast(escale[it]);
    isc +%= cast(u64, (et -% 342), .{}) << 52;
    isc |= cast(u64, sgn, .{}) << 63;
    const cvt2: f64 = @bitCast(isc);
    const c: [8]f64 = .{
        0x1.2319d352ea5d5p-1,  0x1.67ad8ee258d1ap-1,  -0x1.9342edf9cbad9p-2,
        0x1.b6388fc510a75p-3,  -0x1.6002455599e2fp-4, 0x1.7b096936192c4p-6,
        -0x1.e5577187e8bf8p-9, 0x1.169ef81d6c34ep-12,
    };
    const z: f64 = @bitCast(cast(u64, mant, .{}) << 28 | 0x3ff << 52);
    const r0: f64 = -0x1.9931c6c2d19d1p-6 / z;
    const z2: f64 = z * z;
    const z4: f64 = z2 * z2;
    var f: f64 = ((c[0] + z * c[1]) + z2 * (c[2] + z * c[3])) + z4 * ((c[4] + z * c[5]) + z2 * (c[6] + z * c[7])) + r0;
    var r: f64 = f * cvt2;
    var ub: f32 = cast(f32, r, .{});
    const lb: f32 = cast(f32, r - cvt2 * 1.4182e-9, .{});
    if (ub == lb) {
        @branchHint(.likely);
        return ub;
    }

    const U0: f64 = -0x1.ab16ec65d138fp+3;
    const h: f64 = f * f * f - z;
    f -= (f * r0 * U0) * h;
    r = f * cvt2;
    var cvt1: u64 = @bitCast(r);
    ub = cast(f32, r, .{});
    const m0: i64 = @bitCast(cvt1 << 19);
    const m1: i64 = m0 >> 63;
    if ((m0 ^ m1) < (1 << 31)) {
        @branchHint(.unlikely);
        cvt1 = (cvt1 + (1 << 31)) & 0xffffffff00000000;
        ub = cast(f32, @as(f64, @bitCast(cvt1)), .{});
    }
    return ub;
}

fn cbrt64(x: f64) f64 {
    const factor: [5]f64 = .{ 1.0 / 1.5874010519681994748, 1.0 / 1.2599210498948731648, 1.0, 1.2599210498948731648, 1.5874010519681994748 };

    // Reduce X.  XM now is an range 1.0 to 0.5.
    var xe: i32 = undefined;
    const xm = float.frexp(float.abs(x), &xe);

    // If X is not finite or is null return it (with raising exceptions
    // if necessary.
    if (xe == 0 and classify.classify(x) <= classify.ZERO)
        return x + x;

    const u: f64 = (0.354895765043919860 + ((1.50819193781584896 + ((-2.11499494167371287 + ((2.44693122563534430 + ((-1.83469277483613086 + (0.784932344976639262 - 0.145263899385486377 * xm) * xm) * xm)) * xm)) * xm)) * xm));

    const t2: f64 = u * u * u;

    const ym: f64 = u * (t2 + 2.0 * xm) / (2.0 * t2 + xm) * factor[@intCast(2 + @mod(xe, 3))];

    return float.ldexp(if (x > 0) ym else -ym, @divFloor(xe, 3));
}

fn cbrt128(x: f128) f128 {
    const CBRT2: f128 = 1.259921049894873164767210607278228350570251;
    const CBRT4: f128 = 1.587401051968199474751705639272308260391493;
    const CBRT2I: f128 = 0.7937005259840997373758528196361541301957467;
    const CBRT4I: f128 = 0.6299605249474365823836053036391141752851257;

    if (!std.math.isFinite(x))
        return x + x;

    if (x == 0)
        return x;

    var sign: i32 = undefined;
    var xx: f128 = x;
    if (x > 0) {
        sign = 1;
    } else {
        sign = -1;
        xx = -x;
    }

    const z: f128 = xx;
    // extract power of 2, leaving mantissa between 0.5 and 1
    var e: i32 = undefined;
    xx = float.frexp(xx, &e);

    // Approximate cube root of number between .5 and 1,
    // peak relative error = 1.2e-6
    xx = ((((1.3584464340920900529734e-1 * xx - 6.3986917220457538402318e-1) * xx + 1.2875551670318751538055e0) * xx - 1.4897083391357284957891e0) * xx + 1.3304961236013647092521e0) * xx + 3.7568280825958912391243e-1;

    // exponent divided by 3
    if (e >= 0) {
        var rem: i32 = e;
        e = @divFloor(e, 3);
        rem -= 3 * e;
        if (rem == 1) {
            xx *= CBRT2;
        } else if (rem == 2) {
            xx *= CBRT4;
        }
    } else { // argument less than 1
        e = -e;
        var rem: i32 = e;
        e = @divFloor(e, 3);
        rem -= 3 * e;
        if (rem == 1) {
            xx *= CBRT2I;
        } else if (rem == 2) {
            xx *= CBRT4I;
        }
        e = -e;
    }

    // multiply by power of 2
    xx = float.ldexp(xx, e);

    // Newton iteration
    xx -= (xx - (z / (xx * xx))) * 0.3333333333333333333333333333333333333333;
    xx -= (xx - (z / (xx * xx))) * 0.3333333333333333333333333333333333333333;
    xx -= (xx - (z / (xx * xx))) * 0.3333333333333333333333333333333333333333;

    if (sign < 0)
        xx = -xx;

    return xx;
}
