const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const flt32 = @import("flt32.zig");
const dbl64 = @import("dbl64.zig");
const ldbl80 = @import("ldbl80.zig");
const ldbl128 = @import("ldbl128.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn scalbn(x: anytype, n: i32) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return scalbn(cast(EnsureFloat(@TypeOf(x)), x, .{}), n);
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return cast(f16, scalbn32(cast(f32, x), n)),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/s_scalbnf.c
                    return scalbn32(x, n);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/s_scalbn.c
                    return scalbn64(x, n);
                },
                f80 => {
                    // glibc/sysdeps/ieee754/ldbl-96/s_scalbnl.c
                    return scalbn80(x, n);
                },
                f128 => {
                    // glibc/sysdeps/ieee754/ldbl-128/s_scalbnl.c
                    return scalbn128(x, n);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn scalbn32(x: f32, n: i32) f32 {
    const two25: f32 = 3.355443200e+07; // 0x4c000000
    const twom25: f32 = 2.9802322388e-08; // 0x33000000
    const huge: f32 = 1.0e+30;
    const tiny: f32 = 1.0e-30;

    var ix: i32 = undefined;
    flt32.getWord(&ix, x);
    var k: i32 = (ix & 0x7f800000) >> 23; // extract exponent
    var xx: f32 = x;
    if (k == 0) { // 0 or subnormal x
        @branchHint(.unlikely);

        if ((ix & 0x7fffffff) == 0)
            return xx; // +-0

        xx *= two25;
        flt32.getWord(&ix, xx);
        k = ((ix & 0x7f800000) >> 23) - 25;
    }

    if (k == 0xff) {
        @branchHint(.unlikely);
        return xx + xx; // NaN or Inf
    }

    if (n < -50000) {
        @branchHint(.unlikely);
        return tiny * float.copysign(tiny, xx); // underflow
    }

    if (n > 50000 or k + n > 0xfe) {
        @branchHint(.unlikely);
        return huge * float.copysign(huge, xx); // overflow
    }

    // Now k and n are bounded we know that k = k+n does not overflow.
    k = k + n;
    if (k > 0) { // normal result
        @branchHint(.likely);
        flt32.setWord(&xx, (@as(u32, @bitCast(ix)) & 0x807fffff) | (@as(u32, @bitCast(k)) << 23));
        return xx;
    }

    if (k <= -25)
        return tiny * float.copysign(tiny, xx); // underflow

    k += 25; // subnormal result
    flt32.setWord(&xx, (@as(u32, @bitCast(ix)) & 0x807fffff) | (@as(u32, @bitCast(k)) << 23));
    return xx * twom25;
}

fn scalbn64(x: f64, n: i32) f64 {
    const two54: f64 = 1.80143985094819840000e+16; // 0x43500000, 0x00000000
    const twom54: f64 = 5.55111512312578270212e-17; // 0x3C900000, 0x00000000
    const huge: f64 = 1.0e+300;
    const tiny: f64 = 1.0e-300;

    var ix: i64 = undefined;
    dbl64.extractWords64(&ix, x);
    var k: i64 = (ix >> 52) & 0x7ff; // extract exponent
    var xx: f64 = x;
    if (k == 0) { // 0 or subnormal x
        @branchHint(.unlikely);
        if ((ix & 0xfffffffffffff) == 0)
            return xx; // +-0

        xx *= two54;
        dbl64.extractWords64(&ix, xx);
        k = ((ix >> 52) & 0x7ff) - 54;
    }

    if (k == 0x7ff) {
        @branchHint(.unlikely);
        return xx + xx; // NaN or Inf
    }

    if (n < -50000) {
        @branchHint(.unlikely);
        return tiny * float.copysign(tiny, xx); // underflow
    }

    if (n > 50000 or k + n > 0x7fe) {
        @branchHint(.unlikely);
        return huge * float.copysign(huge, xx); // overflow
    }

    // Now k and n are bounded we know that k = k+n does not overflow.
    k = k + n;
    if (k > 0) { // normal result
        @branchHint(.likely);
        dbl64.insertWords64(&xx, (@as(u64, @bitCast(ix)) & 0x800fffffffffffff) | (@as(u64, @bitCast(k)) << 52));
        return xx;
    }

    if (k <= -54)
        return tiny * float.copysign(tiny, xx); // underflow

    k += 54; // subnormal result
    dbl64.insertWords64(&xx, (@as(u64, @bitCast(ix)) & 0x800fffffffffffff) | (@as(u64, @bitCast(k)) << 52));
    return xx * twom54;
}

fn scalbn80(x: f80, n: i32) f80 {
    const two63: f80 = 0x1p63;
    const twom64: f80 = 0x1p-64;
    const huge: f80 = 1.0e+4900;
    const tiny: f80 = 1.0e-4900;

    var es: i32 = undefined;
    var hx: i32 = undefined;
    var lx: i32 = undefined;
    ldbl80.getWords(&es, &hx, &lx, x);
    var k: i32 = es & 0x7fff; // extract exponent
    var xx: f80 = x;
    if (k == 0) { // 0 or subnormal x
        @branchHint(.unlikely);
        if ((lx | (hx & 0x7fffffff)) == 0)
            return xx; // +-0

        xx *= two63;
        ldbl80.getExp(&es, xx);
        k = (es & 0x7fff) - 63;
    }

    if (k == 0x7fff) {
        @branchHint(.unlikely);
        return xx + xx; // NaN or Inf
    }

    if (n < -50000) {
        @branchHint(.unlikely);
        return tiny * float.copysign(tiny, xx);
    }

    if (n > 50000 or k + n > 0x7ffe) {
        @branchHint(.unlikely);
        return huge * float.copysign(huge, xx); // overflow
    }

    // Now k and n are bounded we know that k = k+n does not overflow.
    k = k + n;
    if (k > 0) { // normal result
        @branchHint(.likely);
        ldbl80.setExp(&xx, (es & 0x8000) | k);
        return xx;
    }

    if (k <= -64)
        return tiny * float.copysign(tiny, xx); // underflow

    k += 64; // subnormal result
    ldbl80.setExp(&xx, (es & 0x8000) | k);
    return xx * twom64;
}

fn scalbn128(x: f128, n: i32) f128 {
    const two114: f128 = 2.0769187434139310514121985316880384e+34; // 0x4071000000000000, 0
    const twom114: f128 = 4.8148248609680896326399448564623183e-35; // 0x3F8D000000000000, 0
    const huge: f128 = 1.0e+4900;
    const tiny: f128 = 1.0e-4900;

    var hx: i64 = undefined;
    var lx: i64 = undefined;
    ldbl128.getWords(&hx, &lx, x);
    var k: i64 = (hx >> 48) & 0x7fff; // extract exponent
    var xx: f128 = x;
    if (k == 0) { // 0 or subnormal x
        if (lx | (hx & 0x7fffffffffffffff) == 0)
            return xx; // +-0

        xx *= two114;
        ldbl128.getMsw(&hx, xx);
        k = ((hx >> 48) & 0x7fff) - 114;
    }

    if (k == 0x7fff)
        return xx + xx; // NaN or Inf

    if (n < -50000)
        return tiny * float.copysign(tiny, xx); // underflow

    if (n > 50000 or k + n > 0x7ffe)
        return huge * float.copysign(huge, xx); // overflow

    // Now k and n are bounded we know that k = k+n does not overflow.
    k = k + n;
    if (k > 0) { // normal result
        ldbl128.setMsw(&xx, (@as(u64, @bitCast(hx)) & 0x8000ffffffffffff) | (@as(u64, @bitCast(k)) << 48));
        return xx;
    }

    if (k <= -114)
        return tiny * float.copysign(tiny, xx); // underflow

    k += 114; // subnormal result
    ldbl128.setMsw(&xx, (@as(u64, @bitCast(hx)) & 0x8000ffffffffffff) | (@as(u64, @bitCast(k)) << 48));
    return xx * twom114;
}
