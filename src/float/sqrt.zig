const std = @import("std");
const types = @import("../types.zig");
const dla = @import("dla.zig");
const root = @import("root.zig");
const flt32 = @import("flt32.zig");
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;
const ta = @import("builtin").target.cpu.arch;

fn has_hardware_sqrt(comptime target_arch: std.Target.Cpu.Arch, comptime T: type) bool {
    if (types.numericType(T) != .float)
        @compileError("T must be a float");

    const features = @import("builtin").target.cpu.features;

    return switch (target_arch) {
        .amdgcn => switch (T) {
            f16 => true,
            f32 => true,
            f64 => std.Target.amdgcn.featureSetHas(features, .fp64),
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .arc => false,
        .arm, .armeb, .thumb, .thumbeb => switch (T) {
            f16 => std.Target.arm.featureSetHas(features, .fp16),
            f32 => std.Target.arm.featureSetHas(features, .vfp2),
            f64 => std.Target.arm.featureSetHas(features, .vfp3),
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .aarch64, .aarch64_be => switch (T) {
            f16 => std.Target.aarch64.featureSetHas(features, .fp16),
            f32 => true,
            f64 => true,
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .avr => false,
        .bpfeb, .bpfel => false,
        .csky => false,
        .hexagon => switch (T) {
            f16 => false,
            f32 => std.Target.hexagon.featureSetHas(features, .hvx),
            f64 => false,
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .kalimba => false,
        .lanai => false,
        .loongarch32, .loongarch64 => switch (T) {
            f16 => false,
            f32 => true,
            f64 => true,
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .m68k => false,
        .mips, .mipsel, .mips64, .mips64el => switch (T) {
            f16 => false,
            f32 => std.Target.mips.featureSetHas(features, .fp64),
            f64 => std.Target.mips.featureSetHas(features, .fp64),
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .msp430 => false,
        .nvptx, .nvptx64 => false,
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => switch (T) {
            f16 => false,
            f32 => true,
            f64 => true,
            f80 => false,
            f128 => std.Target.powerpc.featureSetHas(features, .vsx),
            else => unreachable,
        },
        .propeller => false,
        .riscv32, .riscv64 => switch (T) {
            f16 => std.Target.riscv.featureSetHas(features, .zfh),
            f32 => std.Target.riscv.featureSetHas(features, .f),
            f64 => std.Target.riscv.featureSetHas(features, .d),
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .s390x => switch (T) {
            f16 => false,
            f32 => true,
            f64 => true,
            f80 => false,
            f128 => std.Target.s390x.featureSetHas(features, .vector_enhancements_1),
            else => unreachable,
        },
        .sparc, .sparc64 => switch (T) {
            f16 => false,
            f32 => true,
            f64 => true,
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .spirv, .spirv32, .spirv64 => false,
        .ve => false,
        .wasm32, .wasm64 => switch (T) {
            f16 => false,
            f32 => true,
            f64 => true,
            f80 => false,
            f128 => false,
            else => unreachable,
        },
        .x86, .x86_64 => switch (T) {
            f16 => std.Target.x86.featureSetHas(features, .avx512fp16),
            f32 => std.Target.x86.featureSetHas(features, .sse),
            f64 => std.Target.x86.featureSetHas(features, .sse2),
            f80 => std.Target.x86.featureSetHas(features, .x87),
            f128 => false,
            else => unreachable,
        },
        .xcore => false,
        .xtensa => false,
    };
}

pub inline fn sqrt(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return sqrt(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return if (has_hardware_sqrt(ta, f16)) @sqrt(x) else cast(f16, sqrt32(cast(f32, x, .{})), .{}),
                f32 => {
                    // glibc/sysdeps/ieee754/flt-32/e_sqrtf.c
                    return if (has_hardware_sqrt(ta, f32)) @sqrt(x) else sqrt32(x);
                },
                f64 => {
                    // glibc/sysdeps/ieee754/dbl-64/e_sqrt.c
                    return if (has_hardware_sqrt(ta, f64)) @sqrt(x) else sqrt64(x);
                },
                f80 => return if (has_hardware_sqrt(ta, f80)) @sqrt(x) else cast(f80, sqrt128(cast(f128, x, .{})), .{}),
                f128 => {
                    // Adapted from: glibc/sysdeps/ieee754/ldbl-128ibm/e_sqrtl.c
                    return if (has_hardware_sqrt(ta, f128)) @sqrt(x) else sqrt128(x);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn sqrt32(x: f32) f32 {
    var ix: i32 = undefined;
    flt32.getWord(&ix, x);
    const sign: u32 = 0x80000000;

    // take care of Inf and NaN
    if ((ix & 0x7f800000) == 0x7f800000) {
        return x * x + x; // sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN
    }

    // take care of zero
    if (ix <= 0) {
        if ((ix & (~sign)) == 0) {
            return x; // sqrt(+-0) = +-0
        } else if (ix < 0) {
            return (x - x) / (x - x); // sqrt(-ve) = sNaN
        }
    }

    // normalize x
    var m: i32 = (ix >> 23);
    if (m == 0) { // subnormal x
        var i: i32 = 0;
        while ((ix & 0x00800000) == 0) {
            ix <<= 1;

            i += 1;
        }

        m -= i - 1;
    }

    m -= 127; // unbias exponent
    ix = (ix & 0x007fffff) | 0x00800000;

    if ((m & 1) != 0) // odd m, double x to make it even
        ix += ix;

    m >>= 1; // m = [m/2]

    // generate sqrt(x) bit by bit
    ix += ix;
    var s: i32 = 0;
    var q: i32 = 0; // q = sqrt(x)
    var r: i32 = 0x01000000; // r = moving bit from right to left

    while (r != 0) {
        const t: i32 = s + r;
        if (t <= ix) {
            s = t + r;
            ix -= t;
            q += r;
        }

        ix += ix;
        r >>= 1;
    }

    // use floating add to find out rounding direction
    var z: f32 = 0;
    if (ix != 0) {
        z = 0x1p0 - 0x1.4484cp-100; // trigger inexact flag.
        if (z >= 0x1p0) {
            z = 0x1p0 + 0x1.4484cp-100;
            if (z > 0x1p0) {
                q += 2;
            } else {
                q += (q & 1);
            }
        }
    }
    ix = (q >> 1) + 0x3f000000;
    ix += (m << 23);
    return @bitCast(ix);
}

fn sqrt64(x: f64) f64 {
    const rt0: f64 = 9.99999999859990725855365213134618e-01;
    const rt1: f64 = 4.99999999495955425917856814202739e-01;
    const rt2: f64 = 3.75017500867345182581453026130850e-01;
    const rt3: f64 = 3.12523626554518656309172508769531e-01;
    const big: f64 = 134217728.0;

    var a: [2]i32 = @bitCast(x);
    const k: i32 = a[root.HIGH_HALF];
    a[root.HIGH_HALF] = (k & 0x001fffff) | 0x3fe00000;
    var t: f64 = root.inroot[@intCast((k & 0x001fffff) >> 14)];
    const s: f64 = @bitCast(a);
    //----------------- 2^-1022  <= | x |< 2^1024  -----------------
    var c: [2]i32 = .{ 0, 0 };
    if (k > 0x000fffff and k < 0x7ff00000) {
        var y: f64 = 1.0 - t * (t * s);
        t = t * (rt0 + y * (rt1 + y * (rt2 + y * rt3)));
        c[root.HIGH_HALF] = 0x20000000 + ((k & 0x7fe00000) >> 1);
        y = t * s;
        const hy: f64 = (y + big) - big;
        const del: f64 = 0.5 * t * ((s - hy * hy) - (y - hy) * (y + hy));
        var res: f64 = y + del;
        var ret: f64 = undefined;
        if (res == (res + 1.002 * ((y - res) + del))) {
            ret = res * @as(f64, @bitCast(c));
        } else {
            const res1: f64 = res + 1.5 * ((y - res) + del);
            var z: f64 = undefined;
            var zz: f64 = undefined;
            dla.emulv(res, res1, &z, &zz); // (z+zz)=res*res1
            res = if (((z - s) + zz) < 0) @max(res, res1) else @min(res, res1);
            ret = res * @as(f64, @bitCast(c));
        }
        std.mem.doNotOptimizeAway(ret);
        const dret: f64 = x / ret;
        if (dret != ret) {
            const force_inexact: f64 = 1.0 / 3.0;
            std.mem.doNotOptimizeAway(force_inexact);
            // The square root is inexact, ret is the round-to-nearest
            // value which may need adjusting for other rounding
            // modes.
        }
        // Otherwise (x / ret == ret), either the square root was exact or
        // the division was inexact.
        return ret;
    } else {
        if ((k & 0x7ff00000) == 0x7ff00000)
            return x * x + x; // sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN

        if (x == 0)
            return x; // sqrt(+0)=+0, sqrt(-0)=-0

        if (k < 0)
            return (x - x) / (x - x); // sqrt(-ve)=sNaN

        return 0x1p-256 * sqrt64(x * 0x1p512);
    }
}

// Very close to perfect but fails too many tests, although by very little margins
fn sqrt128(x: f128) f128 {
    const t16382: f128 = 0x1p16382;
    const tm8191: f128 = 0x1p-8191;
    const big: f128 = 0x1p120;
    const big1: f128 = 0x1p120 + 1;

    const a: [2]u64 = @bitCast(x);
    const hi_index: u32 = if (@import("builtin").cpu.arch.endian() == .big) 0 else 1;
    const lo_index: u32 = 1 - hi_index;
    // Extract the components of the binary128 number
    const sign_exp_hi: u64 = a[hi_index];
    // Get the exponent bits (15 bits)
    const exp: i64 = @intCast((sign_exp_hi >> 48) & 0x7fff);
    // Mask for the high part mantissa bits
    const mantissa_hi: u64 = sign_exp_hi & 0x0000ffffffffffff;
    // Low part of mantissa
    const mantissa_lo: u64 = a[lo_index];
    //----------------- 2^-16382 <= |x| < 2^16384 -----------------
    if (exp > 0 and exp < 0x7fff) {
        // Normal number
        if (x < 0) return (big1 - big1) / (big - big); // Return NaN for negative input
        // Construct a normalized value with exponent 0x3ffe (bias 0x3fff + 0)
        const adjusted_exp: i64 = exp - 0x3fff;
        const is_odd: bool = (adjusted_exp & 1) != 0;
        const exponent_part: u64 = if (is_odd) 0x4000 else 0x3fff;
        var norm: [2]u64 = undefined;
        norm[hi_index] = (sign_exp_hi & 0x8000000000000000) | (exponent_part << 48) | mantissa_hi;
        norm[lo_index] = mantissa_lo;
        const s: f128 = @bitCast(norm);
        // Compute initial approximation using double precision sqrt
        const d_approx: f64 = sqrt(cast(f64, s, .{}));
        var i: f128 = cast(f128, d_approx, .{});
        // Set the exponent of the result: for sqrt, divide exponent by 2
        const new_exp: i64 = if (is_odd) ((adjusted_exp - 1) >> 1) + 0x3fff else (adjusted_exp >> 1) + 0x3fff;
        var c: [2]u64 = undefined;
        c[hi_index] = (a[hi_index] & 0x8000000000000000) | (@as(u64, @intCast(new_exp & 0x7fff)) << 48);
        c[lo_index] = 0;
        // Newton-Raphson iterations for binary128 precision
        const t: f128 = 0.5 * (i + s / i);
        i = 0.5 * (t + s / t);
        i = 0.5 * (i + s / i); // Extra iteration for quad precision
        return @as(f128, @bitCast(c)) * i;
    } else {
        // Handle special cases
        if (exp == 0x7fff) {
            // sqrt(-Inf) = NaN, sqrt(NaN) = NaN, sqrt(+Inf) = +Inf
            return x * x + x;
        }
        if (x == 0) return x; // sqrt(+0) = +0, sqrt(-0) = -0
        if (x < 0) return (big1 - big1) / (big - big); // Return NaN for negative input
        // Subnormal number: scale up and try again
        return tm8191 * sqrt(x * t16382);
    }
}
