const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const lgamma_data = @import("lgamma_data.zig");
const erf_data = @import("erf_data.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Lgamma(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.lgamma: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

/// Returns the log-gamma function of a float, int or bool  operand. The result
/// type is determined by coercing the operand type to a  float, and the
/// operation is performed by casting the operand to the result type, then
/// computing the log-gamma function at that value.
///
/// The log-gamma function is defined as:
/// $$
/// \log(\Gamma(x)) = \log\left(\int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t\right)
/// $$
///
/// ## Signature
/// ```zig
/// float.lgamma(x: X) float.Lgamma(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the log-gamma function at.
///
/// ## Returns
/// `float.Lgamma(@TypeOf(x))`: The log-gamma function at `x`.
pub fn lgamma(x: anytype) float.Lgamma(@TypeOf(x)) {
    var local_signgam: i32 = undefined;
    switch (float.Lgamma(@TypeOf(x))) {
        f16 => return types.scast(f16, lgamma_r32(types.scast(f32, x), &local_signgam)),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_lgammaf_r.c
            // return lgamma_r32(types.scast(f32, x), &local_signgam);
            return types.scast(f32, lgamma_r128(types.scast(f128, x), &local_signgam));
        },
        f64 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/e_lgamma_r.c
            // return lgamma_r64(types.scast(f64, x), &local_signgam);
            return types.scast(f64, lgamma_r128(types.scast(f128, x), &local_signgam));
        },
        f80 => return types.scast(f80, lgamma_r128(types.scast(f128, x), &local_signgam)),
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_lgammal_r.c
            return lgamma_r128(types.scast(f128, x), &local_signgam);
        },
        else => unreachable,
    }
}

fn lgamma_r32(x: f32, signgamp: *i32) f32 {
    _ = signgamp;
    return x;
}

pub fn lgamma_r64(x: f64, signgamp: *i32) f64 {
    _ = signgamp;
    return x;
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_tgammal.c
//
// Original copyright notice:
// Copyright (c) 2008 Stephen L. Moshier <steve@moshier.net>
//
// Permission to use, copy, modify, and distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
pub fn lgamma_r128(x: f128, signgamp: *i32) f128 {
    signgamp.* = 1;

    if (!std.math.isFinite(x))
        return x * x;

    if (x == 0) {
        if (std.math.signbit(x))
            signgamp.* = -1;

        return 1.0 / float.abs(x);
    }

    if (x < 0) {
        const q: f128 = -x;
        var p: f128 = float.floor(q);
        if (p == q)
            return (1 / (p - p));

        const i: i128 = types.scast(i128, p);
        if ((i & 1) == 0)
            signgamp.* = -1
        else
            signgamp.* = 1;

        var z: f128 = q - p;
        if (z > 0.5) {
            p += 1.0;
            z = p - q;
        }

        z = q * float.sin(float.pi(f128) * z);
        if (z == 0.0)
            return -float.log(q);

        var s: i32 = undefined;
        const w: f128 = lgamma_r128(q, &s);
        z = float.log(float.pi(f128) / z) - w;
        return z;
    }

    if (x < 13.5) {
        var p: f128 = 0;
        const nx: f128 = float.floor(x + 0.5);
        const nn: i32 = types.scast(i32, nx);
        switch (nn) {
            0 => {
                // log gamma (x + 1) = log(x) + log gamma(x)
                if (x <= 0.125) {
                    p = x * erf_data.neval(x, &lgamma_data.RN1_128, 8) / erf_data.deval(x, &lgamma_data.RD1_128, 7);
                } else if (x <= 0.375) {
                    const z: f128 = x - 0.25;
                    p = z * erf_data.neval(z, &lgamma_data.RN1r25_128, 9) / erf_data.deval(z, &lgamma_data.RD1r25_128, 8);
                    p += lgamma_data.lgam1r25b_128;
                    p += lgamma_data.lgam1r25a_128;
                } else if (x <= 0.625) {
                    var z: f128 = x + (1 - lgamma_data.x0a_128);
                    z = z - lgamma_data.x0b_128;
                    p = erf_data.neval(z, &lgamma_data.RN1r5_128, 8) / erf_data.deval(z, &lgamma_data.RD1r5_128, 8);
                    p = p * z * z;
                    p = p + lgamma_data.y0b_128;
                    p = p + lgamma_data.y0a_128;
                } else if (x <= 0.875) {
                    const z: f128 = x - 0.75;
                    p = z * erf_data.neval(z, &lgamma_data.RN1r75_128, 8) / erf_data.deval(z, &lgamma_data.RD1r75_128, 8);
                    p += lgamma_data.lgam1r75b_128;
                    p += lgamma_data.lgam1r75a_128;
                } else {
                    const z: f128 = x - 1;
                    p = z * erf_data.neval(z, &lgamma_data.RN2_128, 9) / erf_data.deval(z, &lgamma_data.RD2_128, 9);
                }
                p = p - float.log(x);
            },
            1 => {
                if (x < 0.875) {
                    if (x <= 0.625) {
                        var z: f128 = x + (1 - lgamma_data.x0a_128);
                        z = z - lgamma_data.x0b_128;
                        p = erf_data.neval(z, &lgamma_data.RN1r5_128, 8) / erf_data.deval(z, &lgamma_data.RD1r5_128, 8);
                        p = p * z * z;
                        p = p + lgamma_data.y0b_128;
                        p = p + lgamma_data.y0a_128;
                    } else if (x <= 0.875) {
                        const z: f128 = x - 0.75;
                        p = z * erf_data.neval(z, &lgamma_data.RN1r75_128, 8) / erf_data.deval(z, &lgamma_data.RD1r75_128, 8);
                        p += lgamma_data.lgam1r75b_128;
                        p += lgamma_data.lgam1r75a_128;
                    } else {
                        const z: f128 = x - 1;
                        p = z * erf_data.neval(z, &lgamma_data.RN2_128, 9) / erf_data.deval(z, &lgamma_data.RD2_128, 9);
                    }
                    p = p - float.log(x);
                } else if (x < 1) {
                    const z: f128 = x - 1;
                    p = z * erf_data.neval(z, &lgamma_data.RNr9_128, 8) / erf_data.deval(z, &lgamma_data.RDr9_128, 8);
                } else if (x == 1) {
                    p = 0;
                } else if (x <= 1.125) {
                    const z: f128 = x - 1;
                    p = z * erf_data.neval(z, &lgamma_data.RN1_128, 8) / erf_data.deval(z, &lgamma_data.RD1_128, 7);
                } else if (x <= 1.375) {
                    const z: f128 = x - 1.25;
                    p = z * erf_data.neval(z, &lgamma_data.RN1r25_128, 9) / erf_data.deval(z, &lgamma_data.RD1r25_128, 8);
                    p += lgamma_data.lgam1r25b_128;
                    p += lgamma_data.lgam1r25a_128;
                } else {
                    // 1.375 <= x+x0 <= 1.625
                    var z: f128 = x - lgamma_data.x0a_128;
                    z = z - lgamma_data.x0b_128;
                    p = erf_data.neval(z, &lgamma_data.RN1r5_128, 8) / erf_data.deval(z, &lgamma_data.RD1r5_128, 8);
                    p = p * z * z;
                    p = p + lgamma_data.y0b_128;
                    p = p + lgamma_data.y0a_128;
                }
            },
            2 => {
                if (x < 1.625) {
                    var z: f128 = x - lgamma_data.x0a_128;
                    z = z - lgamma_data.x0b_128;
                    p = erf_data.neval(z, &lgamma_data.RN1r5_128, 8) / erf_data.deval(z, &lgamma_data.RD1r5_128, 8);
                    p = p * z * z;
                    p = p + lgamma_data.y0b_128;
                    p = p + lgamma_data.y0a_128;
                } else if (x < 1.875) {
                    const z: f128 = x - 1.75;
                    p = z * erf_data.neval(z, &lgamma_data.RN1r75_128, 8) / erf_data.deval(z, &lgamma_data.RD1r75_128, 8);
                    p += lgamma_data.lgam1r75b_128;
                    p += lgamma_data.lgam1r75a_128;
                } else if (x == 2) {
                    p = 0;
                } else if (x < 2.375) {
                    const z: f128 = x - 2;
                    p = z * erf_data.neval(z, &lgamma_data.RN2_128, 9) / erf_data.deval(z, &lgamma_data.RD2_128, 9);
                } else {
                    const z: f128 = x - 2.5;
                    p = z * erf_data.neval(z, &lgamma_data.RN2r5_128, 8) / erf_data.deval(z, &lgamma_data.RD2r5_128, 8);
                    p += lgamma_data.lgam2r5b_128;
                    p += lgamma_data.lgam2r5a_128;
                }
            },
            3 => {
                if (x < 2.75) {
                    const z: f128 = x - 2.5;
                    p = z * erf_data.neval(z, &lgamma_data.RN2r5_128, 8) / erf_data.deval(z, &lgamma_data.RD2r5_128, 8);
                    p += lgamma_data.lgam2r5b_128;
                    p += lgamma_data.lgam2r5a_128;
                } else {
                    const z: f128 = x - 3;
                    p = z * erf_data.neval(z, &lgamma_data.RN3_128, 9) / erf_data.deval(z, &lgamma_data.RD3_128, 9);
                    p += lgamma_data.lgam3b_128;
                    p += lgamma_data.lgam3a_128;
                }
            },
            4 => {
                const z: f128 = x - 4;
                p = z * erf_data.neval(z, &lgamma_data.RN4_128, 9) / erf_data.deval(z, &lgamma_data.RD4_128, 9);
                p += lgamma_data.lgam4b_128;
                p += lgamma_data.lgam4a_128;
            },
            5 => {
                const z: f128 = x - 5;
                p = z * erf_data.neval(z, &lgamma_data.RN5_128, 9) / erf_data.deval(z, &lgamma_data.RD5_128, 8);
                p += lgamma_data.lgam5b_128;
                p += lgamma_data.lgam5a_128;
            },
            6 => {
                const z: f128 = x - 6;
                p = z * erf_data.neval(z, &lgamma_data.RN6_128, 8) / erf_data.deval(z, &lgamma_data.RD6_128, 8);
                p += lgamma_data.lgam6b_128;
                p += lgamma_data.lgam6a_128;
            },
            7 => {
                const z: f128 = x - 7;
                p = z * erf_data.neval(z, &lgamma_data.RN7_128, 8) / erf_data.deval(z, &lgamma_data.RD7_128, 7);
                p += lgamma_data.lgam7b_128;
                p += lgamma_data.lgam7a_128;
            },
            8 => {
                const z: f128 = x - 8;
                p = z * erf_data.neval(z, &lgamma_data.RN8_128, 8) / erf_data.deval(z, &lgamma_data.RD8_128, 7);
                p += lgamma_data.lgam8b_128;
                p += lgamma_data.lgam8a_128;
            },
            9 => {
                const z: f128 = x - 9;
                p = z * erf_data.neval(z, &lgamma_data.RN9_128, 7) / erf_data.deval(z, &lgamma_data.RD9_128, 7);
                p += lgamma_data.lgam9b_128;
                p += lgamma_data.lgam9a_128;
            },
            10 => {
                const z: f128 = x - 10;
                p = z * erf_data.neval(z, &lgamma_data.RN10_128, 7) / erf_data.deval(z, &lgamma_data.RD10_128, 7);
                p += lgamma_data.lgam10b_128;
                p += lgamma_data.lgam10a_128;
            },
            11 => {
                const z: f128 = x - 11;
                p = z * erf_data.neval(z, &lgamma_data.RN11_128, 7) / erf_data.deval(z, &lgamma_data.RD11_128, 6);
                p += lgamma_data.lgam11b_128;
                p += lgamma_data.lgam11a_128;
            },
            12 => {
                const z: f128 = x - 12;
                p = z * erf_data.neval(z, &lgamma_data.RN12_128, 7) / erf_data.deval(z, &lgamma_data.RD12_128, 6);
                p += lgamma_data.lgam12b_128;
                p += lgamma_data.lgam12a_128;
            },
            13 => {
                const z: f128 = x - 13;
                p = z * erf_data.neval(z, &lgamma_data.RN13_128, 7) / erf_data.deval(z, &lgamma_data.RD13_128, 6);
                p += lgamma_data.lgam13b_128;
                p += lgamma_data.lgam13a_128;
            },
            else => unreachable,
        }
        return p;
    }

    if (x > 1.0485738685148938358098967157129705071571e4928)
        return types.scast(f128, signgamp.*) * std.math.inf(f128);

    if (x > 0x1p112)
        return x * (float.log(x) - 1);

    var q: f128 = lgamma_data.ls2pi_128 - x;
    q = (x - 0.5) * float.log(x) + q;
    if (x > 1.0e18)
        return q;

    const p: f128 = 1 / (x * x);
    q += erf_data.neval(p, &lgamma_data.RASY_128, 12) / x;
    return q;
}
