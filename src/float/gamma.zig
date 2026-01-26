const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");

const lgamma = @import("lgamma.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub fn Gamma(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.gamma: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

/// Returns the gamma function of a float, int or bool  operand. The result type
/// is determined by coercing the operand type to a  float, and the operation is
/// performed by casting the operand to the result type, then computing the
/// gamma function at that value.
///
/// The gamma function is defined as:
/// $$
/// \Gamma(x) = \int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t
/// $$
///
/// ## Signature
/// ```zig
/// float.gamma(x: X) float.Gamma(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The value to get the gamma function at.
///
/// ## Returns
/// `float.Gamma(@TypeOf(x))`: The gamma function at `x`.
pub inline fn gamma(x: anytype) float.Gamma(@TypeOf(x)) {
    switch (float.Gamma(@TypeOf(x))) {
        f16 => return types.scast(f16, gamma32(types.scast(f32, x))),
        f32 => {
            // https://github.com/JuliaMath/openlibm/blob/master/src/s_tgammaf.c
            return gamma32(types.scast(f32, x));
        },
        f64 => {
            //
            // return gamma64(types.scast(f64, x));
            return types.scast(f64, gamma128(types.scast(f128, x)));
        },
        f80 => {
            //
            // return gamma80(types.scast(f80, x));
            return types.scast(f80, gamma128(types.scast(f128, x)));
        },
        f128 => {
            // https://github.com/JuliaMath/openlibm/blob/master/ld128/e_tgammal.c
            return gamma128(types.scast(f128, x));
        },
        else => unreachable,
    }
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/src/s_tgammaf.c
//
// Original copyright notice:
// Copyright (c) 2008 David Schultz <das@FreeBSD.ORG>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
fn gamma32(x: f32) f32 {
    return types.scast(f32, gamma64(types.scast(f64, x)));
}

fn gamma64(x: f64) f64 {
    return types.scast(f64, gamma128(types.scast(f128, x)));
}

// Translation of:
// https://github.com/JuliaMath/openlibm/blob/master/ld128/e_tgammal.c
//
// Original copyright notice:
// Copyright (c) 2011 Martynas Venckus <martynas@openbsd.org>
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
fn gamma128(x: f128) f128 {
    const I0: i64 = @bitCast(ldbl128.getHighPart(x));
    const I1: i64 = @bitCast(ldbl128.getLowPart(x));
    if (((I0 & 0x7fffffffffffffff) | I1) == 0)
        return 1.0 / x;

    if (I0 < 0 and @as(u64, @bitCast(I0)) < 0xffff000000000000 and float.rint(x) == x)
        return (x - x) / (x - x);

    if (I0 == 0xffff000000000000 and I1 == 0)
        return x - x;
    var signamp: i32 = undefined;
    const lg: f128 = lgamma.lgamma_r128(x, &signamp);
    return types.scast(f128, signamp) * float.exp(lg);
}
