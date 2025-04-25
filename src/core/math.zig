const std = @import("std");
// Make in-place variants

// Edit all functions to use unified bit extraction

// Add @branchHint to ifs

// Math implementations will be based on glibc and openlibm, choosing the best
// based on:
// Brian Gladman, Vincenzo Innocente, John Mather, Paul Zimmermann. Accuracy of
// Mathematical Functions in Single, Double, Double Extended, and Quadruple
// Precision. 2025. hal-03141101v8
//
// Glibc implementations are based on:
// f32:  glibc/sysdeps/ieee754/flt-32/*.c
// f64:  glibc/sysdeps/ieee754/dbl-64/*.c
// f80:  glibc/sysdeps/ieee754/ldbl-96/*.c
// f128: glibc/sysdeps/ieee754/ldbl-128/*.c
//
// Openlibm implementations are based on:
// f32:  openlibm/src/*f.c
// f64:  openlibm/src/*.c
// f80:  openlibm/ld80/*.c
// f128: openlibm/src/*l.c, or openlibm/ld128/*.c
//
// Each function will state the source of the implementation, and the tests will
// be based on the tests from glibc using tonearest rouding values.

pub const sin = @import("math/sin.zig").sin; // to implement: f80. 20/609 tests fail: 1 for f32, 6 for f64, 13 for f128
pub const cos = @import("math/cos.zig").cos; // to implement: f80. 16/532 tests fail: 2 for f32, 2 for f64, 12 for f128
pub const tan = @import("math/tan.zig").tan; // to implement: f80. 3/546 tests fail: 3 for f128
pub const asin = @import("math/asin.zig").asin; // to implement: f80. 3/385 tests fail: 1 for f64, 2 for f128
pub const acos = @import("math/acos.zig").acos; // to implement: f80. 2/494 tests fail: 1 for f64, 1 for f128
pub const atan = @import("math/atan.zig").atan; // to implement: f80. 3/228 tests fail: 1 for f64, 2 for f128
pub const atan2 = @import("math/atan2.zig").atan2; // to implement: f80. 107/2672 tests fail: 2 for f80, 105 for f128
pub const sincos = @import("math/sincos.zig").sincos; // to implement: f80. 8/674 tests fail: 1 for f64, 7 for f128

pub const sinh = @import("math/sinh.zig").sinh; // to implement: f80. 67/577 tests fail: 29 for f64, 38 for f128
pub const cosh = @import("math/cosh.zig").cosh; // to implement: f80. 52/577 tests fail: 12 for f64, 40 for f128
pub const tanh = @import("math/tanh.zig").tanh; // to implement: f80. 40/452 tests fail: 15 for f64, 25 for f128
pub const asinh = @import("math/asinh.zig").asinh; // to implement: f80. 51/507 test fail: 16 fon f64, 35 for f128
// pub const acosh = @import("math/acosh.zig").acosh; // to implement
// pub const atanh = @import("math/atanh.zig").atanh; // to implement

// pub const sinpi = @import("math/sinpi.zig").sinpi; // to implement
// pub const cospi = @import("math/cospi.zig").cospi; // to implement
// pub const tanpi = @import("math/tanpi.zig").tanpi; // to implement
// pub const asinpi = @import("math/asinpi.zig").asinpi; // to implement
// pub const acospi = @import("math/acospi.zig").acospi; // to implement
// pub const atanpi = @import("math/atanpi.zig").atanpi; // to implement
// pub const atan2pi = @import("math/atan2pi.zig").atan2pi; // to implement

pub const exp = @import("math/exp.zig").exp; // to implement: f80. 5/768 tests fail: 1 for f32, 1 for f64, 3 for f128
// pub const exp2 = @import("math/exp2.zig").exp2; // to implement
pub const expm1 = @import("math/expm1.zig").expm1; // to implement: f80. 28/507 tests fail: 12 for f64, 16 for f128
pub const log = @import("math/log.zig").log; // to implement: f80. 8/256 tests fail: 8 for f128
// pub const log2 = @import("math/log2.zig").log2; // to implement
// pub const log10 = @import("math/log10.zig").log10; // to implement
pub const log1p = @import("math/log1p.zig").log1p; // to implement: f80. 30/421 tests fail: 8 for f64, 22 for f128
// pub const logb = @import("math/logb.zig").logb; // to implement

// pub const pow = @import("math/pow.zig").pow; // to implement
pub const sqrt = @import("math/sqrt.zig").sqrt; // 58/655 tests fail: 1 for f80, 57 for f128
// pub const cbrt = @import("math/cbrt.zig").cbrt; // to implement
pub const hypot = @import("math/hypot.zig").hypot; // to implement: f80. 2/2281 tests fail: 1 for f32, 1 for f128
// pub const fmod = @import("math/fmod.zig").fmod; // to implement
// pub const remainder = @import("math/remainder.zig").remainder; // to implement
// pub const fma = @import("math/fma.zig").fma; // to implement

// pub const erf = @import("math/erf.zig").erf; // to implement
// pub const erfc = @import("math/erfc.zig").erfc; // to implement
// pub const gamma = @import("math/gamma.zig").gamma; // to implement
// pub const lgamma = @import("math/lgamma.zig").lgamma; // to implement
// pub const tgamma = @import("math/tgamma.zig").tgamma; // to implement
// pub const j0 = @import("math/j0.zig").j0; // to implement
// pub const j1 = @import("math/j1.zig").j1; // to implement
// pub const jn = @import("math/jn.zig").jn; // to implement
// pub const y0 = @import("math/y0.zig").y0; // to implement
// pub const y1 = @import("math/y1.zig").y1; // to implement
// pub const yn = @import("math/yn.zig").yn; // to implement

// pub const ceil = @import("math/ceil.zig").ceil; // to implement
pub const floor = @import("math/floor.zig").floor;
// pub const trunc = @import("math/trunc.zig").trunc; // to implement
pub const round = @import("math/round.zig").round; // to implement
// pub const nearbyint = @import("math/nearbyint.zig").nearbyint; // to implement
// pub const rint = @import("math/rint.zig").rint; // to implement
// pub const lrint = @import("math/lrint.zig").lrint; // to implement
// pub const llrint = @import("math/llrint.zig").llrint; // to implement

pub const abs = @import("math/abs.zig").abs;
pub const frexp = @import("math/frexp.zig").frexp;
pub const ldexp = @import("math/ldexp.zig").ldexp;
// pub const modf = @import("math/modf.zig").modf; // to implement
pub const scalbn = @import("math/scalbn.zig").scalbn;
pub const copysign = @import("math/copysign.zig").copysign;
// pub const nextafter = @import("math/nextafter.zig").nextafter; // to implement

test {
    std.testing.refAllDeclsRecursive(@This());
}
