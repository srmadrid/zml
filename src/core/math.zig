const std = @import("std");

// Basic functions
pub const abs = @import("math/abs.zig").abs;
// pub const fmod = @import("math/fmod.zig").fmod; // to implement
// pub const remainder = @import("math/remainder.zig").remainder; // to implement
// pub const remquo = @import("math/remquo.zig").remquo; // to implement
// pub const fdim = @import("math/fdim.zig").fdim; // to implement

// Exponential functions
pub const exp = @import("math/exp.zig").exp; // to implement: f80. 5/768 tests fail: 1 for f32, 1 for f64, 3 for f128
pub const exp10 = @import("math/exp10.zig").exp10; // to implement: f80. 54/650 tests fail: 54 for f128
pub const exp2 = @import("math/exp2.zig").exp2; // to implement: f80. 5/670 tests fail: 1 for f64, 4 for f128
pub const exp10m1 = @import("math/exp10m1.zig").exp10m1; // 235/977 tests fail: 28 for f32, 68 for f64, 55 for f80, 84 for f128
pub const exp2m1 = @import("math/exp2m1.zig").exp2m1; // 189/957 tests fail: 24 for f32, 49 for f64, 48 for f80, 68 for f128
pub const expm1 = @import("math/expm1.zig").expm1; // to implement: f80. 28/507 tests fail: 12 for f64, 16 for f128
pub const log = @import("math/log.zig").log; // to implement: f80. 8/256 tests fail: 8 for f128
pub const log10 = @import("math/log10.zig").log10; // to implement: f80. 33/264 tests fail: 16 for f64, 17 for f128
// pub const log2 = @import("math/log2.zig").log2; // to implement
// pub const log10p1 = @import("math/log10p1.zig").log10p1; // to implement
// pub const log2p1 = @import("math/log2p1.zig").log2p1; // to implement
pub const log1p = @import("math/log1p.zig").log1p; // to implement: f80. 30/421 tests fail: 8 for f64, 22 for f128

// Power functions
// pub const pow = @import("math/pow.zig").pow; // to implement
pub const sqrt = @import("math/sqrt.zig").sqrt; // 58/655 tests fail: 1 for f80, 57 for f128
// pub const cbrt = @import("math/cbrt.zig").cbrt; // to implement
pub const hypot = @import("math/hypot.zig").hypot; // to implement: f80. 2/2281 tests fail: 1 for f32, 1 for f128

// Trigonometric functions
pub const sin = @import("math/sin.zig").sin; // to implement: f80. 20/609 tests fail: 1 for f32, 6 for f64, 13 for f128
pub const cos = @import("math/cos.zig").cos; // to implement: f80. 16/532 tests fail: 2 for f32, 2 for f64, 12 for f128
pub const tan = @import("math/tan.zig").tan; // to implement: f80. 3/546 tests fail: 3 for f128
pub const asin = @import("math/asin.zig").asin; // to implement: f80. 3/385 tests fail: 1 for f64, 2 for f128
pub const acos = @import("math/acos.zig").acos; // to implement: f80. 2/494 tests fail: 1 for f64, 1 for f128
pub const atan = @import("math/atan.zig").atan; // to implement: f80. 3/228 tests fail: 1 for f64, 2 for f128
pub const atan2 = @import("math/atan2.zig").atan2; // to implement: f80. 107/2672 tests fail: 2 for f80, 105 for f128
pub const sincos = @import("math/sincos.zig").sincos; // to implement: f80. 8/674 tests fail: 1 for f64, 7 for f128
pub const sinpi = @import("math/sinpi.zig").sinpi; // 154/1572 tests fail: 28 for f32, 25 for f64, 38 for f80, 63 for f128
pub const cospi = @import("math/cospi.zig").cospi; // 266/1468 tests fail: 43 for f32, 39 for f64, 70 for f80, 114 for f128
pub const tanpi = @import("math/tanpi.zig").tanpi; // 535/1436 tests fail: 58 for f32, 131 for f64, 128 for f80, 218 for f128
pub const asinpi = @import("math/asinpi.zig").asinpi; // 121/437 tests fail: 25 for f32, 20 f64, 30 for f80, 46 for f128
pub const acospi = @import("math/acospi.zig").acospi; // 131/546 tests fail: 17 for f32, 29 for f64, 29 for f80, 56 for f128
pub const atanpi = @import("math/atanpi.zig").atanpi; // 51/272 tests fail: 8 for f32, 11 for f64, 9 for f80, 23 for f128
pub const atan2pi = @import("math/atan2pi.zig").atan2pi; // 177/2413: 7 for f32, 30 for f64, 45 for f80, 95 for f128

// Hyperbolic functions
pub const sinh = @import("math/sinh.zig").sinh; // to implement: f80. 67/577 tests fail: 29 for f64, 38 for f128
pub const cosh = @import("math/cosh.zig").cosh; // to implement: f80. 52/577 tests fail: 12 for f64, 40 for f128
pub const tanh = @import("math/tanh.zig").tanh; // to implement: f80. 40/452 tests fail: 15 for f64, 25 for f128
pub const asinh = @import("math/asinh.zig").asinh; // to implement: f80. 51/507 tests fail: 16 for f64, 35 for f128
pub const acosh = @import("math/acosh.zig").acosh; // to implement: f80. 51/349 tests fail: 23 for f64, 28 for f128
pub const atanh = @import("math/atanh.zig").atanh; // to implement: f80. 48/566 tests fail: 17 for f64, 31 for f128

// Error and gamma functions
// pub const erf = @import("math/erf.zig").erf; // to implement
// pub const erfc = @import("math/erfc.zig").erfc; // to implement
// pub const gamma = @import("math/gamma.zig").gamma; // to implement
// pub const lgamma = @import("math/lgamma.zig").lgamma; // to implement
// pub const tgamma = @import("math/tgamma.zig").tgamma; // to implement

// Bessel functions
// pub const j0 = @import("math/j0.zig").j0; // to implement
// pub const j1 = @import("math/j1.zig").j1; // to implement
// pub const jn = @import("math/jn.zig").jn; // to implement
// pub const y0 = @import("math/y0.zig").y0; // to implement
// pub const y1 = @import("math/y1.zig").y1; // to implement
// pub const yn = @import("math/yn.zig").yn; // to implement

// Nearest integer floating-point operations
// pub const ceil = @import("math/ceil.zig").ceil; // to implement
pub const floor = @import("math/floor.zig").floor;
// pub const trunc = @import("math/trunc.zig").trunc; // to implement
pub const round = @import("math/round.zig").round; // to implement
// pub const nearbyint = @import("math/nearbyint.zig").nearbyint; // to implement
// pub const rint = @import("math/rint.zig").rint; // to implement
// pub const lrint = @import("math/lrint.zig").lrint; // to implement
// pub const llrint = @import("math/llrint.zig").llrint; // to implement. Join the rint's and choose one depending on input integer type?

// Floating-point manipulation functions
pub const frexp = @import("math/frexp.zig").frexp;
pub const ldexp = @import("math/ldexp.zig").ldexp;
// pub const modf = @import("math/modf.zig").modf; // to implement
pub const scalbn = @import("math/scalbn.zig").scalbn;
// pub const ilogb = @import("math/ilogb.zig").ilogb; // to implement
// pub const logb = @import("math/logb.zig").logb; // to implement
// pub const nextafter = @import("math/nextafter.zig").nextafter; // to implement
// pub const nexttoward = @import("math/nexttoward.zig").nexttoward; // to implement
pub const copysign = @import("math/copysign.zig").copysign;

test {
    const test_basic = false;
    const test_exponential = true;
    const test_power = false;
    const test_trigonometric = false;
    const test_hyperbolic = false;
    const test_error_gamma = false;
    const test_bessel = false;
    const test_nearest_integer = false;
    const test_floating_point = false;

    if (test_basic) {
        _ = @import("math/abs.zig");
        //_ = @import("math/fmod.zig");
        _ = @import("math/remainder.zig");
        //_ = @import("math/remquo.zig");
        //_ = @import("math/fdim.zig");
    }

    if (test_exponential) {
        _ = @import("math/exp.zig");
        _ = @import("math/exp10.zig");
        _ = @import("math/exp2.zig");
        _ = @import("math/exp10m1.zig");
        _ = @import("math/exp2m1.zig");
        _ = @import("math/expm1.zig");
        _ = @import("math/log.zig");
        _ = @import("math/log10.zig");
        _ = @import("math/log2.zig");
        //_ = @import("math/log10p1.zig");
        //_ = @import("math/log2p1.zig");
        _ = @import("math/log1p.zig");
    }

    if (test_power) {
        _ = @import("math/pow.zig");
        _ = @import("math/sqrt.zig");
        _ = @import("math/cbrt.zig");
        _ = @import("math/hypot.zig");
    }

    if (test_trigonometric) {
        _ = @import("math/sin.zig");
        _ = @import("math/cos.zig");
        _ = @import("math/tan.zig");
        _ = @import("math/asin.zig");
        _ = @import("math/acos.zig");
        _ = @import("math/atan.zig");
        _ = @import("math/atan2.zig");
        _ = @import("math/sincos.zig");
        _ = @import("math/sinpi.zig");
        _ = @import("math/cospi.zig");
        _ = @import("math/tanpi.zig");
        _ = @import("math/asinpi.zig");
        _ = @import("math/acospi.zig");
        _ = @import("math/atanpi.zig");
        _ = @import("math/atan2pi.zig");
    }

    if (test_hyperbolic) {
        _ = @import("math/sinh.zig");
        _ = @import("math/cosh.zig");
        _ = @import("math/tanh.zig");
        _ = @import("math/asinh.zig");
        _ = @import("math/acosh.zig");
        _ = @import("math/atanh.zig");
    }

    if (test_error_gamma) {
        //_ = @import("math/erf.zig");
        //_ = @import("math/erfc.zig");
        _ = @import("math/gamma.zig");
        _ = @import("math/lgamma.zig");
        //_ = @import("math/tgamma.zig");
    }

    if (test_bessel) {
        //_ = @import("math/j0.zig");
        //_ = @import("math/j1.zig");
        //_ = @import("math/jn.zig");
        //_ = @import("math/y0.zig");
        //_ = @import("math/y1.zig");
        //_ = @import("math/yn.zig");
    }

    if (test_nearest_integer) {
        _ = @import("math/ceil.zig");
        _ = @import("math/floor.zig");
        _ = @import("math/trunc.zig");
        _ = @import("math/round.zig");
        //_ = @import("math/nearbyint.zig");
        //_ = @import("math/rint.zig");
        //_ = @import("math/lrint.zig");
        //_ = @import("math/llrint.zig");
    }

    if (test_floating_point) {
        _ = @import("math/frexp.zig");
        _ = @import("math/ldexp.zig");
        _ = @import("math/modf.zig");
        _ = @import("math/scalbn.zig");
        //_ = @import("math/ilogb.zig");
        //_ = @import("math/logb.zig");
        //_ = @import("math/nextafter.zig");
        //_ = @import("math/nexttoward.zig");
        _ = @import("math/copysign.zig");
    }
}
