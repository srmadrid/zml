test {
    const test_basic = false;
    const test_exponential = true;
    const test_power = true;
    const test_trigonometric = true;
    const test_hyperbolic = true;
    const test_error_gamma = true;
    const test_bessel = false;
    const test_nearest_integer = false;
    const test_floating_point = false;

    if (test_basic) {
        //_ = @import("float/fmod.zig");
        //_ = @import("float/remainder.zig");
        //_ = @import("float/remquo.zig");
        //_ = @import("float/fdim.zig");
    }

    if (test_exponential) {
        _ = @import("float/exp.zig");
        _ = @import("float/exp10.zig");
        _ = @import("float/exp2.zig");
        _ = @import("float/exp10m1.zig");
        _ = @import("float/exp2m1.zig");
        _ = @import("float/expm1.zig");
        _ = @import("float/log.zig");
        _ = @import("float/log10.zig");
        _ = @import("float/log2.zig");
        _ = @import("float/log10p1.zig");
        _ = @import("float/log2p1.zig");
        _ = @import("float/log1p.zig");
    }

    if (test_power) {
        _ = @import("float/pow.zig");
        _ = @import("float/sqrt.zig");
        _ = @import("float/cbrt.zig");
        _ = @import("float/hypot.zig");
    }

    if (test_trigonometric) {
        _ = @import("float/sin.zig");
        _ = @import("float/cos.zig");
        _ = @import("float/tan.zig");
        _ = @import("float/asin.zig");
        _ = @import("float/acos.zig");
        _ = @import("float/atan.zig");
        _ = @import("float/atan2.zig");
        _ = @import("float/sincos.zig");
        _ = @import("float/sinpi.zig");
        _ = @import("float/cospi.zig");
        _ = @import("float/tanpi.zig");
        _ = @import("float/asinpi.zig");
        _ = @import("float/acospi.zig");
        _ = @import("float/atanpi.zig");
        _ = @import("float/atan2pi.zig");
    }

    if (test_hyperbolic) {
        _ = @import("float/sinh.zig");
        _ = @import("float/cosh.zig");
        _ = @import("float/tanh.zig");
        _ = @import("float/asinh.zig");
        _ = @import("float/acosh.zig");
        _ = @import("float/atanh.zig");
    }

    if (test_error_gamma) {
        _ = @import("float/erf.zig");
        _ = @import("float/erfc.zig");
        _ = @import("float/gamma.zig");
        _ = @import("float/lgamma.zig");
    }

    if (test_bessel) {
        //_ = @import("float/j0.zig");
        //_ = @import("float/j1.zig");
        //_ = @import("float/jn.zig");
        //_ = @import("float/y0.zig");
        //_ = @import("float/y1.zig");
        //_ = @import("float/yn.zig");
    }

    if (test_nearest_integer) {
        //_ = @import("float/nearbyint.zig");
        //_ = @import("float/lrint.zig");
        //_ = @import("float/llrint.zig");
    }

    if (test_floating_point) {
        //_ = @import("float/ilogb.zig");
        //_ = @import("float/logb.zig");
        //_ = @import("float/nextafter.zig");
        //_ = @import("float/nexttoward.zig");
    }
}
