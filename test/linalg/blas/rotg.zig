const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const rotg = zml.linalg.blas.rotg;

test rotg {
    var a: f64 = 2;
    var b: f64 = 2;
    var c: f64 = undefined;
    var s: f64 = undefined;

    rotg(f64, &a, &b, &c, &s);

    try std.testing.expectApproxEqAbs(2.8284271247461903, a, 0.0000001);
    try std.testing.expectApproxEqAbs(1.4142135623730951, b, 0.0000001);
    try std.testing.expectApproxEqAbs(0.7071067811865475, c, 0.0000001);
    try std.testing.expectApproxEqAbs(0.7071067811865475, s, 0.0000001);

    var a_c: cf64 = cf64.init(1, 2);
    var b_c: cf64 = cf64.init(3, 4);
    var c_c: f64 = undefined;
    var s_c: cf64 = undefined;

    rotg(cf64, &a_c, &b_c, &c_c, &s_c);

    //try std.testing.expectApproxEqAbs(2.449489742783178, a_c.re, 0.0000001);
    //try std.testing.expectApproxEqAbs(-5.375387381226731e102, a_c.im, 0.0000001);
    //try std.testing.expectApproxEqAbs(3, b_c.re, 0.0000001);
    //try std.testing.expectApproxEqAbs(2, b_c.im, 0.0000001);
    //try std.testing.expectApproxEqAbs(0.6172133998483676, c_c, 0.0000001);
    //try std.testing.expectApproxEqAbs(0.7715167498104596, s_c.re, 0.0000001);
    //try std.testing.expectApproxEqAbs(0.15430334996209188, s_c.im, 0.0000001);
}
