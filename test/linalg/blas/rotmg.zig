const std = @import("std");
const zml = @import("zml");
const rotmg = zml.linalg.blas.rotmg;

test rotmg {
    var d1: f64 = 1;
    var d2: f64 = 2;
    var x1: f64 = 1;
    const y1: f64 = 1;
    var param: [5]f64 = undefined;
    const param_ptr: []f64 = &param;

    rotmg(&d1, &d2, &x1, y1, param_ptr.ptr, .{}) catch unreachable;

    try std.testing.expectApproxEqAbs(1.3333333333333333, d1, 0.0000000001);
    try std.testing.expectApproxEqAbs(0.6666666666666666, d2, 0.0000000001);
    try std.testing.expectApproxEqAbs(1.5, x1, 0.0000000001);
    try std.testing.expectEqual(1, param[0]);
    try std.testing.expectApproxEqAbs(0.5, param[1], 0.0000000001);
    try std.testing.expectApproxEqAbs(0, param[2], 0.0000000001);
    try std.testing.expectApproxEqAbs(0, param[3], 0.0000000001);
    try std.testing.expectApproxEqAbs(1, param[4], 0.0000000001);
}
