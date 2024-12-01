const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

const scalar = core.supported.scalar;

pub inline fn rotmg(comptime T: type, d1: *T, d2: *T, x1: *T, y1: T, param: [*]T) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.rotmg does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            const gam: T = 4096;
            const gamsq: T = 16777216;
            const rgamsq: T = 5.9604645e-8;

            var flag: T = 0;
            var h11: T = 0;
            var h12: T = 0;
            var h21: T = 0;
            var h22: T = 0;
            var p1: T = 0;
            var p2: T = 0;
            var q1: T = 0;
            var q2: T = 0;
            var temp: T = 0;
            var u: T = 0;

            if (d1.* < 0) {
                flag = -1;
                h11 = 0;
                h12 = 0;
                h21 = 0;
                h22 = 0;

                d1.* = 0;
                d2.* = 0;
                x1.* = 0;
            } else {
                p2 = d2.* * y1;
                if (p2 == 0) {
                    flag = -2;
                    param[0] = flag;
                    return;
                }

                p1 = d1.* * x1.*;
                q2 = p2 * y1;
                q1 = p1 * x1.*;

                if (@abs(q1) > @abs(q2)) {
                    h21 = -y1 / x1.*;
                    h12 = p2 / p1;

                    u = 1 - h12 * h21;

                    if (u > 0) {
                        flag = 0;
                        d1.* /= u;
                        d2.* /= u;
                        x1.* *= u;
                    } else {
                        flag = -1;
                        h11 = 0;
                        h12 = 0;
                        h21 = 0;
                        h22 = 0;

                        d1.* = 0;
                        d2.* = 0;
                        x1.* = 0;
                    }
                } else {
                    if (q2 < 0) {
                        flag = -1;
                        h11 = 0;
                        h12 = 0;
                        h21 = 0;
                        h22 = 0;

                        d1.* = 0;
                        d2.* = 0;
                        x1.* = 0;
                    } else {
                        flag = 1;
                        h11 = p1 / p2;
                        h22 = x1.* / y1;
                        u = 1 + h11 * h22;
                        temp = d2.* / u;
                        d2.* = d1.* / u;
                        d1.* = temp;
                        x1.* = y1 * u;
                    }
                }

                if (d1.* != 0) {
                    while ((d1.* <= rgamsq) or (d1.* >= gamsq)) {
                        if (flag == 0) {
                            h11 = 1;
                            h22 = 1;
                            flag = -1;
                        } else {
                            h21 = -1;
                            h12 = 1;
                            flag = -1;
                        }
                        if (d1.* <= rgamsq) {
                            d1.* *= gam * gam;
                            x1.* /= gam;
                            h11 /= gam;
                            h12 /= gam;
                        } else {
                            d1.* /= gam * gam;
                            x1.* *= gam;
                            h11 *= gam;
                            h12 *= gam;
                        }
                    }
                }

                if (d2.* != 0) {
                    while ((@abs(d2.*) <= rgamsq) or (@abs(d2.*) >= gamsq)) {
                        if (flag == 0) {
                            h11 = 1;
                            h22 = 1;
                            flag = -1;
                        } else {
                            h21 = -1;
                            h12 = 1;
                            flag = -1;
                        }
                        if (@abs(d2.*) <= rgamsq) {
                            d2.* *= gam * gam;
                            h21 /= gam;
                            h22 /= gam;
                        } else {
                            d2.* /= gam * gam;
                            h21 *= gam;
                            h22 *= gam;
                        }
                    }
                }
            }

            if (flag < 0) {
                param[1] = h11;
                param[2] = h21;
                param[3] = h12;
                param[4] = h22;
            } else if (flag == 0) {
                param[2] = h21;
                param[3] = h12;
            } else {
                param[1] = h11;
                param[4] = h22;
            }

            param[0] = flag;
        },
        .Complex => @compileError("BLAS.rotmg does not support complex numbers."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.rotmg only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "rotmg" {
    var d1: f64 = 1;
    var d2: f64 = 2;
    var x1: f64 = 1;
    const y1: f64 = 1;
    var param: [5]f64 = undefined;

    @import("BLAS.zig").rotmg(f64, &d1, &d2, &x1, y1, &param);

    try std.testing.expectApproxEqAbs(1.333333, d1, 0.0001);
    try std.testing.expectApproxEqAbs(0.666667, d2, 0.0001);
    try std.testing.expectApproxEqAbs(1.5, x1, 0.0001);
    try std.testing.expectEqual(1, param[0]);
    try std.testing.expectApproxEqAbs(0.5, param[1], 0.0001);
    try std.testing.expectApproxEqAbs(0, param[2], 0.0001);
    try std.testing.expectApproxEqAbs(0, param[3], 0.0001);
    try std.testing.expectApproxEqAbs(1, param[4], 0.0001);
}
