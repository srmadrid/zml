const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

pub inline fn rotmg(comptime T: type, d1: *T, d2: *T, x1: *T, y1: T, param: [*]T) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    switch (numericType) {
        .bool => @compileError("blas.rotmg does not support bool."),
        .int, .float => {
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
        .cfloat => @compileError("blas.rotmg does not support complex numbers."),
        .integer, .rational, .real, .complex, .expression => @compileError("blas.rotmg only supports simple types."),
        .unsupported => unreachable,
    }
}
