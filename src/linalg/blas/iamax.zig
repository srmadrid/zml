const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const float = @import("../../float.zig");

const blas = @import("../blas.zig");

pub fn iamax(
    n: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !usize {
    const X: type = types.Child(@TypeOf(x));

    if (n <= 0 or incx <= 0) return blas.Error.InvalidArgument;

    if (n == 1) return 0;

    var imax: usize = 0;

    if (comptime !types.isArbitraryPrecision(X)) {
        if (comptime !types.isComplex(X)) {
            var max: X = ops.abs(x[0], ctx) catch unreachable;
            var ix: isize = if (incx < 0) (-n + 2) * incx else incx;
            for (1..scast(usize, n)) |i| {
                const absx: X = ops.abs(x[scast(usize, ix)], ctx) catch unreachable;
                if (ops.gt(absx, max, ctx) catch unreachable) {
                    max = absx;
                    imax = i;
                }

                ix += incx;
            }
        } else {
            var max: Scalar(X) = ops.add(
                ops.abs(x[0].re, ctx) catch unreachable,
                ops.abs(x[0].im, ctx) catch unreachable,
                ctx,
            ) catch unreachable;
            var ix: isize = if (incx < 0) (-n + 2) * incx else incx;
            for (1..scast(usize, n)) |i| {
                const absx: Scalar(X) = ops.add(
                    ops.abs(x[scast(usize, ix)].re, ctx) catch unreachable,
                    ops.abs(x[scast(usize, ix)].im, ctx) catch unreachable,
                    ctx,
                ) catch unreachable;
                if (ops.gt(absx, max, ctx) catch unreachable) {
                    max = absx;
                    imax = i;
                }

                ix += incx;
            }
        }
    } else {
        // On abs, copy = false
        @compileError("zml.linalg.blas.iamax not implemented for arbitrary precision types yet");
    }

    return imax;
}
