const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const float = @import("../../float.zig");

const blas = @import("../blas.zig");

pub fn iamin(
    n: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !usize {
    const X: type = types.Child(@TypeOf(x));

    if (n <= 0 or incx <= 0) return blas.Error.InvalidArgument;

    if (n == 1) return 0;

    var imin: usize = 0;

    if (comptime !types.isArbitraryPrecision(X)) {
        if (comptime !types.isComplex(X)) {
            var min: X = ops.abs(x[0], ctx) catch unreachable;
            var ix: isize = if (incx < 0) (-n + 2) * incx else incx;
            for (1..scast(usize, n)) |i| {
                const absx: X = ops.abs(x[scast(usize, ix)], ctx) catch unreachable;
                if (ops.lt(absx, min, ctx) catch unreachable) {
                    min = absx;
                    imin = i;
                }

                ix += incx;
            }
        } else {
            var min: Scalar(X) = ops.add(
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
                if (ops.lt(absx, min, ctx) catch unreachable) {
                    min = absx;
                    imin = i;
                }

                ix += incx;
            }
        }
    } else {
        // On abs, copy = false
        @compileError("zml.linalg.blas.iamin not implemented for arbitrary precision types yet");
    }

    return imin;
}
