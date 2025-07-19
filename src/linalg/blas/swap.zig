const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn swap(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = types.Coerce(X, Y);

    if (n <= 0) return blas.Error.InvalidArgument;

    if (comptime types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.scal not implemented for arbitrary precision types yet");
    } else {
        var temp: C = try ops.init(C, .{});
        var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
        var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
        for (0..scast(usize, n)) |_| {
            ops.set( // temp = x[ix]
                &temp,
                x[scast(usize, ix)],
                ctx,
            ) catch unreachable;

            ops.set( // x[ix] = y[iy]
                &x[scast(usize, ix)],
                y[scast(usize, iy)],
                ctx,
            ) catch unreachable;

            ops.set( // y[iy] = temp
                &y[scast(usize, iy)],
                temp,
                ctx,
            ) catch unreachable;

            ix += incx;
            iy += incy;
        }
    }
}
