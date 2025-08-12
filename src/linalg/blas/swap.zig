const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn swap(
    n: i32,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = types.Coerce(X, Y);

    if (n <= 0) return blas.Error.InvalidArgument;

    if (comptime types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.scal not implemented for arbitrary precision types yet");
    } else {
        var temp: C = ops.init(C, .{}) catch unreachable;

        var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
        var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
        for (0..scast(u32, n)) |_| {
            ops.set( // temp = x[ix]
                &temp,
                x[scast(u32, ix)],
                ctx,
            ) catch unreachable;

            ops.set( // x[ix] = y[iy]
                &x[scast(u32, ix)],
                y[scast(u32, iy)],
                ctx,
            ) catch unreachable;

            ops.set( // y[iy] = temp
                &y[scast(u32, iy)],
                temp,
                ctx,
            ) catch unreachable;

            ix += incx;
            iy += incy;
        }
    }
}
