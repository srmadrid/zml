const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn axpy(
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const C: type = types.Coerce(Al, X);
    const Y: type = types.Child(@TypeOf(y));

    if (n <= 0) return blas.Error.InvalidArgument;

    if (ops.eq(alpha, 0, .{}) catch unreachable) return;

    var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
    var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
    if (comptime types.isArbitraryPrecision(C)) {
        if (comptime types.isArbitraryPrecision(Y)) {
            // Orientative implementation for arbitrary precision types
            var temp: C = try ops.init(C, ctx);
            defer ops.deinit(temp, ctx);
            for (0..scast(usize, n)) |_| {
                try ops.mul_(
                    &temp,
                    alpha,
                    x[scast(usize, ix)],
                    ctx,
                );

                try ops.add_(
                    &y[scast(usize, iy)],
                    y[scast(usize, iy)],
                    temp,
                    ctx,
                );

                ix += incx;
                iy += incy;
            }

            @compileError("zml.linalg.blas.axpy not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.axpy not implemented for arbitrary precision types yet");
        }
    } else {
        if (comptime types.isArbitraryPrecision(Y)) {
            @compileError("zml.linalg.blas.axpy not implemented for arbitrary precision types yet");
        } else {
            for (0..scast(usize, n)) |_| {
                ops.add_( // y[iy] += alpha * x[ix]
                    &y[scast(usize, iy)],
                    y[scast(usize, iy)],
                    ops.mul(
                        alpha,
                        x[scast(usize, ix)],
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ix += incx;
                iy += incy;
            }
        }
    }
}
