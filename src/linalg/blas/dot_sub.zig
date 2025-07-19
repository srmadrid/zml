const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const float = @import("../../float.zig");

const blas = @import("../blas.zig");

pub fn dot_sub(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ret: anytype,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = types.Coerce(X, Y);
    const R: type = types.Child(@TypeOf(ret));

    try ops.set(ret, 0, ctx);

    if (n <= 0) return blas.Error.InvalidArgument;

    var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
    var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
    if (comptime types.isArbitraryPrecision(R)) {
        if (comptime types.isArbitraryPrecision(C)) {
            // Orientative implementation for arbitrary precision types
            var temp: C = try ops.init(C, ctx);
            defer ops.deinit(temp, ctx);
            for (0..scast(usize, n)) |_| {
                try ops.mul_(
                    &temp,
                    x[scast(usize, ix)],
                    y[scast(usize, iy)],
                    ctx,
                );

                try ops.add_(
                    ret,
                    ret.*,
                    temp,
                    ctx,
                );

                ix += incx;
                iy += incy;
            }

            @compileError("zml.linalg.blas.dot_sub not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.dot_sub not implemented for arbitrary precision types yet");
        }
    } else {
        if (comptime types.isArbitraryPrecision(C)) {
            @compileError("zml.linalg.blas.dot_sub not implemented for arbitrary precision types yet");
        } else {
            for (0..scast(usize, n)) |_| {
                ops.add_( // ret += x[ix] * y[iy]
                    ret,
                    ret.*,
                    ops.mul(
                        x[scast(usize, ix)],
                        y[scast(usize, iy)],
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
