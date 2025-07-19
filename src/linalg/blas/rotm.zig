const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const blas = @import("../blas.zig");

pub fn rotm(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    param: anytype,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const P: type = types.Child(@TypeOf(param));
    //const C: type = types.Coerce(X, types.Coerce(Y, P));

    if (n == 0 or param[0] == -2) return blas.Error.InvalidArgument;

    if (comptime types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(P))
    {
        @compileError("zml.linalg.blas.rotm not implemented for arbitrary precision types yet");
    } else {
        const flag = param[0];

        if (flag == -1) {
            const h11: P = param[1];
            const h21: P = param[2];
            const h12: P = param[3];
            const h22: P = param[4];
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            for (0..scast(usize, n)) |_| {
                const x0: X = x[scast(usize, ix)];
                ops.add_( // x[ix] = h11 * x[ix] + h12 * y[iy]
                    &x[scast(usize, ix)],
                    ops.mul(
                        h11,
                        x[scast(usize, ix)],
                        ctx,
                    ) catch unreachable,
                    ops.mul(
                        h12,
                        y[scast(usize, iy)],
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ops.add_( // y[iy] = h21 * x[ix] + h22 * y[iy]
                    &y[scast(usize, iy)],
                    ops.mul(
                        h21,
                        x0,
                        ctx,
                    ) catch unreachable,
                    ops.mul(
                        h22,
                        y[scast(usize, iy)],
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ix += incx;
                iy += incy;
            }
        } else if (flag == 0) {
            const h21: P = param[2];
            const h12: P = param[3];
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            for (0..scast(usize, n)) |_| {
                const x0 = x[scast(usize, ix)];
                ops.add_( // x[ix] = h12 * y[iy] + x[ix]
                    &x[scast(usize, ix)],
                    x[scast(usize, ix)],
                    ops.mul(
                        h12,
                        y[scast(usize, iy)],
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ops.add_( // y[iy] = h21 * x[ix] + y[iy]
                    &y[scast(usize, iy)],
                    y[scast(usize, iy)],
                    ops.mul(
                        h21,
                        x0,
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ix += incx;
                iy += incy;
            }
        } else if (flag == 1) {
            const h11 = param[1];
            const h22 = param[4];
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            for (0..scast(usize, n)) |_| {
                const x0 = x[scast(usize, ix)];
                ops.add_( // x[ix] = h11 * x[ix] - h22 * y[iy]
                    &x[scast(usize, ix)],
                    ops.mul(
                        h11,
                        x[scast(usize, ix)],
                        ctx,
                    ) catch unreachable,
                    y[scast(usize, iy)],
                    ctx,
                ) catch unreachable;

                ops.sub_( // y[iy] = h22 * y[iy] - h11 * x[ix]
                    &y[scast(usize, iy)],
                    ops.mul(
                        h22,
                        y[scast(usize, iy)],
                        ctx,
                    ) catch unreachable,
                    x0,
                    ctx,
                ) catch unreachable;

                ix += incx;
                iy += incy;
            }
        }
    }
}
