const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const blas = @import("../blas.zig");

pub fn rot(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    c: anytype,
    s: anytype,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = @TypeOf(c);
    const S: type = @TypeOf(s);
    const Ca: type = types.Coerce(X, types.Coerce(Y, types.Coerce(C, S)));

    if (n <= 0) return blas.Error.InvalidArgument;

    if (ops.eq(c, 1, .{}) catch unreachable and ops.eq(s, 0, .{}) catch unreachable) return;

    var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
    var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
    if (comptime types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(C) or
        types.isArbitraryPrecision(S))
    {
        // Orientative implementation for arbitrary precision types
        var temp1: Ca = try ops.init(Ca, ctx);
        defer ops.deinit(&temp1, ctx);
        var temp2: Ca = try ops.init(Ca, ctx);
        defer ops.deinit(&temp2, ctx);
        for (0..scast(usize, n)) |_| {
            try ops.mul_(
                &temp1,
                c,
                x[scast(usize, ix)],
                ctx,
            );
            try ops.mul_(
                &temp2,
                s,
                y[scast(usize, iy)],
                ctx,
            );
            try ops.add_(
                &temp1,
                temp1,
                temp2,
                ctx,
            );

            try ops.mul_(
                &temp2,
                s,
                x[scast(usize, ix)],
                ctx,
            );

            try ops.set(
                &x[scast(usize, ix)],
                temp1,
                ctx,
            );

            try ops.mul_(
                &temp1,
                c,
                y[scast(usize, iy)],
                ctx,
            );
            try ops.sub_(
                &temp1,
                temp1,
                temp2,
                ctx,
            );
            try ops.set(
                &y[scast(usize, iy)],
                temp1,
                ctx,
            );

            ix += incx;
            iy += incy;
        }

        @compileError("zml.linalg.blas.rot not implemented for arbitrary precision types yet");
    } else {
        var temp: Ca = try ops.init(Ca, .{});
        for (0..scast(usize, n)) |_| {
            ops.add_( // temp = c * x[ix] + s * y[iy]
                &temp,
                ops.mul(
                    c,
                    x[scast(usize, ix)],
                    ctx,
                ) catch unreachable,
                ops.mul(
                    s,
                    y[scast(usize, iy)],
                    ctx,
                ) catch unreachable,
                ctx,
            ) catch unreachable;

            ops.sub_( // y[iy] = c * y[iy] - s * x[ix]
                &y[scast(usize, iy)],
                ops.mul(
                    c,
                    y[scast(usize, iy)],
                    ctx,
                ) catch unreachable,
                ops.mul(
                    s,
                    x[scast(usize, ix)],
                    ctx,
                ) catch unreachable,
                ctx,
            ) catch unreachable;

            ops.set( // x[ix] = temp
                &x[scast(usize, ix)],
                temp,
                ctx,
            ) catch unreachable;

            ix += incx;
            iy += incy;
        }
    }
}
