const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const blas = @import("../blas.zig");

pub fn rot(
    n: i32,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
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

    var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
    var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
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
        for (0..scast(u32, n)) |_| {
            try ops.mul_(
                &temp1,
                c,
                x[scast(u32, ix)],
                ctx,
            );
            try ops.mul_(
                &temp2,
                s,
                y[scast(u32, iy)],
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
                x[scast(u32, ix)],
                ctx,
            );

            try ops.set(
                &x[scast(u32, ix)],
                temp1,
                ctx,
            );

            try ops.mul_(
                &temp1,
                c,
                y[scast(u32, iy)],
                ctx,
            );
            try ops.sub_(
                &temp1,
                temp1,
                temp2,
                ctx,
            );
            try ops.set(
                &y[scast(u32, iy)],
                temp1,
                ctx,
            );

            ix += incx;
            iy += incy;
        }

        @compileError("zml.linalg.blas.rot not implemented for arbitrary precision types yet");
    } else {
        var temp: Ca = try ops.init(Ca, .{});
        for (0..scast(u32, n)) |_| {
            ops.add_( // temp = c * x[ix] + s * y[iy]
                &temp,
                ops.mul(
                    c,
                    x[scast(u32, ix)],
                    ctx,
                ) catch unreachable,
                ops.mul(
                    s,
                    y[scast(u32, iy)],
                    ctx,
                ) catch unreachable,
                ctx,
            ) catch unreachable;

            ops.sub_( // y[iy] = c * y[iy] - s * x[ix]
                &y[scast(u32, iy)],
                ops.mul(
                    c,
                    y[scast(u32, iy)],
                    ctx,
                ) catch unreachable,
                ops.mul(
                    s,
                    x[scast(u32, ix)],
                    ctx,
                ) catch unreachable,
                ctx,
            ) catch unreachable;

            ops.set( // x[ix] = temp
                &x[scast(u32, ix)],
                temp,
                ctx,
            ) catch unreachable;

            ix += incx;
            iy += incy;
        }
    }
}
