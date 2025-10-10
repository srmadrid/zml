const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const float = @import("../../float.zig");

const blas = @import("../blas.zig");

pub fn dotc_sub(
    n: i32,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    ret: anytype,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = types.Coerce(X, Y);
    const R: type = types.Child(@TypeOf(ret));

    try ops.set(ret, 0, ctx);

    if (n <= 0) return blas.Error.InvalidArgument;

    var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
    var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
    if (comptime types.isArbitraryPrecision(R)) {
        if (comptime types.isArbitraryPrecision(C)) {
            // Orientative implementation for arbitrary precision types
            var temp: C = try ops.init(C, ctx);
            defer ops.deinit(temp, ctx);
            for (0..scast(u32, n)) |_| {
                try ops.mul_(
                    &temp,
                    ops.conj(
                        x[scast(u32, ix)],
                        types.mixStructs(ctx, .{ .copy = false }),
                    ) catch unreachable,
                    y[scast(u32, iy)],
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
        } else {
            @compileError("zml.linalg.blas.dotc_sub not implemented for arbitrary precision types yet");
        }
    } else {
        if (comptime types.isArbitraryPrecision(C)) {
            @compileError("zml.linalg.blas.dotc_sub not implemented for arbitrary precision types yet");
        } else {
            for (0..scast(u32, n)) |_| {
                ops.add_( // ret += conj(x[ix]) * y[iy]
                    ret,
                    ret.*,
                    ops.mul(
                        ops.conj(x[scast(u32, ix)], ctx) catch unreachable,
                        y[scast(u32, iy)],
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
