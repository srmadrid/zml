const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn asum_sub(
    n: i32,
    x: anytype,
    incx: i32,
    ret: anytype,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const R: type = types.Child(@TypeOf(ret));

    try ops.set(ret, 0, ctx);

    if (n <= 0 or incx <= 0)
        return blas.Error.InvalidArgument;

    var ix: i32 = 0;
    if (comptime types.isArbitraryPrecision(R)) {
        if (comptime types.isArbitraryPrecision(X)) {
            // Orientative implementation for arbitrary precision types
            if (comptime types.isComplex(X)) {
                var temp: Scalar(X) = try ops.init(Scalar(X), ctx);
                defer ops.deinit(&temp, ctx);
                for (0..scast(u32, n)) |_| {
                    try ops.add_(
                        &temp,
                        ops.abs(x[scast(u32, ix)].re, types.mixStructs(ctx, .{ .copy = false })) catch unreachable,
                        ops.abs(x[scast(u32, ix)].im, types.mixStructs(ctx, .{ .copy = false })) catch unreachable,
                        ctx,
                    );

                    try ops.add_(ret, ret.*, temp, ctx);

                    ix += incx;
                }
            } else {
                for (0..scast(u32, n)) |_| {
                    try ops.add_(
                        ret,
                        ret.*,
                        ops.abs(x[scast(u32, ix)], types.mixStructs(ctx, .{ .copy = false })) catch unreachable,
                        ctx,
                    );

                    ix += incx;
                }
            }

            @compileError("zml.linalg.blas.asum_sub not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.asum_sub not implemented for arbitrary precision types yet");
        }
    } else {
        if (comptime types.isArbitraryPrecision(X)) {
            @compileError("zml.linalg.blas.asum_sub not implemented for arbitrary precision types yet");
        } else {
            if (comptime types.isComplex(X)) {
                for (0..scast(u32, n)) |_| {
                    ops.add_( // ret += |x[ix].re| + |x[ix].im|
                        ret,
                        ret.*,
                        ops.add(
                            ops.abs(x[scast(u32, ix)].re, ctx) catch unreachable,
                            ops.abs(x[scast(u32, ix)].im, ctx) catch unreachable,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    ix += incx;
                }
            } else {
                for (0..scast(u32, n)) |_| {
                    ops.add_( // ret += |x[ix]|
                        ret,
                        ret.*,
                        ops.abs(x[scast(u32, ix)], ctx) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    ix += incx;
                }
            }
        }
    }
}
