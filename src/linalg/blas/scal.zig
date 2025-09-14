const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn scal(
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const C: type = types.Coerce(Al, X);

    if (n < 0 or incx <= 0) return blas.Error.InvalidArgument;

    if (n == 0) return;

    if (ops.eq(alpha, 1, .{}) catch unreachable) return;

    if (comptime types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.scal not implemented for arbitrary precision types yet");
    } else {
        var ix: i32 = 0;
        for (0..scast(u32, n)) |_| {
            ops.mul_( // x[ix] *= alpha
                &x[scast(u32, ix)],
                x[scast(u32, ix)],
                alpha,
                ctx,
            ) catch unreachable;

            ix += incx;
        }
    }
}
