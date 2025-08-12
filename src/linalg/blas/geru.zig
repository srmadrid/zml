const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;

pub inline fn geru(
    order: Order,
    m: i32,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_geru(m, n, alpha, x, incx, y, incy, a, lda, ctx);
    } else {
        return k_geru(n, m, alpha, y, incy, x, incx, a, lda, ctx);
    }
}

fn k_geru(
    m: i32,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C1: type = types.Coerce(Al, Y);
    const A: type = types.Child(@TypeOf(a));
    const CC: type = types.Coerce(Al, types.Coerce(X, types.Coerce(Y, A)));

    if (m < 0 or n < 0 or lda < int.max(1, m) or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or ops.eq(alpha, 0, ctx) catch unreachable)
        return;

    var jy: i32 = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (incx == 1) {
            var j: i32 = 0;
            while (j < n) : (j += 1) {
                if (ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable) {
                    const temp: C1 = ops.mul( // temp = alpha * y[jy]
                        alpha,
                        y[scast(u32, jy)],
                        ctx,
                    ) catch unreachable;

                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // a[i + j * lda] += x[i] * temp
                            &a[scast(u32, i + j * lda)],
                            a[scast(u32, i + j * lda)],
                            ops.mul(
                                x[scast(u32, i)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }

                jy += incy;
            }
        } else {
            const kx: i32 = if (incx < 0) (-m + 1) * incx else 0;

            var j: i32 = 0;
            while (j < n) : (j += 1) {
                if (ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable) {
                    const temp: C1 = ops.mul( // temp = alpha * y[jy]
                        alpha,
                        y[scast(u32, jy)],
                        ctx,
                    ) catch unreachable;

                    var ix: i32 = kx;
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // a[i + j * lda] += x[ix] * temp
                            &a[scast(u32, i + j * lda)],
                            a[scast(u32, i + j * lda)],
                            ops.mul(
                                x[scast(u32, ix)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ix += incx;
                    }
                }

                jy += incy;
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.geru not implemented for arbitrary precision types yet");
    }

    return;
}
