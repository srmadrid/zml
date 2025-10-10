const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;

pub inline fn gerc(
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
        return k_gerc(m, n, alpha, x, incx, y, incy, a, lda, true, ctx);
    } else {
        return k_gerc(n, m, alpha, y, incy, x, incx, a, lda, false, ctx);
    }
}

fn k_gerc(
    m: i32,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    a: anytype,
    lda: i32,
    noconj: bool,
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
        if (noconj) {
            if (incx == 1) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable) {
                        const temp: C1 = ops.mul( // temp = alpha * conj(y[jy])
                            alpha,
                            ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
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
                        const temp: C1 = ops.mul( // temp = alpha * conj(y[jy])
                            alpha,
                            ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
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
            if (incx == 1) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable) {
                        const temp: C1 = ops.mul( // temp = alpha * conj(y[jy])
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
                                    ops.conj(x[scast(u32, i)], ctx) catch unreachable,
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
                        const temp: C1 = ops.mul( // temp = alpha * conj(y[jy])
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
                                    ops.conj(x[scast(u32, ix)], ctx) catch unreachable,
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
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.gerc not implemented for arbitrary precision types yet");
    }

    return;
}
