const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = linalg.Order;
const Uplo = linalg.Uplo;

pub inline fn spr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ap: anytype,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_spr2(uplo, n, alpha, x, incx, y, incy, ap, ctx);
    } else {
        return k_spr2(uplo.invert(), n, alpha, y, incy, x, incx, ap, ctx);
    }
}

fn k_spr2(
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ap: anytype,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const C2: type = types.Coerce(Al, X);
    const Y: type = types.Child(@TypeOf(y));
    const C1: type = types.Coerce(Al, Y);
    const A: type = types.Child(@TypeOf(ap));
    const CC: type = types.Coerce(Al, types.Coerce(X, types.Coerce(Y, A)));

    if (n < 0 or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or ops.eq(alpha, 0, ctx) catch unreachable)
        return;

    const kx: isize = if (incx < 0) (-n + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        var kk: isize = 0;
        if (uplo == .upper) {
            if (incx == 1 and incy == 1) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable or
                        ops.ne(y[scast(usize, j)], 0, ctx) catch unreachable)
                    {
                        const temp1: C1 = ops.mul( // temp1 = alpha * y[j]
                            y[scast(usize, j)],
                            alpha,
                            ctx,
                        ) catch unreachable;
                        const temp2: C2 = ops.mul( // temp2 = alpha * x[j]
                            x[scast(usize, j)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        var k: isize = kk;
                        var i: isize = 0;
                        while (i < j) : (i += 1) {
                            ops.add_( // ap[k] += x[i] * temp1 + y[i] * temp2
                                &ap[scast(usize, k)],
                                ap[scast(usize, k)],
                                ops.add(
                                    ops.mul(
                                        x[scast(usize, i)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(usize, i)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            k += 1;
                        }

                        ops.add_( // ap[kk + j] += x[j] * temp1 + y[j] * temp2
                            &ap[scast(usize, kk + j)],
                            ap[scast(usize, kk + j)],
                            ops.add(
                                ops.mul(
                                    x[scast(usize, j)],
                                    temp1,
                                    ctx,
                                ) catch unreachable,
                                ops.mul(
                                    y[scast(usize, j)],
                                    temp2,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    kk += j + 1;
                }
            } else {
                var jx: isize = kx;
                var jy: isize = ky;
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable or
                        ops.ne(y[scast(usize, jy)], 0, ctx) catch unreachable)
                    {
                        const temp1: C1 = ops.mul( // temp1 = alpha * y[jy]
                            y[scast(usize, jy)],
                            alpha,
                            ctx,
                        ) catch unreachable;
                        const temp2: C2 = ops.mul( // temp2 = alpha * x[jx]
                            x[scast(usize, jx)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        var ix: isize = kx;
                        var iy: isize = ky;
                        var k: isize = kk;
                        while (k < kk + j) : (k += 1) {
                            ops.add_( // ap[k] += x[ix] * temp1 + y[iy] * temp2
                                &ap[scast(usize, k)],
                                ap[scast(usize, k)],
                                ops.add(
                                    ops.mul(
                                        x[scast(usize, ix)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(usize, iy)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ix += incx;
                            iy += incy;
                        }

                        ops.add_( // ap[kk + j] += x[jx] * temp1 + y[jy] * temp2
                            &ap[scast(usize, kk + j)],
                            ap[scast(usize, kk + j)],
                            ops.add(
                                ops.mul(
                                    x[scast(usize, jx)],
                                    temp1,
                                    ctx,
                                ) catch unreachable,
                                ops.mul(
                                    y[scast(usize, jy)],
                                    temp2,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    jx += incx;
                    jy += incy;
                    kk += j + 1;
                }
            }
        } else {
            if (incx == 1 and incy == 1) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable or
                        ops.ne(y[scast(usize, j)], 0, ctx) catch unreachable)
                    {
                        const temp1: C1 = ops.mul( // temp1 = alpha * y[j]
                            y[scast(usize, j)],
                            alpha,
                            ctx,
                        ) catch unreachable;
                        const temp2: C2 = ops.mul( // temp2 = alpha * x[j]
                            x[scast(usize, j)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // ap[kk] += x[j] * temp1 + y[j] * temp2
                            &ap[scast(usize, kk)],
                            ap[scast(usize, kk)],
                            ops.add(
                                ops.mul(
                                    x[scast(usize, j)],
                                    temp1,
                                    ctx,
                                ) catch unreachable,
                                ops.mul(
                                    y[scast(usize, j)],
                                    temp2,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var k: isize = kk + 1;
                        var i: isize = j + 1;
                        while (i < n) : (i += 1) {
                            ops.add_( // ap[k] += x[i] * temp1 + y[i] * temp2
                                &ap[scast(usize, k)],
                                ap[scast(usize, k)],
                                ops.add(
                                    ops.mul(
                                        x[scast(usize, i)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(usize, i)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            k += 1;
                        }
                    }

                    kk += n - j;
                }
            } else {
                var jx: isize = kx;
                var jy: isize = ky;
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable or
                        ops.ne(y[scast(usize, jy)], 0, ctx) catch unreachable)
                    {
                        const temp1: C1 = ops.mul( // temp1 = alpha * y[jy]
                            y[scast(usize, jy)],
                            alpha,
                            ctx,
                        ) catch unreachable;
                        const temp2: C2 = ops.mul( // temp2 = alpha * x[jx]
                            x[scast(usize, jx)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // ap[kk] += x[jx] * temp1 + y[jy] * temp2
                            &ap[scast(usize, kk)],
                            ap[scast(usize, kk)],
                            ops.add(
                                ops.mul(
                                    x[scast(usize, jx)],
                                    temp1,
                                    ctx,
                                ) catch unreachable,
                                ops.mul(
                                    y[scast(usize, jy)],
                                    temp2,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var ix: isize = jx;
                        var iy: isize = jy;
                        var k: isize = kk + 1;
                        while (k < kk + n - j) : (k += 1) {
                            ix += incx;
                            iy += incy;

                            ops.add_( // ap[k] += x[ix] * temp1 + y[iy] * temp2
                                &ap[scast(usize, k)],
                                ap[scast(usize, k)],
                                ops.add(
                                    ops.mul(
                                        x[scast(usize, ix)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(usize, iy)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }

                    jx += incx;
                    jy += incy;
                    kk += n - j;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.spr2 not implemented for arbitrary precision types yet");
    }

    return;
}
