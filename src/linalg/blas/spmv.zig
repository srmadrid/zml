const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = linalg.Order;
const Uplo = linalg.Uplo;
const Transpose = linalg.Transpose;

pub inline fn spmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    ap: anytype,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_spmv(uplo, n, alpha, ap, x, incx, beta, y, incy, ctx);
    } else {
        return k_spmv(uplo.invert(), n, alpha, ap, x, incx, beta, y, incy, ctx);
    }
}

fn k_spmv(
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    ap: anytype,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(ap));
    const X: type = types.Child(@TypeOf(x));
    const C1: type = types.Coerce(Al, X);
    const C2: type = types.Coerce(A, X);
    const Be: type = @TypeOf(beta);
    const Y: type = types.Child(@TypeOf(y));
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(X, types.Coerce(Be, Y))));

    if (n < 0 or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or
        (ops.eq(alpha, 0, ctx) catch unreachable and ops.eq(beta, 1, ctx) catch unreachable))
        return;

    const kx: isize = if (incx < 0) (-n + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        // First form  y = beta * y.
        if (ops.ne(beta, 1, ctx) catch unreachable) {
            if (incy == 1) {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    for (0..scast(usize, n)) |i| {
                        ops.set( // y[i] = 0
                            &y[i],
                            0,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    for (0..scast(usize, n)) |i| {
                        ops.mul_( // y[i] *= beta
                            &y[i],
                            y[i],
                            beta,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                var iy: isize = ky;
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    for (0..scast(usize, n)) |_| {
                        ops.set( // y[iy] = 0
                            &y[scast(usize, iy)],
                            0,
                            ctx,
                        ) catch unreachable;

                        iy += incy;
                    }
                } else {
                    for (0..scast(usize, n)) |_| {
                        ops.mul_( // y[iy] *= beta
                            &y[scast(usize, iy)],
                            y[scast(usize, iy)],
                            beta,
                            ctx,
                        ) catch unreachable;

                        iy += incy;
                    }
                }
            }
        }

        if (ops.eq(alpha, 0, ctx) catch unreachable) return;

        var kk: isize = 0;
        if (uplo == .upper) {
            if (incx == 1 and incy == 1) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    const temp1: C1 = ops.mul( // temp1 = alpha * x[j]
                        x[scast(usize, j)],
                        alpha,
                        ctx,
                    ) catch unreachable;
                    var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                    var k: isize = kk;
                    var i: isize = 0;
                    while (i < j) : (i += 1) {
                        ops.add_( // y[i] += temp1 * ap[k]
                            &y[scast(usize, i)],
                            y[scast(usize, i)],
                            ops.mul(
                                temp1,
                                ap[scast(usize, k)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // temp2 += ap[k] * x[i]
                            &temp2,
                            temp2,
                            ops.mul(
                                ap[scast(usize, k)],
                                x[scast(usize, i)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        k += 1;
                    }

                    ops.add_( // y[j] += temp1 * ap[kk + j] + alpha * temp2
                        &y[scast(usize, j)],
                        y[scast(usize, j)],
                        ops.add(
                            ops.mul(
                                temp1,
                                ap[scast(usize, kk + j)],
                                ctx,
                            ) catch unreachable,
                            ops.mul(
                                alpha,
                                temp2,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    kk += j + 1;
                }
            } else {
                var jx: isize = kx;
                var jy: isize = ky;
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    const temp1: C1 = ops.mul( // temp1 = alpha * x[jx]
                        x[scast(usize, jx)],
                        alpha,
                        ctx,
                    ) catch unreachable;
                    var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                    var ix: isize = kx;
                    var iy: isize = ky;
                    var k: isize = kk;
                    while (k < kk + j) : (k += 1) {
                        ops.add_( // y[iy] += temp1 * ap[k]
                            &y[scast(usize, iy)],
                            y[scast(usize, iy)],
                            ops.mul(
                                temp1,
                                ap[scast(usize, k)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // temp2 += ap[k] * x[ix]
                            &temp2,
                            temp2,
                            ops.mul(
                                ap[scast(usize, k)],
                                x[scast(usize, ix)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ix += incx;
                        iy += incy;
                    }

                    ops.add_( // y[jy] += temp1 * ap[kk + j] + alpha * temp2
                        &y[scast(usize, jy)],
                        y[scast(usize, jy)],
                        ops.add(
                            ops.mul(
                                temp1,
                                ap[scast(usize, kk + j)],
                                ctx,
                            ) catch unreachable,
                            ops.mul(
                                alpha,
                                temp2,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    jx += incx;
                    jy += incy;
                    kk += j + 1;
                }
            }
        } else {
            if (incx == 1 and incy == 1) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    const temp1: C1 = ops.mul( // temp1 = alpha * x[j]
                        x[scast(usize, j)],
                        alpha,
                        ctx,
                    ) catch unreachable;
                    var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                    ops.add_( // y[j] += temp1 * ap[kk]
                        &y[scast(usize, j)],
                        y[scast(usize, j)],
                        ops.mul(
                            temp1,
                            ap[scast(usize, kk)],
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    var k: isize = kk + 1;
                    var i: isize = j + 1;
                    while (i < n) : (i += 1) {
                        ops.add_( // y[i] += temp1 * ap[k]
                            &y[scast(usize, i)],
                            y[scast(usize, i)],
                            ops.mul(
                                temp1,
                                ap[scast(usize, k)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // temp2 += ap[k] * x[i]
                            &temp2,
                            temp2,
                            ops.mul(
                                ap[scast(usize, k)],
                                x[scast(usize, i)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        k += 1;
                    }

                    ops.add_( // y[j] += alpha * temp2
                        &y[scast(usize, j)],
                        y[scast(usize, j)],
                        ops.mul(
                            alpha,
                            temp2,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    kk += n - j;
                }
            } else {
                var jx: isize = kx;
                var jy: isize = ky;
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    const temp1: C1 = ops.mul( // temp1 = alpha * x[jx]
                        x[scast(usize, jx)],
                        alpha,
                        ctx,
                    ) catch unreachable;
                    var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                    ops.add_( // y[jy] += temp1 * ap[kk]
                        &y[scast(usize, jy)],
                        y[scast(usize, jy)],
                        ops.mul(
                            temp1,
                            ap[scast(usize, kk)],
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

                        ops.add_( // y[iy] += temp1 * ap[k]
                            &y[scast(usize, iy)],
                            y[scast(usize, iy)],
                            ops.mul(
                                temp1,
                                ap[scast(usize, k)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // temp2 += ap[k] * x[ix]
                            &temp2,
                            temp2,
                            ops.mul(
                                ap[scast(usize, k)],
                                x[scast(usize, ix)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    ops.add_( // y[jy] += alpha * temp2
                        &y[scast(usize, jy)],
                        y[scast(usize, jy)],
                        ops.mul(
                            alpha,
                            temp2,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    jx += incx;
                    jy += incy;
                    kk += n - j;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.hemv not implemented for arbitrary precision types yet");
    }

    return;
}
