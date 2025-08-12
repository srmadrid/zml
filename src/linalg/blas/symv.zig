const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;
const Uplo = types.Uplo;

pub inline fn symv(
    order: Order,
    uplo: Uplo,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    x: anytype,
    incx: i32,
    beta: anytype,
    y: anytype,
    incy: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_symv(
            uplo,
            n,
            alpha,
            a,
            lda,
            x,
            incx,
            beta,
            y,
            incy,
            ctx,
        );
    } else {
        return k_symv(
            uplo.invert(),
            n,
            alpha,
            a,
            lda,
            x,
            incx,
            beta,
            y,
            incy,
            ctx,
        );
    }
}

fn k_symv(
    uplo: Uplo,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    x: anytype,
    incx: i32,
    beta: anytype,
    y: anytype,
    incy: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const X: type = types.Child(@TypeOf(x));
    const C1: type = types.Coerce(Al, X);
    const C2: type = types.Coerce(A, X);
    const Be: type = @TypeOf(beta);
    const Y: type = types.Child(@TypeOf(y));
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(X, types.Coerce(Be, Y))));

    if (n < 0 or lda < int.max(1, n) or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or
        (ops.eq(alpha, 0, ctx) catch unreachable and ops.eq(beta, 1, ctx) catch unreachable))
        return;

    const kx: i32 = if (incx < 0) (-n + 1) * incx else 0;
    const ky: i32 = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        // First form  y = beta * y.
        if (ops.ne(beta, 1, ctx) catch unreachable) {
            if (incy == 1) {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    for (0..scast(u32, n)) |i| {
                        ops.set( // y[i] = 0
                            &y[i],
                            0,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    for (0..scast(u32, n)) |i| {
                        ops.mul_( // y[i] *= beta
                            &y[i],
                            y[i],
                            beta,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                var iy: i32 = ky;
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    for (0..scast(u32, n)) |_| {
                        ops.set( // y[iy] = 0
                            &y[scast(u32, iy)],
                            0,
                            ctx,
                        ) catch unreachable;

                        iy += incy;
                    }
                } else {
                    for (0..scast(u32, n)) |_| {
                        ops.mul_( // y[iy] *= beta
                            &y[scast(u32, iy)],
                            y[scast(u32, iy)],
                            beta,
                            ctx,
                        ) catch unreachable;

                        iy += incy;
                    }
                }
            }
        }

        if (ops.eq(alpha, 0, ctx) catch unreachable) return;

        if (uplo == .upper) {
            if (incx == 1 and incy == 1) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    const temp1: C1 = ops.mul( // temp1 = alpha * x[j]
                        x[scast(u32, j)],
                        alpha,
                        ctx,
                    ) catch unreachable;
                    var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                    var i: i32 = 0;
                    while (i < j) : (i += 1) {
                        ops.add_( // y[i] += temp1 * a[i + j * lda]
                            &y[scast(u32, i)],
                            y[scast(u32, i)],
                            ops.mul(
                                a[scast(u32, i + j * lda)],
                                temp1,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // temp2 += a[i + j * lda] * x[i]
                            &temp2,
                            temp2,
                            ops.mul(
                                a[scast(u32, i + j * lda)],
                                x[scast(u32, i)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    ops.add_( // y[j] += temp1 * a[j + j * lda] + alpha * temp2
                        &y[scast(u32, j)],
                        y[scast(u32, j)],
                        ops.add(
                            ops.mul(
                                temp1,
                                a[scast(u32, j + j * lda)],
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
                }
            } else {
                var jx: i32 = kx;
                var jy: i32 = ky;
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    const temp1: C1 = ops.mul( // temp1 = alpha * x[jx]
                        x[scast(u32, jx)],
                        alpha,
                        ctx,
                    ) catch unreachable;
                    var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                    var ix: i32 = kx;
                    var iy: i32 = ky;
                    var i: i32 = 0;
                    while (i < j) : (i += 1) {
                        ops.add_( // y[iy] += temp1 * a[i + j * lda]
                            &y[scast(u32, iy)],
                            y[scast(u32, iy)],
                            ops.mul(
                                a[scast(u32, i + j * lda)],
                                temp1,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // temp2 += a[i + j * lda] * x[ix]
                            &temp2,
                            temp2,
                            ops.mul(
                                a[scast(u32, i + j * lda)],
                                x[scast(u32, ix)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ix += incx;
                        iy += incy;
                    }

                    ops.add_( // y[jy] += temp1 * a[j + j * lda] + alpha * temp2
                        &y[scast(u32, jy)],
                        y[scast(u32, jy)],
                        ops.add(
                            ops.mul(
                                temp1,
                                a[scast(u32, j + j * lda)],
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
                }
            }
        } else {
            if (incx == 1 and incy == 1) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    const temp1: C1 = ops.mul( // temp1 = alpha * x[j]
                        x[scast(u32, j)],
                        alpha,
                        ctx,
                    ) catch unreachable;
                    var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                    ops.add_( // y[j] += temp1 * a[j + j * lda]
                        &y[scast(u32, j)],
                        y[scast(u32, j)],
                        ops.mul(
                            temp1,
                            a[scast(u32, j + j * lda)],
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    var i: i32 = j + 1;
                    while (i < n) : (i += 1) {
                        ops.add_( // y[i] += temp1 * a[i + j * lda]
                            &y[scast(u32, i)],
                            y[scast(u32, i)],
                            ops.mul(
                                a[scast(u32, i + j * lda)],
                                temp1,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // temp2 += a[i + j * lda] * x[i]
                            &temp2,
                            temp2,
                            ops.mul(
                                a[scast(u32, i + j * lda)],
                                x[scast(u32, i)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    ops.add_( // y[j] += alpha * temp2
                        &y[scast(u32, j)],
                        y[scast(u32, j)],
                        ops.mul(
                            alpha,
                            temp2,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;
                }
            } else {
                var jx: i32 = kx;
                var jy: i32 = ky;
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    const temp1: C1 = ops.mul( // temp1 = alpha * x[jx]
                        x[scast(u32, jx)],
                        alpha,
                        ctx,
                    ) catch unreachable;
                    var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                    ops.add_( // y[jy] += temp1 * a[j + j * lda]
                        &y[scast(u32, jy)],
                        y[scast(u32, jy)],
                        ops.mul(
                            temp1,
                            a[scast(u32, j + j * lda)],
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    var ix: i32 = jx;
                    var iy: i32 = jy;
                    var i: i32 = j + 1;
                    while (i < n) : (i += 1) {
                        ix += incx;
                        iy += incy;

                        ops.add_( // y[iy] += temp1 * a[i + j * lda]
                            &y[scast(u32, iy)],
                            y[scast(u32, iy)],
                            ops.mul(
                                a[scast(u32, i + j * lda)],
                                temp1,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // temp2 += a[i + j * lda] * x[ix]
                            &temp2,
                            temp2,
                            ops.mul(
                                a[scast(u32, i + j * lda)],
                                x[scast(u32, ix)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    ops.add_( // y[jy] += alpha * temp2
                        &y[scast(u32, jy)],
                        y[scast(u32, jy)],
                        ops.mul(
                            alpha,
                            temp2,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    jx += incx;
                    jy += incy;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.symv not implemented for arbitrary precision types yet");
    }

    return;
}
