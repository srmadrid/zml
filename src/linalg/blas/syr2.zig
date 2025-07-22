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

pub inline fn syr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_syr2(uplo, n, alpha, x, incx, y, incy, a, lda, ctx);
    } else {
        return k_syr2(uplo.invert(), n, alpha, y, incy, x, incx, a, lda, ctx);
    }
}

pub fn k_syr2(
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const C2: type = types.Coerce(Al, X);
    const Y: type = types.Child(@TypeOf(y));
    const C1: type = types.Coerce(Al, Y);
    const A: type = types.Child(@TypeOf(a));
    const CC: type = types.Coerce(Al, types.Coerce(X, types.Coerce(Y, A)));

    if (n < 0 or lda < int.max(1, n) or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or ops.eq(alpha, 0, ctx) catch unreachable)
        return;

    const kx: isize = if (incx < 0) (-n + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
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

                        var i: isize = 0;
                        while (i < j) : (i += 1) {
                            ops.add_( // a[i + j * lda] += x[i] * temp1 + y[i] * temp2
                                &a[scast(usize, i + j * lda)],
                                a[scast(usize, i + j * lda)],
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
                        }

                        ops.add_( // a[j + j * lda] += x[j] * temp1 + y[j] * temp2
                            &a[scast(usize, j + j * lda)],
                            a[scast(usize, j + j * lda)],
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
                        var i: isize = 0;
                        while (i < j) : (i += 1) {
                            ops.add_( // a[i + j * lda] += x[ix] * temp1 + y[iy] * temp2
                                &a[scast(usize, i + j * lda)],
                                a[scast(usize, i + j * lda)],
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

                        ops.add_( // a[j + j * lda] += x[jx] * temp1 + y[jy] * temp2
                            &a[scast(usize, j + j * lda)],
                            a[scast(usize, j + j * lda)],
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

                        ops.add_( // a[j + j * lda] += x[j] * temp1 + y[j] * temp2
                            &a[scast(usize, j + j * lda)],
                            a[scast(usize, j + j * lda)],
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

                        var i: isize = j + 1;
                        while (i < n) : (i += 1) {
                            ops.add_( // a[i + j * lda] += x[i] * temp1 + y[i] * temp2
                                &a[scast(usize, i + j * lda)],
                                a[scast(usize, i + j * lda)],
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
                        }
                    }
                }
            } else {
                var jx: isize = kx;
                var jy: isize = ky;
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable or
                        ops.ne(y[scast(usize, jy)], 0, ctx) catch unreachable)
                    {
                        const temp1: C1 = ops.mul( // temp1 = alpha * y[jx]
                            y[scast(usize, jy)],
                            alpha,
                            ctx,
                        ) catch unreachable;
                        const temp2: C2 = ops.mul( // temp2 = alpha * x[jx]
                            x[scast(usize, jx)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // a[j + j * lda] += x[jx] * temp1 + y[jy] * temp2
                            &a[scast(usize, j + j * lda)],
                            a[scast(usize, j + j * lda)],
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
                        var i: isize = j + 1;
                        while (i < n) : (i += 1) {
                            ix += incx;
                            iy += incy;

                            ops.add_( // a[i + j * lda] += x[ix] * temp1 + y[iy] * temp2
                                &a[scast(usize, i + j * lda)],
                                a[scast(usize, i + j * lda)],
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
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.ger not implemented for arbitrary precision types yet");
    }

    return;
}
