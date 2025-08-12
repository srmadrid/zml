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

pub inline fn syr(
    order: Order,
    uplo: Uplo,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_syr(uplo, n, alpha, x, incx, a, lda, ctx);
    } else {
        return k_syr(uplo.invert(), n, alpha, x, incx, a, lda, ctx);
    }
}

fn k_syr(
    uplo: Uplo,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const C1: type = types.Coerce(Al, X);
    const A: type = types.Child(@TypeOf(a));
    const CC: type = types.Coerce(Al, types.Coerce(X, A));

    if (n < 0 or lda < int.max(1, n) or incx == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or ops.eq(alpha, 0, ctx) catch unreachable)
        return;

    const kx: i32 = if (incx < 0) (-n + 1) * incx else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (uplo == .upper) {
            if (incx == 1) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable) {
                        const temp: C1 = ops.mul( // temp = alpha * x[j]
                            x[scast(u32, j)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        var i: i32 = 0;
                        while (i < j) : (i += 1) {
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

                        ops.add_( // a[j + j * lda] += x[j] * temp
                            &a[scast(u32, j + j * lda)],
                            a[scast(u32, j + j * lda)],
                            ops.mul(
                                x[scast(u32, j)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                var jx: i32 = kx;
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                        const temp: C1 = ops.mul( // temp = alpha * x[jx]
                            x[scast(u32, jx)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        var ix: i32 = kx;
                        var i: i32 = 0;
                        while (i < j) : (i += 1) {
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

                        ops.add_( // a[j + j * lda] += x[jx] * temp
                            &a[scast(u32, j + j * lda)],
                            a[scast(u32, j + j * lda)],
                            ops.mul(
                                x[scast(u32, jx)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    jx += incx;
                }
            }
        } else {
            if (incx == 1) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable) {
                        const temp: C1 = ops.mul( // temp = alpha * x[j]
                            x[scast(u32, j)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // a[j + j * lda] += x[j] * temp
                            &a[scast(u32, j + j * lda)],
                            a[scast(u32, j + j * lda)],
                            ops.mul(
                                x[scast(u32, j)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var i: i32 = j + 1;
                        while (i < n) : (i += 1) {
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
                }
            } else {
                var jx: i32 = kx;
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                        const temp: C1 = ops.mul( // temp = alpha * x[jx]
                            x[scast(u32, jx)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // a[j + j * lda] += x[jx] * temp
                            &a[scast(u32, j + j * lda)],
                            a[scast(u32, j + j * lda)],
                            ops.mul(
                                x[scast(u32, jx)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var ix: i32 = jx;
                        var i: i32 = j + 1;
                        while (i < n) : (i += 1) {
                            ix += incx;

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
                        }
                    }

                    jx += incx;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.syr not implemented for arbitrary precision types yet");
    }

    return;
}
