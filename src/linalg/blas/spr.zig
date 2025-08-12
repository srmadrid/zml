const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;
const Uplo = types.Uplo;

pub inline fn spr(
    order: Order,
    uplo: Uplo,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    ap: anytype,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_spr(uplo, n, alpha, x, incx, ap, ctx);
    } else {
        return k_spr(uplo.invert(), n, alpha, x, incx, ap, ctx);
    }
}

fn k_spr(
    uplo: Uplo,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    ap: anytype,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const C1: type = types.Coerce(Al, X);
    const A: type = types.Child(@TypeOf(ap));
    const CC: type = types.Coerce(Al, types.Coerce(X, A));

    if (n < 0 or incx == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or ops.eq(alpha, 0, ctx) catch unreachable)
        return;

    const kx: i32 = if (incx < 0) (-n + 1) * incx else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        var kk: i32 = 0;
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

                        var k: i32 = kk;
                        var i: i32 = 0;
                        while (i < j) : (i += 1) {
                            ops.add_( // ap[k] += x[i] * temp
                                &ap[scast(u32, k)],
                                ap[scast(u32, k)],
                                ops.mul(
                                    x[scast(u32, i)],
                                    temp,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            k += 1;
                        }

                        ops.add_( // ap[kk + j] += x[j] * temp
                            &ap[scast(u32, kk + j)],
                            ap[scast(u32, kk + j)],
                            ops.mul(
                                x[scast(u32, j)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    kk += j + 1;
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
                        var k: i32 = kk;
                        while (k < kk + j) : (k += 1) {
                            ops.add_( // ap[k] += x[ix] * temp
                                &ap[scast(u32, k)],
                                ap[scast(u32, k)],
                                ops.mul(
                                    x[scast(u32, ix)],
                                    temp,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ix += incx;
                        }

                        ops.add_( // ap[kk + j] += x[jx] * temp
                            &ap[scast(u32, kk + j)],
                            ap[scast(u32, kk + j)],
                            ops.mul(
                                x[scast(u32, jx)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    jx += incx;
                    kk += j + 1;
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

                        ops.add_( // ap[kk] += x[j] * temp
                            &ap[scast(u32, kk)],
                            ap[scast(u32, kk)],
                            ops.mul(
                                x[scast(u32, j)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var k: i32 = kk + 1;
                        var i: i32 = j + 1;
                        while (i < n) : (i += 1) {
                            ops.add_( // ap[k] += x[i] * temp
                                &ap[scast(u32, k)],
                                ap[scast(u32, k)],
                                ops.mul(
                                    x[scast(u32, i)],
                                    temp,
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
                var jx: i32 = kx;
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                        const temp: C1 = ops.mul( // temp = alpha * x[jx]
                            x[scast(u32, jx)],
                            alpha,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // ap[kk] += x[jx] * temp
                            &ap[scast(u32, kk)],
                            ap[scast(u32, kk)],
                            ops.mul(
                                x[scast(u32, jx)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var ix: i32 = jx;
                        var k: i32 = kk + 1;
                        while (k < kk + n - j) : (k += 1) {
                            ix += incx;

                            ops.add_( // ap[k] += x[ix] * temp
                                &ap[scast(u32, k)],
                                ap[scast(u32, k)],
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
                    kk += n - j;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.spr not implemented for arbitrary precision types yet");
    }

    return;
}
