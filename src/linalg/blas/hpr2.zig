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

pub inline fn hpr2(
    order: Order,
    uplo: Uplo,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    ap: anytype,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_hpr2(
            uplo,
            n,
            alpha,
            x,
            incx,
            y,
            incy,
            ap,
            true,
            ctx,
        );
    } else {
        return k_hpr2(
            uplo.invert(),
            n,
            alpha,
            y,
            incy,
            x,
            incx,
            ap,
            false,
            ctx,
        );
    }
}

fn k_hpr2(
    uplo: Uplo,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    ap: anytype,
    noconj: bool,
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

    const kx: i32 = if (incx < 0) (-n + 1) * incx else 0;
    const ky: i32 = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        var kk: i32 = 0;
        if (uplo == .upper) {
            if (noconj) {
                if (incx == 1 and incy == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, j)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * conj(y[j])
                                ops.conj(y[scast(u32, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * x[j])
                                x[scast(u32, j)],
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var k: i32 = kk;
                            var i: i32 = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // ap[k] += x[i] * temp1 + y[i] * temp2
                                    &ap[scast(u32, k)],
                                    ap[scast(u32, k)],
                                    ops.add(
                                        ops.mul(
                                            x[scast(u32, i)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            y[scast(u32, i)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }

                            ops.add_( // ap[kk + j] = re(ap[kk + j]) + re(x[j] * temp1 + y[j] * temp2)
                                &ap[scast(u32, kk + j)],
                                ops.re(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        x[scast(u32, j)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(u32, j)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // ap[kk + j] = re(ap[kk + j])
                                &ap[scast(u32, kk + j)],
                                ops.re(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        kk += j + 1;
                    }
                } else {
                    var jx: i32 = kx;
                    var jy: i32 = ky;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * conj(y[jy])
                                ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * x[jx])
                                x[scast(u32, jx)],
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var ix: i32 = kx;
                            var iy: i32 = ky;
                            var k: i32 = kk;
                            while (k < kk + j) : (k += 1) {
                                ops.add_( // ap[k] += x[ix] * temp1 + y[iy] * temp2
                                    &ap[scast(u32, k)],
                                    ap[scast(u32, k)],
                                    ops.add(
                                        ops.mul(
                                            x[scast(u32, ix)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            y[scast(u32, iy)],
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

                            ops.add_( // ap[kk + j] = re(ap[kk + j]) + re(x[jx] * temp1 + y[jy] * temp2)
                                &ap[scast(u32, kk + j)],
                                ops.re(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        x[scast(u32, jx)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(u32, jy)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // ap[kk + j] = re(ap[kk + j])
                                &ap[scast(u32, kk + j)],
                                ops.re(ap[scast(u32, kk + j)], ctx) catch unreachable,
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, j)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * y[j]
                                y[scast(u32, j)],
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * conj(x[j]))
                                ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var k: i32 = kk;
                            var i: i32 = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // ap[k] += conj(x[i]) * temp1 + conj(y[i]) * temp2
                                    &ap[scast(u32, k)],
                                    ap[scast(u32, k)],
                                    ops.add(
                                        ops.mul(
                                            ops.conj(x[scast(u32, i)], ctx) catch unreachable,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conj(y[scast(u32, i)], ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }

                            ops.add_( // ap[kk + j] = re(ap[kk + j]) + re(conj(x[j]) * temp1 + conj(y[j]) * temp2)
                                &ap[scast(u32, kk + j)],
                                ops.re(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        ops.conj(y[scast(u32, j)], ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // ap[kk + j] = re(ap[kk + j])
                                &ap[scast(u32, kk + j)],
                                ops.re(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        kk += j + 1;
                    }
                } else {
                    var jx: i32 = kx;
                    var jy: i32 = ky;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * y[jy]
                                y[scast(u32, jy)],
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * conj(x[jx]))
                                ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var ix: i32 = kx;
                            var iy: i32 = ky;
                            var k: i32 = kk;
                            while (k < kk + j) : (k += 1) {
                                ops.add_( // ap[k] += conj(x[ix]) * temp1 + conj(y[iy]) * temp2
                                    &ap[scast(u32, k)],
                                    ap[scast(u32, k)],
                                    ops.add(
                                        ops.mul(
                                            ops.conj(x[scast(u32, ix)], ctx) catch unreachable,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conj(y[scast(u32, iy)], ctx) catch unreachable,
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

                            ops.add_( // ap[kk + j] = re(ap[kk + j]) + re(conj(x[jx]) * temp1 + conj(y[jy]) * temp2)
                                &ap[scast(u32, kk + j)],
                                ops.re(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // ap[kk + j] = re(ap[kk + j])
                                &ap[scast(u32, kk + j)],
                                ops.re(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        jy += incy;
                        kk += j + 1;
                    }
                }
            }
        } else {
            if (noconj) {
                if (incx == 1 and incy == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, j)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * conj(y[j])
                                ops.conj(y[scast(u32, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * x[j])
                                x[scast(u32, j)],
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            ops.add_( // ap[kk] = re(ap[kk]) + re(x[j] * temp1 + y[j] * temp2)
                                &ap[scast(u32, kk)],
                                ops.re(ap[scast(u32, kk)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        x[scast(u32, j)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(u32, j)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var k: i32 = kk + 1;
                            var i: i32 = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // ap[k] += x[i] * temp1 + y[i] * temp2
                                    &ap[scast(u32, k)],
                                    ap[scast(u32, k)],
                                    ops.add(
                                        ops.mul(
                                            x[scast(u32, i)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            y[scast(u32, i)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }
                        } else {
                            ops.set( // ap[kk] = re(ap[kk])
                                &ap[scast(u32, kk)],
                                ops.re(ap[scast(u32, kk)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        kk += n - j;
                    }
                } else {
                    var jx: i32 = kx;
                    var jy: i32 = ky;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * conj(y[jy])
                                ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * x[jx])
                                x[scast(u32, jx)],
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            ops.add_( // ap[kk] = re(ap[kk]) + re(x[jx] * temp1 + y[jy] * temp2)
                                &ap[scast(u32, kk)],
                                ops.re(ap[scast(u32, kk)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        x[scast(u32, jx)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(u32, jy)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var ix: i32 = jx;
                            var iy: i32 = jy;
                            var k: i32 = kk + 1;
                            while (k < kk + n - j) : (k += 1) {
                                ix += incx;
                                iy += incy;

                                ops.add_( // ap[k] += x[ix] * temp1 + y[iy] * temp2
                                    &ap[scast(u32, k)],
                                    ap[scast(u32, k)],
                                    ops.add(
                                        ops.mul(
                                            x[scast(u32, ix)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            y[scast(u32, iy)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // ap[kk] = re(ap[kk])
                                &ap[scast(u32, kk)],
                                ops.re(ap[scast(u32, kk)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        jy += incy;
                        kk += n - j;
                    }
                }
            } else {
                if (incx == 1 and incy == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, j)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * y[j]
                                y[scast(u32, j)],
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * conj(x[j]))
                                ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            ops.add_( // ap[kk] = re(ap[kk]) + re(conj(x[j]) * temp1 + conj(y[j]) * temp2)
                                &ap[scast(u32, kk)],
                                ops.re(ap[scast(u32, kk)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        ops.conj(y[scast(u32, j)], ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var k: i32 = kk + 1;
                            var i: i32 = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // ap[k] += conj(x[i]) * temp1 + conj(y[i]) * temp2
                                    &ap[scast(u32, k)],
                                    ap[scast(u32, k)],
                                    ops.add(
                                        ops.mul(
                                            ops.conj(x[scast(u32, i)], ctx) catch unreachable,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conj(y[scast(u32, i)], ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }
                        } else {
                            ops.set( // ap[kk] = re(ap[kk])
                                &ap[scast(u32, kk)],
                                ops.re(ap[scast(u32, kk)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        kk += n - j;
                    }
                } else {
                    var jx: i32 = kx;
                    var jy: i32 = ky;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * y[jy]
                                y[scast(u32, jy)],
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * conj(x[jx]))
                                ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            ops.add_( // ap[kk] = re(ap[kk]) + re(conj(x[jx]) * temp1 + conj(y[jy]) * temp2)
                                &ap[scast(u32, kk)],
                                ops.re(ap[scast(u32, kk)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var ix: i32 = jx;
                            var iy: i32 = jy;
                            var k: i32 = kk + 1;
                            while (k < kk + n - j) : (k += 1) {
                                ix += incx;
                                iy += incy;

                                ops.add_( // ap[k] += conj(x[ix]) * temp1 + conj(y[iy]) * temp2
                                    &ap[scast(u32, k)],
                                    ap[scast(u32, k)],
                                    ops.add(
                                        ops.mul(
                                            ops.conj(x[scast(u32, ix)], ctx) catch unreachable,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conj(y[scast(u32, iy)], ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // ap[kk] = re(ap[kk])
                                &ap[scast(u32, kk)],
                                ops.re(ap[scast(u32, kk)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        jy += incy;
                        kk += n - j;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.hpr2 not implemented for arbitrary precision types yet");
    }

    return;
}
