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
const Transpose = linalg.Transpose;

pub inline fn gbmv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_gbmv(transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, ctx);
    } else {
        return k_gbmv(transa.invert(), n, m, ku, kl, alpha, a, lda, x, incx, beta, y, incy, ctx);
    }
}

fn k_gbmv(
    transa: Transpose,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
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

    if (m < 0 or n < 0 or kl < 0 or ku < 0 or lda < (kl + ku + 1) or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or
        (ops.eq(alpha, 0, ctx) catch unreachable and ops.eq(beta, 1, ctx) catch unreachable))
        return;

    const noconj: bool = transa == .no_trans or transa == .trans;

    // Set lenx and leny, the lengths of the vectors x and y, and set up the
    // start points in x and y.
    var lenx: isize = 0;
    var leny: isize = 0;
    if (transa == .no_trans or transa == .conj_no_trans) {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    var kx: isize = if (incx < 0) (-lenx + 1) * incx else 0;
    var ky: isize = if (incy < 0) (-leny + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        // First form y = beta * y.
        if (ops.ne(beta, 1, ctx) catch unreachable) {
            if (incy == 1) {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    for (0..scast(usize, leny)) |i| {
                        ops.set( // y[i] = 0
                            &y[i],
                            0,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    for (0..scast(usize, leny)) |i| {
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
                    for (0..scast(usize, leny)) |_| {
                        ops.set( // y[iy] = 0
                            &y[scast(usize, iy)],
                            0,
                            ctx,
                        ) catch unreachable;

                        iy += incy;
                    }
                } else {
                    for (0..scast(usize, leny)) |_| {
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

        if (transa == .no_trans or transa == .conj_no_trans) {
            // Form  y = alpha * A * x + y  or  y = alpha * A*** * x + y.
            var jx: isize = kx;
            if (incy == 1) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    const temp: C1 = ops.mul( // temp = alpha * x[jx]
                        alpha,
                        x[scast(usize, jx)],
                        ctx,
                    ) catch unreachable;

                    const k: isize = ku - j;
                    if (noconj) {
                        var i: isize = int.max(0, j - ku);
                        while (i < int.min(m, j + kl + 1)) : (i += 1) {
                            ops.add_( // y[i] += temp * a[k + i + j * lda]
                                &y[scast(usize, i)],
                                y[scast(usize, i)],
                                ops.mul(
                                    temp,
                                    a[scast(usize, k + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    } else {
                        var i: isize = int.max(0, j - ku);
                        while (i < int.min(m, j + kl + 1)) : (i += 1) {
                            ops.add_( // y[i] += temp * conj(a[k + i + j * lda])
                                &y[scast(usize, i)],
                                y[scast(usize, i)],
                                ops.mul(
                                    temp,
                                    ops.conjugate(a[scast(usize, k + i + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }

                    jx += incx;
                }
            } else {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    const temp: C1 = ops.mul( // temp = alpha * x[jx]
                        alpha,
                        x[scast(usize, jx)],
                        ctx,
                    ) catch unreachable;

                    const k: isize = ku - j;
                    if (noconj) {
                        var iy: isize = ky;
                        var i: isize = int.max(0, j - ku);
                        while (i < int.min(m, j + kl + 1)) : (i += 1) {
                            ops.add_( // y[iy] += temp * a[k + i + j * lda]
                                &y[scast(usize, iy)],
                                y[scast(usize, iy)],
                                ops.mul(
                                    temp,
                                    a[scast(usize, k + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            iy += incy;
                        }
                    } else {
                        var iy: isize = ky;
                        var i: isize = int.max(0, j - ku);
                        while (i < int.min(m, j + kl + 1)) : (i += 1) {
                            ops.add_( // y[iy] += temp * a[k + i + j * lda]
                                &y[scast(usize, iy)],
                                y[scast(usize, iy)],
                                ops.mul(
                                    temp,
                                    ops.conjugate(a[scast(usize, k + i + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            iy += incy;
                        }
                    }

                    jx += incx;

                    if (j >= ku) ky += incy;
                }
            }
        } else {
            // Form  y = alpha * A**T * x + y  or  y = alpha * A**H * x + y.
            var jy: isize = ky;
            if (incx == 1) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    var temp: C2 = constants.zero(C2, ctx) catch unreachable;

                    const k: isize = ku - j;
                    if (noconj) {
                        var i: isize = int.max(0, j - ku);
                        while (i < int.min(m, j + kl)) : (i += 1) {
                            ops.add_( // temp += a[k + i + j * lda] * x[i]
                                &temp,
                                temp,
                                ops.mul(
                                    a[scast(usize, k + i + j * lda)],
                                    x[scast(usize, i)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    } else {
                        var i: isize = int.max(0, j - ku);
                        while (i < int.min(m, j + kl + 1)) : (i += 1) {
                            ops.add_( // temp += a[k + i + j * lda] * x[i]
                                &temp,
                                temp,
                                ops.mul(
                                    ops.conjugate(a[scast(usize, k + i + j * lda)], ctx) catch unreachable,
                                    x[scast(usize, i)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }

                    ops.add_( // y[jy] += alpha * temp
                        &y[scast(usize, jy)],
                        y[scast(usize, jy)],
                        ops.mul(
                            alpha,
                            temp,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    jy += incy;
                }
            } else {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    var temp: C2 = constants.zero(C2, ctx) catch unreachable;

                    const k: isize = ku - j;
                    if (noconj) {
                        var ix: isize = kx;
                        var i: isize = int.max(0, j - ku);
                        while (i < int.min(m, j + kl + 1)) : (i += 1) {
                            ops.add_( // temp += a[k + i + j * lda] * x[ix]
                                &temp,
                                temp,
                                ops.mul(
                                    a[scast(usize, k + i + j * lda)],
                                    x[scast(usize, ix)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ix += incx;
                        }
                    } else {
                        var ix: isize = kx;
                        var i: isize = int.max(0, j - ku);
                        while (i < int.min(m, j + kl + 1)) : (i += 1) {
                            ops.add_( // temp += a[k + i + j * lda] * x[ix]
                                &temp,
                                temp,
                                ops.mul(
                                    ops.conjugate(a[scast(usize, k + i + j * lda)], ctx) catch unreachable,
                                    x[scast(usize, ix)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ix += incx;
                        }
                    }

                    ops.add_( // y[jy] += alpha * temp
                        &y[scast(usize, jy)],
                        y[scast(usize, jy)],
                        ops.mul(
                            alpha,
                            temp,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    jy += incy;

                    if (j >= ku) kx += incx;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.gbmv not implemented for arbitrary precision types yet");
    }

    return;
}
