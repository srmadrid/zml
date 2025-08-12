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

pub inline fn hbmv(
    order: Order,
    uplo: Uplo,
    n: i32,
    k: i32,
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
        return k_hbmv(
            uplo,
            n,
            k,
            alpha,
            a,
            lda,
            x,
            incx,
            beta,
            y,
            incy,
            true,
            ctx,
        );
    } else {
        return k_hbmv(
            uplo.invert(),
            n,
            k,
            ops.conjugate(alpha, ctx) catch unreachable,
            a,
            lda,
            x,
            incx,
            ops.conjugate(beta, ctx) catch unreachable,
            y,
            incy,
            false,
            ctx,
        );
    }
}

fn k_hbmv(
    uplo: Uplo,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    x: anytype,
    incx: i32,
    beta: anytype,
    y: anytype,
    incy: i32,
    noconj: bool,
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

    if (n < 0 or k < 0 or lda < (k + 1) or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or
        (ops.eq(alpha, 0, ctx) catch unreachable and ops.eq(beta, 1, ctx) catch unreachable))
        return;

    var kx: i32 = if (incx < 0) (-n + 1) * incx else 0;
    var ky: i32 = if (incy < 0) (-n + 1) * incy else 0;

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
                    if (noconj) {
                        for (0..scast(u32, n)) |i| {
                            ops.mul_( // y[i] *= beta
                                &y[i],
                                y[i],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    } else {
                        for (0..scast(u32, n)) |i| {
                            ops.mul_( // y[i] = conj(y[i]) * beta
                                &y[i],
                                ops.conjugate(y[i], ctx) catch unreachable,
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
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
                    if (noconj) {
                        for (0..scast(u32, n)) |_| {
                            ops.mul_( // y[iy] *= beta
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            iy += incy;
                        }
                    } else {
                        for (0..scast(u32, n)) |_| {
                            ops.mul_( // y[iy] = conj(y[iy]) * beta
                                &y[scast(u32, iy)],
                                ops.conjugate(y[scast(u32, iy)], ctx) catch unreachable,
                                beta,
                                ctx,
                            ) catch unreachable;

                            iy += incy;
                        }
                    }
                }
            }
        }

        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            if (!noconj) {
                if (incy == 1) {
                    for (0..scast(u32, n)) |i| {
                        ops.conjugate_( // y[i] = conj(y[i])
                            &y[i],
                            y[i],
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
                    for (0..scast(u32, n)) |_| {
                        ops.conjugate_( // y[iy] = conj(y[iy])
                            &y[scast(u32, iy)],
                            y[scast(u32, iy)],
                            ctx,
                        ) catch unreachable;

                        iy += incy;
                    }
                }
            }

            return;
        }

        if (uplo == .upper) {
            // Form  y  when upper triangle of A is stored.
            if (noconj) {
                if (incx == 1 and incy == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * x[j]
                            alpha,
                            x[scast(u32, j)],
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        const l: i32 = k - j;
                        var i: i32 = int.max(0, j - k);
                        while (i < j) : (i += 1) {
                            ops.add_( // y[i] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, i)],
                                y[scast(u32, i)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[i]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    x[scast(u32, i)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.add_( // y[j] += temp1 * re(a[k + j * lda]) + alpha * temp2
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.add(
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, k + j * lda)], ctx) catch unreachable,
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
                            alpha,
                            x[scast(u32, jx)],
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        var ix: i32 = kx;
                        var iy: i32 = ky;
                        const l: i32 = k - j;
                        var i: i32 = int.max(0, j - k);
                        while (i < j) : (i += 1) {
                            ops.add_( // y[iy] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[ix]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    x[scast(u32, ix)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ix += incx;
                            iy += incy;
                        }

                        ops.add_( // y[jy] += temp1 * re(a[k + j * lda]) + alpha * temp2
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.add(
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, k + j * lda)], ctx) catch unreachable,
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

                        if (j >= k) {
                            kx += incx;
                            ky += incy;
                        }
                    }
                }
            } else {
                if (incx == 1 and incy == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * conj(x[j])
                            alpha,
                            ops.conjugate(x[scast(u32, j)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        const l: i32 = k - j;
                        var i: i32 = int.max(0, j - k);
                        while (i < j) : (i += 1) {
                            ops.add_( // y[i] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, i)],
                                y[scast(u32, i)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[i]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    ops.conjugate(x[scast(u32, i)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.add_( // y[j] += temp1 * re(a[k + j * lda]) + alpha * temp2
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.add(
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, k + j * lda)], ctx) catch unreachable,
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
                        const temp1: C1 = ops.mul( // temp1 = alpha * conj(x[jx])
                            alpha,
                            ops.conjugate(x[scast(u32, jx)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        var ix: i32 = kx;
                        var iy: i32 = ky;
                        const l: i32 = k - j;
                        var i: i32 = int.max(0, j - k);
                        while (i < j) : (i += 1) {
                            ops.add_( // y[iy] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * conj(x[ix])
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    ops.conjugate(x[scast(u32, ix)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ix += incx;
                            iy += incy;
                        }

                        ops.add_( // y[jy] += temp1 * re(a[k + j * lda]) + alpha * temp2
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.add(
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, k + j * lda)], ctx) catch unreachable,
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

                        if (j >= k) {
                            kx += incx;
                            ky += incy;
                        }
                    }
                }
            }
        } else {
            // Form  y  when lower triangle of A is stored.
            if (noconj) {
                if (incx == 1 and incy == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * x[j]
                            alpha,
                            x[scast(u32, j)],
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        ops.add_( // y[j] += temp1 * re(a[0 + j * lda])
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.mul(
                                temp1,
                                ops.re(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        const l: i32 = -j;
                        var i: i32 = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            ops.add_( // y[i] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, i)],
                                y[scast(u32, i)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[i]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
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
                            alpha,
                            x[scast(u32, jx)],
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        ops.add_( // y[jy] += temp1 * re(a[0 + j * lda])
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.mul(
                                temp1,
                                ops.re(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        const l: i32 = -j;
                        var ix: i32 = jx;
                        var iy: i32 = jy;
                        var i: i32 = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            ix += incx;
                            iy += incy;

                            ops.add_( // y[iy] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[ix]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
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
            } else {
                if (incx == 1 and incy == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * conj(x[j])
                            alpha,
                            ops.conjugate(x[scast(u32, j)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        ops.add_( // y[j] += temp1 * re(a[0 + j * lda])
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.mul(
                                temp1,
                                ops.re(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        const l: i32 = -j;
                        var i: i32 = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            ops.add_( // y[i] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, i)],
                                y[scast(u32, i)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * conj(x[i])
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    ops.conjugate(x[scast(u32, i)], ctx) catch unreachable,
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
                        const temp1: C1 = ops.mul( // temp1 = alpha * conj(x[jx])
                            alpha,
                            ops.conjugate(x[scast(u32, jx)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        ops.add_( // y[jy] += temp1 * re(a[0 + j * lda])
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.mul(
                                temp1,
                                ops.re(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        const l: i32 = -j;
                        var ix: i32 = jx;
                        var iy: i32 = jy;
                        var i: i32 = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            ix += incx;
                            iy += incy;

                            ops.add_( // y[iy] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * conj(x[ix])
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    ops.conjugate(x[scast(u32, ix)], ctx) catch unreachable,
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
        }

        if (!noconj) {
            if (incy == 1) {
                for (0..scast(u32, n)) |i| {
                    ops.conjugate_( // y[i] = conj(y[i])
                        &y[i],
                        y[i],
                        ctx,
                    ) catch unreachable;
                }
            } else {
                var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
                for (0..scast(u32, n)) |_| {
                    ops.conjugate_( // y[iy] = conj(y[iy])
                        &y[scast(u32, iy)],
                        y[scast(u32, iy)],
                        ctx,
                    ) catch unreachable;

                    iy += incy;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.hbmv not implemented for arbitrary precision types yet");
    }

    return;
}
