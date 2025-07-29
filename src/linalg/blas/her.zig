const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = linalg.Order;
const Uplo = linalg.Uplo;

pub inline fn her(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_her(uplo, n, alpha, x, incx, a, lda, true, ctx);
    } else {
        return k_her(uplo.invert(), n, alpha, x, incx, a, lda, false, ctx);
    }
}

fn k_her(
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    a: anytype,
    lda: isize,
    noconj: bool,
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

    const kx: isize = if (incx < 0) (-n + 1) * incx else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (uplo == .upper) {
            if (noconj) {
                if (incx == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable) {
                            const temp: C1 = ops.mul( // temp = alpha * conj(x[j])
                                ops.conjugate(x[scast(usize, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;

                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // a[i + j * lda] += x[i] * temp
                                    &a[scast(usize, i + j * lda)],
                                    a[scast(usize, i + j * lda)],
                                    ops.mul(
                                        x[scast(usize, i)],
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(x[j] * temp)
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    x[scast(usize, j)],
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // a[j + j * lda] = re(a[j + j * lda])
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            const temp: C1 = ops.mul( // temp = alpha * conj(x[jx])
                                ops.conjugate(x[scast(usize, jx)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;

                            var ix: isize = kx;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // a[i + j * lda] += x[ix] * temp
                                    &a[scast(usize, i + j * lda)],
                                    a[scast(usize, i + j * lda)],
                                    ops.mul(
                                        x[scast(usize, ix)],
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(x[jx] * temp)
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    x[scast(usize, jx)],
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // a[j + j * lda] = re(a[j + j * lda])
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                    }
                }
            } else {
                if (incx == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable) {
                            const temp: C1 = ops.mul( // temp = alpha * x[j]
                                x[scast(usize, j)],
                                alpha,
                                ctx,
                            ) catch unreachable;

                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // a[i + j * lda] += conj(x[i]) * temp
                                    &a[scast(usize, i + j * lda)],
                                    a[scast(usize, i + j * lda)],
                                    ops.mul(
                                        ops.conjugate(x[scast(usize, i)], ctx) catch unreachable,
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(conj(x[j]) * temp)
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    ops.conjugate(x[scast(usize, j)], ctx) catch unreachable,
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // a[j + j * lda] = re(a[j + j * lda])
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            const temp: C1 = ops.mul( // temp = alpha * x[jx]
                                x[scast(usize, jx)],
                                alpha,
                                ctx,
                            ) catch unreachable;

                            var ix: isize = kx;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // a[i + j * lda] += conj(x[ix]) * temp
                                    &a[scast(usize, i + j * lda)],
                                    a[scast(usize, i + j * lda)],
                                    ops.mul(
                                        ops.conjugate(x[scast(usize, ix)], ctx) catch unreachable,
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }
                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(conj(x[jx]) * temp)
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    ops.conjugate(x[scast(usize, jx)], ctx) catch unreachable,
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // a[j + j * lda] = re(a[j + j * lda])
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                    }
                }
            }
        } else {
            if (noconj) {
                if (incx == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable) {
                            const temp: C1 = ops.mul( // temp = alpha * conj(x[j])
                                ops.conjugate(x[scast(usize, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(x[j] * temp)
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    x[scast(usize, j)],
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // a[i + j * lda] += x[i] * temp
                                    &a[scast(usize, i + j * lda)],
                                    a[scast(usize, i + j * lda)],
                                    ops.mul(
                                        x[scast(usize, i)],
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // a[j + j * lda] = re(a[j + j * lda])
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            const temp: C1 = ops.mul( // temp = alpha * conj(x[jx])
                                ops.conjugate(x[scast(usize, jx)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(x[jx] * temp)
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    x[scast(usize, jx)],
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var ix: isize = jx;
                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ix += incx;

                                ops.add_( // a[i + j * lda] += x[ix] * temp
                                    &a[scast(usize, i + j * lda)],
                                    a[scast(usize, i + j * lda)],
                                    ops.mul(
                                        x[scast(usize, ix)],
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // a[j + j * lda] = re(a[j + j * lda])
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                        jx += incx;
                    }
                }
            } else {
                if (incx == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable) {
                            const temp: C1 = ops.mul( // temp = alpha * x[j]
                                x[scast(usize, j)],
                                alpha,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(conj(x[j]) * temp)
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    ops.conjugate(x[scast(usize, j)], ctx) catch unreachable,
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // a[i + j * lda] += conj(x[i]) * temp
                                    &a[scast(usize, i + j * lda)],
                                    a[scast(usize, i + j * lda)],
                                    ops.mul(
                                        ops.conjugate(x[scast(usize, i)], ctx) catch unreachable,
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // a[j + j * lda] = re(a[j + j * lda])
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            const temp: C1 = ops.mul( // temp = alpha * x[jx]
                                x[scast(usize, jx)],
                                alpha,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(conj(x[jx]) * temp)
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    ops.conjugate(x[scast(usize, jx)], ctx) catch unreachable,
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var ix: isize = jx;
                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ix += incx;

                                ops.add_( // a[i + j * lda] += conj(x[ix]) * temp
                                    &a[scast(usize, i + j * lda)],
                                    a[scast(usize, i + j * lda)],
                                    ops.mul(
                                        ops.conjugate(x[scast(usize, ix)], ctx) catch unreachable,
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // a[j + j * lda] = re(a[j + j * lda])
                                &a[scast(usize, j + j * lda)],
                                ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                        jx += incx;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.her not implemented for arbitrary precision types yet");
    }

    return;
}
