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

pub inline fn hpr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    ap: anytype,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_hpr(uplo, n, alpha, x, incx, ap, true, ctx);
    } else {
        return k_hpr(uplo.invert(), n, alpha, x, incx, ap, false, ctx);
    }
}

pub fn k_hpr(
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    ap: anytype,
    noconj: bool,
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

    const kx: isize = if (incx < 0) (-n + 1) * incx else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        var kk: isize = 0;
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

                            var k: isize = kk;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // ap[k] += x[i] * temp
                                    &ap[scast(usize, k)],
                                    ap[scast(usize, k)],
                                    ops.mul(
                                        x[scast(usize, i)],
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }

                            ops.add_( // ap[kk + j] = re(ap[kk + j]) + re(x[j] * temp)
                                &ap[scast(usize, kk + j)],
                                ops.re(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    x[scast(usize, j)],
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // ap[kk + j] = re(ap[kk + j])
                                &ap[scast(usize, kk + j)],
                                ops.re(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        kk += j + 1;
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
                            var k: isize = kk;
                            while (k < kk + j) : (k += 1) {
                                ops.add_( // ap[k] += x[ix] * temp
                                    &ap[scast(usize, k)],
                                    ap[scast(usize, k)],
                                    ops.mul(
                                        x[scast(usize, ix)],
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            ops.add_( // ap[kk + j] = re(ap[kk + j]) + re(x[jx] * temp)
                                &ap[scast(usize, kk + j)],
                                ops.re(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    x[scast(usize, jx)],
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // ap[kk + j] = re(ap[kk + j])
                                &ap[scast(usize, kk + j)],
                                ops.re(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        kk += j + 1;
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

                            var k: isize = kk;
                            std.debug.print("incx: {}\n", .{incx});
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // ap[k] += conj(x[i]) * temp
                                    &ap[scast(usize, k)],
                                    ap[scast(usize, k)],
                                    ops.mul(
                                        ops.conjugate(x[scast(usize, i)], ctx) catch unreachable,
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }

                            ops.add_( // ap[kk + j] = re(ap[kk + j]) + re(conj(x[j]) * temp)
                                &ap[scast(usize, kk + j)],
                                ops.re(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    ops.conjugate(x[scast(usize, j)], ctx) catch unreachable,
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // ap[kk + j] = re(ap[kk + j])
                                &ap[scast(usize, kk + j)],
                                ops.re(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        kk += j + 1;
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
                            var k: isize = kk;
                            while (k < kk + j) : (k += 1) {
                                ops.add_( // ap[k] += conj(x[ix]) * temp
                                    &ap[scast(usize, k)],
                                    ap[scast(usize, k)],
                                    ops.mul(
                                        ops.conjugate(x[scast(usize, ix)], ctx) catch unreachable,
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            ops.add_( // ap[kk + j] = re(ap[kk + j]) + re(conj(x[jx]) * temp)
                                &ap[scast(usize, kk + j)],
                                ops.re(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    ops.conjugate(x[scast(usize, jx)], ctx) catch unreachable,
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // ap[kk + j] = re(ap[kk + j])
                                &ap[scast(usize, kk + j)],
                                ops.re(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        kk += j + 1;
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

                            ops.add_( // ap[kk] = re(ap[kk]) + re(x[j] * temp)
                                &ap[scast(usize, kk)],
                                ops.re(ap[scast(usize, kk)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    x[scast(usize, j)],
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var k: isize = kk + 1;
                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // ap[k] += x[i] * temp
                                    &ap[scast(usize, k)],
                                    ap[scast(usize, k)],
                                    ops.mul(
                                        x[scast(usize, i)],
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }
                        } else {
                            ops.set( // ap[kk] = re(ap[kk])
                                &ap[scast(usize, kk)],
                                ops.re(ap[scast(usize, kk)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        kk += n - j;
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

                            ops.add_( // ap[kk] = re(ap[kk]) + re(x[jx] * temp)
                                &ap[scast(usize, kk)],
                                ops.re(ap[scast(usize, kk)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    x[scast(usize, jx)],
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var ix: isize = jx;
                            var k: isize = kk + 1;
                            while (k < kk + n - j) : (k += 1) {
                                ix += incx;

                                ops.add_( // ap[k] += x[ix] * temp
                                    &ap[scast(usize, k)],
                                    ap[scast(usize, k)],
                                    ops.mul(
                                        x[scast(usize, ix)],
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // ap[kk] = re(ap[kk])
                                &ap[scast(usize, kk)],
                                ops.re(ap[scast(usize, kk)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        kk += n - j;
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

                            ops.add_( // ap[kk] = re(ap[kk]) + re(conj(x[j]) * temp)
                                &ap[scast(usize, kk)],
                                ops.re(ap[scast(usize, kk)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    ops.conjugate(x[scast(usize, j)], ctx) catch unreachable,
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var k: isize = kk + 1;
                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // ap[k] += conj(x[i]) * temp
                                    &ap[scast(usize, k)],
                                    ap[scast(usize, k)],
                                    ops.mul(
                                        ops.conjugate(x[scast(usize, i)], ctx) catch unreachable,
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }
                        } else {
                            ops.set( // ap[kk] = re(ap[kk])
                                &ap[scast(usize, kk)],
                                ops.re(ap[scast(usize, kk)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        kk += n - j;
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

                            ops.add_( // ap[kk] = re(ap[kk]) + re(conj(x[jx]) * temp)
                                &ap[scast(usize, kk)],
                                ops.re(ap[scast(usize, kk)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    ops.conjugate(x[scast(usize, jx)], ctx) catch unreachable,
                                    temp,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var ix: isize = jx;
                            var k: isize = kk + 1;
                            while (k < kk + n - j) : (k += 1) {
                                ix += incx;

                                ops.add_( // ap[k] += conj(x[ix]) * temp
                                    &ap[scast(usize, k)],
                                    ap[scast(usize, k)],
                                    ops.mul(
                                        ops.conjugate(x[scast(usize, ix)], ctx) catch unreachable,
                                        temp,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // ap[kk] = re(ap[kk])
                                &ap[scast(usize, kk)],
                                ops.re(ap[scast(usize, kk)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        kk += n - j;
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
