const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = linalg.Order;
const Side = linalg.Side;
const Uplo = linalg.Uplo;

pub inline fn hemm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_hemm(
            side,
            uplo,
            m,
            n,
            alpha,
            a,
            lda,
            b,
            ldb,
            beta,
            c,
            ldc,
            ctx,
        );
    } else {
        return k_hemm(
            side.invert(),
            uplo.invert(),
            n,
            m,
            alpha,
            a,
            lda,
            b,
            ldb,
            beta,
            c,
            ldc,
            ctx,
        );
    }
}

fn k_hemm(
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const Be: type = @TypeOf(beta);
    const C: type = types.Child(@TypeOf(c));
    const T1: type = types.Coerce(Al, B);
    const T2: type = types.Coerce(A, B);
    const T3: type = types.Coerce(Al, A);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(B, types.Coerce(Be, C))));

    const nrowa: isize = if (side == .left) m else n;

    if (m < 0 or n < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, m) or ldc < int.max(1, m))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or
        (ops.eq(alpha, 0, ctx) catch unreachable and ops.eq(beta, 1, ctx) catch unreachable))
        return;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            if (ops.eq(beta, 0, ctx) catch unreachable) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    var i: isize = 0;
                    while (i < m) : (i += 1) {
                        ops.set( // c[i + j * ldc] = 0
                            &c[scast(usize, i + j * ldc)],
                            0,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    var i: isize = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] *= beta
                            &c[scast(usize, i + j * ldc)],
                            c[scast(usize, i + j * ldc)],
                            beta,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }

            return;
        }

        if (side == .left) {
            if (uplo == .upper) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    var i: isize = 0;
                    while (i < m) : (i += 1) {
                        const temp1: T1 = ops.mul( // temp1 = alpha * b[i + j * ldb];
                            alpha,
                            b[scast(usize, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                        var temp2: T2 = constants.zero(T2, ctx) catch unreachable;

                        var k: isize = 0;
                        while (k < i) : (k += 1) {
                            ops.add_( // c[k + j * ldc] += temp1 * a[k + i * lda];
                                &c[scast(usize, k + j * ldc)],
                                c[scast(usize, k + j * ldc)],
                                ops.mul(
                                    temp1,
                                    a[scast(usize, k + i * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += b[k + j * ldb] * conj(a[k + i * lda]);
                                &temp2,
                                temp2,
                                ops.mul(
                                    b[scast(usize, k + j * ldb)],
                                    ops.conjugate(a[scast(usize, k + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] = temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(usize, i + j * ldc)],
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(usize, i + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ops.mul(
                                    alpha,
                                    temp2,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(usize, i + i * lda)], ctx) catch unreachable,
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
                        } else {
                            ops.mul_( // c[i + j * ldc] *= beta;
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(usize, i + i * lda)], ctx) catch unreachable,
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
                    }
                }
            } else {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    var i: isize = m - 1;
                    while (i >= 0) : (i -= 1) {
                        const temp1: T1 = ops.mul( // temp1 = alpha * b[i + j * ldb];
                            alpha,
                            b[scast(usize, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                        var temp2: T2 = constants.zero(T2, ctx) catch unreachable;

                        var k: isize = i + 1;
                        while (k < m) : (k += 1) {
                            ops.add_( // c[k + j * ldc] += temp1 * a[k + i * lda];
                                &c[scast(usize, k + j * ldc)],
                                c[scast(usize, k + j * ldc)],
                                ops.mul(
                                    temp1,
                                    a[scast(usize, k + i * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += b[k + j * ldb] * conj(a[k + i * lda]);
                                &temp2,
                                temp2,
                                ops.mul(
                                    b[scast(usize, k + j * ldb)],
                                    ops.conjugate(a[scast(usize, k + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] = temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(usize, i + j * ldc)],
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(usize, i + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ops.mul(
                                    alpha,
                                    temp2,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(usize, i + i * lda)], ctx) catch unreachable,
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
                        } else {
                            ops.mul_( // c[i + j * ldc] *= beta;
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(usize, i + i * lda)], ctx) catch unreachable,
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
                    }
                }
            }
        } else {
            var j: isize = 0;
            while (j < n) : (j += 1) {
                var temp1: T3 = ops.mul( // temp1 = alpha * re(a[j + j * lda]);
                    alpha,
                    ops.re(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                    ctx,
                ) catch unreachable;

                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    var i: isize = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] = temp1 * b[i + j * ldb];
                            &c[scast(usize, i + j * ldc)],
                            temp1,
                            b[scast(usize, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                    }
                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                    var i: isize = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + j * ldb];
                            &c[scast(usize, i + j * ldc)],
                            c[scast(usize, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(usize, i + j * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    var i: isize = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] *= beta;
                            &c[scast(usize, i + j * ldc)],
                            c[scast(usize, i + j * ldc)],
                            beta,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // c[i + j * ldc] += temp1 * b[i + j * ldb];
                            &c[scast(usize, i + j * ldc)],
                            c[scast(usize, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(usize, i + j * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }

                var k: isize = 0;
                while (k < j) : (k += 1) {
                    if (uplo == .upper) {
                        ops.mul_( // temp1 = alpha * a[k + j * lda];
                            &temp1,
                            alpha,
                            a[scast(usize, k + j * lda)],
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.mul_( // temp1 = alpha * conj(a[j + k * lda]);
                            &temp1,
                            alpha,
                            ops.conjugate(a[scast(usize, j + k * lda)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    var i: isize = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + k * ldb];
                            &c[scast(usize, i + j * ldc)],
                            c[scast(usize, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(usize, i + k * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }

                k = j + 1;
                while (k < n) : (k += 1) {
                    if (uplo == .upper) {
                        ops.mul_( // temp1 = alpha * conj(a[j + k * lda]);
                            &temp1,
                            alpha,
                            ops.conjugate(a[scast(usize, j + k * lda)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.mul_( // temp1 = alpha * a[k + j * lda];
                            &temp1,
                            alpha,
                            a[scast(usize, k + j * lda)],
                            ctx,
                        ) catch unreachable;
                    }

                    var i: isize = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + k * ldb];
                            &c[scast(usize, i + j * ldc)],
                            c[scast(usize, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(usize, i + k * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.hemm not implemented for arbitrary precision types yet");
    }

    return;
}
