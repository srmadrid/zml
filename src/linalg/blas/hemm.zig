const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;
const Side = linalg.Side;
const Uplo = types.Uplo;

pub inline fn hemm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: i32,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
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
    m: i32,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
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

    const nrowa: i32 = if (side == .left) m else n;

    if (m < 0 or n < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, m) or ldc < int.max(1, m))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or
        (ops.eq(alpha, 0, ctx) catch unreachable and ops.eq(beta, 1, ctx) catch unreachable))
        return;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            if (ops.eq(beta, 0, ctx) catch unreachable) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.set( // c[i + j * ldc] = 0
                            &c[scast(u32, i + j * ldc)],
                            0,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] *= beta
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
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
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        const temp1: T1 = ops.mul( // temp1 = alpha * b[i + j * ldb];
                            alpha,
                            b[scast(u32, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                        var temp2: T2 = constants.zero(T2, ctx) catch unreachable;

                        var k: i32 = 0;
                        while (k < i) : (k += 1) {
                            ops.add_( // c[k + j * ldc] += temp1 * a[k + i * lda];
                                &c[scast(u32, k + j * ldc)],
                                c[scast(u32, k + j * ldc)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, k + i * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += b[k + j * ldb] * conj(a[k + i * lda]);
                                &temp2,
                                temp2,
                                ops.mul(
                                    b[scast(u32, k + j * ldb)],
                                    ops.conj(a[scast(u32, k + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] = temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = m - 1;
                    while (i >= 0) : (i -= 1) {
                        const temp1: T1 = ops.mul( // temp1 = alpha * b[i + j * ldb];
                            alpha,
                            b[scast(u32, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                        var temp2: T2 = constants.zero(T2, ctx) catch unreachable;

                        var k: i32 = i + 1;
                        while (k < m) : (k += 1) {
                            ops.add_( // c[k + j * ldc] += temp1 * a[k + i * lda];
                                &c[scast(u32, k + j * ldc)],
                                c[scast(u32, k + j * ldc)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, k + i * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += b[k + j * ldb] * conj(a[k + i * lda]);
                                &temp2,
                                temp2,
                                ops.mul(
                                    b[scast(u32, k + j * ldb)],
                                    ops.conj(a[scast(u32, k + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] = temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
            var j: i32 = 0;
            while (j < n) : (j += 1) {
                var temp1: T3 = ops.mul( // temp1 = alpha * re(a[j + j * lda]);
                    alpha,
                    ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                    ctx,
                ) catch unreachable;

                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] = temp1 * b[i + j * ldb];
                            &c[scast(u32, i + j * ldc)],
                            temp1,
                            b[scast(u32, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                    }
                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + j * ldb];
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(u32, i + j * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] *= beta;
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            beta,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // c[i + j * ldc] += temp1 * b[i + j * ldb];
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(u32, i + j * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }

                var k: i32 = 0;
                while (k < j) : (k += 1) {
                    if (uplo == .upper) {
                        ops.mul_( // temp1 = alpha * a[k + j * lda];
                            &temp1,
                            alpha,
                            a[scast(u32, k + j * lda)],
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.mul_( // temp1 = alpha * conj(a[j + k * lda]);
                            &temp1,
                            alpha,
                            ops.conj(a[scast(u32, j + k * lda)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + k * ldb];
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(u32, i + k * ldb)],
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
                            ops.conj(a[scast(u32, j + k * lda)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.mul_( // temp1 = alpha * a[k + j * lda];
                            &temp1,
                            alpha,
                            a[scast(u32, k + j * lda)],
                            ctx,
                        ) catch unreachable;
                    }

                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + k * ldb];
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(u32, i + k * ldb)],
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
