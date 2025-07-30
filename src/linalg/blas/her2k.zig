const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = linalg.Order;
const Transpose = linalg.Transpose;
const Uplo = linalg.Uplo;

pub inline fn her2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
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
        return k_her2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
    } else {
        return k_her2k(uplo.invert(), trans.reverse(), n, k, ops.conjugate(alpha, ctx) catch unreachable, b, ldb, a, lda, beta, c, ldc, ctx);
    }
}

fn k_her2k(
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
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
    const T2: type = types.Coerce(Al, A);
    const T3: type = types.Coerce(A, B);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(B, types.Coerce(Be, C))));

    const nrowa: isize = if (trans == .no_trans) n else k;

    if (trans == .conj_no_trans or trans == .trans or
        n < 0 or k < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, nrowa) or ldc < int.max(1, n))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or
        ((ops.eq(alpha, 0, ctx) catch unreachable or k == 0) and
            ops.eq(beta, 1, ctx) catch unreachable))
        return;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            if (uplo == .upper) {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var i: isize = 0;
                        while (i <= j) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(usize, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var i: isize = 0;
                        while (i < j) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(usize, j + j * ldc)],
                            beta,
                            ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var i: isize = j;
                        while (i < n) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(usize, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(usize, j + j * ldc)],
                            beta,
                            ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var i: isize = j + 1;
                        while (i < n) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            }

            return;
        }

        if (trans == .no_trans) {
            if (uplo == .upper) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    if (ops.eq(beta, 0, ctx) catch unreachable) {
                        var i: isize = 0;
                        while (i <= j) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(usize, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                        var i: isize = 0;
                        while (i < j) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(usize, j + j * ldc)],
                            beta,
                            ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.set( // c[j + j * ldc] = re(c[j + j * ldc])
                            &c[scast(usize, j + j * ldc)],
                            ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    var l: isize = 0;
                    while (l < k) : (l += 1) {
                        if (ops.ne(a[scast(usize, j + l * lda)], 0, ctx) catch unreachable or
                            ops.ne(b[scast(usize, j + l * ldb)], 0, ctx) catch unreachable)
                        {
                            const temp1: T1 = ops.mul( // temp1 = alpha * conj(b[j + l * ldb])
                                alpha,
                                ops.conjugate(b[scast(usize, j + l * ldb)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                            const temp2: T2 = ops.conjugate(ops.mul( // temp2 = conj(alpha * a[j + l * lda])
                                alpha,
                                a[scast(usize, j + l * lda)],
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += a[i + l * lda] * temp1 + b[i + l * ldb] * temp2
                                    &c[scast(usize, i + j * ldc)],
                                    c[scast(usize, i + j * ldc)],
                                    ops.add(
                                        ops.mul(
                                            a[scast(usize, i + l * lda)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            b[scast(usize, i + l * ldb)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            ops.add_( // c[j + j * ldc] = re(c[j + j * ldc]) + re(a[j + l * lda] * temp1 + b[j + l * ldb] * temp2)
                                &c[scast(usize, j + j * ldc)],
                                ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        a[scast(usize, j + l * lda)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        b[scast(usize, j + l * ldb)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            } else {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    if (ops.eq(beta, 0, ctx) catch unreachable) {
                        var i: isize = j;
                        while (i < n) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(usize, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                        var i: isize = j + 1;
                        while (i < n) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(usize, j + j * ldc)],
                            beta,
                            ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.set( // c[j + j * ldc] = re(c[j + j * ldc])
                            &c[scast(usize, j + j * ldc)],
                            ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    var l: isize = 0;
                    while (l < k) : (l += 1) {
                        if (ops.ne(a[scast(usize, j + l * lda)], 0, ctx) catch unreachable or
                            ops.ne(b[scast(usize, j + l * ldb)], 0, ctx) catch unreachable)
                        {
                            const temp1: T1 = ops.mul( // temp1 = alpha * conj(b[j + l * ldb])
                                alpha,
                                ops.conjugate(b[scast(usize, j + l * ldb)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                            const temp2: T2 = ops.conjugate(ops.mul( // temp2 = conj(alpha * a[j + l * lda])
                                alpha,
                                a[scast(usize, j + l * lda)],
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += a[i + l * lda] * temp1 + b[i + l * ldb] * temp2
                                    &c[scast(usize, i + j * ldc)],
                                    c[scast(usize, i + j * ldc)],
                                    ops.add(
                                        ops.mul(
                                            a[scast(usize, i + l * lda)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            b[scast(usize, i + l * ldb)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            ops.add_( // c[j + j * ldc] = re(c[j + j * ldc]) + re(a[j + l * lda] * temp1 + b[j + l * ldb] * temp2)
                                &c[scast(usize, j + j * ldc)],
                                ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        a[scast(usize, j + l * lda)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        b[scast(usize, j + l * ldb)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            }
        } else {
            if (uplo == .upper) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    var i: isize = 0;
                    while (i <= j) : (i += 1) {
                        var temp1: T3 = constants.zero(T3, ctx) catch unreachable;
                        var temp2: T3 = constants.zero(T3, ctx) catch unreachable;

                        var l: isize = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp1 += conj(a[l + i * lda]) * b[l + j * ldb]
                                &temp1,
                                temp1,
                                ops.mul(
                                    ops.conjugate(a[scast(usize, l + i * lda)], ctx) catch unreachable,
                                    b[scast(usize, l + j * ldb)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(b[l + i * ldb]) * a[l + j * lda]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(b[scast(usize, l + i * ldb)], ctx) catch unreachable,
                                    a[scast(usize, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (i == j) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                ops.set( // c[j + j * ldc] = re(alpha * temp1 + conj(alpha) * temp2)
                                    &c[scast(usize, j + j * ldc)],
                                    ops.re(ops.add(
                                        ops.mul(
                                            alpha,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conjugate(alpha, ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable, ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                                    &c[scast(usize, j + j * ldc)],
                                    ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                                    beta,
                                    ctx,
                                ) catch unreachable;

                                ops.add_( // c[j + j * ldc] += re(alpha * temp1 + conj(alpha) * temp2)
                                    &c[scast(usize, j + j * ldc)],
                                    c[scast(usize, j + j * ldc)],
                                    ops.re(ops.add(
                                        ops.mul(
                                            alpha,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conjugate(alpha, ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable, ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                ops.set( // c[i + j * ldc] = alpha * temp1 + conj(alpha) * temp2
                                    &c[scast(usize, i + j * ldc)],
                                    ops.add(ops.mul(
                                        alpha,
                                        temp1,
                                        ctx,
                                    ) catch unreachable, ops.mul(
                                        ops.conjugate(alpha, ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable, ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.mul_( // c[i + j * ldc] *= beta
                                    &c[scast(usize, i + j * ldc)],
                                    c[scast(usize, i + j * ldc)],
                                    beta,
                                    ctx,
                                ) catch unreachable;

                                ops.add_( // c[i + j * ldc] += alpha * temp1 + conj(alpha) * temp2
                                    &c[scast(usize, i + j * ldc)],
                                    c[scast(usize, i + j * ldc)],
                                    ops.add(
                                        ops.mul(
                                            alpha,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conjugate(alpha, ctx) catch unreachable,
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
                    var i: isize = j;
                    while (i < n) : (i += 1) {
                        var temp1: T3 = constants.zero(T3, ctx) catch unreachable;
                        var temp2: T3 = constants.zero(T3, ctx) catch unreachable;

                        var l: isize = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp1 += conj(a[l + i * lda]) * b[l + j * ldb]
                                &temp1,
                                temp1,
                                ops.mul(
                                    ops.conjugate(a[scast(usize, l + i * lda)], ctx) catch unreachable,
                                    b[scast(usize, l + j * ldb)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(b[l + i * ldb]) * a[l + j * lda]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conjugate(b[scast(usize, l + i * ldb)], ctx) catch unreachable,
                                    a[scast(usize, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (i == j) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                ops.set( // c[j + j * ldc] = re(alpha * temp1 + conj(alpha) * temp2)
                                    &c[scast(usize, j + j * ldc)],
                                    ops.re(ops.add(
                                        ops.mul(
                                            alpha,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conjugate(alpha, ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable, ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                                    &c[scast(usize, j + j * ldc)],
                                    ops.re(c[scast(usize, j + j * ldc)], ctx) catch unreachable,
                                    beta,
                                    ctx,
                                ) catch unreachable;

                                ops.add_( // c[j + j * ldc] += re(alpha * temp1 + conj(alpha) * temp2)
                                    &c[scast(usize, j + j * ldc)],
                                    c[scast(usize, j + j * ldc)],
                                    ops.re(ops.add(
                                        ops.mul(
                                            alpha,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conjugate(alpha, ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable, ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                ops.set( // c[i + j * ldc] = alpha * temp1 + conj(alpha) * temp2
                                    &c[scast(usize, i + j * ldc)],
                                    ops.add(
                                        ops.mul(
                                            alpha,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conjugate(alpha, ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.mul_( // c[i + j * ldc] *= beta
                                    &c[scast(usize, i + j * ldc)],
                                    c[scast(usize, i + j * ldc)],
                                    beta,
                                    ctx,
                                ) catch unreachable;

                                ops.add_( // c[i + j * ldc] += alpha * temp1 + conj(alpha) * temp2
                                    &c[scast(usize, i + j * ldc)],
                                    c[scast(usize, i + j * ldc)],
                                    ops.add(
                                        ops.mul(
                                            alpha,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conjugate(alpha, ctx) catch unreachable,
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
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.her2k not implemented for arbitrary precision types yet");
    }

    return;
}
