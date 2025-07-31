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

pub inline fn gemm(
    order: Order,
    transa: Transpose,
    transb: Transpose,
    m: isize,
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
        return k_gemm(
            transa,
            transb,
            m,
            n,
            k,
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
        return k_gemm(
            transb,
            transa,
            n,
            m,
            k,
            alpha,
            b,
            ldb,
            a,
            lda,
            beta,
            c,
            ldc,
            ctx,
        );
    }
}

fn k_gemm(
    transa: Transpose,
    transb: Transpose,
    m: isize,
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
    const T2: type = types.Coerce(A, B);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(B, types.Coerce(Be, C))));

    const nota: bool = transa == .no_trans or transa == .conj_no_trans;
    const notb: bool = transb == .no_trans or transb == .conj_no_trans;
    const noconja: bool = transa == .no_trans or transa == .trans;
    const noconjb: bool = transb == .no_trans or transb == .trans;

    const nrowa: isize = if (nota) m else k;
    const nrowb: isize = if (notb) k else n;

    if (m < 0 or n < 0 or k < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, nrowb) or ldc < int.max(1, m))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or
        ((ops.eq(alpha, 0, ctx) catch unreachable or k == 0) and
            ops.eq(beta, 1, ctx) catch unreachable))
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
            } else {
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

        if (notb) {
            if (nota) {
                if (noconjb) {
                    if (noconja) {
                        var j: isize = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(usize, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: isize = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * b[l + j * ldb]
                                    b[scast(usize, l + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            a[scast(usize, i + l * lda)],
                                            temp,
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
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(usize, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: isize = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * b[l + j * ldb]
                                    b[scast(usize, l + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * conj(a[i + l * lda])
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, i + l * lda)], ctx) catch unreachable,
                                            temp,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    if (noconja) {
                        var j: isize = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(usize, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: isize = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * conj(b[l + j * ldb])
                                    ops.conjugate(b[scast(usize, l + j * ldb)], ctx) catch unreachable,
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            a[scast(usize, i + l * lda)],
                                            temp,
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
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(usize, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: isize = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * conj(b[l + j * ldb])
                                    ops.conjugate(b[scast(usize, l + j * ldb)], ctx) catch unreachable,
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * conj(a[i + l * lda])
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, i + l * lda)], ctx) catch unreachable,
                                            temp,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                }
            } else {
                if (noconjb) {
                    if (noconja) {
                        var j: isize = 0;
                        while (j < n) : (j += 1) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: isize = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += a[l + i * lda] * b[l + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(usize, l + i * lda)],
                                            b[scast(usize, l + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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

                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: isize = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += conj(a[l + i * lda]) * b[l + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, l + i * lda)], ctx) catch unreachable,
                                            b[scast(usize, l + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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

                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    if (noconja) {
                        var j: isize = 0;
                        while (j < n) : (j += 1) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: isize = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += a[l + i * lda] * conj(b[l + j * ldb])
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(usize, l + i * lda)],
                                            ops.conjugate(b[scast(usize, l + j * ldb)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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

                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: isize = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += conj(a[l + i * lda]) * conj(b[l + j * ldb])
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, l + i * lda)], ctx) catch unreachable,
                                            ops.conjugate(b[scast(usize, l + j * ldb)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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

                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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
            if (nota) {
                if (noconjb) {
                    if (noconja) {
                        var j: isize = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(usize, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: isize = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * b[j + l * ldb]
                                    b[scast(usize, j + l * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            a[scast(usize, i + l * lda)],
                                            temp,
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
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(usize, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: isize = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * b[j + l * ldb]
                                    b[scast(usize, j + l * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * conj(a[i + l * lda])
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, i + l * lda)], ctx) catch unreachable,
                                            temp,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    if (noconja) {
                        var j: isize = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(usize, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: isize = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * conj(b[j + l * ldb])
                                    ops.conjugate(b[scast(usize, j + l * ldb)], ctx) catch unreachable,
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            a[scast(usize, i + l * lda)],
                                            temp,
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
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(usize, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: isize = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * conj(b[j + l * ldb])
                                    ops.conjugate(b[scast(usize, j + l * ldb)], ctx) catch unreachable,
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * conj(a[i + l * lda])
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, i + l * lda)], ctx) catch unreachable,
                                            temp,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                }
            } else {
                if (noconjb) {
                    if (noconja) {
                        var j: isize = 0;
                        while (j < n) : (j += 1) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: isize = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += a[l + i * lda] * b[j + l * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(usize, l + i * lda)],
                                            b[scast(usize, j + l * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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

                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: isize = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += conj(a[l + i * lda]) * b[j + l * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, l + i * lda)], ctx) catch unreachable,
                                            b[scast(usize, j + l * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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

                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    if (noconja) {
                        var j: isize = 0;
                        while (j < n) : (j += 1) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: isize = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += a[l + i * lda] * conj(b[j + l * ldb])
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(usize, l + i * lda)],
                                            ops.conjugate(b[scast(usize, j + l * ldb)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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

                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: isize = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += conj(a[l + i * lda]) * conj(b[j + l * ldb])
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, l + i * lda)], ctx) catch unreachable,
                                            ops.conjugate(b[scast(usize, j + l * ldb)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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

                                    ops.add_( // c[i + j * ldc] += alpha * temp
                                        &c[scast(usize, i + j * ldc)],
                                        c[scast(usize, i + j * ldc)],
                                        ops.mul(
                                            alpha,
                                            temp,
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
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.gemm not implemented for arbitrary precision types yet");
    }

    return;
}
