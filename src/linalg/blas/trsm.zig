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
const Side = linalg.Side;
const Uplo = linalg.Uplo;
const Diag = linalg.Diag;

pub inline fn trsm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_trsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, ctx);
    } else {
        return k_trsm(side.invert(), uplo.invert(), transa, diag, n, m, alpha, a, lda, b, ldb, ctx);
    }
}

fn k_trsm(
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const T1: type = types.Coerce(Al, types.Coerce(A, B));
    const CC: type = types.Coerce(Al, types.Coerce(A, B));

    const not: bool = transa == .no_trans or transa == .conj_no_trans;
    const noconj: bool = transa == .no_trans or transa == .trans;

    const nrowa: isize = if (side == .left) m else n;

    if (m < 0 or n < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, m))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0)
        return;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            var j: isize = 0;
            while (j < n) : (j += 1) {
                var i: isize = 0;
                while (i < m) : (i += 1) {
                    ops.set( // b[i + j * ldb] = 0
                        &b[scast(usize, i + j * ldb)],
                        0,
                        ctx,
                    ) catch unreachable;
                }
            }

            return;
        }

        if (side == .left) {
            if (not) {
                if (uplo == .upper) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (noconj) {
                            if (ops.ne(alpha, 1, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.mul_( // b[i + j * ldb] *= alpha
                                        &b[scast(usize, i + j * ldb)],
                                        b[scast(usize, i + j * ldb)],
                                        alpha,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            var k: isize = m - 1;
                            while (k >= 0) : (k -= 1) {
                                if (ops.ne(b[scast(usize, k + j * ldb)], 0, ctx) catch unreachable) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= a[k + k * lda]
                                            &b[scast(usize, k + j * ldb)],
                                            b[scast(usize, k + j * lda)],
                                            a[scast(usize, k + k * lda)],
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: isize = 0;
                                    while (i < k) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda]
                                            &b[scast(usize, i + j * ldb)],
                                            b[scast(usize, i + j * ldb)],
                                            ops.mul(
                                                b[scast(usize, k + j * ldb)],
                                                a[scast(usize, i + k * lda)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                }
                            }
                        } else {
                            if (ops.ne(alpha, 1, ctx) catch unreachable) {
                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.mul_( // b[i + j * ldb] *= alpha
                                        &b[scast(usize, i + j * ldb)],
                                        b[scast(usize, i + j * ldb)],
                                        alpha,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            var k: isize = m - 1;
                            while (k >= 0) : (k -= 1) {
                                if (ops.ne(b[scast(usize, k + j * ldb)], 0, ctx) catch unreachable) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= conj(a[k + k * lda])
                                            &b[scast(usize, k + j * ldb)],
                                            b[scast(usize, k + j * lda)],
                                            ops.conjugate(a[scast(usize, k + k * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: isize = 0;
                                    while (i < k) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * conj(a[i + k * lda])
                                            &b[scast(usize, i + j * ldb)],
                                            b[scast(usize, i + j * ldb)],
                                            ops.mul(
                                                b[scast(usize, k + j * ldb)],
                                                ops.conjugate(a[scast(usize, i + k * lda)], ctx) catch unreachable,
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
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(usize, i + j * ldb)],
                                    b[scast(usize, i + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: isize = 0;
                        while (k < m) : (k += 1) {
                            if (ops.ne(b[scast(usize, k + j * ldb)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= a[k + k * lda]
                                            &b[scast(usize, k + j * ldb)],
                                            b[scast(usize, k + j * lda)],
                                            a[scast(usize, k + k * lda)],
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: isize = k + 1;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda]
                                            &b[scast(usize, i + j * ldb)],
                                            b[scast(usize, i + j * ldb)],
                                            ops.mul(
                                                b[scast(usize, k + j * ldb)],
                                                a[scast(usize, i + k * lda)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= conj(a[k + k * lda])
                                            &b[scast(usize, k + j * ldb)],
                                            b[scast(usize, k + j * lda)],
                                            ops.conjugate(a[scast(usize, k + k * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: isize = k + 1;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * conj(a[i + k * lda])
                                            &b[scast(usize, i + j * ldb)],
                                            b[scast(usize, i + j * ldb)],
                                            ops.mul(
                                                b[scast(usize, k + j * ldb)],
                                                ops.conjugate(a[scast(usize, i + k * lda)], ctx) catch unreachable,
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
                if (uplo == .upper) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var i: isize = 0;
                        while (i < m) : (i += 1) {
                            var temp: T1 = scast(T1, ops.mul( // temp = alpha * b[i + j * ldb]
                                b[scast(usize, i + j * ldb)],
                                alpha,
                                ctx,
                            ) catch unreachable);

                            if (noconj) {
                                var k: isize = 0;
                                while (k < i) : (k += 1) {
                                    ops.sub_( // temp -= a[k + i * lda] * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(usize, k + i * lda)],
                                            b[scast(usize, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= a[i + i * lda]
                                        &temp,
                                        temp,
                                        a[scast(usize, i + i * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var k: isize = 0;
                                while (k < i) : (k += 1) {
                                    ops.sub_( // temp -= conj(a[k + i * lda]) * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, k + i * lda)], ctx) catch unreachable,
                                            b[scast(usize, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= conj(a[i + i * lda])
                                        &temp,
                                        temp,
                                        ops.conjugate(a[scast(usize, i + i * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            ops.set( // b[i + j * ldb] = temp
                                &b[scast(usize, i + j * ldb)],
                                temp,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var i: isize = m - 1;
                        while (i >= 0) : (i -= 1) {
                            var temp: T1 = scast(T1, ops.mul( // temp = alpha * b[i + j * ldb]
                                alpha,
                                b[scast(usize, i + j * ldb)],
                                ctx,
                            ) catch unreachable);

                            if (noconj) {
                                var k: isize = i + 1;
                                while (k < m) : (k += 1) {
                                    ops.sub_( // temp -= a[k + i * lda] * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(usize, k + i * lda)],
                                            b[scast(usize, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= a[i + i * lda]
                                        &temp,
                                        temp,
                                        a[scast(usize, i + i * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var k: isize = i + 1;
                                while (k < m) : (k += 1) {
                                    ops.sub_( // temp -= conj(a[k + i * lda]) * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conjugate(a[scast(usize, k + i * lda)], ctx) catch unreachable,
                                            b[scast(usize, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= conj(a[i + i * lda])
                                        &temp,
                                        temp,
                                        ops.conjugate(a[scast(usize, i + i * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            ops.set( // b[i + j * ldb] = temp
                                &b[scast(usize, i + j * ldb)],
                                temp,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            }
        } else {
            if (not) {
                if (uplo == .upper) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(usize, i + j * ldb)],
                                    b[scast(usize, i + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: isize = 0;
                        while (k < j) : (k += 1) {
                            if (ops.ne(a[scast(usize, k + j * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    var i: isize = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb]
                                            &b[scast(usize, i + j * ldb)],
                                            b[scast(usize, i + j * ldb)],
                                            ops.mul(
                                                a[scast(usize, k + j * lda)],
                                                b[scast(usize, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    var i: isize = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= conj(a[k + j * lda]) * b[i + k * ldb]
                                            &b[scast(usize, i + j * ldb)],
                                            b[scast(usize, i + j * ldb)],
                                            ops.mul(
                                                ops.conjugate(a[scast(usize, k + j * lda)], ctx) catch unreachable,
                                                b[scast(usize, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                }
                            }
                        }

                        if (diag == .non_unit) {
                            var temp: A = constants.one(A, ctx) catch unreachable;

                            if (noconj) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= temp
                                    &b[scast(usize, i + j * ldb)],
                                    b[scast(usize, i + j * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                } else {
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(usize, i + j * ldb)],
                                    b[scast(usize, i + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: isize = j + 1;
                        while (k < n) : (k += 1) {
                            if (ops.ne(a[scast(usize, k + j * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    var i: isize = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb]
                                            &b[scast(usize, i + j * ldb)],
                                            b[scast(usize, i + j * ldb)],
                                            ops.mul(
                                                a[scast(usize, k + j * lda)],
                                                b[scast(usize, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    var i: isize = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= conj(a[k + j * lda]) * b[i + k * ldb]
                                            &b[scast(usize, i + j * ldb)],
                                            b[scast(usize, i + j * ldb)],
                                            ops.mul(
                                                ops.conjugate(a[scast(usize, k + j * lda)], ctx) catch unreachable,
                                                b[scast(usize, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                }
                            }
                        }

                        if (diag == .non_unit) {
                            var temp: A = constants.one(A, ctx) catch unreachable;

                            if (noconj) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= temp
                                    &b[scast(usize, i + j * ldb)],
                                    b[scast(usize, i + j * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                if (uplo == .upper) {
                    var k: isize = n - 1;
                    while (k >= 0) : (k -= 1) {
                        var temp: A = constants.one(A, ctx) catch unreachable;

                        if (diag == .non_unit) {
                            if (noconj) {
                                ops.div_( // temp = a[k + k * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, k + k * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp = conj(a[k + k * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, k + k * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= temp
                                    &b[scast(usize, i + k * ldb)],
                                    b[scast(usize, i + k * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var j: isize = 0;
                        while (j < k) : (j += 1) {
                            if (ops.ne(a[scast(usize, j + k * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    ops.set( // temp = a[j + k * lda]
                                        &temp,
                                        a[scast(usize, j + k * lda)],
                                        ctx,
                                    ) catch unreachable;
                                } else {
                                    ops.set( // temp = conj(a[j + k * lda])
                                        &temp,
                                        ops.conjugate(a[scast(usize, j + k * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.sub_( // b[i + j * ldb] -= temp * b[i + k * ldb]
                                        &b[scast(usize, i + j * ldb)],
                                        b[scast(usize, i + j * ldb)],
                                        ops.mul(
                                            temp,
                                            b[scast(usize, i + k * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= alpha
                                    &b[scast(usize, i + k * ldb)],
                                    b[scast(usize, i + k * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                } else {
                    var k: isize = 0;
                    while (k < n) : (k += 1) {
                        var temp: A = constants.one(A, ctx) catch unreachable;

                        if (diag == .non_unit) {
                            if (noconj) {
                                ops.div_( // temp = a[k + k * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, k + k * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp = conj(a[k + k * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, k + k * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= temp
                                    &b[scast(usize, i + k * ldb)],
                                    b[scast(usize, i + k * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var j: isize = k + 1;
                        while (j < n) : (j += 1) {
                            if (ops.ne(a[scast(usize, j + k * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    ops.set( // temp = a[j + k * lda]
                                        &temp,
                                        a[scast(usize, j + k * lda)],
                                        ctx,
                                    ) catch unreachable;
                                } else {
                                    ops.set( // temp = conj(a[j + k * lda])
                                        &temp,
                                        ops.conjugate(a[scast(usize, j + k * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                var i: isize = 0;
                                while (i < m) : (i += 1) {
                                    ops.sub_( // b[i + j * ldb] -= temp * b[i + k * ldb]
                                        &b[scast(usize, i + j * ldb)],
                                        b[scast(usize, i + j * ldb)],
                                        ops.mul(
                                            temp,
                                            b[scast(usize, i + k * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: isize = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= alpha
                                    &b[scast(usize, i + k * ldb)],
                                    b[scast(usize, i + k * ldb)],
                                    alpha,
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
        @compileError("zml.linalg.blas.trsm not implemented for arbitrary precision types yet");
    }

    return;
}
