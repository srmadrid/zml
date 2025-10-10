const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;
const Transpose = linalg.Transpose;
const Side = linalg.Side;
const Uplo = types.Uplo;
const Diag = types.Diag;

pub inline fn trsm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: i32,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_trsm(
            side,
            uplo,
            transa,
            diag,
            m,
            n,
            alpha,
            a,
            lda,
            b,
            ldb,
            ctx,
        );
    } else {
        return k_trsm(
            side.invert(),
            uplo.invert(),
            transa,
            diag,
            n,
            m,
            alpha,
            a,
            lda,
            b,
            ldb,
            ctx,
        );
    }
}

fn k_trsm(
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: i32,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const T1: type = types.Coerce(Al, types.Coerce(A, B));
    const CC: type = types.Coerce(Al, types.Coerce(A, B));

    const not: bool = transa == .no_trans or transa == .conj_no_trans;
    const noconj: bool = transa == .no_trans or transa == .trans;

    const nrowa: i32 = if (side == .left) m else n;

    if (m < 0 or n < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, m))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0)
        return;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            var j: i32 = 0;
            while (j < n) : (j += 1) {
                var i: i32 = 0;
                while (i < m) : (i += 1) {
                    ops.set( // b[i + j * ldb] = 0
                        &b[scast(u32, i + j * ldb)],
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (noconj) {
                            if (ops.ne(alpha, 1, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.mul_( // b[i + j * ldb] *= alpha
                                        &b[scast(u32, i + j * ldb)],
                                        b[scast(u32, i + j * ldb)],
                                        alpha,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            var k: i32 = m - 1;
                            while (k >= 0) : (k -= 1) {
                                if (ops.ne(b[scast(u32, k + j * ldb)], 0, ctx) catch unreachable) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= a[k + k * lda]
                                            &b[scast(u32, k + j * ldb)],
                                            b[scast(u32, k + j * lda)],
                                            a[scast(u32, k + k * lda)],
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: i32 = 0;
                                    while (i < k) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                b[scast(u32, k + j * ldb)],
                                                a[scast(u32, i + k * lda)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                }
                            }
                        } else {
                            if (ops.ne(alpha, 1, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.mul_( // b[i + j * ldb] *= alpha
                                        &b[scast(u32, i + j * ldb)],
                                        b[scast(u32, i + j * ldb)],
                                        alpha,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            var k: i32 = m - 1;
                            while (k >= 0) : (k -= 1) {
                                if (ops.ne(b[scast(u32, k + j * ldb)], 0, ctx) catch unreachable) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= conj(a[k + k * lda])
                                            &b[scast(u32, k + j * ldb)],
                                            b[scast(u32, k + j * lda)],
                                            ops.conj(a[scast(u32, k + k * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: i32 = 0;
                                    while (i < k) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * conj(a[i + k * lda])
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                b[scast(u32, k + j * ldb)],
                                                ops.conj(a[scast(u32, i + k * lda)], ctx) catch unreachable,
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(u32, i + j * ldb)],
                                    alpha,
                                    b[scast(u32, i + j * ldb)],
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: i32 = 0;
                        while (k < m) : (k += 1) {
                            if (ops.ne(b[scast(u32, k + j * ldb)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= a[k + k * lda]
                                            &b[scast(u32, k + j * ldb)],
                                            b[scast(u32, k + j * lda)],
                                            a[scast(u32, k + k * lda)],
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: i32 = k + 1;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                b[scast(u32, k + j * ldb)],
                                                a[scast(u32, i + k * lda)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= conj(a[k + k * lda])
                                            &b[scast(u32, k + j * ldb)],
                                            b[scast(u32, k + j * lda)],
                                            ops.conj(a[scast(u32, k + k * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: i32 = k + 1;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * conj(a[i + k * lda])
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                b[scast(u32, k + j * ldb)],
                                                ops.conj(a[scast(u32, i + k * lda)], ctx) catch unreachable,
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = 0;
                        while (i < m) : (i += 1) {
                            var temp: T1 = scast(T1, ops.mul( // temp = alpha * b[i + j * ldb]
                                b[scast(u32, i + j * ldb)],
                                alpha,
                                ctx,
                            ) catch unreachable);

                            if (noconj) {
                                var k: i32 = 0;
                                while (k < i) : (k += 1) {
                                    ops.sub_( // temp -= a[k + i * lda] * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(u32, k + i * lda)],
                                            b[scast(u32, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= a[i + i * lda]
                                        &temp,
                                        temp,
                                        a[scast(u32, i + i * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var k: i32 = 0;
                                while (k < i) : (k += 1) {
                                    ops.sub_( // temp -= conj(a[k + i * lda]) * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conj(a[scast(u32, k + i * lda)], ctx) catch unreachable,
                                            b[scast(u32, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= conj(a[i + i * lda])
                                        &temp,
                                        temp,
                                        ops.conj(a[scast(u32, i + i * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            ops.set( // b[i + j * ldb] = temp
                                &b[scast(u32, i + j * ldb)],
                                temp,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = m - 1;
                        while (i >= 0) : (i -= 1) {
                            var temp: T1 = scast(T1, ops.mul( // temp = alpha * b[i + j * ldb]
                                alpha,
                                b[scast(u32, i + j * ldb)],
                                ctx,
                            ) catch unreachable);

                            if (noconj) {
                                var k: i32 = i + 1;
                                while (k < m) : (k += 1) {
                                    ops.sub_( // temp -= a[k + i * lda] * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(u32, k + i * lda)],
                                            b[scast(u32, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= a[i + i * lda]
                                        &temp,
                                        temp,
                                        a[scast(u32, i + i * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var k: i32 = i + 1;
                                while (k < m) : (k += 1) {
                                    ops.sub_( // temp -= conj(a[k + i * lda]) * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conj(a[scast(u32, k + i * lda)], ctx) catch unreachable,
                                            b[scast(u32, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= conj(a[i + i * lda])
                                        &temp,
                                        temp,
                                        ops.conj(a[scast(u32, i + i * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            ops.set( // b[i + j * ldb] = temp
                                &b[scast(u32, i + j * ldb)],
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(u32, i + j * ldb)],
                                    b[scast(u32, i + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: i32 = 0;
                        while (k < j) : (k += 1) {
                            if (ops.ne(a[scast(u32, k + j * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    var i: i32 = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                a[scast(u32, k + j * lda)],
                                                b[scast(u32, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    var i: i32 = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= conj(a[k + j * lda]) * b[i + k * ldb]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                ops.conj(a[scast(u32, k + j * lda)], ctx) catch unreachable,
                                                b[scast(u32, i + k * ldb)],
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
                                    a[scast(u32, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conj(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= temp
                                    &b[scast(u32, i + j * ldb)],
                                    b[scast(u32, i + j * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                } else {
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(u32, i + j * ldb)],
                                    b[scast(u32, i + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: i32 = j + 1;
                        while (k < n) : (k += 1) {
                            if (ops.ne(a[scast(u32, k + j * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    var i: i32 = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                a[scast(u32, k + j * lda)],
                                                b[scast(u32, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    var i: i32 = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= conj(a[k + j * lda]) * b[i + k * ldb]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                ops.conj(a[scast(u32, k + j * lda)], ctx) catch unreachable,
                                                b[scast(u32, i + k * ldb)],
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
                                    a[scast(u32, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conj(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= temp
                                    &b[scast(u32, i + j * ldb)],
                                    b[scast(u32, i + j * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                if (uplo == .upper) {
                    var k: i32 = n - 1;
                    while (k >= 0) : (k -= 1) {
                        var temp: A = constants.one(A, ctx) catch unreachable;

                        if (diag == .non_unit) {
                            if (noconj) {
                                ops.div_( // temp = a[k + k * lda]
                                    &temp,
                                    temp,
                                    a[scast(u32, k + k * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp = conj(a[k + k * lda])
                                    &temp,
                                    temp,
                                    ops.conj(a[scast(u32, k + k * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= temp
                                    &b[scast(u32, i + k * ldb)],
                                    b[scast(u32, i + k * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var j: i32 = 0;
                        while (j < k) : (j += 1) {
                            if (ops.ne(a[scast(u32, j + k * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    ops.set( // temp = a[j + k * lda]
                                        &temp,
                                        a[scast(u32, j + k * lda)],
                                        ctx,
                                    ) catch unreachable;
                                } else {
                                    ops.set( // temp = conj(a[j + k * lda])
                                        &temp,
                                        ops.conj(a[scast(u32, j + k * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.sub_( // b[i + j * ldb] -= temp * b[i + k * ldb]
                                        &b[scast(u32, i + j * ldb)],
                                        b[scast(u32, i + j * ldb)],
                                        ops.mul(
                                            temp,
                                            b[scast(u32, i + k * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= alpha
                                    &b[scast(u32, i + k * ldb)],
                                    b[scast(u32, i + k * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                } else {
                    var k: i32 = 0;
                    while (k < n) : (k += 1) {
                        var temp: A = constants.one(A, ctx) catch unreachable;

                        if (diag == .non_unit) {
                            if (noconj) {
                                ops.div_( // temp = a[k + k * lda]
                                    &temp,
                                    temp,
                                    a[scast(u32, k + k * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp = conj(a[k + k * lda])
                                    &temp,
                                    temp,
                                    ops.conj(a[scast(u32, k + k * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= temp
                                    &b[scast(u32, i + k * ldb)],
                                    b[scast(u32, i + k * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var j: i32 = k + 1;
                        while (j < n) : (j += 1) {
                            if (ops.ne(a[scast(u32, j + k * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    ops.set( // temp = a[j + k * lda]
                                        &temp,
                                        a[scast(u32, j + k * lda)],
                                        ctx,
                                    ) catch unreachable;
                                } else {
                                    ops.set( // temp = conj(a[j + k * lda])
                                        &temp,
                                        ops.conj(a[scast(u32, j + k * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.sub_( // b[i + j * ldb] -= temp * b[i + k * ldb]
                                        &b[scast(u32, i + j * ldb)],
                                        b[scast(u32, i + j * ldb)],
                                        ops.mul(
                                            temp,
                                            b[scast(u32, i + k * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= alpha
                                    &b[scast(u32, i + k * ldb)],
                                    b[scast(u32, i + k * ldb)],
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
