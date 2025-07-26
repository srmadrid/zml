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
const Uplo = linalg.Uplo;

pub inline fn syrk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc, ctx);
    } else {
        return k_syrk(uplo.invert(), trans.invert(), n, k, alpha, a, lda, beta, c, ldc, ctx);
    }
}

fn k_syrk(
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const Be: type = @TypeOf(beta);
    const C: type = types.Child(@TypeOf(c));
    const T1: type = types.Coerce(Al, A);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(Be, C)));

    const nrowa: isize = if (trans == .no_trans) n else k;

    if (trans == .conj_no_trans or trans == .conj_trans or
        n < 0 or k < 0 or lda < int.max(1, nrowa) or ldc < int.max(1, n))
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
                        while (i <= j) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(usize, i + j * ldc)],
                                c[scast(usize, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
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
                        var i: isize = j;
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
                        while (i <= j) : (i += 1) {
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
                        if (ops.ne(a[scast(usize, j + l * lda)], 0, ctx) catch unreachable) {
                            const temp: T1 = ops.mul( // temp = alpha * a[j + l * lda]
                                alpha,
                                a[scast(usize, j + l * lda)],
                                ctx,
                            ) catch unreachable;

                            var i: isize = 0;
                            while (i <= j) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                    &c[scast(usize, i + j * ldc)],
                                    c[scast(usize, i + j * ldc)],
                                    ops.mul(
                                        temp,
                                        a[scast(usize, i + l * lda)],
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
                        var i: isize = j;
                        while (i < n) : (i += 1) {
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
                        if (ops.ne(a[scast(usize, j + l * lda)], 0, ctx) catch unreachable) {
                            const temp: T1 = ops.mul( // temp = alpha * a[j + l * lda]
                                alpha,
                                a[scast(usize, j + l * lda)],
                                ctx,
                            ) catch unreachable;

                            var i: isize = j;
                            while (i < n) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                    &c[scast(usize, i + j * ldc)],
                                    c[scast(usize, i + j * ldc)],
                                    ops.mul(
                                        temp,
                                        a[scast(usize, i + l * lda)],
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
            if (uplo == .upper) {
                var j: isize = 0;
                while (j < n) : (j += 1) {
                    var i: isize = 0;
                    while (i <= j) : (i += 1) {
                        var temp: A = constants.zero(A, ctx) catch unreachable;

                        var l: isize = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp += a[l + i * lda] * a[l + j * lda]
                                &temp,
                                temp,
                                ops.mul(
                                    a[scast(usize, l + i * lda)],
                                    a[scast(usize, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.set( // c[i + j * ldc] = alpha * temp
                                &c[scast(usize, i + j * ldc)],
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
                    var i: isize = j;
                    while (i < n) : (i += 1) {
                        var temp: A = constants.zero(A, ctx) catch unreachable;

                        var l: isize = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp += a[l + i * lda] * a[l + j * lda]
                                &temp,
                                temp,
                                ops.mul(
                                    a[scast(usize, l + i * lda)],
                                    a[scast(usize, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.set( // c[i + j * ldc] = alpha * temp
                                &c[scast(usize, i + j * ldc)],
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
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.syrk not implemented for arbitrary precision types yet");
    }

    return;
}
