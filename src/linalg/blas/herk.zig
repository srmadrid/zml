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
const Uplo = types.Uplo;

pub inline fn herk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_herk(
            uplo,
            trans,
            n,
            k,
            alpha,
            a,
            lda,
            beta,
            c,
            ldc,
            ctx,
        );
    } else {
        return k_herk(
            uplo.invert(),
            trans.reverse(),
            n,
            k,
            alpha,
            a,
            lda,
            beta,
            c,
            ldc,
            ctx,
        );
    }
}

fn k_herk(
    uplo: Uplo,
    trans: Transpose,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const Be: type = @TypeOf(beta);
    const C: type = types.Child(@TypeOf(c));
    const T1: type = types.Coerce(Al, A);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(Be, C)));

    const nrowa: i32 = if (trans == .no_trans) n else k;

    if (trans == .conj_no_trans or trans == .trans or
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = 0;
                        while (i <= j) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(u32, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = 0;
                        while (i < j) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(u32, j + j * ldc)],
                            beta,
                            ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = j;
                        while (i < n) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(u32, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(u32, j + j * ldc)],
                            beta,
                            ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var i: i32 = j + 1;
                        while (i < n) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
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
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.eq(beta, 0, ctx) catch unreachable) {
                        var i: i32 = 0;
                        while (i <= j) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(u32, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                        var i: i32 = 0;
                        while (i < j) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(u32, j + j * ldc)],
                            beta,
                            ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.set( // c[j + j * ldc] = re(c[j + j * ldc])
                            &c[scast(u32, j + j * ldc)],
                            ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    var l: i32 = 0;
                    while (l < k) : (l += 1) {
                        if (ops.ne(a[scast(u32, j + l * lda)], 0, ctx) catch unreachable) {
                            const temp: T1 = ops.mul( // temp = alpha * conj(a[j + l * lda])
                                ops.conjugate(a[scast(u32, j + l * lda)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;

                            var i: i32 = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                    &c[scast(u32, i + j * ldc)],
                                    c[scast(u32, i + j * ldc)],
                                    ops.mul(
                                        temp,
                                        a[scast(u32, i + l * lda)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            ops.add_( // c[j + j * ldc] = re(c[j + j * ldc]) + re(temp * a[i + l * lda])
                                &c[scast(u32, j + j * ldc)],
                                ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    temp,
                                    a[scast(u32, i + l * lda)],
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            } else {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.eq(beta, 0, ctx) catch unreachable) {
                        var i: i32 = j;
                        while (i < n) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(u32, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(u32, j + j * ldc)],
                            beta,
                            ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        var i: i32 = j + 1;
                        while (i < n) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    } else {
                        ops.set( // c[j + j * ldc] = re(c[j + j * ldc])
                            &c[scast(u32, j + j * ldc)],
                            ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    var l: i32 = 0;
                    while (l < k) : (l += 1) {
                        if (ops.ne(a[scast(u32, j + l * lda)], 0, ctx) catch unreachable) {
                            const temp: T1 = ops.mul( // temp = alpha * conj(a[j + l * lda])
                                ops.conjugate(a[scast(u32, j + l * lda)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[j + j * ldc] = re(c[j + j * ldc]) + re(temp * a[j + l * lda])
                                &c[scast(u32, j + j * ldc)],
                                ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                                ops.re(ops.mul(
                                    temp,
                                    a[scast(u32, j + l * lda)],
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var i: i32 = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                    &c[scast(u32, i + j * ldc)],
                                    c[scast(u32, i + j * ldc)],
                                    ops.mul(
                                        temp,
                                        a[scast(u32, i + l * lda)],
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
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i < j) : (i += 1) {
                        var temp: A = constants.zero(A, ctx) catch unreachable;

                        var l: i32 = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp += conj(a[l + i * lda]) * a[l + j * lda]
                                &temp,
                                temp,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i * lda)], ctx) catch unreachable,
                                    a[scast(u32, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.mul_( // c[i + j * ldc] = alpha * temp
                                &c[scast(u32, i + j * ldc)],
                                alpha,
                                temp,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += alpha * temp
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    alpha,
                                    temp,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }

                    var rtemp: Scalar(A) = constants.zero(Scalar(A), ctx) catch unreachable;

                    var l: i32 = 0;
                    while (l < k) : (l += 1) {
                        ops.add_( // rtemp += re(conj(a[l + j * lda]) * a[l + j * lda])
                            &rtemp,
                            rtemp,
                            ops.re(ops.mul(
                                ops.conjugate(a[scast(u32, l + j * lda)], ctx) catch unreachable,
                                a[scast(u32, l + j * lda)],
                                ctx,
                            ) catch unreachable, ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    if (ops.eq(beta, 0, ctx) catch unreachable) {
                        ops.mul_( // c[j + j * ldc] = alpha * rtemp
                            &c[scast(u32, j + j * ldc)],
                            alpha,
                            rtemp,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(u32, j + j * ldc)],
                            beta,
                            ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // c[j + j * ldc] += alpha * rtemp
                            &c[scast(u32, j + j * ldc)],
                            c[scast(u32, j + j * ldc)],
                            ops.mul(
                                alpha,
                                rtemp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var rtemp: Scalar(A) = constants.zero(Scalar(A), ctx) catch unreachable;

                    var l: i32 = 0;
                    while (l < k) : (l += 1) {
                        ops.add_( // rtemp += re(conj(a[l + j * lda]) * a[l + j * lda])
                            &rtemp,
                            rtemp,
                            ops.re(ops.mul(
                                ops.conjugate(a[scast(u32, l + j * lda)], ctx) catch unreachable,
                                a[scast(u32, l + j * lda)],
                                ctx,
                            ) catch unreachable, ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    if (ops.eq(beta, 0, ctx) catch unreachable) {
                        ops.mul_( // c[j + j * ldc] = alpha * rtemp
                            &c[scast(u32, j + j * ldc)],
                            alpha,
                            rtemp,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.mul_( // c[j + j * ldc] = beta * re(c[j + j * ldc])
                            &c[scast(u32, j + j * ldc)],
                            beta,
                            ops.re(c[scast(u32, j + j * ldc)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // c[j + j * ldc] += alpha * rtemp
                            &c[scast(u32, j + j * ldc)],
                            c[scast(u32, j + j * ldc)],
                            ops.mul(
                                alpha,
                                rtemp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    var i: i32 = j + 1;
                    while (i < n) : (i += 1) {
                        var temp: A = constants.zero(A, ctx) catch unreachable;

                        l = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp += conj(a[l + i * lda]) * a[l + j * lda]
                                &temp,
                                temp,
                                ops.mul(
                                    ops.conjugate(a[scast(u32, l + i * lda)], ctx) catch unreachable,
                                    a[scast(u32, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.mul_( // c[i + j * ldc] = alpha * temp
                                &c[scast(u32, i + j * ldc)],
                                alpha,
                                temp,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += alpha * temp
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
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
        @compileError("zml.linalg.blas.herk not implemented for arbitrary precision types yet");
    }

    return;
}
