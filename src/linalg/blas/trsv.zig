const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Uplo = linalg.Uplo;
const Diag = linalg.Diag;
const Order = linalg.Order;
const Transpose = linalg.Transpose;

pub inline fn trsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_trsv(uplo, transa, diag, n, a, lda, x, incx, ctx);
    } else {
        return k_trsv(uplo.invert(), transa.invert(), diag, n, a, lda, x, incx, ctx);
    }
}

fn k_trsv(
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    const A: type = types.Child(@TypeOf(a));
    const X: type = types.Child(@TypeOf(x));
    const C1: type = types.Coerce(A, X);
    const CC: type = types.Coerce(A, X);

    if (n < 0 or lda < int.max(1, n) or incx == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0)
        return;

    const noconj: bool = transa == .no_trans or transa == .trans;
    const nounit: bool = diag == .non_unit;

    var kx: isize = if (incx < 0) (-n + 1) * incx else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (transa == .no_trans or transa == .conj_no_trans) {
            if (uplo == .upper) {
                if (incx == 1) {
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= a[j + j * lda]
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        a[scast(usize, j + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var i: isize = j - 1;
                                while (i >= 0) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * a[i + j * lda]
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            a[scast(usize, i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(a[j + j * lda])
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var i: isize = j - 1;
                                while (i >= 0) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * conj(a[i + j * lda])
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(usize, i + j * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    var jx: isize = kx + (n - 1) * incx;
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= a[j + j * lda]
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        a[scast(usize, j + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var ix: isize = jx;
                                var i: isize = j - 1;
                                while (i >= 0) : (i -= 1) {
                                    ix -= incx;

                                    ops.sub_( // x[ix] -= temp * a[i + j * lda]
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            a[scast(usize, i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(a[j + j * lda])
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var ix: isize = jx;
                                var i: isize = j - 1;
                                while (i >= 0) : (i -= 1) {
                                    ix -= incx;

                                    ops.sub_( // x[ix] -= temp * conj(a[i + j * lda])
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(usize, i + j * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        jx -= incx;
                    }
                }
            } else {
                if (incx == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= a[j + j * lda]
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        a[scast(usize, j + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var i: isize = j + 1;
                                while (i < n) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * a[i + j * lda]
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            a[scast(usize, i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(a[j + j * lda])
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var i: isize = j + 1;
                                while (i < n) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * conj(a[i + j * lda])
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(usize, i + j * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= a[j + j * lda]
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        a[scast(usize, j + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var ix: isize = jx;
                                var i: isize = j + 1;
                                while (i < n) : (i += 1) {
                                    ix += incx;

                                    ops.sub_( // x[ix] -= temp * a[i + j * lda]
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            a[scast(usize, i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(a[j + j * lda])
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var ix: isize = jx;
                                var i: isize = j + 1;
                                while (i < n) : (i += 1) {
                                    ix += incx;

                                    ops.sub_( // x[ix] -= temp * conj(a[i + j * lda])
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(usize, i + j * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        jx += incx;
                    }
                }
            }
        } else {
            if (uplo == .upper) {
                if (incx == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(usize, j)]);

                        if (noconj) {
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= a[i + j * lda] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, i + j * lda)],
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= conj(a[i + j * lda]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(usize, i + j * lda)], ctx) catch unreachable,
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, j)] = scast(X, temp);
                    }
                } else {
                    var jx: isize = kx;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(usize, jx)]);

                        if (noconj) {
                            var ix: isize = kx;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= a[i + j * lda] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, i + j * lda)],
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var ix: isize = kx;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= conj(a[i + j * lda]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(usize, i + j * lda)], ctx) catch unreachable,
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, jx)] = scast(X, temp);

                        jx += incx;
                    }
                }
            } else {
                if (incx == 1) {
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(usize, j)]);

                        if (noconj) {
                            var i: isize = n - 1;
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= a[i + j * lda] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, i + j * lda)],
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: isize = n - 1;
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= a[i + j * lda] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, i + j * lda)],
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, j)] = scast(X, temp);
                    }
                } else {
                    kx += (n - 1) * incx;
                    var jx: isize = kx;
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(usize, jx)]);

                        if (noconj) {
                            var ix: isize = kx;
                            var i: isize = n - 1;
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= a[i + j * lda] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, i + j * lda)],
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix -= incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var ix: isize = kx;
                            var i: isize = n - 1;
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= conj(a[i + j * lda]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(usize, i + j * lda)], ctx) catch unreachable,
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix -= incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, jx)] = scast(X, temp);
                        jx -= incx;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.gbmv not implemented for arbitrary precision types yet");
    }

    return;
}
