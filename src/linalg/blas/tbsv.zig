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

pub inline fn tbsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_tbsv(uplo, transa, diag, n, k, a, lda, x, incx, ctx);
    } else {
        return k_tbsv(uplo.invert(), transa.invert(), diag, n, k, a, lda, x, incx, ctx);
    }
}

fn k_tbsv(
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
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

    if (n < 0 or k < 0 or lda < (k + 1) or incx == 0)
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
                            const l: isize = k - j;
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= a[k + j * lda]
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        a[scast(usize, k + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var i: isize = j - 1;
                                while (i >= int.max(0, j - k)) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * a[l + i + j * lda]
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            a[scast(usize, l + i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(a[k + j * lda])
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        ops.conjugate(a[scast(usize, k + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var i: isize = j - 1;
                                while (i >= int.max(0, j - k)) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * conj(a[l + i + j * lda])
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(usize, l + i + j * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    kx += (n - 1) * incx;
                    var jx: isize = kx;
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        kx -= incx;
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            var ix: isize = kx;
                            const l: isize = k - j;
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= a[k + j * lda]
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        a[scast(usize, k + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var i: isize = j - 1;
                                while (i >= int.max(0, j - k)) : (i -= 1) {
                                    ops.sub_( // x[ix] -= temp * a[l + i + j * lda]
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            a[scast(usize, l + i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix -= incx;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(a[k + j * lda])
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        ops.conjugate(a[scast(usize, k + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var i: isize = j - 1;
                                while (i >= int.max(0, j - k)) : (i -= 1) {
                                    ops.sub_( // x[ix] -= temp * conj(a[l + i + j * lda])
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(usize, l + i + j * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix -= incx;
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
                            const l: isize = -j;
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= a[0 + j * lda]
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        a[scast(usize, 0 + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var i: isize = j + 1;
                                while (i < int.min(n, j + k + 1)) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * a[l + i + j * lda]
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            a[scast(usize, l + i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(a[0 + j * lda])
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        ops.conjugate(a[scast(usize, 0 + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var i: isize = j + 1;
                                while (i < int.min(n, j + k + 1)) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * conj(a[l + i + j * lda])
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(usize, l + i + j * lda)], ctx) catch unreachable,
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
                        kx += incx;
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            var ix: isize = kx;
                            const l: isize = -j;
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= a[0 + j * lda]
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        a[scast(usize, 0 + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var i: isize = j + 1;
                                while (i < int.min(n, j + k + 1)) : (i += 1) {
                                    ops.sub_( // x[ix] -= temp * a[l + i + j * lda]
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            a[scast(usize, l + i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix += incx;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(a[0 + j * lda])
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        ops.conjugate(a[scast(usize, 0 + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var i: isize = j + 1;
                                while (i < int.min(n, j + k + 1)) : (i += 1) {
                                    ops.sub_( // x[ix] -= temp * conj(a[l + i + j * lda])
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(usize, l + i + j * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix += incx;
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

                        const l: isize = k - j;
                        if (noconj) {
                            var i: isize = int.max(0, j - k);
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= a[l + i + j * lda] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, l + i + j * lda)],
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[k + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, k + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: isize = int.max(0, j - k);
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= conj(a[l + i + j * lda]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(usize, l + i + j * lda)], ctx) catch unreachable,
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[k + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, k + j * lda)], ctx) catch unreachable,
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

                        var ix: isize = kx;
                        const l: isize = k - j;
                        if (noconj) {
                            var i: isize = int.max(0, j - k);
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= a[l + i + j * lda] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, l + i + j * lda)],
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[k + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, k + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: isize = int.max(0, j - k);
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= conj(a[l + i + j * lda]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(usize, l + i + j * lda)], ctx) catch unreachable,
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[k + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, k + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, jx)] = scast(X, temp);

                        jx += incx;

                        if (j >= k) {
                            kx += incx;
                        }
                    }
                }
            } else {
                if (incx == 1) {
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(usize, j)]);

                        const l: isize = -j;
                        if (noconj) {
                            var i: isize = int.min(n - 1, j + k);
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= a[l + i + j * lda] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, l + i + j * lda)],
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[0 + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, 0 + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: isize = int.min(n - 1, j + k);
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= conj(a[l + i + j * lda]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(usize, l + i + j * lda)], ctx) catch unreachable,
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[0 + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, 0 + j * lda)], ctx) catch unreachable,
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

                        var ix: isize = kx;
                        const l: isize = -j;
                        if (noconj) {
                            var i: isize = int.min(n - 1, j + k);
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= a[l + i + j * lda] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(usize, l + i + j * lda)],
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix -= incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[0 + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(usize, 0 + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: isize = int.min(n - 1, j + k);
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= conj(a[l + i + j * lda]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(usize, l + i + j * lda)], ctx) catch unreachable,
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix -= incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[0 + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(usize, 0 + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, jx)] = scast(X, temp);

                        jx -= incx;

                        if ((n - 1 - j) >= k) {
                            kx -= incx;
                        }
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
