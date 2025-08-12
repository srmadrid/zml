const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Uplo = types.Uplo;
const Diag = types.Diag;
const Order = types.Order;
const Transpose = linalg.Transpose;

pub inline fn tbsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: i32,
    k: i32,
    a: anytype,
    lda: i32,
    x: anytype,
    incx: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_tbsv(
            uplo,
            transa,
            diag,
            n,
            k,
            a,
            lda,
            x,
            incx,
            ctx,
        );
    } else {
        return k_tbsv(
            uplo.invert(),
            transa.invert(),
            diag,
            n,
            k,
            a,
            lda,
            x,
            incx,
            ctx,
        );
    }
}

fn k_tbsv(
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: i32,
    k: i32,
    a: anytype,
    lda: i32,
    x: anytype,
    incx: i32,
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

    var kx: i32 = if (incx < 0) (-n + 1) * incx else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (transa == .no_trans or transa == .conj_no_trans) {
            if (uplo == .upper) {
                if (incx == 1) {
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable) {
                            const l: i32 = k - j;
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= a[k + j * lda]
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        a[scast(u32, k + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, j)];

                                var i: i32 = j - 1;
                                while (i >= int.max(0, j - k)) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * a[l + i + j * lda]
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            a[scast(u32, l + i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(a[k + j * lda])
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ops.conjugate(a[scast(u32, k + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, j)];

                                var i: i32 = j - 1;
                                while (i >= int.max(0, j - k)) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * conj(a[l + i + j * lda])
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
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
                    var jx: i32 = kx;
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        kx -= incx;
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                            var ix: i32 = kx;
                            const l: i32 = k - j;
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= a[k + j * lda]
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        a[scast(u32, k + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, jx)];

                                var i: i32 = j - 1;
                                while (i >= int.max(0, j - k)) : (i -= 1) {
                                    ops.sub_( // x[ix] -= temp * a[l + i + j * lda]
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            a[scast(u32, l + i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix -= incx;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(a[k + j * lda])
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ops.conjugate(a[scast(u32, k + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, jx)];

                                var i: i32 = j - 1;
                                while (i >= int.max(0, j - k)) : (i -= 1) {
                                    ops.sub_( // x[ix] -= temp * conj(a[l + i + j * lda])
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable) {
                            const l: i32 = -j;
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= a[0 + j * lda]
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        a[scast(u32, 0 + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, j)];

                                var i: i32 = j + 1;
                                while (i < int.min(n, j + k + 1)) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * a[l + i + j * lda]
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            a[scast(u32, l + i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(a[0 + j * lda])
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ops.conjugate(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, j)];

                                var i: i32 = j + 1;
                                while (i < int.min(n, j + k + 1)) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * conj(a[l + i + j * lda])
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    var jx: i32 = kx;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        kx += incx;
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                            var ix: i32 = kx;
                            const l: i32 = -j;
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= a[0 + j * lda]
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        a[scast(u32, 0 + j * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, jx)];

                                var i: i32 = j + 1;
                                while (i < int.min(n, j + k + 1)) : (i += 1) {
                                    ops.sub_( // x[ix] -= temp * a[l + i + j * lda]
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            a[scast(u32, l + i + j * lda)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix += incx;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(a[0 + j * lda])
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ops.conjugate(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, jx)];

                                var i: i32 = j + 1;
                                while (i < int.min(n, j + k + 1)) : (i += 1) {
                                    ops.sub_( // x[ix] -= temp * conj(a[l + i + j * lda])
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(u32, j)]);

                        const l: i32 = k - j;
                        if (noconj) {
                            var i: i32 = int.max(0, j - k);
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= a[l + i + j * lda] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(u32, l + i + j * lda)],
                                        x[scast(u32, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[k + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(u32, k + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: i32 = int.max(0, j - k);
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= conj(a[l + i + j * lda]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                        x[scast(u32, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[k + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(u32, k + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, j)] = scast(X, temp);
                    }
                } else {
                    var jx: i32 = kx;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(u32, jx)]);

                        var ix: i32 = kx;
                        const l: i32 = k - j;
                        if (noconj) {
                            var i: i32 = int.max(0, j - k);
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= a[l + i + j * lda] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(u32, l + i + j * lda)],
                                        x[scast(u32, ix)],
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
                                    a[scast(u32, k + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: i32 = int.max(0, j - k);
                            while (i < j) : (i += 1) {
                                ops.sub_( // temp -= conj(a[l + i + j * lda]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                        x[scast(u32, ix)],
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
                                    ops.conjugate(a[scast(u32, k + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, jx)] = scast(X, temp);

                        jx += incx;

                        if (j >= k) {
                            kx += incx;
                        }
                    }
                }
            } else {
                if (incx == 1) {
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(u32, j)]);

                        const l: i32 = -j;
                        if (noconj) {
                            var i: i32 = int.min(n - 1, j + k);
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= a[l + i + j * lda] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(u32, l + i + j * lda)],
                                        x[scast(u32, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= a[0 + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(u32, 0 + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: i32 = int.min(n - 1, j + k);
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= conj(a[l + i + j * lda]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                        x[scast(u32, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(a[0 + j * lda])
                                    &temp,
                                    temp,
                                    ops.conjugate(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, j)] = scast(X, temp);
                    }
                } else {
                    kx += (n - 1) * incx;
                    var jx: i32 = kx;
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(u32, jx)]);

                        var ix: i32 = kx;
                        const l: i32 = -j;
                        if (noconj) {
                            var i: i32 = int.min(n - 1, j + k);
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= a[l + i + j * lda] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        a[scast(u32, l + i + j * lda)],
                                        x[scast(u32, ix)],
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
                                    a[scast(u32, 0 + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: i32 = int.min(n - 1, j + k);
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= conj(a[l + i + j * lda]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                        x[scast(u32, ix)],
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
                                    ops.conjugate(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, jx)] = scast(X, temp);

                        jx -= incx;

                        if (n - 1 - j >= k) {
                            kx -= incx;
                        }
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.tbsv not implemented for arbitrary precision types yet");
    }

    return;
}
