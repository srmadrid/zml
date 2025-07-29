const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Uplo = linalg.Uplo;
const Diag = linalg.Diag;
const Order = linalg.Order;
const Transpose = linalg.Transpose;

pub inline fn tpsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: anytype,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_tpsv(uplo, transa, diag, n, ap, x, incx, ctx);
    } else {
        return k_tpsv(uplo.invert(), transa.invert(), diag, n, ap, x, incx, ctx);
    }
}

fn k_tpsv(
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: anytype,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    const A: type = types.Child(@TypeOf(ap));
    const X: type = types.Child(@TypeOf(x));
    const C1: type = types.Coerce(A, X);
    const CC: type = types.Coerce(A, X);

    if (n < 0 or incx == 0)
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
                var kk: isize = int.div(n * (n + 1), 2) - 1;
                if (incx == 1) {
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= ap[kk]
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        ap[scast(usize, kk)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var k: isize = kk - 1;
                                var i: isize = j - 1;
                                while (i >= 0) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * ap[k]
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            ap[scast(usize, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    k -= 1;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(ap[kk])
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        ops.conjugate(ap[scast(usize, kk)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var k: isize = kk - 1;
                                var i: isize = j - 1;
                                while (i >= 0) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * conj(ap[k])
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(ap[scast(usize, k)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    k -= 1;
                                }
                            }
                        }

                        kk -= j + 1;
                    }
                } else {
                    var jx: isize = kx + (n - 1) * incx;
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= ap[kk]
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        ap[scast(usize, kk)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var ix: isize = jx;
                                var k: isize = kk - 1;
                                while (k >= kk - j) : (k -= 1) {
                                    ix -= incx;

                                    ops.sub_( // x[ix] -= temp * ap[k]
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            ap[scast(usize, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(ap[kk])
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        ops.conjugate(ap[scast(usize, kk)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];

                                var ix: isize = jx;
                                var k: isize = kk - 1;
                                while (k >= kk - j) : (k -= 1) {
                                    ix -= incx;

                                    ops.sub_( // x[ix] -= temp * conj(ap[k])
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(ap[scast(usize, k)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        jx -= incx;
                        kk -= j + 1;
                    }
                }
            } else {
                var kk: isize = 0;
                if (incx == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, j)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= ap[kk]
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        ap[scast(usize, kk)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var k: isize = kk + 1;
                                var i: isize = j + 1;
                                while (i < n) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * ap[k]
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            ap[scast(usize, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    k += 1;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(ap[kk])
                                        &x[scast(usize, j)],
                                        x[scast(usize, j)],
                                        ops.conjugate(ap[scast(usize, kk)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, j)];

                                var k: isize = kk + 1;
                                var i: isize = j + 1;
                                while (i < n) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * conj(ap[k])
                                        &x[scast(usize, i)],
                                        x[scast(usize, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(ap[scast(usize, k)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    k += 1;
                                }
                            }
                        }

                        kk += n - j;
                    }
                } else {
                    var jx: isize = kx;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(usize, jx)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= ap[kk]
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        ap[scast(usize, kk)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];
                                var ix: isize = jx;
                                var k: isize = kk + 1;
                                while (k < kk + n - j) : (k += 1) {
                                    ix += incx;

                                    ops.sub_( // x[ix] -= temp * ap[k]
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            ap[scast(usize, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(ap[kk])
                                        &x[scast(usize, jx)],
                                        x[scast(usize, jx)],
                                        ops.conjugate(ap[scast(usize, kk)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(usize, jx)];
                                var ix: isize = jx;
                                var k: isize = kk + 1;
                                while (k < kk + n - j) : (k += 1) {
                                    ix += incx;

                                    ops.sub_( // x[ix] -= temp * conj(ap[k])
                                        &x[scast(usize, ix)],
                                        x[scast(usize, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(ap[scast(usize, k)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        jx += incx;
                        kk += n - j;
                    }
                }
            }
        } else {
            if (uplo == .upper) {
                var kk: isize = 0;
                if (incx == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(usize, j)]);

                        if (noconj) {
                            var k: isize = kk;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // temp -= ap[k] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(usize, k)],
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }

                            if (nounit) {
                                ops.div_( // temp /= ap[kk + j]
                                    &temp,
                                    temp,
                                    ap[scast(usize, kk + j)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var k: isize = kk;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // temp -= conj(ap[k]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(ap[scast(usize, k)], ctx) catch unreachable,
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(ap[kk + j])
                                    &temp,
                                    temp,
                                    ops.conjugate(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, j)] = scast(X, temp);

                        kk += j + 1;
                    }
                } else {
                    var jx: isize = kx;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(usize, jx)]);

                        var ix: isize = kx;
                        if (noconj) {
                            var k: isize = kk;
                            while (k < kk + j) : (k += 1) {
                                ops.sub_( // temp -= ap[k] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(usize, k)],
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= ap[kk + j]
                                    &temp,
                                    temp,
                                    ap[scast(usize, kk + j)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var k: isize = kk;
                            while (k < kk + j) : (k += 1) {
                                ops.sub_( // temp -= conj(ap[k]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(ap[scast(usize, k)], ctx) catch unreachable,
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(ap[kk + j])
                                    &temp,
                                    temp,
                                    ops.conjugate(ap[scast(usize, kk + j)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, jx)] = scast(X, temp);

                        jx += incx;
                        kk += j + 1;
                    }
                }
            } else {
                var kk: isize = int.div(n * (n + 1), 2) - 1;
                if (incx == 1) {
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(usize, j)]);

                        var k: isize = kk;
                        if (noconj) {
                            var i: isize = n - 1;
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= ap[k] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(usize, k)],
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k -= 1;
                            }

                            if (nounit) {
                                ops.div_( // temp /= ap[kk - (n - 1) + j]
                                    &temp,
                                    temp,
                                    ap[scast(usize, kk - (n - 1) + j)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: isize = n - 1;
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= conj(ap[k]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(ap[scast(usize, k)], ctx) catch unreachable,
                                        x[scast(usize, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k -= 1;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(ap[kk - (n - 1) + j])
                                    &temp,
                                    temp,
                                    ops.conjugate(ap[scast(usize, kk - (n - 1) + j)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, j)] = scast(X, temp);

                        kk -= n - j;
                    }
                } else {
                    kx += (n - 1) * incx;
                    var jx: isize = kx;
                    var j: isize = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(usize, jx)]);

                        if (noconj) {
                            var ix: isize = kx;
                            var k: isize = kk;
                            while (k > kk - (n - (j + 1))) : (k -= 1) {
                                ops.sub_( // temp -= ap[k] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(usize, k)],
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix -= incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= ap[kk - (n - 1) + j]
                                    &temp,
                                    temp,
                                    ap[scast(usize, kk - (n - 1) + j)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var ix: isize = kx;
                            var k: isize = kk;
                            while (k > kk - (n - (j + 1))) : (k -= 1) {
                                ops.sub_( // temp -= conj(ap[k] * x[ix])
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(ap[scast(usize, k)], ctx) catch unreachable,
                                        x[scast(usize, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix -= incx;
                            }

                            if (nounit) {
                                ops.div_( // temp /= conj(ap[kk - (n - 1) + j])
                                    &temp,
                                    temp,
                                    ops.conjugate(ap[scast(usize, kk - (n - 1) + j)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(usize, jx)] = scast(X, temp);

                        jx -= incx;
                        kk -= n - j;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.tpsv not implemented for arbitrary precision types yet");
    }

    return;
}
