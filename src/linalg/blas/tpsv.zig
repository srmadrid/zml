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

pub inline fn tpsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: i32,
    ap: anytype,
    x: anytype,
    incx: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_tpsv(
            uplo,
            transa,
            diag,
            n,
            ap,
            x,
            incx,
            ctx,
        );
    } else {
        return k_tpsv(
            uplo.invert(),
            transa.invert(),
            diag,
            n,
            ap,
            x,
            incx,
            ctx,
        );
    }
}

fn k_tpsv(
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: i32,
    ap: anytype,
    x: anytype,
    incx: i32,
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

    var kx: i32 = if (incx < 0) (-n + 1) * incx else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (transa == .no_trans or transa == .conj_no_trans) {
            if (uplo == .upper) {
                var kk: i32 = int.div(n * (n + 1), 2) - 1;
                if (incx == 1) {
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= ap[kk]
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ap[scast(u32, kk)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, j)];

                                var k: i32 = kk - 1;
                                var i: i32 = j - 1;
                                while (i >= 0) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * ap[k]
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            ap[scast(u32, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    k -= 1;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(ap[kk])
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ops.conj(ap[scast(u32, kk)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, j)];

                                var k: i32 = kk - 1;
                                var i: i32 = j - 1;
                                while (i >= 0) : (i -= 1) {
                                    ops.sub_( // x[i] -= temp * conj(ap[k])
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            ops.conj(ap[scast(u32, k)], ctx) catch unreachable,
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
                    var jx: i32 = kx + (n - 1) * incx;
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= ap[kk]
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ap[scast(u32, kk)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, jx)];

                                var ix: i32 = jx;
                                var k: i32 = kk - 1;
                                while (k >= kk - j) : (k -= 1) {
                                    ix -= incx;

                                    ops.sub_( // x[ix] -= temp * ap[k]
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ap[scast(u32, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(ap[kk])
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ops.conj(ap[scast(u32, kk)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, jx)];

                                var ix: i32 = jx;
                                var k: i32 = kk - 1;
                                while (k >= kk - j) : (k -= 1) {
                                    ix -= incx;

                                    ops.sub_( // x[ix] -= temp * conj(ap[k])
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conj(ap[scast(u32, k)], ctx) catch unreachable,
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
                var kk: i32 = 0;
                if (incx == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[j] /= ap[kk]
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ap[scast(u32, kk)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, j)];

                                var k: i32 = kk + 1;
                                var i: i32 = j + 1;
                                while (i < n) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * ap[k]
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            ap[scast(u32, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    k += 1;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[j] /= conj(ap[kk])
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ops.conj(ap[scast(u32, kk)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, j)];

                                var k: i32 = kk + 1;
                                var i: i32 = j + 1;
                                while (i < n) : (i += 1) {
                                    ops.sub_( // x[i] -= temp * conj(ap[k])
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            ops.conj(ap[scast(u32, k)], ctx) catch unreachable,
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
                    var jx: i32 = kx;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                            if (noconj) {
                                if (nounit) {
                                    ops.div_( // x[jx] /= ap[kk]
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ap[scast(u32, kk)],
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, jx)];
                                var ix: i32 = jx;
                                var k: i32 = kk + 1;
                                while (k < kk + n - j) : (k += 1) {
                                    ix += incx;

                                    ops.sub_( // x[ix] -= temp * ap[k]
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ap[scast(u32, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                if (nounit) {
                                    ops.div_( // x[jx] /= conj(ap[kk])
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ops.conj(ap[scast(u32, kk)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                const temp: X = x[scast(u32, jx)];
                                var ix: i32 = jx;
                                var k: i32 = kk + 1;
                                while (k < kk + n - j) : (k += 1) {
                                    ix += incx;

                                    ops.sub_( // x[ix] -= temp * conj(ap[k])
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conj(ap[scast(u32, k)], ctx) catch unreachable,
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
                var kk: i32 = 0;
                if (incx == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(u32, j)]);

                        if (noconj) {
                            var k: i32 = kk;
                            var i: i32 = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // temp -= ap[k] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(u32, k)],
                                        x[scast(u32, i)],
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
                                    ap[scast(u32, kk + j)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var k: i32 = kk;
                            var i: i32 = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // temp -= conj(ap[k]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conj(ap[scast(u32, k)], ctx) catch unreachable,
                                        x[scast(u32, i)],
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
                                    ops.conj(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, j)] = scast(X, temp);

                        kk += j + 1;
                    }
                } else {
                    var jx: i32 = kx;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(u32, jx)]);

                        var ix: i32 = kx;
                        if (noconj) {
                            var k: i32 = kk;
                            while (k < kk + j) : (k += 1) {
                                ops.sub_( // temp -= ap[k] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(u32, k)],
                                        x[scast(u32, ix)],
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
                                    ap[scast(u32, kk + j)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var k: i32 = kk;
                            while (k < kk + j) : (k += 1) {
                                ops.sub_( // temp -= conj(ap[k]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conj(ap[scast(u32, k)], ctx) catch unreachable,
                                        x[scast(u32, ix)],
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
                                    ops.conj(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, jx)] = scast(X, temp);

                        jx += incx;
                        kk += j + 1;
                    }
                }
            } else {
                var kk: i32 = int.div(n * (n + 1), 2) - 1;
                if (incx == 1) {
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(u32, j)]);

                        var k: i32 = kk;
                        if (noconj) {
                            var i: i32 = n - 1;
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= ap[k] * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(u32, k)],
                                        x[scast(u32, i)],
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
                                    ap[scast(u32, kk - (n - 1) + j)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var i: i32 = n - 1;
                            while (i > j) : (i -= 1) {
                                ops.sub_( // temp -= conj(ap[k]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conj(ap[scast(u32, k)], ctx) catch unreachable,
                                        x[scast(u32, i)],
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
                                    ops.conj(ap[scast(u32, kk - (n - 1) + j)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, j)] = scast(X, temp);

                        kk -= n - j;
                    }
                } else {
                    kx += (n - 1) * incx;
                    var jx: i32 = kx;
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(u32, jx)]);

                        if (noconj) {
                            var ix: i32 = kx;
                            var k: i32 = kk;
                            while (k > kk - (n - (j + 1))) : (k -= 1) {
                                ops.sub_( // temp -= ap[k] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(u32, k)],
                                        x[scast(u32, ix)],
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
                                    ap[scast(u32, kk - (n - 1) + j)],
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            var ix: i32 = kx;
                            var k: i32 = kk;
                            while (k > kk - (n - (j + 1))) : (k -= 1) {
                                ops.sub_( // temp -= conj(ap[k] * x[ix])
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conj(ap[scast(u32, k)], ctx) catch unreachable,
                                        x[scast(u32, ix)],
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
                                    ops.conj(ap[scast(u32, kk - (n - 1) + j)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, jx)] = scast(X, temp);

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
