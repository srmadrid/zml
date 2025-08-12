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

pub inline fn tpmv(
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
        return k_tpmv(
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
        return k_tpmv(
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

fn k_tpmv(
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
                var kk: i32 = 0;
                if (incx == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable) {
                            const temp: X = x[scast(u32, j)];
                            if (noconj) {
                                var k: i32 = kk;
                                var i: i32 = 0;
                                while (i < j) : (i += 1) {
                                    ops.add_( // x[i] += temp * ap[k]
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

                                if (nounit) {
                                    ops.mul_( // x[j] *= ap[kk + j]
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ap[scast(u32, kk + j)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var k: i32 = kk;
                                var i: i32 = 0;
                                while (i < j) : (i += 1) {
                                    ops.add_( // x[i] += temp * conj(ap[k])
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(ap[scast(u32, k)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    k += 1;
                                }

                                if (nounit) {
                                    ops.mul_( // x[j] *= conj(ap[kk + j])
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ops.conjugate(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        kk += j + 1;
                    }
                } else {
                    var jx: i32 = kx;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                            const temp: X = x[scast(u32, jx)];

                            if (noconj) {
                                var ix: i32 = kx;
                                var k: i32 = kk;
                                while (k < kk + j) : (k += 1) {
                                    ops.add_( // x[ix] += temp * ap[k]
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ap[scast(u32, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix += incx;
                                }

                                if (nounit) {
                                    ops.mul_( // x[jx] *= ap[kk + j]
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ap[scast(u32, kk + j)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var ix: i32 = kx;
                                var k: i32 = kk;
                                while (k < kk + j) : (k += 1) {
                                    ops.add_( // x[ix] += temp * conj(ap[k])
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(ap[scast(u32, k)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix += incx;
                                }

                                if (nounit) {
                                    ops.mul_( // x[jx] *= conj(ap[kk + j])
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ops.conjugate(ap[scast(u32, kk + j)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        jx += incx;
                        kk += j + 1;
                    }
                }
            } else {
                var kk: i32 = int.div(n * (n + 1), 2) - 1;
                if (incx == 1) {
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable) {
                            const temp: X = x[scast(u32, j)];

                            if (noconj) {
                                var k: i32 = kk;
                                var i: i32 = n - 1;
                                while (i > j) : (i -= 1) {
                                    ops.add_( // x[i] += temp * ap[k]
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
                                if (nounit) {
                                    ops.mul_( // x[j] *= ap[kk - (n - 1) + j]
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ap[scast(u32, kk - (n - 1) + j)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var k: i32 = kk;
                                var i: i32 = n - 1;
                                while (i > j) : (i -= 1) {
                                    ops.add_( // x[i] += temp * conj(ap[k])
                                        &x[scast(u32, i)],
                                        x[scast(u32, i)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(ap[scast(u32, k)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    k -= 1;
                                }
                                if (nounit) {
                                    ops.mul_( // x[j] *= conj(ap[kk - (n - 1) + j])
                                        &x[scast(u32, j)],
                                        x[scast(u32, j)],
                                        ops.conjugate(ap[scast(u32, kk - (n - 1) + j)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        kk -= n - j;
                    }
                } else {
                    kx += (n - 1) * incx;
                    var jx: i32 = kx;
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable) {
                            const temp: X = x[scast(u32, jx)];

                            if (noconj) {
                                var ix: i32 = kx;
                                var k: i32 = kk;
                                while (k > kk - (n - (j + 1))) : (k -= 1) {
                                    ops.add_( // x[ix] += temp * ap[k]
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ap[scast(u32, k)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix -= incx;
                                }

                                if (nounit) {
                                    ops.mul_( // x[jx] *= ap[kk - (n - 1) + j]
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ap[scast(u32, kk - (n - 1) + j)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var ix: i32 = kx;
                                var k: i32 = kk;
                                while (k > kk - (n - (j + 1))) : (k -= 1) {
                                    ops.add_( // x[ix] += temp * conj(ap[k])
                                        &x[scast(u32, ix)],
                                        x[scast(u32, ix)],
                                        ops.mul(
                                            temp,
                                            ops.conjugate(ap[scast(u32, k)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;

                                    ix -= incx;
                                }

                                if (nounit) {
                                    ops.mul_( // x[jx] *= conj(ap[kk - (n - 1) + j])
                                        &x[scast(u32, jx)],
                                        x[scast(u32, jx)],
                                        ops.conjugate(ap[scast(u32, kk - (n - 1) + j)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        jx -= incx;
                        kk -= n - j;
                    }
                }
            }
        } else {
            if (uplo == .upper) {
                var kk: i32 = int.div(n * (n + 1), 2) - 1;
                if (incx == 1) {
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(u32, j)]);

                        if (noconj) {
                            if (nounit) {
                                ops.mul_( // temp *= ap[kk]
                                    &temp,
                                    temp,
                                    ap[scast(u32, kk)],
                                    ctx,
                                ) catch unreachable;
                            }

                            var k: i32 = kk - 1;
                            var i: i32 = j - 1;
                            while (i >= 0) : (i -= 1) {
                                ops.add_( // temp += ap[k] * x[i]
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
                        } else {
                            if (nounit) {
                                ops.mul_( // temp *= conj(ap[kk])
                                    &temp,
                                    temp,
                                    ops.conjugate(ap[scast(u32, kk)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var k: i32 = kk - 1;
                            var i: i32 = j - 1;
                            while (i >= 0) : (i -= 1) {
                                ops.add_( // temp += conj(ap[k]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(ap[scast(u32, k)], ctx) catch unreachable,
                                        x[scast(u32, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k -= 1;
                            }
                        }

                        x[scast(u32, j)] = scast(X, temp);

                        kk -= j + 1;
                    }
                } else {
                    var jx: i32 = kx + (n - 1) * incx;
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        var temp: C1 = scast(C1, x[scast(u32, jx)]);

                        var ix: i32 = jx;
                        if (noconj) {
                            if (nounit) {
                                ops.mul_( // temp *= ap[kk]
                                    &temp,
                                    temp,
                                    ap[scast(u32, kk)],
                                    ctx,
                                ) catch unreachable;
                            }

                            var k: i32 = kk - 1;
                            while (k >= kk - j) : (k -= 1) {
                                ix -= incx;

                                ops.add_( // temp += ap[k] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(u32, k)],
                                        x[scast(u32, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            if (nounit) {
                                ops.mul_( // temp *= conj(ap[kk])
                                    &temp,
                                    temp,
                                    ops.conjugate(ap[scast(u32, kk)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var k: i32 = kk - 1;
                            while (k >= kk - j) : (k -= 1) {
                                ix -= incx;

                                ops.add_( // temp += conj(ap[k]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(ap[scast(u32, k)], ctx) catch unreachable,
                                        x[scast(u32, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, jx)] = scast(X, temp);

                        jx -= incx;
                        kk -= j + 1;
                    }
                }
            } else {
                var kk: i32 = 0;
                if (incx == 1) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(u32, j)]);

                        var k: i32 = kk + 1;
                        if (noconj) {
                            if (nounit) {
                                ops.mul_( // temp *= ap[kk]
                                    &temp,
                                    temp,
                                    ap[scast(u32, kk)],
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // temp += ap[k] * x[i]
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
                        } else {
                            if (nounit) {
                                ops.mul_( // temp *= conj(ap[kk])
                                    &temp,
                                    temp,
                                    ops.conjugate(ap[scast(u32, kk)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // temp += conj(ap[k]) * x[i]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(ap[scast(u32, k)], ctx) catch unreachable,
                                        x[scast(u32, i)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                k += 1;
                            }
                        }

                        x[scast(u32, j)] = scast(X, temp);

                        kk += n - j;
                    }
                } else {
                    var jx: i32 = kx;
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var temp: C1 = scast(C1, x[scast(u32, jx)]);

                        if (noconj) {
                            if (nounit) {
                                ops.mul_( // temp *= ap[kk]
                                    &temp,
                                    temp,
                                    ap[scast(u32, kk)],
                                    ctx,
                                ) catch unreachable;
                            }

                            var ix: i32 = jx;
                            var k: i32 = kk + 1;
                            while (k < kk + n - j) : (k += 1) {
                                ix += incx;

                                ops.add_( // temp += ap[k] * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ap[scast(u32, k)],
                                        x[scast(u32, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            if (nounit) {
                                ops.mul_( // temp *= conj(ap[kk])
                                    &temp,
                                    temp,
                                    ops.conjugate(ap[scast(u32, kk)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var ix: i32 = jx;
                            var k: i32 = kk + 1;
                            while (k <= kk + n - j - 1) : (k += 1) {
                                ix += incx;

                                ops.add_( // temp += conj(ap[k]) * x[ix]
                                    &temp,
                                    temp,
                                    ops.mul(
                                        ops.conjugate(ap[scast(u32, k)], ctx) catch unreachable,
                                        x[scast(u32, ix)],
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        x[scast(u32, jx)] = scast(X, temp);
                        jx += incx;
                        kk += n - j;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.tpmv not implemented for arbitrary precision types yet");
    }

    return;
}
