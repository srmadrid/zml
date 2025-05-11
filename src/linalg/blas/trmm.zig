const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Side = blas.Side;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;
const Diag = blas.Diag;

pub inline fn trmm(comptime T: type, order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]T, ldb: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0 or m <= 0) return;

    var M = m;
    var N = n;
    var SIDE = side;
    var UPLO = uplo;
    const TRANSA = transA;
    if (order == .RowMajor) {
        M = n;
        N = m;
        SIDE = if (side == .Left) .Right else .Left;
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    if (SIDE == .Left and lda < @max(1, M)) return;
    if (SIDE == .Right and lda < @max(1, N)) return;
    if (ldb < @max(1, M)) return;

    switch (numericType) {
        .bool => @compileError("blas.trmm does not support bool."),
        .int, .float => {
            if (alpha == 0) {
                var j: isize = 0;
                var jbj: isize = 0;
                while (j < N) {
                    var i: isize = 0;
                    var ijbj: isize = jbj;
                    while (i < M) {
                        B[@intCast(ijbj)] = 0;

                        i += 1;
                        ijbj += 1;
                    }

                    j += 1;
                    jbj += ldb;
                }

                return;
            }

            if (SIDE == .Left) {
                if (UPLO == .Upper) {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    const t0 = alpha * B[@intCast(ibkj)];

                                    var i: isize = 0;
                                    var iaik: isize = jak;
                                    var ibij: isize = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)] += t0 * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    B[@intCast(ibkj)] = t0 * A[@intCast(iaik)];

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    const t0 = alpha * B[@intCast(ibkj)];

                                    var i: isize = 0;
                                    var iaik: isize = jak;
                                    var ibij: isize = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)] += t0 * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    B[@intCast(ibkj)] = t0;

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = jbj + M - 1;
                                while (i >= 0) {
                                    var t0 = B[@intCast(ibij)] * A[@intCast(i + jai)];

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0 += A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = alpha * t0;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = jbj + M - 1;
                                while (i >= 0) {
                                    var t0 = B[@intCast(ibij)];

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0 += A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = alpha * t0;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = jbj + M - 1;
                                while (k >= 0) {
                                    const t0 = alpha * B[@intCast(ibkj)];
                                    B[@intCast(ibkj)] = t0 * A[@intCast(k + jak)];

                                    var i: isize = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    var ibij: isize = jbj + k + 1;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = jbj + M - 1;
                                while (k >= 0) {
                                    const t0 = alpha * B[@intCast(ibkj)];
                                    B[@intCast(ibkj)] = t0;

                                    var i: isize = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    var ibij: isize = jbj + k + 1;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0 = B[@intCast(ibij)] * A[@intCast(i + jai)];

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = jbj + i + 1;
                                    while (k < M) {
                                        t0 += A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = alpha * t0;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0 = B[@intCast(ibij)];

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = jbj + i + 1;
                                    while (k < M) {
                                        t0 += A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = alpha * t0;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    }
                }
            } else {
                if (UPLO == .Upper) {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var t0 = alpha * A[@intCast(j + jaj)];

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= t0;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    t0 = alpha * A[@intCast(iakj)];

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        } else {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var t0 = alpha;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= t0;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    t0 = alpha * A[@intCast(iakj)];

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = alpha * A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                const t0 = alpha * A[@intCast(iajk)];

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= t0;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        } else {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = alpha * A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                const t0 = alpha;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= t0;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var t0 = alpha * A[@intCast(j + jaj)];

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= t0;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    t0 = alpha * A[@intCast(iakj)];

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var t0 = alpha;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= t0;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    t0 = alpha * A[@intCast(iakj)];

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = alpha * A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                const t0 = alpha * A[@intCast(k + jak)];

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= t0;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        } else {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = alpha * A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] += t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                const t0 = alpha;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= t0;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        }
                    }
                }
            }
        },
        .cfloat => {
            if (alpha.re == 0 and alpha.im == 0) {
                var j: isize = 0;
                var jbj: isize = 0;
                while (j < N) {
                    var i: isize = 0;
                    var ijbj: isize = jbj;
                    while (i < M) {
                        B[@intCast(ijbj)].re = 0;
                        B[@intCast(ijbj)].im = 0;

                        i += 1;
                        ijbj += 1;
                    }

                    j += 1;
                    jbj += ldb;
                }

                return;
            }

            if (SIDE == .Left) {
                if (UPLO == .Upper) {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * B[@intCast(ibkj)].re - alpha.im * B[@intCast(ibkj)].im;
                                    t0.im = alpha.re * B[@intCast(ibkj)].im + alpha.im * B[@intCast(ibkj)].re;

                                    var i: isize = 0;
                                    var iaik: isize = jak;
                                    var ibij: isize = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re += t0.re * A[@intCast(iaik)].re - t0.im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im += t0.re * A[@intCast(iaik)].im + t0.im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    B[@intCast(ibkj)].re = t0.re * A[@intCast(iaik)].re - t0.im * A[@intCast(iaik)].im;
                                    B[@intCast(ibkj)].im = t0.re * A[@intCast(iaik)].im + t0.im * A[@intCast(iaik)].re;

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * B[@intCast(ibkj)].re - alpha.im * B[@intCast(ibkj)].im;
                                    t0.im = alpha.re * B[@intCast(ibkj)].im + alpha.im * B[@intCast(ibkj)].re;

                                    var i: isize = 0;
                                    var iaik: isize = jak;
                                    var ibij: isize = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re += t0.re * A[@intCast(iaik)].re - t0.im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im += t0.re * A[@intCast(iaik)].im + t0.im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    B[@intCast(ibkj)] = t0;

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * B[@intCast(ibkj)].re - alpha.im * B[@intCast(ibkj)].im;
                                    t0.im = alpha.re * B[@intCast(ibkj)].im + alpha.im * B[@intCast(ibkj)].re;

                                    var i: isize = 0;
                                    var iaik: isize = jak;
                                    var ibij: isize = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re += t0.re * A[@intCast(iaik)].re + t0.im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im += t0.im * A[@intCast(iaik)].re - t0.re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    B[@intCast(ibkj)].re = t0.re * A[@intCast(iaik)].re + t0.im * A[@intCast(iaik)].im;
                                    B[@intCast(ibkj)].im = t0.im * A[@intCast(iaik)].re - t0.re * A[@intCast(iaik)].im;

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * B[@intCast(ibkj)].re - alpha.im * B[@intCast(ibkj)].im;
                                    t0.im = alpha.re * B[@intCast(ibkj)].im + alpha.im * B[@intCast(ibkj)].re;

                                    var i: isize = 0;
                                    var iaik: isize = jak;
                                    var ibij: isize = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re += t0.re * A[@intCast(iaik)].re + t0.im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im += t0.im * A[@intCast(iaik)].re - t0.re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    B[@intCast(ibkj)] = t0;

                                    k += 1;
                                    jak += lda;
                                    ibkj += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .Trans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = jbj + M - 1;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = A[@intCast(i + jai)].re * B[@intCast(ibij)].re - A[@intCast(i + jai)].im * B[@intCast(ibij)].im;
                                    t0.im = A[@intCast(i + jai)].re * B[@intCast(ibij)].im + A[@intCast(i + jai)].im * B[@intCast(ibij)].re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)].re = alpha.re * t0.re - alpha.im * t0.im;
                                    B[@intCast(ibij)].im = alpha.re * t0.im + alpha.im * t0.re;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = jbj + M - 1;
                                while (i >= 0) {
                                    var t0 = B[@intCast(ibij)];

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)].re = alpha.re * t0.re - alpha.im * t0.im;
                                    B[@intCast(ibij)].im = alpha.re * t0.im + alpha.im * t0.re;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = jbj + M - 1;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = A[@intCast(i + jai)].re * B[@intCast(ibij)].re + A[@intCast(i + jai)].im * B[@intCast(ibij)].im;
                                    t0.im = A[@intCast(i + jai)].re * B[@intCast(ibij)].im - A[@intCast(i + jai)].im * B[@intCast(ibij)].re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)].re = alpha.re * t0.re - alpha.im * t0.im;
                                    B[@intCast(ibij)].im = alpha.re * t0.im + alpha.im * t0.re;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = M - 1;
                                var jai: isize = (M - 1) * lda;
                                var ibij: isize = jbj + M - 1;
                                while (i >= 0) {
                                    var t0 = B[@intCast(ibij)];

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)].re = alpha.re * t0.re - alpha.im * t0.im;
                                    B[@intCast(ibij)].im = alpha.re * t0.im + alpha.im * t0.re;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = jbj + M - 1;
                                while (k >= 0) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * B[@intCast(ibkj)].re - alpha.im * B[@intCast(ibkj)].im;
                                    t0.im = alpha.re * B[@intCast(ibkj)].im + alpha.im * B[@intCast(ibkj)].re;
                                    B[@intCast(ibkj)].re = t0.re * A[@intCast(k + jak)].re - t0.im * A[@intCast(k + jak)].im;
                                    B[@intCast(ibkj)].im = t0.re * A[@intCast(k + jak)].im + t0.im * A[@intCast(k + jak)].re;

                                    var i: isize = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    var ibij: isize = jbj + k + 1;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * A[@intCast(iaik)].re - t0.im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im += t0.re * A[@intCast(iaik)].im + t0.im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = jbj + M - 1;
                                while (k >= 0) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * B[@intCast(ibkj)].re - alpha.im * B[@intCast(ibkj)].im;
                                    t0.im = alpha.re * B[@intCast(ibkj)].im + alpha.im * B[@intCast(ibkj)].re;
                                    B[@intCast(ibkj)] = t0;

                                    var i: isize = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    var ibij: isize = jbj + k + 1;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * A[@intCast(iaik)].re - t0.im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im += t0.re * A[@intCast(iaik)].im + t0.im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = jbj + M - 1;
                                while (k >= 0) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * B[@intCast(ibkj)].re - alpha.im * B[@intCast(ibkj)].im;
                                    t0.im = alpha.re * B[@intCast(ibkj)].im + alpha.im * B[@intCast(ibkj)].re;
                                    B[@intCast(ibkj)].re = t0.re * A[@intCast(k + jak)].re + t0.im * A[@intCast(k + jak)].im;
                                    B[@intCast(ibkj)].im = t0.im * A[@intCast(k + jak)].re - t0.re * A[@intCast(k + jak)].im;

                                    var i: isize = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    var ibij: isize = jbj + k + 1;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * A[@intCast(iaik)].re + t0.im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im += t0.im * A[@intCast(iaik)].re - t0.re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = jbj + M - 1;
                                while (k >= 0) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * B[@intCast(ibkj)].re - alpha.im * B[@intCast(ibkj)].im;
                                    t0.im = alpha.re * B[@intCast(ibkj)].im + alpha.im * B[@intCast(ibkj)].re;
                                    B[@intCast(ibkj)] = t0;

                                    var i: isize = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    var ibij: isize = jbj + k + 1;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * A[@intCast(iaik)].re + t0.im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im += t0.im * A[@intCast(iaik)].re - t0.re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

                                    k -= 1;
                                    jak -= lda;
                                    ibkj -= 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .Trans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * A[@intCast(i + jai)].re - B[@intCast(ibij)].im * A[@intCast(i + jai)].im;
                                    t0.im = B[@intCast(ibij)].re * A[@intCast(i + jai)].im + B[@intCast(ibij)].im * A[@intCast(i + jai)].re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = jbj + i + 1;
                                    while (k < M) {
                                        t0.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)].re = alpha.re * t0.re - alpha.im * t0.im;
                                    B[@intCast(ibij)].im = alpha.re * t0.im + alpha.im * t0.re;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0 = B[@intCast(ibij)];

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = jbj + i + 1;
                                    while (k < M) {
                                        t0.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)].re = alpha.re * t0.re - alpha.im * t0.im;
                                    B[@intCast(ibij)].im = alpha.re * t0.im + alpha.im * t0.re;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * A[@intCast(i + jai)].re + B[@intCast(ibij)].im * A[@intCast(i + jai)].im;
                                    t0.im = B[@intCast(ibij)].im * A[@intCast(i + jai)].re - B[@intCast(ibij)].re * A[@intCast(i + jai)].im;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = jbj + i + 1;
                                    while (k < M) {
                                        t0.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)].re = alpha.re * t0.re - alpha.im * t0.im;
                                    B[@intCast(ibij)].im = alpha.re * t0.im + alpha.im * t0.re;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var jai: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    var t0 = B[@intCast(ibij)];

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = jbj + i + 1;
                                    while (k < M) {
                                        t0.re += A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im += A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)].re = alpha.re * t0.re - alpha.im * t0.im;
                                    B[@intCast(ibij)].im = alpha.re * t0.im + alpha.im * t0.re;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
                                }

                                j += 1;
                                jbj += ldb;
                            }
                        }
                    }
                }
            } else {
                if (UPLO == .Upper) {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var t0: T = undefined;
                                t0.re = alpha.re * A[@intCast(j + jaj)].re - alpha.im * A[@intCast(j + jaj)].im;
                                t0.im = alpha.re * A[@intCast(j + jaj)].im + alpha.im * A[@intCast(j + jaj)].re;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * t0.re - B[@intCast(ibij)].im * t0.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * t0.im + B[@intCast(ibij)].im * t0.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    t0.re = alpha.re * A[@intCast(iakj)].re - alpha.im * A[@intCast(iakj)].im;
                                    t0.im = alpha.re * A[@intCast(iakj)].im + alpha.im * A[@intCast(iakj)].re;

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        } else {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var t0 = alpha;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * t0.re - B[@intCast(ibij)].im * t0.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * t0.im + B[@intCast(ibij)].im * t0.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    t0.re = alpha.re * A[@intCast(iakj)].re - alpha.im * A[@intCast(iakj)].im;
                                    t0.im = alpha.re * A[@intCast(iakj)].im + alpha.im * A[@intCast(iakj)].re;

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        }
                    } else if (TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var t0: T = undefined;
                                t0.re = alpha.re * A[@intCast(j + jaj)].re + alpha.im * A[@intCast(j + jaj)].im;
                                t0.im = alpha.im * A[@intCast(j + jaj)].re - alpha.re * A[@intCast(j + jaj)].im;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * t0.re - B[@intCast(ibij)].im * t0.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * t0.im + B[@intCast(ibij)].im * t0.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    t0.re = alpha.re * A[@intCast(iakj)].re + alpha.im * A[@intCast(iakj)].im;
                                    t0.im = alpha.im * A[@intCast(iakj)].re - alpha.re * A[@intCast(iakj)].im;

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        } else {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var t0 = alpha;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * t0.re - B[@intCast(ibij)].im * t0.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * t0.im + B[@intCast(ibij)].im * t0.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    t0.re = alpha.re * A[@intCast(iakj)].re + alpha.im * A[@intCast(iakj)].im;
                                    t0.im = alpha.im * A[@intCast(iakj)].re - alpha.re * A[@intCast(iakj)].im;

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j -= 1;
                                jaj -= lda;
                                jbj -= ldb;
                            }
                        }
                    } else if (TRANSA == .Trans) {
                        if (diag == .NonUnit) {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * A[@intCast(iajk)].re - alpha.im * A[@intCast(iajk)].im;
                                    t0.im = alpha.re * A[@intCast(iajk)].im + alpha.im * A[@intCast(iajk)].re;

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var t0: T = undefined;
                                t0.re = alpha.re * A[@intCast(iajk)].re - alpha.im * A[@intCast(iajk)].im;
                                t0.im = alpha.re * A[@intCast(iajk)].im + alpha.im * A[@intCast(iajk)].re;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                    B[@intCast(ibik)].im = t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        } else {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * A[@intCast(iajk)].re - alpha.im * A[@intCast(iajk)].im;
                                    t0.im = alpha.re * A[@intCast(iajk)].im + alpha.im * A[@intCast(iajk)].re;

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                const t0 = alpha;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                    B[@intCast(ibik)].im = t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * A[@intCast(iajk)].re + alpha.im * A[@intCast(iajk)].im;
                                    t0.im = alpha.im * A[@intCast(iajk)].re - alpha.re * A[@intCast(iajk)].im;

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var t0: T = undefined;
                                t0.re = alpha.re * A[@intCast(iajk)].re + alpha.im * A[@intCast(iajk)].im;
                                t0.im = alpha.im * A[@intCast(iajk)].re - alpha.re * A[@intCast(iajk)].im;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                    B[@intCast(ibik)].im = t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        } else {
                            var k: isize = 0;
                            var jak: isize = 0;
                            var jbk: isize = 0;
                            while (k < N) {
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * A[@intCast(iajk)].re + alpha.im * A[@intCast(iajk)].im;
                                    t0.im = alpha.im * A[@intCast(iajk)].re - alpha.re * A[@intCast(iajk)].im;

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                const t0 = alpha;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                    B[@intCast(ibik)].im = t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var t0: T = undefined;
                                t0.re = alpha.re * A[@intCast(j + jaj)].re - alpha.im * A[@intCast(j + jaj)].im;
                                t0.im = alpha.re * A[@intCast(j + jaj)].im + alpha.im * A[@intCast(j + jaj)].re;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * t0.re - B[@intCast(ibij)].im * t0.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * t0.im + B[@intCast(ibij)].im * t0.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    t0.re = alpha.re * A[@intCast(iakj)].re - alpha.im * A[@intCast(iakj)].im;
                                    t0.im = alpha.re * A[@intCast(iakj)].im + alpha.im * A[@intCast(iakj)].re;

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var t0 = alpha;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * t0.re - B[@intCast(ibij)].im * t0.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * t0.im + B[@intCast(ibij)].im * t0.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    t0.re = alpha.re * A[@intCast(iakj)].re - alpha.im * A[@intCast(iakj)].im;
                                    t0.im = alpha.re * A[@intCast(iakj)].im + alpha.im * A[@intCast(iakj)].re;

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var t0: T = undefined;
                                t0.re = alpha.re * A[@intCast(j + jaj)].re + alpha.im * A[@intCast(j + jaj)].im;
                                t0.im = alpha.im * A[@intCast(j + jaj)].re - alpha.re * A[@intCast(j + jaj)].im;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * t0.re - B[@intCast(ibij)].im * t0.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * t0.im + B[@intCast(ibij)].im * t0.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    t0.re = alpha.re * A[@intCast(iakj)].re + alpha.im * A[@intCast(iakj)].im;
                                    t0.im = alpha.im * A[@intCast(iakj)].re - alpha.re * A[@intCast(iakj)].im;

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        } else {
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var t0 = alpha;

                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * t0.re - B[@intCast(ibij)].im * t0.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * t0.im + B[@intCast(ibij)].im * t0.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    t0.re = alpha.re * A[@intCast(iakj)].re + alpha.im * A[@intCast(iakj)].im;
                                    t0.im = alpha.im * A[@intCast(iakj)].re - alpha.re * A[@intCast(iakj)].im;

                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                j += 1;
                                jaj += lda;
                                jbj += ldb;
                            }
                        }
                    } else if (TRANSA == .Trans) {
                        if (diag == .NonUnit) {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * A[@intCast(iajk)].re - alpha.im * A[@intCast(iajk)].im;
                                    t0.im = alpha.re * A[@intCast(iajk)].im + alpha.im * A[@intCast(iajk)].re;

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var t0: T = undefined;
                                t0.re = alpha.re * A[@intCast(k + jak)].re - alpha.im * A[@intCast(k + jak)].im;
                                t0.im = alpha.re * A[@intCast(k + jak)].im + alpha.im * A[@intCast(k + jak)].re;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                    B[@intCast(ibik)].im = t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        } else {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * A[@intCast(iajk)].re - alpha.im * A[@intCast(iajk)].im;
                                    t0.im = alpha.re * A[@intCast(iajk)].im + alpha.im * A[@intCast(iajk)].re;

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                const t0 = alpha;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                    B[@intCast(ibik)].im = t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        }
                    } else {
                        if (diag == .NonUnit) {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * A[@intCast(iajk)].re + alpha.im * A[@intCast(iajk)].im;
                                    t0.im = alpha.im * A[@intCast(iajk)].re - alpha.re * A[@intCast(iajk)].im;

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var t0: T = undefined;
                                t0.re = alpha.re * A[@intCast(k + jak)].re + alpha.im * A[@intCast(k + jak)].im;
                                t0.im = alpha.im * A[@intCast(k + jak)].re - alpha.re * A[@intCast(k + jak)].im;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                    B[@intCast(ibik)].im = t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        } else {
                            var k: isize = N - 1;
                            var jak: isize = (N - 1) * lda;
                            var jbk: isize = (N - 1) * ldb;
                            while (k >= 0) {
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    var t0: T = undefined;
                                    t0.re = alpha.re * A[@intCast(iajk)].re + alpha.im * A[@intCast(iajk)].im;
                                    t0.im = alpha.im * A[@intCast(iajk)].re - alpha.re * A[@intCast(iajk)].im;

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re += t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im += t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                const t0 = alpha;

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                    B[@intCast(ibik)].im = t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;
                                    B[@intCast(ibik)].re = tmp;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        }
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.trmm only supports simple types."),
        .unsupported => unreachable,
    }
}
