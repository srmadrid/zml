const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Side = blas.Side;
const Uplo = blas.Uplo;
const Transpose = blas.Transpose;
const Diag = blas.Diag;

pub inline fn trsm(comptime T: type, order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]T, ldb: isize) void {
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
        .bool => @compileError("blas.trsm does not support bool."),
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    B[@intCast(ibkj)] /= A[@intCast(k + jak)];

                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)] -= B[@intCast(ibkj)] * A[@intCast(iaik)];

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)] -= B[@intCast(ibkj)] * A[@intCast(iaik)];

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
                                    var t0 = B[@intCast(ibij)] * alpha;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0 -= A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0 / A[@intCast(i + jai)];

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
                                    var t0 = B[@intCast(ibij)] * alpha;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0 -= A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    B[@intCast(ibkj)] /= A[@intCast(k + jak)];

                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= B[@intCast(ibkj)] * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= B[@intCast(ibkj)] * A[@intCast(iaik)];

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

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
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0 = B[@intCast(ibij)] * alpha;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0 -= A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0 / A[@intCast(i + jai)];

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
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0 = B[@intCast(ibij)] * alpha;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0 -= A[@intCast(iaki)] * B[@intCast(ibkj)];

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
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
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= A[@intCast(iakj)] * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] /= A[@intCast(j + jaj)];

                                    i += 1;
                                    ibij += 1;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= A[@intCast(iakj)] * B[@intCast(ibik)];

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
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] /= A[@intCast(k + jak)];

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= alpha;

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
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= alpha;

                                    i += 1;
                                    ibik += 1;
                                }

                                k -= 1;
                                jak -= lda;
                                jbk -= ldb;
                            }
                        }
                    }
                } else {
                    if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= A[@intCast(iakj)] * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] /= A[@intCast(j + jaj)];

                                    i += 1;
                                    ibij += 1;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    B[@intCast(ibij)] *= alpha;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= A[@intCast(iakj)] * B[@intCast(ibik)];

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
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] /= A[@intCast(k + jak)];

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= alpha;

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
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)] -= t0 * B[@intCast(ibik)];

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    B[@intCast(ibik)] *= alpha;

                                    i += 1;
                                    ibik += 1;
                                }

                                k += 1;
                                jak += lda;
                                jbk += ldb;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(i)].im) < @abs(A[@intCast(i)].re)) {
                                        const temp1 = A[@intCast(i)].im / A[@intCast(i)].re;
                                        const temp2 = A[@intCast(i)].re + temp1 * A[@intCast(i)].im;
                                        temp.re = (B[@intCast(ibkj)].re + temp1 * B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (B[@intCast(ibkj)].im - temp1 * B[@intCast(ibkj)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(i)].re / A[@intCast(i)].im;
                                        const temp2 = A[@intCast(i)].im + temp1 * A[@intCast(i)].re;
                                        temp.re = (temp1 * B[@intCast(ibkj)].re + B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibkj)].im - B[@intCast(ibkj)].re) / temp2;
                                    }
                                    B[@intCast(ibkj)] = temp;

                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re - B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].re * A[@intCast(iaik)].im + B[@intCast(ibkj)].im * A[@intCast(iaik)].re;

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re - B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].re * A[@intCast(iaik)].im + B[@intCast(ibkj)].im * A[@intCast(iaik)].re;

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(i)].im) < @abs(A[@intCast(i)].re)) {
                                        const temp1 = -A[@intCast(i)].im / A[@intCast(i)].re;
                                        const temp2 = A[@intCast(i)].re - temp1 * A[@intCast(i)].im;
                                        temp.re = (B[@intCast(ibkj)].re + temp1 * B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (B[@intCast(ibkj)].im - temp1 * B[@intCast(ibkj)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(i)].re / A[@intCast(i)].im;
                                        const temp2 = -A[@intCast(i)].im + temp1 * A[@intCast(i)].re;
                                        temp.re = (temp1 * B[@intCast(ibkj)].re + B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibkj)].im - B[@intCast(ibkj)].re) / temp2;
                                    }
                                    B[@intCast(ibkj)] = temp;

                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re + B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].im * A[@intCast(iaik)].re - B[@intCast(ibkj)].re * A[@intCast(iaik)].im;

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = M - 1;
                                var jak: isize = (M - 1) * lda;
                                var ibkj: isize = M - 1 + jbj;
                                while (k >= 0) {
                                    i = 0;
                                    var iaik: isize = jak;
                                    ibij = jbj;
                                    while (i < k) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re + B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].im * A[@intCast(iaik)].re - B[@intCast(ibkj)].re * A[@intCast(iaik)].im;

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
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    k = i + jai;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (t0.re + temp1 * t0.im) / temp2;
                                        temp.im = (t0.im - temp1 * t0.re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * t0.re + t0.im) / temp2;
                                        temp.im = (temp1 * t0.im - t0.re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

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
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

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
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    k = i + jai;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = -A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re - temp1 * A[@intCast(k)].im;
                                        temp.re = (t0.re + temp1 * t0.im) / temp2;
                                        temp.im = (t0.im - temp1 * t0.re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = -A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * t0.re + t0.im) / temp2;
                                        temp.im = (temp1 * t0.im - t0.re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

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
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = 0;
                                    var iaki: isize = jai;
                                    var ibkj: isize = jbj;
                                    while (k < i) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i += 1;
                                    jai += lda;
                                    ibij += 1;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(i)].im) < @abs(A[@intCast(i)].re)) {
                                        const temp1 = A[@intCast(i)].im / A[@intCast(i)].re;
                                        const temp2 = A[@intCast(i)].re + temp1 * A[@intCast(i)].im;
                                        temp.re = (B[@intCast(ibkj)].re + temp1 * B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (B[@intCast(ibkj)].im - temp1 * B[@intCast(ibkj)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(i)].re / A[@intCast(i)].im;
                                        const temp2 = A[@intCast(i)].im + temp1 * A[@intCast(i)].re;
                                        temp.re = (temp1 * B[@intCast(ibkj)].re + B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibkj)].im - B[@intCast(ibkj)].re) / temp2;
                                    }
                                    B[@intCast(ibkj)] = temp;

                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re - B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].re * A[@intCast(iaik)].im + B[@intCast(ibkj)].im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re - B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].re * A[@intCast(iaik)].im + B[@intCast(ibkj)].im * A[@intCast(iaik)].re;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(i)].im) < @abs(A[@intCast(i)].re)) {
                                        const temp1 = -A[@intCast(i)].im / A[@intCast(i)].re;
                                        const temp2 = A[@intCast(i)].re - temp1 * A[@intCast(i)].im;
                                        temp.re = (B[@intCast(ibkj)].re + temp1 * B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (B[@intCast(ibkj)].im - temp1 * B[@intCast(ibkj)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(i)].re / A[@intCast(i)].im;
                                        const temp2 = -A[@intCast(i)].im + temp1 * A[@intCast(i)].re;
                                        temp.re = (temp1 * B[@intCast(ibkj)].re + B[@intCast(ibkj)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibkj)].im - B[@intCast(ibkj)].re) / temp2;
                                    }
                                    B[@intCast(ibkj)] = temp;

                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re + B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].im * A[@intCast(iaik)].re - B[@intCast(ibkj)].re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var jak: isize = 0;
                                var ibkj: isize = jbj;
                                while (k < M) {
                                    i = k + 1;
                                    var iaik: isize = k + 1 + jak;
                                    ibij = k + 1 + jbj;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= B[@intCast(ibkj)].re * A[@intCast(iaik)].re + B[@intCast(ibkj)].im * A[@intCast(iaik)].im;
                                        B[@intCast(ibij)].im -= B[@intCast(ibkj)].im * A[@intCast(iaik)].re - B[@intCast(ibkj)].re * A[@intCast(iaik)].im;

                                        i += 1;
                                        iaik += 1;
                                        ibij += 1;
                                    }

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
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    k = i + jai;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (t0.re + temp1 * t0.im) / temp2;
                                        temp.im = (t0.im - temp1 * t0.re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * t0.re + t0.im) / temp2;
                                        temp.im = (temp1 * t0.im - t0.re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

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
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re - A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im + A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

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
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    k = i + jai;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = -A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re - temp1 * A[@intCast(k)].im;
                                        temp.re = (t0.re + temp1 * t0.im) / temp2;
                                        temp.im = (t0.im - temp1 * t0.re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = -A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * t0.re + t0.im) / temp2;
                                        temp.im = (temp1 * t0.im - t0.re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

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
                                var ibij: isize = M - 1 + jbj;
                                while (i >= 0) {
                                    var t0: T = undefined;
                                    t0.re = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    t0.im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;

                                    var k: isize = i + 1;
                                    var iaki: isize = i + 1 + jai;
                                    var ibkj: isize = i + 1 + jbj;
                                    while (k < M) {
                                        t0.re -= A[@intCast(iaki)].re * B[@intCast(ibkj)].re + A[@intCast(iaki)].im * B[@intCast(ibkj)].im;
                                        t0.im -= A[@intCast(iaki)].re * B[@intCast(ibkj)].im - A[@intCast(iaki)].im * B[@intCast(ibkj)].re;

                                        k += 1;
                                        iaki += 1;
                                        ibkj += 1;
                                    }

                                    B[@intCast(ibij)] = t0;

                                    i -= 1;
                                    jai -= lda;
                                    ibij -= 1;
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
                            var j: isize = 0;
                            var jaj: isize = 0;
                            var jbj: isize = 0;
                            while (j < N) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re - A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im + A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    k = j + jaj;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (B[@intCast(ibij)].re + temp1 * B[@intCast(ibij)].im) / temp2;
                                        temp.im = (B[@intCast(ibij)].im - temp1 * B[@intCast(ibij)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * B[@intCast(ibij)].re + B[@intCast(ibij)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibij)].im - B[@intCast(ibij)].re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    ibij += 1;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re - A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im + A[@intCast(iakj)].im * B[@intCast(ibik)].re;

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re + A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im - A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    k = j + jaj;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = -A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re - temp1 * A[@intCast(k)].im;
                                        temp.re = (B[@intCast(ibij)].re + temp1 * B[@intCast(ibij)].im) / temp2;
                                        temp.im = (B[@intCast(ibij)].im - temp1 * B[@intCast(ibij)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = -A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * B[@intCast(ibij)].re + B[@intCast(ibij)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibij)].im - B[@intCast(ibij)].re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    ibij += 1;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = 0;
                                var iakj: isize = jaj;
                                var jbk: isize = 0;
                                while (k < j) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re + A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im - A[@intCast(iakj)].im * B[@intCast(ibik)].re;

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
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const j: isize = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(j)].im) < @abs(A[@intCast(j)].re)) {
                                        const temp1 = A[@intCast(j)].im / A[@intCast(j)].re;
                                        const temp2 = A[@intCast(j)].re + temp1 * A[@intCast(j)].im;
                                        temp.re = (B[@intCast(ibik)].re + temp1 * B[@intCast(ibik)].im) / temp2;
                                        temp.im = (B[@intCast(ibik)].im - temp1 * B[@intCast(ibik)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(j)].re / A[@intCast(j)].im;
                                        const temp2 = A[@intCast(j)].im + temp1 * A[@intCast(j)].re;
                                        temp.re = (temp1 * B[@intCast(ibik)].re + B[@intCast(ibik)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibik)].im - B[@intCast(ibik)].re) / temp2;
                                    }
                                    B[@intCast(ibik)] = temp;

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
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
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
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
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const j: isize = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(j)].im) < @abs(A[@intCast(j)].re)) {
                                        const temp1 = -A[@intCast(j)].im / A[@intCast(j)].re;
                                        const temp2 = A[@intCast(j)].re - temp1 * A[@intCast(j)].im;
                                        temp.re = (B[@intCast(ibik)].re + temp1 * B[@intCast(ibik)].im) / temp2;
                                        temp.im = (B[@intCast(ibik)].im - temp1 * B[@intCast(ibik)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(j)].re / A[@intCast(j)].im;
                                        const temp2 = -A[@intCast(j)].im + temp1 * A[@intCast(j)].re;
                                        temp.re = (temp1 * B[@intCast(ibik)].re + B[@intCast(ibik)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibik)].im - B[@intCast(ibik)].re) / temp2;
                                    }
                                    B[@intCast(ibik)] = temp;

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re + t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im - t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
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
                                var j: isize = 0;
                                var iajk: isize = jak;
                                var jbj: isize = 0;
                                while (j < k) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re + t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im - t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
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
                } else {
                    if (TRANSA == .NoTrans) {
                        if (diag == .NonUnit) {
                            var j: isize = N - 1;
                            var jaj: isize = (N - 1) * lda;
                            var jbj: isize = (N - 1) * ldb;
                            while (j >= 0) {
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re - A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im + A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    k = j + jaj;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (B[@intCast(ibij)].re + temp1 * B[@intCast(ibij)].im) / temp2;
                                        temp.im = (B[@intCast(ibij)].im - temp1 * B[@intCast(ibij)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * B[@intCast(ibij)].re + B[@intCast(ibij)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibij)].im - B[@intCast(ibij)].re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    ibij += 1;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re - A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im + A[@intCast(iakj)].im * B[@intCast(ibik)].re;

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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re + A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im - A[@intCast(iakj)].im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    k += 1;
                                    iakj += 1;
                                    jbk += ldb;
                                }

                                i = 0;
                                ibij = jbj;
                                while (i < M) {
                                    k = j + jaj;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(k)].im) < @abs(A[@intCast(k)].re)) {
                                        const temp1 = -A[@intCast(k)].im / A[@intCast(k)].re;
                                        const temp2 = A[@intCast(k)].re + temp1 * A[@intCast(k)].im;
                                        temp.re = (B[@intCast(ibij)].re + temp1 * B[@intCast(ibij)].im) / temp2;
                                        temp.im = (B[@intCast(ibij)].im - temp1 * B[@intCast(ibij)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(k)].re / A[@intCast(k)].im;
                                        const temp2 = -A[@intCast(k)].im + temp1 * A[@intCast(k)].re;
                                        temp.re = (temp1 * B[@intCast(ibij)].re + B[@intCast(ibij)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibij)].im - B[@intCast(ibij)].re) / temp2;
                                    }
                                    B[@intCast(ibij)] = temp;

                                    i += 1;
                                    ibij += 1;
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
                                var i: isize = 0;
                                var ibij: isize = jbj;
                                while (i < M) {
                                    const tmp = B[@intCast(ibij)].re * alpha.re - B[@intCast(ibij)].im * alpha.im;
                                    B[@intCast(ibij)].im = B[@intCast(ibij)].re * alpha.im + B[@intCast(ibij)].im * alpha.re;
                                    B[@intCast(ibij)].re = tmp;

                                    i += 1;
                                    ibij += 1;
                                }

                                var k: isize = j + 1;
                                var iakj: isize = j + 1 + jaj;
                                var jbk: isize = (j + 1) * ldb;
                                while (k < N) {
                                    i = 0;
                                    ibij = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= A[@intCast(iakj)].re * B[@intCast(ibik)].re + A[@intCast(iakj)].im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= A[@intCast(iakj)].re * B[@intCast(ibik)].im - A[@intCast(iakj)].im * B[@intCast(ibik)].re;

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
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const j: isize = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(j)].im) < @abs(A[@intCast(j)].re)) {
                                        const temp1 = A[@intCast(j)].im / A[@intCast(j)].re;
                                        const temp2 = A[@intCast(j)].re + temp1 * A[@intCast(j)].im;
                                        temp.re = (B[@intCast(ibik)].re + temp1 * B[@intCast(ibik)].im) / temp2;
                                        temp.im = (B[@intCast(ibik)].im - temp1 * B[@intCast(ibik)].re) / temp2;
                                    } else {
                                        const temp1 = A[@intCast(j)].re / A[@intCast(j)].im;
                                        const temp2 = A[@intCast(j)].im + temp1 * A[@intCast(j)].re;
                                        temp.re = (temp1 * B[@intCast(ibik)].re + B[@intCast(ibik)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibik)].im - B[@intCast(ibik)].re) / temp2;
                                    }
                                    B[@intCast(ibik)] = temp;

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
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
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re - t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im + t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
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
                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const j: isize = k + jak;
                                    var temp: T = undefined;
                                    if (@abs(A[@intCast(j)].im) < @abs(A[@intCast(j)].re)) {
                                        const temp1 = -A[@intCast(j)].im / A[@intCast(j)].re;
                                        const temp2 = A[@intCast(j)].re - temp1 * A[@intCast(j)].im;
                                        temp.re = (B[@intCast(ibik)].re + temp1 * B[@intCast(ibik)].im) / temp2;
                                        temp.im = (B[@intCast(ibik)].im - temp1 * B[@intCast(ibik)].re) / temp2;
                                    } else {
                                        const temp1 = -A[@intCast(j)].re / A[@intCast(j)].im;
                                        const temp2 = -A[@intCast(j)].im + temp1 * A[@intCast(j)].re;
                                        temp.re = (temp1 * B[@intCast(ibik)].re + B[@intCast(ibik)].im) / temp2;
                                        temp.im = (temp1 * B[@intCast(ibik)].im - B[@intCast(ibik)].re) / temp2;
                                    }
                                    B[@intCast(ibik)] = temp;

                                    i += 1;
                                    ibik += 1;
                                }

                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    i = 0;
                                    var ibij: isize = jbj;
                                    ibik = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re + t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im - t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                i = 0;
                                ibik = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
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
                                var j: isize = k + 1;
                                var iajk: isize = k + 1 + jak;
                                var jbj: isize = (k + 1) * ldb;
                                while (j < N) {
                                    const t0 = A[@intCast(iajk)];

                                    var i: isize = 0;
                                    var ibij: isize = jbj;
                                    var ibik: isize = jbk;
                                    while (i < M) {
                                        B[@intCast(ibij)].re -= t0.re * B[@intCast(ibik)].re + t0.im * B[@intCast(ibik)].im;
                                        B[@intCast(ibij)].im -= t0.re * B[@intCast(ibik)].im - t0.im * B[@intCast(ibik)].re;

                                        i += 1;
                                        ibij += 1;
                                        ibik += 1;
                                    }

                                    j += 1;
                                    iajk += 1;
                                    jbj += ldb;
                                }

                                var i: isize = 0;
                                var ibik: isize = jbk;
                                while (i < M) {
                                    const tmp = B[@intCast(ibik)].re * alpha.re - B[@intCast(ibik)].im * alpha.im;
                                    B[@intCast(ibik)].im = B[@intCast(ibik)].re * alpha.im + B[@intCast(ibik)].im * alpha.re;
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
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.trsm only supports simple types."),
    }
}
