const std = @import("std");
const core = @import("../../core/core.zig");
const Transpose = @import("../ndarray.zig").Transpose;
const Order = @import("../ndarray.zig").Order;
const Uplo = @import("../ndarray.zig").Uplo;
const Diag = @import("../ndarray.zig").Diag;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn tbmv(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, diag: Diag, n: isize, k: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0 or k < 0) return;

    const N = n;
    var UPLO = uplo;
    var TRANS = trans;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        TRANS = if (trans == .NoTrans) .Trans else if (trans == .ConjNoTrans) .ConjTrans else if (trans == .Trans) .NoTrans else .ConjNoTrans;
    }

    if (lda < k + 1) return;

    const LENX = N;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.tbmv does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (UPLO == .Upper) {
                if (TRANS == .NoTrans or TRANS == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)] += t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] *= A[@intCast(iaij)];

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)] += t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0: T = 0;
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0 + A[@intCast(iaij)] * x[@intCast(jx)];

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0: T = 0;
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] += t0;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                }
            } else {
                if (TRANS == .NoTrans or TRANS == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            x[@intCast(jx)] *= A[@intCast(jaj)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)] += t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)] += t0 * A[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            t0 *= A[@intCast(jaj)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                }
            }
        },
        .Complex => {
            if (UPLO == .Upper) {
                if (TRANS == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)].re = t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                            x[@intCast(jx)].im = t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                } else if (TRANS == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)].re = t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                            x[@intCast(jx)].im = t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var kx: isize = 0;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx + kx else kx;
                            while (i < j) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            if (j >= k) kx += incx;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                } else if (TRANS == .Trans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = T.init(0, 0);
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            t0.re += A[@intCast(iaij)].re * x[@intCast(jx)].re - A[@intCast(iaij)].im * x[@intCast(jx)].im;
                            t0.im += A[@intCast(iaij)].re * x[@intCast(jx)].im + A[@intCast(iaij)].im * x[@intCast(jx)].re;
                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = T.init(0, 0);
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)].re += t0.re;
                            x[@intCast(jx)].im += t0.im;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = T.init(0, 0);
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            t0.re += A[@intCast(iaij)].re * x[@intCast(jx)].re + A[@intCast(iaij)].im * x[@intCast(jx)].im;
                            t0.im += A[@intCast(iaij)].re * x[@intCast(jx)].im - A[@intCast(iaij)].im * x[@intCast(jx)].re;
                            x[@intCast(jx)] = t0;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            var t0 = T.init(0, 0);
                            const l: isize = k - j;
                            const I0: isize = if (j - k > 0) j - k else 0;

                            var i: isize = I0;
                            var iaij: isize = jaj + I0 + l;
                            var ix: isize = if (incx < 0) (-LENX + 1 + I0) * incx else I0 * incx;
                            while (i < j) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)].re += t0.re;
                            x[@intCast(jx)].im += t0.im;

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                }
            } else {
                if (TRANS == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            x[@intCast(jx)].re = t0.re * A[@intCast(jaj)].re - t0.im * A[@intCast(jaj)].im;
                            x[@intCast(jx)].im = t0.re * A[@intCast(jaj)].im + t0.im * A[@intCast(jaj)].re;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re - t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.re * A[@intCast(iaij)].im + t0.im * A[@intCast(iaij)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else if (TRANS == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            x[@intCast(jx)].re = t0.re * A[@intCast(jaj)].re + t0.im * A[@intCast(jaj)].im;
                            x[@intCast(jx)].im = t0.im * A[@intCast(jaj)].re - t0.re * A[@intCast(jaj)].im;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = lda * (N - 1);
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        while (j >= 0) {
                            const t0: T = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                x[@intCast(ix)].re += t0.re * A[@intCast(iaij)].re + t0.im * A[@intCast(iaij)].im;
                                x[@intCast(ix)].im += t0.im * A[@intCast(iaij)].re - t0.re * A[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else if (TRANS == .Trans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const temp = t0;
                            t0.re = temp.re * A[@intCast(jaj)].re - temp.im * A[@intCast(jaj)].im;
                            t0.im = temp.re * A[@intCast(jaj)].im + temp.im * A[@intCast(jaj)].re;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re - A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im + A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const temp = t0;
                            t0.re = temp.re * A[@intCast(jaj)].re + temp.im * A[@intCast(jaj)].im;
                            t0.im = temp.im * A[@intCast(jaj)].re - temp.re * A[@intCast(jaj)].im;
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];
                            const I1: isize = if (N - 1 > j + k) j + k else N - 1;

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i <= I1) {
                                t0.re += A[@intCast(iaij)].re * x[@intCast(ix)].re + A[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im += A[@intCast(iaij)].re * x[@intCast(ix)].im - A[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                        }
                    }
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.tbmv only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "tbmv" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const k = 1;

    const A = try a.alloc(f64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);

    @memcpy(A.ptr, &[_]f64{
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    });
    @memcpy(x1.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbmv(f64, .RowMajor, .Upper, .NoTrans, .NonUnit, n, k, A.ptr, k + 1, x1.ptr, 2);

    try std.testing.expectEqual(5, x1[0]);
    try std.testing.expectEqual(18, x1[2]);
    try std.testing.expectEqual(39, x1[4]);
    try std.testing.expectEqual(68, x1[6]);
    try std.testing.expectEqual(45, x1[8]);

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);

    @memcpy(x2.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbmv(f64, .RowMajor, .Upper, .NoTrans, .Unit, n, k, A.ptr, k + 1, x2.ptr, -2);

    try std.testing.expectEqual(5, x2[8]);
    try std.testing.expectEqual(14, x2[6]);
    try std.testing.expectEqual(27, x2[4]);
    try std.testing.expectEqual(44, x2[2]);
    try std.testing.expectEqual(5, x2[0]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);

    @memcpy(x3.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbmv(f64, .ColumnMajor, .Upper, .NoTrans, .NonUnit, n, k, A.ptr, k + 1, x3.ptr, 2);

    try std.testing.expectEqual(8, x3[0]);
    try std.testing.expectEqual(23, x3[2]);
    try std.testing.expectEqual(46, x3[4]);
    try std.testing.expectEqual(77, x3[6]);
    try std.testing.expectEqual(50, x3[8]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);

    @memcpy(x4.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbmv(f64, .ColumnMajor, .Upper, .NoTrans, .Unit, n, k, A.ptr, k + 1, x4.ptr, -2);

    try std.testing.expectEqual(7, x4[8]);
    try std.testing.expectEqual(17, x4[6]);
    try std.testing.expectEqual(31, x4[4]);
    try std.testing.expectEqual(49, x4[2]);
    try std.testing.expectEqual(5, x4[0]);

    const x5 = try a.alloc(f64, 2 * n);
    defer a.free(x5);

    @memcpy(x5.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbmv(f64, .RowMajor, .Upper, .Trans, .NonUnit, n, k, A.ptr, k + 1, x5.ptr, 2);

    try std.testing.expectEqual(1, x5[0]);
    try std.testing.expectEqual(8, x5[2]);
    try std.testing.expectEqual(23, x5[4]);
    try std.testing.expectEqual(46, x5[6]);
    try std.testing.expectEqual(77, x5[8]);

    const x6 = try a.alloc(f64, 2 * n);
    defer a.free(x6);

    @memcpy(x6.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbmv(f64, .RowMajor, .Upper, .Trans, .Unit, n, k, A.ptr, k + 1, x6.ptr, -2);

    try std.testing.expectEqual(1, x6[8]);
    try std.testing.expectEqual(4, x6[6]);
    try std.testing.expectEqual(11, x6[4]);
    try std.testing.expectEqual(22, x6[2]);
    try std.testing.expectEqual(37, x6[0]);

    const x7 = try a.alloc(f64, 2 * n);
    defer a.free(x7);

    @memcpy(x7.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbmv(f64, .ColumnMajor, .Upper, .Trans, .NonUnit, n, k, A.ptr, k + 1, x7.ptr, 2);

    try std.testing.expectEqual(2, x7[0]);
    try std.testing.expectEqual(11, x7[2]);
    try std.testing.expectEqual(28, x7[4]);
    try std.testing.expectEqual(53, x7[6]);
    try std.testing.expectEqual(86, x7[8]);

    const x8 = try a.alloc(f64, 2 * n);
    defer a.free(x8);

    @memcpy(x8.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbmv(f64, .ColumnMajor, .Upper, .Trans, .Unit, n, k, A.ptr, k + 1, x8.ptr, -2);

    try std.testing.expectEqual(1, x8[8]);
    try std.testing.expectEqual(5, x8[6]);
    try std.testing.expectEqual(13, x8[4]);
    try std.testing.expectEqual(25, x8[2]);
    try std.testing.expectEqual(41, x8[0]);

    const x9 = try a.alloc(f64, 2 * n);
    defer a.free(x9);

    @memcpy(x9.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbmv(f64, .RowMajor, .Lower, .NoTrans, .NonUnit, n, k, A.ptr, k + 1, x9.ptr, 2);

    try std.testing.expectEqual(2, x9[0]);
    try std.testing.expectEqual(11, x9[2]);
    try std.testing.expectEqual(28, x9[4]);
    try std.testing.expectEqual(53, x9[6]);
    try std.testing.expectEqual(86, x9[8]);

    const x10 = try a.alloc(f64, 2 * n);
    defer a.free(x10);

    @memcpy(x10.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbmv(f64, .RowMajor, .Lower, .NoTrans, .Unit, n, k, A.ptr, k + 1, x10.ptr, -2);

    try std.testing.expectEqual(1, x10[8]);
    try std.testing.expectEqual(5, x10[6]);
    try std.testing.expectEqual(13, x10[4]);
    try std.testing.expectEqual(25, x10[2]);
    try std.testing.expectEqual(41, x10[0]);

    const x11 = try a.alloc(f64, 2 * n);
    defer a.free(x11);

    @memcpy(x11.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbmv(f64, .ColumnMajor, .Lower, .NoTrans, .NonUnit, n, k, A.ptr, k + 1, x11.ptr, 2);

    try std.testing.expectEqual(1, x11[0]);
    try std.testing.expectEqual(8, x11[2]);
    try std.testing.expectEqual(23, x11[4]);
    try std.testing.expectEqual(46, x11[6]);
    try std.testing.expectEqual(77, x11[8]);

    const x12 = try a.alloc(f64, 2 * n);
    defer a.free(x12);

    @memcpy(x12.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbmv(f64, .ColumnMajor, .Lower, .NoTrans, .Unit, n, k, A.ptr, k + 1, x12.ptr, -2);

    try std.testing.expectEqual(1, x12[8]);
    try std.testing.expectEqual(4, x12[6]);
    try std.testing.expectEqual(11, x12[4]);
    try std.testing.expectEqual(22, x12[2]);
    try std.testing.expectEqual(37, x12[0]);

    const x13 = try a.alloc(f64, 2 * n);
    defer a.free(x13);

    @memcpy(x13.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbmv(f64, .RowMajor, .Lower, .Trans, .NonUnit, n, k, A.ptr, k + 1, x13.ptr, 2);

    try std.testing.expectEqual(8, x13[0]);
    try std.testing.expectEqual(23, x13[2]);
    try std.testing.expectEqual(46, x13[4]);
    try std.testing.expectEqual(77, x13[6]);
    try std.testing.expectEqual(50, x13[8]);

    const x14 = try a.alloc(f64, 2 * n);
    defer a.free(x14);

    @memcpy(x14.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbmv(f64, .RowMajor, .Lower, .Trans, .Unit, n, k, A.ptr, k + 1, x14.ptr, -2);

    try std.testing.expectEqual(7, x14[8]);
    try std.testing.expectEqual(17, x14[6]);
    try std.testing.expectEqual(31, x14[4]);
    try std.testing.expectEqual(49, x14[2]);
    try std.testing.expectEqual(5, x14[0]);

    const x15 = try a.alloc(f64, 2 * n);
    defer a.free(x15);

    @memcpy(x15.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.tbmv(f64, .ColumnMajor, .Lower, .Trans, .NonUnit, n, k, A.ptr, k + 1, x15.ptr, 2);

    try std.testing.expectEqual(5, x15[0]);
    try std.testing.expectEqual(18, x15[2]);
    try std.testing.expectEqual(39, x15[4]);
    try std.testing.expectEqual(68, x15[6]);
    try std.testing.expectEqual(45, x15[8]);

    const x16 = try a.alloc(f64, 2 * n);
    defer a.free(x16);

    @memcpy(x16.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.tbmv(f64, .ColumnMajor, .Lower, .Trans, .Unit, n, k, A.ptr, k + 1, x16.ptr, -2);

    try std.testing.expectEqual(5, x16[8]);
    try std.testing.expectEqual(14, x16[6]);
    try std.testing.expectEqual(27, x16[4]);
    try std.testing.expectEqual(44, x16[2]);
    try std.testing.expectEqual(5, x16[0]);

    const B = try a.alloc(Complex(f64), n * n);
    defer a.free(B);
    const x17 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x17);

    @memcpy(B.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(2, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(4, 4),
        Complex(f64).init(5, 5),
        Complex(f64).init(6, 6),
        Complex(f64).init(7, 7),
        Complex(f64).init(8, 8),
        Complex(f64).init(9, 9),
        Complex(f64).init(10, 10),
        Complex(f64).init(11, 11),
        Complex(f64).init(12, 12),
        Complex(f64).init(13, 13),
        Complex(f64).init(14, 14),
        Complex(f64).init(15, 15),
        Complex(f64).init(16, 16),
        Complex(f64).init(17, 17),
        Complex(f64).init(18, 18),
        Complex(f64).init(19, 19),
        Complex(f64).init(20, 20),
        Complex(f64).init(21, 21),
        Complex(f64).init(22, 22),
        Complex(f64).init(23, 23),
        Complex(f64).init(24, 24),
        Complex(f64).init(25, 25),
    });
    @memcpy(x17.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Upper, .NoTrans, .NonUnit, n, k, B.ptr, k + 1, x17.ptr, 2);

    try std.testing.expectEqual(0, x17[0].re);
    try std.testing.expectEqual(10, x17[0].im);
    try std.testing.expectEqual(0, x17[2].re);
    try std.testing.expectEqual(36, x17[2].im);
    try std.testing.expectEqual(0, x17[4].re);
    try std.testing.expectEqual(78, x17[4].im);
    try std.testing.expectEqual(0, x17[6].re);
    try std.testing.expectEqual(136, x17[6].im);
    try std.testing.expectEqual(0, x17[8].re);
    try std.testing.expectEqual(90, x17[8].im);

    const x18 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x18);

    @memcpy(x18.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Upper, .NoTrans, .Unit, n, k, B.ptr, k + 1, x18.ptr, -2);

    try std.testing.expectEqual(1, x18[8].re);
    try std.testing.expectEqual(9, x18[8].im);
    try std.testing.expectEqual(2, x18[6].re);
    try std.testing.expectEqual(26, x18[6].im);
    try std.testing.expectEqual(3, x18[4].re);
    try std.testing.expectEqual(51, x18[4].im);
    try std.testing.expectEqual(4, x18[2].re);
    try std.testing.expectEqual(84, x18[2].im);
    try std.testing.expectEqual(5, x18[0].re);
    try std.testing.expectEqual(5, x18[0].im);

    const x19 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x19);

    @memcpy(x19.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Upper, .NoTrans, .NonUnit, n, k, B.ptr, k + 1, x19.ptr, 2);

    try std.testing.expectEqual(0, x19[0].re);
    try std.testing.expectEqual(16, x19[0].im);
    try std.testing.expectEqual(0, x19[2].re);
    try std.testing.expectEqual(46, x19[2].im);
    try std.testing.expectEqual(0, x19[4].re);
    try std.testing.expectEqual(92, x19[4].im);
    try std.testing.expectEqual(0, x19[6].re);
    try std.testing.expectEqual(154, x19[6].im);
    try std.testing.expectEqual(0, x19[8].re);
    try std.testing.expectEqual(100, x19[8].im);

    const x20 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x20);

    @memcpy(x20.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Upper, .NoTrans, .Unit, n, k, B.ptr, k + 1, x20.ptr, -2);

    try std.testing.expectEqual(1, x20[8].re);
    try std.testing.expectEqual(13, x20[8].im);
    try std.testing.expectEqual(2, x20[6].re);
    try std.testing.expectEqual(32, x20[6].im);
    try std.testing.expectEqual(3, x20[4].re);
    try std.testing.expectEqual(59, x20[4].im);
    try std.testing.expectEqual(4, x20[2].re);
    try std.testing.expectEqual(94, x20[2].im);
    try std.testing.expectEqual(5, x20[0].re);
    try std.testing.expectEqual(5, x20[0].im);

    const x21 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x21);

    @memcpy(x21.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Upper, .ConjNoTrans, .NonUnit, n, k, B.ptr, k + 1, x21.ptr, 2);

    try std.testing.expectEqual(10, x21[0].re);
    try std.testing.expectEqual(0, x21[0].im);
    try std.testing.expectEqual(36, x21[2].re);
    try std.testing.expectEqual(0, x21[2].im);
    try std.testing.expectEqual(78, x21[4].re);
    try std.testing.expectEqual(0, x21[4].im);
    try std.testing.expectEqual(136, x21[6].re);
    try std.testing.expectEqual(0, x21[6].im);
    try std.testing.expectEqual(90, x21[8].re);
    try std.testing.expectEqual(0, x21[8].im);

    const x22 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x22);

    @memcpy(x22.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Upper, .ConjNoTrans, .Unit, n, k, B.ptr, k + 1, x22.ptr, -2);

    try std.testing.expectEqual(9, x22[8].re);
    try std.testing.expectEqual(1, x22[8].im);
    try std.testing.expectEqual(26, x22[6].re);
    try std.testing.expectEqual(2, x22[6].im);
    try std.testing.expectEqual(51, x22[4].re);
    try std.testing.expectEqual(3, x22[4].im);
    try std.testing.expectEqual(84, x22[2].re);
    try std.testing.expectEqual(4, x22[2].im);
    try std.testing.expectEqual(5, x22[0].re);
    try std.testing.expectEqual(5, x22[0].im);

    const x23 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x23);

    @memcpy(x23.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Upper, .ConjNoTrans, .NonUnit, n, k, B.ptr, k + 1, x23.ptr, 2);

    try std.testing.expectEqual(16, x23[0].re);
    try std.testing.expectEqual(0, x23[0].im);
    try std.testing.expectEqual(46, x23[2].re);
    try std.testing.expectEqual(0, x23[2].im);
    try std.testing.expectEqual(92, x23[4].re);
    try std.testing.expectEqual(0, x23[4].im);
    try std.testing.expectEqual(154, x23[6].re);
    try std.testing.expectEqual(0, x23[6].im);
    try std.testing.expectEqual(100, x23[8].re);
    try std.testing.expectEqual(0, x23[8].im);

    const x24 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x24);

    @memcpy(x24.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Upper, .ConjNoTrans, .Unit, n, k, B.ptr, k + 1, x24.ptr, -2);

    try std.testing.expectEqual(13, x24[8].re);
    try std.testing.expectEqual(1, x24[8].im);
    try std.testing.expectEqual(32, x24[6].re);
    try std.testing.expectEqual(2, x24[6].im);
    try std.testing.expectEqual(59, x24[4].re);
    try std.testing.expectEqual(3, x24[4].im);
    try std.testing.expectEqual(94, x24[2].re);
    try std.testing.expectEqual(4, x24[2].im);
    try std.testing.expectEqual(5, x24[0].re);
    try std.testing.expectEqual(5, x24[0].im);

    const x25 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x25);

    @memcpy(x25.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Upper, .Trans, .NonUnit, n, k, B.ptr, k + 1, x25.ptr, 2);

    try std.testing.expectEqual(0, x25[0].re);
    try std.testing.expectEqual(2, x25[0].im);
    try std.testing.expectEqual(0, x25[2].re);
    try std.testing.expectEqual(16, x25[2].im);
    try std.testing.expectEqual(0, x25[4].re);
    try std.testing.expectEqual(46, x25[4].im);
    try std.testing.expectEqual(0, x25[6].re);
    try std.testing.expectEqual(92, x25[6].im);
    try std.testing.expectEqual(0, x25[8].re);
    try std.testing.expectEqual(154, x25[8].im);

    const x26 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x26);

    @memcpy(x26.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Upper, .Trans, .Unit, n, k, B.ptr, k + 1, x26.ptr, -2);

    try std.testing.expectEqual(1, x26[8].re);
    try std.testing.expectEqual(1, x26[8].im);
    try std.testing.expectEqual(2, x26[6].re);
    try std.testing.expectEqual(6, x26[6].im);
    try std.testing.expectEqual(3, x26[4].re);
    try std.testing.expectEqual(19, x26[4].im);
    try std.testing.expectEqual(4, x26[2].re);
    try std.testing.expectEqual(40, x26[2].im);
    try std.testing.expectEqual(5, x26[0].re);
    try std.testing.expectEqual(69, x26[0].im);

    const x27 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x27);

    @memcpy(x27.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Upper, .Trans, .NonUnit, n, k, B.ptr, k + 1, x27.ptr, 2);

    try std.testing.expectEqual(0, x27[0].re);
    try std.testing.expectEqual(4, x27[0].im);
    try std.testing.expectEqual(0, x27[2].re);
    try std.testing.expectEqual(22, x27[2].im);
    try std.testing.expectEqual(0, x27[4].re);
    try std.testing.expectEqual(56, x27[4].im);
    try std.testing.expectEqual(0, x27[6].re);
    try std.testing.expectEqual(106, x27[6].im);
    try std.testing.expectEqual(0, x27[8].re);
    try std.testing.expectEqual(172, x27[8].im);

    const x28 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x28);

    @memcpy(x28.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Upper, .Trans, .Unit, n, k, B.ptr, k + 1, x28.ptr, -2);

    try std.testing.expectEqual(1, x28[8].re);
    try std.testing.expectEqual(1, x28[8].im);
    try std.testing.expectEqual(2, x28[6].re);
    try std.testing.expectEqual(8, x28[6].im);
    try std.testing.expectEqual(3, x28[4].re);
    try std.testing.expectEqual(23, x28[4].im);
    try std.testing.expectEqual(4, x28[2].re);
    try std.testing.expectEqual(46, x28[2].im);
    try std.testing.expectEqual(5, x28[0].re);
    try std.testing.expectEqual(77, x28[0].im);

    const x29 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x29);

    @memcpy(x29.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Upper, .ConjTrans, .NonUnit, n, k, B.ptr, k + 1, x29.ptr, 2);

    try std.testing.expectEqual(2, x29[0].re);
    try std.testing.expectEqual(0, x29[0].im);
    try std.testing.expectEqual(16, x29[2].re);
    try std.testing.expectEqual(0, x29[2].im);
    try std.testing.expectEqual(46, x29[4].re);
    try std.testing.expectEqual(0, x29[4].im);
    try std.testing.expectEqual(92, x29[6].re);
    try std.testing.expectEqual(0, x29[6].im);
    try std.testing.expectEqual(154, x29[8].re);
    try std.testing.expectEqual(0, x29[8].im);

    const x30 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x30);

    @memcpy(x30.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Upper, .ConjTrans, .Unit, n, k, B.ptr, k + 1, x30.ptr, -2);

    try std.testing.expectEqual(1, x30[8].re);
    try std.testing.expectEqual(1, x30[8].im);
    try std.testing.expectEqual(6, x30[6].re);
    try std.testing.expectEqual(2, x30[6].im);
    try std.testing.expectEqual(19, x30[4].re);
    try std.testing.expectEqual(3, x30[4].im);
    try std.testing.expectEqual(40, x30[2].re);
    try std.testing.expectEqual(4, x30[2].im);
    try std.testing.expectEqual(69, x30[0].re);
    try std.testing.expectEqual(5, x30[0].im);

    const x31 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x31);

    @memcpy(x31.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Upper, .ConjTrans, .NonUnit, n, k, B.ptr, k + 1, x31.ptr, 2);

    try std.testing.expectEqual(4, x31[0].re);
    try std.testing.expectEqual(0, x31[0].im);
    try std.testing.expectEqual(22, x31[2].re);
    try std.testing.expectEqual(0, x31[2].im);
    try std.testing.expectEqual(56, x31[4].re);
    try std.testing.expectEqual(0, x31[4].im);
    try std.testing.expectEqual(106, x31[6].re);
    try std.testing.expectEqual(0, x31[6].im);
    try std.testing.expectEqual(172, x31[8].re);
    try std.testing.expectEqual(0, x31[8].im);

    const x32 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x32);

    @memcpy(x32.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Upper, .ConjTrans, .Unit, n, k, B.ptr, k + 1, x32.ptr, -2);

    try std.testing.expectEqual(1, x32[8].re);
    try std.testing.expectEqual(1, x32[8].im);
    try std.testing.expectEqual(8, x32[6].re);
    try std.testing.expectEqual(2, x32[6].im);
    try std.testing.expectEqual(23, x32[4].re);
    try std.testing.expectEqual(3, x32[4].im);
    try std.testing.expectEqual(46, x32[2].re);
    try std.testing.expectEqual(4, x32[2].im);
    try std.testing.expectEqual(77, x32[0].re);
    try std.testing.expectEqual(5, x32[0].im);

    const x33 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x33);

    @memcpy(x33.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Lower, .NoTrans, .NonUnit, n, k, B.ptr, k + 1, x33.ptr, 2);

    try std.testing.expectEqual(0, x33[0].re);
    try std.testing.expectEqual(4, x33[0].im);
    try std.testing.expectEqual(0, x33[2].re);
    try std.testing.expectEqual(22, x33[2].im);
    try std.testing.expectEqual(0, x33[4].re);
    try std.testing.expectEqual(56, x33[4].im);
    try std.testing.expectEqual(0, x33[6].re);
    try std.testing.expectEqual(106, x33[6].im);
    try std.testing.expectEqual(0, x33[8].re);
    try std.testing.expectEqual(172, x33[8].im);

    const x34 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x34);

    @memcpy(x34.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Lower, .NoTrans, .Unit, n, k, B.ptr, k + 1, x34.ptr, -2);

    try std.testing.expectEqual(1, x34[8].re);
    try std.testing.expectEqual(1, x34[8].im);
    try std.testing.expectEqual(2, x34[6].re);
    try std.testing.expectEqual(8, x34[6].im);
    try std.testing.expectEqual(3, x34[4].re);
    try std.testing.expectEqual(23, x34[4].im);
    try std.testing.expectEqual(4, x34[2].re);
    try std.testing.expectEqual(46, x34[2].im);
    try std.testing.expectEqual(5, x34[0].re);
    try std.testing.expectEqual(77, x34[0].im);

    const x35 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x35);

    @memcpy(x35.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Lower, .NoTrans, .NonUnit, n, k, B.ptr, k + 1, x35.ptr, 2);

    try std.testing.expectEqual(0, x35[0].re);
    try std.testing.expectEqual(2, x35[0].im);
    try std.testing.expectEqual(0, x35[2].re);
    try std.testing.expectEqual(16, x35[2].im);
    try std.testing.expectEqual(0, x35[4].re);
    try std.testing.expectEqual(46, x35[4].im);
    try std.testing.expectEqual(0, x35[6].re);
    try std.testing.expectEqual(92, x35[6].im);
    try std.testing.expectEqual(0, x35[8].re);
    try std.testing.expectEqual(154, x35[8].im);

    const x36 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x36);

    @memcpy(x36.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Lower, .NoTrans, .Unit, n, k, B.ptr, k + 1, x36.ptr, -2);

    try std.testing.expectEqual(1, x36[8].re);
    try std.testing.expectEqual(1, x36[8].im);
    try std.testing.expectEqual(2, x36[6].re);
    try std.testing.expectEqual(6, x36[6].im);
    try std.testing.expectEqual(3, x36[4].re);
    try std.testing.expectEqual(19, x36[4].im);
    try std.testing.expectEqual(4, x36[2].re);
    try std.testing.expectEqual(40, x36[2].im);
    try std.testing.expectEqual(5, x36[0].re);
    try std.testing.expectEqual(69, x36[0].im);

    const x37 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x37);

    @memcpy(x37.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Lower, .ConjNoTrans, .NonUnit, n, k, B.ptr, k + 1, x37.ptr, 2);

    try std.testing.expectEqual(4, x37[0].re);
    try std.testing.expectEqual(0, x37[0].im);
    try std.testing.expectEqual(22, x37[2].re);
    try std.testing.expectEqual(0, x37[2].im);
    try std.testing.expectEqual(56, x37[4].re);
    try std.testing.expectEqual(0, x37[4].im);
    try std.testing.expectEqual(106, x37[6].re);
    try std.testing.expectEqual(0, x37[6].im);
    try std.testing.expectEqual(172, x37[8].re);
    try std.testing.expectEqual(0, x37[8].im);

    const x38 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x38);

    @memcpy(x38.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Lower, .ConjNoTrans, .Unit, n, k, B.ptr, k + 1, x38.ptr, -2);

    try std.testing.expectEqual(1, x38[8].re);
    try std.testing.expectEqual(1, x38[8].im);
    try std.testing.expectEqual(8, x38[6].re);
    try std.testing.expectEqual(2, x38[6].im);
    try std.testing.expectEqual(23, x38[4].re);
    try std.testing.expectEqual(3, x38[4].im);
    try std.testing.expectEqual(46, x38[2].re);
    try std.testing.expectEqual(4, x38[2].im);
    try std.testing.expectEqual(77, x38[0].re);
    try std.testing.expectEqual(5, x38[0].im);

    const x39 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x39);

    @memcpy(x39.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Lower, .ConjNoTrans, .NonUnit, n, k, B.ptr, k + 1, x39.ptr, 2);

    try std.testing.expectEqual(2, x39[0].re);
    try std.testing.expectEqual(0, x39[0].im);
    try std.testing.expectEqual(16, x39[2].re);
    try std.testing.expectEqual(0, x39[2].im);
    try std.testing.expectEqual(46, x39[4].re);
    try std.testing.expectEqual(0, x39[4].im);
    try std.testing.expectEqual(92, x39[6].re);
    try std.testing.expectEqual(0, x39[6].im);
    try std.testing.expectEqual(154, x39[8].re);
    try std.testing.expectEqual(0, x39[8].im);

    const x40 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x40);

    @memcpy(x40.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Lower, .ConjNoTrans, .Unit, n, k, B.ptr, k + 1, x40.ptr, -2);

    try std.testing.expectEqual(1, x40[8].re);
    try std.testing.expectEqual(1, x40[8].im);
    try std.testing.expectEqual(6, x40[6].re);
    try std.testing.expectEqual(2, x40[6].im);
    try std.testing.expectEqual(19, x40[4].re);
    try std.testing.expectEqual(3, x40[4].im);
    try std.testing.expectEqual(40, x40[2].re);
    try std.testing.expectEqual(4, x40[2].im);
    try std.testing.expectEqual(69, x40[0].re);
    try std.testing.expectEqual(5, x40[0].im);

    const x41 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x41);

    @memcpy(x41.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Lower, .Trans, .NonUnit, n, k, B.ptr, k + 1, x41.ptr, 2);

    try std.testing.expectEqual(0, x41[0].re);
    try std.testing.expectEqual(16, x41[0].im);
    try std.testing.expectEqual(0, x41[2].re);
    try std.testing.expectEqual(46, x41[2].im);
    try std.testing.expectEqual(0, x41[4].re);
    try std.testing.expectEqual(92, x41[4].im);
    try std.testing.expectEqual(0, x41[6].re);
    try std.testing.expectEqual(154, x41[6].im);
    try std.testing.expectEqual(0, x41[8].re);
    try std.testing.expectEqual(100, x41[8].im);

    const x42 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x42);

    @memcpy(x42.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Lower, .Trans, .Unit, n, k, B.ptr, k + 1, x42.ptr, -2);

    try std.testing.expectEqual(1, x42[8].re);
    try std.testing.expectEqual(13, x42[8].im);
    try std.testing.expectEqual(2, x42[6].re);
    try std.testing.expectEqual(32, x42[6].im);
    try std.testing.expectEqual(3, x42[4].re);
    try std.testing.expectEqual(59, x42[4].im);
    try std.testing.expectEqual(4, x42[2].re);
    try std.testing.expectEqual(94, x42[2].im);
    try std.testing.expectEqual(5, x42[0].re);
    try std.testing.expectEqual(5, x42[0].im);

    const x43 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x43);

    @memcpy(x43.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Lower, .Trans, .NonUnit, n, k, B.ptr, k + 1, x43.ptr, 2);

    try std.testing.expectEqual(0, x43[0].re);
    try std.testing.expectEqual(10, x43[0].im);
    try std.testing.expectEqual(0, x43[2].re);
    try std.testing.expectEqual(36, x43[2].im);
    try std.testing.expectEqual(0, x43[4].re);
    try std.testing.expectEqual(78, x43[4].im);
    try std.testing.expectEqual(0, x43[6].re);
    try std.testing.expectEqual(136, x43[6].im);
    try std.testing.expectEqual(0, x43[8].re);
    try std.testing.expectEqual(90, x43[8].im);

    const x44 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x44);

    @memcpy(x44.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Lower, .Trans, .Unit, n, k, B.ptr, k + 1, x44.ptr, -2);

    try std.testing.expectEqual(1, x44[8].re);
    try std.testing.expectEqual(9, x44[8].im);
    try std.testing.expectEqual(2, x44[6].re);
    try std.testing.expectEqual(26, x44[6].im);
    try std.testing.expectEqual(3, x44[4].re);
    try std.testing.expectEqual(51, x44[4].im);
    try std.testing.expectEqual(4, x44[2].re);
    try std.testing.expectEqual(84, x44[2].im);
    try std.testing.expectEqual(5, x44[0].re);
    try std.testing.expectEqual(5, x44[0].im);

    const x45 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x45);

    @memcpy(x45.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Lower, .ConjTrans, .NonUnit, n, k, B.ptr, k + 1, x45.ptr, 2);

    try std.testing.expectEqual(16, x45[0].re);
    try std.testing.expectEqual(0, x45[0].im);
    try std.testing.expectEqual(46, x45[2].re);
    try std.testing.expectEqual(0, x45[2].im);
    try std.testing.expectEqual(92, x45[4].re);
    try std.testing.expectEqual(0, x45[4].im);
    try std.testing.expectEqual(154, x45[6].re);
    try std.testing.expectEqual(0, x45[6].im);
    try std.testing.expectEqual(100, x45[8].re);
    try std.testing.expectEqual(0, x45[8].im);

    const x46 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x46);

    @memcpy(x46.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .RowMajor, .Lower, .ConjTrans, .Unit, n, k, B.ptr, k + 1, x46.ptr, -2);

    try std.testing.expectEqual(13, x46[8].re);
    try std.testing.expectEqual(1, x46[8].im);
    try std.testing.expectEqual(32, x46[6].re);
    try std.testing.expectEqual(2, x46[6].im);
    try std.testing.expectEqual(59, x46[4].re);
    try std.testing.expectEqual(3, x46[4].im);
    try std.testing.expectEqual(94, x46[2].re);
    try std.testing.expectEqual(4, x46[2].im);
    try std.testing.expectEqual(5, x46[0].re);
    try std.testing.expectEqual(5, x46[0].im);

    const x47 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x47);

    @memcpy(x47.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Lower, .ConjTrans, .NonUnit, n, k, B.ptr, k + 1, x47.ptr, 2);

    try std.testing.expectEqual(10, x47[0].re);
    try std.testing.expectEqual(0, x47[0].im);
    try std.testing.expectEqual(36, x47[2].re);
    try std.testing.expectEqual(0, x47[2].im);
    try std.testing.expectEqual(78, x47[4].re);
    try std.testing.expectEqual(0, x47[4].im);
    try std.testing.expectEqual(136, x47[6].re);
    try std.testing.expectEqual(0, x47[6].im);
    try std.testing.expectEqual(90, x47[8].re);
    try std.testing.expectEqual(0, x47[8].im);

    const x48 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x48);

    @memcpy(x48.ptr, &[_]Complex(f64){
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    BLAS.tbmv(Complex(f64), .ColumnMajor, .Lower, .ConjTrans, .Unit, n, k, B.ptr, k + 1, x48.ptr, -2);

    try std.testing.expectEqual(9, x48[8].re);
    try std.testing.expectEqual(1, x48[8].im);
    try std.testing.expectEqual(26, x48[6].re);
    try std.testing.expectEqual(2, x48[6].im);
    try std.testing.expectEqual(51, x48[4].re);
    try std.testing.expectEqual(3, x48[4].im);
    try std.testing.expectEqual(84, x48[2].re);
    try std.testing.expectEqual(4, x48[2].im);
    try std.testing.expectEqual(5, x48[0].re);
    try std.testing.expectEqual(5, x48[0].im);
}
