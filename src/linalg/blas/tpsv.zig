const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Transpose = blas.Transpose;
const Order = blas.Order;
const Uplo = blas.Uplo;
const Diag = blas.Diag;

pub inline fn tpsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const T, x: [*]T, incx: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    var TRANSA = transA;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        TRANSA = if (transA == .NoTrans) .Trans else if (transA == .ConjNoTrans) .ConjTrans else if (transA == .Trans) .NoTrans else .ConjNoTrans;
    }

    const LENX = N;

    switch (numericType) {
        .bool => @compileError("blas.tpsv does not support bool."),
        .int, .float => {
            if (UPLO == .Upper) {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = ((N - 1) * N) >> 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = N;
                        while (j >= 0) {
                            x[@intCast(jx)] /= Ap[@intCast(j + jaj)];
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                x[@intCast(ix)] -= t0 * Ap[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            lda -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = ((N - 1) * N) >> 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = N;
                        while (j >= 0) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                x[@intCast(ix)] -= t0 * Ap[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            lda -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = 1;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0 -= Ap[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0 / Ap[@intCast(iaij)];

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda += 1;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = 1;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0 -= Ap[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda += 1;
                        }
                    }
                }
            } else {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = N;
                        while (j < N) {
                            x[@intCast(jx)] /= Ap[@intCast(jaj)];
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)] -= t0 * Ap[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda -= 1;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = N;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)] -= t0 * Ap[@intCast(iaij)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda -= 1;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = (((N + 1) * N) >> 1) - 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = 1;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0 -= Ap[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0 / Ap[@intCast(jaj)];

                            j -= 1;
                            lda += 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = (((N + 1) * N) >> 1) - 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = 1;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0 -= Ap[@intCast(iaij)] * x[@intCast(ix)];

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j -= 1;
                            lda += 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                }
            }
        },
        .cfloat => {
            if (UPLO == .Upper) {
                if (TRANSA == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = ((N - 1) * N) >> 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = N;
                        while (j >= 0) {
                            var temp: T = undefined;
                            if (@abs(Ap[@intCast(j + jaj)].im) < @abs(Ap[@intCast(j + jaj)].re)) {
                                const temp1 = Ap[@intCast(j + jaj)].im / Ap[@intCast(j + jaj)].re;
                                const temp2 = Ap[@intCast(j + jaj)].re + temp1 * Ap[@intCast(j + jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = Ap[@intCast(j + jaj)].re / Ap[@intCast(j + jaj)].im;
                                const temp2 = Ap[@intCast(j + jaj)].im + temp1 * Ap[@intCast(j + jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                x[@intCast(ix)].re -= Ap[@intCast(iaij)].re * t0.re - Ap[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= Ap[@intCast(iaij)].re * t0.im + Ap[@intCast(iaij)].im * t0.re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            lda -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = ((N - 1) * N) >> 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = N;
                        while (j >= 0) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                x[@intCast(ix)].re -= Ap[@intCast(iaij)].re * t0.re - Ap[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= Ap[@intCast(iaij)].re * t0.im + Ap[@intCast(iaij)].im * t0.re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            lda -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = ((N - 1) * N) >> 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = N;
                        while (j >= 0) {
                            var temp: T = undefined;
                            if (@abs(Ap[@intCast(j + jaj)].im) < @abs(Ap[@intCast(j + jaj)].re)) {
                                const temp1 = -Ap[@intCast(j + jaj)].im / Ap[@intCast(j + jaj)].re;
                                const temp2 = Ap[@intCast(j + jaj)].re - temp1 * Ap[@intCast(j + jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = -Ap[@intCast(j + jaj)].re / Ap[@intCast(j + jaj)].im;
                                const temp2 = -Ap[@intCast(j + jaj)].im + temp1 * Ap[@intCast(j + jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                x[@intCast(ix)].re -= Ap[@intCast(iaij)].re * t0.re + Ap[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= Ap[@intCast(iaij)].re * t0.im - Ap[@intCast(iaij)].im * t0.re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            lda -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = ((N - 1) * N) >> 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = N;
                        while (j >= 0) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                x[@intCast(ix)].re -= Ap[@intCast(iaij)].re * t0.re + Ap[@intCast(iaij)].im * t0.im;
                                x[@intCast(ix)].im -= Ap[@intCast(iaij)].re * t0.im - Ap[@intCast(iaij)].im * t0.re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j -= 1;
                            lda -= 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else if (TRANSA == .Trans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = 1;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0.re -= Ap[@intCast(iaij)].re * x[@intCast(ix)].re - Ap[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= Ap[@intCast(iaij)].re * x[@intCast(ix)].im + Ap[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(Ap[@intCast(iaij)].im) < @abs(Ap[@intCast(iaij)].re)) {
                                const temp1 = Ap[@intCast(iaij)].im / Ap[@intCast(iaij)].re;
                                const temp2 = Ap[@intCast(iaij)].re + temp1 * Ap[@intCast(iaij)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = Ap[@intCast(iaij)].re / Ap[@intCast(iaij)].im;
                                const temp2 = Ap[@intCast(iaij)].im + temp1 * Ap[@intCast(iaij)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda += 1;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = 1;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0.re -= Ap[@intCast(iaij)].re * x[@intCast(ix)].re - Ap[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= Ap[@intCast(iaij)].re * x[@intCast(ix)].im + Ap[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda += 1;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = 1;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0.re -= Ap[@intCast(iaij)].re * x[@intCast(ix)].re + Ap[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= Ap[@intCast(iaij)].re * x[@intCast(ix)].im - Ap[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(Ap[@intCast(iaij)].im) < @abs(Ap[@intCast(iaij)].re)) {
                                const temp1 = -Ap[@intCast(iaij)].im / Ap[@intCast(iaij)].re;
                                const temp2 = Ap[@intCast(iaij)].re - temp1 * Ap[@intCast(iaij)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = -Ap[@intCast(iaij)].re / Ap[@intCast(iaij)].im;
                                const temp2 = -Ap[@intCast(iaij)].im + temp1 * Ap[@intCast(iaij)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda += 1;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = 1;
                        while (j < N) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = 0;
                            var iaij: isize = jaj;
                            var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                            while (i < j) {
                                t0.re -= Ap[@intCast(iaij)].re * x[@intCast(ix)].re + Ap[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= Ap[@intCast(iaij)].re * x[@intCast(ix)].im - Ap[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda += 1;
                        }
                    }
                }
            } else {
                if (TRANSA == .NoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = N;
                        while (j < N) {
                            var temp: T = undefined;
                            if (@abs(Ap[@intCast(jaj)].im) < @abs(Ap[@intCast(jaj)].re)) {
                                const temp1 = Ap[@intCast(jaj)].im / Ap[@intCast(jaj)].re;
                                const temp2 = Ap[@intCast(jaj)].re + temp1 * Ap[@intCast(jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = Ap[@intCast(jaj)].re / Ap[@intCast(jaj)].im;
                                const temp2 = Ap[@intCast(jaj)].im + temp1 * Ap[@intCast(jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)].re -= t0.re * Ap[@intCast(iaij)].re - t0.im * Ap[@intCast(iaij)].im;
                                x[@intCast(ix)].im -= t0.im * Ap[@intCast(iaij)].re + t0.re * Ap[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda -= 1;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = N;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)].re -= t0.re * Ap[@intCast(iaij)].re - t0.im * Ap[@intCast(iaij)].im;
                                x[@intCast(ix)].im -= t0.im * Ap[@intCast(iaij)].re + t0.re * Ap[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda -= 1;
                        }
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    if (diag == .NonUnit) {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = N;
                        while (j < N) {
                            var temp: T = undefined;
                            if (@abs(Ap[@intCast(jaj)].im) < @abs(Ap[@intCast(jaj)].re)) {
                                const temp1 = -Ap[@intCast(jaj)].im / Ap[@intCast(jaj)].re;
                                const temp2 = Ap[@intCast(jaj)].re - temp1 * Ap[@intCast(jaj)].im;
                                temp.re = (x[@intCast(jx)].re + temp1 * x[@intCast(jx)].im) / temp2;
                                temp.im = (x[@intCast(jx)].im - temp1 * x[@intCast(jx)].re) / temp2;
                            } else {
                                const temp1 = -Ap[@intCast(jaj)].re / Ap[@intCast(jaj)].im;
                                const temp2 = -Ap[@intCast(jaj)].im + temp1 * Ap[@intCast(jaj)].re;
                                temp.re = (temp1 * x[@intCast(jx)].re + x[@intCast(jx)].im) / temp2;
                                temp.im = (temp1 * x[@intCast(jx)].im - x[@intCast(jx)].re) / temp2;
                            }
                            x[@intCast(jx)] = temp;
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)].re -= t0.re * Ap[@intCast(iaij)].re + t0.im * Ap[@intCast(iaij)].im;
                                x[@intCast(ix)].im -= t0.im * Ap[@intCast(iaij)].re - t0.re * Ap[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda -= 1;
                        }
                    } else {
                        var j: isize = 0;
                        var jaj: isize = 0;
                        var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                        var lda: isize = N;
                        while (j < N) {
                            const t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                x[@intCast(ix)].re -= t0.re * Ap[@intCast(iaij)].re + t0.im * Ap[@intCast(iaij)].im;
                                x[@intCast(ix)].im -= t0.im * Ap[@intCast(iaij)].re - t0.re * Ap[@intCast(iaij)].im;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            j += 1;
                            jaj += lda;
                            jx += incx;
                            lda -= 1;
                        }
                    }
                } else if (TRANSA == .Trans) {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = (((N + 1) * N) >> 1) - 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = 1;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0.re -= Ap[@intCast(iaij)].re * x[@intCast(ix)].re - Ap[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= Ap[@intCast(iaij)].re * x[@intCast(ix)].im + Ap[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(Ap[@intCast(jaj)].im) < @abs(Ap[@intCast(jaj)].re)) {
                                const temp1 = Ap[@intCast(jaj)].im / Ap[@intCast(jaj)].re;
                                const temp2 = Ap[@intCast(jaj)].re + temp1 * Ap[@intCast(jaj)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = Ap[@intCast(jaj)].re / Ap[@intCast(jaj)].im;
                                const temp2 = Ap[@intCast(jaj)].im + temp1 * Ap[@intCast(jaj)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

                            j -= 1;
                            lda += 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = (((N + 1) * N) >> 1) - 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = 1;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0.re -= Ap[@intCast(iaij)].re * x[@intCast(ix)].re - Ap[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= Ap[@intCast(iaij)].re * x[@intCast(ix)].im + Ap[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j -= 1;
                            lda += 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                } else {
                    if (diag == .NonUnit) {
                        var j: isize = N - 1;
                        var jaj: isize = (((N + 1) * N) >> 1) - 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = 1;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0.re -= Ap[@intCast(iaij)].re * x[@intCast(ix)].re + Ap[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= Ap[@intCast(iaij)].re * x[@intCast(ix)].im - Ap[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            var temp: T = undefined;
                            if (@abs(Ap[@intCast(jaj)].im) < @abs(Ap[@intCast(jaj)].re)) {
                                const temp1 = -Ap[@intCast(jaj)].im / Ap[@intCast(jaj)].re;
                                const temp2 = Ap[@intCast(jaj)].re - temp1 * Ap[@intCast(jaj)].im;
                                temp.re = (t0.re + temp1 * t0.im) / temp2;
                                temp.im = (t0.im - temp1 * t0.re) / temp2;
                            } else {
                                const temp1 = -Ap[@intCast(jaj)].re / Ap[@intCast(jaj)].im;
                                const temp2 = -Ap[@intCast(jaj)].im + temp1 * Ap[@intCast(jaj)].re;
                                temp.re = (temp1 * t0.re + t0.im) / temp2;
                                temp.im = (temp1 * t0.im - t0.re) / temp2;
                            }
                            x[@intCast(jx)] = temp;

                            j -= 1;
                            lda += 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    } else {
                        var j: isize = N - 1;
                        var jaj: isize = (((N + 1) * N) >> 1) - 1;
                        var jx: isize = if (incx < 0) 0 else (LENX - 1) * incx;
                        var lda: isize = 1;
                        while (j >= 0) {
                            var t0 = x[@intCast(jx)];

                            var i: isize = j + 1;
                            var iaij: isize = jaj + 1;
                            var ix: isize = jx + incx;
                            while (i < N) {
                                t0.re -= Ap[@intCast(iaij)].re * x[@intCast(ix)].re + Ap[@intCast(iaij)].im * x[@intCast(ix)].im;
                                t0.im -= Ap[@intCast(iaij)].re * x[@intCast(ix)].im - Ap[@intCast(iaij)].im * x[@intCast(ix)].re;

                                i += 1;
                                iaij += 1;
                                ix += incx;
                            }

                            x[@intCast(jx)] = t0;

                            j -= 1;
                            lda += 1;
                            jaj -= lda;
                            jx -= incx;
                        }
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.tpsv only supports simple types."),
        .unsupported => unreachable,
    }
}
