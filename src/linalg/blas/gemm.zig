const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Transpose = blas.Transpose;

pub inline fn gemm(comptime T: type, order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (m <= 0 or n <= 0) return;

    var M = m;
    var N = n;
    var TRANSA = transA;
    var TRANSB = transB;
    var AA = A;
    var BB = B;
    var LDA = lda;
    var LDB = ldb;
    if (order == .RowMajor) {
        M = n;
        N = m;
        TRANSA = transB;
        TRANSB = transA;
        AA = B;
        BB = A;
        LDA = ldb;
        LDB = lda;
    }

    if ((TRANSA == .NoTrans or TRANSA == .ConjNoTrans) and LDA < @max(1, M)) return;
    if ((TRANSA == .Trans or TRANSA == .ConjTrans) and LDA < @max(1, k)) return;
    if ((TRANSB == .NoTrans or TRANSB == .ConjNoTrans) and LDB < @max(1, k)) return;
    if ((TRANSB == .Trans or TRANSB == .ConjTrans) and LDB < @max(1, N)) return;
    if (ldc < @max(1, M)) return;

    switch (numericType) {
        .bool => @compileError("blas.gemm does not support bool."),
        .int, .float => {
            if ((alpha == 0 or k <= 0) and beta == 1) return;

            if (alpha == 0) {
                if (beta == 0) {
                    var j: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var ijcj: isize = jcj;
                        while (i < M) {
                            C[@intCast(ijcj)] = 0;

                            i += 1;
                            ijcj += 1;
                        }

                        j += 1;
                        jcj += ldc;
                    }
                } else if (beta != 1) {
                    var j: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var ijcj: isize = jcj;
                        while (i < M) {
                            C[@intCast(ijcj)] *= beta;

                            i += 1;
                            ijcj += 1;
                        }

                        j += 1;
                        jcj += ldc;
                    }
                }

                return;
            }

            if (TRANSB == .NoTrans or TRANSB == .ConjNoTrans) {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = 0;

                                icij += 1;
                            }
                        } else if (beta != 1) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] *= beta;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var iblj: isize = jbj;
                        while (l < k) {
                            const t0 = alpha * BB[@intCast(iblj)];

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)] += t0 * AA[@intCast(iail)];

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            iblj += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var iai: isize = 0;
                        while (i < M) {
                            var t0: T = 0;

                            var l: isize = 0;
                            var iail: isize = iai;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0 += AA[@intCast(iail)] * BB[@intCast(iblj)];

                                l += 1;
                                iail += 1;
                                iblj += 1;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += alpha * t0;

                            i += 1;
                            iai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                }
            } else {
                if (TRANSA == .NoTrans or TRANSA == .ConjNoTrans) {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = 0;

                                icij += 1;
                            }
                        } else if (beta != 1) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] *= beta;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var ibjl: isize = ibj;
                        while (l < k) {
                            const t0 = alpha * BB[@intCast(ibjl)];

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)] += t0 * AA[@intCast(iail)];

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            ibjl += LDB;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var jai: isize = 0;
                        while (i < M) {
                            var t0: T = 0;

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ibjl: isize = ibj;
                            while (l < k) {
                                t0 += AA[@intCast(iali)] * BB[@intCast(ibjl)];

                                l += 1;
                                iali += 1;
                                ibjl += LDB;
                            }

                            if (beta == 0) {
                                C[@intCast(icij)] = 0;
                            } else if (beta != 1) {
                                C[@intCast(icij)] *= beta;
                            }

                            C[@intCast(icij)] += alpha * t0;

                            i += 1;
                            jai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                }
            }
        },
        .cfloat => {
            if (((alpha.re == 0 and alpha.im == 0) or k <= 0) and (beta.re == 1 and beta.im == 0)) return;

            if (alpha.re == 0 and alpha.im == 0) {
                if (beta.re == 0 and beta.im == 0) {
                    var j: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var ijcj: isize = jcj;
                        while (i < M) {
                            C[@intCast(ijcj)] = T.init(0, 0);

                            i += 1;
                            ijcj += 1;
                        }

                        j += 1;
                        jcj += ldc;
                    }
                } else if (beta.re != 1 or beta.im != 0) {
                    var j: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var ijcj: isize = jcj;
                        while (i < M) {
                            const tmp = C[@intCast(ijcj)].re * beta.re - C[@intCast(ijcj)].im * beta.im;
                            C[@intCast(ijcj)].im = C[@intCast(ijcj)].re * beta.im + C[@intCast(ijcj)].im * beta.re;
                            C[@intCast(ijcj)].re = tmp;

                            i += 1;
                            ijcj += 1;
                        }

                        j += 1;
                        jcj += ldc;
                    }
                }

                return;
            }

            if (TRANSB == .NoTrans) {
                if (TRANSA == .NoTrans) {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = T.init(0, 0);

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var iblj: isize = jbj;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re += alpha.re * BB[@intCast(iblj)].re - alpha.im * BB[@intCast(iblj)].im;
                            t0.im += alpha.re * BB[@intCast(iblj)].im + alpha.im * BB[@intCast(iblj)].re;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * AA[@intCast(iail)].re - t0.im * AA[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.re * AA[@intCast(iail)].im + t0.im * AA[@intCast(iail)].re;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            iblj += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = T.init(0, 0);

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var iblj: isize = jbj;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re = alpha.re * BB[@intCast(iblj)].re - alpha.im * BB[@intCast(iblj)].im;
                            t0.im = alpha.re * BB[@intCast(iblj)].im + alpha.im * BB[@intCast(iblj)].re;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * AA[@intCast(iail)].re + t0.im * AA[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.im * AA[@intCast(iail)].re - t0.re * AA[@intCast(iail)].im;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            iblj += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                } else if (TRANSA == .Trans) {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var iai: isize = 0;
                        while (i < M) {
                            var t0: T = T.init(0, 0);

                            var l: isize = 0;
                            var iail: isize = iai;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0.re += AA[@intCast(iail)].re * BB[@intCast(iblj)].re - AA[@intCast(iail)].im * BB[@intCast(iblj)].im;
                                t0.im += AA[@intCast(iail)].re * BB[@intCast(iblj)].im + AA[@intCast(iail)].im * BB[@intCast(iblj)].re;

                                l += 1;
                                iail += 1;
                                iblj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)] = T.init(0, 0);
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            iai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var iai: isize = 0;
                        while (i < M) {
                            var t0: T = T.init(0, 0);

                            var l: isize = 0;
                            var iail: isize = iai;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0.re += AA[@intCast(iail)].re * BB[@intCast(iblj)].re + AA[@intCast(iail)].im * BB[@intCast(iblj)].im;
                                t0.im += AA[@intCast(iail)].re * BB[@intCast(iblj)].im - AA[@intCast(iail)].im * BB[@intCast(iblj)].re;

                                l += 1;
                                iail += 1;
                                iblj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)] = T.init(0, 0);
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            iai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                }
            } else if (TRANSB == .ConjNoTrans) {
                if (TRANSA == .NoTrans) {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = T.init(0, 0);

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var iblj: isize = jbj;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re += alpha.re * BB[@intCast(iblj)].re + alpha.im * BB[@intCast(iblj)].im;
                            t0.im += alpha.im * BB[@intCast(iblj)].re - alpha.re * BB[@intCast(iblj)].im;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * AA[@intCast(iail)].re - t0.im * AA[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.re * AA[@intCast(iail)].im + t0.im * AA[@intCast(iail)].re;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            iblj += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = T.init(0, 0);

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var iblj: isize = jbj;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re += alpha.re * BB[@intCast(iblj)].re + alpha.im * BB[@intCast(iblj)].im;
                            t0.im += alpha.im * BB[@intCast(iblj)].re - alpha.re * BB[@intCast(iblj)].im;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * AA[@intCast(iail)].re + t0.im * AA[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.im * AA[@intCast(iail)].re - t0.re * AA[@intCast(iail)].im;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            iblj += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                } else if (TRANSA == .Trans) {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var iai: isize = 0;
                        while (i < M) {
                            var t0: T = T.init(0, 0);

                            var l: isize = 0;
                            var iail: isize = iai;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0.re += AA[@intCast(iail)].re * BB[@intCast(iblj)].re + AA[@intCast(iail)].im * BB[@intCast(iblj)].im;
                                t0.im += AA[@intCast(iail)].im * BB[@intCast(iblj)].re - AA[@intCast(iail)].re * BB[@intCast(iblj)].im;

                                l += 1;
                                iail += 1;
                                iblj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)] = T.init(0, 0);
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            iai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var jbj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var iai: isize = 0;
                        while (i < M) {
                            var t0: T = T.init(0, 0);

                            var l: isize = 0;
                            var iail: isize = iai;
                            var iblj: isize = jbj;
                            while (l < k) {
                                t0.re += AA[@intCast(iail)].re * BB[@intCast(iblj)].re - AA[@intCast(iail)].im * BB[@intCast(iblj)].im;
                                t0.im += -AA[@intCast(iail)].re * BB[@intCast(iblj)].im - AA[@intCast(iail)].im * BB[@intCast(iblj)].re;

                                l += 1;
                                iail += 1;
                                iblj += 1;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)] = T.init(0, 0);
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            iai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        jbj += LDB;
                        jcj += ldc;
                    }
                }
            } else if (TRANSB == .Trans) {
                if (TRANSA == .NoTrans) {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = T.init(0, 0);

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var ibjl: isize = ibj;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re += alpha.re * BB[@intCast(ibjl)].re - alpha.im * BB[@intCast(ibjl)].im;
                            t0.im += alpha.re * BB[@intCast(ibjl)].im + alpha.im * BB[@intCast(ibjl)].re;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * AA[@intCast(iail)].re - t0.im * AA[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.re * AA[@intCast(iail)].im + t0.im * AA[@intCast(iail)].re;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            ibjl += LDB;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = T.init(0, 0);

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var ibjl: isize = ibj;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re += alpha.re * BB[@intCast(ibjl)].re - alpha.im * BB[@intCast(ibjl)].im;
                            t0.im += alpha.re * BB[@intCast(ibjl)].im + alpha.im * BB[@intCast(ibjl)].re;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * AA[@intCast(iail)].re + t0.im * AA[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.im * AA[@intCast(iail)].re - t0.re * AA[@intCast(iail)].im;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            ibjl += LDB;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else if (TRANSA == .Trans) {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var jai: isize = 0;
                        while (i < M) {
                            var t0: T = T.init(0, 0);

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ibjl: isize = ibj;
                            while (l < k) {
                                t0.re += AA[@intCast(iali)].re * BB[@intCast(ibjl)].re - AA[@intCast(iali)].im * BB[@intCast(ibjl)].im;
                                t0.im += AA[@intCast(iali)].re * BB[@intCast(ibjl)].im + AA[@intCast(iali)].im * BB[@intCast(ibjl)].re;

                                l += 1;
                                iali += 1;
                                ibjl += LDB;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)] = T.init(0, 0);
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            jai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var jai: isize = 0;
                        while (i < M) {
                            var t0: T = T.init(0, 0);

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ibjl: isize = ibj;
                            while (l < k) {
                                t0.re += AA[@intCast(iali)].re * BB[@intCast(ibjl)].re + AA[@intCast(iali)].im * BB[@intCast(ibjl)].im;
                                t0.im += AA[@intCast(iali)].re * BB[@intCast(ibjl)].im - AA[@intCast(iali)].im * BB[@intCast(ibjl)].re;

                                l += 1;
                                iali += 1;
                                ibjl += LDB;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)] = T.init(0, 0);
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            jai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                }
            } else {
                if (TRANSA == .NoTrans) {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = T.init(0, 0);

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var ibjl: isize = ibj;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re += alpha.re * BB[@intCast(ibjl)].re + alpha.im * BB[@intCast(ibjl)].im;
                            t0.im += alpha.im * BB[@intCast(ibjl)].re - alpha.re * BB[@intCast(ibjl)].im;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * AA[@intCast(iail)].re - t0.im * AA[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.im * AA[@intCast(iail)].re + t0.re * AA[@intCast(iail)].im;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            ibjl += LDB;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else if (TRANSA == .ConjNoTrans) {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        if (beta.re == 0 and beta.im == 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                Cpjcj[@intCast(icij)] = T.init(0, 0);

                                icij += 1;
                            }
                        } else if (beta.re != 1 or beta.im != 0) {
                            var icij: isize = 0;
                            const Cpjcj: [*]T = @ptrCast(&C[@intCast(jcj)]);
                            while (icij < M) {
                                const tmp = Cpjcj[@intCast(icij)].re * beta.re - Cpjcj[@intCast(icij)].im * beta.im;
                                Cpjcj[@intCast(icij)].im = Cpjcj[@intCast(icij)].re * beta.im + Cpjcj[@intCast(icij)].im * beta.re;
                                Cpjcj[@intCast(icij)].re = tmp;

                                icij += 1;
                            }
                        }

                        var l: isize = 0;
                        var jal: isize = 0;
                        var ibjl: isize = ibj;
                        while (l < k) {
                            var t0: T = undefined;
                            t0.re += alpha.re * BB[@intCast(ibjl)].re + alpha.im * BB[@intCast(ibjl)].im;
                            t0.im += alpha.im * BB[@intCast(ibjl)].re - alpha.re * BB[@intCast(ibjl)].im;

                            var i: isize = 0;
                            var iail: isize = jal;
                            var icij: isize = jcj;
                            while (i < M) {
                                C[@intCast(icij)].re += t0.re * AA[@intCast(iail)].re + t0.im * AA[@intCast(iail)].im;
                                C[@intCast(icij)].im += t0.im * AA[@intCast(iail)].re - t0.re * AA[@intCast(iail)].im;

                                i += 1;
                                iail += 1;
                                icij += 1;
                            }

                            l += 1;
                            jal += LDA;
                            ibjl += LDB;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else if (TRANSA == .Trans) {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var jai: isize = 0;
                        while (i < M) {
                            var t0: T = T.init(0, 0);

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ibjl: isize = ibj;
                            while (l < k) {
                                t0.re += AA[@intCast(iali)].re * BB[@intCast(ibjl)].re + AA[@intCast(iali)].im * BB[@intCast(ibjl)].im;
                                t0.im += AA[@intCast(iali)].im * BB[@intCast(ibjl)].re - AA[@intCast(iali)].re * BB[@intCast(ibjl)].im;

                                l += 1;
                                iali += 1;
                                ibjl += LDB;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)] = T.init(0, 0);
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            jai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                } else {
                    var j: isize = 0;
                    var ibj: isize = 0;
                    var jcj: isize = 0;
                    while (j < N) {
                        var i: isize = 0;
                        var icij: isize = jcj;
                        var jai: isize = 0;
                        while (i < M) {
                            var t0: T = T.init(0, 0);

                            var l: isize = 0;
                            var iali: isize = jai;
                            var ibjl: isize = ibj;
                            while (l < k) {
                                t0.re += AA[@intCast(iali)].re * BB[@intCast(ibjl)].re - AA[@intCast(iali)].im * BB[@intCast(ibjl)].im;
                                t0.im += -AA[@intCast(iali)].im * BB[@intCast(ibjl)].re - AA[@intCast(iali)].re * BB[@intCast(ibjl)].im;

                                l += 1;
                                iali += 1;
                                ibjl += LDB;
                            }

                            if (beta.re == 0 and beta.im == 0) {
                                C[@intCast(icij)] = T.init(0, 0);
                            } else if (beta.re != 1 or beta.im != 0) {
                                const tmp = C[@intCast(icij)].re * beta.re - C[@intCast(icij)].im * beta.im;
                                C[@intCast(icij)].im = C[@intCast(icij)].re * beta.im + C[@intCast(icij)].im * beta.re;
                                C[@intCast(icij)].re = tmp;
                            }

                            C[@intCast(icij)].re += alpha.re * t0.re - alpha.im * t0.im;
                            C[@intCast(icij)].im += alpha.re * t0.im + alpha.im * t0.re;

                            i += 1;
                            jai += LDA;
                            icij += 1;
                        }

                        j += 1;
                        ibj += 1;
                        jcj += ldc;
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.gemm only supports simple types."),
        .unsupported => unreachable,
    }
}
