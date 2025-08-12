const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Uplo = types.Uplo;

pub inline fn potrf(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !i32 {
    if (order == .col_major) {
        return k_potrf_c(uplo, n, a, lda, ctx);
    } else {
        return k_potrf_r(uplo, n, a, lda, ctx);
    }
}

fn k_potrf_c(
    uplo: Uplo,
    n: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));

    if (n < 0 or lda < int.max(1, n))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    // Quick return if possible.
    if (n == 0)
        return info;

    // Determine the block size for this environment. Always returns the
    // same number for getrf regardless of 'S', 'D', 'C', or 'Z'.
    const nb: i32 = lapack.ilaenv(1, "DPOTRF", " ", n, -1, -1, -1);
    if (comptime !types.isArbitraryPrecision(A)) {
        if (nb <= 1 or nb >= n) {
            // Use unblocked code.
            info = lapack.potrf2(
                .col_major,
                uplo,
                n,
                a,
                lda,
                ctx,
            ) catch unreachable;
        } else {
            // Use blocked code.
            if (uplo == .upper) {
                // Compute the Cholesky factorization A = U^T * U or A = U^H * U.
                var j: i32 = 0;
                while (j < n) : (j += nb) {
                    // Update and factorize the current diagonal block and test for non-positive-definiteness.
                    const jb: i32 = int.min(nb, n - j);

                    if (comptime !types.isComplex(A)) {
                        blas.syrk(
                            .col_major,
                            .upper,
                            .trans,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j * lda),
                            lda,
                            1,
                            a + scast(u32, j + j * lda),
                            lda,
                            ctx,
                        ) catch unreachable;
                    } else {
                        blas.herk(
                            .col_major,
                            .upper,
                            .conj_trans,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j * lda),
                            lda,
                            1,
                            a + scast(u32, j + j * lda),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }

                    info = lapack.potrf2(
                        .col_major,
                        .upper,
                        jb,
                        a + scast(u32, j + j * lda),
                        lda,
                        ctx,
                    ) catch unreachable;

                    if (info != 0) {
                        info += j;
                        return info;
                    }

                    if (j + jb < n) {
                        // Compute the current block row.
                        blas.gemm(
                            .col_major,
                            .conj_trans,
                            .no_trans,
                            jb,
                            n - j - jb,
                            j,
                            -1,
                            a + scast(u32, j * lda),
                            lda,
                            a + scast(u32, (j + jb) * lda),
                            lda,
                            1,
                            a + scast(u32, j + (j + jb) * lda),
                            lda,
                            ctx,
                        ) catch unreachable;

                        blas.trsm(
                            .col_major,
                            .left,
                            .upper,
                            .conj_trans,
                            .non_unit,
                            jb,
                            n - j - jb,
                            1,
                            a + scast(u32, j + j * lda),
                            lda,
                            a + scast(u32, j + (j + jb) * lda),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                // Compute the Cholesky factorization A = L * L^T or A = L * L^H.
                var j: i32 = 0;
                while (j < n) : (j += nb) {
                    // Update and factorize the current diagonal block and test for non-positive-definiteness.
                    const jb: i32 = int.min(nb, n - j);

                    if (comptime !types.isComplex(A)) {
                        blas.syrk(
                            .col_major,
                            .lower,
                            .no_trans,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j),
                            lda,
                            1,
                            a + scast(u32, j + j * lda),
                            lda,
                            ctx,
                        ) catch unreachable;
                    } else {
                        blas.herk(
                            .col_major,
                            .lower,
                            .no_trans,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j),
                            lda,
                            1,
                            a + scast(u32, j + j * lda),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }

                    info = lapack.potrf2(
                        .col_major,
                        .lower,
                        jb,
                        a + scast(u32, j + j * lda),
                        lda,
                        ctx,
                    ) catch unreachable;

                    if (info != 0) {
                        info += j;
                        return info;
                    }

                    if (j + jb < n) {
                        // Compute the current block column.
                        blas.gemm(
                            .col_major,
                            .no_trans,
                            .conj_trans,
                            n - j - jb,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j + jb),
                            lda,
                            a + scast(u32, j),
                            lda,
                            1,
                            a + scast(u32, (j + jb) + j * lda),
                            lda,
                            ctx,
                        ) catch unreachable;

                        blas.trsm(
                            .col_major,
                            .right,
                            .lower,
                            .conj_trans,
                            .non_unit,
                            n - j - jb,
                            jb,
                            1,
                            a + scast(u32, j + j * lda),
                            lda,
                            a + scast(u32, (j + jb) + j * lda),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.potrf not implemented for arbitrary precision types yet");
    }

    return info;
}

fn k_potrf_r(
    uplo: Uplo,
    n: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));

    if (n < 0 or lda < int.max(1, n))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    // Quick return if possible.
    if (n == 0)
        return info;

    // Determine the block size for this environment. Always returns the
    // same number for getrf regardless of 'S', 'D', 'C', or 'Z'.
    const nb: i32 = lapack.ilaenv(1, "DPOTRF", " ", n, -1, -1, -1);
    if (comptime !types.isArbitraryPrecision(A)) {
        if (nb <= 1 or nb >= n) {
            // Use unblocked code.
            info = lapack.potrf2(
                .row_major,
                uplo,
                n,
                a,
                lda,
                ctx,
            ) catch unreachable;
        } else {
            // Use blocked code.
            if (uplo == .upper) {
                // Compute the Cholesky factorization A = U^T * U or A = U^H * U.
                var j: i32 = 0;
                while (j < n) : (j += nb) {
                    // Update and factorize the current diagonal block and test for non-positive-definiteness.
                    const jb: i32 = int.min(nb, n - j);

                    if (comptime !types.isComplex(A)) {
                        blas.syrk(
                            .row_major,
                            .upper,
                            .trans,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j),
                            lda,
                            1,
                            a + scast(u32, j * lda + j),
                            lda,
                            ctx,
                        ) catch unreachable;
                    } else {
                        blas.herk(
                            .row_major,
                            .upper,
                            .conj_trans,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j),
                            lda,
                            1,
                            a + scast(u32, j * lda + j),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }

                    info = lapack.potrf2(
                        .row_major,
                        .upper,
                        jb,
                        a + scast(u32, j * lda + j),
                        lda,
                        ctx,
                    ) catch unreachable;

                    if (info != 0) {
                        info += j;
                        return info;
                    }

                    if (j + jb < n) {
                        // Compute the current block row.
                        blas.gemm(
                            .row_major,
                            .conj_trans,
                            .no_trans,
                            jb,
                            n - j - jb,
                            j,
                            -1,
                            a + scast(u32, j),
                            lda,
                            a + scast(u32, (j + jb)),
                            lda,
                            1,
                            a + scast(u32, j * lda + (j + jb)),
                            lda,
                            ctx,
                        ) catch unreachable;

                        blas.trsm(
                            .row_major,
                            .left,
                            .upper,
                            .conj_trans,
                            .non_unit,
                            jb,
                            n - j - jb,
                            1,
                            a + scast(u32, j * lda + j),
                            lda,
                            a + scast(u32, j * lda + (j + jb)),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                // Compute the Cholesky factorization A = L * L^T or A = L * L^H.
                var j: i32 = 0;
                while (j < n) : (j += nb) {
                    // Update and factorize the current diagonal block and test for non-positive-definiteness.
                    const jb: i32 = int.min(nb, n - j);

                    if (comptime !types.isComplex(A)) {
                        blas.syrk(
                            .row_major,
                            .lower,
                            .no_trans,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j * lda),
                            lda,
                            1,
                            a + scast(u32, j * lda + j),
                            lda,
                            ctx,
                        ) catch unreachable;
                    } else {
                        blas.herk(
                            .row_major,
                            .lower,
                            .no_trans,
                            jb,
                            j,
                            -1,
                            a + scast(u32, j * lda),
                            lda,
                            1,
                            a + scast(u32, j * lda + j),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }

                    info = lapack.potrf2(
                        .row_major,
                        .lower,
                        jb,
                        a + scast(u32, j * lda + j),
                        lda,
                        ctx,
                    ) catch unreachable;

                    if (info != 0) {
                        info += j;
                        return info;
                    }

                    if (j + jb < n) {
                        // Compute the current block column.
                        blas.gemm(
                            .row_major,
                            .no_trans,
                            .conj_trans,
                            n - j - jb,
                            jb,
                            j,
                            -1,
                            a + scast(u32, (j + jb) * lda),
                            lda,
                            a + scast(u32, j * lda),
                            lda,
                            1,
                            a + scast(u32, (j + jb) * lda + j),
                            lda,
                            ctx,
                        ) catch unreachable;

                        blas.trsm(
                            .row_major,
                            .right,
                            .lower,
                            .conj_trans,
                            .non_unit,
                            n - j - jb,
                            jb,
                            1,
                            a + scast(u32, j * lda + j),
                            lda,
                            a + scast(u32, (j + jb) * lda + j),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.potrf not implemented for arbitrary precision types yet");
    }

    return info;
}
