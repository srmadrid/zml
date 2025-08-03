const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = linalg.Order;
const Uplo = linalg.Uplo;

pub inline fn potrf2(
    order: Order,
    uplo: Uplo,
    n: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !isize {
    if (order == .col_major) {
        return k_potrf2_c(uplo, n, a, lda, ctx);
    } else {
        return k_potrf2_r(uplo, n, a, lda, ctx);
    }
}

fn k_potrf2_c(
    uplo: Uplo,
    n: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !isize {
    const A: type = types.Child(@TypeOf(a));

    if (n < 0 or lda < int.max(1, n))
        return lapack.Error.InvalidArgument;

    var info: isize = 0;

    // Quick return if possible.
    if (n == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(A)) {
        if (n == 1) {
            // Test for non-positive-definiteness
            const ajj: Scalar(A) = ops.re(a[0], ctx) catch unreachable;
            if (ops.le(ajj, 0, ctx) catch unreachable or
                std.math.isNan(ajj))
            {
                info = 1;

                return info;
            }

            // Factor.
            ops.sqrt_( // a[0] = sqrt(ajj)
                &a[0],
                ajj,
                ctx,
            ) catch unreachable;
        } else {
            // Use recursive code.
            const n1: isize = int.div(n, 2);
            const n2: isize = n - n1;

            // Factor A11.
            var iinfo: isize = k_potrf2_c(
                uplo,
                n1,
                a,
                lda,
                ctx,
            ) catch unreachable;

            if (iinfo != 0) {
                info = iinfo;

                return info;
            }

            if (uplo == .upper) {
                // Compute the Cholesky factorization A = U^T * U or A = U^H * U.

                // Update and scale A12.
                blas.trsm(
                    .col_major,
                    .left,
                    .upper,
                    .conj_trans,
                    .non_unit,
                    n1,
                    n2,
                    1,
                    a,
                    lda,
                    a + scast(usize, n1 * lda),
                    lda,
                    ctx,
                ) catch unreachable;

                // Update and factor A22.
                if (comptime !types.isComplex(A)) {
                    blas.syrk(
                        .col_major,
                        uplo,
                        .trans,
                        n2,
                        n1,
                        -1,
                        a + scast(usize, n1 * lda),
                        lda,
                        1,
                        a + scast(usize, n1 + n1 * lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                } else {
                    blas.herk(
                        .col_major,
                        uplo,
                        .conj_trans,
                        n2,
                        n1,
                        -1,
                        a + scast(usize, n1 * lda),
                        lda,
                        1,
                        a + scast(usize, n1 + n1 * lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                }

                iinfo = k_potrf2_c(
                    uplo,
                    n2,
                    a + scast(usize, n1 + n1 * lda),
                    lda,
                    ctx,
                ) catch unreachable;

                if (iinfo != 0) {
                    info = iinfo + n1;

                    return info;
                }
            } else {
                // Compute the Cholesky factorization A = L * L^T or A = L * L^H.

                // Update and scale A21.
                blas.trsm(
                    .col_major,
                    .right,
                    .lower,
                    .conj_trans,
                    .non_unit,
                    n2,
                    n1,
                    1,
                    a,
                    lda,
                    a + scast(usize, n1),
                    lda,
                    ctx,
                ) catch unreachable;

                // Update and factor A22.
                if (comptime !types.isComplex(A)) {
                    blas.syrk(
                        .col_major,
                        uplo,
                        .no_trans,
                        n2,
                        n1,
                        -1,
                        a + scast(usize, n1),
                        lda,
                        1,
                        a + scast(usize, n1 + n1 * lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                } else {
                    blas.herk(
                        .col_major,
                        uplo,
                        .no_trans,
                        n2,
                        n1,
                        -1,
                        a + scast(usize, n1),
                        lda,
                        1,
                        a + scast(usize, n1 + n1 * lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                }

                iinfo = k_potrf2_c(
                    uplo,
                    n2,
                    a + scast(usize, n1 + n1 * lda),
                    lda,
                    ctx,
                ) catch unreachable;

                if (iinfo != 0) {
                    info = iinfo + n1;

                    return info;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.potrf2 not implemented for arbitrary precision types yet");
    }

    return info;
}

fn k_potrf2_r(
    uplo: Uplo,
    n: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !isize {
    const A: type = types.Child(@TypeOf(a));

    if (n < 0 or lda < int.max(1, n))
        return lapack.Error.InvalidArgument;

    var info: isize = 0;

    // Quick return if possible.
    if (n == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(A)) {
        if (n == 1) {
            // Test for non-positive-definiteness
            const ajj: Scalar(A) = ops.re(a[0], ctx) catch unreachable;
            if (ops.le(ajj, 0, ctx) catch unreachable or
                std.math.isNan(ajj))
            {
                info = 1;

                return info;
            }

            // Factor.
            ops.sqrt_( // a[0] = sqrt(ajj)
                &a[0],
                ajj,
                ctx,
            ) catch unreachable;
        } else {
            // Use recursive code.
            const n1: isize = int.div(n, 2);
            const n2: isize = n - n1;

            // Factor A11.
            var iinfo: isize = k_potrf2_r(
                uplo,
                n1,
                a,
                lda,
                ctx,
            ) catch unreachable;

            if (iinfo != 0) {
                info = iinfo;

                return info;
            }

            if (uplo == .upper) {
                // Compute the Cholesky factorization A = U^T * U or A = U^H * U.

                // Update and scale A12.
                blas.trsm(
                    .row_major,
                    .left,
                    .upper,
                    .conj_trans,
                    .non_unit,
                    n1,
                    n2,
                    1,
                    a,
                    lda,
                    a + scast(usize, n1),
                    lda,
                    ctx,
                ) catch unreachable;

                // Update and factor A22.
                if (comptime !types.isComplex(A)) {
                    blas.syrk(
                        .row_major,
                        uplo,
                        .trans,
                        n2,
                        n1,
                        -1,
                        a + scast(usize, n1),
                        lda,
                        1,
                        a + scast(usize, n1 * lda + n1),
                        lda,
                        ctx,
                    ) catch unreachable;
                } else {
                    blas.herk(
                        .row_major,
                        uplo,
                        .conj_trans,
                        n2,
                        n1,
                        -1,
                        a + scast(usize, n1),
                        lda,
                        1,
                        a + scast(usize, n1 * lda + n1),
                        lda,
                        ctx,
                    ) catch unreachable;
                }

                iinfo = k_potrf2_r(
                    uplo,
                    n2,
                    a + scast(usize, n1 * lda + n1),
                    lda,
                    ctx,
                ) catch unreachable;

                if (iinfo != 0) {
                    info = iinfo + n1;

                    return info;
                }
            } else {
                // Compute the Cholesky factorization A = L * L^T or A = L * L^H.

                // Update and scale A21.
                blas.trsm(
                    .row_major,
                    .right,
                    .lower,
                    .conj_trans,
                    .non_unit,
                    n2,
                    n1,
                    1,
                    a,
                    lda,
                    a + scast(usize, n1 * lda),
                    lda,
                    ctx,
                ) catch unreachable;

                // Update and factor A22.
                if (comptime !types.isComplex(A)) {
                    blas.syrk(
                        .row_major,
                        uplo,
                        .no_trans,
                        n2,
                        n1,
                        -1,
                        a + scast(usize, n1 * lda),
                        lda,
                        1,
                        a + scast(usize, n1 * lda + n1),
                        lda,
                        ctx,
                    ) catch unreachable;
                } else {
                    blas.herk(
                        .row_major,
                        uplo,
                        .no_trans,
                        n2,
                        n1,
                        -1,
                        a + scast(usize, n1 * lda),
                        lda,
                        1,
                        a + scast(usize, n1 * lda + n1),
                        lda,
                        ctx,
                    ) catch unreachable;
                }

                iinfo = k_potrf2_r(
                    uplo,
                    n2,
                    a + scast(usize, n1 * lda + n1),
                    lda,
                    ctx,
                ) catch unreachable;

                if (iinfo != 0) {
                    info = iinfo + n1;

                    return info;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.potrf2 not implemented for arbitrary precision types yet");
    }

    return info;
}
