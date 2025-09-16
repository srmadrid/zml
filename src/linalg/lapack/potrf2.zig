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

const utils = @import("../utils.zig");

pub inline fn potrf2(
    order: Order,
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
            const n1: i32 = int.div(n, 2);
            const n2: i32 = n - n1;

            // Factor A11.
            var iinfo: i32 = potrf2(
                order,
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
                    order,
                    .left,
                    .upper,
                    .conj_trans,
                    .non_unit,
                    n1,
                    n2,
                    1,
                    a,
                    lda,
                    a + utils.index(order, 0, n1, lda),
                    lda,
                    ctx,
                ) catch unreachable;

                // Update and factor A22.
                if (comptime !types.isComplex(A)) {
                    blas.syrk(
                        order,
                        uplo,
                        .trans,
                        n2,
                        n1,
                        -1,
                        a + utils.index(order, 0, n1, lda),
                        lda,
                        1,
                        a + utils.index(order, n1, n1, lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                } else {
                    blas.herk(
                        order,
                        uplo,
                        .conj_trans,
                        n2,
                        n1,
                        -1,
                        a + utils.index(order, 0, n1, lda),
                        lda,
                        1,
                        a + utils.index(order, n1, n1, lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                }

                iinfo = potrf2(
                    order,
                    uplo,
                    n2,
                    a + utils.index(order, n1, n1, lda),
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
                    order,
                    .right,
                    .lower,
                    .conj_trans,
                    .non_unit,
                    n2,
                    n1,
                    1,
                    a,
                    lda,
                    a + utils.index(order, n1, 0, lda),
                    lda,
                    ctx,
                ) catch unreachable;

                // Update and factor A22.
                if (comptime !types.isComplex(A)) {
                    blas.syrk(
                        order,
                        uplo,
                        .no_trans,
                        n2,
                        n1,
                        -1,
                        a + utils.index(order, n1, 0, lda),
                        lda,
                        1,
                        a + utils.index(order, n1, n1, lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                } else {
                    blas.herk(
                        order,
                        uplo,
                        .no_trans,
                        n2,
                        n1,
                        -1,
                        a + utils.index(order, n1, 0, lda),
                        lda,
                        1,
                        a + utils.index(order, n1, n1, lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                }

                iinfo = potrf2(
                    uplo,
                    n2,
                    a + utils.index(order, n1, n1, lda),
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
