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

pub inline fn getrf2(
    order: Order,
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    ctx: anytype,
) !i32 {
    if (order == .col_major) {
        return k_getrf2_c(m, n, a, lda, ipiv, ctx);
    } else {
        return k_getrf2_r(m, n, a, lda, ipiv, ctx);
    }
}

fn k_getrf2_c(
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));

    if (m < 0 or n < 0 or lda < int.max(1, m))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    // Quick return if possible.
    if (m == 0 or n == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(A)) {
        if (m == 1) {
            // Use unblocked code for one row case. Just need to handle ipiv and
            // info.
            ipiv[0] = 1;
            if (ops.eq(a[0], 0, ctx) catch unreachable)
                info = 1;
        } else if (n == 1) {
            // Use unblocked code for one column case.

            // Compute machine safe minimum.
            const sfmin: Scalar(A) = lapack.lamch(Scalar(A), .sfmin);

            // Find pivot and test for singularity.
            const ii: u32 = blas.iamax(m, a, 1, ctx) catch unreachable;
            ipiv[0] = scast(i32, ii) + 1;

            if (ops.ne(a[ii], 0, ctx) catch unreachable) {
                // Apply the interchange.
                if (ii != 0) {
                    const temp: A = a[0];
                    a[0] = a[ii];
                    a[ii] = temp;
                }

                // Compute elements 1:m - 1 of the column.
                if (ops.ge(ops.abs(a[0], ctx) catch unreachable, sfmin, ctx) catch unreachable) {
                    blas.scal(
                        m - 1,
                        ops.div(
                            1,
                            a[0],
                            ctx,
                        ) catch unreachable,
                        a + 1,
                        1,
                        ctx,
                    ) catch unreachable;
                } else {
                    var i: i32 = 0;
                    while (i < m - 1) : (i += 1) {
                        ops.div_( // a[i + 1] /= a[0]
                            &a[scast(u32, i + 1)],
                            a[scast(u32, i + 1)],
                            a[0],
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                info = 1;
            }
        } else {
            // Use recursive code.
            const n1: i32 = int.div(int.min(m, n), 2);
            const n2: i32 = n - n1;

            //        [ A11 ]
            // Factor [ --- ]
            //        [ A21 ]
            var iinfo: i32 = k_getrf2_c(
                m,
                n1,
                a,
                lda,
                ipiv,
                ctx,
            ) catch unreachable;

            if (info == 0 and iinfo > 0)
                info = iinfo;

            //                       [ A12 ]
            // Apply interchanges to [ --- ]
            //                       [ A22 ]
            lapack.laswp(
                .col_major,
                n2,
                a + scast(u32, n1 * lda),
                lda,
                1,
                n1,
                ipiv,
                1,
            ) catch unreachable;

            // Solve A12.
            blas.trsm(
                .col_major,
                .left,
                .lower,
                .no_trans,
                .unit,
                n1,
                n2,
                1,
                a,
                lda,
                a + scast(u32, n1 * lda),
                lda,
                ctx,
            ) catch unreachable;

            // Update A22.
            blas.gemm(
                .col_major,
                .no_trans,
                .no_trans,
                m - n1,
                n2,
                n1,
                -1,
                a + scast(u32, n1),
                lda,
                a + scast(u32, n1 * lda),
                lda,
                1,
                a + scast(u32, n1 + n1 * lda),
                lda,
                ctx,
            ) catch unreachable;

            // Factor A22.
            iinfo = k_getrf2_c(
                m - n1,
                n2,
                a + scast(u32, n1 + n1 * lda),
                lda,
                ipiv + scast(u32, n1),
                ctx,
            ) catch unreachable;

            // Adjust info and the pivot indices.
            if (info == 0 and iinfo > 0)
                ops.add_( // info = iinfo + n1
                    &info,
                    iinfo,
                    n1,
                    ctx,
                ) catch unreachable;

            var i: i32 = n1;
            while (i < int.min(m, n)) : (i += 1) {
                ops.add_( // ipiv[i] += n1
                    &ipiv[scast(u32, i)],
                    ipiv[scast(u32, i)],
                    n1,
                    ctx,
                ) catch unreachable;
            }

            // Apply interchanges to A21.
            lapack.laswp(
                .col_major,
                n1,
                a,
                lda,
                n1 + 1,
                int.min(m, n),
                ipiv,
                1,
            ) catch unreachable;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.getrf2 not implemented for arbitrary precision types yet");
    }

    return info;
}

fn k_getrf2_r(
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));

    if (m < 0 or n < 0 or lda < int.max(1, n))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    // Quick return if possible.
    if (m == 0 or n == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(A)) {
        if (m == 1) {
            // Use unblocked code for one row case. Just need to handle ipiv and
            // info.
            ipiv[0] = 1;
            if (ops.eq(a[0], 0, ctx) catch unreachable)
                info = 1;
        } else if (n == 1) {
            // Use unblocked code for one column case.

            // Compute machine safe minimum.
            const sfmin: Scalar(A) = lapack.lamch(Scalar(A), .sfmin);

            // Find pivot and test for singularity.
            const ii: u32 = blas.iamax(m, a, lda, ctx) catch unreachable;
            ipiv[0] = scast(i32, ii) + 1;

            if (ops.ne(a[scast(u32, int.mul(ii, lda, .wrap))], 0, ctx) catch unreachable) {
                // Apply the interchange.
                if (ii != 0) {
                    const temp: A = a[0];
                    a[0] = a[scast(u32, int.mul(ii, lda, .wrap))];
                    a[scast(u32, int.mul(ii, lda, .wrap))] = temp;
                }

                // Compute elements 1:m - 1 of the column.
                if (ops.ge(ops.abs(a[0], ctx) catch unreachable, sfmin, ctx) catch unreachable) {
                    blas.scal(
                        m - 1,
                        ops.div(
                            1,
                            a[0],
                            ctx,
                        ) catch unreachable,
                        a + scast(u32, lda),
                        lda,
                        ctx,
                    ) catch unreachable;
                } else {
                    var i: i32 = 0;
                    while (i < m - 1) : (i += 1) {
                        ops.div_( // a[i + 1] /= a[0]
                            &a[scast(u32, (i + 1) * lda)],
                            a[scast(u32, (i + 1) * lda)],
                            a[0],
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else {
                info = 1;
            }
        } else {
            // Use recursive code.
            const n1: i32 = int.div(int.min(m, n), 2);
            const n2: i32 = n - n1;

            //        [ A11 ]
            // Factor [ --- ]
            //        [ A21 ]
            var iinfo: i32 = k_getrf2_r(
                m,
                n1,
                a,
                lda,
                ipiv,
                ctx,
            ) catch unreachable;

            if (info == 0 and iinfo > 0)
                info = iinfo;

            //                       [ A12 ]
            // Apply interchanges to [ --- ]
            //                       [ A22 ]
            lapack.laswp(
                .row_major,
                n2,
                a + scast(u32, n1),
                lda,
                1,
                n1,
                ipiv,
                1,
            ) catch unreachable;

            // Solve A12.
            blas.trsm(
                .row_major,
                .left,
                .lower,
                .no_trans,
                .unit,
                n1,
                n2,
                1,
                a,
                lda,
                a + scast(u32, n1),
                lda,
                ctx,
            ) catch unreachable;

            // Update A22.
            blas.gemm(
                .row_major,
                .no_trans,
                .no_trans,
                m - n1,
                n2,
                n1,
                -1,
                a + scast(u32, n1 * lda),
                lda,
                a + scast(u32, n1),
                lda,
                1,
                a + scast(u32, n1 + n1 * lda),
                lda,
                ctx,
            ) catch unreachable;

            // Factor A22.
            iinfo = k_getrf2_r(
                m - n1,
                n2,
                a + scast(u32, n1 + n1 * lda),
                lda,
                ipiv + scast(u32, n1),
                ctx,
            ) catch unreachable;

            // Adjust info and the pivot indices.
            if (info == 0 and iinfo > 0)
                ops.add_( // info = iinfo + n1
                    &info,
                    iinfo,
                    n1,
                    ctx,
                ) catch unreachable;

            var i: i32 = n1;
            while (i < int.min(m, n)) : (i += 1) {
                ops.add_( // ipiv[i] += n1
                    &ipiv[scast(u32, i)],
                    ipiv[scast(u32, i)],
                    n1,
                    ctx,
                ) catch unreachable;
            }

            // Apply interchanges to A21.
            lapack.laswp(
                .row_major,
                n1,
                a,
                lda,
                n1 + 1,
                int.min(m, n),
                ipiv,
                1,
            ) catch unreachable;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.getrf2 not implemented for arbitrary precision types yet");
    }

    return info;
}
