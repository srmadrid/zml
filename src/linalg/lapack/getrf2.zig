const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;

const utils = @import("../utils.zig");

pub fn getrf2(
    order: Order,
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));

    if (m < 0 or n < 0 or lda < int.max(1, if (order == .col_major) m else n))
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
            const sfmin: types.Scalar(A) = lapack.lamch(types.Scalar(A), .sfmin);

            // Find pivot and test for singularity.
            const ii: i32 = types.scast(i32, blas.iamax(
                m,
                a,
                utils.col_ld(order, lda),
                ctx,
            ) catch unreachable);
            ipiv[0] = ii + 1;

            if (ops.ne(a[utils.index(order, ii, 0, lda)], 0, ctx) catch unreachable) {
                // Apply the interchange.
                if (ii != 0) {
                    const temp: A = a[0];
                    a[0] = a[utils.index(order, ii, 0, lda)];
                    a[utils.index(order, ii, 0, lda)] = temp;
                }

                // Compute elements 1:m - 1 of the column.
                if (ops.ge(ops.abs1(a[0], ctx) catch unreachable, sfmin, ctx) catch unreachable) {
                    blas.scal(
                        m - 1,
                        ops.div(
                            1,
                            a[0],
                            ctx,
                        ) catch unreachable,
                        a + utils.index(order, 1, 0, lda),
                        utils.col_ld(order, lda),
                        ctx,
                    ) catch unreachable;
                } else {
                    var i: i32 = 0;
                    while (i < m - 1) : (i += 1) {
                        ops.div_( // a[i + 1] /= a[0]
                            &a[utils.index(order, i + 1, 0, lda)],
                            a[utils.index(order, i + 1, 0, lda)],
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
            var iinfo: i32 = getrf2(
                order,
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
                order,
                n2,
                a + utils.index(order, 0, n1, lda),
                lda,
                1,
                n1,
                ipiv,
                1,
            ) catch unreachable;

            // Solve A12.
            blas.trsm(
                order,
                .left,
                .lower,
                .no_trans,
                .unit,
                n1,
                n2,
                1,
                a,
                lda,
                a + utils.index(order, 0, n1, lda),
                lda,
                ctx,
            ) catch unreachable;

            // Update A22.
            blas.gemm(
                order,
                .no_trans,
                .no_trans,
                m - n1,
                n2,
                n1,
                -1,
                a + utils.index(order, n1, 0, lda),
                lda,
                a + utils.index(order, 0, n1, lda),
                lda,
                1,
                a + utils.index(order, n1, n1, lda),
                lda,
                ctx,
            ) catch unreachable;

            // Factor A22.
            iinfo = getrf2(
                order,
                m - n1,
                n2,
                a + utils.index(order, n1, n1, lda),
                lda,
                ipiv + types.scast(u32, n1),
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
                    &ipiv[types.scast(u32, i)],
                    ipiv[types.scast(u32, i)],
                    n1,
                    ctx,
                ) catch unreachable;
            }

            // Apply interchanges to A21.
            lapack.laswp(
                order,
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
