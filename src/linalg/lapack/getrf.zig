const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;

const utils = @import("../utils.zig");

pub fn getrf(
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

    // Determine the block size for this environment. Always returns the
    // same number for getrf regardless of 'S', 'D', 'C', or 'Z'.
    const nb: i32 = lapack.ilaenv(1, "DGETRF", " ", m, n, -1, -1);
    if (comptime !types.isArbitraryPrecision(A)) {
        if (nb <= 1 or nb >= int.min(m, n)) {
            // Use unblocked code.
            info = lapack.getrf2(
                order,
                m,
                n,
                a,
                lda,
                ipiv,
                ctx,
            ) catch unreachable;
        } else {
            // Use blocked code.
            var j: i32 = 0;
            while (j < int.min(m, n)) : (j += nb) {
                const jb: i32 = int.min(int.min(m, n) - j, nb);

                // Factor diagonal and subdiagonal blocks and test for exact singularity.
                //const iinfo: i32 = lapack.getrf2(
                const iinfo: i32 = lapack.getrf2(
                    order,
                    m - j,
                    jb,
                    a + utils.index(order, j, j, lda),
                    lda,
                    ipiv + types.scast(u32, j),
                    ctx,
                ) catch unreachable;

                // Adjust info and the pivot indices.
                if (info == 0 and iinfo > 0)
                    info = iinfo + j;

                var i: i32 = j;
                while (i < int.min(m, j + jb)) : (i += 1) {
                    ops.add_( // ipiv[i] += j;
                        &ipiv[types.scast(u32, i)],
                        ipiv[types.scast(u32, i)],
                        j,
                        ctx,
                    ) catch unreachable;
                }

                // Apply interchanges to columns 1:j - 1.
                lapack.laswp(
                    order,
                    j,
                    a,
                    lda,
                    j + 1,
                    j + jb,
                    ipiv,
                    1,
                ) catch unreachable;

                if (j + jb < n) {
                    // Apply interchanges to columns j + jb:n.
                    lapack.laswp(
                        order,
                        n - j - jb,
                        a + utils.index(order, 0, j + jb, lda),
                        lda,
                        j + 1,
                        j + jb,
                        ipiv,
                        1,
                    ) catch unreachable;

                    // Compute block row of U.
                    blas.trsm(
                        order,
                        .left,
                        .lower,
                        .no_trans,
                        .unit,
                        jb,
                        n - j - jb,
                        1,
                        a + utils.index(order, j, j, lda),
                        lda,
                        a + utils.index(order, j, j + jb, lda),
                        lda,
                        ctx,
                    ) catch unreachable;

                    if (j + jb < m) {
                        // Update trailing submatrix.
                        blas.gemm(
                            order,
                            .no_trans,
                            .no_trans,
                            m - j - jb,
                            n - j - jb,
                            jb,
                            -1,
                            a + utils.index(order, j + jb, j, lda),
                            lda,
                            a + utils.index(order, j, j + jb, lda),
                            lda,
                            1,
                            a + utils.index(order, j + jb, j + jb, lda),
                            lda,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.getrf not implemented for arbitrary precision types yet");
    }

    return info;
}
