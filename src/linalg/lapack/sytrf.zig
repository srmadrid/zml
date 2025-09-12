const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Uplo = types.Uplo;

const utils = @import("utils.zig");

pub fn sytrf(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    work: anytype,
    lwork: i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));

    if (n < 0 or lda < int.max(1, n) or (lwork < 1 and lwork != -1))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    // Determine the block size
    var nb: i32 = lapack.ilaenv(1, "DSYTRF", if (uplo == .upper) "U" else "L", n, -1, -1, -1);
    const lwkopt: i32 = int.max(1, n * nb);
    try ops.set(
        &work[0],
        lwkopt,
        ctx,
    );

    if (lwork == -1)
        return info;

    var nbmin: i32 = 2;
    const ldwork: i32 = if (order == .col_major) n else nb;
    if (nb > 1 and nb < n) {
        if (lwork < ldwork * nb) {
            nb = int.max(int.div(lwork, ldwork), 1);
            nbmin = int.max(2, lapack.ilaenv(2, "DSYTRF", if (uplo == .upper) "U" else "L", n, -1, -1, -1));
        }
    }

    if (nb < nbmin)
        nb = n;

    if (comptime !types.isArbitraryPrecision(A)) {
        if (uplo == .upper) {
            // Factorize a as u**t*d*u using the upper triangle of a
            //
            // k is the main loop index, decreasing from n-1 to 0 in steps of
            // kb, where kb is the number of columns factorized by dlasyf;
            // kb is either nb or nb-1, or k+1 for the last block
            var k: i32 = n - 1;
            while (true) {
                // If k < 0, exit from loop
                if (k < 0)
                    break;

                var iinfo: i32 = 0;
                var kb: i32 = undefined;
                if (k + 1 > nb) {
                    // Factorize columns k-kb+1:k of a and use blocked code to
                    // update columns 0:k-kb
                    iinfo = try lapack.lasyf(
                        order,
                        .upper,
                        k + 1,
                        nb,
                        &kb,
                        a,
                        lda,
                        ipiv,
                        work,
                        ldwork,
                        ctx,
                    );
                } else {
                    // Use unblocked code to factorize columns 0:k of a
                    iinfo = try lapack.sytf2(
                        order,
                        .upper,
                        k + 1,
                        a,
                        lda,
                        ipiv,
                        ctx,
                    );

                    kb = k + 1;
                }

                // Set info on the first occurrence of a zero pivot
                if (info == 0 and iinfo > 0)
                    info = iinfo;

                // Decrease k and return to the start of the main loop
                k -= kb;
            }
        } else {
            // factorize a as l*d*l**t using the lower triangle of a
            //
            // k is the main loop index, increasing from 0 to n-1 in steps of
            // kb, where kb is the number of columns factorized by dlasyf;
            // kb is either nb or nb-1, or n-k for the last block

            var k: i32 = 0;
            while (true) {
                // Ff k >= n, exit from loop
                if (k >= n)
                    break;

                var iinfo: i32 = 0;
                var kb: i32 = undefined;
                if (k <= n - nb - 1) {
                    std.debug.print("Calling lasyf with k = {d}, n - k = {d}\n", .{ k, n - k });
                    // Factorize columns k:k+kb-1 of a and use blocked code to
                    // update columns k+kb:n-1
                    iinfo = try lapack.lasyf(
                        order,
                        .lower,
                        n - k,
                        nb,
                        &kb,
                        a + utils.index(order, k, k, lda),
                        lda,
                        ipiv + types.scast(u32, k),
                        work,
                        ldwork,
                        ctx,
                    );
                } else {
                    std.debug.print("Calling sytf2 with k = {d}, n - k = {d}\n", .{ k, n - k });
                    // use unblocked code to factorize columns k:n-1 of a
                    iinfo = try lapack.sytf2(
                        order,
                        .lower,
                        n - k,
                        a + utils.index(order, k, k, lda),
                        lda,
                        ipiv + types.scast(u32, k),
                        ctx,
                    );

                    kb = n - k;
                }

                // Set info on the first occurrence of a zero pivot
                if (info == 0 and iinfo > 0)
                    info = iinfo + k;

                // Adjust ipiv
                var j: i32 = k;
                while (j < k + kb) : (j += 1) {
                    if (ipiv[types.scast(u32, j)] > 0) {
                        ipiv[types.scast(u32, j)] = ipiv[types.scast(u32, j)] + k;
                    } else {
                        ipiv[types.scast(u32, j)] = ipiv[types.scast(u32, j)] - k;
                    }
                }

                // increase k and return to the start of the main loop
                k += kb;
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.sytrf not implemented for arbitrary precision types yet");
    }

    try ops.set(
        &work[0],
        lwkopt,
        ctx,
    );
    return info;
}
