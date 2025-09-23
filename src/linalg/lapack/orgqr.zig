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
const Direction = lapack.Direction;
const Storage = lapack.Storage;

const utils = @import("../utils.zig");

pub fn orgqr(
    order: Order,
    m: i32,
    n: i32,
    k: i32,
    a: anytype,
    lda: i32,
    tau: anytype,
    work: anytype,
    lwork: i32,
    ctx: anytype,
) !void {
    if (m < 0 or n < 0 or n > m or k < 0 or k > n or lda < int.max(1, if (order == .col_major) m else n) or
        (lwork < int.max(1, n) and lwork != -1))
        return lapack.Error.InvalidArgument;

    var nb = lapack.ilaenv(1, "DORGQR", " ", m, n, k, -1);
    if (lwork == -1) {
        const lwkopt: i32 = int.max(1, n) * nb;

        try ops.set(
            &work[0],
            lwkopt,
            ctx,
        );

        return;
    }

    if (n == 0) {
        try ops.set(
            &work[0],
            1,
            ctx,
        );

        return;
    }

    var nbmin: i32 = 2;
    var nx: i32 = 0;
    var iws: i32 = n;
    var ldwork: i32 = 0;
    if (nb > 1 and nb < k) {
        // Determine when to cross over from blocked to unblocked code.
        nx = int.max(0, lapack.ilaenv(3, "DORGQR", " ", m, n, k, -1));
        if (nx < k) {
            // Determine if workspace is large enough for blocked code.
            ldwork = if (order == .col_major) n else nb;
            iws = ldwork * nb;
            if (lwork < iws) {
                // Not enough workspace to use optimal nb:  reduce nb and
                // determine the minimum value of nb.
                nb = int.div(lwork, ldwork);
                nbmin = int.max(2, lapack.ilaenv(2, "DORGQR", " ", m, n, k, -1));
            }
        }
    }

    var ki: i32 = 0;
    var kk: i32 = 0;
    if (nb >= nbmin and nb < k and nx < k) {
        // Use blocked code after the last block.
        // The first kk columns are handled by the block method.
        ki = int.div(k - nx - 1, nb) * nb;
        kk = int.min(k, ki + nb);

        // Set A(1:kk,kk+1:n) to zero.
        var j: i32 = kk;
        while (j < n) : (j += 1) {
            var i: i32 = 0;
            while (i < kk) : (i += 1) {
                try ops.set(
                    &a[utils.index(order, i, j, lda)],
                    0,
                    ctx,
                );
            }
        }
    }

    // Use unblocked code for the last or only block.
    if (kk < n) {
        try lapack.org2r(
            order,
            m - kk,
            n - kk,
            k - kk,
            a + utils.index(order, kk, kk, lda),
            lda,
            tau + types.scast(u32, kk),
            work,
            ctx,
        );
    }

    if (kk > 0) {
        // Use blocked code
        var i: i32 = ki;
        while (i >= 0) : (i -= nb) {
            const ib: i32 = int.min(k - i, nb);

            if (i + ib < n) {
                // Form the triangular factor of the block reflector
                // H = H(i) H(i+1) . . . H(i+ib-1)
                try lapack.larft(
                    order,
                    .forward,
                    .columnwise,
                    m - i,
                    ib,
                    a + utils.index(order, i, i, lda),
                    lda,
                    tau + types.scast(u32, i),
                    work,
                    ldwork,
                    ctx,
                );

                // Apply H to A(i:m,i+ib:n) from the left
                try lapack.larfb(
                    order,
                    .left,
                    .no_trans,
                    .forward,
                    .columnwise,
                    m - i,
                    n - i - ib,
                    ib,
                    a + utils.index(order, i, i, lda),
                    lda,
                    work,
                    ldwork,
                    a + utils.index(order, i, i + ib, lda),
                    lda,
                    work + utils.index(order, ib, 0, ldwork),
                    ldwork,
                    ctx,
                );
            }

            // Apply H to rows i:m of current block
            try lapack.org2r(
                order,
                m - i,
                ib,
                ib,
                a + utils.index(order, i, i, lda),
                lda,
                tau + types.scast(u32, i),
                work,
                ctx,
            );

            // Set rows 1:i-1 of current block to zero
            var j: i32 = i;
            while (j < i + ib) : (j += 1) {
                var l: i32 = 0;
                while (l < i) : (l += 1) {
                    try ops.set(
                        &a[utils.index(order, l, j, lda)],
                        0,
                        ctx,
                    );
                }
            }
        }
    }

    try ops.set(
        &work[0],
        iws,
        ctx,
    );

    return;
}
